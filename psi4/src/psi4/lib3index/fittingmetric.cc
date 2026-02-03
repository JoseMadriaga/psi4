/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "3index.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/psimath.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/wavefunction.h"
#include <iostream>

// MKL Header
#ifdef USING_LAPACK_MKL
#include <mkl.h>
#endif

// OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

namespace psi {

FittingMetric::FittingMetric(std::shared_ptr<BasisSet> aux, bool force_C1)
    : aux_(aux), is_poisson_(false), is_inverted_(false), force_C1_(force_C1), omega_(0.0) {}
FittingMetric::FittingMetric(std::shared_ptr<BasisSet> aux, double omega, bool force_C1)
    : aux_(aux), is_poisson_(false), is_inverted_(false), force_C1_(force_C1), omega_(omega) {}
FittingMetric::FittingMetric(std::shared_ptr<BasisSet> aux, std::shared_ptr<BasisSet> pois, bool force_C1)
    : aux_(aux), pois_(pois), is_poisson_(true), is_inverted_(false), force_C1_(force_C1), omega_(0.0) {}

FittingMetric::~FittingMetric() {}

void FittingMetric::form_fitting_metric() {
    is_inverted_ = false;
    algorithm_ = "NONE";

    // Sizing/symmetry indexing
    auto auxfact = std::make_shared<IntegralFactory>(aux_, aux_, aux_, aux_);
    auto auxpet = std::make_shared<PetiteList>(aux_, auxfact);
    std::shared_ptr<IntegralFactory> poisfact;
    std::shared_ptr<PetiteList> poispet;
    if (is_poisson_) {
        poisfact = std::make_shared<IntegralFactory>(pois_, pois_, pois_, pois_);
        poispet = std::make_shared<PetiteList>(pois_, poisfact);
    }

    int naux = 0;
    int ngaussian = 0;
    int npoisson = 0;
    Dimension nauxpi(auxpet->nirrep(), "Fitting Metric Dimensions");
    for (int h = 0; h < auxpet->nirrep(); h++) {
        naux += auxpet->SO_basisdim()[h];
        ngaussian += auxpet->SO_basisdim()[h];
        nauxpi[h] = auxpet->SO_basisdim()[h];
        if (is_poisson_) {
            naux += poispet->SO_basisdim()[h];
            npoisson += poispet->SO_basisdim()[h];
            nauxpi[h] += poispet->SO_basisdim()[h];
        }
    }
    Dimension ngauspi = auxpet->SO_basisdim();

    // Build the full DF/Poisson matrix in the AO basis first
    auto AOmetric = std::make_shared<Matrix>("AO Basis DF Metric", naux, naux);
    double** W = AOmetric->pointer(0);
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    // Only thread if not already in parallel (handy for local fitting)
    int nthread = 1;
#ifdef _OPENMP
    if (!omp_in_parallel()) {
        nthread = Process::environment.get_n_threads();
    }
#endif

    // == (A|B) Block == //
    IntegralFactory rifactory_J(aux_, zero, aux_, zero);
    std::vector<const double*> Jbuffer(nthread);  
    std::vector<std::shared_ptr<TwoBodyAOInt>> Jint(nthread);
    for (int Q = 0; Q < nthread; Q++) {
        if (omega_ > 0.0) {
            Jint[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory_J.erf_eri(omega_));
        } else {
            Jint[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory_J.eri());
        }
        if (!Jint[Q]->sieve_initialized()) Jint[Q]->initialize_sieve();
    }

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (int MU = 0; MU < aux_->nshell(); ++MU) {
        int nummu = aux_->shell(MU).nfunction();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (int NU = 0; NU <= MU; ++NU) {
            int numnu = aux_->shell(NU).nfunction();

            Jint[thread]->compute_shell(MU, 0, NU, 0);
            Jbuffer[thread] = Jint[thread]->buffer();

            int index = 0;
            for (int mu = 0; mu < nummu; ++mu) {
                int omu = aux_->shell(MU).function_index() + mu;

                for (int nu = 0; nu < numnu; ++nu, ++index) {
                    int onu = aux_->shell(NU).function_index() + nu;

                    W[omu][onu] = Jbuffer[thread][index];
                    W[onu][omu] = Jbuffer[thread][index];
                }
            }
        }
    }

    if (is_poisson_) {
        // == (AB) Block == //
        IntegralFactory rifactory_RP(pois_, aux_, zero, zero);
        std::vector<const double*> Obuffer(nthread);
        std::vector<std::shared_ptr<OneBodyAOInt>> Oint(nthread);
        for (int Q = 0; Q < nthread; Q++) {
            Oint[Q] = std::shared_ptr<OneBodyAOInt>(rifactory_RP.ao_overlap());
            Obuffer[Q] = Oint[Q]->buffers()[0];
        }

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (int NU = 0; NU < pois_->nshell(); ++NU) {
            int numnu = pois_->shell(NU).nfunction();

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            for (int MU = 0; MU < aux_->nshell(); ++MU) {
                int nummu = aux_->shell(MU).nfunction();

                Oint[thread]->compute_shell(NU, MU);

                int index = 0;
                for (int nu = 0; nu < numnu; ++nu) {
                    int onu = pois_->shell(NU).function_index() + nu;

                    for (int mu = 0; mu < nummu; ++mu, ++index) {
                        int omu = aux_->shell(MU).function_index() + mu;

                        W[omu][onu + ngaussian] = Obuffer[thread][index];
                        W[onu + ngaussian][omu] = Obuffer[thread][index];
                    }
                }
            }
        }

        // == (A|T|B) Block == //
        IntegralFactory rifactory_P(pois_, pois_, zero, zero);
        std::vector<const double*> Tbuffer(nthread);
        std::vector<std::shared_ptr<OneBodyAOInt>> Tint(nthread);
        for (int Q = 0; Q < nthread; Q++) {
            Tint[Q] = std::shared_ptr<OneBodyAOInt>(rifactory_P.ao_kinetic());
            Tbuffer[Q] = Tint[Q]->buffers()[0];
        }

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (int MU = 0; MU < pois_->nshell(); ++MU) {
            int nummu = pois_->shell(MU).nfunction();

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            for (int NU = 0; NU <= MU; ++NU) {
                int numnu = pois_->shell(NU).nfunction();

                Tint[thread]->compute_shell(MU, NU);

                int index = 0;
                for (int mu = 0; mu < nummu; ++mu) {
                    int omu = pois_->shell(MU).function_index() + mu;

                    for (int nu = 0; nu < numnu; ++nu, ++index) {
                        int onu = pois_->shell(NU).function_index() + nu;

                        // These integrals are (A | -1/2 \nabla^2 | B), and should be (A | - 1 / (4 * pi) \nabla^2 |B)
                        // So a factor of 1 / (2 * PI) is the difference
                        W[omu + ngaussian][onu + ngaussian] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                        W[onu + ngaussian][omu + ngaussian] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                    }
                }
            }
        }
    }

    // If C1, form indexing and exit immediately (multiplying by 1 is not so gratifying)
    if (auxpet->nirrep() == 1 || force_C1_ == true) {
        metric_ = AOmetric;
        if (DPC_ or TR_) {
            const double regularizer = 1e-10;
            outfile->Printf("Adding in the regularizer: %.8e\n", regularizer);
            for (int i = 0; i < naux; ++i) {
                metric_->add(i, i, regularizer);
            }
        } 
        metric_->set_name("SO Basis Fitting Metric");
        pivots_ = std::make_shared<IntVector>(naux);
        rev_pivots_ = std::make_shared<IntVector>(naux);
        int* piv = pivots_->pointer();
        int* rpiv = pivots_->pointer();
        for (int Q = 0; Q < naux; Q++) {
            piv[Q] = Q;
            rpiv[Q] = Q;
        }
        return;
    }

    // Get the similarity transform objects
    SharedMatrix auxAO2USO(auxpet->sotoao());
    // auxAO2USO->print();
    SharedMatrix poisAO2USO;
    if (is_poisson_) {
        poisAO2USO = SharedMatrix(poispet->sotoao());
        // poisAO2USO->print();
    }

    // Allocate the fitting metric
    metric_ = std::make_shared<Matrix>("SO Basis Fitting Metric", nauxpi, nauxpi);
    SharedMatrix Temp;
    double** Temp1;

    // Transform AO to SO
    for (int h = 0; h < auxpet->nirrep(); h++) {
        // Gaussian-Gaussian part
        double** J = metric_->pointer(h);
        double** auxU = auxAO2USO->pointer(h);

        if (ngauspi[h] != 0) {
            Temp = std::make_shared<Matrix>("Temp", ngauspi[h], ngaussian);
            Temp1 = Temp->pointer();
            C_DGEMM('N', 'N', ngauspi[h], ngaussian, ngaussian, 1.0, auxU[0], ngaussian, W[0], naux, 0.0, Temp1[0],
                    ngaussian);
            C_DGEMM('N', 'T', ngauspi[h], ngauspi[h], ngaussian, 1.0, Temp1[0], ngaussian, auxU[0], ngaussian, 0.0,
                    J[0], nauxpi[h]);
            Temp.reset();
        }

        if (is_poisson_ && poispet->SO_basisdim()[h] != 0) {
            Dimension npoispi = poispet->SO_basisdim();
            double** poisU = poisAO2USO->pointer(h);

            // Gaussian-Poisson part
            if (ngauspi[h] != 0) {
                Temp = std::make_shared<Matrix>("Temp", ngauspi[h], npoisson);
                Temp1 = Temp->pointer();
                C_DGEMM('N', 'N', ngauspi[h], npoisson, ngaussian, 1.0, auxU[0], ngaussian, &W[0][ngaussian], naux, 0.0,
                        Temp1[0], npoisson);
                C_DGEMM('N', 'T', ngauspi[h], npoispi[h], npoisson, 1.0, Temp1[0], npoisson, poisU[0], npoisson, 0.0,
                        &J[0][ngauspi[h]], nauxpi[h]);
                for (int Q = 0; Q < ngauspi[h]; Q++)
                    for (int P = 0; P < npoispi[h]; P++) J[P + ngauspi[h]][Q] = J[Q][P + ngauspi[h]];
                Temp.reset();
            }

            // Poisson-Poisson part
            size_t AOoffset = ngaussian * (size_t)naux + (size_t)ngaussian;
            size_t SOoffset = ngauspi[h] * (size_t)nauxpi[h] + (size_t)ngauspi[h];
            Temp = std::make_shared<Matrix>("Temp", npoispi[h], npoisson);
            Temp1 = Temp->pointer();
            C_DGEMM('N', 'N', npoispi[h], npoisson, npoisson, 1.0, poisU[0], npoisson, &W[0][AOoffset], naux, 0.0,
                    Temp1[0], npoisson);
            C_DGEMM('N', 'T', npoispi[h], npoispi[h], npoisson, 1.0, Temp1[0], npoisson, poisU[0], npoisson, 0.0,
                    &J[0][SOoffset], nauxpi[h]);
            Temp.reset();
        }
    }

    // Form indexing
    pivots_ = std::make_shared<IntVector>(nauxpi);
    rev_pivots_ = std::make_shared<IntVector>(nauxpi);
    for (int h = 0; h < auxpet->nirrep(); h++) {
        int* piv = pivots_->pointer(h);
        int* rpiv = pivots_->pointer(h);
        for (int Q = 0; Q < nauxpi[h]; Q++) {
            piv[Q] = Q;
            rpiv[Q] = Q;
        }
    }
}
void FittingMetric::form_cholesky_inverse() {
    is_inverted_ = true;
    algorithm_ = "CHOLESKY";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {
        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        int info = C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
        for (int A = 0; A < metric_->colspi()[h]; A++)
            for (int B = 0; B < A; B++) J[A][B] = 0.0;
    }
    metric_->set_name("SO Basis Fitting Inverse (Cholesky)");
}
void FittingMetric::form_QR_inverse(double tol) {
    is_inverted_ = true;
    algorithm_ = "QR";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {
        if (metric_->colspi()[h] == 0) continue;

        //        metric_->print();

        double** J = metric_->pointer(h);
        int n = metric_->colspi()[h];

        // Copy the J matrix to R (actually R')
        auto R = std::make_shared<Matrix>("R", n, n);
        double** Rp = R->pointer();
        C_DCOPY(n * (size_t)n, J[0], 1, Rp[0], 1);

        // QR Decomposition
        std::vector<double> tau(n);

        // First, find out how much workspace to provide
        // Optimal size of work vector is written to work_size
        double work_size;
        C_DGEQRF(n, n, Rp[0], n, tau.data(), &work_size, -1);

        // Now, do the QR decomposition
        int lwork = (int)work_size;
        std::vector<double> work(lwork);
        C_DGEQRF(n, n, Rp[0], n, tau.data(), work.data(), lwork);

        // Copy Jcopy to Q (actually Q')
        auto Q = std::make_shared<Matrix>("Q", n, n);
        double** Qp = Q->pointer();
        C_DCOPY(n * (size_t)n, Rp[0], 1, Qp[0], 1);

        // Put R in the upper triangle where it belongs
        for (int i = 1; i < n; i++)
            for (int j = 0; j < i; j++) {
                Rp[j][i] = 0.0;
            }

        // First, find out how much workspace to provide
        // Optimal size of work vector is written to work_size
        C_DORGQR(n, n, n, Qp[0], n, tau.data(), &work_size, -1);

        // Now, form Q
        lwork = (int)work_size;
        work.resize(lwork);
        C_DORGQR(n, n, n, Qp[0], n, tau.data(), work.data(), lwork);

        // Q->print();
        // R->print();

        // Find the number of significant basis functions
        int nsig = 0;
        double R_max = std::fabs(Rp[0][0]);
        for (int A = 0; A < n; A++) {
            if ((std::fabs(Rp[A][A]) / R_max) < tol) break;
            nsig++;
        }

        // Transform into the reduced basis
        // Just use R's memory, don't need it anymore
        C_DGEMM('N', 'N', nsig, n, n, 1.0, Qp[0], n, J[0], n, 0.0, Rp[0], n);
        C_DGEMM('N', 'T', nsig, nsig, n, 1.0, Rp[0], n, Qp[0], n, 0.0, J[0], nsig);

        // Find the Cholesky factor in the reduced basis
        C_DPOTRF('L', nsig, J[0], nsig);

        // Backsolve the triangular factor against the change of basis matrix
        C_DTRSM('L', 'U', 'N', 'N', nsig, n, 1.0, J[0], nsig, Qp[0], n);

        // Zero out the metric
        memset(static_cast<void*>(J[0]), '\0', n * (size_t)n);

        // Copy the top bit in
        C_DCOPY(n * (size_t)nsig, Qp[0], 1, J[0], 1);

    }
    metric_->set_name("SO Basis Fitting Inverse (QR)");
}
void FittingMetric::form_eig_inverse(double tol) {
    is_inverted_ = true;
    algorithm_ = "EIG";

    form_fitting_metric();
    metric_->power(-0.5, tol);
    metric_->set_name("SO Basis Fitting Inverse (Eig)");
}
void FittingMetric::form_eig_inverse_DPC() {
    // --- Setup DPC ---
    is_inverted_ = true;
    algorithm_ = "DPC";
    DPC_ = true;
    
    std::cout << "prior to form_fitting_metric()" << std::endl;
    form_fitting_metric();
    
    std::cout << "prior to allocating zero, basis, D_ao, naux, nbf" << std::endl;
    auto zero = BasisSet::zero_ao_basis_set();
    auto basis = reference_wavefunction_->basisset();
    //auto D_ao = reference_wavefunction_->Da();
    //int naux = aux_->nbf();
    //int nbf  = basis->nbf();

    auto eri_fact = std::make_shared<IntegralFactory>(aux_, zero, basis, basis);
    auto eri = std::shared_ptr<TwoBodyAOInt>(eri_fact->eri());

    std::cout << "prior to generate RHS" << std::endl;
    
    SharedMatrix D_ao = reference_wavefunction_->Da();
    if (!D_ao) {
        throw std::runtime_error("D_ao is null");
    }
    
    int naux = aux_->nbf();
    int nbf  = basis->nbf();
    
    if (D_ao->nrow() != nbf || D_ao->ncol() != nbf) {
        throw std::runtime_error("D_ao dimension mismatch");
    }
    
    std::vector<double> rhs(naux, 0.0);
    
    // zero basis info
    int n0 = zero->shell(0).nfunction();
    std::cout << "zero basis nfunction = " << n0 << std::endl;
    if (n0 <= 0) {
        throw std::runtime_error("Zero basis has invalid nfunction");
    }
    
    for (int P = 0; P < aux_->nshell(); P++) {
    
        const auto& Pshell = aux_->shell(P);
        int np = Pshell.nfunction();
        int pstart = Pshell.function_index();
    
        if (pstart < 0 || pstart + np > naux) {
            throw std::runtime_error("Aux shell index out of bounds");
        }
    
        std::cout << "P shell " << P
                  << " np=" << np
                  << " pstart=" << pstart << std::endl;
    
        for (int M = 0; M < basis->nshell(); M++) {
    
            const auto& Mshell = basis->shell(M);
            int nm = Mshell.nfunction();
            int mstart = Mshell.function_index();
    
            if (mstart < 0 || mstart + nm > nbf) {
                throw std::runtime_error("Basis shell M index out of bounds");
            }
    
            for (int N = 0; N < basis->nshell(); N++) {
    
                const auto& Nshell = basis->shell(N);
                int nn = Nshell.nfunction();
                int nstart = Nshell.function_index();
    
                if (nstart < 0 || nstart + nn > nbf) {
                    throw std::runtime_error("Basis shell N index out of bounds");
                }
    
                // Compute (P | 0 M N)
                eri->compute_shell(P, 0, M, N);
                const double* buffer = eri->buffer();
    
                if (!buffer) {
                    throw std::runtime_error("ERI buffer is null");
                }
    
                int expected_size = np * n0 * nm * nn;
                int index = 0;
    
                for (int p = 0; p < np; p++) {
                    for (int q = 0; q < n0; q++) {
                        for (int m = 0; m < nm; m++) {
                            for (int n = 0; n < nn; n++) {
    
                                int Pidx = p + pstart;
                                int midx = m + mstart;
                                int nidx = n + nstart;
    
                                rhs[Pidx] += buffer[index]
                                             * (*D_ao)(midx, nidx);
                                index++;
                            }
                        }
                    }
                }
    
                if (index != expected_size) {
                    throw std::runtime_error("ERI buffer index mismatch");
                }
            }
        }
    }
    
    std::cout << "RHS generation completed successfully" << std::endl;


    //std::cout << "prior to generate three-centered and two body integrals" << std::endl;
    //// --- Compute RHS on-the-fly (memory efficient) ---
    //std::vector<double> rhs(naux, 0.0);
    //auto eri_fact = std::make_shared<IntegralFactory>(aux_, zero, basis, basis);
    //auto eri = std::shared_ptr<TwoBodyAOInt>(eri_fact->eri());
    //
    //std::cout << "prior to generate RHS" << std::endl;
    //for (int P = 0; P < aux_->nshell(); P++) {
    //    int np = aux_->shell(P).nfunction();
    //    int pstart = aux_->shell(P).function_index();
    //
    //    for (int M = 0; M < basis->nshell(); M++) {
    //        int nm = basis->shell(M).nfunction();
    //        int mstart = basis->shell(M).function_index();
    //
    //        for (int N = 0; N < basis->nshell(); N++) {
    //            int nn = basis->shell(N).nfunction();
    //            int nstart = basis->shell(N).function_index();
    //
    //            eri->compute_shell(P, 0, M, N);
    //            const double* buffer = eri->buffer();
    //
    //            // accumulate directly into rhs
    //            for (int p = 0, index = 0; p < np; p++) {
    //                for (int m = 0; m < nm; m++) {
    //                    for (int n = 0; n < nn; n++, index++) {
    //                        rhs[p + pstart] += buffer[index] * (*D_ao)(m + mstart, n + nstart);
    //                    }
    //                }
    //            }
    //        }
    //    }
    //} 
    std::cout << "prior to T_flat" << std::endl;
    // --- Flatten metric into contiguous vector for diagonalization ---
    std::vector<double> metric_flat(naux * naux, 0.0);
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            metric_flat[i * naux + j] = (*metric_)(i,j);
    
    std::cout << "prior to diagonalize" << std::endl;
    // --- Diagonalize metric (C_DSYEV requires contiguous memory) ---
    std::vector<double> eigval(naux, 0.0);
    int lwork = naux * 3;
    std::vector<double> work(lwork, 0.0);
    
    int stat = C_DSYEV('v', 'u', naux, metric_flat.data(), naux, eigval.data(), work.data(), lwork);
    if (stat != 0)
        throw std::runtime_error("C_DSYEV failed to diagonalize metric");
    
    work.clear();  // free workspace
    
    std::cout << "prior to picard_coef" << std::endl;
    // --- Compute Picard coefficients ---
    std::vector<double> sigma(naux, 0.0);
    std::vector<double> picard(naux, 0.0);
    for (int i = 0; i < naux; i++) {
        sigma[i] = std::sqrt(std::abs(eigval[i]));
    
        double dotprod = 0.0;
        for (int P = 0; P < naux; P++)
            dotprod += metric_flat[P + i * naux] * rhs[P];  // use flat eigenvectors
    
        picard[i] = std::abs(dotprod) / sigma[i];
    }
    
    rhs.clear();  // free memory
    
    std::cout << "prior to determining eps_opt" << std::endl;
    // --- Find knee of Picard coefficients ---
    double epsilon_opt = 0.0;
    for (int i = 0; i < naux - 1; i++) {
        if (picard[i+1] > picard[i]) {
            epsilon_opt = sigma[i];
            break;
        }
    }
    double tol = epsilon_opt * epsilon_opt;
    std::cout << "tol " << tol << std::endl;
    
    picard.clear();
    sigma.clear();  // free memory
    
    std::cout << "prior to inverse metric" << std::endl;
    // --- Reconstruct inverse metric ---
    for (int r = 0; r < naux; r++)
        for (int c = 0; c < naux; c++)
            (*metric_)(r,c) = 0.0;
    
    for (int i = 0; i < naux; i++) {
        if (eigval[i] > tol) {
            double inv_sqrt = 1.0 / std::sqrt(eigval[i]);
            for (int r = 0; r < naux; r++) {
                for (int c = 0; c < naux; c++) {
                    (*metric_)(r,c) += metric_flat[r + i*naux] * inv_sqrt * metric_flat[c + i*naux];
                }
            }
        }
    }
    
    metric_flat.clear();
    eigval.clear();
    
    metric_->set_name("SO Basis Fitting Inverse (DPC)");

    //is_inverted_ = true;
    //algorithm_ = "DPC";
    //DPC_ = true;

    //form_fitting_metric();

    //std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    //std::shared_ptr<BasisSet> basis_ = reference_wavefunction_->basisset();

    //// This is already squeezed
    //auto eri3c = std::make_shared<IntegralFactory>(aux_, zero, basis_, basis_); //(aux_, zero, aux_, aux_);
    //SharedMatrix D_ao = reference_wavefunction_->Da();
    ////double norm2 = 0.0;
    ////for (int h = 0; h < D_ao->nirrep(); ++h) {
    ////    int rows = D_ao->rowdim(h);
    ////    int cols = D_ao->coldim(h);
    ////    double** Dh = D_ao->pointer(h);
    ////
    ////    for (int i = 0; i < rows; ++i)
    ////        for (int j = 0; j < cols; ++j)
    ////            norm2 += Dh[i][j] * Dh[i][j];
    ////}
    ////double norm_D = std::sqrt(norm2);
    ////std::cout << norm_D << std::endl;

    ////outfile->Printf("||D_ao||_F = %.12e\n", norm_D);
    ////SharedMatrix rhs = mult(D_ao, eri3c, true, false, 1.0, 0.0);
    //int naux = aux_->nbf();       // number of auxiliary functions
    //int nbf  = basis_->nbf();     // number of primary basis functions

    //// Create Bp matrix: naux x (nbf*nbf)
    //std::vector<std::vector<double>> Bp(naux, std::vector<double>(nbf * nbf, 0.0));
    //
    //// Build the integral engine
    //auto fact = std::make_shared<IntegralFactory>(aux_, BasisSet::zero_ao_basis_set(), basis_, basis_);
    //std::shared_ptr<TwoBodyAOInt> eri(fact->eri());
    //
    //// Loop over shells to fill Bp
    //for (int P = 0; P < aux_->nshell(); P++) {
    //    int np = aux_->shell(P).nfunction();
    //    int pstart = aux_->shell(P).function_index();
    //
    //    for (int M = 0; M < basis_->nshell(); M++) {
    //        int nm = basis_->shell(M).nfunction();
    //        int mstart = basis_->shell(M).function_index();
    //
    //        for (int N = 0; N < basis_->nshell(); N++) {
    //            int nn = basis_->shell(N).nfunction();
    //            int nstart = basis_->shell(N).function_index();
    //
    //            // Compute shell (P,M,N) integrals
    //            eri->compute_shell(P, 0, M, N);
    //            const double* buffer = eri->buffer();
    //
    //            // Map shell-local integrals to global AO indices
    //            for (int p = 0, index = 0; p < np; p++) {
    //                for (int m = 0; m < nm; m++) {
    //                    for (int n = 0; n < nn; n++, index++) {
    //                        Bp[p + pstart][(m + mstart) * nbf + (n + nstart)] = buffer[index];
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}
    //
    //// Contract Bp with D_ao to get rhs: rhs_P = sum_uv (P|uv) D_uv
    //std::vector<double> rhs(naux, 0.0);
    //
    //for (int P = 0; P < naux; P++) {
    //    for (int u = 0; u < nbf; u++) {
    //        for (int v = 0; v < nbf; v++) {
    //            rhs[P] += Bp[P][u * nbf + v] * (*D_ao)(u, v);
    //        }
    //    }
    //}
    ////double frobenius_rhs = 0.0;
    ////for (int P = 0; P < naux; P++) {
    ////    frobenius_rhs += rhs[P] * rhs[P];
    ////}
    ////frobenius_rhs = std::sqrt(frobenius_rhs);
    ////
    ////std::cout << "Frobenius norm of rhs: " << frobenius_rhs << std::endl;
    ////// eri3c holds the 3-center integrals (P|uv)
    ////for (int P = 0; P < naux; ++P) {
    ////    double sum = 0.0;
    ////    for (int u = 0; u < nbf; ++u) {
    ////        for (int v = 0; v < nbf; ++v) {
    ////            // eri3c->get(P,u,v) returns (P|uv)
    ////            sum += eri3c->get(P, u, v) * (*D_ao)(u,v);
    ////        }
    ////    }
    ////    (*rhs)(P,0) = sum;
    ////}

    //// Compute eigenvalues and eigenvectors of the metric
    ////std::vector<double> eigvals = metric_->eigenvalues();
    ////Matrix eigvecs = metric_->eigenvectors();
    //
    //// Allocate sigma and Picard coefficient containers
    ////std::vector<double> sigma(eigvals.size());
    ////std::vector<double> picard_coeffs(eigvals.size());
    ////double epsilon_opt = 0.0;
    //
    ////// Loop over eigenvalues
    ////for (size_t i = 0; i < eigvals.size(); ++i) {
    ////    sigma[i] = std::sqrt(std::abs(eigvals[i]));
    ////
    ////    // rhs here is (D_ao * (P|uv)), flatten column access if needed
    ////    picard_coeffs[i] = std::abs(dot(eigvecs.column(i), rhs->column(i))) / sigma[i];
    ////}
    ////
    ////// Find the Picard "knee" to select optimal regularization
    ////for (size_t i = 0; i < sigma.size() - 1; ++i) {
    ////    if (picard_coeffs[i+1] > picard_coeffs[i]) {
    ////        epsilon_opt = sigma[i];
    ////        break;
    ////    }
    ////}
    ////
    ////// Convert to tolerance and regularize the metric
    ////double tol = epsilon_opt * epsilon_opt;
    ////metric_->power(-0.5, tol);
    ////metric_->set_name("SO Basis Fitting Inverse (DPC)");
    //// Step 5: Diagonalize the metric using C_DSYEV (symmetric eigenproblem)
    //// Flatten metric_ matrix into column-major array T
    //std::vector<double> T_flat(naux * naux, 0.0);
    //for (int i = 0; i < naux; i++)
    //    for (int j = 0; j < naux; j++)
    //        T_flat[i * naux + j] = (*metric_)(i, j);

    //std::vector<double> eigval(naux, 0.0);
    //int lwork = naux * 3;
    //std::vector<double> work(lwork, 0.0);

    //int stat = C_DSYEV('v', 'u', naux, T_flat.data(), naux, eigval.data(), work.data(), lwork);
    //if (stat != 0) {
    //    throw std::runtime_error("C_DSYEV failed to diagonalize metric");
    //}

    //// Step 6: Compute Picard coefficients and find optimal tolerance
    //std::vector<double> sigma(naux, 0.0);
    //std::vector<double> picard_coeffs(naux, 0.0);
    //double epsilon_opt = 0.0;

    //// Compute L2 norm of eigenvalues
    ////double eigval_norm = 0.0;
    ////for (int i = 0; i < naux; i++) {
    ////    eigval_norm += eigval[i] * eigval[i];
    ////}
    ////eigval_norm = std::sqrt(eigval_norm);

    ////// Compute Frobenius norm of eigenvectors
    ////double eigvec_norm = 0.0;
    ////for (int i = 0; i < naux * naux; i++) {
    ////    eigvec_norm += T_flat[i] * T_flat[i];
    ////}
    ////eigvec_norm = std::sqrt(eigvec_norm);
    ////std::cout << "Eigenvalue L2 norm: " << eigval_norm << "\n";
    ////std::cout << "Eigenvector Frobenius norm: " << eigvec_norm << "\n";

    //for (int i = 0; i < naux; i++) {
    //    sigma[i] = std::sqrt(std::abs(eigval[i]));

    //    // Dot product of i-th eigenvector with rhs
    //    double dotprod = 0.0;
    //    for (int P = 0; P < naux; P++)
    //        //dotprod += T_flat[P * naux + i] * rhs[P];
    //        dotprod += T_flat[P + i*naux] * rhs[P];
    //    picard_coeffs[i] = std::abs(dotprod) / sigma[i];
    //}

    ////double pc_norm = 0.0;
    ////for (int i = 0; i < naux; i++) {
    ////    pc_norm += picard_coeffs[i] * picard_coeffs[i];
    ////}
    ////pc_norm = std::sqrt(pc_norm);
    ////std::cout << "picard coeffs" << pc_norm << "\n";
    ////std::cout << "Picard coefficients: ";
    ////for (size_t i = 0; i < picard_coeffs.size(); ++i) {
    ////    std::cout << picard_coeffs[i] << " ";
    ////}
    //std::cout << "\n";
    //// Find the "knee" of Picard coefficients
    //for (int i = 0; i < naux - 1; i++) {
    //    if (picard_coeffs[i + 1] > picard_coeffs[i]) {
    //        epsilon_opt = sigma[i];
    //        break;
    //    }
    //}

    //double tol = epsilon_opt * epsilon_opt;
    //std::cout << "tol " << tol << std::endl;
    ////metric_->power(-0.5, tol);
    //// Allocate metric_inv_sqrt
    //// 6. Overwrite metric_ in place with metric_inv_sqrt = U * diag(inv_sqrt(eigvals[mask])) * U^T
    //for (int r = 0; r < naux; r++)
    //    for (int c = 0; c < naux; c++)
    //        (*metric_)(r, c) = 0.0;  // reset
    //
    //for (int i = 0; i < naux; i++) {
    //    if (eigval[i] > tol) {  // mask
    //        double inv_sqrt = 1.0 / std::sqrt(eigval[i]);
    //        for (int r = 0; r < naux; r++) {
    //            for (int c = 0; c < naux; c++) {
    //                (*metric_)(r, c) += T_flat[r + i * naux] * inv_sqrt * T_flat[c + i * naux];
    //            }
    //        }
    //    }
    //}
    //double norm = 0.0;
    //for (int i = 0; i < naux; i++)
    //    for (int j = 0; j < naux; j++)
    //        norm += (*metric_)(i,j) * (*metric_)(i,j);
    //norm = std::sqrt(norm);
    //std::cout << "norm " << norm << std::endl;
    //metric_->set_name("SO Basis Fitting Inverse (DPC)");
}
void FittingMetric::form_eig_inverse_TR() {
    is_inverted_ = true;
    algorithm_ = "TR"; 
    TR_ = true;

    form_fitting_metric();

    auto zero = BasisSet::zero_ao_basis_set();
    auto basis = reference_wavefunction_->basisset();
    auto eri_fact = std::make_shared<IntegralFactory>(aux_, zero, basis, basis);
    auto eri = std::shared_ptr<TwoBodyAOInt>(eri_fact->eri());

    SharedMatrix D_ao = reference_wavefunction_->Da();
    if (!D_ao) {
        throw std::runtime_error("D_ao is null");
    }

    int naux = aux_->nbf();
    int nbf  = basis->nbf();

    if (D_ao->nrow() != nbf || D_ao->ncol() != nbf) {
        throw std::runtime_error("D_ao dimension mismatch");
    }

    std::vector<double> rhs(naux, 0.0);

    // zero basis info
    int n0 = zero->shell(0).nfunction();
    if (n0 <= 0) {
        throw std::runtime_error("Zero basis has invalid nfunction");
    }

    for (int P = 0; P < aux_->nshell(); P++) {

        const auto& Pshell = aux_->shell(P);
        int np = Pshell.nfunction();
        int pstart = Pshell.function_index();

        if (pstart < 0 || pstart + np > naux) {
            throw std::runtime_error("Aux shell index out of bounds");
        }

        for (int M = 0; M < basis->nshell(); M++) {

            const auto& Mshell = basis->shell(M);
            int nm = Mshell.nfunction();
            int mstart = Mshell.function_index();

            if (mstart < 0 || mstart + nm > nbf) {
                throw std::runtime_error("Basis shell M index out of bounds");
	    }

            for (int N = 0; N < basis->nshell(); N++) {

                const auto& Nshell = basis->shell(N);
                int nn = Nshell.nfunction();
                int nstart = Nshell.function_index();

                if (nstart < 0 || nstart + nn > nbf) {
                    throw std::runtime_error("Basis shell N index out of bounds");
                }

                // Compute (P | 0 M N)
                eri->compute_shell(P, 0, M, N);
                const double* buffer = eri->buffer();

                if (!buffer) {
                    throw std::runtime_error("ERI buffer is null");
                }

                int expected_size = np * n0 * nm * nn;
                int index = 0;

                for (int p = 0; p < np; p++) {
                    for (int q = 0; q < n0; q++) {
                        for (int m = 0; m < nm; m++) {
                            for (int n = 0; n < nn; n++) {

                                int Pidx = p + pstart;
                                int midx = m + mstart;
                                int nidx = n + nstart;

                                rhs[Pidx] += buffer[index]
                                             * (*D_ao)(midx, nidx);
                                index++;
                            }
                        }
                    }
                }

                if (index != expected_size) {
                    throw std::runtime_error("ERI buffer index mismatch");
                }
            }
        }
    }

    // --- Flatten metric into contiguous vector for diagonalization ---
    std::vector<double> metric_flat(naux * naux, 0.0);
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            metric_flat[i * naux + j] = (*metric_)(i,j);

    // --- Diagonalize metric (C_DSYEV requires contiguous memory) ---
    std::vector<double> eigval(naux, 0.0);
    int lwork = naux * 3;
    std::vector<double> work(lwork, 0.0);

    int stat = C_DSYEV('v', 'u', naux, metric_flat.data(), naux, eigval.data(), work.data(), lwork);
    if (stat != 0)
        throw std::runtime_error("C_DSYEV failed to diagonalize metric");

    work.clear();  // free workspace

    // --- Compute Picard coefficients ---
    std::vector<double> sigma(naux, 0.0);
    std::vector<double> picard(naux, 0.0);
    for (int i = 0; i < naux; i++) {
        sigma[i] = std::sqrt(std::abs(eigval[i]));

        double dotprod = 0.0;
        for (int P = 0; P < naux; P++)
            dotprod += metric_flat[P + i * naux] * rhs[P];  // use flat eigenvectors

        picard[i] = std::abs(dotprod) / sigma[i];
    }

    rhs.clear();  // free memory

    // --- Find knee of Picard coefficients ---
    double epsilon_opt = 0.0;
    for (int i = 0; i < naux - 1; i++) {
        if (picard[i+1] > picard[i]) {
            epsilon_opt = sigma[i];
            break;
        }
    }
    double tol = epsilon_opt * epsilon_opt;

    picard.clear();

    std::vector<double> d(naux, 0.0);
    for (int i = 0; i < naux; i++) {
        const double s = sigma[i];
        if (s > 0.0) {
            const double s2 = s * s;
            const double f  = s2 / (s2 + tol);
            d[i] = f / s;
        } else {
            d[i] = 0.0;
        }
    }

    sigma.clear();  // free memory

    double d_norm = 0.0;
    for (int i = 0; i < naux; i++) {
        d_norm += d[i] * d[i];
    }
    d_norm = std::sqrt(d_norm);
    std::cout << "\n d_norm" << std::endl; 
    std::cout << d_norm << std::endl; 


    std::cout << "\n metric:" << std::endl;
    metric_->zero();
    // --- Reconstruct inverse metric ---
    for (int i = 0; i < naux; i++) {
	const double di = d[i];
        for (int r = 0; r < naux; r++) {
            const double U_r_i = metric_flat[r + i * naux];
            const double scaled = U_r_i * di;
            for (int c = 0; c < naux; c++) {
                (*metric_)(r,c) += scaled * metric_flat[c + i * naux];
            }
        }
    }

    double norm = 0.0;
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            norm += (*metric_)(i,j) * (*metric_)(i,j);
    norm = std::sqrt(norm);
    std::cout << norm << std::endl;

    //for (int i = 0; i < naux; i++) {
    //    if (eigval[i] > tol) {
    //        double inv_sqrt = 1.0 / std::sqrt(eigval[i]);
    //        for (int r = 0; r < naux; r++) {
    //            for (int c = 0; c < naux; c++) {
    //                (*metric_)(r,c) += metric_flat[r + i*naux] * inv_sqrt * metric_flat[c + i*naux];
    //            }
    //        }
    //    }
    //}
    
    metric_flat.clear();
    eigval.clear();

    metric_->set_name("SO Basis Fitting Inverse (TR)");

}	
void FittingMetric::form_full_eig_inverse(double tol) {
    is_inverted_ = true;
    algorithm_ = "EIG";

    form_fitting_metric();
    metric_->power(-1.0, tol);
    metric_->set_name("SO Basis Fitting Inverse (Eig)");
}
void FittingMetric::form_full_inverse() {
    is_inverted_ = true;
    algorithm_ = "FULL";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {
        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        int info = C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
        // Inverse
        info = C_DPOTRI('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);

        for (int A = 0; A < metric_->colspi()[h]; A++)
            for (int B = 0; B < A; B++) J[A][B] = J[B][A];
    }
    metric_->set_name("SO Basis Fitting Inverse (Full)");
}
void FittingMetric::form_cholesky_factor() {
    is_inverted_ = true;
    algorithm_ = "CHOLESKY";

    form_fitting_metric();

    // pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {
        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        int info = C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
    }
    metric_->set_name("SO Basis Cholesky Factor (Full)");
}
void FittingMetric::pivot() {
    for (int h = 0; h < metric_->nirrep(); h++) {
        if (metric_->colspi()[h] == 0) continue;

        double** J = metric_->pointer(h);
        int* P = pivots_->pointer(h);
        int norbs = metric_->colspi()[h];
        std::vector<double> Temp(norbs);

        // Pivot
        double max;
        int Temp_p;
        int pivot;
        for (int i = 0; i < norbs - 1; i++) {
            max = 0.0;
            // Where's the pivot diagonal?
            for (int j = i; j < norbs; j++)
                if (max <= std::fabs(J[j][j])) {
                    max = std::fabs(J[j][j]);
                    pivot = j;
                }

            // Rows
            C_DCOPY(norbs, &J[pivot][0], 1, Temp.data(), 1);
            C_DCOPY(norbs, &J[i][0], 1, &J[pivot][0], 1);
            C_DCOPY(norbs, Temp.data(), 1, &J[i][0], 1);

            // Columns
            C_DCOPY(norbs, &J[0][pivot], norbs, Temp.data(), 1);
            C_DCOPY(norbs, &J[0][i], norbs, &J[0][pivot], norbs);
            C_DCOPY(norbs, Temp.data(), 1, &J[0][i], norbs);

            Temp_p = P[i];
            P[i] = P[pivot];
            P[pivot] = Temp_p;
        }

        int* R = rev_pivots_->pointer(h);
        for (int i = 0; i < norbs; i++) R[P[i]] = i;
    }
}
}
