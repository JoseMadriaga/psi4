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

#include <fstream>
#include <iomanip>
#include <sstream>

//DPC added
#include <stdexcept>

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

    int naux = aux_->nbf(); 
    // Diagonalize metric BEFORE calling power()
    auto eigvecs = metric_->clone();
    auto eigvals = std::make_shared<Vector>("eigvals", naux);
    
    metric_->diagonalize(eigvecs, eigvals);
    // Count how many eigenvalues survive the tolerance
    int nkept = 0;
    for (int i = 0; i < naux; ++i) {
        if (eigvals->get(i) > tol) {
            ++nkept;
        }
    }
    
    double trunc_ratio = static_cast<double>(nkept) / static_cast<double>(naux);
    
    outfile->Printf("Metric power(-1/2) truncation: kept %d / %d (ratio = %.6f)\n",
                    nkept, naux, trunc_ratio);
    
    metric_->power(-0.5, tol);
    metric_->set_name("SO Basis Fitting Inverse (Eig)");
}
void FittingMetric::form_eig_inverse_DPC() {
    // --- Setup DPC ---
    is_inverted_ = true;
    algorithm_ = "DPC";
    DPC_ = true;
    
    outfile->Printf("DPC-weighted cutoff with TR stability \n");
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

    /// -----------------------------------comment below ---------------
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

    // ----------------------------comment above -----------

    //std::cout << "prior to T_flat" << std::endl;
    //// --- Flatten metric into contiguous vector for diagonalization ---
    //std::vector<double> metric_flat(naux * naux, 0.0);
    //for (int i = 0; i < naux; i++)
    //    for (int j = 0; j < naux; j++)
    //        metric_flat[i * naux + j] = (*metric_)(i,j);
    //
    //std::cout << "prior to diagonalize" << std::endl;
    //// --- Diagonalize metric (C_DSYEV requires contiguous memory) ---
    //std::vector<double> eigval(naux, 0.0);
    //int lwork = naux * 3;
    //std::vector<double> work(lwork, 0.0);
    //
    //int stat = C_DSYEV('v', 'u', naux, metric_flat.data(), naux, eigval.data(), work.data(), lwork);
    //if (stat != 0)
    //    throw std::runtime_error("C_DSYEV failed to diagonalize metric");
    //
    //work.clear();  // free workspace
    //
    //std::cout << "prior to picard_coef" << std::endl;
    //// --- Compute Picard coefficients ---
    //std::vector<double> sigma(naux, 0.0);
    //std::vector<double> picard(naux, 0.0);
    //for (int i = 0; i < naux; i++) {
    //    sigma[i] = std::sqrt(std::abs(eigval[i]));
    //
    //    double dotprod = 0.0;
    //    for (int P = 0; P < naux; P++)
    //        dotprod += metric_flat[P + i * naux] * rhs[P];  // use flat eigenvectors
    //
    //    picard[i] = std::abs(dotprod) / sigma[i];
    //}
    //
    //rhs.clear();  // free memory
    //
    //std::cout << "prior to determining eps_opt" << std::endl;
    //// --- Find knee of Picard coefficients ---
    //double epsilon_opt = 0.0;
    //for (int i = 0; i < naux - 1; i++) {
    //    if (picard[i+1] > picard[i]) {
    //        epsilon_opt = sigma[i];
    //        break;
    //    }
    //}
    //double tol = epsilon_opt * epsilon_opt;
    //std::cout << "tol " << tol << std::endl;
    //
    //picard.clear();
    //sigma.clear();  // free memory
    //
    //std::cout << "prior to inverse metric" << std::endl;
    //// --- Reconstruct inverse metric ---
    //for (int r = 0; r < naux; r++)
    //    for (int c = 0; c < naux; c++)
    //        (*metric_)(r,c) = 0.0;

    //// count kept modes (truncated dimension)
    //int nkept = 0;
    //for (int i = 0; i < naux; ++i) {
    //    if (eigval[i] != 0.0) ++nkept;           // or std::abs(d[i]) > 0.0
    //}

    //const double trunc_ratio = static_cast<double>(nkept) / static_cast<double>(naux);

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

    //outfile->Printf("DPC/TR truncation: kept %d / %d modes (ratio = %.6f)\n",
    //                nkept, naux, trunc_ratio);

    //metric_flat.clear();
    //eigval.clear();
    //
    //metric_->set_name("SO Basis Fitting Inverse (DPC)");


    // ------------------ comment below ------------------
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
    /// ------------------ comment above --------------------  

    /// ------ ad hoc sort of numer and sigma to descending order -----------
    //std::cout << "prior to T_flat" << std::endl;
    // --- Flatten metric into contiguous vector for diagonalization ---
    std::vector<double> metric_flat(naux * naux, 0.0);
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            metric_flat[i * naux + j] = (*metric_)(i, j);

    std::vector<double> A_flat(naux * naux, 0.0);
    A_flat = metric_flat;

    std::cout << "prior to diagonalize" << std::endl;
    // --- Diagonalize metric (C_DSYEV requires contiguous memory) ---
    std::vector<double> eigval(naux, 0.0);
    int lwork = naux * 3;
    std::vector<double> work(lwork, 0.0);

    int stat = C_DSYEV('v', 'u', naux, metric_flat.data(), naux,
                       eigval.data(), work.data(), lwork);
    if (stat != 0)
        throw std::runtime_error("C_DSYEV failed to diagonalize metric");

    work.clear();  // free workspace

    // ------------- adding sort of sigma and uT ---------------------
    // --- Reverse ascending -> descending while preserving correspondence ---
    std::vector<double> eigval_desc(naux, 0.0);
    std::vector<double> metric_desc(naux * naux, 0.0);
    
    for (int new_col = 0; new_col < naux; ++new_col) {
        int old_col = naux - 1 - new_col;
        eigval_desc[new_col] = eigval[old_col];
    
        for (int row = 0; row < naux; ++row) {
            metric_desc[row + new_col * naux] = metric_flat[row + old_col * naux];
        }
    }
    
    //eigval = std::move(eigval_desc);
    //metric_flat = std::move(metric_desc);
    // ------------------------------------------------------------------

    //std::cout << "prior to numer/sigma construction" << std::endl;
    // --- Compute sigma = sqrt(eigval) and numerator = |u_i^T b| ---
    std::vector<double> sigma(naux, 0.0);
    std::vector<double> numer(naux, 0.0);
    std::vector<double> noabs_numer(naux, 0.0);

    for (int i = 0; i < naux; i++) {
        sigma[i] = std::sqrt(std::abs(eigval_desc[i]));

        double dotprod = 0.0;
        for (int P = 0; P < naux; P++)
            dotprod += metric_desc[P + i * naux] * rhs[P];  // eigenvector i

        numer[i] = std::abs(dotprod);
	noabs_numer[i] = dotprod; 
    }

    //std::cout << "prior to descending sort + heuristic picard" << std::endl;

    // ============================================================
    // Heuristic sorting requested by user:
    //   1) sort sigma in descending order
    //   2) sort numerator |u^T b| in descending order
    //   3) form picard_i = numer_desc[i] / sigma_desc[i]
    //
    // WARNING:
    // This breaks the original one-to-one correspondence between
    // sigma_i and its associated numerator |u_i^T b|.
    // Use as a heuristic diagnostic/cutoff, not as the strict DPC.
    // ============================================================

    // ------- redundant code -----------------
    //std::vector<double> sigma_desc = sigma;
    //std::vector<double> numer_desc = numer;
    //std::vector<double> noabs_numer_desc = noabs_numer;

    //std::sort(sigma_desc.begin(), sigma_desc.end(), std::greater<double>());
    //std::sort(numer_desc.begin(), numer_desc.end(), std::greater<double>());

    //std::vector<double> picard_desc(naux, 0.0);
    //for (int i = 0; i < naux; ++i) {
    //    if (sigma[i] > 0.0)
    //        picard_desc[i] = numer[i] / sigma[i];
    //    else
    //        picard_desc[i] = 0.0;
    //}

    // ============================================================
    // Methods of determining the cutoff:
    // (0) DPC inspired criterion:
    //    choose i where the last (picard_i > sigma_i) 
    //   
    // (1) Maximum gap-ratio criterion:
    //    gap_i = picard_i / picard_{i+1}
    //    choose i that maximizes gap_i
    // (2) l-curve 
    //     plot the log of ||\tilde{x}||_{2} vs ||\tilde{r}||_2 
    //     determine i based on elbow method    
    // ============================================================

    int gap_method = 2;
    int cutoff_desc = 0;
    double tol = 0;
    double max_gap_ratio = -1.0;
    double epsilon_sigma = 0;
    if (gap_method == 1) {
        //int cutoff_desc = 0;      // descending-order cutoff index
        //double max_gap_ratio = -1.0;
        //outfile->Printf("\n gap max ratio\n");
        std::vector<double> picard_desc(naux, 0.0);
        for (int i = 0; i < naux; ++i) {
            if (sigma[i] > 0.0)
                picard_desc[i] = numer[i] / sigma[i];
            else
                picard_desc[i] = 0.0;
        }

        for (int i = 0; i < naux - 1; ++i) {
            if (picard_desc[i + 1] > 0.0) {
                double gap = picard_desc[i] / picard_desc[i + 1];
                if (gap > max_gap_ratio) {
                    max_gap_ratio = gap;
                    cutoff_desc = i;
                }
            }
        }
    }
    else if (gap_method == 0) {
	outfile->Printf("\n going into picard > sigma method\n");
        //int cutoff_desc = 0; 
        std::vector<double> picard_desc(naux, 0.0);
        for (int i = 0; i < naux; ++i) {
            if (sigma[i] > 0.0)
                picard_desc[i] = numer[i] / sigma[i];
            else
                picard_desc[i] = 0.0;
        }

	for (int i =  0; i < naux; i++) {
            if (picard_desc[i] > sigma[i]) {
                cutoff_desc = i;
	    }
	}
    }
    else if (gap_method == 2){
	outfile->Printf("\n l-curve method \n");
        // ============================================================
        // MATLAB-style T eigendecomposition L-curve workflow, all inline
        // Assumes:
        //   sigma : descending singular values
        //   U     : eigenvector
        //   (*metric_)    : matrix
        //   rhs     : RHS vector
        // ============================================================
        //
        //int m = static_cast<int>((*metric_).size());
        //outfile->Printf("\n int m %d \n", m);
	//int n = static_cast<int>((*metric_)(0).size());
        
        // ------------------------------------------------------------
        // Full eigendecomposition solution to build reference B1_0 = A * x_full
        // ------------------------------------------------------------
        int rfull = 0;
        for (int i = 0; i < naux; ++i) {
            if (sigma[i] > 0.0) rfull++;
        }

	//outfile->Printf("\n rfull \n");
        //---------uTb -> numer_desc -------------
        //std::vector<double> uty_full(rfull, 0.0);
        //for (int j = 0; j < rfull; ++j) {
        //    for (int i = 0; i < m; ++i) {
        //        uty_full[j] += U[i][j] * y[i];
        //    }
        //}

        //---------- uTb/sigma -> picard_desc-----------
        std::vector<double> zfull(rfull, 0.0);
        for (int i = 0; i < rfull; ++i) {
            zfull[i] = noabs_numer[i] / sigma[i];
	}

        std::vector<double> csvd_full(naux, 0.0);
        for (int i = 0; i < naux; ++i) {
            for (int j = 0; j < rfull; ++j) {
                csvd_full[i] += metric_desc[j + i * naux] * zfull[j];
            }
        }

	//outfile->Printf("\n csvd_full \n");

	//double full_rsum = 0.0;
        //double norm_rhs = 0.0;  
	std::vector<double> B1_0(naux, 0.0);
        for (int i = 0; i < naux; ++i) {
            for (int j = 0; j < naux; ++j) {
                B1_0[i] += (*metric_)(i,j) * csvd_full[j];
            }
	    //full_rsum += B1_0[i] * B1_0[i];
            //norm_rhs += rhs[i] * rhs[i];
        }

	//full_rsum = std::sqrt(full_rsum);
	//norm_rhs = std::sqrt(norm_rhs);
	//outfile->Printf("\n B1_0: %d, rhs: %d \n", full_rsum, norm_rhs);
        // ------------------------------------------------------------
        // Sweep truncations like MATLAB code
        // keep first r singular values, for r = naux-1 down to 1
        // ------------------------------------------------------------
        std::vector<double> soln_norm_matrix;
        std::vector<double> res_norm_matrix;
        std::vector<int>    k_kept_matrix;
        std::vector<double> sigma_cut_matrix;

        soln_norm_matrix.reserve(naux);
        res_norm_matrix.reserve(naux);
        k_kept_matrix.reserve(naux);
        sigma_cut_matrix.reserve(naux);

        for (int r = naux - 1; r >= 1; --r) {

            // zsvdT1 = ShatT1 \ (UhatT1' * y)
            std::vector<double> uty(r, 0.0);
            for (int j = 0; j < r; ++j) {
                for (int i = 0; i < naux; ++i) {
                    uty[j] += metric_desc[i + j *naux] * rhs[i];
                }
            }

            std::vector<double> z(r, 0.0);
            for (int i = 0; i < r; ++i) {
                z[i] = uty[i] / sigma[i];
            }

            // csvdT1 = VhatT1 * zsvdT1
            std::vector<double> csvdT1(naux, 0.0);
            for (int i = 0; i < naux; ++i) {
                for (int j = 0; j < r; ++j) {
                    csvdT1[i] += metric_desc[j + i * naux] * z[j];
                }
            }

            // solution norm
            double xsum = 0.0;
            for (int i = 0; i < naux; ++i) {
                xsum += csvdT1[i] * csvdT1[i];
            }
            double xnorm = std::sqrt(xsum);

            // BprimeSVDT1 = Ak * csvdT1
            std::vector<double> BprimeSVDT1(naux, 0.0);
            for (int i = 0; i < naux; ++i) {
                for (int j = 0; j < naux; ++j) {
                    BprimeSVDT1[i] += (*metric_)(i,j) * csvdT1[j];
                }
            }

            // residual norm: MATLAB-style workflow
            // res = ||B1_0 - BprimeSVDT1||
            double rsum = 0.0;
            for (int i = 0; i < naux; ++i) {
                double diff = B1_0[i] - BprimeSVDT1[i];
                rsum += diff * diff;
            }
            double rnorm = std::sqrt(rsum);

            soln_norm_matrix.push_back(xnorm);
            res_norm_matrix.push_back(rnorm);
            k_kept_matrix.push_back(r);
            sigma_cut_matrix.push_back(sigma[r - 1]);  // smallest kept sigma
        }

        // ------------------------------------------------------------
        // Find L-curve corner by discrete curvature on log-log scale
        // ------------------------------------------------------------
        // using cutoff_desc instead of best_idx
        //int best_idx = 0;
        double best_curv = -1.0;

        if (soln_norm_matrix.size() >= 3) {
            for (int k = 1; k < static_cast<int>(soln_norm_matrix.size()) - 1; ++k) {

                double x1 = std::log(res_norm_matrix[k - 1]  + 1.0e-300);
                double y1 = std::log(soln_norm_matrix[k - 1] + 1.0e-300);
                double x2 = std::log(res_norm_matrix[k]      + 1.0e-300);
                double y2 = std::log(soln_norm_matrix[k]     + 1.0e-300);
                double x3 = std::log(res_norm_matrix[k + 1]  + 1.0e-300);
                double y3 = std::log(soln_norm_matrix[k + 1] + 1.0e-300);

                double ax = x2 - x1;
                double ay = y2 - y1;
                double bx = x3 - x2;
                double by = y3 - y2;
                double cx = x3 - x1;
                double cy = y3 - y1;

                double a = std::sqrt(ax * ax + ay * ay);
                double b = std::sqrt(bx * bx + by * by);
                double c = std::sqrt(cx * cx + cy * cy);

                double area2 = std::abs(ax * cy - ay * cx);
                double curv = 0.0;

                if (a > 0.0 && b > 0.0 && c > 0.0) {
                    curv = 2.0 * area2 / (a * b * c);
                }

                if (curv > best_curv) {
                    best_curv = curv;
                    cutoff_desc = k;
                }
            }
        }

        int k_opt = k_kept_matrix[cutoff_desc];
        // don't need this I think, decided at the different method for hard cutoff or tol 
        epsilon_sigma = sigma_cut_matrix[cutoff_desc];
        //double tol = sigma[cutoff_desc] * sigma[cutoff_desc + 1];

        outfile->Printf("\n=== TSVD L-curve analysis ===\n");
	outfile->Printf("Cutoff desc - %d\n", cutoff_desc);
        outfile->Printf("Selected truncation k_opt = %d\n", k_opt);
        outfile->Printf("Corner sigma_k = %.12e\n", epsilon_sigma);
        //outfile->Printf("Suggested eigval cutoff tol = sigma_k^2 = %.12e\n", tol);
        //outfile->Printf("Discrete curvature at corner = %.12e\n", best_curv);

        //outfile->Printf("\n   idx    k_kept         ||r||               ||x||             sigma_k\n");
        //for (int i = 0; i < static_cast<int>(k_kept_matrix.size()); ++i) {
        //    outfile->Printf("%6d %8d   %16.8e   %16.8e   %16.8e\n",
        //                    i,
        //                    k_kept_matrix[i],
        //                    res_norm_matrix[i],
        //                    soln_norm_matrix[i],
        //                    sigma_cut_matrix[i]);
        //}

        //outfile->Printf("\nNeighborhood around selected corner:\n");
        //outfile->Printf("   idx    k_kept         ||r||               ||x||             sigma_k\n");
        //for (int i = std::max(0, cutoff_desc - 2);
        //     i <= std::min((int)k_kept_matrix.size() - 1, cutoff_desc + 2);
        //     ++i) {
        //    outfile->Printf("%6d %8d   %16.8e   %16.8e   %16.8e\n",
        //                    i,
        //                    k_kept_matrix[i],
        //                    res_norm_matrix[i],
        //                    soln_norm_matrix[i],
        //                    sigma_cut_matrix[i]);
        //}

        // ------------------------------------------------------------
        // Write L-curve data to file
        // ------------------------------------------------------------
        //std::ofstream dat("lcurve_data.txt");
        //dat << std::setprecision(16);
        //dat << "# idx  k_kept  residual_norm  solution_norm  sigma_k\n";
        //for (int i = 0; i < static_cast<int>(k_kept_matrix.size()); ++i) {
        //    dat << i << " "
        //        << k_kept_matrix[i] << " "
        //        << res_norm_matrix[i] << " "
        //        << soln_norm_matrix[i] << " "
        //        << sigma_cut_matrix[i] << "\n";
        //}
        //dat.close();

        //outfile->Printf("\nWrote L-curve data to lcurve_data.txt\n");

        // ------------------------------------------------------------
        // Write Python plotting script
        // ------------------------------------------------------------
       // std::ofstream py("lcurve_plot.py");
       // py << "import numpy as np\n";
       // py << "import matplotlib.pyplot as plt\n\n";
       // py << "data = np.loadtxt('lcurve_data.txt', comments='#')\n";
       // py << "idx = data[:,0].astype(int)\n";
       // py << "k_kept = data[:,1].astype(int)\n";
       // py << "rnorm = data[:,2]\n";
       // py << "xnorm = data[:,3]\n";
       // py << "sigma_k = data[:,4]\n\n";

       // py << "cutoff_desc = 0\n";
       // py << "best_curv = -1.0\n";
       // py << "for k in range(1, len(rnorm)-1):\n";
       // py << "    x1 = np.log(rnorm[k-1] + 1.0e-300)\n";
       // py << "    y1 = np.log(xnorm[k-1] + 1.0e-300)\n";
       // py << "    x2 = np.log(rnorm[k]   + 1.0e-300)\n";
       // py << "    y2 = np.log(xnorm[k]   + 1.0e-300)\n";
       // py << "    x3 = np.log(rnorm[k+1] + 1.0e-300)\n";
       // py << "    y3 = np.log(xnorm[k+1] + 1.0e-300)\n";
       // py << "    ax, ay = x2-x1, y2-y1\n";
       // py << "    bx, by = x3-x2, y3-y2\n";
       // py << "    cx, cy = x3-x1, y3-y1\n";
       // py << "    a = np.hypot(ax, ay)\n";
       // py << "    b = np.hypot(bx, by)\n";
       // py << "    c = np.hypot(cx, cy)\n";
       // py << "    area2 = abs(ax*cy - ay*cx)\n";
       // py << "    curv = 0.0\n";
       // py << "    if a > 0 and b > 0 and c > 0:\n";
       // py << "        curv = 2.0 * area2 / (a*b*c)\n";
       // py << "    if curv > best_curv:\n";
       // py << "        best_curv = curv\n";
       // py << "        cutoff_desc = k\n\n";

       // py << "plt.figure(figsize=(7,5))\n";
       // py << "plt.loglog(rnorm, xnorm, marker='o')\n";
       // py << "plt.loglog(rnorm[cutoff_desc], xnorm[cutoff_desc], marker='s', markersize=9)\n";
       // py << "plt.xlabel(r'$||r_k||$')\n";
       // py << "plt.ylabel(r'$||x_k||$')\n";
       // py << "plt.title('TSVD L-curve')\n";
       // py << "plt.annotate(f'k={k_kept[cutoff_desc]}\\nσ={sigma_k[cutoff_desc]:.2e}',\n";
       // py << "             (rnorm[cutoff_desc], xnorm[cutoff_desc]),\n";
       // py << "             textcoords='offset points', xytext=(10,10))\n";
       // py << "plt.grid(True, which='both', ls='--', alpha=0.5)\n";
       // py << "plt.tight_layout()\n";
       // py << "plt.savefig('lcurve_plot.png', dpi=300)\n";
       // py << "plt.show()\n";
       // py.close();

       // outfile->Printf("Wrote plotting script to lcurve_plot.py\n");
       // outfile->Printf("Run: python lcurve_plot.py\n");
    }

    //double epsilon_sigma = sigma_desc[cutoff_desc];
    //double tol = sigma[cutoff_desc] * sigma[cutoff_desc + 1] ;

    // Choose method:
    //   1 = hard cutoff beyond sigma_desc[cutoff_desc]
    //   2 = threshold eigenvalues by sigma_desc[cutoff_desc]^2
    int cutoff_method = 2;

    outfile->Printf("\nHeuristic sorted-numerator Picard analysis:\n");
    outfile->Printf("  cutoff_method             = %d\n", cutoff_method);
    //outfile->Printf("  cutoff_desc index         = %d\n", cutoff_desc);
    //outfile->Printf("  picard[i]                 = %.12e\n", picard_desc[cutoff_desc]);
    //outfile->Printf("  picard[i+1]               = %.12e\n", picard_desc[cutoff_desc + 1]);
    //outfile->Printf("  max gap ratio             = %.12e\n", max_gap_ratio);
    //outfile->Printf("  sigma_cut                 = %.12e\n", sigma_desc[cutoff_desc]);
   // outfile->Printf("  sigma_cut + 1             = %.12e\n", sigma[cutoff_desc + 1]);
    //outfile->Printf("  sigma_cut^2 (= tol)       = %.12e\n", tol);

    //outfile->Printf("\nTop few heuristic sorted entries around cutoff:\n");
    //outfile->Printf("  idx         numer_desc            sigma_desc            picard_desc\n");
    //for (int i = std::max(0, cutoff_desc - 3); i <= std::min(naux - 1, cutoff_desc + 3); ++i) {
    //    outfile->Printf("%5d   %18.10e   %18.10e   %18.10e\n",
    //                    i, numer_desc[i], sigma_desc[i], picard_desc[i]);
    //}

    // ============================================================
    // Write full sigma / numer / picard table for plotting
    // ============================================================
    //std::ofstream picard_out("picard_data_sorted.txt");
    //if (picard_out) {
    //    picard_out << "# idx sigma_desc numer_desc picard_desc\n";
    //    picard_out << std::setprecision(16);
    //    for (int i = 0; i < naux; ++i) {
    //        picard_out << i << " "
    //                   << sigma_desc[i] << " "
    //                   << numer_desc[i] << " "
    //                   << picard_desc[i] << "\n";
    //    }
    //    picard_out.close();
    //} else {
    //    outfile->Printf("Warning: could not open picard_data_sorted.txt for writing\n");
    //}
    //std::cout << "prior to inverse metric" << std::endl;
    // --- Reconstruct inverse metric ---
    for (int r = 0; r < naux; r++)
        for (int c = 0; c < naux; c++)
            (*metric_)(r, c) = 0.0;

    int nkept = 0;

    if (cutoff_method == 1) {
        // --------------------------------------------------------
        // METHOD 1:
        // Hard cutoff by descending sigma index.
        //
        // Since eigval is ascending, keeping the largest (cutoff_desc+1)
        // sigma values means keeping the largest (cutoff_desc+1)
        // eigenvalues. In ascending eigval indexing, that corresponds to:
        //
        //     i >= naux - 1 - cutoff_desc
        // --------------------------------------------------------
        int keep_from_asc = naux - 1 - cutoff_desc;

        for (int i = 0; i < naux; ++i) {
            if (i >= keep_from_asc) ++nkept;
        }

        const double trunc_ratio =
            static_cast<double>(nkept) / static_cast<double>(naux);
      
        for (int i = 0; i < naux; i++) {
            if (i >= keep_from_asc) {
                if (eigval[i] > 0.0) {
                    double inv_sqrt = 1.0 / std::sqrt(eigval[i]);
                    for (int r = 0; r < naux; r++) {
                        for (int c = 0; c < naux; c++) {
                            (*metric_)(r, c) += metric_flat[r + i * naux] *
                                                inv_sqrt *
                                                metric_flat[c + i * naux];
                        }
                    }
                }
            }
        }

        outfile->Printf("\nDPC Method 1 (hard cutoff beyond sigma_desc[%d]): kept %d / %d modes (ratio = %.6f)\n",
                        cutoff_desc, nkept, naux, trunc_ratio);

    } else if (cutoff_method == 2) {
        // --------------------------------------------------------
        // METHOD 2:
        // Use tol = sigma_desc[cutoff_desc]^2, then keep eigval >= tol
        // --------------------------------------------------------
        outfile->Printf("n epsilon_sigma = %.12e\n", epsilon_sigma);
	double tol = epsilon_sigma * epsilon_sigma;
        outfile->Printf("\n tol = %.12e\n", tol);	
        for (int i = 0; i < naux; ++i) {
            if (eigval[i] >= tol) ++nkept;
	    outfile->Printf("\n count = %d, %.12e\n", i, eigval[i]);  
        }

        const double trunc_ratio =
            static_cast<double>(nkept) / static_cast<double>(naux);

        for (int i = 0; i < naux; i++) {
            if (eigval[i] >= tol) {
                double inv_sqrt = 1.0 / std::sqrt(eigval[i]);
                for (int r = 0; r < naux; r++) {
                    for (int c = 0; c < naux; c++) {
                        (*metric_)(r, c) += metric_flat[r + i * naux] *
                                            inv_sqrt *
                                            metric_flat[c + i * naux];
                    }
                }
            }
        }

        outfile->Printf("\nDPC Method 2 (threshold by sigma_desc[%d]^2): kept %d / %d modes (ratio = %.6f)\n",
                        cutoff_desc, nkept, naux, trunc_ratio);
    } else {
        throw std::runtime_error("Invalid cutoff_method: choose 1 or 2");
    }

    //picard_desc.clear();
    //numer_desc.clear();
    //sigma_desc.clear();
    numer.clear();
    sigma.clear();

    metric_flat.clear();
    eigval.clear();

    metric_->set_name("SO Basis Fitting Inverse (DPC)");
    }
  
    /// ---------------- l curve method --------------------------
    //std::cout << "prior to T_flat" << std::endl;
    //// --- Flatten metric into contiguous vector for diagonalization ---
    //std::vector<double> metric_flat(naux * naux, 0.0);
    //for (int i = 0; i < naux; i++)
    //    for (int j = 0; j < naux; j++)
    //        metric_flat[i * naux + j] = (*metric_)(i,j);
    //
    //std::cout << "prior to diagonalize" << std::endl;
    //// --- Diagonalize metric (C_DSYEV requires contiguous memory) ---
    //std::vector<double> eigval(naux, 0.0);
    //int lwork = naux * 3;
    //std::vector<double> work(lwork, 0.0);
    //
    //int stat = C_DSYEV('v', 'u', naux, metric_flat.data(), naux, eigval.data(), work.data(), lwork);
    //if (stat != 0)
    //    throw std::runtime_error("C_DSYEV failed to diagonalize metric");
    //
    //work.clear();  // free workspace
    //
    //std::cout << "prior to picard_coef" << std::endl;
    //// --- Compute singular values and |u_i^T b| ---
    //std::vector<double> sigma(naux, 0.0);
    //std::vector<double> numer(naux, 0.0);
    //std::vector<double> picard(naux, 0.0);
    //
    //for (int i = 0; i < naux; i++) {
    //    sigma[i] = std::sqrt(std::abs(eigval[i]));
    //
    //    double dotprod = 0.0;
    //    for (int P = 0; P < naux; P++)
    //        dotprod += metric_flat[P + i * naux] * rhs[P];  // eigenvector i
    //
    //    numer[i]  = std::abs(dotprod);
    //    picard[i] = (sigma[i] > 0.0 ? numer[i] / sigma[i] : 0.0);
    //}
    //
    //rhs.clear();  // free memory
    //
    //std::cout << "prior to determining eps_opt" << std::endl;
    //
    //// ------------------------------------------------------------------
    //// --- L-curve construction using descending singular values order ---
    //// C_DSYEV gives ascending eigvals, so largest sigma is at naux-1.
    //// For each k = 1..naux:
    ////   ||x_k||^2 = sum_{kept modes} (numer/sigma)^2
    ////   ||r_k||^2 = sum_{discarded modes} numer^2
    //// ------------------------------------------------------------------
    //std::vector<double> xnorm(naux, 0.0);
    //std::vector<double> rnorm(naux, 0.0);
    //
    //double xsum = 0.0;
    //
    //// residual for k=1 starts as sum over all but largest mode, etc.
    //for (int k = 0; k < naux; ++k) {
    //    int idx = naux - 1 - k;  // descending-order index
    //
    //    if (sigma[idx] > 0.0)
    //        xsum += (numer[idx] / sigma[idx]) * (numer[idx] / sigma[idx]);
    //    xnorm[k] = std::sqrt(xsum);
    //
    //    double rsum = 0.0;
    //    for (int j = 0; j < idx; ++j)  // smaller singular values not yet kept
    //        rsum += numer[j] * numer[j];
    //    rnorm[k] = std::sqrt(rsum);
    //}
    //
    //// --- Find corner of log-scale L-curve with simple discrete curvature ---
    //int k_opt = 1;   // 1-based truncation index
    //double best_curv = -1.0;
    //
    //for (int k = 1; k < naux - 1; ++k) {
    //    // points: (log ||r_k||, log ||x_k||)
    //    double x1 = std::log(rnorm[k - 1] + 1.0e-300);
    //    double y1 = std::log(xnorm[k - 1] + 1.0e-300);
    //    double x2 = std::log(rnorm[k]     + 1.0e-300);
    //    double y2 = std::log(xnorm[k]     + 1.0e-300);
    //    double x3 = std::log(rnorm[k + 1] + 1.0e-300);
    //    double y3 = std::log(xnorm[k + 1] + 1.0e-300);
    //
    //    double ax = x2 - x1, ay = y2 - y1;
    //    double bx = x3 - x2, by = y3 - y2;
    //    double cx = x3 - x1, cy = y3 - y1;
    //
    //    double a = std::sqrt(ax * ax + ay * ay);
    //    double b = std::sqrt(bx * bx + by * by);
    //    double c = std::sqrt(cx * cx + cy * cy);
    //
    //    double area2 = std::abs(ax * cy - ay * cx);  // 2 * triangle area
    //    double curv = 0.0;
    //    if (a > 0.0 && b > 0.0 && c > 0.0)
    //        curv = 2.0 * area2 / (a * b * c);
    //
    //    if (curv > best_curv) {
    //        best_curv = curv;
    //        k_opt = k + 1;  // convert to 1-based k
    //    }
    //}
    //
    //// --- map k_opt back to sigma_k in descending order ---
    //int idx_opt = naux - k_opt;
    //double epsilon_sigma = sigma[idx_opt];
    //double tol = epsilon_sigma * epsilon_sigma;
    //
    //outfile->Printf("L-curve selected k = %d\n", k_opt);
    //outfile->Printf("Corresponding sigma_k = %.12e\n", epsilon_sigma);
    //outfile->Printf("Using eigval cutoff tol = sigma_k^2 = %.12e\n", tol);
    //
    //// Optional: print a few nearby L-curve points
    //outfile->Printf("\n   k           ||r_k||            ||x_k||            sigma_k\n");
    //for (int kk = std::max(1, k_opt - 2); kk <= std::min(naux, k_opt + 2); ++kk) {
    //    int id = naux - kk;
    //    outfile->Printf("%4d   %16.8e   %16.8e   %16.8e\n",
    //                    kk, rnorm[kk - 1], xnorm[kk - 1], sigma[id]);
    //}
    //
    //picard.clear();
    //xnorm.clear();
    //rnorm.clear();
    //numer.clear();
    //sigma.clear();
    //
    //std::cout << "prior to inverse metric" << std::endl;
    //// --- Reconstruct inverse metric ---
    //for (int r = 0; r < naux; r++)
    //    for (int c = 0; c < naux; c++)
    //        (*metric_)(r,c) = 0.0;
    //
    //// count kept modes using tol
    //int nkept = 0;
    //for (int i = 0; i < naux; ++i) {
    //    if (eigval[i] > tol) ++nkept;
    //}
    //
    //const double trunc_ratio = static_cast<double>(nkept) / static_cast<double>(naux);
    //
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
    //
    //outfile->Printf("DPC/L-curve truncation: kept %d / %d modes (ratio = %.6f)\n",
    //                nkept, naux, trunc_ratio);
    //
    //metric_flat.clear();
    //eigval.clear();
    //
    //metric_->set_name("SO Basis Fitting Inverse (DPC)");
    //}
    ///// ------------- l curve method above ----------------
    
    void FittingMetric::form_eig_inverse_TR() {
        is_inverted_ = true;
        algorithm_ = "TR"; 
        TR_ = true;
    
    
        outfile->Printf("DPC-weighted filter for TR \n");
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
    
        // count kept modes (truncated dimension)
        int nkept = 0;
        for (int i = 0; i < naux; ++i) {
            if (d[i] != 0.0) ++nkept;           // or std::abs(d[i]) > 0.0
        }
    
        const double trunc_ratio = static_cast<double>(nkept) / static_cast<double>(naux);
    
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
    
        outfile->Printf("DPC/TR truncation: kept %d / %d modes (ratio = %.6f)\n",
                    nkept, naux, trunc_ratio);
    
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
