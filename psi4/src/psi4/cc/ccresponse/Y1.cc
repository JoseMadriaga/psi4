/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"


namespace psi {
namespace ccresponse {

void local_filter_T1(dpdfile2 *);
void lambda_residuals();
void L2HX2(const char *pert, int irrep, double omega);

void Y1_inhomogenous_build(const char *pert, int irrep, double omega) {
    dpdfile2 newYIA, newYia, YIA, Yia;
    dpdfile2 dIA, dia, Fme, FME;
    dpdfile2 YFaet2, YFAEt2, YFmit2, YFMIt2;
    dpdfile2 GMI, Gmi, Gae, XIA, Xia;
    dpdfile2 GAE;
    dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ;
    dpdbuf4 WMBIJ, Wmbij, WMbIj, WmBiJ;
    dpdbuf4 YIJAB, Yijab, YIjAb, YiJaB, Y2;
    dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE;
    dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF; //W;
    dpdbuf4 Z;
    dpdfile2 YLD;
    int Gim, Gi, Gm, Ga, Gam, nrows, ncols, A, a, am;
    int Gei, ei, e, i, Gef, Ge, Gf, E, I, af, fa, f;
    double *X;

    //Added
    dpdfile2 F, z1, z2;
    dpdbuf4 W, WL ,D, X2, Z2, Z3, lx_iajb, X2test, L2test, LIjAb; 
    dpdfile2 Y1, Y1new, mu1, L1, lt, lx, lx_AB, X1;
    dpdbuf4 L2, mu2, Hx_ijab, lx_ijab;
    char lbl[32];
    //int L_irr;
    //L_irr = irrep;
    double Y1_norm;
    double *Y;
    dpdfile2 test;

    outfile->Printf("\tY1 Build...");

    //sprintf(lbl, "New Y_%s_IA (%5.3f)", pert, omega);
    sprintf(lbl, "Inhomo Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);

/*
    // Homogenous terms
    sprintf(lbl, "New Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl);

    global_dpd_->file2_axpy(&Y1, &Y1new, omega, 0);

    // L1 RHS += Yie*Fea
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&Y1, &F, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F);

    // L1 RHS += -Yma*Fim 
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&F, &Y1, &Y1new, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F);

    // L1 RHS += Yme*Wieam  
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->contract422(&W, &Y1, &Y1new, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&Y1);

    // L1 RHS += 1/2 Yimef*Wefam 
    // Y(i,a) += [ 2 Y(im,ef) - Y(im,fe) ] * W(am,ef) 
    // Note: W(am,ef) is really Wabei (ei,ab) //
    sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
    //global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "2 YIjAb - YIjBa");
        //    dpd_contract442(&Y2, &W, &newYIA, 0, 0, 1.0, 1.0); 
        global_dpd_->file2_mat_init(&Y1new);
        global_dpd_->file2_mat_rd(&Y1new);
        for (Gam = 0; Gam < moinfo.nirreps; Gam++) {
            Gef = Gam; // W is totally symmetric 
            Gim = Gef ^ irrep;
            global_dpd_->buf4_mat_irrep_init(&Y2, Gim);
            global_dpd_->buf4_mat_irrep_rd(&Y2, Gim);
            global_dpd_->buf4_mat_irrep_shift13(&Y2, Gim);

            for (Gi = 0; Gi < moinfo.nirreps; Gi++) {
                Ga = Gi ^ irrep;
                Gm = Ga ^ Gam;
                W.matrix[Gam] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gm], W.params->coltot[Gam]);

                nrows = moinfo.occpi[Gi];
                ncols = moinfo.occpi[Gm] * W.params->coltot[Gam];

                for (A = 0; A < moinfo.virtpi[Ga]; A++) {
                    a = moinfo.vir_off[Ga] + A;
                    am = W.row_offset[Gam][a];

                    global_dpd_->buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);

                    if (nrows && ncols && moinfo.virtpi[Ga])
                        C_DGEMV('n', nrows, ncols, 1, Y2.shift.matrix[Gim][Gi][0], ncols, W.matrix[Gam][0], 1, 1,
                                &(Y1new.matrix[Gi][0][A]), moinfo.virtpi[Ga]);
                }

                global_dpd_->free_dpd_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gam]);
            }
            global_dpd_->buf4_mat_irrep_close(&Y2, Gim);
        }
        global_dpd_->file2_mat_wrt(&Y1new);
        global_dpd_->file2_mat_close(&Y1new);
        global_dpd_->buf4_close(&Y2);
        global_dpd_->buf4_close(&W);

    // Y1 RHS += -1/2 Ymnae*Wiemn 
    sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract442(&WMbIj, &Y2, &Y1new, 0, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&WMbIj);

    //L1 RHS += -Gef*Weifa //
    //sprintf(lbl, "G_%s_AE (%5.3f)", pert, omega);
    //global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); 
    //global_dpd_->dot13(&GAE,&W,&Y1new, 0, 0, -1.0, 1.0);
    //global_dpd_->buf4_close(&W); 
     //global_dpd_->file2_close(&GAE); 
    // Above code replaced to remove disk-space and memory bottlenecks 7/26/05, -TDC /
    sprintf(lbl, "G_%s_AE (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->file2_mat_init(&GAE);
    global_dpd_->file2_mat_rd(&GAE);
    global_dpd_->file2_mat_init(&Y1new);
    global_dpd_->file2_mat_rd(&Y1new);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");

    for (Gei = 0; Gei < moinfo.nirreps; Gei++) {
        global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
        X = init_array(W.params->coltot[Gei]);
        for (ei = 0; ei < W.params->rowtot[Gei]; ei++) {
            global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
            e = W.params->roworb[Gei][ei][0];
            i = W.params->roworb[Gei][ei][1];
            Ge = W.params->psym[e];
            Gf = Ge ^ irrep;
            Gi = Ge ^ Gei;
            Ga = Gi ^ irrep;
            E = e - moinfo.vir_off[Ge];
            I = i - moinfo.occ_off[Gi];

            zero_arr(X, W.params->coltot[Gei]);

            for (fa = 0; fa < W.params->coltot[Gei]; fa++) {
                f = W.params->colorb[Gei][fa][0];
                a = W.params->colorb[Gei][fa][1];
                af = W.params->colidx[a][f];
                X[fa] = 2.0 * W.matrix[Gei][0][fa] - W.matrix[Gei][0][af];
            }

            nrows = moinfo.virtpi[Gf];
            ncols = moinfo.virtpi[Ga];
            if (nrows && ncols)
                  C_DGEMV('t', nrows, ncols, -1, &X[W.col_offset[Gei][Gf]], ncols, GAE.matrix[Ge][E], 1, 1,
                         Y1new.matrix[Gi][I], 1);
        }
        global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
        free(X);

    }
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_mat_wrt(&Y1new);
    global_dpd_->file2_mat_close(&Y1new);
    global_dpd_->file2_mat_close(&GAE);
    global_dpd_->file2_close(&GAE);

    // Y1 RHS += -Gmn*Wmina 
    sprintf(lbl, "G_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->dot14(&GMI, &WmNiE, &Y1new, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&WmNiE);
    global_dpd_->file2_close(&GMI);

*/

    // Inhomogenous terms

    /*** Mu * L1 ***/ 

    sprintf(lbl, "%s_IA", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_axpy(&mu1, &Y1new, 2, 0);  
    global_dpd_->file2_close(&mu1);

    /*** L1 * MuBAR + L2 * MuBAR ***/

    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    sprintf(lbl, "%sBAR_MI", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->contract222(&mu1, &L1, &Y1new, 0, 1, -1, 1.0);
    global_dpd_->file2_close(&mu1);

    sprintf(lbl, "%sBAR_AE", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->contract222(&L1, &mu1, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&mu1);
    global_dpd_->file2_close(&L1);

    sprintf(lbl, "%s_IA", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_init(&lt, PSIF_CC_OEI, 0, 0, 0, "Lt_IJ");
    global_dpd_->contract222(&lt, &mu1, &Y1new, 0, 1, 1, 1.0);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    sprintf(lbl, "%sBAR_MbIj", pert, omega);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    global_dpd_->contract442(&mu2, &L2, &Y1new, 0, 2, -0.5, 1.0);

    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&mu2);


    //Here Sort L
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, pqsr, 0, 5, "(2 LIjAb - LIjBa) (ij|ba)");
    global_dpd_->buf4_close(&L2);    

    sprintf(lbl, "%sBAR_MbIj",pert, omega);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    sprintf(lbl, "%sBAR_MbjI", pert, omega);    
    global_dpd_->buf4_sort(&mu2, PSIF_CC_LR, pqsr, 10, 0, lbl);
    global_dpd_->buf4_close(&mu2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "(2 LIjAb - LIjBa) (ij|ba)");
    sprintf(lbl, "%sBAR_MbjI", pert);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    global_dpd_->contract442(&mu2, &L2, &Y1new, 0, 2, -0.5, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&mu2);
    //END HERE

     //-----------------------

    // <O|[Hbar(0), X1]|0>
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract422(&D, &X1, &Y1new, 0, 0, 2, 1); 
    global_dpd_->buf4_close(&D); 
    global_dpd_->file2_close(&X1);   


    // <O|L1(0)|[Hbar(0), X1]|0>
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->contract424(&W, &L1, &Z, 3, 0, 0, -1, 1);

    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 0, 11, "2WMnIe - WnMIe (nM,eI)");   //sort
    global_dpd_->buf4_close(&W);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    //global_dpd_->dot14(&X1, &Z, &Y1new, 0, 0, -1, 1);
 
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (nM,eI)");
    global_dpd_->contract424(&W, &L1, &Z2, 3, 0, 0, -1, 1);   
    global_dpd_->buf4_close(&W);
    //global_dpd_->file2_close(&L1);
    //global_dpd_->dot13(&X1, &Z2, &test, 0, 0, -1, 0);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, pqsr, 0, 5, "Z2 (ij,ba)");   //sort
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ba)");
    global_dpd_->buf4_axpy(&Z2, &Z, 1);  
    //global_dpd_->dot14(&X1, &Z, &Y1new, 0, 0, -1, 1);
    //global_dpd_->dot14(&X1, &Z2, &test, 0, 0, -1, 0);
    //global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z2);
    //global_dpd_->buf4_close(&Z); 

//Here Compute out of core!!!!!!!!

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_scmcopy(&W, PSIF_CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 10, 5, "WAmEf 2(mA,Ef) - (mA,fE)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z3 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf 2(mA,Ef) - (mA,fE)");
    //global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    //global_dpd_->contract424(&W, &L1, &Z2, 1, 1, 1, 1, 0); 
    global_dpd_->contract424(&W, &L1, &Z2, 1, 1, 1, 1, 0);        //What is the difference between contract424 and contract244 

    global_dpd_->buf4_axpy(&Z2, &Z, 1);

    //global_dpd_->dot14(&X1, &Z, &Y1new, 0, 0, 1, 1);	   


    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W); 

  
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z4 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)");
    global_dpd_->contract244(&L1, &W, &Z2, 1, 0, 0, 1, 1);
   
    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->dot14(&X1, &Z, &Y1new, 0, 0, 1, 1);

    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&L1); 
    global_dpd_->buf4_close(&Z);
 
 
    // <O|L2(0)|[Hbar(0), X1]|0>
    //************This part of the code take straight from L2.cc******************
    /* Generate spin-adapted Type-I Lambda-residual contributions to LHX1Y1 */

    //******** COmbine these two!! ***************************
    //Make the intermediate "WMnIj + WnMIj"
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Znew(ij,ab)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->contract444(&W, &L2, &Z, 0, 1, 0.5, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&W);

    //sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    //global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    //global_dpd_->dot13(&X1, &Z, &Y1new, 0, 0, 1, 1);

    //global_dpd_->file2_close(&X1); 
    //global_dpd_->buf4_close(&Z);

    //Sort WMnIj -> WnMIj
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 0, 0, "WMnIj (nM,Ij)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Znew(ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj (nM,Ij)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract444(&W, &L2, &Z2, 0, 1, 0.5, 0);
    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    //sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    //global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    //global_dpd_->dot23(&X1, &Z, &Y1new, 0, 0, 1, 1);
    //global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z2);
    //global_dpd_->buf4_close(&Z);

    //******** COmbine these two!! ***************************
    //Hvvvv x L2    
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2");

    global_dpd_->buf4_scmcopy(&WL, PSIF_CC_LAMPS, "WefabL2 2(ij,ab) - (ij,ba)", 2);
    global_dpd_->buf4_sort_axpy(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ij,ab) - (ij,ba)", -1);
    global_dpd_->buf4_close(&WL);

    //Sort
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ab) - (ij,ba)");
    global_dpd_->buf4_sort(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ji,ab) - (ji,ba)");
    global_dpd_->buf4_close(&WL);
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ji,ab) - (ji,ba)"); 
    //sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    //global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_axpy(&WL, &Z, 0.5);
    global_dpd_->buf4_close(&WL); 
    //global_dpd_->file2_close(&X1); 


    //Sort WefabL2 2(ij,ab) - (ij,ba) ->WefabL2 2(ij,ba) - (ij,ab)
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ab) - (ij,ba)");
    global_dpd_->buf4_sort(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ij,ba) - (ij,ab)"); 
    global_dpd_->buf4_close(&WL);

    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ba) - (ij,ab)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_axpy(&WL, &Z, 0.5);
    global_dpd_->dot23(&X1, &Z, &Y1new, 0, 0, 1, 1);	
    global_dpd_->buf4_close(&WL);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z);

    //Improve this part....
    //Build Gae*.Loovv -> Take it straight from L2???
    global_dpd_->file2_init(&z1, PSIF_CC_OEI, 0, 0, 1, "Z_ia");
    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, 0, 1, 1, "GAE");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract222(&X1, &GAE, &z1, 0, 1, 1, 0);
    global_dpd_->dot24(&z1, &D, &Y1new, 0, 0, 1, 1);
    global_dpd_->file2_close(&z1);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->contract424(&D, &GAE, &Z2, 3, 1, 0, 1, 0); 
    global_dpd_->dot13(&X1, &Z2, &Y1new, 0, 0, 1, 1);
    global_dpd_->file2_close(&GAE);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->file2_init(&z1, PSIF_CC_OEI, 0, 1, 0, "Z_ij");
    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, 0, 0, 0, "GMI");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract222(&X1, &GMI, &z1, 1, 0, 1, 0);
    global_dpd_->dot24(&z1, &D, &Y1new, 1, 0, -1, 1);
    global_dpd_->file2_close(&z1);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->contract424(&D, &GMI, &Z2, 1, 0, 1, 1, 0);  
    global_dpd_->dot13(&X1, &Z2, &Y1new, 0, 0, -1, 1);
    global_dpd_->file2_close(&GMI);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z2);


/*
    // Type-I L2 residual 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "LHX1Y1 I (2 Lijab - Lijba)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);

    global_dpd_->dot14(&X1, &L2, &Y1new, 0, 0, 1, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&X1);    
*/ 

    // Type-II L2 residual 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "LHX1Y1 Residual II");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract422(&L2, &X1, &Y1new, 0, 0, -1, 1);  
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "LHX1Y1 Residual IItest");
    global_dpd_->contract422(&L2, &X1, &Y1new, 0, 0, -1, 1);
    global_dpd_->file2_close(&X1);


    //# <O|L1(0)|[Hbar(0), X2]|phi^a_i>
    sprintf(lbl, "Z_%s_ME", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 0, 1, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    global_dpd_->dot24(&L1, &X2, &z1, 0, 0, 2, 0); 
    
    sprintf(lbl, "X_%s_IjbA (%5.3f)", pert, omega);
    global_dpd_->buf4_sort(&X2, PSIF_CC_LR, pqsr, 0, 5, lbl);
    global_dpd_->buf4_close(&X2);

    sprintf(lbl, "X_%s_IjbA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);  
    sprintf(lbl, "Z2_%s_ME", pert);
    global_dpd_->file2_init(&z2, PSIF_CC_TMP0, irrep, 0, 1, lbl);
    global_dpd_->dot24(&L1, &X2, &z2, 0, 0, -1, 0);
    global_dpd_->file2_axpy(&z2, &z1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z2);
 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot24(&z1, &D, &Y1new, 0, 0, 1, 1); 
    global_dpd_->buf4_close(&D); 
    global_dpd_->file2_close(&z1);

    sprintf(lbl, "Z_%s_MN", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 0, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z1, 0, 0, 1.0, 0.0);
    //global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); //Remove?
    global_dpd_->contract222(&z1, &L1, &Y1new, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&z1);

    sprintf(lbl, "Z_%s_AE", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 1, 1, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z1, 2, 2, -1.0, 0.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);
    global_dpd_->contract222(&L1, &z1, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&z1);
    global_dpd_->file2_close(&L1);

    // <O|L2(0)|[Hbar(0), X2]|0> 
    // Lijab * Xijab -> Lx_IJ //
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 0, 0, "Lx_IJ");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&L2, &X2, &lx, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 0, 0, "Lx_IJ");    //Lx or Lt??
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME"); 
    global_dpd_->contract222(&lx, &F, &Y1new, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->file2_close(&F);  


    // Lijab * Xijab -> Lx_AB 
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 1, 1, "Lx_AB");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&X2, &L2, &lx, 2, 2, 1.0, 0.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 1, 1, "Lx_AB");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");	
    global_dpd_->contract222(&F, &lx, &Y1new, 0, 1, -1.0, 1.0);   
    global_dpd_->file2_close(&lx);
    global_dpd_->file2_close(&F);

/*
    //Sorted  L (ib|ja))
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, psqr, 10, 10, "(2 LIjAb - LIjBa) (ib|ja)");
    global_dpd_->buf4_close(&L2);

    sprintf(lbl, "Lx_%s_ijab", pert);
    global_dpd_->buf4_init(&lx_ijab, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
    //global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    //sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    //global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);

    Y1_norm = 0;
    Y1_norm = global_dpd_->buf4_dot_self(&X2);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of the X2sorted.... %20.15f\n", Y1_norm);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "(2 LIjAb - LIjBa) (ib|ja)");

    Y1_norm = 0;
    Y1_norm = global_dpd_->buf4_dot_self(&L2);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of the L2sorted.... %20.15f\n", Y1_norm);


    //global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 1, 1.0, 0.0);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2test, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
 outfile->Printf("\tTest1.... \n");

    global_dpd_->contract444(&X2, &X2test, &lx_ijab, 1, 0, 1.0, 0.0);

 outfile->Printf("\tTest2.... \n");

    global_dpd_->buf4_close(&X2test);
 
    Y1_norm = 0;
    Y1_norm = global_dpd_->buf4_dot_self(&lx_ijab);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of the lx_iajb.... %20.15f\n", Y1_norm);

 outfile->Printf("\tTest1.... %20.15f\n");
*/

/*
    global_dpd_->buf4_init(&lx_ijab, PSIF_CC_LR, 0, 0, 0, 0, 0, 0, "LXijab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2test, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    global_dpd_->contract444(&X2test, &X2, &lx_ijab, 0, 0, 1, 1);

    Y1_norm = 0;
    Y1_norm = global_dpd_->buf4_dot_self(&lx_ijab);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of the lx_ijab.... %20.15f\n", Y1_norm);
*/


    //----------------------------------------------------------------------------------------------------- 

    //Sorted  L (ib|ja))
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, psqr, 10, 10, "(2 LIjAb - LIjBa) (ib|ja)");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, prqs, 10, 10, "(2 LIjAb - LIjBa) (ia|jb)");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, qspr, 10, 10, "(2 LIjAb - LIjBa) (jb|ia)");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");

    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 1, 1, 0);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); //Compute this part out of core
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, -1, 1);


    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);


    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)"); //Check ths one!!!

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_2");

    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 1, 1, 0);
 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); //Compute this part out of core
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf (am,fe)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (am,fe)"); //Compute this part out of core
    
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, -1, 1);  
 

    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);    
    global_dpd_->buf4_close(&W);

   //Here

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_init(&Hx_ijab, PSIF_CC_LR, irrep, 0, 11, 0, 11, 0, "Hx_ijab");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); //Compute this part out of core
   
    global_dpd_->contract444(&X2, &W, &Hx_ijab, 0, 0, 1, 0);
    global_dpd_->contract442(&Hx_ijab, &L2, &Y1new, 3, 3, -1, 1);

    global_dpd_->buf4_close(&Hx_ijab);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&W);

    //------------------------------------------------------------------------------------------------------

    //This is computed ealier
    /*
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); 
    global_dpd_->buf4_scmcopy(&W, PSIF_CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);
    */


    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, irrep, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl); 
    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, 0, 10, 10, 10, 10, 0, "LXiajb_new");

    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 0, 1, 0);
    
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx_iajb);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, 0, 10, 10, 10, 10, 0, "LXiajb_new");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)");	

    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, 1, 1);

    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); //Make it out of core
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep, 1, 1, "Lx_AB");
   
    global_dpd_->dot13(&lx, &W, &Y1new, 1, 0, 1.0, 1.0); 

    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&lx);    

  
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_init(&lx_ijab, PSIF_CC_LR, irrep, 0, 0, 0, 0, 0, "Lx_ijkl");
    
    global_dpd_->contract444(&L2, &X2, &lx_ijab, 0, 0, 1, 0);

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx_ijab); 


    global_dpd_->buf4_init(&lx_ijab, PSIF_CC_LR, irrep, 0, 0, 0, 0, 0, "Lx_ijkl");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    	    	    
    global_dpd_->contract442(&lx_ijab, &W, &Y1new, 1, 3, 1, 1);

    global_dpd_->buf4_close(&lx_ijab);
    global_dpd_->buf4_close(&W); 


    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_2"); //I am reusing here
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");

    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, 1, 1);

    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_3");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");

    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 0, 1, 0);
    
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx_iajb);

    //Sort WMnIE (Mn,eI) -> WMnIE (nM,eI)
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 0, 11, "WMnIe (nM,eI)");

    //global_dpd_->file2_init(&test, PSIF_CC_OEI, irrep, 0, 1, "test");
    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_3");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (nM,eI)");
    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep, 0, 0, "Lx_IJ");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->dot14(&lx, &W, &Y1new, 1, 0, -1.0, 1.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->buf4_close(&W);
   
    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, irrep, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 1, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx_iajb);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)"); 
    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, -1, 1);   
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);


/*
    Y1_norm = 0;
    Y1_norm = global_dpd_->file2_dot_self(&Y1new);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of the Y1new part1i_Final.... %20.15f\n", Y1_norm);
*/


    global_dpd_->file2_close(&Y1new);


    return;
}

}  // namespace ccresponse
}  // namespace psi