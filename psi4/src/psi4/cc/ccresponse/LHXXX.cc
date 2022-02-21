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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

double LHXXX_p1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, z, z1, l1, xl;
    dpdbuf4 D; 
    char lbl[32];
    double Y1_norm;

//outfile->Printf("\n\tirrep_x %d", irrep_x);
//outfile->Printf("\n\tirrep_y %d", irrep_y);
//outfile->Printf("\n\tirrep_z %d", irrep_z);

/*
    sprintf(lbl, "z_%s_IA", pert_x);
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_x, 0, 1, lbl);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->contract422(&D, &X1, &z, 0, 0, 1, 0); 
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);
*/
    sprintf(lbl, "z_IJ%s%s", pert_x, pert_z);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep_x^irrep_z, 0, 0, lbl);
    global_dpd_->file2_scm(&z1, 0);	

    sprintf(lbl, "DX_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&z, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    //global_dpd_->contract222(&X1, &z, &z1, 0, 0, 1, 0);
    global_dpd_->contract222(&z, &X1, &z1, 0, 0, 1, 0);
    global_dpd_->file2_close(&z);
    global_dpd_->file2_close(&X1);

    sprintf(lbl, "xl_%s_mi (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&xl, PSIF_CC_OEI, irrep_y, 0, 0, lbl);

/*
    sprintf(lbl, "XL_%s_IJ", pert_y);
    global_dpd_->file2_init(&lx, PSIF_CC_TMP0, irrep_y, 0, 0, lbl);

    global_dpd_->file2_init(&l1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&l1, &X1, &lx, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&l1);
*/
    result = global_dpd_->file2_dot(&xl, &z1);

    global_dpd_->file2_close(&z1);
    global_dpd_->file2_close(&xl);
    
    return result;
}

double LHX1X1X1_p1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1; 
    dpdbuf4 l2, W, WX, LX, Z; 
    char lbl[32];
    double Y1_norm;

/*
    sprintf(lbl, "XW_%s (ij,kl) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract424(&W, &X1, &XW, 3, 1, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    sprintf(lbl, "XW_%s (kl,ij) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&XW, PSIF_CC_TMP0, rspq, 0, 0, lbl);
    global_dpd_->buf4_close(&XW);
*/

/*
    sprintf(lbl, "XL_%s (ij,al) (%5.3f)", pert_z, omega_z);    
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 11, 0, 11, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract424(&l2, &X1, &XL, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/

    sprintf(lbl, "Z(ij,kl)%s%s", pert_y, pert_z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_z, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "LX_%s_jiak (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 11, 0, 11, 0, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &LX, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&LX);
    global_dpd_->file2_close(&X1);

//    sprintf(lbl, "XW_%s (kl,ij) (%5.3f)", pert_x, omega_x);
//    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "WX_%s_ijkl (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&WX, PSIF_CC_LR, irrep_x, 0, 0, 0, 0, 0, lbl);

    result = global_dpd_->buf4_dot(&Z, &WX);

    global_dpd_->buf4_close(&WX);
    global_dpd_->buf4_close(&Z);
    
    return result; 
}

double LHX1X1X1_p2(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1; 
    dpdbuf4 l2, Z, WX, LX, W;
    char lbl[32];

/*
    sprintf(lbl, "XW_%s (ij,kl) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    //global_dpd_->contract244(&X1, &W, &XW, 1, 3, 0, 1, 0);
    global_dpd_->contract424(&W, &X1, &XW, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&XW);
*/

    //FInd intermediate for XL!!!!

/*
    sprintf(lbl, "XL_%s (ij,ka) (%5.3f)", pert_z, omega_z); 
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&X1, &l2, &XL, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/

    sprintf(lbl, "Z(ij,kl)%s%s", pert_z, pert_y);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_z^irrep_y, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 10, 0, 10, 0, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    //global_dpd_->contract244(&X1, &XL, &Z, 1, 3, 0, 1, 0);
    global_dpd_->contract424(&LX, &X1, &Z, 3, 1, 1, 1, 0);
    global_dpd_->buf4_close(&LX);
    global_dpd_->file2_close(&X1);

    //sprintf(lbl, "XW_%s (ij,kl) (%5.3f)", pert_x, omega_x); 
    //global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);
    sprintf(lbl, "WX_%s_ijkl (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&WX, PSIF_CC_LR, irrep_x, 0, 0, 0, 0, 0, lbl);

    result = global_dpd_->buf4_dot(&Z, &WX);

    global_dpd_->buf4_close(&WX);
    global_dpd_->buf4_close(&Z);

    return result;
}

double LHX1X1X1_p3(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1; 
    dpdbuf4 l2, Z, WX, LX, W;
    char lbl[32];
  
/*  
    sprintf(lbl, "XW_%s (aj,bl) (%5.3f)", pert_x, omega_x); 
    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, irrep_x, 11, 11, 11, 11, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract424(&W, &X1, &XW, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);
*/

    sprintf(lbl, "Z(al,kl)%s%s", pert_y, pert_x);	
    //sprintf(lbl, "Z_%s (al,kj) (%5.3f)", pert_y, omega_y);    
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_x, 11, 0, 11, 0, 0, lbl);

    sprintf(lbl, "WX_%s_fibm (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&WX, PSIF_CC_LR, irrep_x, 11, 11, 11, 11, 0, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &WX, &Z, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&WX);
    global_dpd_->file2_close(&X1);

    sprintf(lbl, "Z(jk,la)%s%s", pert_y, pert_x); 
    //sprintf(lbl, "Z_%s (jk,la) (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, srqp, 0, 10, lbl);
    global_dpd_->buf4_close(&Z);

/*
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "XL(ij,la)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&X1, &l2, &XL, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 10, 0, 10, 0, lbl);

    //sprintf(lbl, "Z_%s (jk,la) (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "Z(jk,la)%s%s", pert_y, pert_x); 
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_x, 0, 10, 0, 10, 0, lbl);

    result = global_dpd_->buf4_dot(&LX, &Z);

    global_dpd_->buf4_close(&LX);
    global_dpd_->buf4_close(&Z);

    return result;
}

double LHX1X1X1_p4(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1; 
    dpdbuf4 l2, Z, XW, LX, W;
    char lbl[32];

/*
    sprintf(lbl, "XW_%s (aj,lb) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, irrep_x, 11, 10, 11, 10, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract244(&X1, &W, &XW, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);
*/
    sprintf(lbl, "Z(al,kl)%s%s", pert_x, pert_y);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_x^irrep_y, 11, 0, 11, 0, 0, lbl);

    sprintf(lbl, "WAmEfX1_%s (fj,ma) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XW, PSIF_CC_LR, irrep_x, 11, 10, 11, 10, 0, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract424(&XW, &X1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&XW);
    global_dpd_->file2_close(&X1);

    sprintf(lbl, "Z(jk,al)%s%s", pert_x, pert_y);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, rspq, 0, 11, lbl);
    global_dpd_->buf4_close(&Z);

/*
    sprintf(lbl, "XL_%s (ij,al) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 11, 0, 11, 0, "XL(ij,al)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract424(&l2, &X1, &XL, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/

    sprintf(lbl, "LX_%s_jiak (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 11, 0, 11, 0, lbl);

    //sprintf(lbl, "Z_%s (jk,al) (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "Z(jk,al)%s%s", pert_x, pert_y);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_x^irrep_y, 0, 11, 0, 11, 0, lbl);

    result = global_dpd_->buf4_dot(&LX, &Z);

    global_dpd_->buf4_close(&LX);
    global_dpd_->buf4_close(&Z);

    return result;
}

double LHX2X1X1_p1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1; 
    dpdbuf4 l2, X2, Z, XD, LX, D;
    char lbl[32];
    double Y2_norm;

/*
    sprintf(lbl, "XL_%s (ij,al) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_x, 0, 11, 0, 11, 0, lbl);


    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract424(&l2, &X1, &XL, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
    global_dpd_->buf4_close(&XL);
*/

/*
    sprintf(lbl, "XD_%s (ij,al) (%5.3f)", pert_y, omega_y); 
    global_dpd_->buf4_init(&XD, PSIF_CC_TMP0, omega_y, 0, 11, 0, 11, 0, lbl);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract424(&D, &X1, &XD, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);
    //global_dpd_->buf4_close(&XD);
*/
    sprintf(lbl, "XD_%s_ijam (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&XD, PSIF_CC_LR, irrep_y, 0, 11, 0, 11, 0, lbl);

    sprintf(lbl, "XD_%s (ia,jl) (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_sort(&XD, PSIF_CC_TMP1, prqs, 10, 0, lbl);
    global_dpd_->buf4_close(&XD);

    //sprintf(lbl, "Z_%s (ja,kl) (%5.3f)", pert_z, omega_z);
    sprintf(lbl, "Z(ja,kl)%s%s", pert_z, pert_y); 
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, irrep_z^irrep_y, 10, 0, 10, 0, 0, lbl);
    global_dpd_->buf4_scm(&Z, 0);

    sprintf(lbl, "XD_%s (ia,jl) (%5.3f)", pert_y, omega_y);  
    global_dpd_->buf4_init(&XD, PSIF_CC_TMP1, irrep_y, 10, 0, 10, 0, 0, lbl);
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &XD, &Z, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&XD);

    //sprintf(lbl, "Z_%s (jl,ak) (%5.3f)", pert_z, omega_z);
    sprintf(lbl, "Z(jl,ak)%s%s", pert_z, pert_y);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, psqr, 0, 11, lbl);
    global_dpd_->buf4_close(&Z);

    //sprintf(lbl, "Z_%s (jl,ak) (%5.3f)", pert_z, omega_z); 
    sprintf(lbl, "Z(jl,ak)%s%s", pert_z, pert_y);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, irrep_z^irrep_y, 0, 11, 0, 11, 0, lbl);
    //sprintf(lbl, "XL_%s (ij,al) (%5.3f)", pert_x, omega_x);
    //global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_x, 0, 11, 0, 11, 0, lbl);  

    sprintf(lbl, "LX_%s_jiak (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_x, 0, 11, 0, 11, 0, lbl);

    result = global_dpd_->buf4_dot(&Z, &LX);
     
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&LX);

    return result;	
}

double LHX2X1X1_p2(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, lx, z, dx;
    dpdbuf4 l2, X2, D;
    char lbl[32];

//    sprintf(lbl, "xl_%s_AB (%5.3f)", pert_x, omega_x);
//    global_dpd_->file2_init(&xl_ab, PSIF_CC_OEI, irrep_x, 1, 1, lbl);

/*
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&X2, &l2, &xl_ab, 2, 2, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&l2);
*/

    //sprintf(lbl, "z_%s_IA (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "z_IA%s%s", pert_y, pert_x);
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_y^irrep_x, 0, 1, lbl);

    sprintf(lbl, "LX_%s_AB (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep_x, 1, 1, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    //global_dpd_->contract222(&X1, &xl_ab, &z, 0, 0, 1, 0);    
    global_dpd_->contract222(&X1, &lx, &z, 0, 1, 1, 0); 
    global_dpd_->file2_close(&X1); 
    global_dpd_->file2_close(&lx);	   

    sprintf(lbl, "DX_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&dx, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    
    /*
    sprintf(lbl, "dx_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&dx, PSIF_CC_OEI, irrep_z, 0, 1, "dx_IA");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)"); 
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->contract422(&D, &X1, &dx, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);
    */

    result = global_dpd_->file2_dot(&z, &dx);

    global_dpd_->file2_close(&z);
    global_dpd_->file2_close(&dx);

    return result;
}

double LHX2X1X1_p3(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, lx, dx, z;
    dpdbuf4 l2, X2, D;
    char lbl[32];

/*
    sprintf(lbl, "xl_%s_IJ (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&xl_ij, PSIF_CC_OEI, irrep_x, 0, 0, "xl_IJ");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&X2, &l2, &xl_ij, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&l2);
*/

/*  
    sprintf(lbl, "dx_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&dx, PSIF_CC_OEI, irrep_y, 0, 1, "dx_IA");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract422(&D, &X1, &dx, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);
*/

    //sprintf(lbl, "xl_%s_IJ (%5.3f)", pert_z, omega_z);
    sprintf(lbl, "xl_IJ%s%s", pert_z, pert_y);
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_z^irrep_y, 0, 0, lbl);

    sprintf(lbl, "DX_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&dx, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->contract222(&X1, &dx, &z, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&dx);

    sprintf(lbl, "Lx_%s_IJ (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep_x, 0, 0, lbl);

    result = global_dpd_->file2_dot(&lx, &z);

    global_dpd_->file2_close(&z);
    global_dpd_->file2_close(&lx);

    return result;
}

double LHX2X1X1_p4(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1;
    dpdbuf4 l2, X2, Z, XD, LX, D;
    char lbl[32];
    double Y2_norm;

/*
    sprintf(lbl, "XD_%s (ij,kl) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XD, PSIF_CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&X2, &D, &XD, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);
*/

/*
    sprintf(lbl, "XL_%s (ij,al) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 11, 0, 11, 0, "XL(ij,al)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract424(&l2, &X1, &XL, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/

    //sprintf(lbl, "Z_%s (ij,kl) (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "Z(ij,kl)%s%s", pert_y, pert_z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_z, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "LX_%s_jiak (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 11, 0, 11, 0, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &LX, &Z, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&LX);
    global_dpd_->file2_close(&X1);

    sprintf(lbl, "X2D_%s_ijlk (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XD, PSIF_CC_LR, irrep_x, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "X2D_%s_ijkl (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&XD, PSIF_CC_LR, pqsr, 0, 0, lbl);
    global_dpd_->buf4_close(&XD);

    sprintf(lbl, "X2D_%s_ijkl (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XD, PSIF_CC_LR, irrep_x, 0, 0, 0, 0, 0, lbl);

    result = global_dpd_->buf4_dot(&Z, &XD);

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&XD);

    return result;
}

double LHX2X1X1_p5(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1;
    dpdbuf4 l2, X2, Z, XD, LX, D;
    char lbl[32];

    sprintf(lbl, "XD_%s (ia,jk) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XD, PSIF_CC_LR, irrep_x, 10, 0, 10, 0, 0, lbl);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D (ia|jb)"); 
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->contract424(&D, &X1, &XD, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);

    //sprintf(lbl, "Z_%s (ja,kl) (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "Z(ja,kl)%s%s", pert_y, pert_x);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_x, 10, 0, 10, 0, 0, lbl);

    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &XD, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&XD);

    //sprintf(lbl, "Z_%s (jk,al) (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "Z(jk,al)%s%s", pert_y, pert_x);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, prqs, 0, 11, lbl);
    global_dpd_->buf4_close(&Z);

/*
    sprintf(lbl, "XL_%s (ij,al) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 11, 0, 11, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract424(&l2, &X1, &XL, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);

    sprintf(lbl, "XL_%s (il,aj) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_sort(&XL, PSIF_CC_TMP0, psrq, 0, 11, lbl);
    global_dpd_->buf4_close(&XL);
    
    sprintf(lbl, "XL_%s (il,aj) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 11, 0, 11, 0, lbl);
*/

    sprintf(lbl, "LX_%s_jiak (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 11, 0, 11, 0, lbl);

    sprintf(lbl, "XL_%s (il,aj) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_sort(&LX, PSIF_CC_TMP0, psrq, 0, 11, lbl);
    global_dpd_->buf4_close(&LX);

    sprintf(lbl, "XL_%s (il,aj) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_TMP0, irrep_z, 0, 11, 0, 11, 0, lbl);

    //sprintf(lbl, "Z_%s (jk,al) (%5.3f)", pert_y, omega_y); 
    sprintf(lbl, "Z(jk,al)%s%s", pert_y, pert_x);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_x, 0, 11, 0, 11, 0, lbl);

    result = global_dpd_->buf4_dot(&Z, &LX);

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&LX);

    return result;
}

double LHX2X1X1_p6(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1;
    dpdbuf4 l2, X2, Z, XD, LX, D;
    char lbl[32];

/*
    sprintf(lbl, "XD_%s (ia,jl) (%5.3f)", pert_x, omega_x); 
    global_dpd_->buf4_init(&XD, PSIF_CC_TMP0, irrep_x, 10, 0, 10, 0, 0, lbl);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D (ia|jb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->contract424(&D, &X1, &XD, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);

*/

    //sprintf(lbl, "Z_%s (ja,kl) (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "Z(ja,kl)%s%s", pert_y, pert_x); 
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_x, 10, 0, 10, 0, 0, lbl);

    sprintf(lbl, "XD_%s (ia,jk) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XD, PSIF_CC_LR, irrep_x, 10, 0, 10, 0, 0, lbl);
    sprintf(lbl, "XD_%s (ja,jl) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&XD, PSIF_CC_TMP0, rqps, 10, 0, lbl);
    global_dpd_->buf4_close(&XD);

    sprintf(lbl, "XD_%s (ja,jl) (%5.3f)", pert_x, omega_x); 
    global_dpd_->buf4_init(&XD, PSIF_CC_TMP0, irrep_x, 10, 0, 10, 0, 0, lbl);
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &XD, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&XD);

    //sprintf(lbl, "Z_%s (jk,al) (%5.3f)", pert_y, omega_y); 
    sprintf(lbl, "Z(jk,al)%s%s", pert_y, pert_x); 
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, prqs, 0, 11, lbl);
    global_dpd_->buf4_close(&Z);

/*
    sprintf(lbl, "XL_%s (ij,la) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&X1, &l2, &XL, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "LX_%s (il,aj) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_sort(&LX, PSIF_CC_TMP0, prsq, 0, 11, lbl);
    global_dpd_->buf4_close(&LX);

    sprintf(lbl, "LX_%s (il,aj) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_TMP0, irrep_z, 0, 11, 0, 11, 0, lbl);
    sprintf(lbl, "Z(jk,al)%s%s", pert_y, pert_x);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_x, 0, 11, 0, 11, 0, lbl);

    result = global_dpd_->buf4_dot(&Z, &LX);

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&LX);

    return result;
}

double LHX2X1X1_p7(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1;
    dpdbuf4 l2, X2, Z, XD, LX, D;
    char lbl[32];

    sprintf(lbl, "XD_%s (ij,al) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&XD, PSIF_CC_TMP0, irrep_x, 0, 11, 0, 11, 0, lbl);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->contract424(&D, &X1, &XD, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);

    //sprintf(lbl, "Z_%s (ij,kl) (%5.3f)", pert_y, omega_y);
    sprintf(lbl, "Z(ij,kl)%s%s", pert_y, pert_x);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y^irrep_x, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl); 
    //global_dpd_->contract244(&X1, &XD, &Z, 1, 2, 1, 1, 0);
    global_dpd_->contract244(&X1, &XD, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&XD);
    global_dpd_->file2_close(&X1); 

/*
    sprintf(lbl, "XL_%s (ij,kl) (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, irrep_z, 0, 0, 0, 0, 0, lbl);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract444(&X2, &l2, &XL, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&l2);
*/

    sprintf(lbl, "LX_%s_ijkl (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep_z, 0, 0, 0, 0, 0, lbl);
  
    result = 2.0 * global_dpd_->buf4_dot(&Z, &LX);

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&LX);

    return result;
}

double LHXXX(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		     const char *pert_z, int irrep_z, double omega_z) {

    double hyper = 0.0;

    /*** <L1(0)|[[[H_bar,X1(A)],X1(B)],X1(C)]|0> ***/ 

    hyper -= LHXXX_p1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z); 
    hyper -= LHXXX_p1(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y); 
    hyper -= LHXXX_p1(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x); 
    hyper -= LHXXX_p1(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z); 
    hyper -= LHXXX_p1(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);
    hyper -= LHXXX_p1(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x); 

   // <L2(0)|[[[H_bar,X1(A)],X1(B)],X1(C)]|0> 

    hyper += LHX1X1X1_p1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper += LHX1X1X1_p1(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);
    hyper += LHX1X1X1_p1(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);

    hyper += LHX1X1X1_p2(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper += LHX1X1X1_p2(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);
    hyper += LHX1X1X1_p2(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);    

    hyper -= LHX1X1X1_p3(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper -= LHX1X1X1_p3(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);
    hyper -= LHX1X1X1_p3(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);

    hyper -= LHX1X1X1_p4(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper -= LHX1X1X1_p4(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);
    hyper -= LHX1X1X1_p4(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);

    // <L2(0)|[[[H_bar,X2(A)],X1(B)],X1(C)]|0> 
   
    hyper -= LHX2X1X1_p1(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x); 
    hyper -= LHX2X1X1_p1(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);

    hyper -= LHX2X1X1_p2(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);
    hyper -= LHX2X1X1_p2(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z); 

    hyper -= LHX2X1X1_p3(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper -= LHX2X1X1_p3(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);

    hyper += LHX2X1X1_p4(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z); 

    hyper += LHX2X1X1_p5(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z);
    hyper += LHX2X1X1_p5(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);   

    hyper += LHX2X1X1_p6(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z); 
    hyper += LHX2X1X1_p6(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);  

    hyper += LHX2X1X1_p7(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);

    // <L2(0)|[[[H_bar,X1(A)],X2(B)],X1(C)]|0> 

    hyper -= LHX2X1X1_p1(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);
    hyper -= LHX2X1X1_p1(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);
    
    hyper -= LHX2X1X1_p2(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);
    hyper -= LHX2X1X1_p2(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z);

    hyper -= LHX2X1X1_p3(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z);
    hyper -= LHX2X1X1_p3(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);

    hyper += LHX2X1X1_p4(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z);

    hyper += LHX2X1X1_p5(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper += LHX2X1X1_p5(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);

    hyper += LHX2X1X1_p6(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper += LHX2X1X1_p6(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);

    hyper += LHX2X1X1_p7(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);

    // <L2(0)|[[[H_bar,X1(A)],X1(B)],X2(C)]|0> 

    hyper -= LHX2X1X1_p1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper -= LHX2X1X1_p1(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z);

    hyper -= LHX2X1X1_p2(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);
    hyper -= LHX2X1X1_p2(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);

    hyper -= LHX2X1X1_p3(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);
    hyper -= LHX2X1X1_p3(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);

    hyper += LHX2X1X1_p4(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);

    hyper += LHX2X1X1_p5(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);
    hyper += LHX2X1X1_p5(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);

    hyper += LHX2X1X1_p6(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);
    hyper += LHX2X1X1_p6(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);

    hyper += LHX2X1X1_p7(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);

    return hyper;
}

}  // namespace ccresponse
}  // namespace psi
