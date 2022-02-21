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
      
    Computes LCXX contributions to Quadratic Response Functions.
    
    Author: Monika Kodrycka
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

double LCXX(const char *pert_c, int irrep_c, double omega_c, const char *pert_x, int irrep_x, double omega_x,
		     const char *pert_y, int irrep_y, double omega_y) {

    double hyper = 0.0;
    dpdfile2 Y1, yt, LX;
    dpdfile2 X1, mu1, z, z1, l1, mu, lt, xc, xl;
    dpdbuf4 X2, Y2, l2, mu2, z2, Z, t2;
    char lbl[32];
    double Y1_norm;

    /*** L1 * MuBAR * X1 * X1 + L2 * MuBAR * X1 * X2
                + L2 * MuBAR * X2 * X1 ***/

/*
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, 0, 0, 0, "z_IJ");

    global_dpd_->file2_init(&l1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&l1, &X1, &z, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&l1);
*/

/*
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, 0, 0, 0, "z_IJ");

    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    //global_dpd_->contract442(&l2, &X2, &z, 0, 0, 1.0, 1.0);
    global_dpd_->contract442(&X2, &l2, &z, 0, 0, 1.0, 0.0);
*/

    sprintf(lbl, "z_IJ%s", pert_y); 
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_y, 0, 0, lbl);
    global_dpd_->file2_scm(&z, 0); 

    sprintf(lbl, "Lx_%s_IJ (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&LX, PSIF_CC_OEI, irrep_y, 0, 0, lbl);

    global_dpd_->file2_axpy(&LX, &z, 1, 1);

    global_dpd_->file2_close(&LX);

    sprintf(lbl, "xl_%s_mi (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&xl, PSIF_CC_OEI, irrep_y, 0, 0, lbl);

    global_dpd_->file2_axpy(&xl, &z, 1, 0);
    global_dpd_->file2_close(&xl);

    sprintf(lbl, "z_IJ%s%s", pert_c, pert_x);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP1, irrep_c^irrep_x, 0, 0, lbl);

    sprintf(lbl, "%s_IA", pert_c);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    //global_dpd_->contract222(&X1, &mu1, &z1, 0, 0, 1, 0);
    global_dpd_->contract222(&mu1, &X1, &z1, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&mu1);  

    hyper -= global_dpd_->file2_dot(&z, &z1);

    global_dpd_->file2_close(&z); 
    global_dpd_->file2_close(&z1); 

    sprintf(lbl, "z_IJ%s%s", pert_c, pert_y); 
    global_dpd_->file2_init(&z, PSIF_CC_TMP2, irrep_c^irrep_y, 0, 0, lbl);

    sprintf(lbl, "%s_IA", pert_c);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&mu1, &X1, &z, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&mu1);

    //global_dpd_->file2_init(&z1, PSIF_CC_TMP1, 0, 0, 0, "z_IJ");

/*
    global_dpd_->file2_init(&l1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->contract222(&X1, &l1, &z1, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&l1);
*/

    sprintf(lbl, "xl_%s_mi (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&xl, PSIF_CC_OEI, irrep_x, 0, 0, lbl);

    hyper -= global_dpd_->file2_dot(&z, &xl);

    global_dpd_->file2_close(&z);
    global_dpd_->file2_close(&xl);

/*
    global_dpd_->buf4_init(&z2, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Ij,kb)");

    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &l2, &z2, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/


    sprintf(lbl, "z_IA%s%s", pert_x, pert_y);	
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_x^irrep_y, 0, 1, lbl);

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert_x, omega_x);  
    global_dpd_->buf4_init(&z2, PSIF_CC_LR, irrep_x, 0, 10, 0, 10, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&z2, &X2, &z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&X2); 
    global_dpd_->buf4_close(&z2);

    sprintf(lbl, "%s_IA", pert_c);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);

    hyper -= global_dpd_->file2_dot(&z, &mu1);
     
    global_dpd_->file2_close(&z);
    global_dpd_->file2_close(&mu1); 

    sprintf(lbl, "z_IJ%s%s", pert_c, pert_y);
    global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_c^irrep_y, 0, 0, lbl); ///I changed here....

    sprintf(lbl, "%s_IA", pert_c);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&X1, &mu1, &z, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&mu1);

    sprintf(lbl, "Lx_%s_IJ (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&z1, PSIF_CC_OEI, irrep_x, 0, 0, lbl);

/*
    global_dpd_->file2_init(&z1, PSIF_CC_TMP1, 0, 0, 0, "z_IJ");

    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&X2, &l2, &z1, 0, 0, 1.0, 0.0);
*/

    hyper -= global_dpd_->file2_dot(&z, &z1);

    global_dpd_->file2_close(&z);
    global_dpd_->file2_close(&z1);


/*
    global_dpd_->buf4_init(&z2, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Ij,kb)");

    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &l2, &z2, 1, 2, 1, 1, 0);

    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);
*/

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&z2, PSIF_CC_LR, irrep_y, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "z_IA%s%s", pert_y, pert_x); 
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep_y^irrep_x, 0, 1, lbl); 
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&z2, &X2, &z1, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&z2);
    global_dpd_->buf4_close(&X2);

    sprintf(lbl, "%s_IA", pert_c);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);    
    hyper -= global_dpd_->file2_dot(&z1, &mu1);

    global_dpd_->file2_close(&z1);
    global_dpd_->file2_close(&mu1);

    return hyper;
}

}  // namespace ccresponse
}  // namespace psi
