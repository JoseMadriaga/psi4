/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 Copyright (c) 2007-2021 The Psi4 Developers.
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
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

double pseudopolar_Y(const char *pert, int irrep, double omega) {
    dpdfile2 Y1, mubar1;
    dpdbuf4 Y2, mubar2;
    char lbl[32];
    double polar1, polar2, Y1_norm;

    sprintf(lbl, "%sBAR_IA", pert);
    global_dpd_->file2_init(&mubar1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl);

    polar1 = 2.0 * global_dpd_->file2_dot(&mubar1, &Y1);
    global_dpd_->file2_close(&mubar1);
    global_dpd_->file2_close(&Y1);

    //Change to Y2
    sprintf(lbl, "%sBAR_IjAb", pert);
    global_dpd_->buf4_init(&mubar2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    polar2 = global_dpd_->buf4_dot(&mubar2, &Y2);
    global_dpd_->buf4_close(&mubar2);
    global_dpd_->buf4_close(&Y2);

    /*   if(params.print & 2) { */
    /*     outfile->Printf( "\tpolar1 = %20.12f\n", -2.0*polar1); */
    /*     outfile->Printf( "\tpolar2 = %20.12f\n", -2.0*polar2); */
    /*   } */

    return polar1 + polar2;
}

}  // namespace ccresponse
}  // namespace psi