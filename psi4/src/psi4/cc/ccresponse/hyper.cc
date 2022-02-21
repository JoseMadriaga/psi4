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
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <sstream>

#include "psi4/libpsi4util/process.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"

#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

void pertbar(const char *pert, int irrep, int anti);
void compute_X(const char *pert, int irrep, double omega);
void compute_Y(const char *pert, int irrep, double omega);
void quadresp(double *tensor, double A, double B, const char *pert_x, int x_irrep, double omega_x, const char *pert_y,
             int y_irrep, double omega_y, const char *pert_z, int z_irrep, double omega_z);

void hyper() {
    double ***tensor;
    char **cartcomp, pert[32], pert_x[32], pert_y[32], pert_z[32];
    int alpha, beta, gamma, i;
    double omega_nm, omega_ev, omega_cm, *trace;
    char lbl[32];
    double hyper_ave;

    cartcomp = (char **)malloc(3 * sizeof(char *));
    cartcomp[0] = strdup("X");
    cartcomp[1] = strdup("Y");
    cartcomp[2] = strdup("Z");

    //tensor = (double ***)malloc(params.nomega * sizeof(double **));
    //for (i = 0; i < params.nomega; i++) tensor[i] = block_matrix(3, 3, 3);

     tensor = (double ***)malloc(3 * sizeof(double **));
        for (i = 0; i < 3; i++) tensor[i] = block_matrix(3, 3);

//Do I need to assume that user provides 3 frequencies?
    for (i = 0; i < params.nomega; i++) {
    	sprintf(lbl, "<<Mu;Mu;Mu>>_(%5.3f)", params.omega[i]);
        outfile->Printf("\t-----------------------------------------\n");
        if (!params.restart || !psio_tocscan(PSIF_CC_INFO, lbl)) {
            for (alpha = 0; alpha < 3; alpha++) {

		zero_mat(tensor[alpha], 3, 3);

            	sprintf(pert, "Mu_%1s", cartcomp[alpha]);
                pertbar(pert, moinfo.mu_irreps[alpha], 0);
                compute_X(pert, moinfo.mu_irreps[alpha], params.omega[i]);
                compute_Y(pert, moinfo.mu_irreps[alpha], params.omega[i]);

		//This block is super redundant!!!!!!!!
                /*
		if (params.omega[i] != 0.0){
		    compute_X(pert, moinfo.mu_irreps[alpha], -params.omega[i]);
		    compute_Y(pert, moinfo.mu_irreps[alpha], -params.omega[i]);
		
		}
		*/
		
            }
        }

     }

     //outfile->Printf("\n\tComputing %s tensor.\n", lbl);
     for (alpha = 0; alpha < 3; alpha++) {
     	for (beta = 0; beta < 3; beta++) {
            for (gamma = 0; gamma < 3; gamma++) {
                sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
                sprintf(pert_y, "Mu_%1s", cartcomp[beta]);
                sprintf(pert_z, "Mu_%1s", cartcomp[gamma]);
		    
		quadresp(&tensor[alpha][beta][gamma], 0.5, 0.0, pert_x, moinfo.mu_irreps[alpha], params.omega[0],
		    pert_y, moinfo.mu_irreps[beta], params.omega[1],
                    pert_z, moinfo.mu_irreps[gamma], params.omega[2]);

	 	if (params.omega[1] == 0.0 && params.omega[0] != 0.0 && params.omega[2] != 0.0) {
                    quadresp(&tensor[alpha][beta][gamma], 0.5, 1.0, pert_x, moinfo.mu_irreps[alpha], -params.omega[0],
                    pert_y, moinfo.mu_irreps[beta], params.omega[1],
                    pert_z, moinfo.mu_irreps[gamma], -params.omega[2]);
	        }

                else if (params.omega[0] == 0.0 && params.omega[1] != 0.0 && params.omega[2] != 0.0) {
                    quadresp(&tensor[alpha][beta][gamma], 0.5, 1.0, pert_x, moinfo.mu_irreps[alpha], params.omega[0],
                    pert_y, moinfo.mu_irreps[beta], -params.omega[1],
                    pert_z, moinfo.mu_irreps[gamma], -params.omega[2]);
                }

		else if (params.omega[0] != 0.0 && params.omega[1] != 0.0 && params.omega[2] == 0.0) {
                    quadresp(&tensor[alpha][beta][gamma], 0.5, 1.0, pert_x, moinfo.mu_irreps[alpha], -params.omega[0],
                    pert_y, moinfo.mu_irreps[beta], -params.omega[1],
                    pert_z, moinfo.mu_irreps[gamma], params.omega[2]);
                }

		else {
                    quadresp(&tensor[alpha][beta][gamma], 0.5, 1.0, pert_x, moinfo.mu_irreps[alpha], -params.omega[0],
                    pert_y, moinfo.mu_irreps[beta], -params.omega[1],
                    pert_z, moinfo.mu_irreps[gamma], -params.omega[2]);
		}

		outfile->Printf("\tMu_%1s, Mu_%1s, Mu_%1s, %5.10f\n", cartcomp[alpha], cartcomp[beta], cartcomp[gamma], tensor[alpha][beta][gamma]);
            } 
        }
    }

	outfile->Printf("\n\n \t      CCSD Hyperpolarizability\n");
	outfile->Printf("  -------------------------------------------------------------------------\n");
	outfile->Printf("   Evaluated at: \n");
        outfile->Printf("   \t omega_1 = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[0],
                            omega_nm, omega_ev, omega_cm);
	outfile->Printf("   \t omega_2 = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[1],
                            omega_nm, omega_ev, omega_cm);	
        outfile->Printf("   \t omega_3 = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[2],
                            omega_nm, omega_ev, omega_cm);

        for (i = 0; i < params.nomega; i++) {
        	outfile->Printf("  -------------------------------------------------------------------------\n");
        	mat_print(tensor[i], 3, 3, "outfile");
	}

	// Compute the hyperpolarizability average
	hyper_ave = 0.0;
	for (i=0; i<3; i++){
		hyper_ave += (tensor[i][i][2] + tensor[i][2][i] + tensor[2][i][i]) / 5;
	}

	outfile->Printf("\n\tbeta_ave = %20.12f a.u.\n", hyper_ave);         
	   
    //	    psio_write_entry(PSIF_CC_INFO, lbl, (char *)tensor[i][0], 27 * sizeof(double));

    //        psio_close(PSIF_CC_LR, 0);
    //        psio_open(PSIF_CC_LR, 0);
       // } else {
       //     outfile->Printf("Using %s tensor found on disk.\n", lbl);
       //     psio_read_entry(PSIF_CC_INFO, lbl, (char *)tensor[i], 27 * sizeof(double));  
       // }
  //}

    for (i = 0; i < params.nomega; i++) free_block(tensor[i]);
    free(tensor);

    free(cartcomp[0]);
    free(cartcomp[1]);
    free(cartcomp[2]);
    free(cartcomp);
}

}  // namespace ccresponse
}  // namespace psi
