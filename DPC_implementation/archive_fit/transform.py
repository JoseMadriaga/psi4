import numpy as np
import libinv
import os
import sys

coord = np.zeros((3,3))
with open("water.xyz", "r") as f:
     line = f.readline()
     parts = line.split()
     coord[0,0] = parts[1]
     coord[0,1] = parts[2]
     coord[0,2] = parts[3]

     line = f.readline()
     parts = line.split()
     coord[1,0] = parts[1]
     coord[1,1] = parts[2]
     coord[1,2] = parts[3]

     line = f.readline()
     parts = line.split()
     coord[2,0] = parts[1]
     coord[2,1] = parts[2]
     coord[2,2] = parts[3]

flabel = ["O","H1","H2"]

coord = coord/0.52917720859

par = open("inputA4.prm", "r")
aux = open("input.aux", "r")

#for i in range(3):
#    print(flabel[i], sys.argv[1], coord[i,0], coord[i,1], coord[i,2])

libinv.AHBASE_load_auxbasis(par,aux)
libinv.AHBASE_calcnorms()
libinv.AHTYPE_load_site_info(par,aux)
libinv.AHCOEFF_copy_norms_from_basis()
libinv.AHFRAME_get_molframe(par,flabel,coord)
libinv.AHFRAME_load_deflist()
libinv.AHFRAME_build_frames(coord)

#localcoeffs =  sys.argv[1]
#local_herm_coeffs=np.loadtxt(localcoeffs)

origin = sys.argv[1]
aux_coefs = np.loadtxt(origin)
aux_coefs = libinv.reordenar_inverso(aux_coefs)
global_herm_coeffs = libinv.cart_to_gherm(aux_coefs)
local_herm_coeffs = libinv.AHFRAME_global_to_local(global_herm_coeffs)

os.system('rm -rf __pycache__')
