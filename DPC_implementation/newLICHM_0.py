import numpy as np
import sys
import math
from psi4.driver import constants
psi4.set_memory('24 gb')
molA = psi4.geometry("""
 0 1
 O -1.4721136585713 -0.0000963760323 2.13908682505816
 H -1.2604340970020 1.43032236532654 3.22583052477304
 H -1.2793408068782 -1.5039271978783 3.12562212783960
 symmetry c1
 no_reorient
 units bohr
 no_com
""")

molB = psi4.geometry("""
 0 1
 O 1.65149969039789 0.00000000000000 -2.3565583971992
 H 0.64636192364623 0.00000000000000 -3.8604289032938
 H 0.54428647730172 0.00000000000000 -0.9261774503628
 symmetry c1
 no_reorient
 units bohr
 no_com
""")

methodname = 'wB97X-D'
basisname = 'aug-cc-pVDZ'
gembasisname = 'gemherm_old'
K_exchange = 7.07566

epsilon = 1e-12
regularizer = 1e-8
constrain = True

psi4.set_options({'puream' : False, 'print' : 1,'scf_type' : 'pk'})

eA, wfnA = psi4.energy(methodname+'/'+basisname, molecule=molA, return_wfn=True)
eB, wfnB = psi4.energy(methodname+'/'+basisname, molecule=molB, return_wfn=True)
print("eA", eA * constants.hartree2kcalmol)

def invert(mat):
    if epsilon == 0:
        c = np.linalg.inv(np.linalg.cholesky(mat))
        return np.dot(c.T,c)
    else:
        evals, evecs = np.linalg.eigh(mat)
        evals = np.where(evals < epsilon, 0.0, 1.0/evals)
        return np.einsum('ik,k,jk->ij', evecs, evals, evecs)

def build_gem_field(wfn, add_coulomb, add_exchange, fit):
    print("building gem field")
    orbital_basis = wfn.basisset()
    molecule = orbital_basis.molecule()
    aux_basis = psi4.core.BasisSet.build(molecule, 'ORBITAL', gembasisname)
    aux_basis.apply_hermite_normalization()
    zero_basis = psi4.core.BasisSet.zero_ao_basis_set()
    mints = psi4.core.MintsHelper(orbital_basis)
    field = psi4.QMMMbohr().extern
    print("going into the fit")
    if fit:
        J_PQ = np.squeeze(mints.ao_eri(aux_basis, zero_basis, aux_basis, zero_basis))
        J_PQinv = invert(J_PQ + regularizer*np.eye(aux_basis.nbf()))
        J_Pmn = np.squeeze(mints.ao_eri(aux_basis, zero_basis, orbital_basis, orbital_basis))
        d =  2 * np.einsum('Qmn,mn->Q', J_Pmn, wfn.Da())
        fit_coefficients = np.einsum('PQ,Q->P', J_PQinv, d)
        print("gone through the inverse metric")
        if constrain:
            def dfact(n):
                if n <= 1:
                    return 1.0
                if n%2 == 1:
                    k = (n+1)//2
                    return math.factorial(2*k) / (2**k * math.factorial(k))
                else:
                    k = n//2
                    return 2**k * math.factorial(k)
            nshell = aux_basis.nshell()
            q = []
            print("after the constraint")
            for sh in range(nshell):
                shell = aux_basis.shell(sh)
                nprim = shell.nprimitive
                if nprim > 1:
                    raise('This code currently assumes uncontracted GEM basis sets')
                L = shell.am
                ex = shell.exp(0)
                coef = shell.coef(0)
                prefac = coef * np.sqrt((np.pi**3)/(2**L * ex**(L+3)))
                for lx in range(L,-1,-1):
                    for lz in range(L-lx+1):
                        ly = L - lx - lz
                        if (lx%2==0) and (ly%2==0) and (lz%2==0):
                            q.append(prefac*dfact(lx-1)*dfact(ly-1)*dfact(lz-1))
                        else:
                            q.append(0.0)
            natom = molecule.natom()
            Q = 0
            print("adding in charge")
            for i in range(natom):
                Q += molecule.Z(i)
            Q -= molecule.molecular_charge()
            old_nelec = np.dot(fit_coefficients, q)
            psi4.core.print_out(f'Number of electrons before constraining: {old_nelec:.6f}')
            qdot = np.dot(J_PQinv, q)
            lam = (Q - np.dot(qdot,  d)) / np.dot(qdot, q)
            fit_coefficients += lam * qdot
            new_nelec = np.dot(fit_coefficients, q)
            psi4.core.print_out(f'\nNumber of electrons after constraining: {new_nelec:.6f}')
            np.savetxt('fitcoeffs_test.dat',fit_coefficients)
            print("should have printed fitcoefs")
    else:
        fit_coefficients = np.loadtxt('fitcoeffs.dat')
    coefs_vec = psi4.core.Vector(aux_basis.nbf())
    for n, coef in enumerate(fit_coefficients): coefs_vec.set(n, coef)
    field.addBasis(aux_basis, coefs_vec)
    if add_exchange:
        K_coefs_vec = psi4.core.Vector(aux_basis.nbf())
        for n, coef in enumerate(fit_coefficients): K_coefs_vec.set(n, K_exchange*coef)
        field.addExchangeBasis(aux_basis, K_coefs_vec)
    for n in range(molecule.natom()):
        field.addCharge(molecule.Z(n), molecule.x(n), molecule.y(n), molecule.z(n))
    print("field", field)
    return field

field = build_gem_field(wfnB, add_coulomb=True, add_exchange=False, fit=True)
psi4.core.set_global_option_python('EXTERN', field)
eA_Jperturbed = psi4.energy(methodname+'/'+basisname, molecule=molA)
print("EA_JPerturbed", eA_Jperturbed*constants.hartree2kcalmol)

field = build_gem_field(wfnB, add_coulomb=True, add_exchange=True, fit=True)
psi4.core.set_global_option_python('EXTERN', field)
eA_JKperturbed = psi4.energy(methodname+'/'+basisname, molecule=molA)
print("EA_JKperturbed", eA_JKperturbed* constants.hartree2kcalmol)

psi4.core.print_out('\n Coulomb Energy: {:6} kcal/mol'.format(constants.hartree2kcalmol*(eA_Jperturbed-eA)))
psi4.core.print_out('\n Exchange Energy: {:6} kcal/mol\n'.format(constants.hartree2kcalmol*(eA_JKperturbed-eA_Jperturbed)))
psi4.core.print_out('\n QM+GEM Energy: {:6} \n'.format(eA_JKperturbed))
