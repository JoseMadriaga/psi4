import numpy as np
import psi4
import math
from psi4.driver import constants
psi4.core.clean()
psi4.set_memory('8gb')
# ===== System Setup =====
molA = psi4.geometry("""
 O   1.65152       0.000000     -2.35691
 H   0.646298      0.000000     -3.86023
 H   0.544274      0.000000     -0.926141
 symmetry c1
 no_reorient
 units bohr
 no_com
""")

molB = psi4.geometry("""
 O -0.108659252162 0.0 6.43678512549
 H 1.58963761598 -0.000188972612457 5.73890926769
 H 0.135682335744 0.0 8.25659138345

 symmetry c1
 no_reorient
 units bohr
 no_com
""")

# ===== Calculation Parameters =====
methodname = 'B3LYP'
basisname = 'aug-cc-pvtz'
gembasisname = 'gemherm'
K_exchange = 7.075662

# DPC and regularization parameters
use_dpc = True
regularizer = 0
epsilon = 1e-12

psi4.set_options({
    'puream': False,
    'print': 1,
    'scf_type': 'df'
})

#===== Helper Functions =====
def invert(mat, rhs=None):
    """Matrix inversion with DPC-based thresholding"""
    evals, evecs = np.linalg.eigh(mat)
     
    idx = np.argsort(evals)[::-1]
    evals_sort = evals[idx]
    evecs_sort = evecs[:, idx] 

    if use_dpc and rhs is not None:
        sigma = np.sqrt(np.abs(evals_sort))
        picard_coeffs = np.abs(np.dot(evecs_sort.T, rhs)) / sigma
        print("sigma", sigma)
        print("numer", np.abs(np.dot(evecs_sort.T, rhs)))
        print("picard_coeffs", picard_coeffs)
        epsilon_opt = sigma[-1]
        for i in range(len(sigma)-1):
            #if picard_coeffs[i+1] > picard_coeffs[i]:
            #    epsilon_opt = sigma[i]
            #    break
            if picard_coeffs[i] > sigma[i]:
                index_cutoff = i; 
            
        epsilon_opt = sigma[i]
        psi4.core.print_out(f"\nDPC-selected epsilon: {epsilon_opt**2:.2e}\n")
        evals = np.where(evals < epsilon_opt**2, 0.0, 1.0/evals)
    else:
        evals = np.where(evals < epsilon, 0.0, 1.0/evals)
    
    return np.einsum('ik,k,jk->ij', evecs, evals, evecs)

def invert_svd(mat, rhs=None):
    """
    Pseudoinverse via SVD, with optional DPC-based cutoff selection.

    SVD: mat = U diag(s) V^T
    For symmetric PSD metrics, U ~ V and s ~ eigenvalues, but this works generally.

    DPC rule (your exact logic):
        picard = |U^T rhs| / sigma
        epsilon_opt starts at sigma[-1]
        first i where picard[i+1] > picard[i] -> epsilon_opt = sigma[i]
    where sigma = sqrt(s) (to mirror your eig-based code).

    Notes:
      - If sort_desc=True, we sort singular values s from largest->smallest,
        and reorder U,V consistently. This is recommended for Picard logic.
      - If you *want* to reproduce your "unsorted" behavior, set sort_desc=False.
    """
    mat = np.asarray(mat, dtype=float)

    # Full SVD (thin is fine too, but keep full for clarity)
    U, s, Vt = np.linalg.svd(mat, full_matrices=False)  # s >= 0, usually sorted desc by NumPy

    # Optional: enforce sorting (NumPy already returns sorted descending, but keep explicit)
    #if sort_desc:
    #    idx = np.argsort(s)[::-1]
    #    s = s[idx]
    #    U = U[:, idx]
    #    Vt = Vt[idx, :]

    # Build sigma analogous to your eig-code: sigma = sqrt(|evals|)
    sigma = np.sqrt(s)

    if rhs is not None:
        rhs = np.asarray(rhs, dtype=float)

        # numer = |U^T rhs|
        numer = np.abs(U.T @ rhs)

        # avoid divide by zero safely
        picard_coeffs = numer / np.where(sigma == 0.0, np.inf, sigma)

        print("sigma", sigma)
        print("numer", numer)
        print("picard_coeffs", picard_coeffs)

        epsilon_opt = sigma[-1]
        for i in range(len(sigma) - 1):
            if picard_coeffs[i + 1] > picard_coeffs[i]:
                epsilon_opt = sigma[i]
                break

        # threshold is on s (since s ~ evals for SPD); your print uses epsilon_opt^2

        print(f"\nDPC-selected epsilon: {epsilon_opt**2:.2e}\n")

        # Pseudoinverse of singular values: keep only s >= epsilon_opt^2
        s_inv = np.where(s < epsilon_opt**2, 0.0, 1.0 / s)
    else:
        # fixed threshold on s (analogous to evals < epsilon)
        s_inv = np.where(s < epsilon, 0.0, 1.0 / s)

    # mat^+ = V diag(s_inv) U^T
    inv = (Vt.T * s_inv) @ U.T
    return inv

def build_gem_field(wfn, add_coulomb=True, add_exchange=False):
    """Builds GEM external potential field without using factorial2"""
    psi4.core.print_out("\nComputing GEM representation\n")
    
    psi4.core.print_out("\norbital_basis\n")
    orbital_basis = wfn.basisset()
    psi4.core.print_out("\nmolecule\n")
    molecule = orbital_basis.molecule()
    psi4.core.print_out("\naux_basis\n")
    aux_basis = psi4.core.BasisSet.build(molecule, "ORBITAL", gembasisname)
    psi4.core.print_out("\napply_hermite_normalization\n")
    aux_basis.apply_hermite_normalization()
    psi4.core.print_out("\nzero_basis\n")
    zero_basis = psi4.core.BasisSet.zero_ao_basis_set()
    psi4.core.print_out("\nmints\n")
    mints = psi4.core.MintsHelper(orbital_basis)
    psi4.core.print_out("\nfield\n")
    field = psi4.QMMMbohr().extern

    # Coulomb terms
    J_PQ = np.squeeze(mints.ao_eri(aux_basis, zero_basis, aux_basis, zero_basis))
    J_Pmn = np.squeeze(mints.ao_eri(aux_basis, zero_basis, orbital_basis, orbital_basis))
    
    d = 2 * np.einsum('Qmn,mn->Q', J_Pmn, wfn.Da())
    J_PQinv = invert(J_PQ + regularizer*np.eye(aux_basis.nbf()), d if use_dpc else None)
    fit_coefficients = np.einsum('PQ,Q->P', J_PQinv, d)
    print(fit_coefficients)
    # Charge constraint
    nshell = aux_basis.nshell()
    q = []
    for sh in range(nshell):
        shell = aux_basis.shell(sh)
        L = shell.am
        ex = shell.exp(0)
        coef = shell.coef(0)
        prefac = coef * np.sqrt((np.pi**3)/(2**L * ex**(L+3)))
        
        fact_terms = [1]
        for n in range(1, L+1):
            fact_terms.append(fact_terms[-1] * n if n%2 == 1 else fact_terms[-2] * n)
        
        for lx in range(L, -1, -1):
            for lz in range(L - lx + 1):
                ly = L - lx - lz
                if (lx%2 == 0) and (ly%2 == 0) and (lz%2 == 0):
                    term = 1
                    for val in [lx-1, ly-1, lz-1]:
                        if val >= 0:
                            term *= fact_terms[val]
                    q.append(prefac * term)
                else:
                    q.append(0.0)

    Q = sum(molecule.Z(i) for i in range(molecule.natom())) - molecule.molecular_charge()
    q = np.array(q)
    qdot = np.dot(J_PQinv, q)
    lam = (Q - np.dot(qdot, d)) / np.dot(qdot, q)
    fit_coefficients += lam * qdot

    # Add potential to field - CORRECTED VERSION
    coefs_vec = psi4.core.Vector.from_array(fit_coefficients)
    if add_coulomb:
        field.addBasis(aux_basis, coefs_vec)
    
    if add_exchange:
        exchange_coefs = psi4.core.Vector(aux_basis.nbf())
        for i in range(aux_basis.nbf()):
            exchange_coefs.set(i, fit_coefficients[i] * K_exchange)
        field.addExchangeBasis(aux_basis, exchange_coefs)
    
    for n in range(molecule.natom()):
        field.addCharge(molecule.Z(n), molecule.x(n), molecule.y(n), molecule.z(n))
    
    return field
# ===== Main Calculation =====
psi4.core.print_out("\n===== Computing Monomer A =====\n")
eA, wfnA = psi4.energy(methodname+'/'+basisname, molecule=molA, return_wfn=True)

psi4.core.print_out("\n===== Computing Monomer B =====\n")
eB, wfnB = psi4.energy(methodname+'/'+basisname, molecule=molB, return_wfn=True)

# Coulomb perturbation
psi4.core.print_out("\n===== Computing Coulomb Perturbation =====\n")
field = build_gem_field(wfnA, add_coulomb=True, add_exchange=False)
psi4.core.set_global_option_python('EXTERN', field)
eB_Jperturbed = psi4.energy(methodname+'/'+basisname, molecule=molB)

# Exchange perturbation
psi4.core.print_out("\n===== Computing Exchange Perturbation =====\n")
field = build_gem_field(wfnA, add_coulomb=True, add_exchange=True)
psi4.core.set_global_option_python('EXTERN', field)
eB_JKperturbed = psi4.energy(methodname+'/'+basisname, molecule=molB)

# ===== Results =====
psi4.core.print_out("\n===== GEM Results =====\n")
psi4.core.print_out(f"Coulomb Energy:     {constants.hartree2kcalmol*(eB_Jperturbed-eB):.4f} kcal/mol\n")
psi4.core.print_out(f"Exchange Energy:    {constants.hartree2kcalmol*(eB_JKperturbed-eB_Jperturbed):.4f} kcal/mol\n")
