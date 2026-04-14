This README file is to use DPC to obtain the coefficients using the GEM implementation of Psi4. To obtain proper coefficients, please go through the following steps:
-Obtain your optimized water dimer. In the first part of the code, where the co-ordinates have to be entered as the example provided below:

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

Here, molA is the molecule whose density will be calculated, using which the coefficients will be obtained. !!MAKE SURE THE CO-ORDIANTES ARE IN BOHR AND NOT IN ANGSTROMS(COMMON UNIT IN MOST FILES)!! 1 angstorm = 1.88973 bohr.

-Select appropriate level of theory and the value of kexch like the example here:

# ===== Calculation Parameters =====
methodname = 'B3LYP'
basisname = 'aug-cc-pvtz'
gembasisname = 'gemherm'
K_exchange = 7.075662

!!The "gembasisname" file contains information regarding the auxiliary basis set to be used. You can find all the available ABSs in the "input.aux" file. 

-Change the "use_dpc" flag to True if you want to use DPC. If not, please make sure to input appropriate truncation value for "epsilon".

-Execute the file using "python psi4_input_with_dpc.py > output.dat". This will ensure that the "fit_coefficients" are printed to the file. (They will be present at the very end of the file.)

-Put the coefficients in the GEM format in a new file: 
#Number of coefficients
coeff#1
coeff#2
.
.
.
.
.
coeff#x

-Use the "transform.py" file to transform these cofficients into gem compatible coefficients
To use this, three files are needed: an input.prm, input.aux and the file with the coefficients in the GEM format.
An example of the input.prm for the A4 basis set is present in the directory. Currently, the transformation only works for the A4 ABS.

-Use these coefficients in GEM_site_site to calculate inter-molecular interactions
