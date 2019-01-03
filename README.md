# x_band

## Project Description
This code is written to extract the orbital projections from a VASP PROCAR file containing s, p, d and f projections for specific atoms and allows also to extract the total spin orientations (mx, my mz).  
Furthermore the states can be ordered by the eigenenergies, which turns out to be useful for hybrid functional computations, here the states are sometimes not energy-ordered in the PROCAR.  
As this small project is still in its infancy, the code is not well structured yet and can only read the PROCARs from VASP.5.4.1 and VASP.5.4.4 as the format of the PROCAR is version-dependent.

## Compilation
The repository contains a Makefile for ifortran.

## Required input files
x_input.dat - Input options  
PROCAR  
CONTCAR - Needed to compute reciprocal lattice vectors  
