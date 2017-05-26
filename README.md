# monomer_dimer_tilings

Code for paper "Power series expansions for the planar monomer-dimer problem"

The code consists of two Sage scripts and one C program.
It was tested with Sage 6.8 and gcc 4.9.2.
If you want to launch the code, you should compile C program with make.

**Compute the free energy modulo prime**
Launch 
sage tilings.sage p n
where p is a prime number not exceeding 2^(31) - 1, and n is the number of correct terms to be computed.

**Reconstruct from values modulo prime**
Outputs of tilings.sage should have names of the form Mask1, Mask2, ..., MaskN.
If you want to reconstruct first n terms of the result, launch
sage collect.sage Mask n
