# N3p_PESs
Potential energy surfaces for the ground state of N3+

**Requirements**

(1) RKHS toolkit : Download from https://github.com/MeuwlyGroup/RKHS

**Compile a particular PES**

A test program file (pes_test.f90) is given and it can be compiled

`gfortran RKHS.f90 N3+_PES.f90 pes_test.f90`

**Running the executable**

Before running the executable make sure that the asymp.dat and pes1.kernel files for that PES present in the current directory (or change the file path in the fortran program).

**Cite as**
 Debasish Koner, Max Schwilk, Sarbani Patra, Evan J. Bieske and Markus Meuwly,  	arXiv:2004.12404 [physics.chem-ph]
