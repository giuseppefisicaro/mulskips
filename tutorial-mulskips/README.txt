!> @Tutorial
!!
!!
!!   Copyright (C) 2019-2020
!!   @authors: Ioannis Deretzis, Giuseppe Fisicaro and Antonino La Magna
!!   This file is part of the MulSKIPS code.
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/licenses/gpl-3.0.txt .
!!   For the list of contributors, see ~/AUTHORS

!!   MulSKIPS is a free software: you can redistribute it and/or modify
!!   it under the terms of the GNU General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.

!!   MulSKIPS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details

This is a tutorial for the MulSKIPS code.

With MulSKIPS you can simulate the epitaxial growth of cubic silicon carbide (3C-SiC)
starting from various input structures.
See references [1,2] for details about the code and its application.

Each folder of the tutorial, that are "sphere, cube, surface, and inverted-pyramid",
contains a start.dat file which has been generated with the jupyter notebook "run-mulskips.ipynb".
To run the code you just need an empty folder containing the input file "start.dat".

The "run-mulskips.ipynb" notebook allows you to set the input file "start.dat", to compile and run the mulskips
Superlattice Kinetic Monte Carlo (KMC) code.

In particular the notebook allows you to (in order):
- generate the input file `start.dat` of mulskips code;
- make a copy of the code in the folder `mulskips-source` and compile it;
- make a folder where mulskips runs and all output files are written;
- analize the output.

In case you want to generate your own start.dat input file, you can easily follow
all steps in "run-mulskips.ipynb". For example, to modify the process temperature
or the deposition frequencies.

This is an example of the start.dat file:

######################## start.dat #########################################################################
1)    0.000458642486241 ! PtransE_Si[0,1]     evaporation frequency for one Si bound to one C
2)    0.00107413320314 ! PtransE_Si[1,0]      evaporation frequency for one Si bound to one Si
3)    0.000654677677311 ! PtransE_Si[1,1]     evaporation frequency for one Si bound to one Si and one C
4)    9.03452031205e-05 ! PtransE_Si[0,2]     evaporation frequency for one Si bound to two C
5)    0.000717033435226 ! PtransE_Si[2,0]     evaporation frequency for one Si bound to two Si
6)    0.000142385075632 ! PtransE_Si[2,1]     evaporation frequency for one Si bound to two Si and one C
7)    2.91989021579e-05 ! PtransE_Si[1,2]     evaporation frequency for one Si bound to one Si and two C
8)    3.18575387165e-07 ! PtransE_Si[0,3]     evaporation frequency for one Si bound to three C
9)    0.000518948689858 ! PtransE_Si[3,0]     evaporation frequency for one Si bound to three Si
10)   10693067.5485 ! PtransE_C[0,1]          evaporation frequency for one C bound to one C
11)   336.433275007 ! PtransE_C[1,0]          evaporation frequency for one C bound to one Si
12)   9.14549727245 ! PtransE_C[1,1]          evaporation frequency for one C bound to one Si and one C
13)   327556.955526 ! PtransE_C[0,2]          evaporation frequency for one C bound to two C
14)   4.99391159317e-06 ! PtransE_C[2,0]      evaporation frequency for one C bound to two Si
15)   2.81716126678e-06 ! PtransE_C[2,1]      evaporation frequency for one C bound to two Si and one C
16)   0.288419678521 ! PtransE_C[1,2]         evaporation frequency for one C bound to one Si and two C
17)   20148.1158068 ! PtransE_C[0,3]          evaporation frequency for one C bound to three C
18)   2.85295646127e-07 ! PtransE_C[3,0]      evaporation frequency for one C bound to three Si
19)   9.98730082416e-05 ! PtransD[1]          deposition frequency a one coordinated Si
20)   0.000998782318633 ! PtransD[2]          deposition frequency a one coordinated C
21)   2.49682520604e-05 ! PtransD[3]          deposition frequency a two coordinated Si
22)   1.24847789829e-05 ! PtransD[4]          deposition frequency a two coordinated C
23)   0.00399492032966 ! PtransD[5]           deposition frequency a three coordinated Si
24)   0.00199756463727 ! PtransD[6]           deposition frequency a three coordinated C
25)   1.0 ! PtransZig
26)   C ! Initstat: S Sphere, C Parellelepipid, F Flat (100) surface, A APB, I inverted pyramid, D inverted pyramid of C, Z inverted pyramid of Si, J inverted pyramid with APB
27)   120 120 120 ! Len1 Len2 Len3
28)   100000  ! OutMolMol -> output frequency
29)   5000000 ! IterMax -> Max number of iterations
30)   10000   ! exit_zeta strategy after lenz=500
######################## end start.dat #########################################################################

Lines 1-9 refer to the evaporation frequencies for one/two/three coordinated Si atoms.
Lines 10-18 refer to the evaporation frequencies for one/two/three coordinated C atoms.
Lines 19-24 refer to the deposition frequencies for one/two/three coordinated Si and C atoms.
Line 25 sets the probability for stacking fault formation (1, no stacking fault with growth of the perfect cubic SiV. See Ref. [1] for more details).
Line 26 sets the input structure.
Line 27 sets sizes of the input structure.
Line 28 sets output frequency.
Line 29 sets the max number of KMC iterations.
Line 30 exit_zeta strategy after lenz=500.

You can directily compile the mulskips code and runs it on each folder.
To compile mulskips, touch the makefile specifying the "DESTDIR" and the "FC" Fortran compiler.
Then type "make" from command line.

Once mulskips is compiled you can run it on each folder "sphere, cube, surface, and inverted-pyramid"
which contains the file "start.dat".
Type from command line the command:

<path-of-the-compiled-code>/mulskips.e

Folders refer to:
sphere            -> epitaxial growth starting from a sphere
cube              -> epitaxial growth starting from a cube
surface           -> epitaxial growth starting from a flat (001) surface
inverted-pyramid  -> epitaxial growth starting from a patterned substrate with an inverted pyramid

In these input files start.dat we just modified line 26.

After the mulskips runs, several output files are generated.
The files of interest are the *.xyz files which store atomic coordinates of KMC particles.
To visualize the evolving KMC structures, like in a molecular dynamics simulation,
you have to open with any software for atomistic visualization (VESTA, VMD, V_sim) the *.xyz files:

VESTA I000000*.xyz   -> These files store all undercoordinated atoms (atoms with coordination lower than four)
VESTA I000000*_d.xyz -> These files store atoms which are not epitaxially ordered
VESTA I000000*_v.xyz -> These files store all vacancy point defects generated during the epitaxial growth

References:
[1] A. La Magna, A. Alberti, E. Barbagiovanni, C. Bongiorno, M. Cascio, I. Deretzis, F. La Via, and E. Smecca, "Simulation of the growth kinetics in group iv compound semiconductors" Phys. Status Solidi A 216, 1800597 (2019).
[2] G. Fisicaro, C. Bongiorno, I. Deretzis, F. Giannazzo, F. La Via, F. Roccaforte, M. Zielinski, M. Zimbone, and A. La Magna, "Genesis and evolution of extended defects: The role of evolving interface instabilities in cubic SiC" Applied Physics Reviews, 7(2), 021402 (2020).
