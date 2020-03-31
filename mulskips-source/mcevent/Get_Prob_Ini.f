!> @file
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
**************************************************************************8

      SUBROUTINE Get_Prob_Ini(Occ,NSiC,Coor,Site,NextN,  ! Only for pure 3C Configurations
     >                    Index_Event,Prob_Event) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Occ,NdiC,Coor,CoorNN
      INTEGER Site(3),NextN(3,4)
      INTEGER i,Index_Event,SiOrC,NSiC,N_C,N_Si
      REAL(8) Prob_Event
      REAL(8):: Ptot,Pcurr
      REAL(8) :: random
      EXTERNAL random

      IF (Occ.EQ.0)THEN ! Depositions
         IF (Coor.EQ.1)THEN
            CoorNN=IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosCoor,LenCoor)
            IF(CoorNN.GE.2)THEN ! in the case Coor=1 the single first neighbor should have CoorNN=3
              Ptot=PtransD(1)+PtransD(2)
              Pcurr=random(idum)*Ptot
              IF(Pcurr.LE.PtransD(1))THEN
                Index_Event = 1
              ELSE
                Index_Event = 2
              END IF
              Prob_Event = Ptot
            ELSE ! Event to be excluded but MC particle is present in the list
              Index_Event = 1
              Prob_Event = 1.e-6
            END IF
         ELSE IF (Coor.EQ.2)THEN
            Ptot=PtransD(3)+PtransD(4)
            Pcurr=random(idum)*Ptot
            IF(Pcurr.LE.PtransD(3))THEN
               Index_Event = 1
            ELSE
               Index_Event = 2
            END IF
            Prob_Event = Ptot
         ELSE IF (Coor.EQ.3)THEN
            Ptot=PtransD(5)+PtransD(6)
            Pcurr=random(idum)*Ptot
            IF(Pcurr.LE.PtransD(5))THEN
               Index_Event = 1
            ELSE
               Index_Event = 2
            END IF
            Prob_Event = Ptot
         ELSE IF(Coor.EQ.4)THEN! Vacancy Site?  Error
            Index_Event = 0
            Prob_Event = 0.d0
            write(*,*)Coor,'Get Prob Vacancy Generation'
         END IF
      ELSE ! Occ=1 Evaporation
         N_Si=0
         N_C =0
         !write(*,*)'Get Prob Site',Site,Occ,Coor
         DO i=1,4
           SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
           IF(SiOrC.EQ.1)N_Si=N_Si+1
           IF(SiOrC.EQ.2)N_C=N_C+1
           !write(*,*)'Get Prob SiOrC',NextN(1:3,i),SiOrC
         END DO
         IF(N_Si+N_C.NE.Coor)THEN
            write(*,*)Coor,N_Si,N_C,'Get Prob Error N_Si+N_C.NE.Coor'
            STOP
         END IF
         IF (N_Si+N_C.GT.3)THEN
            Index_Event = 0
            Prob_Event = 0.d0
            !write(*,*)Coor,'Get Prob Bulk Site'
         ELSE IF(N_Si+N_C.EQ.0)THEN
            Index_Event = 3
            Prob_Event = 100000.*Tree(1)  ! it makes its evaporation a fast event
            write(*,*)Coor,'Get Prob Isolated Atom'
         ELSE
            Index_Event = 3
            Prob_Event = PTransE(NSic,N_Si,N_C)
         END IF
      END IF
      END SUBROUTINE Get_Prob_Ini
