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
**************************************************************************************
**
**
**************************************************************************************
      PROGRAM KLMCSiC3C
      USE Definitions
      USE DefSystem
      IMPLICIT NONE
      INTEGER :: i,j,k,l,in,irot,K1,Itrans,Iter=1
      INTEGER(8) :: MaxIter=2100000000
      INTEGER :: Ind,IndTest
      INTEGER :: MGigaCycle=400,IGigaCycle=0,Kilo=1000,Mega=1000000
      INTEGER :: exit_zeta=10000
      INTEGER :: Len=LenZ-1,Len1=120,Len2=120,Len3=120
      INTEGER :: OutMolMol=100000000,IterOut=0,MolmOutpSequ=0
! Random NumBer
      REAL(8) :: Pi=3.14159265358979323
      REAL(8)  :: Time=0.,Dt,Patom
      INTEGER :: Site(3),Coord,Nbon,Njum,Jump,ijump
      INTEGER :: Site1(3),Site2(3),SiteC(3),Newsite(3,2)
      INTEGER :: SNSite(3,3),SNZigZag(3,3),SNArmChair(3,3)
      LOGICAL   MolmOutpStat
      REAL(8) Random
      EXTERNAL Random
! Executable statments
8     FORMAT(I4,I5,' ',8ES10.3)
9     FORMAT(I4,' ',8ES10.3)
10    FORMAT('Jump: ',3I3,'   from: ',3I4)
11    FORMAT(I7,I5,I3,' ',ES10.3)
12    FORMAT(I8,ES10.3,' ',2F8.3)
13    FORMAT('o_site:',3I4,'  bon:',I2,'  jum: ',I2, ' nn:',8I2)
14    FORMAT('n_site:',3I4,'  bon:',I2,'  jum: ',I2, ' nn:',8I2)
115   FORMAT('JumpSt:',12I9)

      LattCoo=0
      LattInd=0
      IDUM = -9117116
      WRITE(*,*)'IDUM',IDUM
      CALL Fileopen()
      CALL SetProbability()
      CALL AllocateArrays()
      NumAtoms=0
      NumVoids=0
      NumAdAtom=0
      READ(IPF,*)InitSt
      READ(IPF,*)Len1,Len2,Len3
      READ(IPF,*)OutMolMol
      READ(IPF,*)MaxIter
      READ(IPF,*)exit_zeta
     
      IF(InitSt.EQ.'F')THEN
         Len1=Len1-MOD(Len1,24)
         CALL SetSiC3C(Len1)
      ELSE IF(InitSt.EQ.'A')THEN
         CALL SetAFBSymZimb1(Len1)
      ELSE IF(InitSt.EQ.'I')THEN
         CALL SetInvPOff(Len1,Len2)
      ELSE IF(InitSt.EQ.'D')THEN
         CALL SetInvPC(Len1,Len2)
      ELSE IF(InitSt.EQ.'Z')THEN
         CALL SetInvPSi(Len1,Len2)
      ELSE IF(InitSt.EQ.'J')THEN
         CALL SetInvPwAPB(Len1,Len2)
      ELSE IF(InitSt.EQ.'K')THEN
         CALL SetInvPwOnlySiorC(Len1,Len2)
      ELSE IF(InitSt.EQ.'L')THEN
         CALL SetInvPwOnlyCorSi(Len1,Len2)         
      ELSE IF(InitSt.EQ.'S')THEN
         CALL SetNCSp(Len1)
      ELSE IF(InitSt.EQ.'C')THEN
         CALL SetNC(Len1,Len2,Len3)
      ELSE IF(InitSt.EQ.'T')THEN
         CALL SetTrench(Len1,Len2,Len3)
      ELSE
         write(*,*)'initialization not_implemented'
         STOP
      END IF
      write(*,*)NumAtoms,NumVoids,NumAdAtom
      CALL WriteMolMolSource(Time, IterOut,1)
      DO WHILE(Iter.LE.MaxIter)
       !IterOut=Iter/10
       IterOut=Iter
       IF(NumAdAtom.LE.20)THEN
         write(*,*)"few MC particles: stop MC",Iter,NumAdAtom
         EXIT
       END IF
       CALL PickTreeEvent(Ind,Patom)
       Time=Time+1.0/Tree(1)
       Itrans=ListAdAtom(Ind) % Ind_Event
       IF(MOD(Iter,OutMolMol).EQ.0)THEN
          CALL WriteMolMolSource(Time, IterOut,OutMolMol)
       END IF
       !IF(MOD(Iter,10000000).EQ.0)THEN
       IF(MOD(Iter,OutMolMol).EQ.0)THEN
          Site=ListAdAtom(Ind) % AtomXYZ
          write(*,*)'Iter',Iter,Time,Site
          write(*,*)'NatNvoidNad',NumAtoms,NumVoids,NumAdAtom
       END IF

!       Site=ListAdAtom(Ind) % AtomXYZ
!       Coord=IBITS(LattCoo(Site(1),Site(2),Site(3)),
!     >                  PosCoor,LenCoor)
       IF(Itrans.EQ.3)THEN
!         write(*,*)"Evap",Iter,Ind,Coord
!         write(*,*)ListAdAtom(Ind) % AtomXYZ
!         write(*,*)ListAdAtom(Ind) % NextNXYZ
         CALL Evaporation(Ind)
!         Coord=IBITS(LattCoo(Site(1),Site(2),Site(3)),
!     >                  PosCoor,LenCoor)
!         write(*,*)Coord
!         write(*,*)ListAdAtom(Ind) % NextNXYZ
       ELSE
!         write(*,*)"Depo",Iter,Ind,Coord
!         write(*,*)ListAdAtom(Ind) % AtomXYZ
!         write(*,*)ListAdAtom(Ind) % NextNXYZ
         CALL Deposition(Ind,Itrans)
!         Coord=IBITS(LattCoo(Site(1),Site(2),Site(3)),
!     >                  PosCoor,LenCoor)
!         write(*,*)Coord
!         write(*,*)ListAdAtom(Ind) % NextNXYZ
       END IF
       IF (Site(3).GT.exit_zeta) exit
       Iter = Iter +1
      END DO
      CALL WriteMolMolSource(Time, IterOut,OutMolMol)
      CALL Fileclose()

      END PROGRAM KLMCSiC3C
