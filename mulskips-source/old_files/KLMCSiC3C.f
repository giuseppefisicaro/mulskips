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

**************************************************************************************
**    This subroutine find remanent two next-neighbour site of SiteC when
**    two next-neighbour are given
**************************************************************************************
      SUBROUTINE FIND_NN(Site1,Site2,SiteC,NewSite,LX,LY)
      IMPLICIT NONE
      REAL(8) :: Pi=3.14159265358979323
      REAL(8)  :: Cth,Sth,snorm,Rot(3,3),ax(3)
      REAL(8)  :: RNewsite(3)
      INTEGER :: i,k,l,irot,LX,LY
      INTEGER :: Site1(3),Site2(3),SiteC(3),Newsite(3,2)
      INTEGER :: Site11(3),Site22(3)
      Site11=Site1
      Site22=Site2
! rotation axis as versor of the direction SiteC-Site1
      IF(IABS(SiteC(1)-Site11(1)).GT.LX/2)
     >   Site11(1)=Site11(1)+ISIGN(LX,SiteC(1)-Site11(1))
      IF(IABS(SiteC(2)-Site11(2)).GT.LY/2)
     >   Site11(2)=Site11(2)+ISIGN(LY,SiteC(2)-Site11(2))

      IF(IABS(SiteC(1)-Site22(1)).GT.LX/2)
     >   Site22(1)=Site22(1)+ISIGN(LX,SiteC(1)-Site22(1))
      IF(IABS(SiteC(2)-Site22(2)).GT.LY/2)
     >   Site22(2)=Site22(2)+ISIGN(LY,SiteC(2)-Site22(2))

      ax=DFLOAT(SiteC)-DFLOAT(Site11)
      snorm=DSQRT(ax(1)**2+ax(2)**2+ax(3)**2)
      ax=ax/snorm
      DO irot=1,2
! rotation matrix for an angle 2/3 Pi or 4/3 Pi
         Cth=DCOS(DFLOAT(2*irot)*pi/3.)
         Sth=DSIN(DFLOAT(2*irot)*pi/3.)
         DO i=1,3
            Rot(i,i)=ax(i)**2*(1-Cth)+Cth
         END DO
         Rot(2,1)=ax(1)*ax(2)*(1-Cth)-ax(3)*Sth
         Rot(1,2)=ax(1)*ax(2)*(1-Cth)+ax(3)*Sth
         Rot(3,1)=ax(1)*ax(3)*(1-Cth)+ax(2)*Sth
         Rot(1,3)=ax(1)*ax(3)*(1-Cth)-ax(2)*Sth
         Rot(3,2)=ax(2)*ax(3)*(1-Cth)-ax(1)*Sth
         Rot(2,3)=ax(2)*ax(3)*(1-Cth)+ax(1)*Sth
! rotate Site2
         DO k=1,3
            RNewSite(k)=0.0
            DO l=1,3
               RNewSite(k)=RNewSite(k)+Rot(l,k)*
     >                   (DFLOAT(Site22(l))-DFLOAT(SiteC(l)))
            END DO
         END DO
         NewSite(1:3,irot)=NINT(RNewSite)+SiteC
         IF(NewSite(1,irot).GE.LX)NewSite(1,irot)=NewSite(1,irot)-LX
         IF(NewSite(1,irot).LT.0)NewSite(1,irot)=NewSite(1,irot)+LX
         IF(NewSite(2,irot).GE.LY)NewSite(2,irot)=NewSite(2,irot)-LY
         IF(NewSite(2,irot).LT.0)NewSite(2,irot)=NewSite(2,irot)+LY
!         write(*,*)'FindNN',NewSite(1:3,irot),RNewSite
      END DO
      END SUBROUTINE FIND_NN


**************************************************************************************
**    This subroutine find Arm-Chair and Zig-Zag configurations which next-neighbour
**    of the site of Site and Second Neighbour of the Site SiteC. Site C has a
**    fixed Coordination
**************************************************************************************
      SUBROUTINE FIND_SN(Site,SiteC_inp,SNSite_inp,SNZiGZag,
     >                   SNArmChair,LX,LY)
      IMPLICIT NONE
      REAL(8)  :: snorm,ax(3),Cent(3)
      REAL(8)  :: RSite(3),RSiteC(3)
      REAL(8)  :: RSNZiGZag(3),RSNArmChair(3)
      INTEGER :: i,LX,LY,Site(3),SiteC_inp(3)
      INTEGER :: SNSite_inp(3,3),SNZiGZag(3,3),SNArmChair(3,3)
      INTEGER :: SNSite(3,3),SiteC(3)
      LOGICAL debug_bc

      SNSite=SNSite_inp
      SiteC=SiteC_inp
!      write(*,*)'standard',Site,SiteC,SNSite
      debug_bc=.FALSE.
      IF(IABS(Site(1)-SiteC(1)).GT.LX/2)THEN
         SiteC(1)=SiteC(1)+ISIGN(LX,Site(1)-SiteC(1))
         debug_bc=.TRUE.
      END IF
      IF(IABS(Site(2)-SiteC(2)).GT.LY/2)THEN
         SiteC(2)=SiteC(2)+ISIGN(LY,Site(2)-SiteC(2))
         debug_bc=.TRUE.
      END IF
      DO i=1,3
        IF(IABS(Site(1)-SNSite(1,i)).GT.LX/2)THEN
           SNSite(1,i)=SNSite(1,i)+ISIGN(LX,Site(1)-SNSite(1,i))
           debug_bc=.TRUE.
        END IF
        IF(IABS(Site(2)-SNSite(2,i)).GT.LY/2)THEN
           SNSite(2,i)=SNSite(2,i)+ISIGN(LY,Site(2)-SNSite(2,i))
           debug_bc=.TRUE.
        END IF
      END DO
!      if(debug_bc)write(*,*)'debug_bc',Site,SiteC,SNSite
      Cent=0.5*(DFLOAT(SiteC)+FLOAT(Site))
      ax=DFLOAT(Site)-DFLOAT(SiteC)
      snorm=DSQRT(ax(1)**2+ax(2)**2+ax(3)**2)
      ax=ax/snorm
!      write(*,*)ax
      DO i=1,3
         RSNZigZag=2*Cent-DFLOAT(SNsite(1:3,i))
         RSNArmChair=DFLOAT(SNsite(1:3,i))+8.6602254037*ax
         SNZiGZag(1:3,i)=NINT(RSNZigZag)
         SNArmChair(1:3,i)=NINT(RSNArmChair)
         IF(SNZiGZag(1,i).GE.LX)SNZiGZag(1,i)=SNZiGZag(1,i)-LX
         IF(SNZiGZag(1,i).LT.0)SNZiGZag(1,i)=SNZiGZag(1,i)+LX
         IF(SNZiGZag(2,i).GE.LY)SNZiGZag(2,i)=SNZiGZag(2,i)-LY
         IF(SNZiGZag(2,i).LT.0)SNZiGZag(2,i)=SNZiGZag(2,i)+LY
         IF(SNArmChair(1,i).GE.LX)SNArmChair(1,i)=SNArmChair(1,i)-LX
         IF(SNArmChair(1,i).LT.0)SNArmChair(1,i)=SNArmChair(1,i)+LX
         IF(SNArmChair(2,i).GE.LY)SNArmChair(2,i)=SNArmChair(2,i)-LY
         IF(SNArmChair(2,i).LT.0)SNArmChair(2,i)=SNArmChair(2,i)+LY
      END DO
!      write(*,*)'standard Arm',SNArmChair
!      write(*,*)'standard Zig',SNZiGZag
!      if(debug_bc)write(*,*)'debug_bc Arm',SNArmChair
!      if(debug_bc)write(*,*)'debug_bc Zig',SNZiGZag

      END SUBROUTINE FIND_SN


**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetSiC3C(Len)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),Buf(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      DO z=0,Len-4
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
!  Filled  Sites in 100 SiC lattice MC particles..NumAdatoms
      DO z=Len-3,Len-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si Atoms
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >            Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
            Prob = PtransE(1,2,0)! Si atom on SiC Surface
            Index_Event = 3
            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C Atoms
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >           Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
            Prob = PtransE(2,2,0)! C atom on SiC Surface
            Index_Event = 3
            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAdAtom",NumAdAtom,Len
      END DO
!  Empty  Sites in 100 SiC lattice (to be also consider MC particles in NumAdatoms)
      DO z=Len,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
            CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Buf=0
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               write(*,*)NextN(1:3,i)
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
            CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAtoms,NumAdAtom",NumAtoms,NumAdAtom
      END DO
      END SUBROUTINE SetSiC3C

**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetAFB(Len)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      DO z=0,Len-4
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF(x+y.LE.LenY)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF(x+y.LE.LenY)THEN  ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Filled  Sites in 100 SiC lattice MC particles..NumAdatoms
      DO z=Len-3,Len-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si Atoms
            IF(x+y.LE.LenY)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >            Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(1,2,0)! Si atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF(x+y.LE.LenY)THEN  ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >           Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(2,2,0)! C atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAdAtom",NumAdAtom,Len
      END DO
!  Empty  Sites in 100 SiC lattice (to be also consider MC particles in NumAdatoms)
      DO z=Len,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
 !           CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
 !           CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               write(*,*)NextN(1:3,i)
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
!            CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAtoms,NumAdAtom",NumAtoms,NumAdAtom
      END DO

! Set the probability
      DO z=Len-3,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! Si sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! C sites
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE ! Si Sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

      END SUBROUTINE SetAFB

**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetAFBSym(Len)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      Lenh=LenX/2-1+24
      DO z=0,Len-4
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Filled  Sites in 100 SiC lattice MC particles..NumAdatoms
      DO z=Len-3,Len-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si Atoms
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >            Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(1,2,0)! Si atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >           Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(2,2,0)! C atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAdAtom",NumAdAtom,Len
      END DO
!  Empty  Sites in 100 SiC lattice (to be also consider MC particles in NumAdatoms)
      DO z=Len,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
 !           CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
 !           CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               write(*,*)NextN(1:3,i)
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
!            CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAtoms,NumAdAtom",NumAtoms,NumAdAtom
      END DO

! Set the probability
      DO z=Len-3,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! Si sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! C sites
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE ! Si Sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

      END SUBROUTINE SetAFBSym


**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetAFBSymZimb(Len)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)

      Lenh=LenX/2-1+24
      DO z=0,Len/2 - 4 ! epi substrate
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

      DO z=Len/2-3, Len-4 
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=4 ! coordination
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=4
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Filled  Sites in 100 SiC lattice MC particles..NumAdatoms
      DO z=Len-3,Len-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si Atoms
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >            Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(1,2,0)! Si atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >           Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(2,2,0)! C atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAdAtom",NumAdAtom,Len
      END DO
!  Empty  Sites in 100 SiC lattice (to be also consider MC particles in NumAdatoms)
      DO z=Len,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
!           CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
!           CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               write(*,*)NextN(1:3,i)
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
!            CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAtoms,NumAdAtom",NumAtoms,NumAdAtom
      END DO

! Set the coordination 
      DO z=Len/2-3, Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      
! Set the probability
      DO z=21,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! Si sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! C sites
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE ! Si Sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

      END SUBROUTINE SetAFBSymZimb


**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetAFBSymZimb1(Len)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)

      Lenh=LenX/2-1+24
      DO z=0,Len/2 - 4 ! epi substrate
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            IF(z.eq.0)THEN
                Bon=4 ! coordination
                CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            END IF
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LE.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LE.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LE.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
            IF(z.EQ.0)THEN
                Bon=4
                CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            END IF
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

      DO z=Len/2-3, Len-4 
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF((ABS(x)+ABS(y).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LT.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LT.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE IF((ABS(x)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x)+ABS(y-LenY).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y-LenY).GT.Lenh+3))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=4 ! coordination
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LT.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LT.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE IF((ABS(x)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x)+ABS(y-LenY).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y-LenY).GT.Lenh+3))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=4
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Filled  Sites in 100 SiC lattice MC particles..NumAdatoms
      DO z=Len-3,Len-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si Atoms
            IF((ABS(x)+ABS(y).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LT.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LT.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE IF((ABS(x)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x)+ABS(y-LenY).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y-LenY).GT.Lenh+3))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >            Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(1,2,0)! Si atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((ABS(x)+ABS(y).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y).LT.Lenh).OR.
     >         (ABS(x)+ABS(y-LenY).LT.Lenh).OR.
     >         (ABS(x-LenX)+ABS(y-LenY).LT.Lenh))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE IF((ABS(x)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y).GT.Lenh+3).AND.
     >         (ABS(x)+ABS(y-LenY).GT.Lenh+3).AND.
     >         (ABS(x-LenX)+ABS(y-LenY).GT.Lenh+3))THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >           Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
!            Prob = PtransE(2,2,0)! C atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAdAtom",NumAdAtom,Len
      END DO
!  Empty  Sites in 100 SiC lattice (to be also consider MC particles in NumAdatoms)
      DO z=Len,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
!           CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
!           CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
!            Bon=2
!            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               write(*,*)NextN(1:3,i)
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
!            CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAtoms,NumAdAtom",NumAtoms,NumAdAtom
      END DO

! Set the coordination 
      DO z=1, Len  !Len/2, Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      
! Set the probability
      DO z=21,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! Si sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! C sites
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE ! Si Sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

      END SUBROUTINE SetAFBSymZimb1

**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetNC(Len1,Len2,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len1,Len2,Len3
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob
      write(*,*)LenZ/2-Len3/2,LenZ/2+Len3/2+1
      write(*,*)LenX/2-Len1/2,LenX/2+Len1/2+1

!  Fill Sites of an NC in 100 SiC lattice
      DO z=LenZ/2-Len3/2,LenZ/2+Len3/2-1
       DO y=LenY/2-Len2/2,LenY/2+Len2/2-1
        DO x=LenX/2-Len1/2,LenX/2+Len1/2-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=LenZ/2-Len3/2-10,LenZ/2+Len3/2+9
       DO y=LenY/2-Len2/2-10,LenY/2+Len2/2+9
        DO x=LenX/2-Len1/2-10,LenX/2+Len1/2+9
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=LenZ/2-Len3/2-10,LenZ/2+Len3/2+9
       DO y=LenY/2-Len2/2-10,LenY/2+Len2/2+9
        DO x=LenX/2-Len1/2-10,LenX/2+Len1/2+9
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetNC


**************************************************************************************
*
      SUBROUTINE SetInvPwAPB(Len,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF(x+y.LE.LenY)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF           
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF(x+y.LE.LenY)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF           
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(.NOT.(((ABS(x-Lenh)+ABS(y-Lenh)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y-LenY)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y-LenY)-(z-Len3)).LE.Len)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
              IF(x+y.LE.LenY)THEN ! Si sites
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              ELSE
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              END IF           
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
              IF(x+y.LE.LenY)THEN ! Si sites
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              ELSE
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              END IF           
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetInvPwAPB


**************************************************************************************
*
      SUBROUTINE SetInvPwOnlySiorC(Len,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF           
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF           
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(.NOT.(((ABS(x-Lenh)+ABS(y-Lenh)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y-LenY)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y-LenY)-(z-Len3)).LE.Len)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              ELSE
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              END IF           
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              ELSE
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              END IF           
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetInvPwOnlySiorC


**************************************************************************************
*
      SUBROUTINE SetInvPwOnlyCorSi(Len,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF           
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF           
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(.NOT.(((ABS(x-Lenh)+ABS(y-Lenh)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y-LenY)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y-LenY)-(z-Len3)).LE.Len)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              ELSE
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              END IF           
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
            IF((x.LT.LenX/2.AND.y.GE.LenY/2).OR.
     >         (x.GE.LenX/2.AND.y.LT.LenY/2))THEN ! Si sites
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              ELSE
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              END IF           
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetInvPwOnlyCorSi

      
**
**************************************************************************************
      SUBROUTINE SetInvP(Len,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(.NOT.(((ABS(x-Lenh)+ABS(y-Lenh)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y-LenY)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y-LenY)-(z-Len3)).LE.Len)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetInvP


**
**************************************************************************************
      SUBROUTINE SetTrench(Len,Len2,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len2,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
      IF(.NOT.((ABS(x-Lenh).LE.Len).AND.(z.GT.Len2.AND.z.LT.Len3)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetTrench

**************************************************************************************
      SUBROUTINE SetInvPOff(Len,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(.NOT.(((ABS(x-Lenh)+ABS(y-Lenh)-1.10*(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y)-1.10*(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y)-1.10*(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y-LenY)-1.10*(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y-LenY)-1.10*(z-Len3)).LE.Len)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetInvPOff


      
**************************************************************************************
**
**************************************************************************************
      SUBROUTINE SetInvPC(Len,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(.NOT.(((ABS(x-Lenh)+ABS(y-Lenh)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y-LenY)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y-LenY)-(z-Len3)).LE.Len)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetInvPC


**************************************************************************************
**
**************************************************************************************
      SUBROUTINE SetInvPSi(Len,Len3)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Len3,Lenh
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

      Lenh=LenX/2-1
      write(*,*)lenh,len,len3

      DO z=0,20
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
            END DO
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

!  Fill Sites of an inverted pyramid in 100 SiC lattice
      DO z=21,Len3-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(.NOT.(((ABS(x-Lenh)+ABS(y-Lenh)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x)+ABS(y-LenY)-(z-Len3)).LE.Len).OR.
     >         ((ABS(x-LenX)+ABS(y-LenY)-(z-Len3)).LE.Len)))THEN
             IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=21,Len3+9
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetInvPSi


**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetNCSp(Len1)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len1,d2
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)
      INTEGER :: Index_Event
      REAL(8) :: Prob

!  Fill Sites of an NC in 100 SiC lattice
      DO z=LenZ/2-Len1/2,LenZ/2+Len1/2-1
       DO y=LenY/2-Len1/2,LenY/2+Len1/2-1
        DO x=LenX/2-Len1/2,LenX/2+Len1/2-1
         d2=(x-LenX/2)*(x-LenX/2)+(y-LenY/2)*(y-LenY/2)+
     >(z-LenZ/2)*(z-LenZ/2)
         IF(d2.LE.Len1*Len1/4)THEN
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
           END IF
          END IF
         END IF
         END IF
        END DO
       END DO
      END DO

! Set the coordination
      DO z=LenZ/2-Len1/2-10,LenZ/2+Len1/2+9
       DO y=LenY/2-Len1/2-10,LenY/2+Len1/2+9
        DO x=LenX/2-Len1/2-10,LenX/2+Len1/2+9
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Coor=0
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)Coor=Coor+1
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
               LattCoo(x,y,z)=0
               DO i=1,4
                 CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                 CoorNN=CoorNN-1
                 IF(CoorNN.LT.0)CoorNN=0
                 CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
               END DO
             ELSE
               CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

! Set the probability
      DO z=LenZ/2-Len1/2-10,LenZ/2+Len1/2+9
       DO y=LenY/2-Len1/2-10,LenY/2+Len1/2+9
        DO x=LenX/2-Len1/2-10,LenX/2+Len1/2+9
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetNCSp




!     Deposition SUBROUTINE
!     it updates the systems's state after the deposition of an atom of the type IAtom
!     the IndVoid points to a empty site which have one of its NN already filled
      SUBROUTINE Deposition(IndVoid,IAtom) ! IAtom 1=01 Si 2=10 C
      USE DefDerType
      USE DefSystem
      USE Definitions

      IMPLICIT NONE
      INTEGER IndVoid,IAtom,IndNN,IndNNN,Index_Event
      INTEGER DisNN2,Dx,Dy
      INTEGER i,j,jj,Coor,CoorOld,Occ,OccN,CoorNN,CoorNNN,SiOrC
      INTEGER Site(3),SiteC(3),NextN(3,4),NextNN(3,4),NextNNN(3,4)
      INTEGER Newsite(3,2),NNArm(3,3),NNZig(3,3),SNSite(3,3)
      REAL(8):: Prob
      LOGICAL IF_Arm,IF_Arm_Choice,Recipr
      EXTERNAL IF_Arm_Choice

      Site=ListAdAtom(IndVoid) % AtomXYZ
      !write(*,*)'Deposition Site',Site
      Coor=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosCoor,LenCoor) ! Get the previous coordination
      CoorOld=Coor
      CALL MVBITS(1,0,LenOcc,
     >            LattCoo(Site(1),Site(2),Site(3)),PosOcc) ! Set occupancy
      CALL MVBITS(IAtom,0,LenSiC,
     >            LattCoo(Site(1),Site(2),Site(3)),PosSiC)  ! Set atom kind
      !write(*,*)" "
      !write(*,*)"Deposition Coor IndVoid",Coor
      IF(Coor.EQ.1)THEN ! the new full site has not get yet a fixed coordination
        SiteC = ListAdAtom(IndVoid) % NextNXYZ(1:3,1)
        IndNN = LattInd(SiteC(1),SiteC(2),SiteC(3))
        CoorNN=IBITS(LattCoo(SiteC(1),SiteC(2),SiteC(3)),
     >                  PosCoor,LenCoor)  ! not used
        SNSite=0
        j=0
        DO i=1,4
         IF(.NOT.ALL(ListAdAtom(IndNN) % NextNXYZ(1:3,i).EQ.
     >      Site(1:3)))THEN
            j=j+1  ! should reaches the value = 3
            SNSite(1:3,j)= ListAdAtom(IndNN) % NextNXYZ(1:3,i)
         END IF
        END DO
        !write(*,*)'Dep. Call SN Routine Site,SiteC',Site,SiteC
        !write(*,*)SNSite,IndNN
        CALL FIND_SN(Site,SiteC,SNSite,NNZig,NNArm,LenX,LenY)
        !write(*,*)NNZiG
        !write(*,*)NNArm
        IF_Arm = IF_Arm_Choice(Site,SiteC,NNZig,NNArm)  ! Fix the coordination type and store the NN positions
        IF(IF_Arm)THEN
          ListAdAtom(IndVoid) % NextNXYZ(1:3,2:4)=NNArm
        ELSE
          ListAdAtom(IndVoid) % NextNXYZ(1:3,2:4)=NNZig
        END IF
      END IF

      NextN = ListAdAtom(IndVoid) % NextNXYZ
      !write(*,*)"Deposition NextN",NextN
      DO i=1,4
        Occ = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosOcc,LenOcc)
        IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
        CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)
        !write(*,*)"Deposition IndNN,CoorNN,Occ",IndNN,CoorNN,Occ
        NextNN=0
        IF (Occ.EQ.0)THEN  ! empty site
            IF (CoorNN.EQ.0) THEN  ! This empty site is a new MC particle its Coordination is now 1 the type is not fixed
              CoorNN=1
              CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
              NextNN(1:3,1)=Site ! Only one NN is fixed
              CALL Get_Prob(0,0,CoorNN,NextN(1:3,i),NextNN,
     >                      Index_Event,Prob) ! Occ=1/0, NSiC 1=01 Si, 2=10 C, If Occ=0 NSiC=0
              IF(Index_Event.NE.0)THEN
                 !write(*,*)"Dep. Coor1 Ind_Ev,Prob",Index_Event,Prob
                 CALL AddMC(NextN(1:3,i),NextNN,Index_Event,Prob)
              END IF
            ELSE IF (CoorNN.EQ.1) THEN ! This empty site now fixes its next neighbours if a tetraedral configuration can be find
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Dx=ABS(NextNN(1,1)-Site(1))
              Dx=MIN(Dx,LenX-Dx)
              Dy=ABS(NextNN(2,1)-Site(2))
              Dy=MIN(Dy,LenX-Dy)
              DisNN2=Dx*Dx+Dy*Dy+(NextNN(3,1)-Site(3))**2
              IF(DisNN2.EQ.72)THEN  ! 72 is the second neighbour distance in tetrahedral configuration
                CoorNN = 2
                NextNN(1:3,2) = Site
                CALL FIND_NN(NextNN(1:3,1),NextNN(1:3,2),
     >                     NextN(1:3,i),NewSite,LenX,LenY) ! find remanent nn
                DO j=3,4
                   NextNN(1:3,j)=NewSite(1:3,j-2)
                   OccN = IBITS(LattCoo(NextNN(1,j),NextNN(2,j),
     >                          NextNN(3,j)),PosOcc,LenOcc)
                   IF(OccN.NE.0)THEN
                     !write(*,*)'NN occupied',NextN(1:3,i)
                     !write(*,*)'NN occupied',NextNN(1:3,j)
                     CoorNNN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j)
     >                  ,NextNN(3,j)),PosCoor,LenCoor)
                     IF(CoorNNN.LE.3)THEN !
                        IndNNN=LattInd(NextNN(1,j),NextNN(2,j),
     >                                 NextNN(3,j))
                        NextNNN = ListAdAtom(IndNNN) % NextNXYZ
                        Recipr=.FALSE.
                        DO jj=1,4 ! check is the site belongs to the NN of its NN
                           IF(ALL(NextNNN(1:3,jj).EQ.NextN(1:3,i)))
     >                        Recipr=.TRUE.
                        END DO
                        IF(Recipr)THEN
                          CoorNN = CoorNN+1
                        END IF
                     ELSE
                        Recipr=.FALSE.
                     END IF
                     !write(*,*)Recipr
                   END IF
                END DO
                CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
                ListAdAtom(IndNN) % NextNXYZ = NextNN
                CALL Get_Prob(0,0,CoorNN,NextN(1:3,i),NextNN,
     >                      Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                IF(Index_Event.EQ.0)THEN
                  CALL EraseMC(NextN(1:3,i))
                  CALL AddVoid(NextN(1:3,i),NextNN)
                ELSE ! 72 is the NN in tetrahedral configuration Site can belong to NextNN
                  ListAdAtom(IndNN) % Ind_Event = Index_Event
                  ListAdAtom(IndNN) % ProbTrans = Prob
                  CALL Updatetree(IndNN,Prob) !
                END IF
              !ELSE ! distance between NextNN and Site is not seconf NN one Site cannot be added to NextNN and CoorNN=1
                  !write(*,*)'DisNN2',DisNN2
                  !write(*,*)NextNN(1:3,1)
                  !write(*,*)Site
              !    STOP
              END IF
            ELSE IF(CoorNN.EQ.2.OR.CoorNN.EQ.3)THEN ! This empty site now increases the coordination if Site is in its NN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the site belongs to the NN of its NN
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
              END DO
              !write(*,*)'Deposition Recipr',Recipr
              IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                CoorNN = CoorNN + 1
                CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
                CALL Get_Prob(0,0,CoorNN,NextN(1:3,i),NextNN,
     >                      Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                IF(Index_Event.EQ.0)THEN
                  !write(*,*)"Dep. Vacancy Generation",NextN(1:3,i)
                  CALL EraseMC(NextN(1:3,i))
                  CALL AddVoid(NextN(1:3,i),NextNN)
                ELSE ! 72 is the NN in tetrahedral configuration Site can belong to NextNN
                  ListAdAtom(IndNN) % Ind_Event = Index_Event
                  ListAdAtom(IndNN) % ProbTrans = Prob
                  !write(*,*)"Dep.C3 Index_Event,Prob",Index_Event,Prob
                  CALL Updatetree(IndNN,Prob) !
                END IF
              ELSE ! Site in not in the NN list of Next Coord cannot increase
                !write(*,*)NextNN
                !write(*,*)Site
               ! STOP
              END IF
            END IF
        ELSE ! full site
          IF(CoorNN.LE.3)THEN ! Only if belong to MC particles
            NextNN = ListAdAtom(IndNN) % NextNXYZ
            Recipr=.FALSE.
            DO j=1,4 ! check is the site belongs to the NN of its NN
              IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
            END DO
            !write(*,*)'Deposition Recipr',Recipr
            IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
              CoorNN=CoorNN + 1
              CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
              SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
              CALL Get_Prob(1,SiOrC,CoorNN,NextN(1:3,i),NextNN,
     >                    Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
              IF(Index_Event.EQ.0)THEN
                !write(*,*)"Dep. Bulk Site Generation",Index_Event,Prob
                CALL EraseMC(NextN(1:3,i))
                CALL AddBulk(NextN(1:3,i),NextNN)
              ELSE
                ListAdAtom(IndNN) % Ind_Event = Index_Event
                ListAdAtom(IndNN) % ProbTrans = Prob
                !write(*,*)"Dep.full NN",Index_Event,Prob,SiOrC,CoorNN
                CALL Updatetree(IndNN,Prob) !
              END IF
              IF(CoorOld.EQ.1.AND.i.GE.2)THEN
                !write(*,*)'Increase Coor',Site,Coor
                Coor=Coor+1
              END IF
            END IF
          END IF
        END IF
      END DO
      IndVoid=LattInd(Site(1),Site(2),Site(3)) ! it can change due to the Next Neigh. update
      IF(Coor.NE.CoorOld)CALL MVBITS(Coor,0,LenCoor,LattCoo(Site(1),
     >                               Site(2),Site(3)),PosCoor)
      CALL Get_Prob(1,IAtom,Coor,Site,NextN,Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
      IF(Index_Event.EQ.0)THEN
                !write(*,*)"Dep. Bulk Site Generation",Index_Event,Prob
        CALL EraseMC(Site(1:3))
        CALL AddBulk(Site(1:3),NextN)
      ELSE
        ListAdAtom(IndVoid) % Ind_Event = Index_Event
        ListAdAtom(IndVoid) % ProbTrans = Prob
      !write(*,*)"Dep. IndVoid Index_Event,Prob",Index_Event,Prob
        CALL Updatetree(IndVoid,Prob)
      END IF
      END SUBROUTINE Deposition

      FUNCTION IF_Arm_ChoiceOld(SiteOld,SiteCOld,NZig,NArm)
      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE

      INTEGER i,j,IndArm,IndZig,NumArm,NumZig,CoorArm,CoorZig
      INTEGER NArm(3,3),NZig(3,3)
      INTEGER NNArm(3,3),NNZig(3,3),SNSite(3,3)
      INTEGER SiteOld(3),SiteCOld(3),Site(3),SiteC(3)
      LOGICAL IF_Arm_ChoiceOld,IF_Arm_Choice
      REAL(8) random
      EXTERNAL random

      CoorArm=0
      CoorZig=0
      SiteC=SiteOld
      DO i=1,3
       SNSite=NArm
       SNSite(1:3,i)=SiteCOld
       Site=NArm(1:3,i)
       CALL FIND_SN(Site,SiteC,SNSite,NNZig,NNArm,LenX,LenY)
       DO j=1,3
          CoorArm=CoorArm+IBITS(LattCoo(NNZig(1,j)
     >    ,NNZig(2,j),NNZig(3,j)),PosCoor,LenCoor) ! ???
       END DO
       SNSite=NZig
       SNSite(1:3,i)=SiteCOld
       Site=NZig(1:3,i)
       CALL FIND_SN(Site,SiteC,SNSite,NNZig,NNArm,LenX,LenY)
       DO j=1,3
          CoorZig=CoorZig+IBITS(LattCoo(NNZig(1,j)
     >    ,NNZig(2,j),NNZig(3,j)),PosCoor,LenCoor) ! ???
       END DO
      END DO
      IF(CoorZig.GT.CoorArm)THEN
         IF_Arm_Choice=.FALSE.
      ELSE IF(CoorZig.EQ.CoorArm)THEN
         IF(random(idum).le.PTransZig)THEN
           IF_Arm_Choice=.FALSE.
         ELSE
           IF_Arm_Choice=.TRUE.
         END IF
      ELSE
         IF_Arm_Choice=.TRUE.
      END IF
      IF_Arm_ChoiceOld=IF_Arm_Choice
      END FUNCTION IF_Arm_ChoiceOld


      FUNCTION IF_Arm_Choice(SiteOld,SiteCOld,NZig,NArm)
      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE

      INTEGER i,j,IndArm,IndZig,NumArm,NumZig,CoorArm,CoorZig
      INTEGER NArm(3,3),NZig(3,3)
      INTEGER NNArm(3,3),NNZig(3,3),SNSite(3,3)
      INTEGER SiteOld(3),SiteCOld(3),Site(3),SiteC(3)
      LOGICAL IF_Arm_Choice
      REAL(8) random
      EXTERNAL random

      CoorArm=0
      CoorZig=0
      SiteC=SiteOld
      DO i=1,3
       CoorArm=CoorArm+IBITS(LattCoo(NArm(1,i)
     >          ,NArm(2,i),NArm(3,i)),PosCoor,LenCoor) ! ???
       CoorZig=CoorZig+IBITS(LattCoo(NZig(1,i)
     >          ,NZig(2,i),NZig(3,i)),PosCoor,LenCoor) ! ???
      END DO
      IF(CoorZig.GT.CoorArm)THEN
         IF_Arm_Choice=.FALSE.
      ELSE IF(CoorZig.EQ.CoorArm)THEN
         IF(random(idum).le.PTransZig)THEN
           IF_Arm_Choice=.FALSE.
         ELSE
           IF_Arm_Choice=.TRUE.
         END IF
      ELSE
         IF_Arm_Choice=.TRUE.
      END IF
      !write(*,*)'IF_Arm_choice',CoorZig,CoorArm,IF_Arm_Choice
      END FUNCTION IF_Arm_Choice

!     Evaporation SUBROUTINE
!     it updates the systems's state after the evaporation of an atom
!     the IndAtom points to a full site before the update procedure
      SUBROUTINE Evaporation(IndAtom) ! IAtom 1=01 Si 2=10 C
      USE DefDerType
      USE DefSystem
      USE Definitions

      IMPLICIT NONE
      INTEGER IndAtom,IAtom,IndNN,IndNNN,Index_Event
      INTEGER i,j,jj,Coor,Occ,OccN,CoorNN,CoorNNN,SiOrC
      INTEGER Site(3),SiteC(3),NextN(3,4),NextNN(3,4),SNSite(3,3)
      INTEGER NNArm(3,3),NNZig(3,3),NBuff(3),NextNNN(3,4)
      REAL(8)Prob
      LOGICAL Recipr

      Site=ListAdAtom(IndAtom) % AtomXYZ
      IAtom=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosSiC,LenSiC)
      Coor=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosCoor,LenCoor) ! not used
      !write(*,*)'Evaporation Site',Site,Coor
      CALL MVBITS(0,0,LenOcc,
     >            LattCoo(Site(1),Site(2),Site(3)),PosOcc)  ! No occupancy
      CALL MVBITS(0,0,LenSiC,
     >            LattCoo(Site(1),Site(2),Site(3)),PosSiC)  ! No atom type

      NextN = ListAdAtom(IndAtom) % NextNXYZ ! All the NextN are fixed for
      DO i=1,4  ! The coordination of the next neighbour sites decreases of 1
        Occ = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosOcc,LenOcc)
        CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)
        IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
        !write(*,*)"Ev. NextN,Occ,CoorNN",NextN(1:3,i),Occ,CoorNN
        NextNN=0

        IF (Occ.EQ.0)THEN  ! empty site
            IF (CoorNN.EQ.1) THEN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
              END DO
              IF(Recipr)THEN ! this empty site has not now any NN is must be deleted from the MC list
                  CALL EraseMC(NextN(1:3,i)) ! out of the MC particles
                  LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))=0
              END IF
            ELSE IF (CoorNN.EQ.2) THEN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                 CoorNN=1
                 CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                 DO j=1,4
                   OccN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j),
     >                  NextNN(3,j)),PosOcc,LenOcc)
                   IF(OccN.EQ.0)THEN
                     ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                   ELSE
                     CoorNNN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j)
     >                  ,NextNN(3,j)),PosCoor,LenCoor)
                     IF(CoorNNN.LE.3)THEN !
                        IndNNN=LattInd(NextNN(1,j),NextNN(2,j),
     >                                 NextNN(3,j))
                        NextNNN = ListAdAtom(IndNNN) % NextNXYZ
                        Recipr=.FALSE.
                        DO jj=1,4
                          IF(ALL(NextNNN(1:3,jj).EQ.NextN(1:3,i)))
     >                      Recipr=.TRUE.
                        END DO
                        IF(Recipr)THEN
                          NBuff(1:3)=
     >                    ListAdAtom(IndNN) % NextNXYZ(1:3,j)
                          ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                        ELSE
                          ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                        END IF
                     ELSE ! If CoorNNN=4 NextNN Cannot be NN of Next
                        ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                     END IF
                   END IF
                 END DO
                 ListAdAtom(IndNN) % NextNXYZ(1:3,1)=NBuff(1:3)
                 CALL Get_Prob(0,0,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ListAdAtom(IndNN) % Ind_Event = Index_Event
                 ListAdAtom(IndNN) % ProbTrans = Prob
                 CALL Updatetree(IndNN,Prob) !
              END IF
            ELSE IF (CoorNN.EQ.3) THEN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                CoorNN=2
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                CALL Get_Prob(0,0,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                ListAdAtom(IndNN) % Ind_Event = Index_Event
                ListAdAtom(IndNN) % ProbTrans = Prob
                CALL Updatetree(IndNN,Prob) !
              END IF
            ELSE ! CoorNN.EQ.4 A previously Vacancy site now is a MC particle
              NextNN = ListVoid(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                CoorNN=3
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                CALL Get_Prob(0,0,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL EraseVoid(NextN(1:3,i))
                CALL AddMC(NextN(1:3,i),NextNN,Index_Event,Prob)
              END IF
            END IF
          ELSE
            IF (CoorNN.EQ.4)THEN! full site
              NextNN = ListAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                CoorNN=3
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
                CALL Get_Prob(1,SiOrC,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL EraseBulk(NextN(1:3,i))
                CALL AddMC(NextN(1:3,i),NextNN,Index_Event,Prob)
              END IF
            ELSE
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                CoorNN=CoorNN-1
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
                CALL Get_Prob(1,SiOrC,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                ListAdAtom(IndNN) % Ind_Event = Index_Event
                ListAdAtom(IndNN) % ProbTrans = Prob
                CALL Updatetree(IndNN,Prob) !
              END IF
            END IF
        END IF
      END DO
      IndAtom=LattInd(Site(1),Site(2),Site(3)) ! it can change due to the Next Neigh. update
      IF(Coor.EQ.0)THEN  ! Isolated Atom Evaporation
         CALL EraseMC(Site)
      ELSE
        CALL Get_Prob(0,0,Coor,Site,NextN,Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
        ListAdAtom(IndAtom) % Ind_Event = Index_Event
        ListAdAtom(IndAtom) % ProbTrans = Prob
        IF(Coor.EQ.1)THEN ! the new empty site will not have a fixed coordination
          NextN = ListAdAtom(IndAtom) % NextNXYZ
          DO i=1,4
           OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosOcc,LenOcc)
           IF(OccN.EQ.0)THEN
             ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
           ELSE ! Note that Coor
             CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)
             IF(CoorNN.LE.3)THEN
               IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
               NextNN = ListAdAtom(IndNN) % NextNXYZ
               Recipr=.FALSE.
               DO j=1,4 ! check is the Site belongs to the NN of its NN
                 IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
               END DO
             !write(*,*)'Evap 2 Recipr',Recipr
               IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                 NBuff(1:3)=
     >           ListAdAtom(IndAtom) % NextNXYZ(1:3,i)
                 ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
               ELSE
                 ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
               END IF
             ELSE
               !write(*,*)'Evap Coordination = 4'
               ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
             END IF
           END IF
          END DO
          ListAdAtom(IndAtom) % NextNXYZ(1:3,1)=NBuff(1:3)
        END IF
        CALL Updatetree(IndAtom,Prob)
      END IF
      END SUBROUTINE Evaporation

      SUBROUTINE Get_Prob(Occ,NSiC,Coor,Site,NextN,
     >                    Index_Event,Prob_Event) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Occ,OccN,NdiC,Coor,CoorNN,IndNN
      INTEGER Site(3),NextN(3,4),NextNN(3,4)
      INTEGER i,j,Index_Event,SiOrC,NSiC,N_C,N_Si
      REAL(8) Prob_Event
      LOGICAL Recipr
      REAL(8):: Ptot,Pcurr
      REAL(8) :: random
      EXTERNAL random

      IF (Occ.EQ.0)THEN ! Depositions
         IF (Coor.EQ.1)THEN
            CoorNN=IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosCoor,LenCoor)
            IF(CoorNN.GE.2)THEN ! in the case Coor=1 the single first neighbor should have CoorNN>=2
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
              Prob_Event = 1.e-8
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
         ELSE IF(Coor.GE.4)THEN! Vacancy Site?  Error
            Index_Event = 0
            Prob_Event = 0.d0
            !write(*,*)Coor,'Get Prob Vacancy Generation'
            IF(Coor.GT.4)THEN
               write(*,*)'anomalous coordination'
            END IF
         ELSE
            write(*,*)'anomalous coordination'
            write(*,*)Coor
         END IF
      ELSE ! Occ=1 Evaporation
         N_Si=0
         N_C =0
         !write(*,*)'Get Prob Site',Site,LattInd(Site(1),Site(2),Site(3))
         DO i=1,4
           OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                PosOcc,LenOcc)
           IF(OccN.NE.0)THEN
             CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)
             IF(CoorNN.GE.1.AND.CoorNN.LE.3)THEN
               IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
               NextNN = ListAdAtom(IndNN) % NextNXYZ
               Recipr=.FALSE.
               DO j=1,4 ! check is the site belongs to the NN of its NN
                 IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
               END DO
!           write(*,*)'Get Prob Recipr',Recipr
               IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                 SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                 ,PosSiC,LenSiC)
                 IF(SiOrC.EQ.1)N_Si=N_Si+1
                 IF(SiOrC.EQ.2)N_C=N_C+1
!            write(*,*)'Get Prob SiOrC',NextN(1:3,i),SiOrC
               END IF
             ELSE
               IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
               NextNN = ListAtom(IndNN) % NextNXYZ
               Recipr=.FALSE.
               DO j=1,4 ! check is the site belongs to the NN of its NN
                 IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
               END DO
!           write(*,*)'Get Prob Recipr',Recipr
               IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                 SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                 NextN(3,i)),PosSiC,LenSiC)
                 IF(SiOrC.EQ.1)N_Si=N_Si+1
                 IF(SiOrC.EQ.2)N_C=N_C+1
               END IF
             END IF
           END IF
         END DO
         IF(N_Si+N_C.NE.Coor)THEN
            write(*,*)Coor,N_Si,N_C,'Get Prob Error N_Si+N_C.NE.Coor'
            write(*,*)Site,NextN
            STOP
         END IF
         IF (N_Si+N_C.GT.3)THEN
            Index_Event = 0
            Prob_Event = 0.d0
            !write(*,*)Coor,'Get Prob Bulk Site'
         ELSE IF(N_Si+N_C.EQ.0)THEN
            Index_Event = 3
            Prob_Event = 100000.*Tree(1)  ! it makes its evaporation a fast event
            !write(*,*)Coor,'Get Prob Isolated Atom'
         ELSE
            Index_Event = 3
            Prob_Event = PTransE(NSic,N_Si,N_C)
         END IF
      END IF
      IF(Index_Event.NE.0.AND.Prob_Event.LT.1e-12)THEN
        write(*,*)Index_Event,Prob_Event,Coor
      END IF
      END SUBROUTINE Get_Prob

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


      SUBROUTINE SetProbability()

      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
! Ptrans (1,N_Si,N_C) Silicon evaporation
      PtransE=0.d0 ! not allowed case (e.g. 1,0,0) or (1,3,1) the default is 0
      READ(IPF,*)PtransE(1,0,1)
      READ(IPF,*)PtransE(1,1,0)
      READ(IPF,*)PtransE(1,1,1)
      READ(IPF,*)PtransE(1,0,2)
      READ(IPF,*)PtransE(1,2,0)
      READ(IPF,*)PtransE(1,2,1)
      READ(IPF,*)PtransE(1,1,2)
      READ(IPF,*)PtransE(1,0,3)
      READ(IPF,*)PtransE(1,3,0)
! Ptrans (2,N_Si,N_C) Carbon evaporation
      READ(IPF,*)PtransE(2,0,1)
      READ(IPF,*)PtransE(2,1,0)
      READ(IPF,*)PtransE(2,1,1)
      READ(IPF,*)PtransE(2,0,2)
      READ(IPF,*)PtransE(2,2,0)
      READ(IPF,*)PtransE(2,2,1)
      READ(IPF,*)PtransE(2,1,2)
      READ(IPF,*)PtransE(2,0,3)
      READ(IPF,*)PtransE(2,3,0)
!
      READ(IPF,*)PtransD(1)
      READ(IPF,*)PtransD(2)
      READ(IPF,*)PtransD(3)
      READ(IPF,*)PtransD(4)
      READ(IPF,*)PtransD(5)
      READ(IPF,*)PtransD(6)

      READ(IPF,*)PtransZig

      END SUBROUTINE SetProbability

      SUBROUTINE SetProbabilityOld()

      USE DefDerType
      USE DefSystem
      USE Definitions
! Ptrans (1,N_Si,N_C) Silicon evaporation
      PtransE=0.d0 ! not allowed case (e.g. 1,0,0) or (1,3,1) the default is 0
      PtransE(1,0,1)=3.0
      PtransE(1,1,0)=50.0
      PtransE(1,1,1)=20.0
      PtransE(1,0,2)=2.
      PtransE(1,2,0)=20.0
      PtransE(1,2,1)=10.
      PtransE(1,1,2)=10.
      PtransE(1,0,3)=0.1
      PtransE(1,3,0)=10.
! Ptrans (2,N_Si,N_C) Carbon evaporation
      PtransE(2,0,1)=50.0
      PtransE(2,1,0)=3.0
      PtransE(2,1,1)=20.
      PtransE(2,0,2)=20.
      PtransE(2,2,0)=2.
      PtransE(2,2,1)=10.
      PtransE(2,1,2)=10.
      PtransE(2,0,3)=10.
      PtransE(2,3,0)=0.1
!
      PtransD(1)=1.0
      PtransD(2)=1.0
      PtransD(3)=2.0
      PtransD(4)=2.0
      PtransD(5)=3.0
      PtransD(6)=3.0

      PtransZig=1.0

      END SUBROUTINE SetProbabilityOld


******|****************************************************************

!      SUBROUTINE UpdateDepos(Itrans,Site,NewSite)

!      OldCoo=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosCoor,LenCoor)
!      IF(OldCoo.EQ.2.OR.OldCoo.EQ.3)THEN
!       FIND_NN(                         ! find coordinate or remaining NN
!      ELSE IF(OldCoo
!      END SUBROUTINE UpdateDepos

***********************************************************************
       SUBROUTINE UpdateCoor_1(PosFrom,LenFrom,PosTo,LenTo,CoorNN)

        USE Definitions

        IMPLICIT NONE
        INTEGER :: PosFrom,LenFrom,PosTo,LenTo,CoorNN

        INTEGER    :: CoorFrom,CoorTo

        CoorTo   = IBITS(CoorNN,PosTo,LenTo) + 1
        CoorFrom = IBITS(CoorNN,PosFrom,LenFrom) - 1
        CALL MVBITS(CoorFrom,0,LenFrom,CoorNN,PosFrom)
        CALL MVBITS(CoorTo,0,LenTo,CoorNN,PosTo)
       END SUBROUTINE UpdateCoor_1
******|****************************************************************

***********************************************************************
       SUBROUTINE UpdateCoor_2(PosFrom,LenFrom,CoorNN)

        USE Definitions

        IMPLICIT NONE
        INTEGER :: PosFrom,LenFrom,CoorNN

        INTEGER    :: CoorFrom

        CoorFrom = IBITS(CoorNN,PosFrom,LenFrom) - 1
        CALL MVBITS(CoorFrom,0,LenFrom,CoorNN,PosFrom)
       END SUBROUTINE UpdateCoor_2
******|****************************************************************

***********************************************************************
       SUBROUTINE UpdateCoor_3(PosTo,LenTo,CoorNN)

        USE Definitions

        IMPLICIT NONE
        INTEGER :: PosTo,LenTo,CoorNN

        INTEGER   :: CoorTo

        CoorTo   = IBITS(CoorNN,PosTo,LenTo) + 1
        CALL MVBITS(CoorTo,0,LenTo,CoorNN,PosTo)
       END SUBROUTINE UpdateCoor_3
******|***************************************************************



***********************************************************************
       SUBROUTINE ExtrCoor(Site,Coor)

        USE DefSystem

        IMPLICIT NONE
        INTEGER :: Site,Coor(3)
        INTEGER :: res

        Coor(3) =Site/LenX/LenY
        res = Site - Coor(3)*LenX*LenY
        Coor(2) =res/LenX
	    Coor(1) = Res - Coor(2)*LenX

       END SUBROUTINE ExtrCoor
******|****************************************************************

*******************************************************
*  Allocate Binary Tree for the Transition probability
*******************************************************
       SUBROUTINE AllocateArrays()
       USE DefSystem
       IMPLICIT NONE
       INTEGER ind

!!!!!!it calculates necessary number of levels !!!!!!!
       DO ind=0,31
         IF(BTEST(NumMcMax,ind))Levels=ind ! level necessary in the tree
       END DO
       IF(NumMcMax.GT.ISHFT(1,Levels))THEN
         Levels = Levels + 2
         SizeTree = ISHFT(1,Levels) - 1
       ELSE
         Levels = Levels + 1
         SizeTree = ISHFT(1,Levels) - 1
       END IF
       WRITE(*,*) 'numParticelle=',NumMcMax,
     >     'Levels =',Levels,', SizeTree =',SizeTree

       ALLOCATE(Tree(SizeTree))

       END SUBROUTINE AllocateArrays
******************************************************

********************************************************************
*     It picks the leaf of the current even in the tree and returns
*     index in the list and related event probability
********************************************************************
      SUBROUTINE PickTreeEvent(Index,Patom)

!     Index - output
!     Patom - output
!     Levels,SizeTree,Tree - input, info su albero da usare
       USE Definitions
       USE DefSystem
       IMPLICIT NONE
       INTEGER :: Index
       REAL(8) :: Patom,Ptot
       INTEGER Parent,i
       REAL(8) Random
       EXTERNAL Random
       Parent = 1
       Index=1
       DO
	   IF (Parent.GT.Levels-1) EXIT
	   Patom=Random(idum)*Tree(Index)
	   Index = ISHFT(Index,1)
	   IF(Tree(Index).LT.Patom) Index=Index+1
	   Parent=Parent+1
       END DO
       Patom = Tree(Index)

       IF (Index-ISHFT(1,Levels-1)+1.GT.NumAdAtom)THEN
         Ptot=0
         DO i=ISHFT(1,Levels-1),ISHFT(1,Levels-1)+NumAdAtom-1
            Ptot=Ptot+Tree(i)
            write(*,*)i,i - ISHFT(1,Levels-1)+1,Tree(i)
         END DO
         write(*,*)'write Index> NumAdAtom',Ptot,Tree(1)
         STOP
       END IF

       Index = Index - ISHFT(1,Levels-1)+1

      END SUBROUTINE PickTreeEvent

**********************************************************************
*     It modifies to Patom the probability related to the particle
*     and recursively update the tree
***********************************************************************
      SUBROUTINE UpdateTree(Index,Patom)

!       modifica la particella Index, assegnadogli la probabilita
!       Patom, ATTENTO che Patom rappresenta prob totale
!       aggiorna poi l'albero di conseguenza

	USE DefSystem
	USE Definitions
	IMPLICIT NONE
	INTEGER :: Index
	Real(8)  Patom
	INTEGER Pos,Parent,NewPos,NeighPos,ControlBit

	IF(Patom.EQ.0.and.Index.lt.NumAdAtom) THEN
		write(*,*) Index,ListAdatom(Index) % AtomXYZ,Patom
		!CALL ERR("stabilizzato  da dentro Update!")
		write(*,*)"error stabilizzato  da dentro Update!"
		STOP
	ELSE if (Index.gt.NumAdAtom)THEN
		write(*,*) Index,NumAdAtom
		!CALL ERR("stabilizzato  da dentro Update!")
		write(*,*)"Index > NumAdAtom"
                STOP
	ENDIF
	Pos = Index + ISHFT(1,Levels-1)-1
	Tree(Pos) = Patom
	Parent=Levels-1
	DO
	   IF (Parent.LT.1) EXIT
	   ControlBit = IAND(Pos,1)
	   SELECT CASE (ControlBit)
	   CASE (0)
	      NeighPos = Pos + 1
	   CASE (1)
	      NeighPos = Pos - 1
	   END SELECT
	   NewPos = ISHFT(Pos,-1)
	   Tree(NewPos) = Tree(Pos) + Tree(NeighPos)
	   Pos = NewPos
	   Parent = Parent - 1
	END DO
       END SUBROUTINE UpdateTree

********************************************************************
*  It adds a MC particle to the lists and in the lattice
********************************************************************
       SUBROUTINE AddMC(Site,NextN,Index_Event,Prob)
!    Questo metodo aggiorna l'array delle cinetiche: AtprobList
!    inserisce la particella nell'albero e ne aggiorna le probabilita
!    MODIFICA: AtprobList, Tree (attraverso UpdateTree), LattInd, Latt

	USE DefDerType
	USE DefSystem
	IMPLICIT NONE
	Real*8 :: Prob
	INTEGER Index_Event
	INTEGER, DIMENSION(3) :: Site
	INTEGER, DIMENSION(3,4) :: NextN

	IF(Prob.EQ.0) THEN
           WRITE(*,*)"provo ad aggiungere prob nulla"
           STOP
        END IF
	NumAdAtom = NumAdAtom + 1
	IF(NumAdAtom.GE.NumMCMax) then
	  write(*,*) NumAdAtom,NumMCMax
	  STOP "troppe particelle cinetiche"
	endif
	LattInd(Site(1),Site(2),Site(3)) = NumAdAtom
	ListAdAtom(NumAdAtom) % AtomXYZ = Site
	ListAdAtom(NumAdAtom) % NextNXYZ = NextN
	ListAdAtom(NumAdAtom) % Ind_Event = Index_Event
	ListAdAtom(NumAdAtom) % ProbTrans = Prob
	CALL UpdateTree(NumAdAtom,Prob)
	END SUBROUTINE AddMC

********************************************************************
*  It adds a bulk particle to the lists and in the lattice
********************************************************************
       SUBROUTINE AddBulk(Site,NextN)
!    Questo metodo aggiorna l'array delle cinetiche: AtprobList
!    inserisce la particella nell'albero e ne aggiorna le probabilita
!    MODIFICA: AtprobList, Tree (attraverso UpdateTree), LattInd, Latt

	USE DefDerType
	USE DefSystem
	IMPLICIT NONE
	INTEGER, DIMENSION(3) :: Site
	INTEGER, DIMENSION(3,4) :: NextN

	NumAtoms = NumAtoms + 1
	IF(NumAtoms.GE.NumAtMax) then
	  write(*,*) NumAtoms,NumAtMax
	  STOP "troppe particelle bulk"
	endif
	LattInd(Site(1),Site(2),Site(3)) = NumAtoms
	ListAtom(NumAtoms) % AtomXYZ = Site
	ListAtom(NumAtoms) % NextNXYZ = NextN
	END SUBROUTINE AddBulk


********************************************************************
*  It adds a Vacancy to the lists and in the lattice
********************************************************************
       SUBROUTINE AddVoid(Site,NextN)
!    Questo metodo aggiorna l'array delle cinetiche: AtprobList
!    inserisce la particella nell'albero e ne aggiorna le probabilita
!    MODIFICA: AtprobList, Tree (attraverso UpdateTree), LattInd, Latt

	USE DefDerType
	USE DefSystem
	IMPLICIT NONE
	INTEGER, DIMENSION(3) :: Site
	INTEGER, DIMENSION(3,4) :: NextN

	NumVoids = NumVoids + 1
	IF(NumVoids.GE.NumVdMax) then
	  write(*,*) NumVoids,NumVdMax
	  STOP "troppe vacanze"
	endif
	LattInd(Site(1),Site(2),Site(3)) = NumVoids
	ListVoid(NumVoids) % AtomXYZ = Site
	ListVoid(NumVoids) % NextNXYZ = NextN
	END SUBROUTINE AddVoid



***********************************************************************
*    Erase a particle from the MC list
************************************************************************
      SUBROUTINE EraseMC(Site)
!	It modifies:
!	AtProbList(LattInd(x,y,z)), LattInd(x,y,z)


      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Site(3),Itype
      INTEGER, DIMENSION(3)   :: SiteTemp
      INTEGER, DIMENSION(3,4) :: NextNTemp
      INTEGER x,y,z,IndNN
      INTEGER n,ind
      x=Site(1)
      y=Site(2)
      z=Site(3)
      IF(z.eq.-1)THEN
       !CALL ERR("accesso a latt fuori bond")
       write(*,*)"accesso a latt fuori bond"
       STOP
      END IF
      indNN=LattInd(x,y,z)
      IF(IndNN.GT.NumMcMax.or.IndNN.EQ.0) THEN
        write(*,*) "err: EraseMC",IndNN,x,y,z
        STOP
      ENDIF

      IF(IndNN.GE.NumAdAtom)THEN
            ListAdAtom(IndNN) % AtomXYZ = 0
            ListAdAtom(IndNN) % NextNXYZ = 0
            ListAdAtom(IndNN) % ProbTrans = 0.d0
            ListAdAtom(IndNN) % Ind_Event = 0
            LattInd(Site(1),Site(2),Site(3))=0
            CALL UpdateTree(IndNN,0.d0)
            NumAdAtom = NumAdAtom - 1
      ELSE
            ListAdAtom(IndNN) = ListAdAtom(NumAdAtom)
       	    ind = NumAdAtom + ISHFT(1,Levels-1)-1
       	    CALL UpdateTree(IndNN,Tree(ind))
       	    CALL UpdateTree(NumAdAtom,0.d0)
       	    SiteTemp =  ListAdAtom(NumAdAtom) % AtomXYZ
!            NextNTemp = AtprobList(NumAdAtom) % NextNXYZ
       	    LattInd(SiteTemp(1),SiteTemp(2),SiteTemp(3)) = IndNN
            ListAdAtom(NumAdAtom) % ProbTrans = 0.d0
            ListAdAtom(NumAdAtom) % AtomXYZ = 0
            ListAdAtom(NumAdAtom) % NextNXYZ = 0
            ListAdAtom(NumAdAtom) % Ind_Event = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumAdAtom = NumAdAtom - 1
      END IF
c	      IF(NumMc.EQ.2) STOP "sono finite le particelle cinetiche...."
      END SUBROUTINE EraseMC
***********************************************************************

***********************************************************************
*    Erase a particle from the bulk particle list
************************************************************************
      SUBROUTINE EraseBulk(Site)
!	It modifies:
!	AtProbList(LattInd(x,y,z)), LattInd(x,y,z)


      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Site(3),Itype
      INTEGER, DIMENSION(3)   :: SiteTemp
      INTEGER, DIMENSION(3,4) :: NextNTemp
      INTEGER x,y,z,IndNN
      INTEGER n,ind
      x=Site(1)
      y=Site(2)
      z=Site(3)
      IF(z.eq.-1)THEN
       !CALL ERR("accesso a latt fuori bond")
       write(*,*)"accesso a latt fuori bond"
       STOP
      END IF
      indNN=LattInd(x,y,z)
      IF(IndNN.GT.NumAtMax.or.IndNN.EQ.0) THEN
        write(*,*) "err: EraseBulk",IndNN,x,y,z
        STOP
      ENDIF
      !write(*,*)'erase bulk indnn site',indnn,site
      IF(IndNN.EQ.NumAtoms)THEN
            ListAtom(NumAtoms) % AtomXYZ = 0
            ListAtom(NumAtoms) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumAtoms = NumAtoms - 1
      ELSE
            ListAtom(IndNN) = ListAtom(NumAtoms)
       	    SiteTemp =  ListAtom(NumAtoms) % AtomXYZ
       	    !write(*,*)SiteTemp
!            NextNTemp = AtprobList(NumAdAtom) % NextNXYZ
       	    LattInd(SiteTemp(1),SiteTemp(2),SiteTemp(3)) = IndNN
            ListAtom(NumAtoms) % AtomXYZ = 0
            ListAtom(NumAtoms) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumAtoms = NumAtoms - 1
      END IF
c	      IF(NumMc.EQ.2) STOP "sono finite le particelle cinetiche...."
      END SUBROUTINE EraseBulk
***********************************************************************


***********************************************************************
*    Erase a particle from the bulk particle list
************************************************************************
      SUBROUTINE EraseVoid(Site)
!	It modifies:
!	AtProbList(LattInd(x,y,z)), LattInd(x,y,z)


      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Site(3),Itype
      INTEGER, DIMENSION(3)   :: SiteTemp
      INTEGER, DIMENSION(3,4) :: NextNTemp
      INTEGER x,y,z,IndNN
      INTEGER n,ind
      x=Site(1)
      y=Site(2)
      z=Site(3)
      IF(z.eq.-1)THEN
       !CALL ERR("accesso a latt fuori bond")
       write(*,*)"accesso a latt fuori bond"
       STOP
      END IF
      indNN=LattInd(x,y,z)
      IF(IndNN.GT.NumVdMax.or.IndNN.EQ.0) THEN
        write(*,*) "err: EraseVoid",IndNN,x,y,z
        STOP
      ENDIF
      !write(*,*)'erase void indnn site',indnn,site
      IF(IndNN.EQ.NumVoids)THEN
            ListVoid(NumVoids) % AtomXYZ = 0
            ListVoid(NumVoids) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumVoids = NumVoids - 1
      ELSE
            ListVoid(IndNN) = ListVoid(NumVoids)
       	    SiteTemp =  ListVoid(NumVoids) % AtomXYZ
       	    !write(*,*)SiteTemp
!            NextNTemp = AtprobList(NumAdAtom) % NextNXYZ
       	    LattInd(SiteTemp(1),SiteTemp(2),SiteTemp(3)) = IndNN
            ListVoid(NumVoids) % AtomXYZ = 0
            ListVoid(NumVoids) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumVoids = NumVoids - 1
      END IF
c	      IF(NumMc.EQ.2) STOP "sono finite le particelle cinetiche...."
      END SUBROUTINE EraseVoid
***********************************************************************


***********************************************************************
      SUBROUTINE FileOpen()

      USE Definitions
      IMPLICIT    NONE

C   File descriptor 18 used in 'SUBROUTINE ReadDefectData'
      IPF = 19
      OPF0 = 20
      OPF1 = 21
      OPF2 = 22
      OPF3 = 23
      OPF4 = 24
      OPF5 = 25
      OPF6 = 26
      OPF7 = 27
      OPF8 = 28
      OPF9 = 29
      OPF10 = 30
      OPF11 = 31
      OPF12 = 32
      OPF13 = 33
      OPF14 = 34
      OPF15 = 35

      OPEN(IPF,FILE='start.dat',STATUS='OLD')
      OPEN(OPF0,FILE='Run_time.dat',STATUS='REPLACE')
!      OPEN(OPF1,FILE='check01.dat',STATUS='REPLACE')
!      OPEN(OPF2,FILE='check02.dat',STATUS='REPLACE')
!      OPEN(OPF3,FILE='check03.dat',STATUS='REPLACE')
!      OPEN(OPF4,FILE='check04.dat',STATUS='REPLACE')
!      OPEN(OPF5,FILE='check05.dat',STATUS='REPLACE')
!      OPEN(OPF6,FILE='check06.dat',STATUS='REPLACE')
!      OPEN(OPF7,FILE='check07.dat',STATUS='REPLACE')
!      OPEN(OPF8,FILE='check08.dat',STATUS='REPLACE')
!      OPEN(OPF9,FILE='check09.dat',STATUS='REPLACE')
!      OPEN(OPF10,FILE='check10.dat',STATUS='REPLACE')
!      OPEN(OPF11,FILE='check11.dat',STATUS='REPLACE')
!      OPEN(OPF12,FILE='check12.dat',STATUS='REPLACE')
!      OPEN(OPF13,FILE='check13.dat',STATUS='REPLACE')
!      OPEN(OPF14,FILE='check14.dat',STATUS='REPLACE')
!      OPEN(OPF15,FILE='check15.dat',STATUS='REPLACE')

      END SUBROUTINE FileOpen
******|****************************************************************

***********************************************************************
      SUBROUTINE FileClose()

      USE Definitions

      IMPLICIT    NONE

      CLOSE(IPF)
      CLOSE(OPF0)
!      CLOSE(OPF1)
!      CLOSE(OPF2)
!      CLOSE(OPF3)
!      CLOSE(OPF4)
!      CLOSE(OPF5)
!      CLOSE(OPF6)
!      CLOSE(OPF7)
!      CLOSE(OPF8)
!      CLOSE(OPF9)
!      CLOSE(OPF10)
!      CLOSE(OPF11)
!      CLOSE(OPF12)
!      CLOSE(OPF13)
!      CLOSE(OPF14)
!      CLOSE(OPF15)

      END SUBROUTINE FileClose
******|****************************************************************

***********************************************************************
      real(8) function Random(idum)

       Real(8) :: Fac
       INTEGER :: Idum
       INTEGER :: Mbig,Mseed,Mz,Iff,Mj,Mk,I,II,Ma,K,Inext,Inextp


         Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)

         Dimension MA(55)  !the dimension of MA is  special
         save
!       inizialize MA(55) using idum and mseed
         if (idum.lt.0.or.iff.eq.0) then
          iff=1
          mj=mseed-iabs(idum)
          mj=mod(mj,mbig)
          ma(55)=mj
          mk=1
!        inizialize the rest of MA in slightly random order with
!        Number that are not especially random
          do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if (mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
11        continue
!         randomize MA NumBer by the "warming up generator"
          do 13 k=1,4
            do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if (ma(i).lt.mz)ma(i)=ma(i)+mbig
12          continue
13        continue
          inext=0
          inextp=31     ! 31 is special
          idum=1
        end if
        inext=inext+1
        if (inext.eq.56) inext=1
         inextp=inextp+1
        if (inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)     ! generate a new random NumBer
        if (mj.lt.mz)mj=mj+mbig     ! Be sure that it is in the range
        ma(inext)=mj                ! Store it
        random=mj*fac                 ! outpup ran3 belongs (0,1)

      end function random


***********************************************************************
      SUBROUTINE WriteMolMolSource(Time,Iter,Iout)

      USE  Definitions
      USE  DefSystem
      IMPLICIT    NONE
      INTEGER     Iter
      REAL(8)     Time

      INTEGER FN
      CHARACTER(Len=13) :: XYZFileName,SourceFileName,GifFileName
      CHARACTER(Len=15) :: XYZ_dFileName,XYZ_vFileName,XYZ_wFileName
      CHARACTER(Len=9)  :: FileNameBase
      CHARACTER(Len=4)  :: XYZExt='.xyz',SourceExt='.src'
      CHARACTER(Len=6)  :: DefectExt='_d.xyz',VacancyExt='_v.xyz'
      CHARACTER(Len=6)  :: WrongExt='_w.xyz'
      CHARACTER(Len=4)  :: GifExt='.gif'
      INTEGER :: Iout

      FN = 97
      CALL GetOutputFileName(Time, Iter/Iout, FileNameBase)
      write(*,*)FileNameBase
      XYZFileName = FileNameBase//XYZExt
      XYZ_dFileName = FileNameBase//DefectExt
      XYZ_vFileName = FileNameBase//VacancyExt
      XYZ_wFileName = FileNameBase//WrongExt
      SourceFileName = FileNameBase//SourceExt
      GifFileName = FileNameBase//GifExt
      OPEN(FN+6,FILE=XYZ_wFileName(:LEN_TRIM(XYZ_wFileName)))
      OPEN(FN+5,FILE=XYZ_vFileName(:LEN_TRIM(XYZ_vFileName)))
      OPEN(FN+4,FILE=XYZ_dFileName(:LEN_TRIM(XYZ_dFileName)))
      OPEN(FN,FILE=SourceFileName(:LEN_TRIM(SourceFileName)))
      OPEN(FN+1,FILE=XYZFileName(:LEN_TRIM(XYZFileName)))
      OPEN(FN+2,FILE='movie.mos')
      OPEN(FN+3,FILE='DefAtEnergy.dat')
      CALL WriteMolMolXYZFile(FN+1)
      write(OPF10,*)time

      WRITE(FN,*) 'load xyz ', XYZFileName
      WRITE(FN,*) 'color temperature'
      WRITE(FN,*) 'set background white'

      WRITE(FN,*) 'set boundbox true'
      WRITE(FN,*) 'zoom 140'
      WRITE(FN,*) 'rotate x 41'
      WRITE(FN,*) 'rotate y 4'
C      WRITE(FN,*) 'spacefill 45'
      WRITE(FN,*) 'spacefill 2'
      WRITE(FN,*) 'color axes black'
      WRITE(FN,*) 'wireframe 2'
!      WRITE(FN,*) 'select atomno<',NumS
      WRITE(FN,*) 'color atom green'
!      WRITE(FN,*) 'select atomno>=',NumS
      WRITE(FN,*) 'wireframe 12'
      WRITE(FN,*) 'color atom purple'
      WRITE(FN,*) 'set bondmode and'

32    FORMAT('select atomno = ',I5)
33    FORMAT('select atomno = ',I5,', atomno = ',I5)
34    FORMAT(I5,2X,5I5,2X,F8.4)

      CLOSE(FN)
      WRITE(FN+3,*)''
      CLOSE(FN+3)
      WRITE(FN+2,*)'source    ',SourceFileName
      WRITE(FN+2,*)'write gif ',GifFileName
      WRITE(FN+2,*)'zap'
      CLOSE(FN+2)
      OPEN(FN+2,FILE='movie.mos',STATUS='OLD',POSITION='APPEND')
      OPEN(FN+3,FILE='DefAtEnergy.dat',STATUS='OLD',POSITION='APPEND')

      END SUBROUTINE WriteMolMolSource
******|****************************************************************

***********************************************************************
      SUBROUTINE WriteMolMolXYZFile(FN)

      USE  DefSystem
      USE  Definitions
      IMPLICIT NONE
      INTEGER IfSi,IfC,IcellSiC,JcellSiC,IGr,JGr,indAtom
      INTEGER x,y,z,x1,y1,z1,IfOcc,Coo,NCoorCorr
      INTEGER NumAtTotal,NoEpiBulk,NoEpiSurf,NvoidMerged,NumAtWrong
      INTEGER NumIsWrong,ind_atotype
      INTEGER FN,NumSurfAt,Site(3),i,ic,sitec(3),NextN(3,4)
      REAL(8) rad,sf,a1(3),a2(3),a3(3),cr(3)
      LOGICAL noepi
      REAL(8), dimension(3) :: alat

!      DATA a1 / 1,  1, -1/  !  BCC
!      DATA a2 /-1,  1,  1/
!      DATA a3 / 1, -1,  1/

!      DATA a1 / 1,  1,  0/   ! FCC
!      DATA a2 / 0,  1,  1/
!      DATA a3 / 1,  0,  1/

       DATA a1 / 1.,  0.,   0./   ! HEX
       DATA a2 / 0.5,  0.866025037844386,  0./
!       DATA a3 / 0.,  0.,   1./
       DATA a3 / 0.,  0.,   1.414213562373/


      sf = 4.63/12.0

      alat(1)=sf*real(LenX,kind=8)
      alat(2)=sf*real(LenY,kind=8)
      alat(3)=sf*real(LenZ,kind=8)

100   FORMAT('Si ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
101   FORMAT('B ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
102   FORMAT('C ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
103   FORMAT('Ge ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
104   FORMAT('O ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)

      NumAtTotal=0
      NumAtWrong=0
      NumIsWrong=0
      DO z=1,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF(IfOcc.EQ.1.AND.Coo.LE.3)THEN
                NumAtTotal=NumAtTotal+1
            ELSE IF (IfOcc.EQ.1.AND.Coo.EQ.4)THEN
                indAtom=LattInd(x,y,z)
                IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                NextN=ListAtom(indAtom) % NextNXYZ
                NCoorCorr=0
                IF(IfSi.EQ.1)THEN
                    DO i=1,4
                        NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosC,LenC)
                    END DO
                ELSE IF(IfC.EQ.1)THEN
                    DO i=1,4
                        NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosSi,LenSi)
                    END DO
                END IF
                IF(NCoorCorr.EQ.0)THEN ! only isolated 
                    NumIsWrong=NumIsWrong+1
                    IF(IfSi.EQ.1)THEN
                        LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosSi)
                        LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                        ListAtBuffer(NumIsWrong) % AtomXYZ(1)=x
                        ListAtBuffer(NumIsWrong) % AtomXYZ(2)=y
                        ListAtBuffer(NumIsWrong) % AtomXYZ(3)=z
                        ListAtBuffer(NumIsWrong) % Ind_Atype = 1
                    ELSE 
                        LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosC)
                        LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                        ListAtBuffer(NumIsWrong) % AtomXYZ(1)=x
                        ListAtBuffer(NumIsWrong) % AtomXYZ(2)=y
                        ListAtBuffer(NumIsWrong) % AtomXYZ(3)=z
                        ListAtBuffer(NumIsWrong) % Ind_Atype = 2
                    END IF
                END IF
            END IF
         END DO
       END DO
      END DO

      DO z=1,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF (IfOcc.EQ.1.AND.Coo.EQ.4)THEN
                indAtom=LattInd(x,y,z)
                IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                NextN=ListAtom(indAtom) % NextNXYZ
                NCoorCorr=0
                IF(IfSi.EQ.1)THEN
                    DO i=1,4
                        NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosC,LenC)
                    END DO
                ELSE IF(IfC.EQ.1)THEN
                    DO i=1,4
                        NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosSi,LenSi)
                    END DO
                END IF
                IF(NCoorCorr.EQ.3)NumAtWrong=NumAtWrong+1
            END IF
         END DO
        END DO
      END DO

      write(*,*)'NumAtWrong: ', NumAtWrong 

      NoEpiBulk=0
      DO z=0,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF(IfOcc.EQ.1.AND.Coo.EQ.4)THEN
              noepi=.TRUE.
              IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.
     >           AND.MOD(z,3).EQ.0)THEN
                x1=x/3
                y1=y/3
                z1=z/3
                IF((MOD(x1,2).EQ.MOD(y1,2)).AND.
     >             (MOD(x1,2).EQ.MOD(z1,2)))THEN
                  IF(MOD(x1+y1+z1,4).EQ.0.OR.
     >               MOD(x1+y1+z1,4).EQ.3)THEN ! Si sites
                     noepi=.FALSE.
                  END IF
                END IF
              END IF
              IF(noepi)NoEpiBulk=NoEpiBulk+1
            END IF
         END DO
        END DO
      END DO

!      write(*,*)NoEpiBulk
!      WRITE(FN+3,*) NoEpiBulk  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
!      WRITE(FN+3,*) ''

      write(*,*)NoEpiBulk
      WRITE(FN+3,'(1x,i8,1x,a)') NoEpiBulk,'angstroem'  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
      WRITE(FN+3,'(1x,a,3(1x,e15.8))')'periodic',alat(1),alat(2),alat(3)

      DO z=0,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF(IfOcc.EQ.1.AND.Coo.EQ.4)THEN
              noepi=.TRUE.
              IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.
     >           AND.MOD(z,3).EQ.0)THEN
                x1=x/3
                y1=y/3
                z1=z/3
                IF((MOD(x1,2).EQ.MOD(y1,2)).AND.
     >             (MOD(x1,2).EQ.MOD(z1,2)))THEN
                  IF(MOD(x1+y1+z1,4).EQ.0.OR.
     >               MOD(x1+y1+z1,4).EQ.3)THEN ! Si sites
                     noepi=.FALSE.
                  END IF
                END IF
              END IF
              IF(noepi)THEN
                IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                IF(IfSi.EQ.1)WRITE(FN+3,103)sf*x,sf*y,sf*z
                IF(IfC.EQ.1)WRITE(FN+3,101)sf*x,sf*y,sf*z
              END IF
            END IF
         END DO
        END DO
      END DO

!      NvoidMerged=NumVoids
!      WRITE(FN,*) NumAtTotal+NvoidMerged  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
!      WRITE(FN,*) ''

      NvoidMerged=NumVoids
      WRITE(FN,'(1x,i8,1x,a)') NumAtTotal+NvoidMerged,'angstroem'  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
      WRITE(FN,'(1x,a,3(1x,e15.8))')'periodic',alat(1),alat(2),alat(3)

      DO z=0,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
            IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF(IfSi.EQ.1.AND.Coo.LE.3)WRITE(FN,100)sf*x,sf*y,sf*z !, rad
            IF(IfC.EQ.1.AND.Coo.LE.3)WRITE(FN,102)sf*x,sf*y,sf*z !, rad
         END DO
        END DO
      END DO
      DO i=1,NvoidMerged
        x = ListVoid (i) % AtomXYZ(1)
        y = ListVoid (i) % AtomXYZ(2)
        z = ListVoid (i) % AtomXYZ(3)
        WRITE(FN,104)sf*x,sf*y,sf*z !, rad
      END DO


!      WRITE(FN+4,*) NumVoids  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
!      WRITE(FN+4,*) ''

      WRITE(FN+4,'(1x,i8,1x,a)') NumVoids,'angstroem'  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
      WRITE(FN+4,'(1x,a,3(1x,e15.8))')'periodic',alat(1),alat(2),alat(3)

      DO i=1,NumVoids
        x = ListVoid (i) % AtomXYZ(1)
        y = ListVoid (i) % AtomXYZ(2)
        z = ListVoid (i) % AtomXYZ(3)
        WRITE(FN+4,104)sf*x,sf*y,sf*z !, rad
      END DO

      
      WRITE(FN+5,'(1x,i8,1x,a)') NumAtWrong,'angstroem'  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
      WRITE(FN+5,'(1x,a,3(1x,e15.8))')'periodic',alat(1),alat(2),alat(3)
      DO z=1,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF(IfOcc.EQ.1.AND.Coo.LE.3)THEN
                NumAtTotal=NumAtTotal+1
            ELSE IF (IfOcc.EQ.1.AND.Coo.EQ.4)THEN
                indAtom=LattInd(x,y,z)
                IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                NextN=ListAtom(indAtom) % NextNXYZ
                NCoorCorr=0
                IF(IfSi.EQ.1)THEN
                    DO i=1,4
                        NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosC,LenC)
                    END DO
                ELSE IF(IfC.EQ.1)THEN
                    DO i=1,4
                        NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosSi,LenSi)
                    END DO
                END IF
                IF(NCoorCorr.EQ.3)THEN 
                    IF(IfSi.EQ.1)WRITE(FN+5,100)sf*x,sf*y,sf*z
                    IF(IfC.EQ.1)WRITE(FN+5,102)sf*x,sf*y,sf*z
                END IF 
            END IF
         END DO
        END DO
      END DO

      DO i=1,NumIsWrong
        x=ListAtBuffer(i) % AtomXYZ(1)
        y=ListAtBuffer(i) % AtomXYZ(2)
        z=ListAtBuffer(i) % AtomXYZ(3)
        ind_atotype=ListAtBuffer(NumIsWrong) % Ind_Atype 
        IF(ind_atotype.EQ.1)THEN
            LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
        ELSE
            LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
        END IF
      END DO 
      
      CLOSE(FN)
      CLOSE(FN+3)
      CLOSE(FN+4)
      CLOSE(FN+5)


      END SUBROUTINE WriteMolMolXYZFile
******|****************************************************************




***********************************************************************
      SUBROUTINE GetOutputFileName(Time, Iter, FileNameBase)

      IMPLICIT    NONE
      INTEGER     Iter,T
      REAL(8)     Time
      CHARACTER(Len=10)   FileNameBase
      CHARACTER(LEN=1)   C1,C1A,C1B,C1C,C1D,C1E,C1F,C1G,C1H

      Iter=Iter*10

      IF (Iter.GE.1000000000)
     &  STOP 'GetOutputFileName:: Iter can not be represented'
C    Character representation offset: 48 !!!
      T = ITER
      IF (T.GE.100000000) THEN
        C1 = CHAR(INT(T/100000000) + 48)
        T = MOD(T,100000000)
      ELSE
        C1 = '0'
      ENDIF
      IF (T.GE.10000000) THEN
        C1A = CHAR(INT(T/10000000) + 48)
        T = MOD(T,10000000)
      ELSE
        C1A = '0'
      ENDIF
      IF (T.GE.1000000) THEN
        C1B = CHAR(INT(T/1000000) + 48)
        T = MOD(T,1000000)
      ELSE
        C1B = '0'
      ENDIF
      IF (T.GE.100000) THEN
        C1C = CHAR(INT(T/100000) + 48)
        T = MOD(T,100000)
      ELSE
        C1C = '0'
      ENDIF
      IF (T.GE.10000) THEN
        C1D = CHAR(INT(T/10000) + 48)
        T = MOD(T,10000)
      ELSE
        C1D = '0'
      ENDIF
      IF (T.GE.1000) THEN
        C1E = CHAR(INT(T/1000) + 48)
        T = MOD(T,1000)
      ELSE
        C1E = '0'
      ENDIF
      IF (T.GE.100) THEN
        C1F = CHAR(INT(T/100) + 48)
        T = MOD(T,100)
      ELSE
        C1F = '0'
      ENDIF
      C1G = CHAR(INT(T/10) + 48)
      C1H = CHAR(MOD(T,10) + 48)

      FileNameBase = 'I'//C1//C1A//C1B//C1C//C1D//C1E//C1F//C1G//C1H

      END SUBROUTINE GetOutputFileName




******|****************************************************************

!      RSite(1)=3.0
!      RSite(2)=3.0
!      RSite(3)=3.0
!      RSiteC(1)=8.0
!      RSiteC(2)=2.0
!      RSiteC(3)=2.0
!      ax(1)=RSiteC(1)-RSite(1)
!      ax(2)=RSiteC(2)-RSite(2)
!      ax(3)=RSiteC(3)-RSite(3)
!      snorm=DSQRT(ax(1)**2+ax(2)**2+ax(3)**2)
!      ax(1)=ax(1)/snorm
!      ax(2)=ax(2)/snorm
!      ax(3)=ax(3)/snorm
!      DO irot=1,5
!         Cth=DCOS(DFLOAT(irot)*pi/3.)
!         Sth=DSIN(DFLOAT(irot)*pi/3.)
!         DO i=1,3
!            Rot(i,i)=ax(i)**2*(1-Cth)+Cth
!         END DO
!         Rot(2,1)=ax(1)*ax(2)*(1-Cth)-ax(3)*Sth
!         Rot(1,2)=ax(1)*ax(2)*(1-Cth)+ax(3)*Sth
!         Rot(3,1)=ax(1)*ax(3)*(1-Cth)+ax(2)*Sth
!         Rot(1,3)=ax(1)*ax(3)*(1-Cth)-ax(2)*Sth
!         Rot(3,2)=ax(2)*ax(3)*(1-Cth)-ax(1)*Sth
!         Rot(2,3)=ax(2)*ax(3)*(1-Cth)+ax(1)*Sth
!         RoldSite(1)=1.
!         RoldSite(2)=-5.
!         RoldSite(3)=1.
!         DO k=1,3
!            RNewSite(k)=0.0
!            DO l=1,3
!               RNewSite(k)=RNewSite(k)+Rot(l,k)*RoldSite(l)
!            END DO
!         END DO
!         NewSite=NINT(RNewSite)
!         write(*,*)NewSite,RNewSite
!      END DO

!      NumAtoms=0
!      DO i=1,8
!         DO j=1,3
!            in=NComplDia(j,i)
!            write(*,*)i,in,JumpDia(1:3,in)
!            DO k=1,3
!               NewSite(k)=0
!               DO l=1,3
!                  NewSite(k)=NewSite(k)+IMatRotat(l,k,i)*JumpDia(l,in)/3
!               END DO
!            END DO
!            write(*,*)NewSite
!         END DO
!         write(*,*)' '
!      END DO
