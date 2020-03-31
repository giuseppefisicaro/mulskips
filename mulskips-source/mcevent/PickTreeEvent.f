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
