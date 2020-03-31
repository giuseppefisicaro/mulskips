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
