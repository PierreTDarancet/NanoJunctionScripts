!
! Copyright (C) 2009 Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET

!***********************************************
 MODULE Model_parameters_module
   !***********************************************
   IMPLICIT NONE 
   !
   ! Declarations
   !
   INTEGER                :: nbias        ! dimension of the energy grid
   REAL(dbl)              :: biasmin      !
   REAL(dbl)              :: biasmax      ! egrid extrema 
   REAL(dbl)              :: delta     ! i\delta for GFs
   INTEGER                :: ne             ! max dimension of the energy grid

   CHARACTER(nstrx)       :: Field_level_alignment_variation

   !
   ! Variable privacy
   !

   PUBLIC                :: nbias        ! dimension of the energy grid
   PUBLIC                :: biasmin      !
   PUBLIC                :: biasmax      ! egrid extrema 
   PUBLIC                :: delta     ! i\delta for GFs
   PUBLIC                :: ne             ! max dimension of the energy grid


   !
   ! Definitions 
   !
   nbias=100
   biasmin=-4.00 !!!!!!!!!!!!!!! biasmin = - biasmax
   biasmax=4.00 !!!!!!!!!!!!!!!!! otherwise RR is false
   delta=0.00001 !idelta for gf
   ne= 1001



END MODULE
