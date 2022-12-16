!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE dim_variable_module
!*********************************************
   USE kinds,      ONLY : dbl

   IMPLICIT NONE
   PRIVATE 
   SAVE
!
! Contains GLOBAL DIM variables
! 
  INTEGER :: n_orb
! dimension of the initial state for the recursion
!
   REAL(dbl) :: cell_size
! useful to shift the cell

   REAL(dbl) :: limit_0
! lower limit
   REAL(dbl) :: limit_1
! limit for first neightbor
   REAL(dbl) :: limit_2
! limit for second neightbor
   REAL(dbl) :: limit_3
! limit for third neightbor


! 
! 

!
! end delcarations
!
!
   PUBLIC :: cell_size
!
   PUBLIC :: n_orb
!
   PUBLIC :: limit_0
!
   PUBLIC :: limit_1
!
   PUBLIC :: limit_2
!
   PUBLIC :: limit_3
!

END MODULE dim_variable_module

