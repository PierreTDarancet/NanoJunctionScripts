!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE dim_variable_module
!*********************************************

   IMPLICIT NONE
   PRIVATE 
   SAVE
!
! Contains GLOBAL DIM variables
! 

! 
! 
!
   INTEGER :: dim_subspace
! dimension of the initial state for the recursion

!
   INTEGER :: dim_recursion
! dimension of the hamiltonian used for the recursion

!
   INTEGER :: dimwan
! dimension of the total hamiltonian, given in inputs

!
! end delcarations
!
!
   PUBLIC ::  dim_subspace
   PUBLIC ::  dim_recursion
   PUBLIC ::  dimwan

END MODULE dim_variable_module

