!
! Copyright (C) 2006 LEPES-CNRS/2007 Institut Néel Grenoble
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
   USE kinds

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
   INTEGER :: n_orb
! dimension of the total hamiltonian, given in inputs
!
   INTEGER :: nb_max_first
! max number of first neighbor
   INTEGER :: nb_max_second
!
   INTEGER :: nb_max_third
!
   REAL(dbl) :: limit_0
!
   REAL(dbl) :: limit_1
!
   REAL(dbl) :: limit_2
!
   REAL(dbl) :: limit_3
!

!
! end delcarations
!
!
   PUBLIC ::  dim_subspace
!
   PUBLIC ::  dim_recursion
!
   PUBLIC ::  n_orb
!
   PUBLIC ::  limit_0
!
   PUBLIC ::  limit_1
!
   PUBLIC ::  limit_2
!
   PUBLIC ::  limit_3
!
   PUBLIC ::  nb_max_first
!
   PUBLIC ::  nb_max_second
!
   PUBLIC ::  nb_max_third
!


END MODULE dim_variable_module

