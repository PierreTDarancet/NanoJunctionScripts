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
   MODULE control_variable_module
!*********************************************
   USE kinds,      ONLY : dbl
   USE parameters, ONLY : nstrx
   IMPLICIT NONE
   PRIVATE 
   SAVE
!
! Contains GLOBAL CONTROL variables
! 

!
   CHARACTER(nstrx) :: calculation_mode
!  method for recusrsion states' calculation 
! 
   CHARACTER(nstrx) :: variation_mode
!  how calculate the variation between 2 iterations
! 
!
   INTEGER :: max_iter
! max number of iteration
!
   INTEGER :: min_iter
! min number of iteration

!
   REAL(dbl) :: conv_criteria
! convergency criteria

!
   REAL(dbl) :: numerical_cut_off
! convergency criteria

!
   LOGICAL :: use_second
!
   LOGICAL :: use_third
!
   LOGICAL :: print_hamiltonian
!
   LOGICAL :: print_center
!
   LOGICAL :: print_state
!
   LOGICAL :: print_variation
!
   LOGICAL :: print_A_sum
!
   LOGICAL :: print_B_sum
!
   LOGICAL :: print_matrix
!
   LOGICAL :: print_overlap
!
   LOGICAL :: debug_mode
!
   LOGICAL :: print_orbital
!
   LOGICAL :: print_distance
!

!
! end delcarations
!
!
!
   PUBLIC ::  use_second
!
   PUBLIC ::  use_third
!
   PUBLIC ::  calculation_mode
!
   PUBLIC ::  variation_mode
!
   PUBLIC ::  max_iter
!
   PUBLIC ::  min_iter
!
   PUBLIC :: conv_criteria
!
   PUBLIC :: print_orbital
!
   PUBLIC :: print_distance
!
   PUBLIC :: print_hamiltonian
!
   PUBLIC :: print_center
!
   PUBLIC :: print_state
!
   PUBLIC :: print_variation
!
   PUBLIC :: print_A_sum
!
   PUBLIC :: print_B_sum
!
   PUBLIC :: print_matrix
!
   PUBLIC :: print_overlap
!
   PUBLIC :: debug_mode
!
   PUBLIC :: numerical_cut_off
!
END MODULE control_variable_module

