!
! Copyright (C) 2005 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*********************************************
   MODULE subspace_variable_module
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
   REAL(dbl) :: x_limit
! begin of the recursion chain
!
   REAL(dbl) :: cell_size_C
!
   REAL(dbl) :: cell_size_R
!
   REAL(dbl) :: cell_size_L
!
   LOGICAL :: perform_recursion_x_greater
!
   LOGICAL :: perform_recursion_x_lesser
! take hamiltonian data with x >/< x0 to generate the recursion hamiltonian 
!
! end delcarations
!
!
   PUBLIC ::  x_limit
!
   PUBLIC :: perform_recursion_x_greater
!
   PUBLIC :: perform_recursion_x_lesser
!
   PUBLIC :: cell_size_C
!
   PUBLIC :: cell_size_R
!
   PUBLIC :: cell_size_L

END MODULE subspace_variable_module

