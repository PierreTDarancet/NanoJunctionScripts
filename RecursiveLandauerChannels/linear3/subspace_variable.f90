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

END MODULE subspace_variable_module

