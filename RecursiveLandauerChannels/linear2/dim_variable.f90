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
   INTEGER :: dimwan_total
!
   INTEGER :: dim_L
   INTEGER :: dim_R
! dimension of the l/R hamiltonian
!
!   INTEGER :: dim_extension
!
   INTEGER :: nb_replica_L
   INTEGER :: nb_replica_R
!   INTEGER :: nb_replica_extension
!  number of replica 



!
! end delcarations
!
!
   PUBLIC ::  dim_subspace
   PUBLIC ::  dim_recursion
   PUBLIC ::  dimwan
   PUBLIC ::  dimwan_total

   PUBLIC ::  dim_L
   PUBLIC ::  dim_R
!   PUBLIC ::  dim_extension
!
   PUBLIC ::  nb_replica_R
   PUBLIC ::  nb_replica_L
!   PUBLIC ::  nb_replica_extension
!

END MODULE dim_variable_module

