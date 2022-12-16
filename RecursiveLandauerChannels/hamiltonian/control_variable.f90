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

   CHARACTER(nstrx) :: datafile_H
!
   LOGICAL :: use_second 
!
   LOGICAL :: use_third 
!
   LOGICAL :: print_distance 
!
   LOGICAL :: print_orbital 
!
   LOGICAL :: print_hamiltonian
!

! end delcarations
!
!
   PUBLIC ::  datafile_H
!
   PUBLIC :: use_second
!
   PUBLIC :: use_third
!
   PUBLIC :: print_distance 
!
   PUBLIC :: print_orbital 
!
   PUBLIC :: print_hamiltonian

!

END MODULE control_variable_module

