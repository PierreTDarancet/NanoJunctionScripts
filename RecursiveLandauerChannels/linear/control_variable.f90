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
   CHARACTER(nstrx) :: datafile_H
!  datafile for hamiltonian

!
   INTEGER :: max_iter
! max number of iteration

!
   REAL(dbl) :: conv_criteria
! convergency criteria

!
! end delcarations
!
!
   PUBLIC ::  datafile_H
!
   PUBLIC ::  max_iter
!
   PUBLIC :: conv_criteria

END MODULE control_variable_module

