!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!********************************************
   MODULE output_module
!********************************************
   USE kinds, ONLY : dbl
   USE io_module, ONLY : stdin, stdout
   USE constants, ONLY : ZERO
   USE iotk_module


   IMPLICIT NONE
   PRIVATE
! routines in this module:
! SUBROUTINE output_manager()


   PUBLIC :: output_manager


CONTAINS

!
! subroutines
!


!**********************************************************
   SUBROUTINE output_manager()
   !**********************************************************
      USE hamiltonian_module,       ONLY : print_hamiltonian
      USE orbital_module,           ONLY : print_orbital
      USE distance_module,          ONLY : print_distance
      USE control_variable_module,  ONLY : log_print_hamiltonian => print_hamiltonian, &
                                           log_print_orbital     => print_orbital, &
                                           log_print_distance    => print_distance

      INTEGER :: ierr
      CHARACTER(13) :: subname='output_manager'
!


   !CALL timing( 'output_manager', OPR='start' )

!
!--------------------------------------------
! ... Startup
!--------------------------------------------
!
!
    IF (log_print_hamiltonian)  CALL print_hamiltonian()
!
    IF (log_print_orbital)      CALL print_orbital()
!
    IF (log_print_distance)     CALL print_distance()
!
!      CALL timing_upto_now(stdout)
   END SUBROUTINE output_manager

 END MODULE output_module

