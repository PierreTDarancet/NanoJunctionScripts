!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
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
      USE recursion_module, ONLY : print_hamiltonian, print_matrix, &
                                   print_states, print_wan_center, print_variation, &
                                   print_A_sum, print_B_sum, print_state_overlap
      USE distance_module,  ONLY : print_distance
      USE orbital_module,   ONLY : print_orbital
      USE control_variable_module,  ONLY : log_print_hamiltonian => print_hamiltonian, &
                                           log_print_matrix => print_matrix, &
                                           log_print_state => print_state, &
                                           log_print_center => print_center, &
                                           log_print_variation => print_variation, &
                                           log_print_A_sum => print_A_sum, &
                                           log_print_B_sum => print_B_sum, &
                                           log_print_overlap => print_overlap, &
                                           log_print_orbital => print_orbital, &
                                           log_print_distance => print_distance

      INTEGER :: ierr
      CHARACTER(13) :: subname='output_manager'
!

  !    INTEGER :: ik, m, n, ierr
      !
  !    CHARACTER( LEN=nstrx )  :: filename
      !
      ! ... end of declarations
      !

   !CALL timing( 'output_manager', OPR='start' )

!
!--------------------------------------------
! ... Startup
!--------------------------------------------
!
!
!--------------------------------------------
!...  Init Wannier Functions localization procedure
!--------------------------------------------
! 

!
    IF (log_print_hamiltonian)  CALL print_hamiltonian()
!
    IF (log_print_matrix)       CALL print_matrix()
!
    IF (log_print_state)        CALL print_states()
!
    IF (log_print_center)       CALL print_wan_center()
!
    IF (log_print_variation)    CALL print_variation()
!
    IF (log_print_A_sum)        CALL print_A_sum()
!
    IF (log_print_B_sum)        CALL print_B_sum()
!
    IF (log_print_overlap)      CALL print_state_overlap()
!
    IF (log_print_distance)     CALL print_distance()
!
    IF (log_print_orbital)      CALL print_orbital()

!      CALL timing_upto_now(stdout)
   END SUBROUTINE output_manager
 END MODULE output_module

