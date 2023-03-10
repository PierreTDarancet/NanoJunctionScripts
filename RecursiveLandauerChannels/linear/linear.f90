!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!=====================================================
   PROGRAM linear
   !=====================================================
 ! This program allows to linearize an arbitrary 3D real space Hamiltonian thanks to the matricial
 ! recursion method.
 ! inputs are an Hamiltonian given by hamiltonian.x


       USE kinds
       USE parameters, ONLY : nstrx
       USE input_module, ONLY : input_manager
       USE output_module, ONLY : output_manager
       USE io_global_module, ONLY : stdout
       USE summary_module, ONLY : summary_input
       USE timing_module, ONLY : timing, timing_upto_now, timing_overview, global_list
       USE version_module, ONLY : version_number

!
!
      IMPLICIT NONE

!       INTEGER :: ik, m, n, ierr

!
!--------------------------------------------
! ... Startup
!--------------------------------------------
!
      CALL startup(version_number,'wannier')

      !
      ! ... Read input parameter
      !
      CALL input_manager()
      ! ... Check input
      !      CALL test_input()
      ! ...Summarize the inputs
      CALL summary_input( stdout )
      !      CALL timing_upto_now(stdout)
      ! ... main part
      CALL recursion()
      !... output
      CALL output_manager()
      !      cleanup
      CALL cleanup()

END PROGRAM linear

