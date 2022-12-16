!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!=====================================================
   PROGRAM hamiltonian
   !=====================================================
       USE kinds
       USE parameters, ONLY : nstrx
       USE input_module, ONLY : input_manager
       USE data_module, ONLY : orb_init
       USE output_module, ONLY : output_manager
       USE io_module, ONLY : stdout
       USE summary_module, ONLY : summary_input
!       USE timing_module, ONLY : timing, timing_upto_now, timing_overview, global_list
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

      CALL orb_init()
      !
      ! ...Summarize the inputs
      CALL summary_input( stdout )
      !      CALL timing_upto_now(stdout)

      ! ... main part
      CALL construct_hamiltonian()
      !... output
      CALL output_manager()
      !      cleanup
      CALL cleanup()

END PROGRAM hamiltonian

