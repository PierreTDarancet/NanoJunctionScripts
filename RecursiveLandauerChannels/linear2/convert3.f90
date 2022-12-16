!
! Copyright (C) 2007 Institut Néel Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!=====================================================
   PROGRAM convert3
   !=====================================================



    !   USE parameters, ONLY : nstrx
       USE convert3_input_module, ONLY : convert3_input_manager
       USE convert3_module, ONLY : matrix3_allocate, &
                                  read3_matrix, &
                                  build3_matrix, &
                                  build3_hamiltonian, &
                                  print3_want, &
                                  matrix3_deallocate, &
                                  convert3_summary_input

       !USE io_global_module, ONLY : stdout
       !USE summary_module, ONLY : convert_summary_input
      ! USE timing_module, ONLY : timing, timing_upto_now, timing_overview, global_list
    !   USE version_module, ONLY : version_number

!
!
      IMPLICIT NONE

!       INTEGER :: ik, m, n, ierr

!
!--------------------------------------------
! ... Startup
!--------------------------------------------
!
     ! CALL startup(version_number,'wannier')
      !
      ! ... Read input parameter
      !
      CALL convert3_input_manager()
      CALL convert3_summary_input()
      CALL matrix3_allocate()
      CALL read3_matrix()
      CALL build3_matrix()
      CALL build3_hamiltonian()
      CALL print3_want()
      CALL matrix3_deallocate()

END PROGRAM convert3

