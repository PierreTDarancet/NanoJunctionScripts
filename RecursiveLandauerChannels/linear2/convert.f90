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
   PROGRAM convert
   !=====================================================



    !   USE parameters, ONLY : nstrx
       USE convert_input_module, ONLY : convert_input_manager
       USE convert_module, ONLY : matrix_allocate, &
                                  read_matrix, &
                                  build_matrix, &
                                  build_hamiltonian, &
                                  print_want, &
                                  matrix_deallocate, &
                                  convert_summary_input

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
      CALL convert_input_manager()
      CALL convert_summary_input()
      CALL matrix_allocate()
      CALL read_matrix()
      CALL build_matrix()
      CALL build_hamiltonian()
      CALL print_want()
      CALL matrix_deallocate()

END PROGRAM convert

