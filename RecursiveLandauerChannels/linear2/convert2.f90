!
! Copyright (C) 2007 Institut N?el Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!=====================================================
   PROGRAM convert2
   !=====================================================



    !   USE parameters, ONLY : nstrx
       USE convert2_input_module, ONLY : convert2_input_manager
       USE convert2_module, ONLY : matrix2_allocate, &
                                  read2_matrix, &
                                  build2_matrix, &
                                  build2_hamiltonian, &
                                  print2_want, &
                                  matrix2_deallocate, &
                                  convert2_summary_input

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
      CALL convert2_input_manager()
      CALL convert2_summary_input()
      CALL matrix2_allocate()
      CALL read2_matrix()
      CALL build2_matrix()
      CALL build2_hamiltonian()
      CALL print2_want()
      CALL matrix2_deallocate()

END PROGRAM convert2

