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
   PROGRAM convert4
   !=====================================================



    !   USE parameters, ONLY : nstrx
       USE convert4_input_module, ONLY : convert4_input_manager
       USE convert4_module, ONLY : matrix4_allocate, &
                                  read4_matrix, &
                                  build4_matrix, &
                                  build4_hamiltonian, &
                                  print4_want, &
                                  matrix4_deallocate, &
                                  convert4_summary_input

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
      CALL convert4_input_manager()
      CALL convert4_summary_input()
      CALL matrix4_allocate()
      CALL read4_matrix()
      CALL build4_matrix()
      CALL build4_hamiltonian()
      CALL print4_want()
      CALL matrix4_deallocate()

END PROGRAM convert4

