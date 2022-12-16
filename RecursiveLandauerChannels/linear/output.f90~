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
      USE recursion_module,  ONLY : print_hamiltonian, print_matrix, &
                                   print_states, print_wan_center

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

      CALL print_hamiltonian()
      CALL print_matrix()
      CALL print_states()
      CALL print_wan_center()

!      CALL timing_upto_now(stdout)
   END SUBROUTINE output_manager
 END MODULE output_module

