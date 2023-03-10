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
   SUBROUTINE recursion
   !=====================================================
 ! This routine applies the recursion method to a semi infinite hamiltonian
 ! Initial states are assumed to be orthogonal


       USE kinds
       USE constants, ONLY: CZERO
       USE parameters, ONLY : nstrx
       USE recursion_module, ONLY : first_iteration, recursion_loop

!
!
      IMPLICIT NONE

       INTEGER :: iter, ierr

!
!--------------------------------------------
! ... Startup
!--------------------------------------------
!

    ! Allocate Haimltonian, initialize hamiltonian and trial states
 CALL recursion_init()

  !first iteration
 CALL first_iteration()
  ! recursion
 CALL recursion_loop()

END SUBROUTINE recursion

