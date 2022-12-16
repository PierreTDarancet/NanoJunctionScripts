!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*******************************************************************
   SUBROUTINE recursion_init
   !*******************************************************************

   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, EPS_m4
   USE io_global_module,     ONLY : stdin, stdout
   USE io_module,     ONLY :  aux_unit
   USE iotk_module
   USE subspace_variable_module, ONLY : x_limit, perform_recursion_x_greater, perform_recursion_x_lesser

   USE recursion_module, ONLY : recursion_allocate, initial_value, init_wan_num,  &
                                select_wf, init_first_recursion_state
   USE timing_module
   USE distance_module,      ONLY : distance_allocate, distance_alloc, &
                                    init_metric



   IMPLICIT NONE
   !
   ! local variables
   !
   CHARACTER(14) :: subname="recursion_init"
   INTEGER       :: i, ierr

   !
   ! end of declarations
   !

!
!----------------------------------------
! main Body
!----------------------------------------
!
    CALL timing( 'recursion_init', OPR='start' )

 !
   ! allocations
   !
    CALL recursion_allocate()
   !
    CALL initial_value()


    ! geometry finder(x_limit)
    CALL select_wf(x_limit, perform_recursion_x_lesser,perform_recursion_x_greater)


    distance_alloc=.FALSE.

    CALL distance_allocate()
    !
    CALL init_metric() 

   !  total_hamiltonian(:,:) => recursion_hamiltonian
     ! select hamiltonian components
   ! CALL reduce_hamiltonian()
    !
    CALL init_wan_num()
    !
    CALL init_first_recursion_state()
    !
    !
    ! local cleaning
    !

    CALL timing( 'recursion_init', OPR='stop' )


END SUBROUTINE recursion_init
