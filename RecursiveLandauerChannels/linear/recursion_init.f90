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
   USE constants,            ONLY : ZERO, CZERO, CONE
   USE io_global_module,     ONLY : stdin, stdout
   USE io_module,     ONLY :  aux_unit
   USE iotk_module
   USE control_variable_module,     ONLY : datafile_H
   USE dim_variable_module, ONLY : dimwan, dim_recursion, dim_subspace
   USE subspace_variable_module, ONLY : x_limit, perform_recursion_x_greater, perform_recursion_x_lesser
   USE recursion_module, ONLY : recursion_allocate, initial_value, init_wan_num, total_hamiltonian, &
                                build_ham_and_coor, select_wf, reduce_hamiltonian, init_first_recursion_state
   USE timing_module



   IMPLICIT NONE
   !
   ! local variables
   !
   CHARACTER(14) :: subname="recursion_init"
   COMPLEX(dbl), ALLOCATABLE :: aux(:,:)
   REAL(dbl), ALLOCATABLE :: wannier_coordinates_aux(:,:)
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

   ALLOCATE( aux(dimwan,dimwan), STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'allocating aux',ABS(ierr))
   ALLOCATE( wannier_coordinates_aux(3,dimwan), STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'allocating wannier_coordinates_aux',ABS(ierr))

   aux(:,:)=CZERO
   wannier_coordinates_aux(:,:)=ZERO
        !return total hamiltonian components and wannier center coordinates, contained in the datafile_H
   CALL recursion_read_matrix( datafile_H, aux, dimwan, dimwan, wannier_coordinates_aux, 3, dimwan)
   CALL build_ham_and_coor(aux, wannier_coordinates_aux)

    ! geometry finder(x_limit)
    CALL select_wf(x_limit, perform_recursion_x_lesser,perform_recursion_x_greater)

   !  total_hamiltonian(:,:) => recursion_hamiltonian
     ! select hamiltonian components
    CALL reduce_hamiltonian()

    CALL init_wan_num()
    CALL init_first_recursion_state()

   !
   ! local cleaning
   !
   DEALLOCATE( aux, STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'deallocating aux',ABS(ierr))
   DEALLOCATE( wannier_coordinates_aux, STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'deallocating aux',ABS(ierr))

   CALL timing( 'recursion_init', OPR='stop' )


END SUBROUTINE recursion_init
