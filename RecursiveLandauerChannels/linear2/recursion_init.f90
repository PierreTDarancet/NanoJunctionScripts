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
   USE control_variable_module,     ONLY : datafile_H, datafile_R, datafile_L, use_extension_L, use_extension_R
   USE dim_variable_module, ONLY : dimwan, dim_recursion, dim_subspace, dimwan_total, dim_L, dim_R, nb_replica_L, nb_replica_R
   USE subspace_variable_module, ONLY : x_limit, perform_recursion_x_greater, perform_recursion_x_lesser, &
                                        cell_size_C, cell_size_L, cell_size_R
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

   ALLOCATE( aux(dimwan_total,dimwan_total), STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'allocating aux',ABS(ierr))
   ALLOCATE( wannier_coordinates_aux(3,dimwan_total), STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'allocating wannier_coordinates_aux',ABS(ierr))

   aux(:,:)=CZERO
   wannier_coordinates_aux(:,:)=ZERO
        !return total hamiltonian components and wannier center coordinates, contained in the datafile_H
   CALL recursion_read_matrix( datafile_H, 'H_C' , aux, dimwan_total, dimwan_total, wannier_coordinates_aux, 3, dimwan_total, dimwan, 0, 0, ZERO, ZERO)
   CALL build_ham_and_coor( aux,  wannier_coordinates_aux)
   IF (use_extension_L) THEN
         aux(:,:)=CZERO
         wannier_coordinates_aux(:,:)=ZERO
         ! attention, ne pas upgrader les wannier coordinates mais seulement l'Hamiltonien!
         CALL recursion_read_matrix( datafile_H, 'H_CL', aux, dimwan_total, dimwan_total, wannier_coordinates_aux, 3, dimwan_total, dimwan, dim_L, 0, cell_size_C, cell_size_L)
         CALL build_ham_and_coor(aux, wannier_coordinates_aux)

         aux(:,:)=CZERO
         wannier_coordinates_aux(:,:)=ZERO
         CALL recursion_read_matrix( datafile_L, 'H_L', aux, dimwan_total, dimwan_total, wannier_coordinates_aux, 3, dimwan_total,  dimwan, dim_L, nb_replica_L, cell_size_C, cell_size_L)
         CALL build_ham_and_coor(aux, wannier_coordinates_aux)
   ELSE IF (use_extension_R) THEN
         aux(:,:)=CZERO
         wannier_coordinates_aux(:,:)=ZERO
         ! attention, ne pas upgrader les wannier coordinates mais seulement l'Hamiltonien!
         CALL recursion_read_matrix( datafile_H, 'H_CR', aux, dimwan_total, dimwan_total, wannier_coordinates_aux, 3, dimwan_total,  dimwan, dim_R, 0, cell_size_C, cell_size_R)
         CALL build_ham_and_coor(aux, wannier_coordinates_aux)

         aux(:,:)=CZERO
         wannier_coordinates_aux(:,:)=ZERO
         CALL recursion_read_matrix( datafile_R, 'H_R', aux, dimwan_total, dimwan_total, wannier_coordinates_aux, 3, dimwan_total,  dimwan, dim_R, nb_replica_R, cell_size_C, cell_size_R)
         CALL build_ham_and_coor( aux, wannier_coordinates_aux)

   ENDIF

    ! geometry finder(x_limit)
    CALL select_wf(x_limit, perform_recursion_x_lesser,perform_recursion_x_greater)

   !  total_hamiltonian(:,:) => recursion_hamiltonian
     ! select hamiltonian components
    CALL reduce_hamiltonian()
     !
    CALL init_wan_num()
     !
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
