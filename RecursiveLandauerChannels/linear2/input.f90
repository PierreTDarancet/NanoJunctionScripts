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
   MODULE input_module
!********************************************
   USE kinds, ONLY : dbl
   USE io_module, ONLY : stdin, stdout
   USE constants, ONLY : ZERO
   USE iotk_module

   IMPLICIT NONE
   PRIVATE

!
! This module handles the reading of input data
!
! routines in this module:
! SUBROUTINE input_manager()
! SUBROUTINE setup_control()
! SUBROUTINE setup_dim()
! SUBROUTINE setup_subspace()

!

!
   PUBLIC :: input_manager

!

CONTAINS

!
! subroutines
!


!**********************************************************
   SUBROUTINE input_manager()
   !**********************************************************
      USE input_parameters_module,  ONLY : read_namelist_control

     IMPLICIT NONE
      INTEGER :: ierr
      CHARACTER(13) :: subname='input_manager'

      !
      ! reading and checking namelists
      !

      !!!!!!!!!!!!!!!!!!!!!!!!have to be improved : for the moment all 
      !!!!!!!!!!!!!!!!!!!!!!!!the variables are supposed to be read in the control 
      !!!!!!!!!!!!!!!!!!!!!!!!namelist even if they are next stored in different variables module

      CALL read_namelist_control(stdin)



      !
      ! scattering data in their own modules
      !
      CALL setup_control()
      CALL setup_io()
      CALL setup_dim()
      CALL setup_subspace()


   END SUBROUTINE input_manager


!**********************************************************
   SUBROUTINE setup_control()
   !**********************************************************
      USE control_variable_module,     ONLY :    datafile_H, datafile_L, datafile_R, &
                                                 min_iter, max_iter, conv_criteria, &
                                                 use_extension_R, use_extension_L, &
                                                 variation_mode, calculation_mode, &
                                                 print_hamiltonian, print_center, &
                                                 print_state,  print_variation, &
                                                 print_A_sum, print_B_sum, &
                                                 print_matrix, print_overlap,&
                                                 debug_mode, numerical_cut_off

      USE input_parameters_module, ONLY :        datafile_H_ => datafile_H, &
                                                 datafile_L_ => datafile_L, &
                                                 datafile_R_ => datafile_R,  &
                                                 min_iter_ => min_iter,     &
                                                 max_iter_ => max_iter,     &
                                                 conv_criteria_ => conv_criteria, &
                                                 use_extension_R_ => use_extension_R, &
                                                 use_extension_L_ => use_extension_L, &
                                                 variation_mode_ => variation_mode, &
                                                 calculation_mode_ => calculation_mode, &
                                                 print_hamiltonian_ => print_hamiltonian, &
                                                 print_center_ => print_center, &
                                                 print_state_ => print_state, &
                                                 print_variation_ => print_variation, &
                                                 print_A_sum_ => print_A_sum, &
                                                 print_B_sum_ => print_B_sum, &
                                                 print_matrix_ => print_matrix, &
                                                 print_overlap_ => print_overlap, &
                                                 debug_mode_ => debug_mode, &
                                                 numerical_cut_off_ => numerical_cut_off



      IMPLICIT NONE

     datafile_H = datafile_H_
     max_iter = max_iter_
     min_iter = min_iter_
     conv_criteria = conv_criteria_
     use_extension_R = use_extension_R_
     use_extension_L = use_extension_L_
     datafile_L = datafile_L_
     datafile_R = datafile_R_
     variation_mode = variation_mode_
     calculation_mode = calculation_mode_
     print_hamiltonian = print_hamiltonian_
     print_center = print_center_
     print_state = print_state_
     print_variation = print_variation_
     print_A_sum = print_A_sum_
     print_B_sum = print_B_sum_
     print_matrix = print_matrix_
     print_overlap = print_overlap_
     debug_mode = debug_mode_
     numerical_cut_off = numerical_cut_off_


   END SUBROUTINE setup_control

!**********************************************************
   SUBROUTINE setup_io()
   !**********************************************************
      USE io_module,                ONLY : prefix, postfix, work_dir, title
      USE input_parameters_module,  ONLY : prefix_    => prefix,       &
                                           postfix_   => postfix,      &
                                           work_dir_  => work_dir,     &
                                           title_     => title
      IMPLICIT NONE
      prefix   = prefix_
      postfix  = postfix_
      work_dir = work_dir_
      title    = title_
   END SUBROUTINE setup_io

!**********************************************************
   SUBROUTINE setup_dim()
   !**********************************************************
      USE dim_variable_module,    ONLY :  nb_replica_L, nb_replica_R, dimwan, dim_recursion, dim_subspace, dim_L, dim_R, dimwan_total
      USE input_parameters_module, ONLY :  dimwan_ => dimwan, &
                                           dim_recursion_ => dim_recursion, &
                                           dim_subspace_ => dim_subspace, &
                                           nb_replica_L_ => nb_replica_L, &
                                           nb_replica_R_ => nb_replica_R, &
                                           dim_L_ => dim_L, &
                                           dim_R_ => dim_R, &
                                           use_extension_R_ => use_extension_R, &
                                           use_extension_L_ => use_extension_L

      IMPLICIT NONE
      CHARACTER(9) :: subname='setup_dim'


     dimwan = dimwan_
     dim_recursion = dim_recursion_
     dim_subspace = dim_subspace_
     dim_L = dim_L_
     dim_R = dim_R_
     nb_replica_L = nb_replica_L_
     nb_replica_R = nb_replica_R_
     IF (use_extension_R_) THEN
        dimwan_total= dimwan + (nb_replica_R + 1) * dim_R
        IF ( dim_recursion >  dimwan_total ) &
               CALL errore(subname,'dim_recursion > dimwan_total',1)

     ELSE IF (use_extension_L_) THEN
        dimwan_total= dimwan + (nb_replica_L + 1) * dim_L
        IF ( dim_recursion >  dimwan_total ) &
               CALL errore(subname,'dim_recursion > dimwan_total',2)

     ELSE
        dimwan_total = dimwan
        IF ( dim_recursion >  dimwan_total ) &
               CALL errore(subname,'dim_recursion > dimwan_total',3)

     ENDIF

   END SUBROUTINE setup_dim


!**********************************************************
   SUBROUTINE setup_subspace()
   !**********************************************************
      USE subspace_variable_module,    ONLY :    x_limit, perform_recursion_x_greater, perform_recursion_x_lesser, &
                                                 cell_size_C, cell_size_L, cell_size_R
      USE input_parameters_module, ONLY :        x_limit_ => x_limit, &
                                                 perform_recursion_x_greater_ => perform_recursion_x_greater, &
                                                 perform_recursion_x_lesser_  => perform_recursion_x_lesser, &
                                                 cell_size_C_ => cell_size_C, &
                                                 cell_size_L_ => cell_size_L, &
                                                 cell_size_R_ => cell_size_R
      IMPLICIT NONE

     x_limit = x_limit_
     perform_recursion_x_greater = perform_recursion_x_greater_
     perform_recursion_x_lesser  = perform_recursion_x_lesser_
     cell_size_C = cell_size_C_
     cell_size_L = cell_size_L_
     cell_size_R = cell_size_R_

   END SUBROUTINE setup_subspace


END MODULE input_module

