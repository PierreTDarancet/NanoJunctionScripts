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
      USE orbital_module,           ONLY : orbital_allocate, orbital_init, orbital_alloc
      USE input_base_module,        ONLY : read_cards, card_orb

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

      orbital_alloc=.FALSE.
      CALL orbital_allocate()
      CALL orbital_init()


      CALL read_cards(stdin)


   END SUBROUTINE input_manager


!**********************************************************
   SUBROUTINE setup_control()
   !**********************************************************
      USE control_variable_module,     ONLY :    min_iter, max_iter, conv_criteria, &
                                                 variation_mode, calculation_mode, &
                                                 print_hamiltonian, print_center, &
                                                 print_state,  print_variation, &
                                                 print_A_sum, print_B_sum, &
                                                 print_matrix, print_overlap, &
                                                 print_orbital, print_distance, &
                                                 use_second, use_third, &
                                                 debug_mode, numerical_cut_off


      USE input_parameters_module, ONLY :        min_iter_ => min_iter,     &
                                                 max_iter_ => max_iter,     &
                                                 use_second_    => use_second, &
                                                 use_third_     => use_third, &
                                                 conv_criteria_ => conv_criteria, &
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
                                                 print_orbital_ => print_orbital, &
                                                 print_distance_ => print_distance, &
                                                 debug_mode_ => debug_mode, &
                                                 numerical_cut_off_ => numerical_cut_off

      IMPLICIT NONE


     use_second = use_second_
     use_third  = use_third_
     max_iter = max_iter_
     min_iter = min_iter_
     conv_criteria = conv_criteria_
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
     print_orbital = print_orbital_
     print_distance = print_distance_
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
      USE dim_variable_module,    ONLY :  limit_0, limit_1, limit_2, &
                                          limit_3, n_orb, &
                                          dim_recursion, dim_subspace, &
                                          nb_max_first, &
                                          nb_max_second, &
                                          nb_max_third

      USE input_parameters_module, ONLY :     n_orb_      => n_orb, &
                                              limit_0_    => limit_0, &
                                              limit_1_    => limit_1, &
                                              limit_2_    => limit_2, &
                                              limit_3_    => limit_3, &
                                              dim_recursion_ => dim_recursion, &
                                              dim_subspace_ => dim_subspace, &
                                              nb_max_first_  => nb_max_first, &
                                              nb_max_second_  => nb_max_second, &
                                              nb_max_third_  => nb_max_third


      IMPLICIT NONE
      CHARACTER(9) :: subname='setup_dim'

     limit_0    = limit_0_
     limit_1    = limit_1_
     limit_2    = limit_2_
     limit_3    = limit_3_
     n_orb      = n_orb_
     dim_recursion = dim_recursion_
     dim_subspace = dim_subspace_
     nb_max_first  = nb_max_first_
     nb_max_second  = nb_max_second_
     nb_max_third  = nb_max_third_



   END SUBROUTINE setup_dim


!**********************************************************
   SUBROUTINE setup_subspace()
   !**********************************************************
      USE subspace_variable_module,    ONLY :    x_limit, &
                                                 perform_recursion_x_greater, &
                                                 perform_recursion_x_lesser

      USE input_parameters_module, ONLY :        x_limit_ => x_limit, &
                                                 perform_recursion_x_greater_ => perform_recursion_x_greater, &
                                                 perform_recursion_x_lesser_  => perform_recursion_x_lesser

      IMPLICIT NONE

     x_limit = x_limit_
     perform_recursion_x_greater = perform_recursion_x_greater_
     perform_recursion_x_lesser  = perform_recursion_x_lesser_


   END SUBROUTINE setup_subspace


END MODULE input_module

