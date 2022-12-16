!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
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

      orbital_alloc=.FALSE.
      CALL orbital_allocate()
      CALL orbital_init()


      CALL read_cards(stdin)


   END SUBROUTINE input_manager


!**********************************************************
   SUBROUTINE setup_control()
   !**********************************************************
      USE control_variable_module,     ONLY :    datafile_H, &
                                                 use_second, use_third, &
                                                 print_orbital, &
                                                 print_hamiltonian, &
                                                 print_distance
      USE input_parameters_module, ONLY :        datafile_H_        => datafile_H, &
                                                 use_second_        => use_second, &
                                                 use_third_         => use_third, &
                                                 print_orbital_     => print_orbital, &
                                                 print_hamiltonian_ => print_hamiltonian, &
                                                 print_distance_    => print_distance


      IMPLICIT NONE



     use_second = use_second_
     use_third  = use_third_
     datafile_H = datafile_H_
     print_orbital     = print_orbital_
     print_hamiltonian = print_hamiltonian_
     print_distance    = print_distance_


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
      USE dim_variable_module,     ONLY :   cell_size, limit_0, limit_1, limit_2, &
                                            limit_3, n_orb
      USE input_parameters_module, ONLY :   n_orb_      => n_orb, &
                                            limit_0_    => limit_0, &
                                            cell_size_  => cell_size, &
                                            limit_1_    => limit_1, &
                                            limit_2_    => limit_2, &
                                            limit_3_    => limit_3

      IMPLICIT NONE

     cell_size  = cell_size_
     limit_0    = limit_0_
     limit_1    = limit_1_
     limit_2    = limit_2_
     limit_3    = limit_3_
     n_orb      = n_orb_

   END SUBROUTINE setup_dim
END MODULE input_module

