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


   PUBLIC :: input_manager


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
      USE control_variable_module,     ONLY :    datafile_H, max_iter, conv_criteria
      USE input_parameters_module, ONLY :        datafile_H_ => datafile_H, &
                                                 max_iter_ => max_iter,     &
                                                 conv_criteria_ => conv_criteria

      IMPLICIT NONE

     datafile_H = datafile_H_
     max_iter = max_iter_
     conv_criteria = conv_criteria_

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
      USE dim_variable_module,         ONLY :    dimwan, dim_recursion, dim_subspace
      USE input_parameters_module, ONLY :        dimwan_ => dimwan, &
                                                 dim_recursion_ => dim_recursion, &
                                                 dim_subspace_ => dim_subspace
      IMPLICIT NONE
     dimwan = dimwan_
     dim_recursion = dim_recursion_
     dim_subspace = dim_subspace_

   END SUBROUTINE setup_dim


!**********************************************************
   SUBROUTINE setup_subspace()
   !**********************************************************
      USE subspace_variable_module,    ONLY :    x_limit, perform_recursion_x_greater, perform_recursion_x_lesser
      USE input_parameters_module, ONLY :        x_limit_ => x_limit, &
                                                 perform_recursion_x_greater_ => perform_recursion_x_greater, &
                                                 perform_recursion_x_lesser_  => perform_recursion_x_lesser

      IMPLICIT NONE

     x_limit = x_limit_
     perform_recursion_x_greater = perform_recursion_x_greater_
     perform_recursion_x_lesser  = perform_recursion_x_lesser_

   END SUBROUTINE setup_subspace

END MODULE input_module

