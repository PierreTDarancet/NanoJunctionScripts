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
   MODULE convert3_input_module
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
! SUBROUTINE convert_input_manager()
! SUBROUTINE convert_setup_control()
! SUBROUTINE convert_setup_io()

!

!
   PUBLIC :: convert3_input_manager

!

CONTAINS

!
! subroutines
!


!**********************************************************
   SUBROUTINE convert3_input_manager()
   !**********************************************************
      USE convert3_input_parameters_module,  ONLY : convert3_read_namelist_control

     IMPLICIT NONE
      INTEGER :: ierr
      CHARACTER(22) :: subname='convert3_input_manager'

      !
      ! reading and checking namelists
      !

      !!!!!!!!!!!!!!!!!!!!!!!!have to be improved : for the moment all 
      !!!!!!!!!!!!!!!!!!!!!!!!the variables are supposed to be read in the control 
      !!!!!!!!!!!!!!!!!!!!!!!!namelist even if they are next stored in different variables module

      CALL convert3_read_namelist_control(stdin)



      !
      ! scattering data in their own modules
      !
      CALL convert3_setup_control()
      CALL convert3_setup_io()
      !CALL convert_setup_dim()
      !CALL convert_setup_subspace()


   END SUBROUTINE convert3_input_manager

!**********************************************************
   SUBROUTINE convert3_setup_control()
   !**********************************************************
      USE convert3_control_variable_module,     ONLY :    out_datafile_C, out_datafile_L, out_datafile_R, &
                                                         in_datafile_L, in_datafile_R, &
                                                         dim_subspace, &
                                                         n_iter_R, n_iter_L, &
                                                         out_iter_R, out_iter_L, &
                                                         out_iter_C
      USE convert3_input_parameters_module, ONLY :   out_datafile_C_ => out_datafile_C, &
                                                     out_datafile_L_ => out_datafile_L, &
                                                     out_datafile_R_ => out_datafile_R, &
                                                     in_datafile_L_ => in_datafile_L, &
                                                     in_datafile_R_ => in_datafile_R, &
                                                     dim_subspace_ => dim_subspace, &
                                                     n_iter_R_ => n_iter_R,     &
                                                     n_iter_L_ => n_iter_L,     &
                                                     out_iter_R_ => out_iter_R, &
                                                     out_iter_L_ => out_iter_L, &
                                                     out_iter_C_ => out_iter_C

      IMPLICIT NONE

     out_datafile_C = out_datafile_C_
     out_datafile_L = out_datafile_L_
     out_datafile_R = out_datafile_R_
     out_iter_R = out_iter_R_
     out_iter_L = out_iter_L_
     out_iter_C = out_iter_C_

     in_datafile_L = in_datafile_L_
     in_datafile_R = in_datafile_R_
     dim_subspace = dim_subspace_
     n_iter_R = n_iter_R_
     n_iter_L = n_iter_L_

   END SUBROUTINE convert3_setup_control
! 
! !**********************************************************
!    SUBROUTINE convert_setup_dim()
!    !**********************************************************
!       USE convert_dim_variable_module,    ONLY :  nb_replica_L, nb_replica_R, dimwan, dim_recursion, dim_subspace, dim_L, dim_R, dimwan_total
!       USE convert_input_parameters_module, ONLY :  dimwan_ => dimwan, &
!                                            dim_recursion_ => dim_recursion, &
!                                            dim_subspace_ => dim_subspace, &
!                                            nb_replica_L_ => nb_replica_L, &
!                                            nb_replica_R_ => nb_replica_R, &
!                                            dim_L_ => dim_L, &
!                                            dim_R_ => dim_R, &
!                                            use_extension_R_ => use_extension_R, &
!                                            use_extension_L_ => use_extension_L
! 
!       IMPLICIT NONE
!       CHARACTER(9) :: subname='setup_dim'
! 
! 
!      dimwan = dimwan_
!      dim_recursion = dim_recursion_
!      dim_subspace = dim_subspace_
!      dim_L = dim_L_
!      dim_R = dim_R_
!      nb_replica_L = nb_replica_L_
!      nb_replica_R = nb_replica_R_
!      IF (use_extension_R_) THEN
!         dimwan_total= dimwan + (nb_replica_R + 1) * dim_R
!         IF ( dim_recursion >  dimwan_total ) &
!                CALL errore(subname,'dim_recursion > dimwan_total',1)
! 
!      ELSE IF (use_extension_L_) THEN
!         dimwan_total= dimwan + (nb_replica_L + 1) * dim_L
!         IF ( dim_recursion >  dimwan_total ) &
!                CALL errore(subname,'dim_recursion > dimwan_total',2)
! 
!      ELSE
!         dimwan_total = dimwan
!         IF ( dim_recursion >  dimwan_total ) &
!                CALL errore(subname,'dim_recursion > dimwan_total',3)
! 
!      ENDIF
! 
!    END SUBROUTINE convert_setup_dim

! 
! !**********************************************************
!    SUBROUTINE control_setup_subspace()
!    !**********************************************************
!       USE convert_subspace_variable_module,    ONLY :    x_limit, perform_recursion_x_greater, perform_recursion_x_lesser, &
!                                                  cell_size_C, cell_size_L, cell_size_R
!       USE convert_input_parameters_module, ONLY :        x_limit_ => x_limit, &
!                                                  perform_recursion_x_greater_ => perform_recursion_x_greater, &
!                                                  perform_recursion_x_lesser_  => perform_recursion_x_lesser, &
!                                                  cell_size_C_ => cell_size_C, &
!                                                  cell_size_L_ => cell_size_L, &
!                                                  cell_size_R_ => cell_size_R
!       IMPLICIT NONE
! 
!      x_limit = x_limit_
!      cell_size_C = cell_size_C_
!      cell_size_L = cell_size_L_
!      cell_size_R = cell_size_R_
! 
!    END SUBROUTINE convert_setup_subspace

!**********************************************************
   SUBROUTINE convert3_setup_io()
   !**********************************************************
      USE io_module,                ONLY : prefix, postfix, work_dir, title
      USE convert_input_parameters_module,  ONLY : prefix_    => prefix,       &
                                           postfix_   => postfix,      &
                                           work_dir_  => work_dir,     &
                                           title_     => title
      IMPLICIT NONE
      prefix   = prefix_
      postfix  = postfix_
      work_dir = work_dir_
      title    = title_
   END SUBROUTINE convert3_setup_io


END MODULE convert3_input_module

