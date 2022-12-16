! 
! Copyright (C) 2005 WanT Group
! 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
!********************************************
   MODULE T_input_module
!********************************************
   USE kinds, ONLY : dbl
   USE io_module, ONLY : stdin, stdout
   USE constants, ONLY : ZERO
   IMPLICIT NONE
   PRIVATE
!
! This module handles the reading of input data
!
! routines in this module:
! SUBROUTINE input_manager()
! SUBROUTINE setup_control()
! SUBROUTINE setup_hamiltonian()
! SUBROUTINE setup_egrid()
! 


   PUBLIC :: input_manager


CONTAINS

!
! subroutines
!


!**********************************************************
   SUBROUTINE input_manager()
   !**********************************************************
      USE T_input_parameters_module,  ONLY : read_namelist_input_conductor
      IMPLICIT NONE

      !
      ! reading and checking namelists
      !
      CALL read_namelist_input_conductor(stdin)

      !
      ! scattering data in their own modules
      !
      CALL setup_control()
      CALL setup_egrid()
!       CALL setup_hamiltonian()

      !
      ! reading further input data
      !

   END SUBROUTINE input_manager


!**********************************************************
   SUBROUTINE setup_control()
   !**********************************************************

      USE T_control_module,         ONLY : nprint,            &
                                           in_datafile_C,     &
                                           in_max_iter_C,     &
                                           dim_subspace,      &
                                           print_gamma,       &
                                           imp_gamma,         &
                                           dimC,              &
                                           in_datafile_gamma_R, &
                                           in_datafile_gamma_L, &
                                           in_datafile_sigma_R, &
                                           in_datafile_sigma_L, &
                                           datafile_C_form



      USE T_input_parameters_module,ONLY : nprint_            => nprint,        &
                                           in_datafile_C_     => in_datafile_C, &
                                           in_max_iter_C_     => in_max_iter_C, &
                                           dim_subspace_      => dim_subspace,  &
                                           print_gamma_       => print_gamma,   &
                                           imp_gamma_         => imp_gamma,     &
                                           dimC_              => dimC,          &
                                           in_datafile_gamma_R_  => in_datafile_gamma_R, &
                                           in_datafile_gamma_L_  => in_datafile_gamma_L, &
                                           in_datafile_sigma_R_  => in_datafile_sigma_R, &
                                           in_datafile_sigma_L_  => in_datafile_sigma_L, &
                                           datafile_C_form_   => datafile_C_form

      IMPLICIT NONE

      nprint              = nprint_
      in_datafile_C       = in_datafile_C_
      dim_subspace        = dim_subspace_
      in_max_iter_C       = in_max_iter_C_
      print_gamma         = print_gamma_
      imp_gamma           = imp_gamma_
      dimC                = dimC_
      in_datafile_gamma_R = in_datafile_gamma_R_
      in_datafile_gamma_L = in_datafile_gamma_L_
      in_datafile_sigma_R = in_datafile_sigma_R_
      in_datafile_sigma_L = in_datafile_sigma_L_
      datafile_C_form     = datafile_C_form_


   END SUBROUTINE setup_control


!**********************************************************
   SUBROUTINE setup_egrid()
   !**********************************************************

      USE T_egrid_module,           ONLY : ne,           &
                                           emin, emax,   &
                                           delta
      USE T_input_parameters_module,ONLY : ne_     => ne,   &
                                           emin_   => emin, &
                                           emax_   => emax, &
                                           delta_  => delta
      IMPLICIT NONE

      ne     = ne_
      emin   = emin_
      emax   = emax_
      delta  = delta_
      !

   END SUBROUTINE setup_egrid


END MODULE T_input_module

