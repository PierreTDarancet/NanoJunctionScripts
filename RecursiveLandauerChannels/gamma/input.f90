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
                                           in_datafile,       &
                                           in_max_iter,       &
                                           max_iter_renorm,   &
                                           max_iter_term,     &
                                           print_gamma,       &
                                           print_gamma_tilde, &
                                           method_sigma,      &
                                           dim_subspace,      &
                                           a_analytique,      &
                                           b_analytique

      USE T_input_parameters_module,ONLY : nprint_            => nprint, &
                                           in_datafile_       => in_datafile, &
                                           in_max_iter_       => in_max_iter, &
                                           max_iter_renorm_   => max_iter_renorm,    &
                                           max_iter_term_     => max_iter_term,    &
                                           dim_subspace_      => dim_subspace,  &
                                           print_gamma_       => print_gamma,   &
                                           print_gamma_tilde_ => print_gamma_tilde,  &
                                           method_sigma_      => method_sigma,    &
                                           a_analytique_      => a_analytique,    &
                                           b_analytique_      => b_analytique


      IMPLICIT NONE

      nprint              = nprint_
      in_datafile         = in_datafile_
      dim_subspace        = dim_subspace_
      in_max_iter         = in_max_iter_
      max_iter_renorm     = max_iter_renorm_
      max_iter_term       = max_iter_term_
      print_gamma         = print_gamma_
      print_gamma_tilde   = print_gamma_tilde_
      method_sigma        = method_sigma_
      a_analytique        = a_analytique_
      b_analytique        = b_analytique_

   END SUBROUTINE setup_control


!**********************************************************
   SUBROUTINE setup_egrid()
   !**********************************************************
      USE T_egrid_module,           ONLY : ne,           &
                                           emin, emax,   &
                                           delta_lead
      USE T_input_parameters_module,ONLY : ne_     => ne,   &
                                           emin_   => emin, &
                                           emax_   => emax, &
                                           delta_lead_  => delta_lead
      IMPLICIT NONE

      ne     = ne_
      emin   = emin_
      emax   = emax_
      delta_lead  = delta_lead_

   END SUBROUTINE setup_egrid
    !
END MODULE T_input_module

