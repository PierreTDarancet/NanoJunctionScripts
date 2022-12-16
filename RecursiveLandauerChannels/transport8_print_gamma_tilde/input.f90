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
                                           in_datafile_L,     &
                                           in_datafile_R,     &
                                           in_max_iter_L,     &
                                           in_max_iter_R,     &
                                           max_iter_L,        &
                                           max_iter_R,        &
                                           dim_subspace,      &
                                           conv_criterion,    &
                                           numerical_cut_off, &
                                           max_iter_CL,       &
                                           max_iter_CR,       &
                                           max_iter_CC,       &
                                           cut_chain,         &
                                           print_gamma

      USE T_input_parameters_module,ONLY : nprint_            => nprint, &
                                           in_datafile_L_     => in_datafile_L, &
                                           in_datafile_R_     => in_datafile_R, &
                                           in_max_iter_L_     => in_max_iter_L, &
                                           in_max_iter_R_     => in_max_iter_R, &
                                           max_iter_L_        => max_iter_L,    &
                                           max_iter_R_        => max_iter_R,    &
                                           dim_subspace_      => dim_subspace,  &
                                           conv_criterion_    => conv_criterion, &
                                           numerical_cut_off_ => numerical_cut_off, &
                                           max_iter_CL_       => max_iter_CL,   &
                                           max_iter_CR_       => max_iter_CR,   &
                                           max_iter_CC_       => max_iter_CC,   &
                                           cut_chain_         => cut_chain,     &
                                           print_gamma_       => print_gamma
      IMPLICIT NONE

      nprint              = nprint_
      in_datafile_L       = in_datafile_L_
      in_datafile_R       = in_datafile_R_
      dim_subspace        = dim_subspace_
      in_max_iter_R       = in_max_iter_R_
      in_max_iter_L       = in_max_iter_L_
      max_iter_R          = max_iter_R_
      max_iter_L          = max_iter_L_
      conv_criterion      = conv_criterion_
      numerical_cut_off   = numerical_cut_off_
      max_iter_CR         = max_iter_CR_
      max_iter_CL         = max_iter_CL_
      max_iter_CC         = max_iter_CC_
      cut_chain           = cut_chain_
      print_gamma         = print_gamma_

   END SUBROUTINE setup_control


!**********************************************************
   SUBROUTINE setup_egrid()
   !**********************************************************
      USE T_egrid_module,           ONLY : ne,           &
                                           emin, emax,   &
                                           delta, delta_lead
      USE T_input_parameters_module,ONLY : ne_     => ne,   &
                                           emin_   => emin, &
                                           emax_   => emax, &
                                           delta_  => delta, &
                                           delta_lead_  => delta_lead
      IMPLICIT NONE

      ne     = ne_
      emin   = emin_
      emax   = emax_
      delta  = delta_
      delta_lead  = delta_lead_

   END SUBROUTINE setup_egrid


!    !**********************************************************
!       SUBROUTINE setup_hamiltonian()
!       !**********************************************************
!          USE T_hamiltonian_module, ONLY :     dimC, &
!                                              dimR, &
!                                              dimL
!    
!          USE T_input_parameters_module,ONLY : dimC_              => dimC, &
!                                              dimR_              => dimR, &
!                                              dimL_              => dimL
!    
!          IMPLICIT NONE
!    
!          dimL                = dimL_
!          dimR                = dimR_
!          dimC                = dimC_
!    
!    
!       END SUBROUTINE setup_hamiltonian
!    
END MODULE T_input_module

