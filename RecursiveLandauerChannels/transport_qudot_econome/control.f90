!
! Copyright (C) 2005 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*********************************************
   MODULE T_control_module
!*********************************************
   USE kinds,      ONLY : dbl
   USE parameters, ONLY : nstrx
   IMPLICIT NONE
   PRIVATE 
   SAVE
!
! Contains GLOBAL CONTROL variables for transport calculations
! 
   !
   CHARACTER(nstrx)          :: in_datafile_C
   !
   CHARACTER(nstrx)          :: in_datafile_gamma_R
   !
   CHARACTER(nstrx)          :: in_datafile_gamma_L
   !
   CHARACTER(nstrx)          :: in_datafile_sigma_R
   !
   CHARACTER(nstrx)          :: in_datafile_sigma_L
   !
   CHARACTER(nstrx)          :: datafile_C_form
   !
!   INTEGER                   :: transport_dir
   !
!   INTEGER                   :: niterx
   !
   INTEGER                   :: nprint
   !
   INTEGER                   :: dim_subspace
   !
   INTEGER                   :: in_max_iter_C
   !
   INTEGER                   :: dimC
   !
   LOGICAL                   :: print_gamma
   !
   LOGICAL                   :: imp_gamma
   !
!
! end delcarations
!
   !
   PUBLIC :: in_datafile_C
   !
   PUBLIC :: in_datafile_sigma_R
   !
   PUBLIC :: in_datafile_gamma_R
   !
   PUBLIC :: in_datafile_sigma_L
   !
   PUBLIC :: in_datafile_gamma_L
   !
   PUBLIC :: datafile_C_form
   !
   PUBLIC :: nprint
   !
   PUBLIC :: dim_subspace
   !
   PUBLIC :: in_max_iter_C
   !
   PUBLIC :: dimC
   !
   PUBLIC :: print_gamma
   !
   PUBLIC :: imp_gamma
   !

END MODULE T_control_module

