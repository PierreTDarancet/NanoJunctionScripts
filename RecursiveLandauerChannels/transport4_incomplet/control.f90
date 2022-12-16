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
   
   CHARACTER(nstrx)          :: calculation_type
   CHARACTER(nstrx)          :: conduct_formula
   !
   CHARACTER(nstrx)          :: in_datafile_L, in_datafile_R
   !
   INTEGER                   :: transport_dir
   !
   !
!   INTEGER                   :: niterx
   !
   INTEGER                   :: nprint
   !
   !
   INTEGER                   :: dim_subspace
   !
   INTEGER                   :: max_iter_R
   !
   INTEGER                   :: max_iter_L
   !
   INTEGER                   :: max_iter_CR
   !
   INTEGER                   :: max_iter_LC
   !
   INTEGER                   :: in_max_iter_R
   !
   INTEGER                   :: in_max_iter_L
   !
   LOGICAL                   :: debug_mode
   !
   LOGICAL                   :: cut_chain
   !
   REAL(dbl)                 :: conv_criterion
   !
   REAL(dbl)                 :: numerical_cut_off

!

! end delcarations
!

   PUBLIC :: calculation_type
   !
   PUBLIC :: conduct_formula
   !
   PUBLIC :: in_datafile_L,  in_datafile_R
   !
   PUBLIC :: transport_dir
   !
   !
!   PUBLIC :: niterx
   !
   PUBLIC :: nprint
   !
   !
   PUBLIC :: dim_subspace
   !
   PUBLIC :: max_iter_R
   !
   PUBLIC :: max_iter_L
   !
   PUBLIC :: max_iter_CR
   !
   PUBLIC :: max_iter_LC
   !
   PUBLIC :: in_max_iter_R
   !
   PUBLIC :: in_max_iter_L
   !
   PUBLIC :: debug_mode
   !
   PUBLIC :: cut_chain
   !
   PUBLIC :: conv_criterion
   !
   PUBLIC :: numerical_cut_off
   !

END MODULE T_control_module

