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
   CHARACTER(nstrx)          :: in_datafile
   !
   CHARACTER(nstrx)          :: method_sigma
   !

!   INTEGER                   :: transport_dir
   !
   !
!   INTEGER                   :: niterx
   !
   INTEGER                   :: nprint
   !
   !
   INTEGER                   :: dim_subspace
   !
   INTEGER                   :: max_iter_term
   !
   INTEGER                   :: max_iter_renorm
   !
   INTEGER                   :: in_max_iter
   !
   LOGICAL                   :: debug_mode=.FALSE.
   !
   LOGICAL                   :: print_gamma
   !
   LOGICAL                   :: print_gamma_tilde
   !
   REAL(dbl)                 :: a_analytique
   !
   REAL(dbl)                 :: b_analytique
!

! end delcarations
!

   !
   PUBLIC :: in_datafile
   !
!   PUBLIC :: transport_dir
   !
   !
!   PUBLIC :: niterx
   !
   PUBLIC :: nprint
   !
   PUBLIC :: dim_subspace
   !
   PUBLIC :: max_iter_term
   !
   PUBLIC :: max_iter_renorm
   !
   PUBLIC :: in_max_iter
   !
   PUBLIC :: print_gamma
   !
   PUBLIC :: print_gamma_tilde
   !
   PUBLIC :: method_sigma
   !
   PUBLIC :: a_analytique
   !
   PUBLIC :: b_analytique
   !

END MODULE T_control_module

