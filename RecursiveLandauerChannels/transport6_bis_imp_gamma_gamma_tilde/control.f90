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
   INTEGER                   :: in_max_iter_C

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
   !
   PUBLIC :: nprint
   !
   !
   PUBLIC :: dim_subspace
   !
   PUBLIC :: in_max_iter_C
   !
   PUBLIC :: print_gamma
   !
   PUBLIC :: imp_gamma
   !
   !

END MODULE T_control_module

