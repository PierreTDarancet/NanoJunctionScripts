!
! Copyright (C) 2009 Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!***********************************************
   MODULE  T_control_module
   !***********************************************
   IMPLICIT NONE
   PRIVATE 
   SAVE



   ! Public
!
! Contains control variables
! 
   ! Printed    iterations
   INTEGER   ::  nprint
   ! Printed  PDOS at bias  iterations
   INTEGER   ::  nprint_PDOS
   ! Printed  transmission at bias  iterations
   INTEGER   ::  nprint_transmission
   ! Calculate rectification ratio Y/N
   LOGICAL :: calculate_RR
   ! To be fully implemented 
   ! 0 normal mode
   ! 1-5 debug modes with progressive verbosity
   INTEGER :: debug_level
   ! Printed    Output files
   ! Rectification Ratio
   LOGICAL :: print_RR
   ! Photocurrent
   LOGICAL :: print_Photocurrent
   ! Current
   LOGICAL :: print_current
   ! Phtotocurrent for omega = ELUMO - EHOMO
   LOGICAL :: print_Photocurrentmax

   ! Private variables

   ! Public variables
   PUBLIC :: nprint
   PUBLIC :: nprint_PDOS
   PUBLIC :: nprint_transmission
   PUBLIC :: calculate_RR
   PUBLIC :: debug_level
   PUBLIC :: print_RR
   PUBLIC :: print_Photocurrent
   PUBLIC :: print_current
   PUBLIC :: print_Photocurrentmax

!
! PUBLIC ROUTINES
! 



END MODULE T_control_module
