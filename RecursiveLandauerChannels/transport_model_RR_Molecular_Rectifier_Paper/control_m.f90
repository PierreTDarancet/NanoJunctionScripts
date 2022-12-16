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
  USE kinds, ONLY : dbl
   IMPLICIT NONE
   PRIVATE 
   SAVE



   ! Public
!
! Contains control variables
! 

   ! To be fully implemented 
   ! 0 normal mode
   ! 1-5 debug modes with progressive verbosity
   INTEGER :: debug_level
   REAL(dbl) :: gamma_big
   REAL(dbl) :: gamma_small
   ! Printed    Output files
   ! Rectification Ratio
    !PUBLIC :: nprint
   !PUBLIC :: nprint_PDOS
   !PUBLIC :: nprint_transmission

   PUBLIC :: debug_level

   PUBLIC :: gamma_big
   PUBLIC :: gamma_small

   !PUBLIC :: print_Photocurrentmax

!
! PUBLIC ROUTINES
! 



END MODULE T_control_module
