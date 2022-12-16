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
   ! Printed    iterations
   !INTEGER   ::  nprint
   ! Printed  PDOS at bias  iterations
   !INTEGER   ::  nprint_PDOS
   ! Printed  transmission at bias  iterations
   !INTEGER   ::  nprint_transmission
   ! Calculate rectification ratio Y/N
   !LOGICAL :: calculate_RR
   ! To be fully implemented 
   ! 0 normal mode
   ! 1-5 debug modes with progressive verbosity
   INTEGER :: debug_level
   INTEGER :: nmaxiter_gf
   CHARACTER(3) :: ee_method !DFT or SIG or EXC
   CHARACTER(2) :: ep_method !AM Aeberhard Morf PRB 77 (2008) - SCBA
                             !GN Galperin Nitzan PRL 2005 with real dipole transition
                             !G0 Galperin Nitzan PRL 2005 without real dipole transition
                             ! NO No electron-photon contribution
  
   ! Printed    Output files
   ! Rectification Ratio
   !LOGICAL :: print_RR
   ! Photocurrent
   !LOGICAL :: print_Photocurrent
   ! Current
   !LOGICAL :: print_current
   ! Phtotocurrent for omega = ELUMO - EHOMO
   !LOGICAL :: print_Photocurrentmax
   REAL(dbl) :: amplitude
   REAL(dbl) :: photonfrequency
   ! Private variables

   ! Public variables
   PUBLIC :: debug_level
   PUBLIC :: nmaxiter_gf
   PUBLIC :: ee_method 
   PUBLIC :: ep_method 
   PUBLIC ::  amplitude
   PUBLIC :: photonfrequency
 
!
! PUBLIC ROUTINES
! 



END MODULE T_control_module
