!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE data_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE
   USE identity_module,    ONLY : orbital_def

  IMPLICIT NONE
   PRIVATE
   SAVE
   !
   !


   INTEGER, PARAMETER :: n_id=2
   !
   TYPE (orbital_def) :: orb_id(n_id)


   !
   PUBLIC :: orb_id
   PUBLIC :: n_id
   !

   PUBLIC :: orb_init


    CONTAINS
!*********************************************
   SUBROUTINE orb_init
!*********************************************
  CHARACTER(8) :: subname="orb_init"
!

   orb_id(1)%name   = "pz"
   orb_id(1)%onsite = CZERO
   orb_id(1)%s(:)   = CZERO
   orb_id(1)%px(:)  = CZERO
   orb_id(1)%py(:)  = CZERO
   orb_id(1)%pz(1)  = 3.000*CONE
   orb_id(1)%pz(2)  = 0.160*CONE
   orb_id(1)%pz(3)  = CZERO
   orb_id(1)%d(:)   = CZERO
   orb_id(1)%f(:)   = CZERO

   orb_id(2)%name   = "px"
   orb_id(2)%onsite = -4.000*CONE
   orb_id(2)%s(:)   = CZERO
   orb_id(2)%px(:)  = CZERO
   orb_id(2)%py(:)  = CZERO
   orb_id(2)%pz(1)  = 6.000*CONE
   orb_id(2)%pz(2)  = CZERO
   orb_id(2)%pz(3)  = CZERO
   orb_id(2)%d(:)   = CZERO
   orb_id(2)%f(:)   = CZERO

   END SUBROUTINE orb_init


END MODULE data_module