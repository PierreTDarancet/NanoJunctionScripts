!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!**********************************************************
   SUBROUTINE cleanup()
   !**********************************************************
   !
   ! This module contains the routine CLEANUP that
   ! that deallocates all the data stored in modules
   ! If data is not allocated the routine goes through 
   ! and nothing happens.
   !
   USE hamiltonian_module,     ONLY : hamiltonian_deallocate,  hamiltonian_alloc
   USE orbital_module,         ONLY : orbital_deallocate, orbital_alloc
   USE distance_module,        ONLY : distance_deallocate, distance_alloc
   IMPLICIT NONE
      !
       IF ( orbital_alloc )      CALL orbital_deallocate()
       IF ( distance_alloc )     CALL distance_deallocate()
       IF ( hamiltonian_alloc )  CALL hamiltonian_deallocate()

END SUBROUTINE cleanup


