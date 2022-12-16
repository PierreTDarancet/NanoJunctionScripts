! 
! Copyright (C) 2004 WanT Group
! 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 

!**********************************************************
   SUBROUTINE cleanup()
   !**********************************************************
   !
   ! This module contains the routine CLEANUP that
   ! that deallocates all the data stored in modules
   ! in the WanT-transport code. 
   ! If data is not allocated the routine goes through 
   ! and nothing happens.
   !
   USE timing_module,        ONLY : timing_deallocate, timing_alloc => alloc 
   USE T_egrid_module,       ONLY : egrid_deallocate, egrid_alloc => alloc
   USE in_matrix_module,     ONLY : in_matrix_deallocate, in_mat_alloc => alloc
   USE matrix_module,        ONLY : matrix_deallocate, mat_alloc => alloc
   USE distance_module,      ONLY : distance_deallocate, distance_alloc


   IMPLICIT NONE
      !
      IF ( egrid_alloc )    CALL egrid_deallocate()
      IF ( timing_alloc )   CALL timing_deallocate()
      IF ( in_mat_alloc )   CALL in_matrix_deallocate()
      IF ( mat_alloc )      CALL matrix_deallocate()
      IF ( distance_alloc ) CALL distance_deallocate()

      !

END SUBROUTINE cleanup


