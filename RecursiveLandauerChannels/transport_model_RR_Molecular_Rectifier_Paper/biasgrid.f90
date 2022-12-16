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
   MODULE T_biasgrid_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public

   INTEGER                :: nbias        ! dimension of the energy grid for each bias
   REAL(dbl)              :: dbias        

   REAL(dbl), ALLOCATABLE :: biasgrid(:)  ! grid values

   REAL(dbl)              :: biasmin      !  Minimum energy (input) 
                                       ! (Should be automatized) like H00_C(i,i) min - biasmin - photon_frequency - a bit more
   REAL(dbl)              :: biasmax      !  Maximum energy (input) 
                                       ! (Should be automatized) like H00_C(i,i) max - biasmax + photon_frequency + a bit more

   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: nbias
   PUBLIC                 :: dbias
   PUBLIC                 :: biasgrid

   PUBLIC                 :: biasmin
   PUBLIC                 :: biasmax

   ! Public routines:
   PUBLIC                 :: biasgrid_allocate
   PUBLIC                 :: biasgrid_deallocate
   PUBLIC                 :: biasgrid_init

   CONTAINS 
!***********************************************
   SUBROUTINE biasgrid_allocate
   !***********************************************
   IMPLICIT NONE

       CHARACTER(14)      :: subname="biasgrid_allocate"
       INTEGER  :: ie, ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( nbias <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nbiasmax definition", " nbiasmax = ", nbias
            STOP
         ENDIF 

       ALLOCATE( biasgrid(nbias), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating biasgrid"
            STOP
       ENDIF 
    


  END SUBROUTINE biasgrid_allocate
!***********************************************
   SUBROUTINE biasgrid_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(16)      :: subname="biasgrid_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(biasgrid) ) THEN
            DEALLOCATE(biasgrid, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating biasgrid"
               STOP
           ENDIF 
       ENDIF
       alloc = .FALSE.

   END SUBROUTINE biasgrid_deallocate

!***********************************************
   SUBROUTINE biasgrid_init( )
   !***********************************************
  IMPLICIT NONE


   CHARACTER(10)      :: subname="biasgrid_init"
   INTEGER      :: ie, ierr



          dbias = (biasmax - biasmin) / REAL(nbias-1, dbl)
          DO ie = 1, nbias
              biasgrid(ie) = biasmin + REAL(ie -1, dbl) * dbias
          ENDDO
 
          ! check this part another time
   END SUBROUTINE biasgrid_init

  END  MODULE T_biasgrid_module

 
