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
   MODULE T_gapgrid_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE
   PRIVATE 
   SAVE
!



   ! Public
   INTEGER                :: ngap     ! dimension of the gap grid (input)
   REAL(dbl)              :: dgap        
   REAL(dbl), ALLOCATABLE :: gapgrid(:)  ! grid values
   REAL(dbl)              :: gapmax ! max gap (input) 
   REAL(dbl)              :: gapmin ! min gap (input). At this stage of the code and for rectification ratio purpose, the gap is intended to be symmetric
                                     ! A non symmetric gap grid while calculating RR should have only bad consequences on this part of the code
                                     ! This point could be easily improved, without stopping to calculate the RR, by calculating the RR in this module
                                     ! rather than in the main source

   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: ngap
   PUBLIC                 :: dgap
   PUBLIC                 :: gapgrid
   PUBLIC                 :: gapmax
   PUBLIC                 :: gapmin


   ! Public routines:
   PUBLIC                 :: gapgrid_allocate
   PUBLIC                 :: gapgrid_deallocate
   PUBLIC                 :: gapgrid_init

   CONTAINS 
!***********************************************
   SUBROUTINE gapgrid_allocate
   !***********************************************
   IMPLICIT NONE
       CHARACTER(17)      :: subname="gapgrid_allocate"
       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( ngap <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in ngap definition", " ngap = ", ngap
            STOP
         ENDIF 

       ALLOCATE( gapgrid(ngap), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating gapgrid"
            STOP
       ENDIF 
       gapgrid(:) = ZERO
       alloc = .TRUE.


  END SUBROUTINE gapgrid_allocate
!***********************************************
   SUBROUTINE gapgrid_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(19)      :: subname="gapgrid_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(gapgrid) ) THEN
            DEALLOCATE(gapgrid, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating gapgrid"
               STOP
           ENDIF 
       ENDIF
       alloc = .FALSE.

   END SUBROUTINE gapgrid_deallocate

!***********************************************
   SUBROUTINE gapgrid_init
   !***********************************************
  IMPLICIT NONE
    CHARACTER(13)      :: subname="gapgrid_init"
       INTEGER :: ierr, igap

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
       !
       IF ( ngap > 1 ) THEN
          !
          dgap = (gapmax - gapmin) / REAL(ngap -1, dbl)
          !
          DO igap = 1, ngap
             gapgrid(igap) = gapmin + REAL(igap -1, dbl) * dgap
          ENDDO
          !
       ELSE ! if ngap = 1 
          ! 
          dgap = ZERO
          gapgrid(ngap)= gapmax
          !
       ENDIF


   END SUBROUTINE gapgrid_init


  END  MODULE T_gapgrid_module


