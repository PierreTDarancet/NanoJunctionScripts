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
   MODULE T_alphagrid_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   INTEGER                :: nalpha     ! dimension of the alpha grid (input)
   REAL(dbl)              :: dalpha        
   REAL(dbl), ALLOCATABLE :: alphagrid(:)  ! grid values
   REAL(dbl)              :: alphamax ! max alpha (input) 
   REAL(dbl)              :: alphamin ! min alpha (input). 
   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: nalpha
   PUBLIC                 :: dalpha
   PUBLIC                 :: alphagrid
   PUBLIC                 :: alphamax
   PUBLIC                 :: alphamin


   ! Public routines:
   PUBLIC                 :: alphagrid_allocate
   PUBLIC                 :: alphagrid_deallocate
   PUBLIC                 :: alphagrid_init

   CONTAINS 
!***********************************************
   SUBROUTINE alphagrid_allocate
   !***********************************************
   IMPLICIT NONE
       CHARACTER(18)      :: subname="alphagrid_allocate"
       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( nalpha <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nalpha definition", " nalpha = ", nalpha
            STOP
         ENDIF 

       ALLOCATE( alphagrid(nalpha), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating alphagrid"
            STOP
       ENDIF 
      alphagrid(:) = ZERO
       alloc = .TRUE.


  END SUBROUTINE alphagrid_allocate
!***********************************************
   SUBROUTINE alphagrid_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(20)      :: subname="alphagrid_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(alphagrid) ) THEN
            DEALLOCATE(alphagrid, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating alphagrid"
               STOP
           ENDIF 
       ENDIF
       alloc = .FALSE.

   END SUBROUTINE alphagrid_deallocate

!***********************************************
   SUBROUTINE alphagrid_init
   !***********************************************
  IMPLICIT NONE
    CHARACTER(14)      :: subname="alphagrid_init"
       INTEGER :: ierr, ialpha

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
       !
       IF ( nalpha > 1 ) THEN
          !
          dalpha = (alphamax - alphamin) / REAL(nalpha -1, dbl)
          !
          DO ialpha = 1, nalpha
             alphagrid(ialpha) = alphamin + REAL(ialpha -1, dbl) * dalpha
          ENDDO
          !
       ELSE ! if nalpha = 1 
          ! 
          dalpha = ZERO
          alphagrid(nalpha)= alphamax
          !
       ENDIF


   END SUBROUTINE alphagrid_init


  END  MODULE T_alphagrid_module


