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
   MODULE T_gammagrid_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   INTEGER                :: ngamma     ! dimension of the gamma grid (input)
   REAL(dbl)              :: d_gamma        
   REAL(dbl), ALLOCATABLE :: gammagrid(:)  ! grid values
   REAL(dbl)              :: gammamax ! max gamma (input) 
   REAL(dbl)              :: gammamin ! min gamma (input). 
   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: ngamma
   PUBLIC                 :: d_gamma
   PUBLIC                 :: gammagrid
   PUBLIC                 :: gammamax
   PUBLIC                 :: gammamin


   ! Public routines:
   PUBLIC                 :: gammagrid_allocate
   PUBLIC                 :: gammagrid_deallocate
   PUBLIC                 :: gammagrid_init

   CONTAINS 
!***********************************************
   SUBROUTINE gammagrid_allocate
   !***********************************************
   IMPLICIT NONE
       CHARACTER(18)      :: subname="gammagrid_allocate"
       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( ngamma <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in ngamma definition", " ngamma = ", ngamma
            STOP
         ENDIF 

       ALLOCATE( gammagrid(ngamma), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating gammagrid"
            STOP
       ENDIF 
      gammagrid(:) = ZERO
       alloc = .TRUE.


  END SUBROUTINE gammagrid_allocate
!***********************************************
   SUBROUTINE gammagrid_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(20)      :: subname="gammagrid_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(gammagrid) ) THEN
            DEALLOCATE(gammagrid, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating gammagrid"
               STOP
           ENDIF 
       ENDIF
       alloc = .FALSE.

   END SUBROUTINE gammagrid_deallocate

!***********************************************
   SUBROUTINE gammagrid_init
   !***********************************************
  IMPLICIT NONE
    CHARACTER(14)      :: subname="gammagrid_init"
       INTEGER :: ierr, igamma

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
       !
       IF (ngamma >1) THEN
          !
          d_gamma = (log(gammamax) - log(gammamin)) / REAL(ngamma -1, dbl)
          !
          DO igamma = 1, ngamma
             gammagrid(igamma) = exp( log(gammamin) + REAL(igamma -1, dbl) * d_gamma )
          ENDDO
          !
       ELSE ! if ngamma = 1 
          ! 
          d_gamma = ZERO
          gammagrid(ngamma)= gammamax
          !
       ENDIF


   END SUBROUTINE gammagrid_init


  END  MODULE T_gammagrid_module


