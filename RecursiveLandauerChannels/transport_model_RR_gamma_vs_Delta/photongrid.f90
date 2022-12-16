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
   MODULE photon_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE 
   PRIVATE 
   SAVE
!


   ! Public
   INTEGER                :: nphotonfrequency     ! dimension of the photon grid (input)
   REAL(dbl)              :: dephoton        
   REAL(dbl), ALLOCATABLE :: photonfrequency(:)  ! grid values
   REAL(dbl)              :: photonmax ! max frequency (input)
   REAL(dbl)              :: photonmin ! min frequency (input or 0) 
   REAL(dbl)              :: photoncouplingconstant ! Vo
   REAL(dbl)              :: FGR_broadening ! broadening for Fermi golden rule (lorentzian)

   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: nphotonfrequency
   PUBLIC                 :: dephoton
   PUBLIC                 :: photonfrequency
   PUBLIC                 :: photonmax
   PUBLIC                 :: photonmin
   PUBLIC                 :: photoncouplingconstant
   PUBLIC                 :: FGR_broadening
   ! Public routines:
   PUBLIC                 :: photongrid_allocate
   PUBLIC                 :: photongrid_deallocate
   PUBLIC                 :: photongrid_init

   CONTAINS 
!***********************************************
   SUBROUTINE photongrid_allocate
   !***********************************************
   IMPLICIT NONE
       CHARACTER(19)      :: subname="photongrid_allocate"

       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( nphotonfrequency <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nphotonfrequency definition", " nphotonfrequency = ", nphotonfrequency
            STOP
         ENDIF 

       ALLOCATE( photonfrequency(nphotonfrequency), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating Photonfrequency"
            STOP
       ENDIF 
       photonfrequency(:) = ZERO
       alloc = .TRUE.


  END SUBROUTINE photongrid_allocate
!***********************************************
   SUBROUTINE photongrid_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(21)      :: subname="photongrid_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(photonfrequency) ) THEN
            DEALLOCATE(photonfrequency, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating photonfreqeuncy"
               STOP
           ENDIF 
       ENDIF
       alloc = .FALSE.

   END SUBROUTINE photongrid_deallocate

!***********************************************
   SUBROUTINE photongrid_init
   !***********************************************
  IMPLICIT NONE
    CHARACTER(15)      :: subname="photongrid_init"
       INTEGER :: ierr, ie

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
       !
       IF ( nphotonfrequency > 1 ) THEN
          !
          dephoton = (photonmax - photonmin) / REAL(nphotonfrequency-1, dbl)
          DO ie = 1, nphotonfrequency
              photonfrequency(ie) = photonmin + REAL(ie -1, dbl) * dephoton
          ENDDO
          !
       ELSE
          ! 
          dephoton = ZERO
          photonfrequency(nphotonfrequency)= photonmax
          !
       ENDIF


   END SUBROUTINE photongrid_init


  END  MODULE photon_module

 
