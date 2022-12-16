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
   MODULE kpoint_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI

   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   REAL(dbl), ALLOCATABLE ::  kweight(:) ! 
   INTEGER :: nk
 
   ! Private
   LOGICAL :: alloc = .FALSE.


   ! Public variables
   PUBLIC                 :: nk
   PUBLIC                 :: kweight

   ! Public routines:

   PUBLIC                 :: kgrid_allocate
   PUBLIC                 :: kgrid_deallocate
   PUBLIC                 :: kgrid_init

   CONTAINS 
!***********************************************
   SUBROUTINE kgrid_allocate
   !***********************************************
   IMPLICIT NONE
       CHARACTER(14)      :: subname="kgrid_allocate"

       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( nk <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nk definition", " nk ", nk
            STOP
         ENDIF 


       ALLOCATE( kweight(nk), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating kgrid"
            STOP
       ENDIF 



       kweight(:) = ZERO

       alloc = .TRUE.


  END SUBROUTINE kgrid_allocate
!***********************************************
   SUBROUTINE kgrid_deallocate
   !***********************************************
  IMPLICIT NONE
       CHARACTER(16)      :: subname="kgrid_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(kweight) ) THEN
            DEALLOCATE(kweight, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating kweight"
               STOP
           ENDIF 
       ENDIF


       alloc = .FALSE.

   END SUBROUTINE kgrid_deallocate

!***********************************************
   SUBROUTINE kgrid_init()
   !***********************************************
  IMPLICIT NONE
    CHARACTER(10)         ::     fname_in="KPOINTS.in"
    CHARACTER(10)      :: subname="kgrid_init"
    INTEGER :: ierr, ik

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
        !

  	PRINT*, "     ... OPENING File:", fname_in
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!
  	PRINT*, "     ...... READING kweight"
        DO ik=1,nk
           !
           READ(100,"(F10.10)")  kweight(ik)
           PRINT*, "          ik=", ik, " weight= ", kweight(ik)
           !
        ENDDO
        !
  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)
   END SUBROUTINE kgrid_init


  END  MODULE kpoint_module

 
