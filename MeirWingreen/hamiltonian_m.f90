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
   MODULE hamiltonian_workspace_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE
   PRIVATE 
   SAVE



   ! Public
   INTEGER   ::   dimC
!
! Contains Hamiltonian data
! 
   
   COMPLEX(dbl), ALLOCATABLE :: h00_C(:,:)
   ! Private
   LOGICAL :: alloc = .FALSE.


   ! Public variables
   PUBLIC ::  dimC
   PUBLIC ::  h00_C
 
! PUBLIC ROUTINES
! 

   PUBLIC :: hamiltonian_allocate
   PUBLIC :: hamiltonian_deallocate
   PUBLIC :: hamiltonian_init
 CONTAINS 

!***********************************************
   SUBROUTINE hamiltonian_allocate
    !***********************************************
   IMPLICIT NONE
       CHARACTER(20)      :: subname="hamiltonian_allocate"
       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( dimC <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in dimC definition", " dimC = ", dimC
            STOP
         ENDIF 

       ALLOCATE( h00_C(dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h00_C"
            STOP
       ENDIF 
      h00_C(:,:) = CZERO
       alloc = .TRUE.

   END SUBROUTINE hamiltonian_allocate

!***********************************************
   SUBROUTINE hamiltonian_deallocate
   !***********************************************
      CHARACTER(22)      :: subname="hamiltonian_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED( h00_C ) ) THEN
            DEALLOCATE( h00_C , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h00_C "
               STOP
           ENDIF 
       ENDIF

 
       alloc = .FALSE.
   
   END SUBROUTINE hamiltonian_deallocate
!***********************************************
   SUBROUTINE hamiltonian_init
    !***********************************************
   IMPLICIT NONE
       CHARACTER(16)      :: subname="hamiltonian_init"
       INTEGER  :: ierr, idiag
       CHARACTER(14)         ::     fname_in="EIGENVALUES.in"



        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
        !

        h00_C(:,:)=CZERO
  	PRINT*, "     ... OPENING File:", fname_in
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!

  	PRINT*, "     ...... READING EIGENVALUES"
        DO idiag=1,dimC
           !
           READ(100,"(a11,I10)") h00_C(idiag,idiag)
           PRINT*, "          n=", idiag, " Energy= ",  h00_C(idiag,idiag)
           !
        ENDDO
        !
  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)

       alloc = .TRUE.

   END SUBROUTINE hamiltonian_init

END MODULE hamiltonian_workspace_module
