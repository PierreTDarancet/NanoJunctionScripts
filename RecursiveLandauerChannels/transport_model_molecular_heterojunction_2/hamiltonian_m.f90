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
  INTEGER   ::   dimL
  INTEGER   ::   dimR
  INTEGER   ::   dimC
!
! Contains Hamiltonian data
! 
   
   COMPLEX(dbl), ALLOCATABLE :: h00_C(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h_CR(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h_RC(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h_LC(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h_CL(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h00_R(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h01_R(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h00_L(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h01_L(:,:)
   ! Private
   LOGICAL :: alloc = .FALSE.


   ! Public variables
   PUBLIC ::  dimC
   PUBLIC ::  dimR
   PUBLIC ::  dimL
   PUBLIC ::  h00_C
   PUBLIC ::  h_CR
   PUBLIC ::  h_LC

   PUBLIC ::  h_RC
   PUBLIC ::  h_CL
   PUBLIC ::  h00_R
   PUBLIC ::  h01_R
   PUBLIC ::  h00_L
   PUBLIC ::  h01_L

!
! PUBLIC ROUTINES
! 

   PUBLIC :: hamiltonian_allocate
   PUBLIC :: hamiltonian_deallocate

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

        IF( dimL <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in dimL definition", " dimL = ", dimL
            STOP
         ENDIF 

        IF( dimR <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in dimR definition", " dimR = ", dimR
            STOP
         ENDIF 

       ALLOCATE( h00_C(dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h00_C"
            STOP
       ENDIF 

       ALLOCATE( h00_R(dimR,dimR)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h00_R"
            STOP
       ENDIF 
       ALLOCATE( h00_L(dimL,dimL)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h00_L"
            STOP
       ENDIF 
       ALLOCATE( h01_L(dimL,dimL)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h01_L"
            STOP
       ENDIF 
       ALLOCATE( h01_R(dimR,dimR)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h01_R"
            STOP
       ENDIF 
       ALLOCATE( h_LC(dimL,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h_LC"
            STOP
       ENDIF 
       ALLOCATE( h_CR(dimC,dimR)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h_CR"
            STOP
       ENDIF 
       ALLOCATE( h_CL(dimC,dimL)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h_CL"
            STOP
       ENDIF 
       ALLOCATE( h_RC(dimR,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h_RC"
            STOP
       ENDIF 
       h_RC(:,:)  = CZERO
       h_CL(:,:)  = CZERO
       h_LC(:,:)  = CZERO 
       h00_L(:,:) = CZERO
       h00_R(:,:) = CZERO 
       h01_L(:,:) = CZERO 
       h01_R(:,:) = CZERO 
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

       IF ( ALLOCATED( h_CR ) ) THEN
            DEALLOCATE( h_CR , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h_CR"
               STOP
           ENDIF 
       ENDIF

      IF ( ALLOCATED( h_RC ) ) THEN
            DEALLOCATE( h_RC , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h_RC"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( h_LC   ) ) THEN
            DEALLOCATE( h_LC   , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  h_LC "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( h_CL   ) ) THEN
            DEALLOCATE( h_CL   , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  h_CL "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(h00_L ) ) THEN
            DEALLOCATE( h00_L, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h00_L"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(h00_R  ) ) THEN
            DEALLOCATE( h00_R , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h00_R "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(h01_L ) ) THEN
            DEALLOCATE( h01_L, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h01_L"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( h01_R  ) ) THEN
            DEALLOCATE( h01_R , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h01_R "
               STOP
           ENDIF 
       ENDIF

       alloc = .FALSE.
   
   END SUBROUTINE hamiltonian_deallocate

END MODULE hamiltonian_workspace_module
