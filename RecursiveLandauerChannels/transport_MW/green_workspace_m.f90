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
   MODULE  green_workspace_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI

   IMPLICIT NONE
   PRIVATE 
   SAVE



   ! Public
!
! Contains green functions data
!                                                   

   COMPLEX(dbl), ALLOCATABLE :: g_C(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: glesser(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: ggreater(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: gamma_L(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: gamma_R(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_L(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_R(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_L_lesser(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_R_lesser(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_L_greater(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_R_greater(:,:,:)

   COMPLEX(dbl), ALLOCATABLE :: sigma_C_r(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_C_lesser(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_C_greater(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_ep_lesser(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_ep_greater(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_ep_r(:,:,:)

   COMPLEX(dbl), ALLOCATABLE :: sigma_ee_contact_lesser(:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_ee_contact_r(:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_ee_mol_r(:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_ee_mol_lesser(:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_ee_mol_corr(:)
   REAL(dbl), ALLOCATABLE :: A_C(:,:,:)
   REAL(dbl), ALLOCATABLE :: dipole_matrix(:,:)

   ! Private
   LOGICAL :: alloc = .FALSE.


   ! Public variables
   PUBLIC :: g_C
   PUBLIC :: glesser
   PUBLIC :: gamma_L
   PUBLIC :: gamma_R
   PUBLIC :: sigma_L
   PUBLIC :: sigma_R
   PUBLIC :: sigma_L_lesser
   PUBLIC :: sigma_R_lesser
   PUBLIC :: sigma_C_r
   PUBLIC :: sigma_C_lesser

   PUBLIC :: sigma_ep_lesser
   PUBLIC :: sigma_ep_r

   PUBLIC :: sigma_ee_contact_lesser
   PUBLIC :: sigma_ee_contact_r
   PUBLIC :: sigma_ee_mol_r
   PUBLIC :: sigma_ee_mol_lesser
   PUBLIC :: sigma_ee_mol_corr
   PUBLIC :: dipole_matrix
   PUBLIC :: A_C

   PUBLIC :: ggreater
   PUBLIC :: sigma_L_greater
   PUBLIC :: sigma_R_greater
   PUBLIC :: sigma_C_greater
   PUBLIC :: sigma_ep_greater
 
! PUBLIC ROUTINES
! 

   PUBLIC :: green_allocate
   PUBLIC :: green_deallocate
 

 CONTAINS 

!***********************************************
   SUBROUTINE green_allocate(ne, dimC)
    !***********************************************
   IMPLICIT NONE
      INTEGER, INTENT(in)   :: ne
      INTEGER, INTENT(in)   :: dimC
      ! Local
       CHARACTER(14)      :: subname="green_allocate"
       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( ne <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in ne definition", " ne = ", ne
            STOP
         ENDIF 

        IF( dimC <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in dimC definition", " dimC = ", dimC
            STOP
         ENDIF 

       ALLOCATE( g_C(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating g_C"
            STOP
       ENDIF 

     ALLOCATE( ggreater(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating ggreater"
            STOP
       ENDIF 

       ALLOCATE( glesser(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating glesser"
            STOP
       ENDIF 

        ALLOCATE( gamma_L(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating gamma_L"
            STOP
       ENDIF 

       ALLOCATE( gamma_R(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating gamma_R"
            STOP
       ENDIF 

       ALLOCATE( sigma_R(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_R"
            STOP
       ENDIF 

       ALLOCATE( sigma_L(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_L"
            STOP
       ENDIF 
      ALLOCATE( sigma_R_lesser(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_R"
            STOP
       ENDIF 

       ALLOCATE( sigma_L_lesser(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_L"
            STOP
       ENDIF 
  
 
    ALLOCATE( sigma_R_greater(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_R"
            STOP
       ENDIF 

       ALLOCATE( sigma_L_greater(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_L"
            STOP
       ENDIF 
  
       
       ALLOCATE(sigma_C_r(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_C_r "
            STOP
       ENDIF 
      ALLOCATE(sigma_C_lesser(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_C_lesser"
            STOP
       ENDIF 

      ALLOCATE(sigma_C_greater(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating sigma_C_greater"
            STOP
       ENDIF 


      ALLOCATE(  sigma_ep_greater(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating  sigma_ep_greater"
            STOP
       ENDIF 

   
      ALLOCATE(  sigma_ep_lesser(ne,dimC,dimC)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating  sigma_ep_lessserr"
            STOP
       ENDIF 
   ALLOCATE(   sigma_ep_r(ne,dimC,dimC) , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating  sigma_ep_r"
            STOP
       ENDIF 
      ALLOCATE(    sigma_ee_contact_lesser(dimC,dimC), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating   sigma_ee_contact_lesser"
            STOP
       ENDIF 
      ALLOCATE(  sigma_ee_contact_r(dimC,dimC) , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating  sigma_ee_contact_r"
            STOP
       ENDIF 
       ALLOCATE(    sigma_ee_mol_r(dimC,dimC), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating   sigma_ee_mol_r"
            STOP
       ENDIF 
      ALLOCATE(   sigma_ee_mol_lesser(dimC,dimC), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating  sigma_ee_mol_lesser"
            STOP
       ENDIF 
      ALLOCATE(   sigma_ee_mol_corr(dimC), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating  sigma_ee_mol_corr"
            STOP
       ENDIF 
      ALLOCATE(  A_C(ne,dimC,dimC) , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating  A_C"
            STOP
       ENDIF 
     ALLOCATE( dipole_matrix(dimC,dimC) , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating   dipole_matrix"
            STOP
       ENDIF 
       
       g_C(:,:,:) = CZERO
       glesser(:,:,:) = CZERO
       ggreater(:,:,:) = CZERO
       gamma_L(:,:,:) = CZERO
       gamma_R(:,:,:) = CZERO
       sigma_L(:,:,:) = CZERO
       sigma_R(:,:,:) = CZERO
       sigma_L_lesser(:,:,:) = CZERO
       sigma_R_lesser(:,:,:) = CZERO
       sigma_L_greater(:,:,:) = CZERO
       sigma_R_greater(:,:,:) = CZERO
       sigma_C_greater(:,:,:)  = CZERO
       sigma_ep_greater(:,:,:)  = CZERO
       sigma_C_r(:,:,:)  = CZERO
       sigma_C_lesser(:,:,:)  = CZERO
       sigma_ep_lesser(:,:,:)  = CZERO
       sigma_ep_r(:,:,:)  = CZERO
       sigma_ee_contact_lesser(:,:)  = CZERO
       sigma_ee_contact_r(:,:)  = CZERO
       sigma_ee_mol_r(:,:)  = CZERO
       sigma_ee_mol_lesser(:,:)  = CZERO
       sigma_ee_mol_corr(:)  = CZERO
       A_C(:,:,:)  = ZERO
       dipole_matrix(:,:)  = ZERO

       alloc = .TRUE.

   END SUBROUTINE green_allocate

!***********************************************
   SUBROUTINE green_deallocate
   !***********************************************
      CHARACTER(17)      :: subname="green_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED( g_C ) ) THEN
            DEALLOCATE( g_C , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating gC_ "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( glesser) ) THEN
            DEALLOCATE( glesser , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating glesser"
               STOP
           ENDIF 
       ENDIF

      IF ( ALLOCATED( ggreater) ) THEN
            DEALLOCATE( ggreater , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating ggreater"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( sigma_L ) ) THEN
            DEALLOCATE( sigma_L, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_L"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( sigma_R ) ) THEN
            DEALLOCATE( sigma_R , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_R "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( gamma_L ) ) THEN
            DEALLOCATE( gamma_L, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating gamma_L"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( gamma_R ) ) THEN
            DEALLOCATE( gamma_R, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating gamma_R"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( dipole_matrix ) ) THEN
            DEALLOCATE( dipole_matrix, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating dipole_matrix"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( sigma_L_lesser ) ) THEN
            DEALLOCATE(sigma_L_lesser, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_L_lesser"
               STOP
           ENDIF 
       ENDIF

       

       IF ( ALLOCATED( sigma_R_greater  ) ) THEN
            DEALLOCATE(  sigma_R_greater, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_R_grer"
               STOP
           ENDIF 
       ENDIF
      
       IF ( ALLOCATED( sigma_L_greater) ) THEN
            DEALLOCATE(sigma_L_greater, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_L_ger"
               STOP
           ENDIF 
       ENDIF

       

       IF ( ALLOCATED( sigma_R_lesser  ) ) THEN
            DEALLOCATE(  sigma_R_lesser, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_R_lesser"
               STOP
           ENDIF 
       ENDIF
 
       IF ( ALLOCATED( sigma_C_r  ) ) THEN
            DEALLOCATE( sigma_C_r, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_C_r"
               STOP
           ENDIF 
       ENDIF
      

       IF ( ALLOCATED(  sigma_C_lesser ) ) THEN
            DEALLOCATE( sigma_C_lesser, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_C_lesser"
               STOP
           ENDIF 
       ENDIF
      IF ( ALLOCATED(  sigma_C_greater ) ) THEN
            DEALLOCATE( sigma_C_greater, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_Cgreaterer"
               STOP
           ENDIF 
       ENDIF
       

       IF ( ALLOCATED( sigma_ep_lesser ) ) THEN
            DEALLOCATE(sigma_ep_lesser, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_ep_lesser"
               STOP
           ENDIF 
       ENDIF
       
    IF ( ALLOCATED( sigma_ep_greater ) ) THEN
            DEALLOCATE(sigma_ep_greater, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_ep_greater"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(  sigma_ep_r ) ) THEN
            DEALLOCATE( sigma_ep_r, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_ep_r"
               STOP
           ENDIF 
       ENDIF
      

       IF ( ALLOCATED( sigma_ee_contact_lesser ) ) THEN
            DEALLOCATE( sigma_ee_contact_lesser, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_ee_contact_lesser"
               STOP
           ENDIF 
       ENDIF
      

       IF ( ALLOCATED( sigma_ee_contact_r ) ) THEN
            DEALLOCATE(sigma_ee_contact_r, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_ee_contact_r"
               STOP
           ENDIF 
       ENDIF
       

       IF ( ALLOCATED( sigma_ee_mol_r ) ) THEN
            DEALLOCATE(sigma_ee_mol_r, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating sigma_ee_mol_r"
               STOP
           ENDIF 
       ENDIF
       

       IF ( ALLOCATED(  sigma_ee_mol_lesser ) ) THEN
            DEALLOCATE( sigma_ee_mol_lesser, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_ee_mol_lesser"
               STOP
           ENDIF 
       ENDIF
      

       IF ( ALLOCATED( sigma_ee_mol_corr  ) ) THEN
            DEALLOCATE( sigma_ee_mol_corr, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  sigma_ee_mol_corr"
               STOP
           ENDIF 
       ENDIF
      

       IF ( ALLOCATED( A_C  ) ) THEN
            DEALLOCATE( A_C, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  A_C"
               STOP
           ENDIF 
       ENDIF
      
       

       alloc = .FALSE.
   
   END SUBROUTINE green_deallocate

END MODULE green_workspace_module
