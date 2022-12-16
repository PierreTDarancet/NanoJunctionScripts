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
   MODULE hamiltonian_construction_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE
   PRIVATE 
   SAVE



   ! Public
  INTEGER  ::   dimL_ref=1
  INTEGER  ::   dimR_ref=1
  INTEGER   ::   dimC_ref
  CHARACTER(100) :: HC_file
  CHARACTER(100) :: HL_file
  CHARACTER(100) :: HR_file
   ! Private

!
! Contains reference Hamiltonian data
! 
  
   COMPLEX(dbl), ALLOCATABLE :: h00_C_ref(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h_CR_ref(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h_LC_ref(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h00_R_ref(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h01_R_ref(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h00_L_ref(:,:)
   COMPLEX(dbl), ALLOCATABLE :: h01_L_ref(:,:)

   LOGICAL :: alloc = .FALSE.
   LOGICAL :: init  = .FALSE.
  CHARACTER(31) :: modulename="hamiltonian_construction_module"
   ! Public variables
   PUBLIC ::  dimC_ref
   PUBLIC ::  dimR_ref
   PUBLIC ::  dimL_ref
   PUBLIC ::  HC_file
   PUBLIC ::  HL_file
   PUBLIC ::  HR_file

!
! PUBLIC ROUTINES
! 

   PUBLIC :: hamiltonian_reference_allocate
   PUBLIC :: hamiltonian_reference_deallocate
   PUBLIC :: Molecular_hamiltonian_init
   PUBLIC :: hamiltonian_reference_init

 CONTAINS 

!***********************************************
   SUBROUTINE hamiltonian_reference_allocate
    !***********************************************
   IMPLICIT NONE
       CHARACTER(30)      ::  subname="hamiltonian_reference_allocate"
       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( dimC_ref <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in dimC definition", " dimC_ref = ", dimC_ref
            STOP
         ENDIF 

        IF( dimL_ref <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in dimL definition", " dimL_ref = ", dimL_ref
            STOP
         ENDIF 

        IF( dimR_ref <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in dimR definition", " dimR_ref = ", dimR_ref
            STOP
         ENDIF 

       ALLOCATE( h00_C_ref(dimC_ref,dimC_ref)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h00_C_ref"
            STOP
       ENDIF 

       ALLOCATE( h00_R_ref(dimR_ref,dimR_ref)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h00_R_ref"
            STOP
       ENDIF 
       ALLOCATE( h00_L_ref(dimL_ref,dimL_ref)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h00_L_ref"
            STOP
       ENDIF 
       ALLOCATE( h01_L_ref(dimL_ref,dimL_ref)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h01_L_ref"
            STOP
       ENDIF 
       ALLOCATE( h01_R_ref(dimR_ref,dimR_ref)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h01_R_ref"
            STOP
       ENDIF 
       ALLOCATE( h_LC_ref(dimL_ref,dimC_ref)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h_LC_ref"
            STOP
       ENDIF 
       ALLOCATE( h_CR_ref(dimC_ref,dimR_ref)  , STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating h_CR_ref"
            STOP
       ENDIF 

       h_CR_ref(:,:)  = CZERO
       h_LC_ref(:,:)  = CZERO 
       h00_L_ref(:,:) = CZERO
       h00_R_ref(:,:) = CZERO 
       h01_L_ref(:,:) = CZERO 
       h01_R_ref(:,:) = CZERO 
       h00_C_ref(:,:) = CZERO
       alloc = .TRUE.

   END SUBROUTINE hamiltonian_reference_allocate

!***********************************************
   SUBROUTINE hamiltonian_reference_deallocate
   !***********************************************
      CHARACTER(32)      :: subname="hamiltonian_reference_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED( h00_C_ref ) ) THEN
            DEALLOCATE( h00_C_ref , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h00_C_ref "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( h_CR_ref ) ) THEN
            DEALLOCATE( h_CR_ref , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h_CR_ref"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( h_LC_ref   ) ) THEN
            DEALLOCATE( h_LC_ref   , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  h_LC_ref "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(h00_L_ref ) ) THEN
            DEALLOCATE( h00_L_ref, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h00_L_ref"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(h00_R_ref  ) ) THEN
            DEALLOCATE( h00_R_ref , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h00_R_ref "
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(h01_L_ref ) ) THEN
            DEALLOCATE( h01_L_ref, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h01_L_ref"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED( h01_R_ref ) ) THEN
            DEALLOCATE( h01_R_ref , STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating h01_R_ref "
               STOP
           ENDIF 
       ENDIF

       alloc = .FALSE.
   
   END SUBROUTINE hamiltonian_reference_deallocate


!***********************************************
   SUBROUTINE Molecular_hamiltonian_init(h00_C, h_CR, h_LC, h_RC, h_CL,  h00_L, h01_L, h00_R, h01_R, bias, gamma_ratio, gap, dimCaux, debug_level)  
   !***********************************************
           !        CALL Molecular_hamiltonian_init(h00_C, h_CR, h_LC, h_RC, h_CL, h00_L, h01_L, h00_R, h01_R, biasaux, gammagrid(igamma), gapgrid(igap), dimC, debug_level)

      INTEGER, PARAMETER  ::  dimLaux=1
      INTEGER, PARAMETER  ::  dimRaux=1
      INTEGER, INTENT(in)  ::  dimCaux
      COMPLEX(dbl), INTENT(out) :: h00_C(dimCaux,dimCaux) 
      COMPLEX(dbl), INTENT(out) :: h_CR(dimCaux,dimRaux) 
      COMPLEX(dbl), INTENT(out) :: h_RC(dimRaux,dimCaux) 
      COMPLEX(dbl), INTENT(out) :: h_LC(dimLaux,dimCaux) 
      COMPLEX(dbl), INTENT(out) :: h_CL(dimCaux,dimLaux)
      COMPLEX(dbl), INTENT(out) :: h00_L(dimLaux,dimLaux) 
      COMPLEX(dbl), INTENT(out) :: h01_L(dimLaux,dimLaux) 
      COMPLEX(dbl), INTENT(out) :: h00_R(dimRaux,dimRaux) 
      COMPLEX(dbl), INTENT(out) :: h01_R(dimRaux,dimRaux) 

      REAL(dbl), INTENT(in)  :: bias !Notused
      REAL(dbl), INTENT(in)  :: gamma_ratio
      REAL(dbl), INTENT(in)  :: gap
      REAL(dbl) :: realaux

      CHARACTER(26)      :: subname="Molecular_hamiltonian_init"
      INTEGER :: ierr, irow, icol, icountl, icounth, i_dim
      INTEGER, INTENT(in)  ::  debug_level


        !
        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Hamiltonian not allocated"
            STOP
         ENDIF 

       IF( .NOT. init ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Reference Hamiltonian not initialized"
            STOP
         ENDIF 

      realaux = h_CR_ref(1,1)**2 +  h_LC_ref(1,1)**2
      
      h_CR(1,1)  =  h_CR_ref(1,1)             !sqrt( realaux / (ONE +  gamma_ratio) )
      h_LC(1,1)  =  h_LC_ref(1,1)          !sqrt( realaux / (ONE +  (ONE/gamma_ratio) ) )
      h_CR(2,1)  =  h_CR_ref(2,1)          !sqrt( realaux / (ONE +  (ONE/gamma_ratio) ) )
      h_LC(1,2)  =  h_LC_ref(1,2)          !sqrt( realaux / (ONE + gamma_ratio) )

      h_RC(:,:) = CONJG( TRANSPOSE( h_CR(:,:) ) )
      h_CL(:,:) = CONJG( TRANSPOSE( h_LC(:,:) ) )


      h00_L(:,:) =  h00_L_ref(1:dimLaux,1:dimLaux) 
      h01_L(:,:) =  h01_L_ref(1:dimLaux,1:dimLaux) 
        ! 
      h00_R(:,:) =  h00_R_ref(1:dimRaux,1:dimRaux) 
      h01_R(:,:) =  h01_R_ref(1:dimRaux,1:dimRaux) 

       
      h00_C(1:dimCaux,1:dimCaux)  =   h00_C_ref(1:dimCaux,1:dimCaux)  
  !!!!!!!!!!!!!!!!!!!!WARNING CHANGE!!!!!!
      h00_C(1,1)=  h00_C(1,1) +  gap
      h00_C(2,2)=  h00_C(2,2) +  gap
   !   h00_C(2,2)  =  h00_C(1,1) + gap    CHANGE THAT BACK Trick for quick Plot restore for having the real gap loop
   END SUBROUTINE Molecular_hamiltonian_init





!***********************************************
   SUBROUTINE hamiltonian_reference_init(debug_level)
   !***********************************************
      INTEGER, INTENT(in)  ::  debug_level
      INTEGER, PARAMETER  :: Fil1=201
      CHARACTER(26)      :: subname="hamiltonian_reference_init"
      INTEGER :: ierr, dim1, dim2, irow, icol
      REAL(dbl) :: rowaux(dimC_ref)
      REAL(dbl) :: rowauxl(dimL_ref)
      REAL(dbl) :: rowauxr(dimR_ref)
      INTEGER, PARAMETER :: Hcfil=101 , Hlfil=102, Hrfil=103
      CHARACTER(100)  :: chr1

        !
        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Hamiltonian not allocated"
            STOP
         ENDIF 

       IF( init ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Hamiltonian already initialized"
            STOP
         ENDIF 


   PRINT*, "... READING Hamiltonian Input Files"
   !
   !   
   !
   !  reading H_C
   !
   !
   !

   !
   PRINT*, "     ... READING HC_file:",   TRIM(HC_file)
   OPEN(Hcfil, FILE=TRIM(HC_file),FORM='formatted')
   PRINT*, "          ...... READING H00C"
         ! Tests...
         READ(Hcfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<H00_C>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <H00_C>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hcfil,"(a3,I5,I5)") chr1, dim1, dim2
           IF ( TRIM(chr1) /= "dim" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dim of <H00_C>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
           IF ( dim1 /= dimC_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H00_C>", " Input = ", dim1
            STOP
         ENDIF 
           IF ( dim2 /= dimC_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H00_C>", " Input = ", dim2
            STOP
         ENDIF 
         READ(Hcfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<BEGIN>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <BEGIN>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         ! Reading 

         DO irow=1, dim1
             !
             READ(Hcfil,*) (rowaux(icol) , icol=1,dim2) 
             DO icol=1, dim2
                !
                h00_C_ref(irow,icol)= CMPLX (rowaux(icol), ZERO)
                !
             ENDDO
             !
         ENDDO 
         READ(Hcfil,"(a5)") chr1
           IF ( TRIM(chr1) /= "<END>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <END>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hcfil,"(a5)") chr1 ! Line just empty


   PRINT*, "          ...... READING HLC"
        ! Tests...
         READ(Hcfil,"(a6)") chr1
           IF ( TRIM(chr1) /= "<H_LC>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <H_LC>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hcfil,"(a3,I5,I5)") chr1, dim1, dim2
           IF ( TRIM(chr1) /= "dim" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dim of <H_LC>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
           IF ( dim1 /= dimL_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H_LC>", " Input = ", dim1
            STOP
         ENDIF 
           IF ( dim2 /= dimC_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H_LC>", " Input = ", dim2
            STOP
         ENDIF 
         READ(Hcfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<BEGIN>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <BEGIN>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         ! Reading 

         DO irow=1, dim1
             !
             READ(Hcfil,*) (rowaux(icol) , icol=1,dim2) 
             DO icol=1, dim2
                !
                h_LC_ref(irow,icol)= CMPLX (rowaux(icol), ZERO)
                !
             ENDDO
             !
         ENDDO 
         READ(Hcfil,"(a5)") chr1
           IF ( TRIM(chr1) /= "<END>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <END>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hcfil,"(a5)") chr1 ! Line just empty
   PRINT*, "          ...... READING HCR"
        ! Tests...
         READ(Hcfil,"(a6)") chr1
           IF ( TRIM(chr1) /= "<H_CR>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <H_CR>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hcfil,"(a3,I5,I5)") chr1, dim1, dim2
           IF ( TRIM(chr1) /= "dim" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dim of <H_CR>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
           IF ( dim1 /= dimC_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H_CR>", " Input = ", dim1
            STOP
         ENDIF 
           IF ( dim2 /= dimR_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H_CR>", " Input = ", dim2
            STOP
         ENDIF 
         READ(Hcfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<BEGIN>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <BEGIN>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         !
         ! Reading 
         !
         DO irow=1, dim1
             !
             READ(Hcfil,*) (rowauxr(icol) , icol=1,dim2) 
             DO icol=1, dim2
                !
                h_CR_ref(irow,icol)= CMPLX (rowauxr(icol), ZERO)
                !
             ENDDO
             !
         ENDDO 
         READ(Hcfil,"(a5)") chr1
           IF ( TRIM(chr1) /= "<END>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <END>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         !
   CLOSE(Hcfil)
   PRINT*, "     ... READING HC_file finished"
   !
   !   
   !
   !  End of reading for H_C
   !
   !
   !


   !
   !   
   !
   !  read H_l
   !
   !
   !

   PRINT*, "     ... READING HL_file:",   TRIM(HL_file)
   OPEN(Hlfil, FILE=TRIM(HL_file),FORM='formatted')
   PRINT*, "          ...... READING H00L"
         ! Tests...
         READ(Hlfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<H00_L>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <H00_L>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hlfil,"(a3,I5,I5)") chr1, dim1, dim2
           IF ( TRIM(chr1) /= "dim" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dim of <H00_L>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
           IF ( dim1 /= dimL_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H00_L>", " Input = ", dim1
            STOP
         ENDIF 
           IF ( dim2 /= dimL_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H00_L>", " Input = ", dim2
            STOP
         ENDIF 
         READ(Hlfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<BEGIN>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <BEGIN>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         ! Reading 

         DO irow=1, dim1
             !
             READ(Hlfil,*) (rowauxl(icol) , icol=1,dim2) 
             DO icol=1, dim2
                !
                h00_L_ref(irow,icol)= CMPLX (rowauxl(icol), ZERO)
                !
             ENDDO
             !
         ENDDO 
         READ(Hlfil,"(a5)") chr1
           IF ( TRIM(chr1) /= "<END>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <END>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hlfil,"(a5)") chr1 ! Line just empty

   PRINT*, "          ...... READING H01L"
         ! Tests...
         READ(Hlfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<H01_L>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <H01_L>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hlfil,"(a3,I5,I5)") chr1, dim1, dim2
           IF ( TRIM(chr1) /= "dim" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dim of <H01_L>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
           IF ( dim1 /= dimL_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H01_L>", " Input = ", dim1
            STOP
         ENDIF 
           IF ( dim2 /= dimL_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H01_L>", " Input = ", dim2
            STOP
         ENDIF 
         READ(Hlfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<BEGIN>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <BEGIN>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         ! Reading 

         DO irow=1, dim1
             !
             READ(Hlfil,"(16F5.10)") (rowauxl(icol) , icol=1,dim2) 
             DO icol=1, dim2
                !
                h01_L_ref(irow,icol)= CMPLX (rowauxl(icol), ZERO)
                !
             ENDDO
             !
         ENDDO 
         READ(Hlfil,"(a5)") chr1
           IF ( TRIM(chr1) /= "<END>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <END>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
   CLOSE(Hlfil)
   PRINT*, "     ... READING HL_file finished"

   !
   !   
   !
   !  End of reading for H_l
   !
   !
   !


   !
   !   
   !
   !  read H_r
   !
   !
   !

   PRINT*, "     ... READING HR_file:",   TRIM(HR_file)
   OPEN(Hrfil, FILE=TRIM(HR_file),FORM='formatted')
   PRINT*, "          ...... READING H00R"
         ! Tests...
         READ(Hrfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<H00_R>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <H00_R>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hrfil,"(a3,I5,I5)") chr1, dim1, dim2
           IF ( TRIM(chr1) /= "dim" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dim of <H00_R>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
           IF ( dim1 /= dimR_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H00_R>", " Input = ", dim1
            STOP
         ENDIF 
           IF ( dim2 /= dimR_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H00_R>", " Input = ", dim2
            STOP
         ENDIF 
         READ(Hrfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<BEGIN>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <BEGIN>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         ! Reading 

         DO irow=1, dim1
             !
             READ(Hrfil,*) (rowauxr(icol) , icol=1,dim2) 
             DO icol=1, dim2
                !
                h00_R_ref(irow,icol)= CMPLX (rowauxr(icol), ZERO)
                !
             ENDDO
             !
         ENDDO 
         READ(Hrfil,"(a5)") chr1
           IF ( TRIM(chr1) /= "<END>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <END>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hrfil,"(a5)") chr1 ! Line just empty

   PRINT*, "          ...... READING H01R"
         ! Tests...
         READ(Hrfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<H01_R>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <H01_R>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         READ(Hrfil,"(a3,I5,I5)") chr1, dim1, dim2
           IF ( TRIM(chr1) /= "dim" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dim of <H01_R>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
           IF ( dim1 /= dimR_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H01_R>", " Input = ", dim1
            STOP
         ENDIF 
           IF ( dim2 /= dimR_ref ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Error in dim of <H01_R>", " Input = ", dim2
            STOP
         ENDIF 
         READ(Hrfil,"(a7)") chr1
           IF ( TRIM(chr1) /= "<BEGIN>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <BEGIN>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
         ! Reading 

         DO irow=1, dim1
             !
             READ(Hrfil,*) (rowauxr(icol) , icol=1,dim2) 
             DO icol=1, dim2
                !
                h01_R_ref(irow,icol)= CMPLX (rowauxr(icol), ZERO)
                !
             ENDDO
             !
         ENDDO 
         READ(Hrfil,"(a5)") chr1
           IF ( TRIM(chr1) /= "<END>" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for <END>", " Input = ", TRIM(chr1)
            STOP
         ENDIF 
   CLOSE(Hrfil)
   PRINT*, "     ... READING HR_file finished"
   !


  IF ( debug_level > 5) THEN
          !PRINT READHAMILTONIAN 
         OPEN ( Fil1, FILE='READhamiltonian.dat', FORM='formatted' )
         WRITE(Fil1,*)  "<H00_C>" 
         WRITE(Fil1,*)  "<BEGIN>" 
         DO irow=1, dimC_ref
             WRITE (Fil1, * )  ( h00_C_ref(irow,icol) , icol = 1, dimC_ref) 
         ENDDO 
         WRITE( Fil1 ,*)  "<END>"
         !
         WRITE(Fil1,*)  "<H_LC>" 
         WRITE(Fil1,*)  "<BEGIN>" 
         DO irow=1, dimL_ref
             WRITE (Fil1, * )  ( h_LC_ref(irow,icol) , icol = 1, dimC_ref) 
         ENDDO 
         WRITE( Fil1 ,*)  "<END>"
         !
         WRITE(Fil1,*)  "<H_CR>" 
         WRITE(Fil1,*)  "<BEGIN>" 
         DO irow=1, dimC_ref
             WRITE (Fil1, * )  ( h_CR_ref(irow,icol) , icol = 1, dimR_ref) 
         ENDDO 
         WRITE( Fil1 ,*)  "<END>"
         !
         WRITE(Fil1,*)  "<H00_L>" 
         WRITE(Fil1,*)  "<BEGIN>" 
         DO irow=1, dimL_ref
             WRITE (Fil1, * )  ( h00_L_ref(irow,icol) , icol = 1, dimL_ref) 
         ENDDO 
         WRITE( Fil1 ,*)  "<END>"
         !
         WRITE(Fil1,*)  "<H01_L>" 
         WRITE(Fil1,*)  "<BEGIN>" 
         DO irow=1, dimL_ref
             WRITE (Fil1, * )  ( h01_L_ref(irow,icol) , icol = 1, dimL_ref) 
         ENDDO 
         WRITE( Fil1 ,*)  "<END>"
         !
         WRITE(Fil1,*)  "<H00_R>" 
         WRITE(Fil1,*)  "<BEGIN>" 
         DO irow=1, dimR_ref
             WRITE (Fil1, * )  ( h00_R_ref(irow,icol) , icol = 1, dimR_ref) 
         ENDDO 
         WRITE( Fil1 ,*)  "<END>"
         !
         WRITE(Fil1,*)  "<H01_R>" 
         WRITE(Fil1,*)  "<BEGIN>" 
         DO irow=1, dimR_ref
             WRITE (Fil1, * )  ( h01_R_ref(irow,icol) , icol = 1, dimR_ref) 
         ENDDO 
         WRITE( Fil1 ,*)  "<END>"

         CLOSE( Fil1 )
   ENDIF

   PRINT*, "... END of READING HAMILTONIAN DATA Files"

    init=.TRUE.

   END SUBROUTINE hamiltonian_reference_init

END MODULE hamiltonian_construction_module
