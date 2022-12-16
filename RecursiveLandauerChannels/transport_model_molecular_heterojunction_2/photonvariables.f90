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
   MODULE photon_variable_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE 
   PRIVATE 
   SAVE
!


   ! Public
   INTEGER, ALLOCATABLE   :: ihomo(:)     ! index of homo states
   INTEGER, ALLOCATABLE   :: ilumo(:)     ! index of lumo states
   REAL(dbl), ALLOCATABLE :: ehomo(:)     ! energy values of the HOMO
   REAL(dbl), ALLOCATABLE :: elumo(:)     ! energy values of the HOMO
   REAL(dbl), ALLOCATABLE :: photoncouplingconstant(:) ! Contains energy of the transmission homo -  lumo
   REAL(dbl)              :: Field_Intensity ! Vo (input)
   REAL(dbl)              :: omegamax ! Vo (input)
   REAL(dbl)              :: FGR_broadening ! broadening for Fermi golden rule (lorentzian) from input

   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: ihomo
   PUBLIC                 :: ilumo
   PUBLIC                 :: ehomo
   PUBLIC                 :: elumo
   PUBLIC                 :: Field_Intensity 
   PUBLIC                 :: photoncouplingconstant
   PUBLIC                 :: FGR_broadening
   PUBLIC                 :: omegamax
   ! Public routines:
   PUBLIC                 :: photonvariable_allocate
   PUBLIC                 :: photonvariable_deallocate
   PUBLIC                 :: photonvariable_print

   CONTAINS 
!***********************************************
   SUBROUTINE photonvariable_allocate(ncrosssection,dimC)
   !***********************************************
   IMPLICIT NONE
       INTEGER, INTENT(in)  :: ncrosssection

       INTEGER, INTENT(in)  :: dimC
       CHARACTER(23)      :: subname="photonvariable_allocate"
      
       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 


        IF(ncrosssection <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nalpha definition", " nalpha = ", ncrosssection
            STOP
         ENDIF 


   

        IF( dimC <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nalpha definition", " dimC = ", dimC
            STOP
         ENDIF 

     ALLOCATE ( ehomo(dimC), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating ehomo"
            STOP
         ENDIF 
     ALLOCATE ( elumo(dimC), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating elumo"
            STOP
         ENDIF 
     ALLOCATE ( ihomo(dimC), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating ihomo"
            STOP
         ENDIF 
     ALLOCATE ( ilumo(dimC), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating ilumo"
            STOP
         ENDIF 
    ALLOCATE ( photoncouplingconstant(ncrosssection), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating ilumo"
            STOP
         ENDIF 

 !    ALLOCATE ( ehomo(nalpha, nbias, (INT(dimC/2))), STAT=ierr )
 !        IF( ierr /=0 ) THEN 
 !           PRINT*, "Error in Routine ",  TRIM(subname) 
 !           PRINT*, "Error allocating ehomo"
 !           STOP
 !        ENDIF 
 !    ALLOCATE ( elumo(nalpha, nbias, (INT(dimC/2))), STAT=ierr )
 !        IF( ierr /=0 ) THEN 
 !           PRINT*, "Error in Routine ",  TRIM(subname) 
 !           PRINT*, "Error allocating elumo"
 !           STOP
 !        ENDIF 
 !    ALLOCATE ( ihomo(nalpha, nbias, (INT(dimC/2))), STAT=ierr )
 !        IF( ierr /=0 ) THEN 
 !           PRINT*, "Error in Routine ",  TRIM(subname) 
 !           PRINT*, "Error allocating ihomo"
 !           STOP
 !        ENDIF 
 !    ALLOCATE ( ilumo(nalpha, nbias, (INT(dimC/2))), STAT=ierr )
 !        IF( ierr /=0 ) THEN 
 !           PRINT*, "Error in Routine ",  TRIM(subname) 
 !           PRINT*, "Error allocating ilumo"
 !           STOP
 !        ENDIF 



       alloc = .TRUE.


  END SUBROUTINE photonvariable_allocate
!***********************************************
   SUBROUTINE photonvariable_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(25)      :: subname="photonvariable_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(photoncouplingconstant) ) THEN
            DEALLOCATE(photoncouplingconstant, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating omegamax"
               STOP
           ENDIF 
       ENDIF



    IF ( ALLOCATED(ehomo) ) THEN
            DEALLOCATE(ehomo, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine",  TRIM(subname) 
               PRINT*, "Error deallocating ehomo"
               STOP
           ENDIF 
       ENDIF


    IF ( ALLOCATED(elumo) ) THEN
            DEALLOCATE(elumo, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating elumo"
               STOP
           ENDIF 
       ENDIF


    IF ( ALLOCATED( ihomo) ) THEN
            DEALLOCATE( ihomo, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating ihomo"
               STOP
           ENDIF 
       ENDIF


    IF ( ALLOCATED(ilumo) ) THEN
            DEALLOCATE(ilumo, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating ilumo"
               STOP
           ENDIF 
       ENDIF


      alloc = .FALSE.


   END SUBROUTINE photonvariable_deallocate

!***********************************************
   SUBROUTINE photonvariable_print(Fil, nalpha,nbias,dimC, alphagrid, biasgrid)
   !***********************************************
   IMPLICIT NONE
       INTEGER, INTENT(in)  :: Fil
       INTEGER, INTENT(in)  :: nalpha
       INTEGER, INTENT(in)  :: nbias
       INTEGER, INTENT(in)  :: dimC
       REAL(dbl), INTENT(in)  :: alphagrid(nalpha)
       REAL(dbl), INTENT(in)  :: biasgrid(nbias)

       CHARACTER(20)      :: subname="photonvariable_print"
       INTEGER  :: ierr, ialpha, ibias, i_dim, i_dim2
       !
        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 

       IF( nalpha <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nalpha definition", " nalpha = ", nalpha
            STOP
         ENDIF 


        IF( nbias <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nalpha definition", " nbias = ", nbias
            STOP
         ENDIF 

        IF( dimC <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nalpha definition", " dimC = ", dimC
            STOP
         ENDIF 




	 OPEN ( Fil, FILE='ehomo.dat', FORM='formatted' )

	   DO ialpha = 1, nalpha
	      DO ibias = 1, nbias
	           WRITE (Fil, '(30(f15.9))' )  alphagrid(ialpha) , biasgrid(ibias) , ( ehomo(i_dim2) , i_dim2=1, dimC )
	      ENDDO
              WRITE (Fil, '(2(e15.9))' )
	   ENDDO
        CLOSE( Fil )

	 OPEN ( Fil, FILE='elumo.dat', FORM='formatted' )

	   DO ialpha = 1, nalpha
	      DO ibias = 1, nbias
	           WRITE (Fil, '(30(f15.9))' )  alphagrid(ialpha) , biasgrid(ibias) , ( elumo( i_dim2) , i_dim2=1, dimC)
	      ENDDO
              WRITE (Fil, '(2(e15.9))' )
	   ENDDO
        CLOSE( Fil )
	 OPEN ( Fil, FILE='ihomo.dat', FORM='formatted' )

	   DO ialpha = 1, nalpha
	      DO ibias = 1, nbias
	           WRITE (Fil, '(30(I10))' )  ialpha ,ibias , (  ihomo(i_dim), i_dim=1, dimC )
	      ENDDO
              WRITE (Fil, '(2(e15.9))' )
	   ENDDO
        CLOSE( Fil )

	 OPEN ( Fil, FILE='ilumo.dat', FORM='formatted' )

	   DO ialpha = 1, nalpha
	      DO ibias = 1, nbias
	           WRITE (Fil, '(30(I10))' )  ialpha , ibias , (  ilumo( i_dim), i_dim=1, dimC )
	      ENDDO
              WRITE (Fil, '(2(e15.9))' )
	   ENDDO
        CLOSE( Fil )

	 OPEN ( Fil, FILE='Crosssectiongrid.dat', FORM='formatted' )

	   DO ibias = 1, nbias
	           WRITE (Fil, '(2(f15.9))' )  biasgrid(ibias) , photoncouplingconstant(ibias)
	   ENDDO
        CLOSE( Fil )


   END SUBROUTINE photonvariable_print

  END  MODULE photon_variable_module

 
