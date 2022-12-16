!
! Copyright (C) 2011 Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET

!***********************************************
   PROGRAM conductor
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok

   USE T_control_module,                   ONLY : debug_level,  gamma_big,  gamma_small


   USE T_biasgrid_module,                      ONLY : biasgrid_allocate, biasgrid_deallocate, biasgrid_init, nbias, biasgrid !ok
   USE T_gapgrid_module,                   ONLY : gapgrid_allocate, gapgrid_deallocate, gapgrid_init, ngap, gapgrid, gapmax, gapmin !ok
   USE T_alphagrid_module,                 ONLY : alphagrid_allocate, alphagrid_deallocate, alphagrid_init, nalpha, alphagrid !ok


   USE T_input_module,                     ONLY : input_manager








   IMPLICIT NONE

   !
   ! local variables
   !

   INTEGER, PARAMETER :: Fil1 = 11 
   INTEGER, PARAMETER :: Fil2 = 12 
   INTEGER, PARAMETER :: Fil3 = 13
   INTEGER, PARAMETER :: Fil4 = 14
   INTEGER, PARAMETER :: Fil5 = 15
   INTEGER, PARAMETER :: Fil6 = 16
   INTEGER, PARAMETER :: Fil7 = 17
   REAL(dbl), PARAMETER :: G0 = 77.480917*ONE*1000.0  ! NanoS

 
   CHARACTER(100)  :: outputfile_name, chr
   INTEGER          :: ie, ierr, ncount, igap, ialpha, ibias
   !

! Contains transport fundamental variables
! 

   REAL(dbl), ALLOCATABLE              :: CurrentPositive(:,:,:) ! NanoS
   REAL(dbl), ALLOCATABLE              :: CurrentNegative(:,:,:) ! NanoS
   REAL(dbl), ALLOCATABLE              :: RR(:,:,:)
   REAL(dbl), ALLOCATABLE              :: PrefactorPositive(:,:,:)
   REAL(dbl), ALLOCATABLE              :: PrefactorNegative(:,:,:)
   REAL(dbl)             :: Ene_resonance
!
!------------------------------
! main body
!------------------------------
!   !
   ! read input file
   !


   CALL input_manager()
   !
   ! init
   !
      !
  


     ALLOCATE ( CurrentPositive(nalpha, ngap,nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating current"
            STOP
         ENDIF 

     ALLOCATE ( CurrentNegative(nalpha, ngap,nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating current"
            STOP
         ENDIF 


     ALLOCATE ( RR(nalpha, ngap,nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating rectification ratio"
            STOP
         ENDIF 
    ALLOCATE ( PrefactorPositive(nalpha, ngap,nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating rectification ratio2"
            STOP
         ENDIF 
    ALLOCATE ( PrefactorNegative(nalpha, ngap,nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating rectification ratio2"
            STOP
         ENDIF 

       RR(:,:,:) = ZERO
       CurrentPositive(:,:,:) = ZERO
       CurrentNegative(:,:,:) = ZERO
       PrefactorPositive(:,:,:) = ZERO
       PrefactorNegative(:,:,:) = ZERO

    PRINT*, 'Beginning grid allocations'

   !
    CALL biasgrid_allocate()
    CALL alphagrid_allocate()
    CALL gapgrid_allocate()
  

    PRINT*, 'Beginning Main Loop'

    CALL alphagrid_init()
    CALL gapgrid_init()
    CALL biasgrid_init()
   !

   DO ibias = 1,nbias 


   alpha_loop: &
   DO ialpha=1, nalpha

	   gap_loop: &
	   DO igap=1, ngap

                   Ene_resonance = (biasgrid(ibias)/2.0) * alphagrid(ialpha) + gapgrid(igap)

                   
		   !Forward bias
		   CurrentPositive(ialpha,igap,ibias)=ZERO
                   PrefactorPositive(ialpha,igap,ibias) = 2.0*PI*  gamma_small / (ONE +  (gamma_small / gamma_big))
                   CurrentPositive(ialpha,igap,ibias)=G0 * (PrefactorPositive(ialpha,igap,ibias) /PI) * (  ATAN(( (biasgrid(ibias)/2.0) + Ene_resonance )/( (ONE +  (gamma_small / gamma_big))*gamma_big/2.0 ) ) +  ATAN(( (biasgrid(ibias)/2.0) - Ene_resonance )/( (ONE +  (gamma_small / gamma_big))*gamma_big/2.0 ) ) )

	           !Reverse bias
                   Ene_resonance = (-biasgrid(ibias)/2.0) * alphagrid(ialpha) + gapgrid(igap)

                   PrefactorNegative(ialpha,igap,ibias) = 2.0*PI*  gamma_small / (ONE +  (gamma_small / gamma_big))
       		   CurrentNegative(ialpha,igap,ibias)=ZERO
                   CurrentNegative(ialpha,igap,ibias)=G0 * (PrefactorNegative(ialpha,igap,ibias) /PI) * (  ATAN(( (-biasgrid(ibias)/2.0) + Ene_resonance )/( (ONE +  (gamma_small / gamma_big))*gamma_big/2.0 ) ) +  ATAN(( (-biasgrid(ibias)/2.0) - Ene_resonance )/( (ONE +  (gamma_small / gamma_big))*gamma_big/2.0 ) ) )
  
  
         	    RR(ialpha, igap,ibias)=-CurrentPositive(ialpha,igap,ibias) / CurrentNegative(ialpha,igap,ibias)
 
                   !
	   ENDDO gap_loop  

   ENDDO alpha_loop

   ENDDO ! bias loop


!
! ... writedata on files
!

  PRINT*, "Writing output files"



  PRINT*, "  ...Writing results"    
   
   !

   DO ibias = 1,nbias 
        WRITE (outputfile_name,*) TRIM("RectificationRatio"),  biasgrid(ibias), TRIM(".dat")   
        OPEN ( Fil4, FILE=TRIM(outputfile_name), FORM='formatted' )
	   DO ialpha = 1, nalpha
		DO igap = 1, ngap
            		WRITE (Fil4, '(3(f15.9))' ) alphagrid(ialpha)*100.0 ,  gapgrid(igap) ,  RR(ialpha,igap,ibias)
        	ENDDO
		WRITE (Fil4, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil4)
   ENDDO 

   DO ibias = 1,nbias 
        WRITE (outputfile_name,*) TRIM("RectificationRatioInv"),  biasgrid(ibias), TRIM(".dat")   
        OPEN ( Fil4, FILE=TRIM(outputfile_name), FORM='formatted' )
	   DO ialpha = 1, nalpha
		DO igap = 1, ngap
            		WRITE (Fil4, '(3(f15.9))' ) alphagrid(ialpha)*100.0 ,  gapgrid(igap) ,  ONE/RR(ialpha,igap,ibias)
        	ENDDO
		WRITE (Fil4, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil4)
   ENDDO 

   IF ( debug_level > 2) THEN

          !PRINT GRIDS
	   OPEN ( Fil5, FILE='alphagrid.dat', FORM='formatted' )
	   DO ialpha = 1, nalpha
	           WRITE (Fil5, '(I10,f15.9)' )  ialpha, alphagrid(ialpha) 
	   ENDDO
	   CLOSE( Fil5 )
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='gapgrid.dat', FORM='formatted' )
	   DO igap = 1, ngap
	           WRITE (Fil5, '(I10,f15.9)' )  igap, gapgrid(igap) 
	   ENDDO
	   CLOSE( Fil5 )

          !PRINT GRIDS
	   OPEN ( Fil5, FILE='biasgrid.dat', FORM='formatted' )
	   DO ibias = 1, nbias
	           WRITE (Fil5, '(I10,f15.9)' )  ibias, biasgrid(ibias) 
	   ENDDO
	   CLOSE( Fil5 )

   ENDIF


!
!...  free memory
!

  PRINT*, "Deallocate grids"
  CALL alphagrid_deallocate()
  CALL biasgrid_deallocate()
  CALL gapgrid_deallocate()
  


  PRINT*, "Deallocate transport variables"
    DEALLOCATE ( CurrentPositive, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating current"
            STOP
         ENDIF 
    DEALLOCATE ( CurrentNegative, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating current"
            STOP
         ENDIF 
    DEALLOCATE ( PrefactorPositive, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating prefactor"
            STOP
         ENDIF 
  DEALLOCATE ( PrefactorNegative, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating prefactor"
            STOP
         ENDIF 

     DEALLOCATE ( RR, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating rectification ratio"
            STOP
         ENDIF 

 
 PRINT*, "End"


END PROGRAM conductor
