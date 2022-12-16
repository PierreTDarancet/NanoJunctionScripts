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
   PROGRAM conductor
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   USE T_control_module,                   ONLY : calculate_RR, debug_level, print_RR, print_current, bias, electronic_temperature


   USE T_egrid_module,                     ONLY : egrid_allocate, egrid_deallocate, egrid_init, ne, de, egrid, delta !ok
   USE T_gapgrid_module,                   ONLY : gapgrid_allocate, gapgrid_deallocate, gapgrid_init, ngap, gapgrid, gapmax, gapmin !ok
   USE T_gammagrid_module,                 ONLY : gammagrid_allocate, gammagrid_deallocate, gammagrid_init, ngamma, gammagrid !ok

   USE hamiltonian_workspace_module,       ONLY : hamiltonian_allocate, hamiltonian_deallocate, dimL, dimC, dimR, h00_L, h01_L, h00_R, h01_R, h00_C, h_LC, h_CR, h_CL, H_RC !ok
   USE occupation_module,                  ONLY : occupation_allocate, occupation_deallocate, occupation_function_init, f_L, f_R, occupation_function_onoff_init
   USE T_input_module,                     ONLY : input_manager !ok

   USE T_green_operation_module,           ONLY : lead_green_function, calcul_gamma, calcul_g_C !ok
   USE green_workspace_module,             ONLY : green_allocate, green_deallocate, gL, gR, gC, gamma_L , gamma_R, sigma_L, sigma_R, dos_L, dos_R !ok

   USE hamiltonian_construction_module,    ONLY : Molecular_hamiltonian_init, hamiltonian_reference_allocate, hamiltonian_reference_deallocate, hamiltonian_reference_init !ok

   USE transport_formula_module,           ONLY : transmittance_calculation, current_calculation, lorentzian_transmittance_calculation  !ok




   IMPLICIT NONE

   !
   ! local variables
   !
   COMPLEX(dbl)     :: ene
   INTEGER, PARAMETER :: Fil1 = 11 
   INTEGER, PARAMETER :: Fil2 = 12 
   INTEGER, PARAMETER :: Fil3 = 13
   INTEGER, PARAMETER :: Fil4 = 14
   INTEGER, PARAMETER :: Fil5 = 15
   INTEGER, PARAMETER :: Fil6 = 16
   INTEGER, PARAMETER :: Fil7 = 17
   CHARACTER(100)  :: outputfile_name, chr
   INTEGER          :: ie, ierr, ncount, igap, igamma
   !

! Contains transport fundamental variables
! 

   REAL(dbl), ALLOCATABLE              :: Current(:,:,:)
   REAL(dbl), ALLOCATABLE              :: RR(:,:)
 
   REAL(dbl), ALLOCATABLE              :: transmittance(:)
   REAL(dbl), ALLOCATABLE              :: conductanceperchannel_aux(:,:)

   REAL(dbl)             :: cond_aux, current_aux, biasaux, eta
!
!------------------------------
! main body
!------------------------------
!   !
   ! read input file
   !


   CALL input_manager()

   CALL hamiltonian_reference_allocate()
   CALL hamiltonian_reference_init(debug_level)  

   !
   ! init
   !
      !
  

     ALLOCATE ( transmittance(ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating transmittance"
            STOP
         ENDIF 

     ALLOCATE ( Current(ngamma, ngap,2), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating current"
            STOP
         ENDIF 
     ALLOCATE ( RR(ngamma, ngap), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating rectification ratio"
            STOP
         ENDIF 
     ALLOCATE ( conductanceperchannel_aux(ne, dimC), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating conductanceperchannel_aux"
            STOP
         ENDIF 

       RR(:,:) = ZERO
       Current(:,:,:) = ZERO

    PRINT*, 'Beginning grid allocations'

   !
    CALL egrid_allocate()
    CALL gammagrid_allocate()
    CALL gapgrid_allocate()
  
    PRINT*, 'Beginning Hamiltonian allocations'
    CALL occupation_allocate()
    CALL hamiltonian_allocate()
    CALL green_allocate(ne, dimC, dimL, dimR)

    PRINT*, 'Beginning Main Loop'

    CALL gammagrid_init()
    CALL gapgrid_init()
    CALL egrid_init()
   !
   gamma_loop: &
   DO igamma=1, ngamma

           PRINT*, "#gamma =", igamma, "of", ngamma,  "RATIO =", gammagrid(igamma)
	   gap_loop: &
	   DO igap=1, ngap
                   biasaux =ZERO


                   CALL Molecular_hamiltonian_init(h00_C, h_CR, h_LC, h_RC, h_CL, h00_L, h01_L, h00_R, h01_R, biasaux, gammagrid(igamma), gapgrid(igap), dimC, debug_level)
                   !Initialization
                   conductanceperchannel_aux(:,:)=ZERO
                   f_L(:)=ZERO
                   f_R(:)=ZERO
                   gamma_L(:,:,:)=CZERO
                   gamma_R(:,:,:)=CZERO
                   transmittance(:) = ZERO
   
        	   energy_loop: &
		   DO ie = 1, ne
		      !
		      ! grids and misc
		      !
		      ene =  egrid(ie) ! + delta * CI
                      cond_aux=ZERO

		      ! 
		      ! construct leads self-energies 
		      ! 
			 !
			 CALL lead_green_function(  gL(:,:), h00_L, h01_L, egrid(ie), dimL)
			 CALL lead_green_function(  gR(:,:), h00_R, h01_R, egrid(ie), dimR)


			 CALL calcul_gamma( gamma_L(ie,:,:), sigma_L(:,:) ,  gL(:,:), h_LC(:,:), h_CL(:,:), dimL, dimC )
			 CALL calcul_gamma( gamma_R(ie,:,:), sigma_R(:,:) ,  gR(:,:), h_RC(:,:), h_CR(:,:), dimR, dimC )
    
                         ! gC is assumed to the retarded Green Function
			 !CALL calcul_g_C( gC(ie,:,:) , h00_C(:,:), sigma_L(:,:) , sigma_R(:,:), ene, dimC)

                         !CALL transmittance_calculation(cond_aux, conductanceperchannel_aux(ie,:), gC(ie,:,:), gamma_R(ie,:,:), gamma_L(ie,:,:), dimC)
                         CALL lorentzian_transmittance_calculation(conductanceperchannel_aux(ie,:), h00_C(:,:),  gamma_R(ie,:,:), gamma_L(ie,:,:), ene, dimC)
                         transmittance(ie)= cond_aux
			 !

                         !  
 
       
		   ENDDO energy_loop

                   IF ( print_current  .OR. (debug_level > 0)) THEN
		           biasaux = bias
				   !Forward bias
				   Current(igamma,igap,1)=ZERO
				   !Transmission for resonnance 1 shifted by eta
                                   ! Left Electrode at the real potential E_Fermi (+0eV)  + bias/2
                                   ! Right Electrode at the real potential E_Fermi (+0eV) - bias/2

				   eta = -biasaux * (ONE - ( ONE/gammagrid(igamma) )  ) / (2.0*ONE)

                                   ! Left Electrode at the effective potential E_Fermi (+0eV)  + bias/2 - shift
                                   ! Right Electrode at the effective potential E_Fermi (+0eV) - bias/2 - shift

				   DO ie = 1, ne
				         !CALL occupation_function_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta, electronic_temperature) 
				         CALL occupation_function_onoff_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta,  electronic_temperature)  
				         CALL current_calculation(current_aux, conductanceperchannel_aux(ie,1),f_L(ie), f_R(ie))
				         Current(igamma,igap,1)= Current(igamma,igap,1) + (current_aux *de )
				   ENDDO

				    IF ( debug_level > 4 ) THEN
					   PRINT*, "      ...Writing Transmission times occ"    
					   WRITE(outputfile_name, '(a24,i3,i3,a4)') "Transmission_withocc_d_a", igamma, igap, ".dat"
					   OPEN ( Fil1, FILE=TRIM(outputfile_name), FORM='formatted' )
						   DO ie = 1, ne       
				 			WRITE (Fil1, '(2(f15.9))' )  egrid(ie), (conductanceperchannel_aux(ie,1)*(f_L(ie) - f_R(ie)))
				       		ENDDO
					   CLOSE(Fil1)
	
				   END IF


				   !Transmission for resonnance 2 shifted by  - eta
                                   ! Left Electrode at the effective potential E_Fermi (+0eV)  + bias/2 + shift
                                   ! Right Electrode at the effective potential E_Fermi (+0eV) - bias/2 + shift

				   eta = + biasaux * (ONE - ( ONE/gammagrid(igamma) )  ) / (2.0*ONE)
				   DO ie = 1, ne
				         !CALL occupation_function_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta,  electronic_temperature) 
				         CALL occupation_function_onoff_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta,  electronic_temperature) 
				         CALL current_calculation(current_aux, conductanceperchannel_aux(ie,2),f_L(ie), f_R(ie))
				         Current(igamma,igap,1)= Current(igamma,igap,1) + (current_aux *de )
				   ENDDO

			           IF ( debug_level > 4 ) THEN
					   PRINT*, "      ...Writing Transmission times occ"    
					   WRITE(outputfile_name, '(a24,i3,i3,a4)') "Transmission_withocc_d_b", igamma, igap, ".dat"
					   OPEN ( Fil1, FILE=TRIM(outputfile_name), FORM='formatted' )
						   DO ie = 1, ne       
				 			WRITE (Fil1, '(2(f15.9))' )  egrid(ie), (conductanceperchannel_aux(ie,2)*(f_L(ie) - f_R(ie)))
				       		ENDDO
					   CLOSE(Fil1)
	
				   END IF


		           !Reverse bias
		           biasaux = -bias
				   Current(igamma,igap,2)=ZERO
				   !Transmission for resonnance 1 shifted by eta
				   eta = - biasaux * (ONE - ( ONE/gammagrid(igamma) )  ) / (2.0*ONE)
				   DO ie = 1, ne
				         CALL occupation_function_onoff_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta,  electronic_temperature)  
				         !CALL occupation_function_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta, electronic_temperature) 
				         CALL current_calculation(current_aux, conductanceperchannel_aux(ie,1),f_L(ie), f_R(ie))
				         Current(igamma,igap,2)= Current(igamma,igap,2) + (current_aux *de )
				   ENDDO

			          IF ( debug_level > 4 ) THEN
					   PRINT*, "      ...Writing Transmission times occ"    
					   WRITE(outputfile_name, '(a24,i3,i3,a4)') "Transmission_withocc_i_1", igamma, igap, ".dat"
					   OPEN ( Fil1, FILE=TRIM(outputfile_name), FORM='formatted' )
						   DO ie = 1, ne       
				 			WRITE (Fil1, '(2(f15.9))' )  egrid(ie), (conductanceperchannel_aux(ie,1)*(f_L(ie) - f_R(ie)) )
				       		ENDDO
					   CLOSE(Fil1)
	
				   END IF

				   !Transmission for resonnance 2 shifted by  - eta
				   eta =  + biasaux * (ONE - ( ONE/gammagrid(igamma) )  ) / (2.0*ONE)
				   DO ie = 1, ne
				         !CALL occupation_function_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta, electronic_temperature) 
				         CALL occupation_function_onoff_init(f_L(ie) , f_R(ie), egrid(ie), biasaux, eta,  electronic_temperature)  
				         CALL current_calculation(current_aux, conductanceperchannel_aux(ie,2),f_L(ie), f_R(ie))
				         Current(igamma,igap,2)= Current(igamma,igap,2) + (current_aux *de )
				   ENDDO

	         	          IF ( debug_level > 4 ) THEN
					   PRINT*, "      ...Writing Transmission times occ", de, ne    
					   WRITE(outputfile_name, '(a24,i3,i3,a4)') "Transmission_withocc_i_2", igamma, igap, ".dat"
					   OPEN ( Fil1, FILE=TRIM(outputfile_name), FORM='formatted' )
						   DO ie = 1, ne       
				 			WRITE (Fil1, '(2(f15.9))' )  egrid(ie), (conductanceperchannel_aux(ie,2)*(f_L(ie) - f_R(ie)) )
				       		ENDDO
					   CLOSE(Fil1)
	
				   END IF

                   ENDIF 
  
  

                   IF ( calculate_RR ) THEN
                           !PRINT*, "Calculating RR"
 	        	    RR(igamma, igap)=ABS( Current(igamma,igap,1) / Current(igamma,igap,2) )
                   END IF

                   IF ( debug_level > 4 ) THEN
			   PRINT*, "      ...Writing Transmission"    
			   WRITE(outputfile_name, '(a12,i3,i3,a4)') "Transmission", igamma, igap, ".dat"
			   OPEN ( Fil1, FILE=TRIM(outputfile_name), FORM='formatted' )
				   DO ie = 1, ne       
		 			WRITE (Fil1, '(3(f15.9))' )  egrid(ie), conductanceperchannel_aux(ie,1), conductanceperchannel_aux(ie,2)
		       		ENDDO
			   CLOSE(Fil1)
			   PRINT*, "      ...Writing Gamma"    
			   WRITE(outputfile_name, '(a6,i3,i3,a4)') "GAMMA", igamma, igap, ".dat"
			   OPEN ( Fil1, FILE=TRIM(outputfile_name), FORM='formatted' )
				   DO ie = 1, ne       
		 			WRITE (Fil1, '(5(f15.9))' )  egrid(ie),  REAL(gamma_L(ie,1,1)),  REAL(gamma_L(ie,2,2)),  REAL(gamma_R(ie,1,1)),  REAL(gamma_R(ie,2,2))
		       		ENDDO
			   CLOSE(Fil1)

		   END IF

                   !
	   ENDDO gap_loop  

   ENDDO gamma_loop
   PRINT*, "End of the GAMMA loop"

!
! ... writedata on files
!

  PRINT*, "Writing output files"



  PRINT*, "  ...Writing results"    
   
   !
   IF ( ( print_current ) .OR. ( debug_level > 0 )) THEN
           PRINT*, "      ...Writing Current"    
	   OPEN ( Fil1, FILE='Current.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	
			WRITE (Fil1, '(4(f15.9))' ) gammagrid(igamma) , gapgrid(igap),   Current(igamma,igap,1) *77.480917, Current(igamma,igap,2) *77.480917
		ENDDO
	   ENDDO
           CLOSE(Fil1)
         PRINT*, "      ...Writing Current +"    
	   OPEN ( Fil1, FILE='Current1.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	
			WRITE (Fil1, '(3(f15.9))' ) gammagrid(igamma) , gapgrid(igap),   Current(igamma,igap,1) *77.480917*1000
		ENDDO
		WRITE (Fil1, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil1)
          PRINT*, "      ...Writing Current -"    
	   OPEN ( Fil1, FILE='Current2.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	
			WRITE (Fil1, '(3(f15.9))' ) gammagrid(igamma) , gapgrid(igap),   ABS(Current(igamma,igap,2) *77.480917*1000)
		ENDDO
		WRITE (Fil1, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil1)
 

        PRINT*, "      ...Writing Current + Small"    
	   OPEN ( Fil1, FILE='Current1Small.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	
			!IF ( Current(igamma,igap,1) *77.480917*1000 < 1000.00) THEN
                           WRITE (Fil1, '(3(f15.9))' ) gammagrid(igamma) , gapgrid(igap),    log10 (Current(igamma,igap,1) *77.480917*1000 )
                        !ELSE
                        !    WRITE (Fil1, '(3(f15.9))' ) gammagrid(igamma) , gapgrid(igap),   ONE*1000.00
                        !ENDIF
		ENDDO
		WRITE (Fil1, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil1)
          PRINT*, "      ...Writing Current - Small"    
	   OPEN ( Fil1, FILE='Current2Small.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	

			!IF (Current(igamma,igap,2) *77.480917*1000 <  1000.00) THEN
			   WRITE (Fil1, '(3(f15.9))' ) gammagrid(igamma) , gapgrid(igap),    log10( ABS(Current(igamma,igap,2) *77.480917*1000) )

                        !ELSE
                        !    WRITE (Fil1, '(3(f15.9))' ) gammagrid(igamma) , gapgrid(igap),   ONE*1000.00
                        !ENDIF

		ENDDO
		WRITE (Fil1, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil1)

   ENDIF 

  
   IF ( ( calculate_RR .AND. print_RR ) .OR. ( debug_level > 0 ) )  THEN 
	   OPEN ( Fil4, FILE='RectificationRatio.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	                IF ( RR(igamma,igap) < 10000.0) THEN

               			WRITE (Fil4, '(3(f15.9))' ) gammagrid(igamma) ,  (REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ),  RR(igamma,igap)
                        ELSE
               			WRITE (Fil4, '(3(f15.9))' ) gammagrid(igamma) ,(REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ),  9999.99
                        ENDIF
		ENDDO
		WRITE (Fil4, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil4)
   ENDIF

  IF ( ( calculate_RR .AND. print_RR ) .OR. ( debug_level > 0 ) )  THEN 
	   OPEN ( Fil4, FILE='RectificationRatioSmall.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	                IF ( RR(igamma,igap) < 20.0) THEN

               			WRITE (Fil4, '(3(f15.9))' ) gammagrid(igamma) ,(REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ),  RR(igamma,igap)
                        ELSE
               			WRITE (Fil4, '(3(f15.9))' ) gammagrid(igamma) ,(REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ), 20.000
                        ENDIF
		ENDDO
		WRITE (Fil4, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil4)
   ENDIF
   IF ( ( calculate_RR .AND. print_RR ) .OR. ( debug_level > 0 ) )  THEN 
	   OPEN ( Fil4, FILE='RectificationRatioInv.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	                IF ( RR(igamma,igap) < 10000.0) THEN
               			WRITE (Fil4, '(3(f15.9))' ) (1.00/gammagrid(igamma)) ,(REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ),  RR(igamma,igap)
                       ELSE
               			WRITE (Fil4, '(3(f15.9))' ) (1.00/gammagrid(igamma)) ,(REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ),  9999.99
                        ENDIF
		ENDDO
		WRITE (Fil4, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil4)
   ENDIF
   IF ( ( calculate_RR .AND. print_RR ) .OR. ( debug_level > 0 ) )  THEN 
	   OPEN ( Fil4, FILE='RectificationRatioInvLog.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
		DO igap = 1, ngap
	                IF ( RR(igamma,igap) < 10000.0) THEN
               			WRITE (Fil4, '(3(f15.9))' ) (1.00/gammagrid(igamma)) ,(REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ),  log10(RR(igamma,igap))
                       ELSE
               			WRITE (Fil4, '(3(f15.9))' ) (1.00/gammagrid(igamma)) ,(REAL(h00_C(1,1)) + gapgrid(igap) - gapgrid(ngap) ),  log10(RR(igamma,igap))
                        ENDIF
		ENDDO
		WRITE (Fil4, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil4)
   ENDIF

   IF ( debug_level > 2) THEN

          !PRINT GRIDS
	   OPEN ( Fil5, FILE='gammagrid.dat', FORM='formatted' )
	   DO igamma = 1, ngamma
	           WRITE (Fil5, '(I10,f15.9)' )  igamma, gammagrid(igamma) 
	   ENDDO
	   CLOSE( Fil5 )
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='gapgrid.dat', FORM='formatted' )
	   DO igap = 1, ngap
	           WRITE (Fil5, '(I10,f15.9)' )  igap, gapgrid(igap) 
	   ENDDO
	   CLOSE( Fil5 )

   ENDIF


!
!...  free memory
!

  PRINT*, "Deallocate Hamiltonian"
  CALL hamiltonian_deallocate()
  CALL green_deallocate()
  CALL occupation_deallocate()
  CALL hamiltonian_reference_deallocate()
  PRINT*, "Deallocate grids"
  CALL gammagrid_deallocate()
  CALL egrid_deallocate()
  CALL gapgrid_deallocate()
  

  PRINT*, "Deallocate transport variables"
    DEALLOCATE ( Current, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating current"
            STOP
         ENDIF 
     DEALLOCATE ( RR, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating rectification ratio"
            STOP
         ENDIF 

    DEALLOCATE ( transmittance, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating transmittance"
            STOP
         ENDIF 
     DEALLOCATE ( conductanceperchannel_aux, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating conductanceperchannel_aux"
            STOP
         ENDIF 

 PRINT*, "End"


END PROGRAM conductor
