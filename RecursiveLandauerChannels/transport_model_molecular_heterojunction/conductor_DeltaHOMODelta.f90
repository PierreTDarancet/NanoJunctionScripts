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
   USE T_control_module,                   ONLY : nprint, nprint_PDOS, nprint_transmission, calculate_RR, debug_level, print_RR, print_Photocurrent, print_current, print_Photocurrentmax !ok
   USE photon_module,                      ONLY : photongrid_allocate, photongrid_deallocate, photongrid_init, photonfrequency, nphotonfrequency, photonmax, photonmin
   USE photon_variable_module,             ONLY : photoncouplingconstant,  FGR_broadening, photonvariable_allocate, photonvariable_deallocate, photonvariable_print, elumo, ehomo, ilumo, ihomo, omegamax
   USE T_egrid_module,                     ONLY : egrid_allocate, egrid_deallocate, egrid_init, ne, de, egrid, delta, nemax !ok
   USE T_biasgrid_module,                  ONLY : biasgrid_allocate, biasgrid_deallocate, biasgrid_init, nbias, biasgrid, biasmax, biasmin !ok
   USE T_alphagrid_module,                 ONLY : alphagrid_allocate, alphagrid_deallocate, alphagrid_init, nalpha, alphagrid !ok
   USE hamiltonian_workspace_module,       ONLY : hamiltonian_allocate, hamiltonian_deallocate, dimL, dimC, dimR, h00_L, h01_L, h00_R, h01_R, h00_C, h_LC, h_CR, h_CL, H_RC !ok
   USE occupation_module,                  ONLY : occupation_allocate, occupation_deallocate, occupation_function_init, f_L, f_R, f_shifted_L, f_shifted_R, f_shifted_max_L, f_shifted_max_R !ok
   USE T_input_module,                     ONLY : input_manager !ok
   !USE T_output_module,                    ONLY : write_output
   USE T_green_operation_module,           ONLY : lead_green_function, calcul_gamma, calcul_g_C !ok
   USE green_workspace_module,             ONLY : green_allocate, green_deallocate, gL, gR, gC, gamma_L , gamma_R, sigma_L, sigma_R, dos_L, dos_R !ok
   USE hamiltonian_construction_module,    ONLY : Molecular_hamiltonian_init, hamiltonian_reference_allocate, hamiltonian_reference_deallocate, hamiltonian_reference_init !ok
   USE transport_formula_module,           ONLY : transmittance_calculation, photocurrent_calculation, current_calculation !ok
   USE print_module,                       ONLY : print_PDOS, print_transmission
 !+ include emitted light - to be implemented in an upcoming version 


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
   INTEGER          :: ie, ierr, ncount, ibias, ialpha, iph, i_R, i_L
   !

! Contains transport fundamental variables
! 

   REAL(dbl), ALLOCATABLE              :: Current(:,:)
   REAL(dbl), ALLOCATABLE              :: RR(:,:)
   REAL(dbl), ALLOCATABLE              :: Photocurrent(:,:,:)
!   REAL(dbl), ALLOCATABLE              :: Emittedlight(:,:)
   REAL(dbl), ALLOCATABLE              :: Photocurrentmax(:,:)

   REAL(dbl), ALLOCATABLE              :: transmittance(:)
   REAL(dbl), ALLOCATABLE              :: Photocurrent_aux(:)
   REAL(dbl)             :: cond_aux, current_aux, Photocurrentmax_aux
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
     ALLOCATE ( Photocurrent(nphotonfrequency, nalpha, nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating Photocurrent"
            STOP
         ENDIF 
     ALLOCATE ( Photocurrent_aux(nphotonfrequency), STAT=ierr)
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating Photocurrent"
            STOP
         ENDIF 

     ALLOCATE ( Photocurrentmax(nalpha, nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating Photocurrentmax"
            STOP
         ENDIF 

     ALLOCATE ( transmittance(nemax), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating transmittance"
            STOP
         ENDIF 

     ALLOCATE ( Current(nalpha, nbias), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating current"
            STOP
         ENDIF 
     ALLOCATE ( RR(nalpha, (INT(nbias/2))), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating rectification ratio"
            STOP
         ENDIF 

       Photocurrent(:,:,:)= ZERO
       Photocurrentmax(:,:) = ZERO
       RR(:,:) = ZERO
       Current(:,:) = ZERO
!     ALLOCATE ( Emittedlight(nalpha, nbias), STAT=ierr )
!         IF( ierr /=0 ) THEN 
!            PRINT*, "Error allocating emitted light"
!            STOP
!         ENDIF 
    PRINT*, 'Beginning grid allocations'

   !
    CALL egrid_allocate()
    CALL alphagrid_allocate()
    CALL biasgrid_allocate()
    CALL photongrid_allocate()
    CALL photonvariable_allocate(nalpha,nbias,dimC)

   !conduct(:,:) = ZERO
   !dos(:,:) = ZERO

    PRINT*, 'Beginning Hamiltonian allocations'
    CALL occupation_allocate()
    CALL hamiltonian_allocate()
    CALL green_allocate(nemax, dimC, dimL, dimR)

    PRINT*, 'Beginning Main Loop'
    CALL photongrid_init()
    CALL alphagrid_init()
    CALL biasgrid_init()
   !
   alpha_loop: &
   DO ialpha=1, nalpha
           !ncountalpha=ialpha
           PRINT*, "#Alpha =", ialpha, "of", nalpha,  "Alpha =", alphagrid(ialpha)
	   bias_loop: &
	   DO ibias=1, nbias
                   ncount = ibias
                   CALL Molecular_hamiltonian_init(h00_C, h_CR, h_LC, h_RC, h_CL, h00_L, h01_L, h00_R, h01_R, omegamax(ialpha,ibias), ihomo(ialpha,ibias,:), ilumo(ialpha,ibias,:), ehomo(ialpha,ibias,:), elumo(ialpha,ibias,:), biasgrid(ibias), alphagrid(ialpha), dimC, dimL, dimR, debug_level) 
                   !PRINT*,  ihomo(ialpha,ibias,:), ilumo(ialpha,ibias,:)

  		   CALL egrid_init(biasmax, biasmin, biasgrid(ibias)) ! To be optimized

                   !Initialization
                   Current(ialpha,ibias)=ZERO
                   f_L(:)=ZERO
                   f_R(:)=ZERO
                   f_shifted_R(:,:)=ZERO
                   f_shifted_L(:,:)=ZERO
                   f_shifted_max_R(:)=ZERO
                   f_shifted_max_L(:)=ZERO
                   gamma_L(:,:,:)=CZERO
                   gamma_R(:,:,:)=CZERO
                   transmittance(:) = ZERO
                   dos_L(:) = ZERO
                   dos_R(:) = ZERO
                   Photocurrent_aux(:) = ZERO

          	   IF ( MOD( ncount, nprint) == 0 .OR. ncount == 1 .OR. ( debug_level > 2 )) THEN
                        PRINT*, "  #Bias =  ", ibias, " of ", nbias
                        PRINT*, "         Bias = ", biasgrid(ibias)
                        PRINT*, "         Number of Energy points in the bias window = ", ne
		   ENDIF

  		   energy_loop: &
		   DO ie = 1, ne
		      !
		      ! grids and misc
		      !
		      ene =  egrid(ie)  + delta * CI
                      cond_aux=ZERO
                      current_aux=ZERO
		      ! 
		      ! construct leads self-energies 
		      ! 
			 !
			 CALL lead_green_function(  gL(:,:), h00_L, h01_L, egrid(ie), dimL)
			 CALL lead_green_function(  gR(:,:), h00_R, h01_R, egrid(ie), dimR)


			 CALL calcul_gamma( gamma_L(ie,:,:), sigma_L(:,:) ,  gL(:,:), h_LC(:,:), h_CL(:,:), dimL, dimC )
			 CALL calcul_gamma( gamma_R(ie,:,:), sigma_R(:,:) ,  gR(:,:), h_RC(:,:), h_CR(:,:), dimR, dimC )
                         DO i_L = 1, dimC
                           dos_L(ie) = dos_L(ie) + ONE* gamma_L(ie,i_L,i_L)  / PI
                         ENDDO
                         DO i_R = 1, dimC
                           dos_R(ie) = dos_R(ie) + ONE* gamma_R(ie,i_R,i_R)  / PI
                         ENDDO

                         ! gC is assumed to the retarded Green Function
			 CALL calcul_g_C( gC(ie,:,:) , h00_C(:,:), sigma_L(:,:) , sigma_R(:,:), ene, dimC)

                         CALL transmittance_calculation(cond_aux, gC(ie,:,:), gamma_R(ie,:,:), gamma_L(ie,:,:), dimC)
                         transmittance(ie)= cond_aux
			 !
                         CALL occupation_function_init(f_L(ie) , f_R(ie), f_shifted_L(:,ie) , f_shifted_R(:,ie), f_shifted_max_L(ie) , f_shifted_max_R(ie), egrid(ie), biasgrid(ibias), photonfrequency(:), omegamax(ialpha, ibias), nphotonfrequency) 
                         !
                         IF ( print_current  .OR. (debug_level > 0)) THEN
                              CALL current_calculation(current_aux, cond_aux,f_L(ie), f_R(ie))
                              Current(ialpha,ibias)=Current(ialpha,ibias) + (current_aux/REAL(ne,dbl)) *de
                         ENDIF 
                         !
                         IF ( print_Photocurrent .OR. print_Photocurrentmax .OR. (debug_level > 0)) THEN
                             !
                             !CALL photocurrent_calculation_alt(Photocurrentmax_aux, egrid(ie), gamma_L(ie,:,:), gamma_R(ie,:,:), f_L(ie), f_R(ie), f_shifted_max_L(ie), f_shifted_max_R(ie), omegamax(ialpha, ibias), photoncouplingconstant, dimC, ihomo(ialpha, ibias,1), ilumo(ialpha, ibias,1), ehomo(ialpha, ibias,1), elumo(ialpha, ibias,1), 0.00)
                             !
                             !Photocurrentmax(ialpha,ibias) = Photocurrentmax(ialpha,ibias) + (Photocurrentmax_aux/REAL(ne,dbl)) *de
                             ! 
                             !DO iph=1, nphotonfrequency
		                     !
		             !        Photocurrent(iph,ialpha,ibias) =  Photocurrent(iph,ialpha,ibias) + (Photocurrent_aux(iph)/REAL(ne,dbl)) *de
		                     !
                             !ENDDO

                             ! Previous one

                             CALL photocurrent_calculation( Photocurrent_aux(:), Photocurrentmax_aux, egrid(ie), gamma_L(ie,:,:), gamma_R(ie,:,:), f_L(ie), f_R(ie),  f_shifted_L(:,ie), f_shifted_R(:,ie), f_shifted_max_L(ie), f_shifted_max_R(ie), photonfrequency(:), photoncouplingconstant, dimC, nphotonfrequency, omegamax(ialpha, ibias), ihomo(ialpha, ibias,:), ilumo(ialpha, ibias, :), ehomo(ialpha, ibias,:), elumo(ialpha, ibias,:), FGR_broadening)
                             DO iph=1, nphotonfrequency
		                     !
		                     Photocurrent(iph,ialpha,ibias) =  Photocurrent(iph,ialpha,ibias) + (Photocurrent_aux(iph)/REAL(ne,dbl)) *de
		                     !
    
                             ENDDO
                            Photocurrentmax(ialpha,ibias) = Photocurrentmax(ialpha,ibias) + (Photocurrentmax_aux/REAL(ne,dbl)) *de
                         END IF

		   ENDDO energy_loop
                     
          	   !IF ( ( ( MOD( ncount, nprint_PDOS) == 0 ) .AND. (MOD( ncountalpha, nprintalpha_PDOS) == 0) ).OR. ncount == 1 .OR. ( debug_level > 4 )) THEN
          	   IF ( (MOD( ncount, nprint_PDOS) == 0 ) .OR. ncount == 1 .OR. ( debug_level > 4 )) THEN
                        PRINT*, " PRINTING PDOS"
                        WRITE (outputfile_name, *) "PDOS_", ialpha, ibias, ".dat"
                        CALL print_PDOS(Fil7, h00_C(:,:),  gamma_L(1:ne,:,:),  gamma_R(1:ne,:,:), f_L(1:ne), f_R(1:ne), dos_L(1:ne), dos_R(1:ne),  dimC, ne, egrid(1:ne), outputfile_name)
		   ENDIF
         	   !IF ( (( MOD( ncount, nprint_transmission) == 0 ) .AND. (MOD( ncountalpha, nprintalpha_PDOS) == 0) ) .OR. ncount == 1 .OR. ( debug_level > 4 )) THEN
         	   IF ( (MOD( ncount, nprint_transmission) == 0 ).OR. ncount == 1 .OR. ( debug_level > 4 )) THEN
                        PRINT*, " PRINTING transmittance"
                        WRITE (outputfile_name, *) "Transmittance_", ialpha, ibias, ".dat"
                        CALL print_transmission(Fil6, transmittance(1:ne), egrid(1:ne),ne, outputfile_name)
		   ENDIF
                   !
	   ENDDO bias_loop  
           PRINT*, "Calculating RR"
           IF ( calculate_RR ) THEN
		   DO ibias=1, INT(nbias/2)  !!!!!!!1bias is assumed to be symmetric !!!!!!
		       RR(ialpha, ibias)=ABS( Current(ialpha,(nbias - ibias + 1)) / Current(ialpha,ibias) )
		   ENDDO 
           END IF
   ENDDO alpha_loop
   PRINT*, "End of the Alpha loop"

!
! ... writedata on files
!

  PRINT*, "Writing output files"

!  CALL write_output(Current(:,:), Photocurrent(:,:,:), Photocurrentmax(:,:), RR(:,:),  nalpha, alphagrid(:),  nbias, biasgrid(:), nphotonfrequency, photonfrequency(:))

  PRINT*, "  ...Writing results"    
   
   !
   IF ( ( print_current ) .OR. ( debug_level > 0 )) THEN
           PRINT*, "      ...Writing Current"    
	   OPEN ( Fil1, FILE='Current.dat', FORM='formatted' )
	   DO ialpha = 1, nalpha
		DO ibias = 1, nbias
                        ! Atomic units
			!Current(ialpha,ibias) = Current(ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
			WRITE (Fil1, '(3(f15.9))' ) alphagrid(ialpha) , biasgrid(ibias),   Current(ialpha,ibias)
		ENDDO
		WRITE (Fil1, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil1)
   ENDIF 

   IF ( ( ( photoncouplingconstant  /= ZERO ) .AND. print_Photocurrentmax )  .OR. ( debug_level > 0 ) ) THEN
           PRINT*, "      ...Writing PhotoCurrentmax"    
	   OPEN ( Fil2, FILE='PhotoCurrentmax.dat', FORM='formatted' )
	   DO ialpha = 1, nalpha
		DO ibias = 1, nbias
                        ! Comparison with max value of the table
                        Photocurrentmax_aux = MAXVAL(Photocurrent(:,ialpha,ibias) ) * 6623617.82 / (27.2113845) 
			Photocurrentmax(ialpha,ibias) = Photocurrentmax(ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
			WRITE (Fil2, '(5(f15.9))' ) alphagrid(ialpha) , biasgrid(ibias),   Photocurrentmax(ialpha,ibias), omegamax(ialpha,ibias),  Photocurrentmax_aux
		ENDDO
		WRITE (Fil2, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil2)
   ENDIF     

   IF ( ( ( photoncouplingconstant  /= ZERO ) .AND. print_Photocurrentmax .AND. print_current)  .OR. ( debug_level > 0 ) ) THEN
           PRINT*, "      ...Writing PhotoCurrentmax"    
	   OPEN ( Fil2, FILE='RatioPhotoCurrentmaxCurrent.dat', FORM='formatted' )
	   DO ialpha = 1, nalpha
		DO ibias = 1, nbias
                        ! Comparison with max value of the table
                        Photocurrentmax_aux = ABS (Photocurrentmax(ialpha,ibias) /  Current(ialpha,ibias))
			!Photocurrentmax(ialpha,ibias) = Photocurrentmax(ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
			!Current(ialpha,ibias) = Current(ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
			WRITE (Fil2, '(6(f15.9))' ) alphagrid(ialpha) , biasgrid(ibias), Photocurrentmax_aux,   Photocurrentmax(ialpha,ibias), Current(ialpha,ibias), omegamax(ialpha,ibias)
		ENDDO
		WRITE (Fil2, '(2(e15.9))' )
	   ENDDO
           CLOSE(Fil2)
   ENDIF     


   IF ( ( ( photoncouplingconstant  /= ZERO ) .AND. print_Photocurrent ) .OR. ( debug_level > 0 ) ) THEN
           PRINT*, "      ...Writing PhotoCurrent"
           DO iph=1, nphotonfrequency
                   WRITE (outputfile_name, *) "Photocurrent_", iph, ".dat"
		   PRINT*, "             ...For Photon_frequency = ", iph
	    	   OPEN ( Fil3, FILE=TRIM(outputfile_name), FORM='formatted' )
		   DO ialpha = 1, nalpha
			DO ibias = 1, nbias
				Photocurrent(iph,ialpha,ibias) = Photocurrent(iph,ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
				WRITE (Fil3, '(4(f15.9))' ) alphagrid(ialpha) , biasgrid(ibias),   Photocurrent(iph,ialpha,ibias), photonfrequency(iph) 
			ENDDO
			WRITE (Fil3, '(2(e15.9))' )
		   ENDDO
		   CLOSE(Fil3)
           ENDDO
   ENDIF     

   IF ( ( ( photoncouplingconstant  /= ZERO ) .AND. print_Photocurrent ) .OR. ( debug_level > 0 ) ) THEN
           PRINT*, "      ...Writing PhotoCurrent f(photonfrequency)"
           DO ialpha = 1, nalpha
                   WRITE (outputfile_name, *) "Photocurrentofphotofrequency_", ialpha, ".dat"
		   PRINT*, "             ...For ialpha = ", ialpha
	    	   OPEN ( Fil3, FILE=TRIM(outputfile_name), FORM='formatted' )
                    DO iph=1, nphotonfrequency
        			DO ibias = 1, nbias
				!Photocurrent(iph,ialpha,ibias) = Photocurrent(iph,ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
				WRITE (Fil3, '(4(f15.9))' )  photonfrequency(iph)  , biasgrid(ibias),   Photocurrent(iph,ialpha,ibias), alphagrid(ialpha) 
			ENDDO
			WRITE (Fil3, '(2(e15.9))' )
		   ENDDO
		   CLOSE(Fil3)
           ENDDO
   ENDIF     


   IF ( ( calculate_RR .AND. print_RR ) .OR. ( debug_level > 0 ) )  THEN 
	   OPEN ( Fil4, FILE='RectificationRatio.dat', FORM='formatted' )
	   DO ialpha = 1, nalpha
	      DO ibias = 1, INT(nbias/2)
	           WRITE (Fil4, '(3(f15.9))' )  alphagrid(ialpha) , ABS(biasgrid(ibias)),   RR(ialpha,ibias)
	      ENDDO
              WRITE (Fil4, '(2(e15.9))' )
	   ENDDO
	   CLOSE( Fil4 )
   ENDIF


   IF ( debug_level > 2) THEN
          !PRINT PHOTON VARIaBLE
           CALL  photonvariable_print(Fil5, nalpha, nbias, dimC, alphagrid(:), biasgrid(:) )
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='alphagrid.dat', FORM='formatted' )
	   DO ialpha = 1, nalpha
	           WRITE (Fil5, '(I100,f15.9)' )  ialpha, alphagrid(ialpha) 
	   ENDDO
	   CLOSE( Fil5 )
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='biasgrid.dat', FORM='formatted' )
	   DO ibias = 1, nbias
	           WRITE (Fil5, '(I100,f15.9)' )  ibias, biasgrid(ibias) 
	   ENDDO
	   CLOSE( Fil5 )
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='photongrid.dat', FORM='formatted' )
	   DO iph = 1, nphotonfrequency
	           WRITE (Fil5, '(I100,f15.9)' )  iph,  photonfrequency(iph) 
	   ENDDO
	   CLOSE( Fil5 )
	   !OPEN ( Fil5, FILE='occupation.dat', FORM='formatted' )
           !DO ie = 1, ne
	   !        WRITE (Fil5, '(3(f15.9))' )  egrid(ie) , f_L(ie), f_R(ie)
           !ENDDO
	   !CLOSE( Fil5 )

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
  CALL alphagrid_deallocate()
  CALL egrid_deallocate()
  CALL biasgrid_deallocate()
  CALL photonvariable_deallocate()
  CALL photongrid_deallocate()
  !

  PRINT*, "Deallocate transport variables"

    DEALLOCATE ( Photocurrentmax, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating Photocurrentmax"
            STOP
         ENDIF 

    DEALLOCATE ( Photocurrent, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating Photocurrent"
            STOP
         ENDIF 
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
  
    DEALLOCATE ( Photocurrent_aux, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating Photocurrent_aux"
            STOP
         ENDIF 

   

 !    DEALLOCATE ( Emittedlight, STAT=ierr )
 !        IF( ierr /=0 ) THEN 
 !           PRINT*, "Error deallocating emitted light"
 !           STOP
 !        ENDIF 
 PRINT*, "End"


END PROGRAM conductor
