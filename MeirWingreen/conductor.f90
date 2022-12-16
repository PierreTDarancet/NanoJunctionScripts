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
   USE T_control_module,                   ONLY : debug_level, nmaxiter_gf, amplitude, photonfrequency, ee_method, ep_method ! = Galperin Nitzan or SCBA  !ok  !for the moment not energy dependent method = dft+sigma and  dft+sigma+EXC
   USE T_egrid_module,                     ONLY : egrid_allocate, egrid_deallocate, egrid_init, ne, de, egrid, delta, fermi_energy, electronic_temperature,  green_egrid_allocate, green_egrid_deallocate, green_egrid_init, ne_green, green_egrid, egridplusomega, egridminusomega, egridshifted_init, grid_condtogreen, evaluate_negreen !ok

   USE kpoint_module,                      ONLY : kgrid_allocate, kgrid_deallocate,  nk, kweight, kgrid_init !ok
   USE hamiltonian_workspace_module,       ONLY : hamiltonian_allocate, hamiltonian_deallocate, dimC, h00_C !ok

   USE occupation_module,                  ONLY : occupation_allocate, occupation_deallocate, occupation_function_init, f_L, f_R !ok


   USE T_green_operation_module,           ONLY :  calcul_spectralfunction, calcul_g_lesser !ok
   USE self_energy_formula_module,         ONLY :calcul_electron_photon_sigma, calcul_sigma_ee_retarded !OK
   USE transport_formula_module,           ONLY : transmittance_mw_calculation, transmittance_calculation !OK
 
 
   USE read_module,             ONLY : read_mol_selfenergy, read_imagecharge_selfenergy, read_lead_selfenergy,  read_lead_gamma, read_electron_photon !ok
   USE green_workspace_module,             ONLY : g_C, glesser, gamma_L, gamma_R, sigma_L, sigma_R, &
                                                  sigma_L_lesser, sigma_R_lesser, sigma_C_r, sigma_C_lesser, &
                                                  sigma_ep_lesser, sigma_ep_r, sigma_ee_contact_lesser, &
                                                  sigma_ee_contact_r, sigma_ee_mol_r, sigma_ee_mol_lesser, &
                                                  sigma_ee_mol_corr, dipole_matrix, A_C, &
                                                  green_allocate, green_deallocate !OK
 

   USE T_input_module,                     ONLY : input_manager




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
!   CHARACTER(100)  :: outputfile_name, chr
 
   INTEGER          :: ie, ierr, ik, ncount, i_dim
   !

! Contains transport fundamental variables
! 
   ! Non-Eq conductance
   REAL(dbl), ALLOCATABLE   :: transmittance_tot(:), transmittance(:,:), transmittance_perchannel(:,:,:)
   ! Landauer conductance
   REAL(dbl), ALLOCATABLE   :: transmittance_FisherLee_tot(:), transmittance_FisherLee(:,:)
   ! Density of States
   REAL(dbl), ALLOCATABLE   ::  dosC(:,:), dos_tot(:)
   ! G lesser
   REAL(dbl), ALLOCATABLE   ::  generalized_dosC_tot(:) generalized_dosC_tot(:) 
   ! Out of equilibrium distribution
   REAL(dbl), ALLOCATABLE   ::  generalized_distribution_tot(:), generalized_distribution(:,:) 


   REAL(dbl)              :: cond_aux
!------------------------------
! main body
!------------------------------
!   !
   ! read input file
   !


   CALL input_manager()

    PRINT*, 'Beginning grid allocations'
    PRINT*, ' ... Energy grid allocations'
   !
    CALL egrid_allocate()
    CALL evaluate_negreen(photon_frequency)
    CALL green_egrid_allocate()
    CALL green_allocate

    PRINT*, ' ... Kpoint grid allocations'
    CALL kgrid_allocate()

   !conduct(:,:) = ZERO
   !dos(:,:) = ZERO

    PRINT*, 'Beginning Hamiltonian allocations'
    CALL occupation_allocate()
    CALL green_allocate(ne_green, dimC)
    CALL hamiltonian_allocate()
    CALL hamiltonian_init()  

    CALL egrid_init()

    CALL green_egrid_init()
    CALL occupation_init(green_egrid)
    CALL egridshifted_init(photon_frequency)
    CALL kgrid_init()



   !
   ! init
   !
      !
     ALLOCATE ( transmittance_tot(ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating transmittance_tot"
            STOP
         ENDIF 

     ALLOCATE ( transmittance(nk,ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating transmittance"
            STOP
         ENDIF 
    ALLOCATE ( transmittance_perchannel(nk,ne,dimC), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating transmittance_per channel"
            STOP
         ENDIF 
    ALLOCATE (   transmittance_FisherLee_tot(ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating transmittance_FisherLee_tot"
            STOP
         ENDIF 

    ALLOCATE (    transmittance_FisherLee(nk,ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating transmittance_FisherLee"
            STOP
         ENDIF 
         
    ALLOCATE (    dos_tot(ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating  dos_tot"
            STOP
         ENDIF 

    ALLOCATE (     dosC(nk,ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating  dosC"
            STOP
         ENDIF 

    ALLOCATE (     generalized_distribution_tot(ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating  dosC"
            STOP
         ENDIF 
    ALLOCATE (     generalized_distribution(nk,ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating  dosC"
            STOP
         ENDIF 
    ALLOCATE (     generalized_dosC_tot(ne), STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating  generalized_dosC_tot"
            STOP
         ENDIF 
    
    ALLOCATE (   generalized_dosC(nk,ne) , STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error allocating  generalized_dosC"
            STOP
         ENDIF 

                 


   ! BEGINS


   CALL read_electron_photon( dipole_matrix(:,:), dimC, ep_method)
   ! For the moment the electron electron sigma is entirely static. This will (hopefully) change in the future, as well as the definition of sigma lesser, in that case it will move in the negreen loop
   !

                 read_mol_selfenergy(sigma_ee_mol_r(:,:),dimC)
                 read_imagecharge_selfenergy(sigma_ee_contact_r(:,:),dimC)
   sigma_ee_contact_lesser(:,:) = CZERO
   sigma_ee_mol_lesser(:,:)  = CZERO


   transmittance_tot(:) = ZERO
   transmittance_FisherLee_tot(:) = ZERO
   dos_tot(:) = ZERO
   generalized_distribution_tot(:) = ZERO
   generalized_dosC_tot(:) = ZERO
    PRINT*, 'Beginning Main Loop'

   k_loop: &
   DO ik=1, nk
      PRINT*, "#ik =", ik, "of", nk,  "kweight =", kweight(ik)

                   !Initialization
                   transmittance(ik,:) = ZERO
                   transmittance_per_channel(ik,:,:) = ZERO
                   transmittance_FisherLee(ik,:) = ZERO
                   dosC(ik,:) = ZERO
                   generalized_distribution(ik,:) = ZERO
                   generalized_dosC(ik,:) = ZERO

                
                   !dos_L(ik,:) = ZERO
                   !dos_R(ik,:) = ZERO


                 !Initialization for ik
                 ! Lead
                 CALL read_lead_selfenergy(  sigma_L(:,:,:), sigma_R(:,:,:), ik, ne_green, dimC)
                 CALL read_lead_gamma( gamma_L(:,:,:), gamma_R(:,:,:), ik,  ne_green,dimC)

                 DO ie = 1, ne_green 
                         sigma_L_lesser(ie,:,:) = CI * f_L(ie) * gamma_L(ie,:,:)
                         sigma_R_lesser(ie,:,:) = CI * f_R(ie) * gamma_R(ie,:,:)
                 ENDDO

                 ! Init
  

                 ! Electron-photon

                 DO ie = 1, ne_green 
                         ! Spectral function
                         sigma_C_r(ie,:,:) = sigma_ee_contact_r(:,:) +   sigma_ee_mol_r(:,:) + sigma_ep_r(ie,:,:) 
		         calcul_spectralfunction( A_C(ie,:,:) , g_C(ie,:,:), h00_C(:,:), sigma_L(ie,:,:), sigma_R(ie,:,:) , sigma_C_r(ie,:,:) , ene, dimC)
                         !init glesser
                         glesser(ie,:,:) = CI* f_L(ie) * A_C(ie,:,:)

                         ! G_lesser
                         CALL calcul_electron_photon_sigma(sigma_ep_r(ie,:,:), sigma_ep_lesser(ie,:,:),  g_C(egridplusomega(ie),:,:), glesser(egridplusomega(ie),:,:), g_C(egridminusomega(ie),:,:), glesser(egridminusomega(ie),:,:), amplitude, photonfrequency, dipole_matrix(:,:), dimC, ep_method)          
                         sigma_C_lesser(ie,:,:) = sigma_ee_contact_lesser(:,:) +   sigma_ee_mol_lesser(:,:) + sigma_ep_lesser(ie,:,:) 
                 	 calcul_g_lesser( glesser(ie,:,:) , sigma_L_lesser(ie,:,:), sigma_R_lesser(ie,:,:) , sigma_C_lesser(ie,:,:) , g_C(ie,:,:), dimC)
                 ENDDO
 


                 IF (nmaxiter_gf>=1) THEN
 		         DO iter=1, nmaxiter_gf
                           !
	 		   green_function_energy_loop: &
			   DO ie = 1, ne_green
			      !
			      ! grids and misc
			      !
			      ene =   green_egrid(ie)  + delta * CI
	  		      ! 
			      ! construct self-energies 
			      ! 
			      !         , egridminusomega,

                              CALL calcul_electron_photon_sigma(sigma_ep_r(ie,:,:), sigma_ep_lesser(ie,:,:),  g_C(egridplusomega(ie),:,:), glesser(egridplusomega(ie),:,:), g_C(egridminusomega(ie),:,:), glesser(egridminusomega(ie),:,:), amplitude, photonfrequency, dipole_matrix(:,:), dimC, ep_method)          
                              !CALL calcul_electron_electron_sigma(sigma_ee_contact_r(:,:),  sigma_ee_mol_r(ie,:,:), sigma_ee_contact_lesser(:,:),  sigma_ee_mol_lesser(ie,:,:), ik, dimC, ee_method)        !  static for the moment

                              sigma_C_lesser(ie,:,:) = sigma_ee_contact_lesser(:,:) +   sigma_ee_mol_lesser(:,:) + sigma_ep_lesser(ie,:,:) 
                              sigma_C_r(ie,:,:) = sigma_ee_contact_r(:,:) +   sigma_ee_mol_r(:,:) + sigma_ep_r(ie,:,:) 

 
		              ! green function
		              calcul_spectralfunction( A_C(ie,:,:) , g_C(ie,:,:), h00_C(:,:), sigma_L(ie,:,:), sigma_R(ie,:,:) , sigma_C_r(ie,:,:) , ene, dimC)
		              calcul_g_lesser( glesser(ie,:,:) , sigma_L_lesser(ie,:,:), sigma_R_lesser(ie,:,:) , sigma_C_lesser(ie,:,:) , g_C(ie,:,:), dimC)

		           ENDDO
                              ! Correction to the band structure due to Excitonic effects
                              CALL calcul_sigma_ee_retarded( glesser(:,:,:),  A_C(:,:,:), sigma_ee_mol_corr(:), fermi_energy,  electronic_temperature, green_egrid(:) ,  dimC, ne_green, ee_method)
                              DO i_dim = 1, dimC
                                   sigma_C_r(ie,i_dim,i_dim) =  sigma_C_r(ie,i_dim,i_dim)  + sigma_ee_mol_corr(i_dim) 
                              ENDDO

		          ENDDO
                  ENDIF


                 DO ie=1, ne
                      !              
                      cond_aux=ZERO
                      CALL  transmittance_mw_calculation( cond_aux, transmittance_perchannel(ik,ie,:) ,  f_L(grid_condtogreen(ie)), f_R(grid_condtogreen(ie)), A_C(grid_condtogreen(ie),:,:) , g_lesser(grid_condtogreen(ie),:.:), gamma_R(grid_condtogreen(ie),:,:), gamma_L(grid_condtogreen(ie),:,:), dimC)
 		      !   SUBROUTINE  transmittance_mw_calculation(conductance, conduct(dimx),  f_L, f_R, A_C(dimx,dimx), gC(dimx,dimx), gamR(dimx,dimx), gamL(dimx,dimx), dimx)
                      transmittance(ik,ie) = cond_aux


                      generalized_dosC(ik,ie) = ZERO
                      DO i_dim = 1, dimC
                          generalized_dosC(ik,ie) = generalized_dosC(ik,ie) - (ONE/PI) *  AIMAG(g_lesser(grid_condtogreen(ie),i_dim,i_dim) )
                      ENDDO 


                      dosC(ik,ie) = ZERO
                      DO i_dim = 1, dimC
                          dosC(ik,ie) = dosC(ik,ie) + (ONE/PI) *  A_C(grid_condtogreen(ie),i_dim,i_dim)
                      ENDDO 

                      ! Attempt not rigorous
                      generalized_distribution(ik,ie) =    ( generalized_dosC(ik,ie) )  / ( dosC(ik,ie) + EPS_m2 )





                      cond_aux=ZERO
                      CALL  transmittance_calculation( cond_aux,  g_C(grid_condtogreen(ie),:,:),   gamma_R(grid_condtogreen(ie),:,:), gamma_L(grid_condtogreen(ie),:,:), dimC)
                      ! transmittance_calculation(conductance, gC(dimx,dimx), gamR(dimx,dimx), gamL(dimx,dimx), dimx)
                      transmittance_FisherLee(ik,ie) = cond_aux
                      !


          	 ENDDO energy_loop

                 transmittance_tot(:) = transmittance_tot(:) + transmittance(ik,:)* kweight(ik)
                 transmittance_FisherLee_tot(:) = transmittance_FisherLee_tot(:) + transmittance_FisherLee(ik,:)* kweight(ik)
                 dos_tot(:) = dos_tot(:) + dosC(ik,:)* kweight(ik) 
                 generalized_distribution_tot(:) = generalized_distribution_tot(:) + generalized_distribution(ik,:) * kweight(ik)
                 generalized_dosC_tot(:) = generalized_dosC_tot(:) +  generalized_dosC(ik,:)  * kweight(ik)
                 !
                 !
   ENDDO    k_loop
   PRINT*, "End of the k loop"

!
! ... writedata on files
!

  PRINT*, "Writing output files"

	   OPEN ( Fil1, FILE='Transmittance.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), transmittance_tot(ie)
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='Transmittance_FisherLee.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), transmittance_FisherLee_tot(ie)
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='dos.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), dos_tot(ie) 
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='distribution.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), generalized_distribution_tot(ie)  
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='glesser.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), generalized_dosC_tot(ie)  
	   ENDDO
	   CLOSE( Fil1 )


  
  PRINT*, "Writing optional output files (if required)"
   IF ( debug_level > 2) THEN

	   OPEN ( Fil1, FILE='Transmittance_per_k.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), (transmittance(ik,ie), ik=1, nk)
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='Transmittance_FisherLee_per_k.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), ( transmittance_FisherLee(ik,ie), ik=1, nk)
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='dos_per_k.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie),  (dosC(ik,ie), ik=1, nk)
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='distribution_per_k.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), (generalized_distribution(ik,ie), ik=1, nk)
	   ENDDO
	   CLOSE( Fil1 )

	   OPEN ( Fil1, FILE='glesser_per_k.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), (generalized_dosC(ik,ie), ik=1, nk)
	   ENDDO
	   CLOSE( Fil1 )

   ENDIF


   IF ( debug_level > 3) THEN

	   OPEN ( Fil1, FILE='Transmittance_per_k_per_channel.dat', FORM='formatted' )
           DO ik=1, nk
	   WRITE (Fil1, '(I100,f15.9)' ) ik
	   DO ie = 1, ne
	           WRITE (Fil1, '(I100,f15.9)' ) egrid(ie), (transmittance_perchannel(ik,ie,i_dim), i_dim=1, dimC)
	   ENDDO
	   ENDDO
	   CLOSE( Fil1 )
   ENDIF
                 

  PRINT*, "  ...Writing results"    
   
   !
 

   IF ( debug_level > 2) THEN
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='kparallelgrid.dat', FORM='formatted' )
	   DO ik = 1, nk
	           WRITE (Fil5, '(I100,f15.9)' )  ik, kweight(ik) 
	   ENDDO
	   CLOSE( Fil5 )
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='energygrid.dat', FORM='formatted' )
	   DO ie = 1, ne
	           WRITE (Fil5, '(I100,f15.9)' )  ie, egrid(ie),  green_egrid(grid_condtogreen(ie)) 
	   ENDDO
	   CLOSE( Fil5 )
          !PRINT GRIDS
	   OPEN ( Fil5, FILE='energygrid_greenfunction.dat', FORM='formatted' )
	   DO ie = 1, ne_green
	           WRITE (Fil5, '(I100,f15.9)' )  ie,  green_egrid(ie) 
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

  CALL egrid_deallocate()
  CALL green_egrid_deallocate()
  !

  PRINT*, "Deallocate transport variables"
 
    DEALLOCATE ( transmittance, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating transmittance"
            STOP
         ENDIF 
    DEALLOCATE ( transmittance_tot, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating transmittance_tot"
            STOP
         ENDIF 
    DEALLOCATE (  transmittance_FisherLee_tot, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
    DEALLOCATE (  transmittance_FisherLee, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
    DEALLOCATE (dos_tot , STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
    DEALLOCATE ( dosC, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
    DEALLOCATE (  generalized_distribution_tot, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
    DEALLOCATE (  generalized_distribution, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
    DEALLOCATE ( generalized_dosC_tot, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
    DEALLOCATE (generalized_dosC , STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating "
            STOP
         ENDIF 
   
 PRINT*, "End"


END PROGRAM conductor
