    PROGRAM pdos_plot
        IMPLICIT NONE


        ! 
        ! Parameters
        !
 
         ! Input Parameters
         !Number of energies in the Siesta PDOS file
         INTEGER, PARAMETER :: Nb_energy_in=2000
         !Number of orbitals in the Siesta PDOS file
         INTEGER, PARAMETER :: norb=2248
         ! Name of the input PDOS file
         CHARACTER*50, PARAMETER :: fname_in="BAPA_junc_final.PDOS"
         !

         ! Output Parameters
		 ! Name of the output PDOS file
		 CHARACTER*50, PARAMETER :: fname_out="LDOS_JuncBAPA_"
         !Energy Grid
		 !Number of energies in the Output PDOS file
		 INTEGER, PARAMETER :: Nb_energy_out=2000 
		 !Minimum energy of the energy spectrum of the output PDOS file
		 REAL, PARAMETER :: Output_Energy_grid_min=-15.0
		 !Maximum energy of the energy spectrum of the output PDOS file
		 REAL, PARAMETER :: Output_Energy_grid_max=5.0
                 ! Fermi energy (Just used as an origin in energy in the output grid, can be set to 0)
                 REAL, PARAMETER :: Fermi_energy=-4.9880
         !Position Grid
                 ! Number of points on the spatial grid
                 INTEGER, PARAMETER :: nspatialgrid=200
		 !Maximum z of the position grid  - 
                 !Please note that Spatial_grid_max can be calculated automatically 
                 !Starting from the positions of the atomic orbitals
                 !For consistency between different systems though, I am keeping that as a parameter 
                 ! Just comment this line if needed and uncomment later in the code - indicated
		 REAL, PARAMETER ::Spatial_grid_max =60.0
		 !Minimum z of the position grid 
                 !Please note that Spatial_grid_min can be calculated automatically 
                 !Starting from the positions of the atomic orbitals
                 !For consistency between different systems though, I am keeping that as a parameter 
                 ! Just comment this line if needed and uncomment later in the code - indicated
		 REAL, PARAMETER ::Spatial_grid_min =-25.0
                 !Lowest z printed in the output file
                 REAL, PARAMETER :: print_left=-2.7
                  !Highest z printed in the output file                
                 REAL, PARAMETER :: print_right=37.0


         ! Broadening/Convolution/Filtering Parameters
                 ! Geometry:
                 !         !          !          !
                 !  Left   !          !  Right   !
                 !Electrode! Molecule !Electrode !
                 !         !          !          !
                 !         !          !          !
                 !Left electrode limits (z) (for left electrode broadening Scheme) in Bohr
                 REAL, PARAMETER :: electrodeL_left=-50.0
                 REAL, PARAMETER :: electrodeL_right=2.8 
                 !Molecule limits (z) (for molecule broadening Scheme) in Bohr
                 REAL, PARAMETER :: molecule_left=2.9
                 REAL, PARAMETER :: molecule_right=32.0
                 !Right electrode limits (z) (for right electrode broadening Scheme) in Bohr
                 REAL, PARAMETER :: electrodeR_left=32.1
                 REAL, PARAMETER :: electrodeR_right=100.00
                 ! Left Electrode Gaussian smearing parameter(\sigma) 
                 !      For energy
                 REAL, PARAMETER :: energy_broadening_L=0.1
                 !      For position
                 REAL, PARAMETER :: spatial_broadening_L=5.0
                 ! Left Electrode Number of points used in the gaussian summation: 
                 ! In principle, this could be set to the number of energy points in the output grid
                 ! However it would make the calculation reaaaally slow.
                 ! In practice, I just set up this parameter to sum over 10*broadening for energy, less for position  (2* for electrode in general)
                 !       For energy: Resolution of the grid 0.01eV. Broadening*10 = 1eV = 100 energy points
                 INTEGER, PARAMETER :: energyinterpol_L=100
                 !      For position:Resolution of the grid 0.425Bohr.  Broadening*2 = 10Bohr = 23 Position points
                 INTEGER, PARAMETER :: spatialinterpol_L=23
                 !Right Electrode Gaussian smearing parameter(\sigma) 
                 !      For energy 
                 REAL, PARAMETER :: energy_broadening_R=0.1
                 !      For position
                 REAL, PARAMETER :: spatial_broadening_R=5.0
                 ! Right Electrode Number of points used in the gaussian summation: 
                 ! In principle, this could be set to the number of energy points in the output grid
                 ! However it would make the calculation reaaaally slow.
                 ! In practice, I just set up this parameter to sum over 10*broadening for energy, less for position (2* for electrode in general)
                 !       For energy: Resolution of the grid 0.01eV. Broadening*10 = 1eV = 100 energy points
                 INTEGER, PARAMETER :: energyinterpol_R=100
                 !      For position:Resolution of the grid 0.425Bohr.  Broadening*2 = 10Bohr = 23 Position points
                 INTEGER, PARAMETER :: spatialinterpol_R=23
                 ! Molecule Gaussian smearing parameter(\sigma) 
                 !      For energy
                 REAL, PARAMETER :: energy_broadening_C= 0.01 ! No Broadening
                 !      For position
                 REAL, PARAMETER :: spatial_broadening_C=1.0
                 ! Left Electrode Number of points used in the gaussian summation: 
                 ! In principle, this could be set to the number of energy points in the output grid
                 ! However it would make the calculation reaaaally slow.
                 ! In practice, I just set up this parameter to sum over 10*broadening for energy, less for position (4* for molecule in general)
                 !       For energy: No Broadening
                 INTEGER, PARAMETER :: energyinterpol_C=1
                 !      For position:Resolution of the grid 0.425Bohr.  Broadening*4 = 4Bohr = 8 Position points
                 INTEGER, PARAMETER :: spatialinterpol_C=8

        ! 
        ! Parameters End
        !
 
        ! ==================================================================================!
        ! ==================================================================================!
        ! ==================================================================================!

        ! 
        ! Variables Definitions 
        !


        ! Local Variables 

        ! INPUT related VARIABLES   
	        !READING PDOS File     
	        ! Useless things (browsing the file)
	        CHARACTER*100  :: chr
	        CHARACTER*1    :: chr11
	        ! Important things
	        ! Input energy grid
	        REAL, ALLOCATABLE :: In_energygrid(:)
	        ! energy grid mesh
 	        REAL :: Delta_energy_input_energy_grid
	        ! Input PDOS
	        REAL, ALLOCATABLE :: Input_orb_PDOS(:,:)
	        ! Input Orbital center positions
	        REAL, ALLOCATABLE :: Orb_center_coordinates(:,:)
        ! OUTPUT related VARIABLES     
                ! Spatial   
	        ! Spatial Grid on the z direction
	        REAL, ALLOCATABLE :: Spatial_grid(:) 
		 !Maximum z of the position grid  - Uncomment for automatical choice
	        ! Delta z of the Spatial Grid 
	        REAL :: dspace
                ! OUTPUT Energy grid (Nb_energy_out)
                REAL, ALLOCATABLE :: Out_energygrid(:)
                ! Delta energy of the Output energy Grid
                 REAL :: Delta_energy_output_energy_grid
	        ! OUTPUT PDOS(z,Nb_energy_out)
	        REAL, ALLOCATABLE :: PDOS_raw_smeared(:,:) 
        ! Variables used in the calculations: 
                ! PDOS arrays
	        REAL, ALLOCATABLE :: PDOS_L(:,:) 
	        REAL, ALLOCATABLE :: PDOS_C(:,:)
	        REAL, ALLOCATABLE :: PDOS_R(:,:)  
                REAL, ALLOCATABLE :: PDOS_raw_smeared_L(:,:)
                REAL, ALLOCATABLE :: PDOS_raw_smeared_C(:,:)
                REAL, ALLOCATABLE :: PDOS_raw_smeared_R(:,:)
               ! Convolution/Broadening Variables
	         REAL ::   summation_raw
	         REAL :: spatial_broadening
	         REAL :: energy_broadening
	         INTEGER :: spatialinterpol
	         INTEGER :: energyinterpol
 	       ! Test Variables
	        LOGICAL :: inthegrid
 	       ! LOOP Variables
 	       INTEGER :: ie, igrid, iat, ie2, inter, inter2

        ! 
        ! End of Definitions 
        !

        ! ==================================================================================!
        ! ==================================================================================!
        ! ==================================================================================!

        ! 
        ! Main Body
        !

        PRINT*, "Beginning..."

        ! Allocation part       
        PRINT*, "Beginning Memory Allocations"
        ! Input grids
        ALLOCATE (Spatial_grid(nspatialgrid) ) 
        ALLOCATE (In_energygrid(Nb_energy_in)  ) 
        ALLOCATE (Input_orb_PDOS(norb,Nb_energy_in))
        ! Output Grids
        ALLOCATE (PDOS_L(nspatialgrid,Nb_energy_out)  ) 
        ALLOCATE (PDOS_C(nspatialgrid,Nb_energy_out)  ) 
        ALLOCATE (PDOS_R(nspatialgrid,Nb_energy_out)  ) 
        ALLOCATE (Out_energygrid(Nb_energy_out)  ) 
        ALLOCATE (Orb_center_coordinates(3,norb)) 

        PRINT*, "End of Allocations"
 
        ! Reading Input File

        PRINT*, "READING Incoming PDOS File"
      open(100,file=trim(fname_in),form='formatted')
      rewind(100)
      read(100,*)chr
      read(100,*)chr
      read(100,*)chr
      read(100,*)chr
      do ie = 1,Nb_energy_in
        read(100,*)In_energygrid(ie)
      end do
      read(100,*)chr
      do iat = 1,norb
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,"(a11,3f11.6,a1)")chr, Orb_center_coordinates(1,iat), Orb_center_coordinates(2,iat), &
     &  Orb_center_coordinates(3,iat),chr11
      ! print*,trim(chr),trim(chr11),Orb_center_coordinates(1:3,iat)
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       do ie = 1, Nb_energy_in
           read(100,*)Input_orb_PDOS(iat,ie)
       end do
       read(100,*)chr
       read(100,*)chr
      end do
      close(100)
        PRINT*, "END of READING"

   ! End of reading


   ! Grid Constuction
        PRINT*, "Constructing Output Grids:"
        PRINT*, "...Spatial Grid Initialization"

   ! Spatial Grid
       ! Define spatial grid extrema
       ! Uncomment for automatical choice
       !     Spatial_grid_max =  MAXVAL( (Orb_center_coordinates(3,:)) )  +2.0
       !     Spatial_grid_min =  MINVAL( (Orb_center_coordinates(3,:))  ) -2.0
       ! Spatial separation of two elements of the spatial grid
     dspace = (Spatial_grid_max - Spatial_grid_min) / REAL(nspatialgrid -1)
       ! Spatial Grid calculation
     DO igrid = 1, nspatialgrid
       Spatial_grid(igrid) = Spatial_grid_min + REAL(igrid -1) * dspace
     ENDDO 

    ! Energy Grid 
        PRINT*, "...Energy Grid Initialization"
       ! Energy separation on the output energy grid
     Delta_energy_output_energy_grid = (Output_Energy_grid_max - Output_Energy_grid_min) / REAL(Nb_energy_out-1)
       ! Energy Grid calculation   
     DO ie = 1, Nb_energy_out
       Out_energygrid(ie) =Output_Energy_grid_min + REAL(ie -1) * Delta_energy_output_energy_grid
     ENDDO    
     Delta_energy_input_energy_grid = ABS ( In_energygrid(1) - In_energygrid(2) )

        PRINT*, "End of Output Grids Initialization"
    ! End of Grids Construction

    ! OUTPUT PDOS INIT
        PRINT*, "PDOS Initialization"
    PDOS_L(:,:)= 0.00000000 
    PDOS_C(:,:)= 0.00000000 
    PDOS_R(:,:)= 0.00000000 
    ! OUTPUT PDOS calculation

    ! Orbital Loop
        PRINT*, "PDOS Calculation"

    DO iat=1, norb
       ! Test variable
       inthegrid=.FALSE.
       ! Spatial Grid loop
       grid_loop: &
       DO igrid = 1, nspatialgrid-1
            ! Is Orbital center in the range igrid - igrid+1 ?
            IF ( (Orb_center_coordinates(3,iat) <= Spatial_grid(igrid+1) ).AND. ( Orb_center_coordinates(3,iat) > Spatial_grid(igrid)) )THEN

               IF  ( (Orb_center_coordinates(3,iat) <= electrodeL_right).AND. ( Orb_center_coordinates(3,iat) >electrodeL_left) )THEN
		       ! Yes
		       ! Energy loop
		        DO ie = 1, Nb_energy_in
		              ! Is energy on the output energy grid ?
		              IF ( ( In_energygrid(ie) <= Output_Energy_grid_max ).AND. ( In_energygrid(ie) > Output_Energy_grid_min) ) THEN
		                 ! Yes, Enter contribution
		                 DO ie2 = 1, Nb_energy_out-1
		                     IF ( ( In_energygrid(ie) <= (Out_energygrid(ie2)+Delta_energy_output_energy_grid) ).AND. ( In_energygrid(ie) > Out_energygrid(ie2)) )THEN
		                        PDOS_L(igrid,ie2) = PDOS_L(igrid,ie2) + ( Input_orb_PDOS(iat,ie) * Delta_energy_input_energy_grid/ Delta_energy_output_energy_grid  ) 
		                     ENDIF
		                 ENDDO
		                 !
		             ENDIF
		             ! No, nothing happens  
		        ENDDO
               ELSE IF  ( (Orb_center_coordinates(3,iat) <= electrodeR_right).AND. ( Orb_center_coordinates(3,iat) > electrodeR_left) )THEN
		       ! Yes
		       ! Energy loop
		        DO ie = 1, Nb_energy_in
		              ! Is energy on the output energy grid ?
		              IF ( ( In_energygrid(ie) <= Output_Energy_grid_max ).AND. ( In_energygrid(ie) > Output_Energy_grid_min) ) THEN
		                 ! Yes, Enter contribution
		                 DO ie2 = 1, Nb_energy_out-1
		                     IF ( ( In_energygrid(ie) <= (Out_energygrid(ie2)+Delta_energy_output_energy_grid) ).AND. ( In_energygrid(ie) > Out_energygrid(ie2)) )THEN

		                        PDOS_R(igrid,ie2) = PDOS_R(igrid,ie2) + ( Input_orb_PDOS(iat,ie) * Delta_energy_input_energy_grid/ Delta_energy_output_energy_grid  ) 

		                     ENDIF
		                 ENDDO
		                 !
		             ENDIF
		             ! No, nothing happens  
		        ENDDO

               ELSE IF  ( (Orb_center_coordinates(3,iat) <= molecule_right).AND. ( Orb_center_coordinates(3,iat) >=molecule_left) )THEN
		      ! Yes
		       ! Energy loop
		        DO ie = 1, Nb_energy_in
		              ! Is energy on the output energy grid ?
		              IF ( ( In_energygrid(ie) <= Output_Energy_grid_max ).AND. ( In_energygrid(ie) > Output_Energy_grid_min) ) THEN
		                 ! Yes, Enter contribution
		                 DO ie2 = 1, Nb_energy_out-1
		                     IF ( ( In_energygrid(ie) <= (Out_energygrid(ie2)+Delta_energy_output_energy_grid) ).AND. ( In_energygrid(ie) > Out_energygrid(ie2)) )THEN
		                        PDOS_C(igrid,ie2) = PDOS_C(igrid,ie2) + ( Input_orb_PDOS(iat,ie) * Delta_energy_input_energy_grid/ Delta_energy_output_energy_grid  ) 
		                     ENDIF
		                 ENDDO
		                 !
		             ENDIF
		             ! No, nothing happens  
		        ENDDO

		ENDIF
		       ! Set the set variable to true and exit the llop
		       inthegrid=.TRUE.
		       EXIT grid_loop
           ENDIF

           ! No, nothing happens 
      ENDDO grid_loop
      IF  (inthegrid==.FALSE.) THEN
         PRINT*, "================================================================"
         PRINT*, "Error, atomic orbital number ", iat, "is not in the spatial grid"
         PRINT*, "Position ", Orb_center_coordinates(3,iat), "doesn't belong to [", Spatial_grid_min ," ", Spatial_grid_max, "]"
         PRINT*, "================================================================"
         STOP
       ENDIF
    END DO
    ! End of the Orbital loop
        PRINT*, "End of PDOS Calculation"



        PRINT*, "Fancy Convolutions and Filtering"
   ! ==================================================================================!
   ! PUT ALL THE FANCY FILTERING, CONVOLUTION, SMEARING HERE !!!!!!!!!!!!!!!!!!!!!!!!
   ! ==================================================================================!


   ALLOCATE (PDOS_raw_smeared(nspatialgrid,Nb_energy_out)) 
   ALLOCATE (PDOS_raw_smeared_C(nspatialgrid,Nb_energy_out)) 
   ALLOCATE (PDOS_raw_smeared_L(nspatialgrid,Nb_energy_out)) 
   ALLOCATE (PDOS_raw_smeared_R(nspatialgrid,Nb_energy_out)) 

   PDOS_raw_smeared(:,:) = 0.00
   PDOS_raw_smeared_C(:,:) = 0.00
   PDOS_raw_smeared_L(:,:) = 0.00
   PDOS_raw_smeared_R(:,:) = 0.00
   DO igrid = 1+spatialinterpol, nspatialgrid-spatialinterpol

          

         IF ( ( spatial_grid(igrid) <=   electrodeL_right ) .AND. ( spatial_grid(igrid) >=  electrodeL_left) ) THEN
            spatialinterpol = spatialinterpol_L
            spatial_broadening = spatial_broadening_L
            energyinterpol = energyinterpol_L
            energy_broadening = energy_broadening_L 
            PRINT*, "In L"
		 DO ie = 1, Nb_energy_out
		     summation_raw=0.0000000
		     DO inter = 1,  spatialinterpol
		       IF (igrid-inter > 0) THEN  
		        DO inter2 = 1, energyinterpol                  
		          IF (ie-inter2 > 0) THEN  
		            summation_raw  = summation_raw + PDOS_L((igrid-inter),(ie-inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		         IF (ie+inter2 <= Nb_energy_out) THEN  
                         	summation_raw  = summation_raw + PDOS_L((igrid-inter),(ie+inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		        ENDDO
		       ENDIF
		       IF (igrid+inter <= nspatialgrid) THEN  
		       DO inter2 = 1, energyinterpol                  
		          IF (ie-inter2 > 0) THEN  
			summation_raw  = summation_raw + PDOS_L((igrid+inter),(ie-inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		         IF (ie+inter2 <= Nb_energy_out) THEN  
			summation_raw  = summation_raw + PDOS_L((igrid+inter),(ie+inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		        ENDDO

		       ENDIF
		     ENDDO
		     PDOS_raw_smeared_L(igrid,ie) = PDOS_L(igrid,ie) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) ))   + summation_raw
	  
		 ENDDO

         ELSE IF ( ( spatial_grid(igrid) <=   electrodeR_right ) .AND. ( spatial_grid(igrid) >=  electrodeR_left ) ) THEN
            spatialinterpol = spatialinterpol_R
            spatial_broadening = spatial_broadening_R
            energyinterpol = energyinterpol_R
            energy_broadening = energy_broadening_R 
            PRINT*, "In R"
		 DO ie = 1, Nb_energy_out
		     summation_raw=0.0000000
		     DO inter = 1,  spatialinterpol
		       IF (igrid-inter > 0) THEN  
		        DO inter2 = 1, energyinterpol                  
		          IF (ie-inter2 > 0) THEN  
		            summation_raw  = summation_raw + PDOS_R((igrid-inter),(ie-inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		         IF (ie+inter2 <= Nb_energy_out) THEN  
                         	summation_raw  = summation_raw + PDOS_R((igrid-inter),(ie+inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		        ENDDO
		       ENDIF
		       IF (igrid+inter <= nspatialgrid) THEN  
		       DO inter2 = 1, energyinterpol                  
		          IF (ie-inter2 > 0) THEN  
			summation_raw  = summation_raw + PDOS_R((igrid+inter),(ie-inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		         IF (ie+inter2 <= Nb_energy_out) THEN  
			summation_raw  = summation_raw + PDOS_R((igrid+inter),(ie+inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		        ENDDO

		       ENDIF
		     ENDDO
		     PDOS_raw_smeared_R(igrid,ie) = PDOS_R(igrid,ie) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) ))   + summation_raw
	  
		 ENDDO


         ELSE IF ( ( spatial_grid(igrid) <=  molecule_right  ) .AND. ( spatial_grid(igrid) >= molecule_left) ) THEN
            spatialinterpol = spatialinterpol_C
            spatial_broadening = spatial_broadening_C
            energyinterpol = energyinterpol_C
            energy_broadening = energy_broadening_C 
            PRINT*, "In C"
		 DO ie = 1, Nb_energy_out
		     summation_raw=0.0000000
		     DO inter = 1,  spatialinterpol
		       IF (igrid-inter > 0) THEN  
		        DO inter2 = 1, energyinterpol                  
		          IF (ie-inter2 > 0) THEN  
		            summation_raw  = summation_raw + PDOS_C((igrid-inter),(ie-inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		         IF (ie+inter2 <= Nb_energy_out) THEN  
                         	summation_raw  = summation_raw + PDOS_C((igrid-inter),(ie+inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		        ENDDO
		       ENDIF
		       IF (igrid+inter <= nspatialgrid) THEN  
		       DO inter2 = 1, energyinterpol                  
		          IF (ie-inter2 > 0) THEN  
			summation_raw  = summation_raw + PDOS_C((igrid+inter),(ie-inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		         IF (ie+inter2 <= Nb_energy_out) THEN  
			summation_raw  = summation_raw + PDOS_C((igrid+inter),(ie+inter2)) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * exp(-  ((inter*dspace)**2) /( 2* spatial_broadening**2)) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
		         ENDIF
		        ENDDO

		       ENDIF
		     ENDDO
		     PDOS_raw_smeared_C(igrid,ie) = PDOS_C(igrid,ie) * ( 1.00/ sqrt(2 * 3.1415 * ( spatial_broadening**2) )) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) ))   + summation_raw
	  
		 ENDDO





         ELSE 
            spatialinterpol = 1
            spatial_broadening = 1000.000
            energyinterpol = 1 
            energy_broadening = 1000.000
         ENDIF 

 
   ENDDO



        PRINT*, "Merging"
        DO igrid = 1, nspatialgrid
                  IF ( (Spatial_grid(igrid) <= molecule_right ).AND. ( Spatial_grid(igrid) >= molecule_left) ) THEN
                     PDOS_raw_smeared(igrid,:) = PDOS_raw_smeared_C(igrid,:)
                  ELSE IF ( (Spatial_grid(igrid) <= electrodeL_right ).AND. ( Spatial_grid(igrid) >= electrodeL_left) )THEN
                     PDOS_raw_smeared(igrid,:) = PDOS_raw_smeared_L(igrid,:)
                  ELSE IF ( (Spatial_grid(igrid) <= electrodeR_right ).AND. ( Spatial_grid(igrid) >= electrodeR_left) )THEN
                       PDOS_raw_smeared(igrid,:) = PDOS_raw_smeared_R(igrid,:)
                ENDIF
        ENDDO 

   !OUTPUT PART
    PRINT*, "PRINTING OUTPUT"
    
   ! ==================================================================================!
   ! PRINTED in the gnuplot 3Dformat
   ! For visualization just launch gnuplot and type: splot 'name of the output file' with pm3d
   ! Example of Script:
   !set pm3d map                                  
   !set xrange [-17.5:17.5] 
   !set yrange [-5.0:5.0]   
   !set palette rgbformulae -23,-23,-3 
   !set palette rgbformulae -31,-31,-22
   !set cbrange [0.0:7.0]                         
   !set xlabel 'Position [Bohr]' font "Times-Roman, 26" 
   !set ylabel 'Energy [eV]' font "Times-Roman, 26" 
   !set title 'BAPA: Local Density of States' font "Times-Roman,40" 
   !set xtics font "Times-Roman, 26"                                   
   !set ytics font "Times-Roman, 26" 
   !set ylabel offset -6 
   !set xlabel offset 0,-3 
   !set title offset 0,3
   !set ytics 1                                   
   !set mytics 4 
   !set cblabel 'Log[LDOS]' font "Times-Roman, 26"            
   !set cbtics  font "Times-Roman, 26" 
   !set cblabel  offset 4,0 
   !splot 'LDOS_MolStilAcetyl_rawlog' with pm3d 
   ! ==================================================================================!

   OPEN ( 12, FILE=TRIM(fname_out)//"raw", FORM='formatted' )

   DO ie = 1, Nb_energy_out
       DO igrid = 1, nspatialgrid
                       IF ((Spatial_grid(igrid) >= print_left).AND.(Spatial_grid(igrid) <= print_right)) &
                       WRITE (12, '(3(f15.9))' )  (Spatial_grid(igrid) -((electrodeL_right + electrodeR_left)/2.0)  ) , (Out_energygrid(ie)- Fermi_energy) ,  ABS(PDOS_raw_smeared(igrid,ie))
       ENDDO
       WRITE (12, '(2(e15.9))' )
   ENDDO
   CLOSE( 12 )

       

   OPEN ( 12, FILE=TRIM(fname_out)//"rawlog", FORM='formatted' )

   DO ie = 1, Nb_energy_out
       DO igrid = 1, nspatialgrid
                       IF ((Spatial_grid(igrid) >= print_left).AND.(Spatial_grid(igrid) <= print_right)) &
                       WRITE (12, '(3(f15.9))' )  (Spatial_grid(igrid) -((electrodeL_right + electrodeR_left)/2.0)  ),  (Out_energygrid(ie)- Fermi_energy) , log(ABS(PDOS_raw_smeared(igrid,ie)) + 1.00)
       ENDDO
       WRITE (12, '(2(e15.9))' )
   ENDDO
   CLOSE( 12 )



    PRINT*, "OUTPUT printed"



 ! De-Allocation part       
        PRINT*, "Deallocating"
        ! Input grids
        DEALLOCATE (Spatial_grid) 
        DEALLOCATE (In_energygrid) 
        DEALLOCATE (Input_orb_PDOS)
        ! Output Grids
        DEALLOCATE (PDOS_C, PDOS_L, PDOS_R) 
        DEALLOCATE (PDOS_raw_smeared) 
        DEALLOCATE (PDOS_raw_smeared_C, PDOS_raw_smeared_L, PDOS_raw_smeared_R) 
        DEALLOCATE (Out_energygrid) 
        DEALLOCATE (Orb_center_coordinates) 
        PRINT*, "End of Deallocations"

        PRINT*, "End of Program"
END  program pdos_plot

