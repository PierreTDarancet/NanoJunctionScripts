
! Copyright (C) 2009 Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET

!***********************************************
   PROGRAM dosbroadening
   !*********************************************** 
        IMPLICIT NONE


        ! 
        ! Parameters
        !
 
         ! Input Parameters
         !Number of energies in the Siesta PDOS file
         INTEGER, PARAMETER :: Nb_energy_in=3000
         ! Name of the input PDOS file
         CHARACTER*50, PARAMETER :: fname_in="" 
         !

         ! Output Parameters
         !Number of energies in the Output PDOS file
         INTEGER, PARAMETER :: Nb_energy_out=3000 
         !Minimum energy of the energy spectrum of the output PDOS file
         REAL, PARAMETER :: Output_Energy_grid_min=-20.0
         !Maximum energy of the energy spectrum of the output PDOS file
         REAL, PARAMETER :: Output_Energy_grid_max=10.0
         ! Name of the output PDOS file
         CHARACTER*50, PARAMETER :: fname_out=""
         ! Number of points on the spatial grid
         INTEGER, PARAMETER :: energyinterpol=100
         REAL, PARAMETER :: energy_broadening=0.01
         REAL, PARAMETER :: ZERO=0.000000000000000

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
        ! Useless things
        CHARACTER*100  :: chr
        CHARACTER*1    :: chr11
        ! Important things
        ! Input energy grid
        REAL, ALLOCATABLE :: In_energygrid(:)
        ! Input PDOS
        REAL, ALLOCATABLE :: Input_orb_PDOS(:)
        ! OUTPUT related VARIABLES       
        REAL, ALLOCATABLE :: PDOS_out(:) 
        ! OUTPUT PDOS(z,Nb_energy_out)
        REAL, ALLOCATABLE :: PDOS_out_smeared(:)

        ! OUTPUT Energy grid (Nb_energy_out)
        REAL, ALLOCATABLE :: Out_energygrid(:)
        ! Delta energy of the Output energy Grid
        REAL :: Delta_energy_output_energy_grid
        REAL :: Delta_energy_input_energy_grid
        ! Interpolation variable
        REAL :: summation_out 
        ! Test Variables

        ! LOOP Variables
        INTEGER :: ie, ie2, inter2, nbpoint

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
        ALLOCATE (In_energygrid(Nb_energy_in)  ) 
        ALLOCATE (Input_orb_PDOS(Nb_energy_in))
        ! Output Grids
        ALLOCATE (PDOS_out(Nb_energy_out)  ) 
        ALLOCATE (Out_energygrid(Nb_energy_out)  ) 
        ALLOCATE (PDOS_out_smeared(Nb_energy_out))
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
      !do iat = 1,norb
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
      ! read(100,"(a11,3f11.6,a1)")chr, Orb_center_coordinates(1,iat), Orb_center_coordinates(2,iat), &
     !&  Orb_center_coordinates(3,iat),chr11
      ! print*,trim(chr),trim(chr11),Orb_center_coordinates(1:3,iat)
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       do ie = 1, Nb_energy_in
           read(100,*)Input_orb_PDOS(ie)
       end do
       read(100,*)chr
       read(100,*)chr
     ! end do
      close(100)
        PRINT*, "END of READING"

   ! End of reading


   ! Grid Constuction
        PRINT*, "Constructing Output Grids:"
 
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
    PDOS_out(:)= ZERO
    PDOS_out_smeared(:)= ZERO

    ! OUTPUT PDOS calculation

        PRINT*, "PDOS Calculation"

                 ! Energy loop
                DO ie = 1, Nb_energy_in
                      ! Is energy on the output energy grid ?
                      IF ( ( In_energygrid(ie) <= Output_Energy_grid_max ).AND. ( In_energygrid(ie) > Output_Energy_grid_min) ) THEN
                         ! Yes, Enter contribution
                         DO ie2 = 1, Nb_energy_out-1
                             IF ( ( In_energygrid(ie) <= (Out_energygrid(ie2)+Delta_energy_output_energy_grid) ).AND. ( In_energygrid(ie) > Out_energygrid(ie2)) )THEN
                                PDOS_out(ie2) = PDOS_out(ie2) + ( Input_orb_PDOS(ie) * Delta_energy_input_energy_grid/ Delta_energy_output_energy_grid  ) 
                             ENDIF
                         ENDDO
                         !
                     ENDIF
                     ! No, nothing happens  
                ENDDO
        PRINT*, "End of PDOS Calculation"






        PRINT*, "Fancy Convolutions and Filtering"
   ! ==================================================================================!
   ! PUT ALL THE FANCY FILTERING, CONVOLUTION, SMEARING HERE !!!!!!!!!!!!!!!!!!!!!!!!
   ! ==================================================================================!

   PDOS_out_smeared(:) = ZERO

         DO ie = 1, Nb_energy_out
             summation_out= ZERO
             nbpoint=0
             DO inter2 = 1, energyinterpol                  
                  IF (ie-inter2 > 0) THEN  
                     summation_out  = summation_out + PDOS_out(ie-inter2) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2*  energy_broadening**2) )  
                      nbpoint=nbpoint+1
                 ENDIF
                 IF (ie+inter2 <= Nb_energy_out) THEN  
                     summation_out  = summation_out + PDOS_out(ie+inter2) *( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) * exp(- ( (inter2*Delta_energy_output_energy_grid)**2) /( 2* energy_broadening**2) )  
                      nbpoint=nbpoint+1
                 ENDIF
            ENDDO
            PDOS_out_smeared(ie) = ( PDOS_out(ie) * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) ))   + summation_out ) * ( REAL( (2*energyinterpol+1) / nbpoint ))
 
         ENDDO




   !OUTPUT PART
    PRINT*, "PRINTING OUTPUT"
    
   ! ==================================================================================!
   ! PRINTED in the gnuplot 3Dformat
   ! For visualization just launch gnuplot and type: splot 'name of the output file' with pm3d
   ! ==================================================================================!

   OPEN ( 12, FILE=TRIM(fname_out), FORM='formatted' )
   DO ie = 1, Nb_energy_out
                       WRITE (12, '(3(f15.9))' ) Out_energygrid(ie) , PDOS_out_smeared(ie),  PDOS_out(ie)
   ENDDO
   CLOSE( 12 )



    PRINT*, "OUTPUT printed"



 ! De-Allocation part       
        PRINT*, "Deallocating"
        ! Input grids
        DEALLOCATE (In_energygrid) 
        DEALLOCATE (Input_orb_PDOS)
        ! Output Grids
        DEALLOCATE (PDOS_out) 
        DEALLOCATE (PDOS_out_smeared) 
        DEALLOCATE (Out_energygrid) 
        PRINT*, "End of Deallocations"

        PRINT*, "End of Program"
END  program    dosbroadening

