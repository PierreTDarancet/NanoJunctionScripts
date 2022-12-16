    PROGRAM shift_electrodes
        IMPLICIT NONE


        ! 
        ! Parameters
        !

        ! Input Parameters
         !Number of atoms
         INTEGER, PARAMETER :: norb=198
         ! Name of the input PDOS file
         CHARACTER*26, PARAMETER :: fname_in="Electrodes2.xyz"
         ! Output Parameters
         ! Name of the output PDOS file
         CHARACTER*21, PARAMETER :: fname_out="testshift.out"
         ! Number of points on the spatial grid
        CHARACTER*2    :: orbtype="Au"
        INTEGER, PARAMETER :: Input_orb_type=6
        REAL, PARAMETER :: Shift_Left_x = 1.514869
        REAL, PARAMETER :: Shift_Left_y = 4.42232193
        REAL, PARAMETER :: Shift_Left_z = 22.17995

        REAL, PARAMETER :: Shift_Right_x = -5.48557 !15.043663
        REAL, PARAMETER :: Shift_Right_y =-4.59620907
        REAL, PARAMETER :: Shift_Right_z =14.939312

        ! Local Variables 

        ! INPUT related VARIABLES        
        ! Useless things
        CHARACTER*100  :: chr
        CHARACTER*2  :: chr11

        ! Important things
        ! Input PDOS

        ! Input Orbital center positions
        REAL, ALLOCATABLE :: Orb_Left_center_coordinates(:,:)
        REAL, ALLOCATABLE :: Orb_Right_center_coordinates(:,:)

        REAL :: orbx, orby, orbz
        ! Test Variables
        ! LOOP Variables
        INTEGER :: icoord, iat, iatl, iatr

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
        ALLOCATE ( Orb_Left_center_coordinates(3,  (norb/2) ) )
        ALLOCATE ( Orb_Right_center_coordinates(3,  (norb/2) ) )

        PRINT*, "End of Allocations"
 
        ! Reading Input File

        PRINT*, "READING Incoming PDOS File"
      open(100,file=trim(fname_in),form='formatted')
      rewind(100)
      iatl = 1
      iatr = 1
      do iat = 1, (norb/2)
       read(100,*)chr, Orb_Left_center_coordinates(1,iatl), Orb_Left_center_coordinates(2,iatl), &
     &  Orb_Left_center_coordinates(3,iatl),chr11
       iatl =       iatl + 1  

       read(100,*)chr, Orb_Right_center_coordinates(1,iatr), Orb_Right_center_coordinates(2,iatr), &
     &  Orb_Right_center_coordinates(3,iatr),chr11
       iatr =       iatr + 1  
       end do
       close(100)
        PRINT*, "END of READING"

   !OUTPUT PART
    PRINT*, "PRINTING OUTPUT"
    
   OPEN ( 12, FILE=TRIM(fname_out), FORM='formatted' )
      iatl = 1
      iatr = 1

   DO iat = 1, norb/2
     orbx =  Orb_Left_center_coordinates(1,iatl) + Shift_Left_x 
     orby =  Orb_Left_center_coordinates(2,iatl) + Shift_Left_y 
     orbz =  Orb_Left_center_coordinates(3,iatl) + Shift_Left_z 
     WRITE (12, '(a11,3(f15.9), i4)' ) orbtype,  orbx  , orby, orbz , Input_orb_type
     orbx =  Orb_Right_center_coordinates(1,iatr) + Shift_Right_x 
     orby =  Orb_Right_center_coordinates(2,iatr) + Shift_Right_y 
     orbz =  Orb_Right_center_coordinates(3,iatr) + Shift_Right_z 
     WRITE (12, '(a11,3(f15.9), i4)' ) orbtype,  orbx  , orby, orbz , Input_orb_type
     iatl =    iatl + 1
     iatr =    iatr + 1
   ENDDO
   CLOSE( 12 )
   ! 
   PRINT*, "OUTPUT printed"


     ! De-Allocation part       
        PRINT*, "Deallocating"
        ! Input grids
         ! Output Grids
        DEALLOCATE (Orb_Left_center_coordinates) 
        DEALLOCATE (Orb_Right_center_coordinates) 
        PRINT*, "End of Deallocations"

        PRINT*, "End of Program"
END  program shift_electrodes

