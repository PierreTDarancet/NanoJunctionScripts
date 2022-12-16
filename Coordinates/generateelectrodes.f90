    PROGRAM generate_electrodes
        IMPLICIT NONE


        ! 
        ! Parameters
        !

        ! Input Parameters
         !Number of atoms
         INTEGER, PARAMETER :: nb_Vect1=4
         INTEGER, PARAMETER :: nb_Vect2=4
         INTEGER, PARAMETER :: nb_Vect3=3
         INTEGER, PARAMETER :: norb=nb_Vect1*nb_Vect2*nb_Vect3*2

         CHARACTER*26, PARAMETER :: fname_in="Lead.xyz"

         ! Name of the input PDOS file
         ! Output Parameters
         ! Name of the output PDOS file
         CHARACTER*50, PARAMETER :: fname_out="Enforced_Leads.xyz"
         ! Number of points on the spatial grid
        CHARACTER*2    :: orbtype="Au"
        INTEGER, PARAMETER :: Input_orb_type=4

        REAL, PARAMETER :: Shift_Left_x = 3.030346  
        REAL, PARAMETER :: Shift_Left_y = 4.057368
        REAL, PARAMETER :: Shift_Left_z = -8.160021
        REAL, PARAMETER :: Shift_Right_x = 5.376043
        REAL, PARAMETER :: Shift_Right_y = 17.816868
        REAL, PARAMETER :: Shift_Right_z = 26.858720
        REAL, PARAMETER :: Vector_Left_1(2) =(/2.983990617,0.0/)
        REAL, PARAMETER :: Vector_Left_2(2) =(/1.491995308,2.584211679/)
        REAL, PARAMETER :: Vector_Right_1(2) =(/-2.983990617, 0.0/)
        REAL, PARAMETER :: Vector_Right_2(2) =(/-1.491995308,-2.584211679/)
        REAL, PARAMETER :: zshift = 2.436418136
        REAL, PARAMETER :: yshift =-1.722807785
        REAL, PARAMETER :: delta=0.5
  
         

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
        REAL, ALLOCATABLE :: Orb_Left_center_coordinates2(:,:)
        REAL, ALLOCATABLE :: Orb_Right_center_coordinates2(:,:)
        REAL, ALLOCATABLE :: Orb_Left_center_coordinates3(:,:)
        REAL, ALLOCATABLE :: Orb_Right_center_coordinates3(:,:)

        REAL ::      Maxd_left
        REAL ::     rmsleft

        REAL ::      Maxd_right
        REAL ::      rmsright
        INTEGER ::     iVect3 
        INTEGER ::       iVect1 
        INTEGER ::          iVect2
 
        LOGICAL :: FOUND 
        REAL :: orbx, orby, orbz
        ! Test Variables
        ! LOOP Variables
        INTEGER :: icoord, iat, iatl, iatr, iatcheck

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
        ALLOCATE ( Orb_Left_center_coordinates2(3,  (norb/2) ) )
        ALLOCATE ( Orb_Right_center_coordinates2(3,  (norb/2) ) )
        ALLOCATE ( Orb_Left_center_coordinates3(3,  (norb/2) ) )
        ALLOCATE ( Orb_Right_center_coordinates3(3,  (norb/2) ) )

        PRINT*, "End of Allocations"
 
        ! Reading Input File

        PRINT*, "Calculating Incoming POsition File"
 
        iat=1
       do iVect3 = 1, nb_Vect3
       do iVect1 = 1, nb_Vect1
          do iVect2 = 1, nb_Vect2
            Orb_Left_center_coordinates(1,iat)= ((iVect1-1) *  Vector_Left_1(1)) + ((iVect2-1) *  Vector_Left_2(1)) + Shift_Left_x 
            Orb_Left_center_coordinates(2,iat)= ((iVect1-1) *  Vector_Left_1(2)) + ((iVect2-1) *  Vector_Left_2(2)) + Shift_Left_y + ( yshift* (iVect3-1) ) 
            Orb_Left_center_coordinates(3,iat)= ((iVect3-1) *  (-zshift))  + Shift_Left_z 
            Orb_Right_center_coordinates(1,iat)= ((iVect1-1) *  Vector_Right_1(1)) + ( (iVect2-1) *  Vector_Right_2(1) ) + Shift_Right_x 
            Orb_Right_center_coordinates(2,iat)= ((iVect1-1) *  Vector_Right_1(2)) + ((iVect2-1) *  Vector_Right_2(2))  + Shift_Right_y - ( yshift* (iVect3-1) ) 
            Orb_Right_center_coordinates(3,iat)= ((iVect3-1) *  (zshift) )+ Shift_Right_z 
            IF (iVect2 > 2) THEN
                    Orb_Left_center_coordinates(1,iat)=Orb_Left_center_coordinates(1,iat) -  Vector_Left_1(1) 
                    Orb_Left_center_coordinates(2,iat)=Orb_Left_center_coordinates(2,iat) -  Vector_Left_1(2) 
                    Orb_Right_center_coordinates(1,iat)=Orb_Right_center_coordinates(1,iat) - Vector_Right_1(1)
                    Orb_Right_center_coordinates(2,iat)=Orb_Right_center_coordinates(2,iat) - Vector_Right_1(2)
            ENDIF
            iat=iat+1
     ENDDO
     ENDDO      
     ENDDO
   !OUTPUT PART
    PRINT*, "PRINTING OUTPUT"
    
   OPEN ( 12, FILE=TRIM(fname_out), FORM='formatted' )
      iatl = 1
      iatr = 1

   DO iat = 1, norb/2
     orbx =  Orb_Left_center_coordinates(1,iatl) 
     orby =  Orb_Left_center_coordinates(2,iatl) 
     orbz =  Orb_Left_center_coordinates(3,iatl)
     WRITE (12, '(a11,3(f15.9), i4)' ) orbtype,  orbx  , orby, orbz , Input_orb_type
     orbx =  Orb_Right_center_coordinates(1,iatr) 
     orby =  Orb_Right_center_coordinates(2,iatr)
     orbz =  Orb_Right_center_coordinates(3,iatr) 
     WRITE (12, '(a11,3(f15.9), i4)' ) orbtype,  orbx  , orby, orbz , Input_orb_type
     iatl =    iatl + 1
     iatr =    iatr + 1
   ENDDO


   CLOSE( 12 )
   ! 
   PRINT*, "OUTPUT printed"


        PRINT*, "READING Incoming POsition File"
      open(100,file=trim(fname_in),form='formatted')
      rewind(100)
      iatl = 1
      iatr = 1
      do iat = 1, (norb/2)
       read(100,*)chr, Orb_Left_center_coordinates2(1,iatl), Orb_Left_center_coordinates2(2,iatl), &
     &  Orb_Left_center_coordinates2(3,iatl)
       iatl =       iatl + 1  
        PRINT*, iat
       read(100,*)chr, Orb_Right_center_coordinates2(1,iatr), Orb_Right_center_coordinates2(2,iatr), &
     &  Orb_Right_center_coordinates2(3,iatr)
       iatr =       iatr + 1  
      end do
       close(100)
        PRINT*, "END of READING"



        PRINT*, "Check_left"

      do iatcheck=1,norb/2
        FOUND = .FALSE.
        do iat = 1, norb/2
          IF ( (Orb_Left_center_coordinates2(1,iat) <= Orb_Left_center_coordinates(1,iatcheck)+delta ).AND.(Orb_Left_center_coordinates2(1,iat) >= Orb_Left_center_coordinates(1,iatcheck)- delta  ) ) THEN
            IF ( (Orb_Left_center_coordinates2(2,iat) <= Orb_Left_center_coordinates(2,iatcheck)+delta ).AND.(Orb_Left_center_coordinates2(2,iat) >= Orb_Left_center_coordinates(2,iatcheck)- delta  ) ) THEN
               IF ( (Orb_Left_center_coordinates2(3,iat) <= Orb_Left_center_coordinates(3,iatcheck)+delta ).AND.(Orb_Left_center_coordinates2(3,iat) >= Orb_Left_center_coordinates(3,iatcheck)- delta  ) ) THEN   
                   Orb_Left_center_coordinates3(1,iatcheck)= Orb_Left_center_coordinates2(1,iat)
                   Orb_Left_center_coordinates3(2,iatcheck)= Orb_Left_center_coordinates2(2,iat)
                   Orb_Left_center_coordinates3(3,iatcheck)= Orb_Left_center_coordinates2(3,iat)
                   FOUND = .TRUE.
               ENDIF
            ENDIF
         ENDIF
        enddo
        IF (.NOT.FOUND) PRINT*, "WARNING!!!!!!!!!!! No correspondance found for left electrode iat=", iatcheck
      end do

      PRINT*, "Check_right"

      do iatcheck=1,norb/2
        FOUND = .FALSE.
        do iat = 1, norb/2
          IF ( (Orb_Right_center_coordinates2(1,iat) <= Orb_Right_center_coordinates(1,iatcheck)+delta ).AND.(Orb_Right_center_coordinates2(1,iat) >= Orb_Right_center_coordinates(1,iatcheck)- delta  ) ) THEN
            IF ( (Orb_Right_center_coordinates2(2,iat) <= Orb_Right_center_coordinates(2,iatcheck)+delta ).AND.(Orb_Right_center_coordinates2(2,iat) >= Orb_Right_center_coordinates(2,iatcheck)- delta  ) ) THEN
               IF ( (Orb_Right_center_coordinates2(3,iat) <= Orb_Right_center_coordinates(3,iatcheck)+delta ).AND.(Orb_Right_center_coordinates2(3,iat) >= Orb_Right_center_coordinates(3,iatcheck)- delta  ) ) THEN   
                   Orb_Right_center_coordinates3(1,iatcheck)= Orb_Right_center_coordinates2(1,iat)
                   Orb_Right_center_coordinates3(2,iatcheck)= Orb_Right_center_coordinates2(2,iat)
                   Orb_Right_center_coordinates3(3,iatcheck)= Orb_Right_center_coordinates2(3,iat)
                   FOUND = .TRUE.
               ENDIF
            ENDIF
         ENDIF
        enddo
        IF (.NOT.FOUND) PRINT*, "WARNING!!!!!!!!!!! No correspondance found for right electrode iat=", iatcheck
      end do

      PRINT*, "Check Differences"
      Maxd_right=0.0
      Maxd_left=0.0
      rmsright=0.0
      rmsleft=0.0

      do iatcheck=1,norb/2
           rmsright= rmsright+ (Orb_Right_center_coordinates3(1,iatcheck) - Orb_Right_center_coordinates(1,iatcheck) )**2  + & 
                           (Orb_Right_center_coordinates3(2,iatcheck) - Orb_Right_center_coordinates(2,iatcheck) )**2  + & 
                           (Orb_Right_center_coordinates3(1,iatcheck) - Orb_Right_center_coordinates(1,iatcheck) )**2  

           rmsleft= rmsleft  +  (Orb_Left_center_coordinates3(1,iatcheck) - Orb_Left_center_coordinates(1,iatcheck) )**2  + & 
                           (Orb_Left_center_coordinates3(2,iatcheck) - Orb_Left_center_coordinates(2,iatcheck) )**2  + & 
                           (Orb_Left_center_coordinates3(1,iatcheck) - Orb_Left_center_coordinates(1,iatcheck) )**2  

           Maxd_right = MAX(Maxd_right, ABS(Orb_Right_center_coordinates3(1,iatcheck) - Orb_Right_center_coordinates(1,iatcheck) ) )
           Maxd_right = MAX(Maxd_right, ABS(Orb_Right_center_coordinates3(2,iatcheck) - Orb_Right_center_coordinates(2,iatcheck) ) )
           Maxd_right = MAX(Maxd_right, ABS(Orb_Right_center_coordinates3(3,iatcheck) - Orb_Right_center_coordinates(3,iatcheck) ) )
           Maxd_left = MAX(Maxd_left, ABS(Orb_Left_center_coordinates3(1,iatcheck) - Orb_Left_center_coordinates(1,iatcheck) ) )
           Maxd_left = MAX(Maxd_left, ABS(Orb_Left_center_coordinates3(2,iatcheck) - Orb_Left_center_coordinates(2,iatcheck) ) )
           Maxd_left = MAX(Maxd_left, ABS(Orb_Left_center_coordinates3(3,iatcheck) - Orb_Left_center_coordinates(3,iatcheck) ) )

      enddo
     PRINT*, "Max Diff Left", Maxd_left
     PRINT*, "RMS Left", sqrt(rmsleft)

     PRINT*, "Max Diff Right", Maxd_right
     PRINT*, "RMS Right", sqrt(rmsright)

     ! De-Allocation part       
        PRINT*, "Deallocating"
        ! Input grids
         ! Output Grids
        DEALLOCATE (Orb_Left_center_coordinates,Orb_Left_center_coordinates2,Orb_Left_center_coordinates3) 
        DEALLOCATE (Orb_Right_center_coordinates,Orb_Right_center_coordinates3,Orb_Right_center_coordinates2) 
        PRINT*, "End of Deallocations"

        PRINT*, "End of Program"
END  program generate_electrodes

