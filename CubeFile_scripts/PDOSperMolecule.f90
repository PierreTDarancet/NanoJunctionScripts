    PROGRAM pdos_plot
        IMPLICIT NONE
   INTEGER, PARAMETER :: dimensions = 3    
   INTEGER, PARAMETER :: iufile=100
   CHARACTER*100, PARAMETER :: datafile="pdosplot.in"


   INTEGER :: nbenergies, nbmolecules, nborbitals, nbatoms !,nbenergiesout
   INTEGER, ALLOCATABLE :: atomindex(:,:)
   REAL*8, ALLOCATABLE :: energybroadening(:), In_energygrid(:)
   REAL*8, ALLOCATABLE :: Input_orb_PDOS(:,:), PDOSpermolecule(:,:), PDOSpermolecule_int(:,:)  ! Out_energygrid(:), 
   REAL*8 ::  energy_broadening, deltaE

   CHARACTER*100 :: FileInput, FileOutput
  
    CHARACTER*100,        ::  chr, chr2
     INTEGER :: imol, iint, ipos,iatom, imol2, energyinterpol, ie, nbpoints, iat
        
        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInput
        PRINT*, "...READING Incoming Data File: Input file :", Trim(FileInput)
        READ(iufile,*) nbenergies 
        PRINT*, "...READING Incoming Data File: nbenergies:", nbenergies
        READ(iufile,*) nbatoms
        PRINT*, "...READING Incoming Data File: nbatoms:", nbatoms

        READ(iufile,*) nbmolecules
        PRINT*, "...READING Incoming Data File: nbmolecules:", nbmolecules
        PRINT*, "...ALLOCATE: atomindex:"
        
        ALLOCATE(atomindex(nbmolecules, 2) ) !Starting indices and ending indices of atoms 
        PRINT*, "...ALLOCATE: energybroadening:"

        ALLOCATE(energybroadening(nbmolecules))
        PRINT*, "...ALLOCATE: In_energygrid :"

        ALLOCATE(In_energygrid(nbenergies))
        PRINT*, "...READING Incoming Data File: atomindex :"

        do imol = 1,nbmolecules
                READ(iufile,*) atomindex(imol,1), atomindex(imol,2), energybroadening(imol)
                PRINT*, "...indices of the atoms for molecule:", imol, atomindex(imol,:), energybroadening(imol)
            
        end do
        READ(iufile,*) FileOutput
        PRINT*, "...READING Incoming Data File: FileOutput                      :", FileOutput

        !READ(iufile,*) nbenergiesout
        !PRINT*, "...READING Incoming Data File: nbenergiesout                   :", nbenergiesout
       !ALLOCATE (    PDOSpermolecule(nbmolecules,nbenergiesout) ) 
       ALLOCATE (    PDOSpermolecule(nbmolecules,nbenergies) ) 
       ALLOCATE (     PDOSpermolecule_int(nbmolecules,nbenergies) ) 
       PRINT*, "...READING Incoming Data File: Closing    "

        CLOSE(iufile)
        PRINT*, "...Done"
        ! 
        ! Parameters
        !
        
        
      ! Reading Input DOS File

        PRINT*, "READING Incoming PDOS File", trim(FileInput)
      open(100,file=trim(FileInput),form='formatted')
      rewind(100)
      read(100,'(A)')chr
      read(100,'(A)')chr

      read(100,'(A)')chr !      <norbitals>6576</norbitals>
      !PRINT*, chr, LEN_TRIM(chr)
      !PRINT*, TRIM(chr)
      chr2='</norbitals>'
      ipos=index(chr, TRIM(chr2))
      IF ( ipos < 2 ) THEN
         PRINT*, "Problem in reading norbital", ipos, chr, chr2
         STOP
      ENDIF
      chr2=TRIM(chr(12:ipos-1))
      read ( chr2,'(I10)') nborbitals
      PRINT*, "READING Incoming PDOS File: nborbitals", nborbitals
      PRINT*, "Allocating PDOS"
      ALLOCATE (Input_orb_PDOS(nbatoms,nbenergies))
      read(100,'(A)') chr
      PRINT*, "READING Incoming PDOS File: energy grid"

      do ie = 1,nbenergies
        read(100,*) In_energygrid(ie)
      end do
      read(100,'(A)')chr
      
      PRINT*, "READING Incoming PDOS File: energy grid ... done!"
      
      do iat = 1,nborbitals
       read(100,'(A)')chr !<orbital
       read(100,'(A)')chr ! index="                        1" ORBITAL INDEX
       read(100,'(A)')chr ! atom_index="                        1"
       ipos=LEN_TRIM(chr)
       !PRINT*, ipos, chr(14:ipos-1), chr
       chr2=chr(14:ipos-1)
       chr=TRIM(adjustl(chr2))
       read ( chr,'(I10)') iatom
       PRINT*, "READING Incoming PDOS File: atom index", iatom
       read(100,'(A)')chr ! species="C"
       read(100,'(A)')chr !! position=" -17.299902 -12.124205  37.578357"
       !       read(100,"(a11,3f11.6,a1)")chr, Orb_center_coordinates(1,iat), Orb_center_coordinates(2,iat), &
       !     &  Orb_center_coordinates(3,iat),chr11
       ! print*,trim(chr),trim(chr11),Orb_center_coordinates(1:3,iat)
       read(100,'(A)')chr !! n="                        2"
       read(100,'(A)')chr !! l="                        0"
       read(100,'(A)')chr !! m="                        0"
       read(100,'(A)')chr !! z="                        1"
       read(100,'(A)')chr !!>
       read(100,'(A)')chr !!<data>
       do ie = 1, nbenergies
           read(100,*) Input_orb_PDOS(iatom,ie)
       end do
       read(100,'(A)')chr        !</data>
       read(100,'(A)')chr !</orbital>
      end do
      close(100)
        PRINT*, "END of READING"
   ! End of reading

    ! Energy Grid 
    !    PRINT*, "...Energy Grid Initialization"
    !   ! Energy separation on the output energy grid
    ! Delta_energy_output_energy_grid = (In_energygrid(nbenergies) - In_energygrid(1)) / REAL(nbenergiesout-1)
    ! ALLOCATE(Out_energygrid(nbenergiesout))
    !   ! Energy Grid calculation   
    ! DO ie = 1, Nb_energy_out
    !   Out_energygrid(ie) =Output_Energy_grid_min + REAL(ie -1) * Delta_energy_output_energy_grid
    ! END DO    

     !   PRINT*, "End of Output Grids Initialization"
    ! End of Grids Construction

    ! OUTPUT PDOS INIT
    PRINT*, "PDOS Initialization"
    PDOSpermolecule(:,:)= 0.00000000   
    
    PRINT*, "PDOS ... Atoms to molecules"
    
     DO  iatom = 1, nbatoms
          imol=0
          DO imol2 = 1, nbmolecules
              IF ((atomindex(imol2,2) >= iatom).AND.(atomindex(imol2,1) <= iatom) )  THEN
                 imol=imol2
              END IF
          END DO
          PDOSpermolecule(imol,:) = PDOSpermolecule(imol,:) +  Input_orb_PDOS(iatom,:)
    END DO
    
    
    PRINT*, "PDOS ... broadening"
 
    PDOSpermolecule_int(:,:)= 0.00000000   
    DO imol=1, nbmolecules
         PRINT*, "PDOS ... broadening, mol", imol

         energy_broadening = energybroadening(imol) 
         energyinterpol = 16*INT(energy_broadening/(In_energygrid(2) - In_energygrid(1))) ! integrate over a subset of the energy values. 
         PRINT*, energyinterpol
         DO ie=1, nbenergies
            !PRINT*, "ie", ie
            nbpoints=0 ! counter to renormalize if points are missing
            DO iint=1, energyinterpol
               IF ( (ie-iint) > 1)  THEN 
                   !PRINT*, "-", (ie-iint)
                   deltaE=(In_energygrid(ie-iint) - In_energygrid(ie))
                   PDOSpermolecule_int(imol,ie) = PDOSpermolecule_int(imol,ie) + (  PDOSpermolecule(imol,ie-iint)   * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) *  exp(- ( (deltaE**2) /( 2* energy_broadening**2) )    ))  
               !ELSE 
                   nbpoints=nbpoints+1    
               END IF
               IF ( (ie+iint) < nbenergies + 1 )  THEN 
                   !PRINT*, "+", (ie+iint), (nbenergies+1)
                   deltaE=(In_energygrid(ie+iint) - In_energygrid(ie))
                   PDOSpermolecule_int(imol,ie) = PDOSpermolecule_int(imol,ie) + (  PDOSpermolecule(imol,ie+iint)   * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) *  exp(- ( (deltaE**2) /( 2* energy_broadening**2) )    ))  
                   nbpoints=nbpoints+1      
               END IF
           END DO
               
           PDOSpermolecule_int(imol,ie) =  PDOSpermolecule_int(imol,ie) +  ( PDOSpermolecule(imol,ie)   * ( 1.00/ sqrt(2 * 3.1415 * ( energy_broadening**2) )) )
           ! ADD RENORMALIZATION FACTOR on the exponential there to correct for missing terms
               
         END DO
    END DO
         
   
    ! PRINT

      
      OPEN ( 12, FILE=TRIM(FileOutput)//".grid", FORM='formatted' )
      DO ie = 1, nbenergies
          WRITE (12, * ) In_energygrid(ie)
      END DO     
      CLOSE(12)
      
      DO imol=1, nbmolecules
         WRITE(chr2,'(I10)') imol
         PRINT*, TRIM(ADJUSTL(chr2))
         WRITE(chr,*) TRIM(FileOutput)//"_",TRIM(ADJUSTL(chr2))
         OPEN ( 12, FILE=TRIM(chr), FORM='formatted' )
         DO ie = 1, nbenergies
              WRITE (12, '(2(f15.9))' )   PDOSpermolecule(imol,ie), PDOSpermolecule_int(imol,ie)
         END DO
         CLOSE(12)
      END DO
      PRINT*, "DEALLOCATING"
      DEALLOCATE (    PDOSpermolecule ) 
      DEALLOCATE (     PDOSpermolecule_int ) 
      DEALLOCATE (Input_orb_PDOS)
      DEALLOCATE(atomindex ) !Starting indices and ending indices of atoms 
      DEALLOCATE(energybroadening)
      DEALLOCATE(In_energygrid )


 
      END PROGRAM
   
   