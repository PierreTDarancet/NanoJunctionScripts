!
! Copyright (C) 2011 The Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!***********************************************
   PROGRAM WaveFunctionInPotentialDifference
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   USE IO_module,                              ONLY : gnuplotwrite, datwrite, ReadInputFileWFinPotential, GaussianFilereadHeader, GaussianFilewrite, GaussianFileread, PotentialFilereadHeader,  PotentialFileread,  XYZFilereadHeader, XYZFileread
   USE calculatequantities_Module,         ONLY : CalculateFirstOrderEnergy, CalculateVolume, CalculateCenterOfCharge, CalculateDipoleOperatorWF

   IMPLICIT NONE
   INTEGER, PARAMETER :: dimensions = 3    
   INTEGER, PARAMETER :: iufile=100
   REAL(dbl), PARAMETER ::  Cutoffdistance=20004.5
   REAL(dbl), PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   REAL(dbl), PARAMETER ::  Ryd_to_eV = 13.605698066*ONE
   CHARACTER*100, PARAMETER :: datafile="DATA.in"


   !
   ! Input variables
   !
   CHARACTER*100 ::  FileInput0bias
   CHARACTER*100 ::  FileInputFinitebias
   CHARACTER*100 ::  FileInputGrid0bias
   CHARACTER*100 ::  FileInputGridFinitebias
   CHARACTER*100 ::  FileInputxyz
   CHARACTER*100 ::  FileInputWF
   CHARACTER*100 ::  FileOutputGaussian
   CHARACTER*100 ::  FileOutputGaussian0bias="rho0.cube"
   CHARACTER*100 ::  FileOutputGaussianFiniteBias="rhofinite.cube"
   CHARACTER*100 ::  FileOutputGnuplotx
   CHARACTER*100 ::  FileOutputGnuploty
   CHARACTER*100 ::  FileOutputDat
   CHARACTER*100 ::  FileOutputDat0bias
   CHARACTER*100 ::  FileOutputDatFinitebias

   REAL(dbl) ::   cell(dimensions,dimensions)
   REAL(dbl) ::   cell_check(dimensions,dimensions)
   REAL(dbl) ::   X0Cell(dimensions)
   REAL(dbl) ::   X0Cell_check(dimensions)
   INTEGER   ::   mesh(dimensions)
   INTEGER   ::   mesh_check(dimensions)
   INTEGER   ::   Nb_atoms
   INTEGER   ::   Nb_atoms_check

   REAL(dbl), ALLOCATABLE ::   field_0bias(:,:,:)
   REAL(dbl), ALLOCATABLE ::   field_Finitebias(:,:,:)
   REAL(dbl), ALLOCATABLE ::   Wavefunction(:,:,:)
   INTEGER, ALLOCATABLE   ::   atoms_type(:)
   REAL(dbl), ALLOCATABLE ::   atoms_position(:,:)
   INTEGER, ALLOCATABLE   ::   atoms_type_check(:)
   REAL(dbl), ALLOCATABLE ::   atoms_position_check(:,:)


   !
   ! Output variables
   !
   REAL(dbl), ALLOCATABLE ::   target_field(:,:,:)
   REAL(dbl), ALLOCATABLE ::   Integrated_field1(:,:)
   REAL(dbl), ALLOCATABLE ::   Integrated_field2(:,:)
   REAL(dbl), ALLOCATABLE ::   Integrated_field_xy(:) 
   REAL(dbl) :: FirstOrderCorrection 
   REAL(dbl) :: DipoleOperator(dimensions)  ! dipole operator associ
   REAL(dbl) :: Integrated_Density  ! dipole operator associ
   !REAL(dbl)              :: PositionXYZ(dimensions)  ! position in the cell
   REAL(dbl)              :: CenterOfCharge(dimensions)  ! center of the charges
   REAL(dbl)              :: Integrated_volume, dV
   REAL(dbl)              :: Check_Volume
   REAL(dbl) ::   Nb_electrons ! Valence electrons only
   ! 
   ! Local variables
   ! 
   INTEGER :: ia, ib, ic, idimension
   REAL(dbl), ALLOCATABLE  :: vect1(:), vect2(:)
!***********************************************
!
!------------------------------
! main body
!------------------------------
!  !
   ! read input General Information
   !
   CALL ReadInputFileWFinPotential(iufile, datafile, FileInputGrid0bias, FileInputGridFinitebias, FileInput0bias, FileInputFinitebias, FileInputxyz , FileInputWF,  FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty,FileOutputDat,FileOutputDat0bias,FileOutputDatFinitebias)



!  !
   ! read input file xyz
   !

   CALL XYZFilereadHeader(FileInputxyz, iufile, Nb_atoms)

   ALLOCATE( atoms_type( Nb_atoms) )   
   ALLOCATE(atoms_position( Nb_atoms, dimensions) )   

   CALL XYZFileread( FileInputxyz, iufile, dimensions, Nb_atoms, atoms_type, atoms_position)


!  !
   ! read input file 0 bias
   !
   CALL PotentialFilereadHeader( FileInputGrid0bias, iufile, Cutoffdistance, dimensions, mesh, X0Cell, cell)

   ALLOCATE( field_0bias(mesh(1),mesh(2),mesh(3)) )   
   field_0bias(:,:,:) = ZERO 
   CALL PotentialFileread( FileInput0bias, iufile,Cutoffdistance, dimensions,  mesh, field_0bias)
   field_0bias(:,:,:) =    field_0bias(:,:,:)
   CALL GaussianFilewrite(FileOutputGaussian0bias, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, field_0bias)

!  !
   ! read input file Finite bias
   !
   CALL PotentialFilereadHeader(FileInputGridFinitebias, iufile, Cutoffdistance, dimensions, mesh_check, X0Cell, cell)
   IF ( ( ABS(mesh_check(1)- mesh(1)) >  EPS_m5 ).OR.( ABS(mesh_check(2) - mesh(2)) >  EPS_m5 ).OR.( ABS(mesh_check(3)-mesh(3))>  EPS_m5 ) ) THEN
       PRINT*, "Discrepancies in the meshes" 
       PRINT*, mesh_check(1), mesh(1)
       PRINT*, mesh_check(2), mesh(2) 
       PRINT*, mesh_check(3), mesh(3) 
      STOP
   ENDIF

   ALLOCATE( field_Finitebias(mesh(1),mesh(2),mesh(3)) )   
   field_Finitebias(:,:,:) = ZERO 
   CALL PotentialFileread( FileInputFinitebias, iufile,Cutoffdistance, dimensions,  mesh,  field_Finitebias)

   field_Finitebias(:,:,:) =    field_Finitebias(:,:,:)
   CALL GaussianFilewrite(FileOutputGaussianFinitebias, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, field_Finitebias)


!  !
   ! Calculate the Difference
   !
   ALLOCATE( target_field(mesh(1),mesh(2),mesh(3)) )   
   target_field(:,:,:) = ZERO
   target_field(:,:,:) = field_Finitebias(:,:,:) - field_0bias(:,:,:)

!  !
   ! Print the Difference in a Gaussian file
   !

   CALL GaussianFilewrite(FileOutputGaussian, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, target_field)

!  !
   !Difference  Average   on x and y lines
   !
   ALLOCATE( Integrated_field1(mesh(2),mesh(3)) )
   ALLOCATE( Integrated_field2(mesh(1),mesh(3)) )
   Integrated_field1(:,:) = ZERO
   Integrated_field2(:,:) = ZERO
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_field1(ib,ic) =  Integrated_field1(ib,ic) +  target_field(ia,ib,ic) / mesh(1)
		   ENDDO
	   ENDDO
   ENDDO

  DO ic=1,mesh(3)
       	  DO  ia=1,mesh(1)
	           DO ib=1,mesh(2)
		     Integrated_field2(ia,ic) =  Integrated_field2(ia,ic) +  target_field(ia,ib,ic) / mesh(2)
		   ENDDO
	  ENDDO
   ENDDO
!  !
   ! Average the Difference on x and y lines
   !

!  !
   ! Print the Difference in a Gnuplot file
   !
   ALLOCATE ( vect1(2) )
   vect1(:) = ZERO
   ALLOCATE ( vect2(2) )
   vect2(:) = ZERO
   vect1(1) = cell(2,2)
   vect1(2) = cell(2,3)
   vect2(1) = cell(3,2)
   vect2(2) = cell(3,3)

   CALL gnuplotwrite(iufile, FileOutputGnuplotx, 2,mesh(2),mesh(3),vect1, vect2, Integrated_field1(:,:))

   vect1(:) = ZERO
   vect2(:) = ZERO
   vect1(1) = cell(1,1)
   vect1(2) = cell(1,3)
   vect2(1) = cell(3,1)
   vect2(2) = cell(3,3)

!  !
   ! Print the Difference in a Gnuplot file
   !
   CALL gnuplotwrite(iufile, FileOutputGnuploty, 2,mesh(1),mesh(3),vect1, vect2,  Integrated_field2(:,:))


!  !
   ! Integrate on both x and y lines
   !

   ALLOCATE( Integrated_field_xy(mesh(3)) )
   Integrated_field_xy(:) = ZERO 
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_field_xy(ic) =  Integrated_field_xy(ic) +  target_field(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_field_xy(ic) = Integrated_field_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO

!  !
   ! Print the Difference in a .dat file
   !
   CALL datwrite(iufile, FileOutputDat, mesh(3), cell(3,3), Integrated_field_xy(:))

   Integrated_field_xy(:) = ZERO 
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_field_xy(ic) =  Integrated_field_xy(ic) +   field_0bias(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_field_xy(ic) = Integrated_field_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO


!  !
   ! Print the 0bias file averaged in a .dat file
   !
   CALL datwrite(iufile, FileOutputDat0bias, mesh(3), cell(3,3), Integrated_field_xy(:))

   Integrated_field_xy(:) = ZERO 
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_field_xy(ic) =  Integrated_field_xy(ic) +   field_Finitebias(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_field_xy(ic) = Integrated_field_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO

!  !
   ! Print the Finite Bias averaged in a .dat file
   !
   CALL datwrite(iufile, FileOutputDatFinitebias, mesh(3), cell(3,3), Integrated_field_xy(:))


!  !
   ! Load Wavefunction File in the .cube format
   !

   CALL  GaussianFilereadHeader(FileInputWF, iufile, Cutoffdistance, dimensions, Nb_atoms_check, mesh_check)

   IF (Nb_atoms /=  Nb_atoms_check) THEN 
       PRINT*, "Discrepancies in the atom numbers" 
       PRINT*, Nb_atoms_check, Nb_atoms
      STOP
   ENDIF
   IF ( ( ABS(mesh_check(1)- mesh(1)) >  EPS_m5 ).OR.( ABS(mesh_check(2) - mesh(2)) >  EPS_m5 ).OR.( ABS(mesh_check(3)-mesh(3))>  EPS_m5 ) ) THEN
       PRINT*, "Discrepancies in the meshes" 
       PRINT*, mesh_check(1), mesh(1)
       PRINT*, mesh_check(2), mesh(2) 
       PRINT*, mesh_check(3), mesh(3) 
      STOP
   ENDIF

   ALLOCATE( Wavefunction(mesh(1),mesh(2),mesh(3)) )   
   Wavefunction(:,:,:) = ZERO 
   ALLOCATE( atoms_type_check( Nb_atoms) )   
   ALLOCATE(atoms_position_check( Nb_atoms, dimensions) )   

   CALL  GaussianFileread(FileInputWF, iufile, Cutoffdistance, dimensions, Nb_atoms_check,mesh_check, X0Cell_check, cell_check, atoms_type_check, atoms_position_check, Wavefunction)

   DO ia =1, dimensions
        IF ( ( ABS(mesh(ia)*cell_check(ia,1)- cell(ia,1)) >  EPS_m5 ).OR.( ABS(mesh(ia)*cell_check(ia,2)- cell(ia,2)) >  EPS_m5 ).OR.( ABS(mesh(ia)*cell_check(ia,3)- cell(ia,3))>  EPS_m5 ) ) THEN
	        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...Test vector     : ",  mesh(ia)*cell_check(ia,1), " | " ,  mesh(ia)*cell_check(ia,2), " | " ,  mesh(ia)*cell_check(ia,3)
	        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...Reference vector: ",  cell(ia,1), " | " , cell(ia,2), " | " ,  cell(ia,3)
      		STOP
        ELSE IF ( ( ABS( X0Cell_check(ia)-X0Cell(ia) ) >  EPS_m5 )  ) THEN
                PRINT*, "Discrepancies in the origins" 
      		STOP
        ENDIF 
   ENDDO

   Check_Volume = ZERO
   dV = ZERO

   CALL CalculateVolume(Check_Volume, dV, dimensions, cell(:,:), mesh(:))



   CALL  CalculateFirstOrderEnergy(FirstOrderCorrection, Integrated_volume, Integrated_density, dimensions,  X0Cell, dV, cell, mesh,Wavefunction, target_field)
   PRINT*, "Integrated Volume [Ang^3]   : ", Integrated_volume 
   PRINT*, " Box Volume [Ang^3]         : ", Check_Volume
   PRINT*, "Integrated Density [e]      : ", Integrated_density 
   PRINT*, "First Order Correction [eV] : ", FirstOrderCorrection


   CALL CalculateCenterOfCharge(CenterOfCharge, Nb_electrons, dimensions, Nb_atoms, atoms_type, atoms_position)
   CALL CalculateDipoleOperatorWF(DipoleOperator(:), Integrated_volume, Integrated_density, dimensions, X0Cell(:), dV, cell(:,:), mesh(:),Wavefunction(:,:,:) )
   PRINT*, "Wavefunction Dipole Operator [eA-Coc]     :  "
   PRINT*,  (-DipoleOperator(idimension)+CenterOfCharge(idimension), idimension=1,dimensions)

!  !
   ! Deallocate
   !
   PRINT*, "Deallocate"
   DEALLOCATE(field_0bias, field_Finitebias, Wavefunction)
   DEALLOCATE(atoms_type, atoms_position,atoms_type_check, atoms_position_check)
   DEALLOCATE(target_field, Integrated_field1,Integrated_field2, Integrated_field_xy)
   DEALLOCATE ( vect1,vect2 )



   PRINT*, "END"
END PROGRAM WaveFunctionInPotentialDifference

