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
   PROGRAM DipoleFromChargeDensity
   !***********************************************
   ! Conventions: All the distances in this file are in ANGSTROMS
   !              All the energies  in this file are in eV
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   USE calculatequantities_Module,         ONLY : CalculateVolume, CalculateCenterOfMass, CalculateCenterOfCharge, CalculateDipoleOperator
   USE IO_module,                          ONLY : gnuplotwrite, datwrite, ReadInputFileDipole, GaussianFilewrite, GaussianFileread, RhoFilereadHeader,  RhoFileread,  XYZFilereadHeader, XYZFileread
   IMPLICIT NONE
   INTEGER, PARAMETER :: dimensions = 3    
   INTEGER, PARAMETER :: iufile=100
   REAL(dbl), PARAMETER ::  Cutoffdistance=20000000*ONE
   REAL(dbl), PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   CHARACTER*100, PARAMETER :: datafile="DATA.in"


   !
   ! Input variables
   !
   CHARACTER*100 ::  FileInput ! For the moment, just one density.RHO
   CHARACTER*100 ::  FileInputxyz ! For the ionic charges
   CHARACTER*100 ::  FileOutputGaussian ! For plotting the input density in a .cube file
   CHARACTER*100 ::  FileOutputGnuplotx ! For plotting the average of the input density in a .gp file in the yz plan
   CHARACTER*100 ::  FileOutputGnuploty ! For plotting the average of the input density in a .gp file in the xz plan
   CHARACTER*100 ::  FileOutputDat      ! For plotting the averaged  input density in .dat file along the z axis
   REAL(dbl) ::   cell(dimensions,dimensions), norm
   REAL(dbl) ::   X0Cell(dimensions)
   INTEGER ::   mesh(dimensions)
   INTEGER ::   mesh_test(dimensions)
   INTEGER ::   Nb_atoms 
   REAL(dbl) ::   Nb_electrons ! Valence electrons only

   REAL(dbl), ALLOCATABLE ::   density(:,:,:)
   INTEGER, ALLOCATABLE ::   atoms_type(:)
   REAL(dbl), ALLOCATABLE ::   atoms_position(:,:)
   !
   ! Output variables
   !
! Output
   REAL(dbl)              :: DipoleOperator(dimensions)  ! position operator

   !REAL(dbl)              :: PositionXYZ(dimensions)  ! position in the cell
   REAL(dbl)              :: CenterOfMass(dimensions)  ! barycenter
   REAL(dbl)              :: CenterOfCharge(dimensions)  ! center of the charges
   REAL(dbl)              :: Integrated_volume, dV
   REAL(dbl)              :: Check_Volume
   REAL(dbl), ALLOCATABLE :: Integrated_density1(:,:)
   REAL(dbl), ALLOCATABLE :: Integrated_density2(:,:)
   REAL(dbl), ALLOCATABLE :: Integrated_density_xy(:) 
   REAL(dbl)              :: Integrated_density  ! Check for the number of electrons
   ! 
   ! Local variables
   ! 
   INTEGER :: ia, ib, ic, idimension
   REAL(dbl), ALLOCATABLE :: vect1(:), vect2(:)
!***********************************************
!
!------------------------------
! main body
!------------------------------
!  !
   ! read input General Information
   !
   CALL ReadInputFileDipole(iufile, datafile, FileInput, FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat)

!  !
   ! read input file xyz
   !

   CALL XYZFilereadHeader(FileInputxyz, iufile, Nb_atoms)

   ALLOCATE( atoms_type( Nb_atoms) )   
   ALLOCATE(atoms_position( Nb_atoms, dimensions) )   

   CALL XYZFileread( FileInputxyz, iufile, dimensions, Nb_atoms, atoms_type, atoms_position)

   CALL CalculateCenterOfMass(CenterOfMass, dimensions, Nb_atoms, atoms_type, atoms_position)

   CALL CalculateCenterOfCharge(CenterOfCharge, Nb_electrons, dimensions, Nb_atoms, atoms_type, atoms_position)

!  !
   ! read input file 0 bias
   !
   CALL RhoFilereadHeader( FileInput, iufile, Cutoffdistance, dimensions, mesh)

   ALLOCATE( density(mesh(1),mesh(2),mesh(3)) )   
   density(:,:,:) = ZERO 
   CALL RhoFileread( FileInput, iufile,Cutoffdistance, dimensions,  mesh, X0Cell, cell, density, norm)
   density(:,:,:) =  Nb_electrons*density(:,:,:)/norm
   ! Output the re-scaled density in a Gaussian Cube File.
   CALL GaussianFilewrite(FileOutputGaussian, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, density)

!  !
   ! Average   on x and y lines
   !
   ALLOCATE( Integrated_density1(mesh(2),mesh(3)) )
   ALLOCATE( Integrated_density2(mesh(1),mesh(3)) )
   Integrated_density1(:,:) = ZERO
   Integrated_density2(:,:) = ZERO
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_density1(ib,ic) =  Integrated_density1(ib,ic) +  density(ia,ib,ic) / mesh(1)
		   ENDDO
	   ENDDO
   ENDDO

  DO ic=1,mesh(3)
       	  DO  ia=1,mesh(1)
	           DO ib=1,mesh(2)
		     Integrated_density2(ia,ic) =  Integrated_density2(ia,ic) +  density(ia,ib,ic) / mesh(2)
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

   CALL gnuplotwrite(iufile, FileOutputGnuplotx, 2,mesh(2),mesh(3),vect1, vect2, Integrated_density1(:,:))

   vect1(:) = ZERO
   vect2(:) = ZERO
   vect1(1) = cell(1,1)
   vect1(2) = cell(1,3)
   vect2(1) = cell(3,1)
   vect2(2) = cell(3,3)

!  !
   ! Print the Difference in a Gnuplot file
   !
   CALL gnuplotwrite(iufile, FileOutputGnuploty, 2,mesh(1),mesh(3),vect1, vect2,  Integrated_density2(:,:))


!  !
   ! Integrate on both x and y lines
   !

   ALLOCATE( Integrated_density_xy(mesh(3)) )
   Integrated_density_xy(:) = ZERO 
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_density_xy(ic) =  Integrated_density_xy(ic) + density(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_density_xy(ic) = Integrated_density_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO

!  !
   ! Print the Difference in a .dat file
   !
   CALL datwrite(iufile, FileOutputDat, mesh(3), cell(3,3), Integrated_density_xy(:))

!  !
   ! Calculate Dipole Moment through Position operator \int\int\int dxdydz r rho(x,y,z)
   !

   !PositionXYZ(:) = ZERO

   Check_Volume = ZERO
   dV = ZERO

   CALL CalculateVolume(Check_Volume, dV, dimensions, cell(:,:), mesh(:))

   CALL CalculateDipoleOperator(DipoleOperator(:), Integrated_volume, Integrated_density, dimensions, X0Cell(:), dV, cell(:,:), mesh(:),density(:,:,:) )

   PRINT*, "Integrated Volume [Ang^3]   : ", Integrated_volume 
   PRINT*, " Box Volume [Ang^3]         : ", Check_Volume
   PRINT*, "Integrated Density [e]      : ", Integrated_density 
   PRINT*, "Number of Electrons [Ang^3] : ", Nb_electrons
   PRINT*, "Ionic Dipole Operator [eA-Origin] :  "
   PRINT*,  (Nb_electrons*CenterOfCharge(idimension), idimension=1,dimensions)
   PRINT*, "Electonic Dipole Operator [eA-Origin]     :  "
   PRINT*,  (-DipoleOperator(idimension), idimension=1,dimensions)
   PRINT*, "Total Dipole Operator [eA]      :  "
   PRINT*,  (((Nb_electrons*CenterOfCharge(idimension))-DipoleOperator(idimension)) , idimension=1,dimensions)



!  !
   ! Deallocate
   !
   PRINT*, "Deallocate"
   DEALLOCATE(density)
   DEALLOCATE(atoms_type, atoms_position)
   DEALLOCATE(Integrated_density1,Integrated_density2, Integrated_density_xy)
   DEALLOCATE ( vect1,vect2 )



   PRINT*, "END"
END PROGRAM DipoleFromChargeDensity

