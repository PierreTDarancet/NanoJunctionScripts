!
! Copyright (C) 2012 The Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Sahar Sharifzadeh, Pierre Darancet
!***********************************************
!
! This program calculates the electron-hole
! 2-particle correlation function from the 
! output of Berkeley GW.
!

   PROGRAM Calculate2PartCorr
   !***********************************************
    ! USES SUBROUTINES:
    !invert3x3
    !_
    !GaussianFileread
    !ReadInput
    !gnuplotwrite
    !GaussianFilereadHeader1
    !GaussianFilereadHeader2
    !GaussianFilewrite
    !GaussianFileread
  
   !***********************************************
   IMPLICIT NONE
   !***********************************************
   !** Parameters
   Integer, parameter :: maxchar_in = 55
   Integer, parameter :: maxchar_out = 100
   INTEGER, PARAMETER :: dimensions = 3    
   INTEGER, PARAMETER :: iufile=100
   !REAL*8, PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   CHARACTER*100, PARAMETER :: datafile="DATA.in"
   CHARACTER*18, PARAMETER :: subroutinename="Calculate2PartCorr"

! ...   Numerical constants
        
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::     TWO = 2.0d0
        REAL*8, PARAMETER ::   THREE = 3.0d0
        REAL*8, PARAMETER ::    FOUR = 4.0d0
        COMPLEX*16, PARAMETER:: CZERO = (0.0d0, 0.0d0)
        COMPLEX*16, PARAMETER::  CONE = (1.0d0, 0.0d0)
        COMPLEX*16, PARAMETER::    CI = (0.0d0, 1.0d0)

        REAL*8, PARAMETER ::      PI = 3.14159265358979323846d0
        REAL*8, PARAMETER ::     TPI = 2.0d0 * 3.14159265358979323846d0
        REAL*8, PARAMETER ::     FPI = 4.0d0 * 3.14159265358979323846d0
        REAL*8, PARAMETER ::   SQRT2 = 1.41421356237309504880d0
        REAL*8, PARAMETER ::   SQRT3 = 1.73205080756887729353d0
        REAL*8, PARAMETER ::  SQRTPI = 1.77245385090551602729d0
        REAL*8, PARAMETER :: SQRTPM1 = 1.0d0 / 1.77245385090551602729d0

        REAL*8, PARAMETER ::      EPS_m1  = 0.1d0
        REAL*8, PARAMETER ::      EPS_m2  = 0.01d0
        REAL*8, PARAMETER ::      EPS_m3  = 0.001d0
        REAL*8, PARAMETER ::      EPS_m4  = 0.0001d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0
        REAL*8, PARAMETER ::      EPS_m6  = 0.000001d0
        REAL*8, PARAMETER ::      EPS_m7  = 0.0000001d0
        REAL*8, PARAMETER ::      EPS_m8  = 0.00000001d0
        REAL*8, PARAMETER ::      EPS_m9  = 0.000000001d0
        REAL*8, PARAMETER ::      EPS_m10 = 0.0000000001d0
        REAL*8, PARAMETER ::      EPS_m11 = 0.00000000001d0
        REAL*8, PARAMETER ::      EPS_m12 = 0.000000000001d0
        REAL*8, PARAMETER ::      EPS_m13 = 0.0000000000001d0
        REAL*8, PARAMETER ::      EPS_m14 = 0.00000000000001d0

! ...   Physical constants
        REAL*8, PARAMETER :: K_BOLTZMAN_SI    = 1.38066D-23       ! J K^-1 
        REAL*8, PARAMETER :: K_BOLTZMAN_AU    = 3.1667D-6         ! Hartree K^-1 
        REAL*8, PARAMETER :: K_BOLTZMAN_M1_AU = 315795.26D0       ! Hartree^-1 K 
        REAL*8, PARAMETER :: FACTEM           = 315795.26D0       ! Hartree^-1 K 

! ...   Physical constants defining the Atomic Units System
        REAL*8, PARAMETER :: BOHR_RADIUS_SI   = 0.529177D-10      ! m
        REAL*8, PARAMETER :: BOHR_RADIUS_CM   = 0.529177D-8       ! cm
        REAL*8, PARAMETER :: BOHR_RADIUS_ANGS = 0.529177D0        ! angstrom
        REAL*8, PARAMETER :: ELECTRONMASS_SI  = 9.10953D-31       ! Kg
        REAL*8, PARAMETER :: ELECTRONMASS_UMA = 5.4858D-4         ! uma

! ...   Units conversion factors
        REAL*8, PARAMETER :: ELECTRONVOLT_SI  = 1.6021892D-19     ! J  
        REAL*8, PARAMETER :: UMA_SI           = 1.66057D-27       ! Kg
        REAL*8, PARAMETER :: ANGSTROM_AU      = 1.889727D0        ! au
        REAL*8, PARAMETER :: AU_TO_OHMCMM1    = 46000.0D0         ! (ohm cm)^-1
        REAL*8, PARAMETER :: AU_KB            = 294210.0D0        ! Kbar
        REAL*8, PARAMETER :: KB_AU            = 1.0D0/294210.0D0  ! au
        REAL*8, PARAMETER :: AU               = 27.211652d0       ! eV
        REAL*8, PARAMETER :: RYD              = 13.605826d0       ! eV 
        REAL*8, PARAMETER :: SCMASS           = 1822.89D0         ! uma to au
        REAL*8, PARAMETER :: UMA_AU           = 1822.89D0         ! au
        REAL*8, PARAMETER :: AU_TERAHERTZ     = 2.418D-5          ! THz
        REAL*8, PARAMETER :: TERAHERTZ        = 2.418D-5          ! from au to THz
        REAL*8, PARAMETER :: AU_SEC           = 2.4189D-17        ! sec
       

        REAL*8, PARAMETER :: rhothr = 1.0d-5 ! tolerance
        REAL*8, PARAMETER :: gsmall = 1.0d-12

        REAL*8, PARAMETER :: e2 = 2.d0           ! the square of the electron charge
        REAL*8, PARAMETER :: degspin = 2.d0      ! the number of spins per level
        REAL*8, PARAMETER :: rytoev=13.6058d0    ! conversion from Ry to eV

        !  mass conversion: a.m.u to a.u. (Ry)
        REAL*8, PARAMETER :: amconv= 1.66042d-24/9.1095d-28*0.5d0 
        !  pressure conversion from Ry/(a.u)^3 to K
        REAL*8, PARAMETER :: uakbar= 147105.d0



   !***********************************************
   !** Output Variables
   REAL*8, ALLOCATABLE ::   twopartcorr(:,:,:)
   REAL*8, ALLOCATABLE ::   twopartcorr_distance(:)
   REAL*8, ALLOCATABLE ::   HolePositions(:,:)
   REAL*8, ALLOCATABLE ::   cdf(:,:,:)
   REAL*8, ALLOCATABLE ::   cdf_distance(:)
   REAL*8, ALLOCATABLE ::   cdf_distance_renormalized(:)
   INTEGER, ALLOCATABLE ::  counter(:,:,:)
   INTEGER, ALLOCATABLE ::  counter_distance(:)

   !** Local Variables
   ! Calculated Quantities
   !   
   INTEGER ::   Nb_Distances

   REAL*8  ::   DistanceStep
   REAL*8 ::   ElectronPosition(dimensions)
   REAL*8 ::   AverageDistance
   REAL*8 ::   Dipole(dimensions)
   REAL*8 ::   Quadrupole(dimensions,dimensions)

   REAL*8 ::   diagonal
   REAL*8 ::   distance
   REAL*8 ::   VolCel, dV
   REAL*8 ::   latticetoxyztransfermatrix(dimensions,dimensions)
   REAL*8 ::    xyztolatticetransfermatrix(dimensions,dimensions)
   REAL*8,ALLOCATABLE ::   HolePositionLattice(:,:)
   INTEGER, ALLOCATABLE ::   HolePositioniLattice(:,:)
   ! Input Quantities
   !   
   REAL*8, ALLOCATABLE :: norm(:)
   REAL*8 ::   cell(dimensions,dimensions), sphereradius
   REAL*8 ::   X0Cell(dimensions), reducedgridforplottingX0Cell(dimensions)
   INTEGER ::   grid(dimensions)
   INTEGER ::   grid_aux(dimensions)
   INTEGER ::   reducedgridforplotting(dimensions)
   INTEGER ::   holegrid(dimensions)
   INTEGER ::   grid_test(dimensions)
   INTEGER ::   Nb_atoms, Nb_atoms_aux
   REAL*8, ALLOCATABLE  ::   ElectronicDensity(:,:,:), AvgElectronicDensity(:,:,:)
   REAL*8, ALLOCATABLE  ::   atoms_position(:,:)
   INTEGER, ALLOCATABLE ::  atoms_type(:), atoms_type_aux(:)
   CHARACTER*(maxchar_in) ::fname=' ', GenericNameInputFile=' '
   CHARACTER*(maxchar_out) ::GenericNameOutputFile=' ', outputfname=' ', outputtype=' '

   ! Loop variables
   ! 
   INTEGER :: ia, ib, ic, id, id2, iah, ibh, ich, ia_aux, ib_aux, ic_aux, ih, iatom

   ! Auxiliary variables
   ! 
   REAL*8, ALLOCATABLE  :: vect1(:), vect2(:), TWODFunction_aux1(:,:), TWODFunction_aux2(:,:), TWODFunction_aux3(:,:),ONEDFunction_aux1(:), ONEDFunction_aux2(:), ONEDFunction_aux3(:), atoms_position_aux1(:,:),  atoms_position_aux2(:,:)
   REAL*8 :: Maximum_norm, vect_norm(dimensions), mat3x3_aux(dimensions,dimensions), holeposition_aux(dimensions), holeposition_aux2(dimensions), vect3_aux(dimensions), vect3_aux2(dimensions), real_aux, real_aux2
!***********************************************
!
!------------------------------
! main body
!------------------------------



   !***********************************************
   !
   !------------------------------
   ! Input
   !------------------------------
   PRINT*, "***************************************************"
   PRINT*, " Copyright (C) 2012 The Molecular Foundry Berkeley"
   PRINT*, " This file is distributed under the terms of the"
   PRINT*, " GNU General Public License. See the file `License'"
   PRINT*, " in the root directory of the present distribution,"
   PRINT*, " or http://www.gnu.org/copyleft/gpl.txt ."
   PRINT*, " Contributors  : Sahar Sharifzadeh, Pierre Darancet"
   PRINT*, "**************************************************"
   PRINT*, " This program calculates the electron-hole"
   PRINT*, " 2-particle correlation function from the"
   PRINT*, " output of Berkeley GW."
   PRINT*, "**************************************************"
   PRINT*, " This general input file is: ",  TRIM(datafile)
   PRINT*, "**************************************************"
   PRINT*, " Beginning..."
!  !
   ! read input General Information
   !
   PRINT*, " Reading generic input file..."
   CALL ReadInputFile(iufile, datafile, GenericNameInputFile, holegrid(1), holegrid(2), holegrid(3),GenericNameOutputFile )
   PRINT*, " ...done"
   !------------------------------
   ! Input
   !------------------------------
   ! Done
   !***********************************************

   !***********************************************
   !
   !------------------------------
   ! Allocations
   !------------------------------
   !  HOLE GRID
   PRINT*, "Allocating Hole Grid..."
   PRINT*, "... Number of Input Files:", (holegrid(1)*holegrid(2)*holegrid(3))
   PRINT*, "... The electronic densities will be read from files:"
   DO iah=1, holegrid(1)
      DO ibh=1, holegrid(2)
         DO ich=1, holegrid(3)
            write(fname, '(A50,3i1)') TRIM(GenericNameInputFile), (iah-1), (ibh-1), (ich-1)
            PRINT*,  TRIM(fname)
         ENDDO
      ENDDO
   ENDDO
   PRINT*, " ...done"
   !  Electronic GRID
   PRINT*, " Reading first electronic density file for allocation..."
   iah=1
   ibh=1
   ich=1
   fname=' '
   WRITE(fname, '(A55)')  "                                                                                                    "
   WRITE(fname, '(A50,3I1)') TRIM(GenericNameInputFile), (iah-1), (ibh-1), (ich-1)
   PRINT*, " File...", TRIM(fname)
   CALL  GaussianFilereadHeader1(trim(fname), iufile, dimensions, Nb_atoms, grid, X0Cell, cell, VolCel, dV)
   PRINT*, " ...done"
   !  Atomic GRID
   PRINT*, " Reading first electronic density file for atoms position..." 
   ALLOCATE(atoms_type(Nb_atoms))
   ALLOCATE(atoms_position(Nb_atoms,dimensions))
   CALL  GaussianFilereadHeader2(trim(fname), iufile, dimensions, Nb_atoms, grid, X0Cell, cell,atoms_type,atoms_position)
   !  Distance GRID
   PRINT*, " Calculating the distance grid..."
   diagonal=sqrt((abs(grid(1)*cell(1,1))  +  abs(grid(2)*cell(2,1)) +  abs(grid(3)*cell(3,1)))**2 +(abs(grid(1)*cell(1,2))  +  abs(grid(2)*cell(2,2)) +  abs(grid(3)*cell(3,2)))**2 +(abs(grid(1)*cell(1,3))  +  abs(grid(2)*cell(2,3)) +  abs(grid(3)*cell(3,3)))**2)
   Maximum_norm=ZERO
   vect_norm(:)=ZERO
   DO id=1,dimensions
   	DO id2=1,dimensions
   		vect_norm(id)=vect_norm(id) + ((cell(id,id2))**2 )
   	ENDDO
        vect_norm(id)=sqrt(vect_norm(id))
        IF (vect_norm(id)>Maximum_norm ) Maximum_norm=vect_norm(id)
   ENDDO
   DO id=1,(grid(1)*grid(2)*grid(3))
	   DistanceStep = diagonal/REAL(id)
          IF (DistanceStep > Maximum_norm) Nb_Distances=id
   ENDDO
   DistanceStep = diagonal/Nb_Distances
   PRINT*, "Allocating Variables..."
   ALLOCATE(HolePositions((holegrid(1)*holegrid(2)*holegrid(3)),dimensions))
   ALLOCATE(norm( (holegrid(1)*holegrid(2)*holegrid(3)) ))
   ALLOCATE( HolePositionLattice((holegrid(1)*holegrid(2)*holegrid(3)),dimensions))
   ALLOCATE(HolePositioniLattice((holegrid(1)*holegrid(2)*holegrid(3)),dimensions))
   ALLOCATE(twopartcorr(((2*grid(1))-1),((2*grid(2))-1),((2*grid(3))-1)))
   ALLOCATE(cdf(((2*grid(1))-1),((2*grid(2))-1),((2*grid(3))-1)))
   ALLOCATE(counter(((2*grid(1))-1),((2*grid(2))-1),((2*grid(3))-1)))
   ALLOCATE(ElectronicDensity(grid(1),grid(2),grid(3)))
   ALLOCATE(AvgElectronicDensity(grid(1),grid(2),grid(3)))
   ALLOCATE(twopartcorr_distance(Nb_Distances))
   ALLOCATE(cdf_distance(Nb_Distances))
   ALLOCATE(cdf_distance_renormalized(Nb_Distances))
   ALLOCATE(counter_distance(Nb_Distances))
   PRINT*, " ...done"
   !  Printing Summary
   PRINT*, "... Largest distance in the supercell [Bohr]:", diagonal
   PRINT'(A35,3F12.4)', "... Norm of Lattice Vector [Bohr]: ", (vect_norm(id), id=1,dimensions)
   PRINT*, "... Maximal Norm of Lattice Vector [Bohr]:  ", Maximum_norm
   PRINT*, "... Step in distance [Bohr]:  ", DistanceStep
   PRINT*, "... Number of steps        :  ", Nb_Distances
   PRINT*, " ...done"
   !------------------------------
   ! Allocations
   !------------------------------
   ! Done
   !***********************************************


   !***********************************************
   !
   !------------------------------
   ! Initializing Variables
   !------------------------------
   PRINT*, "Initializing Variables..."
   HolePositions(:,:)=ZERO
   HolePositionLattice(:,:)=ZERO
   twopartcorr(:,:,:)=ZERO
   twopartcorr_distance(:)=ZERO
   cdf(:,:,:)=ZERO
   cdf_distance(:)=ZERO
   counter(:,:,:)=0
   counter_distance(:)=0
   HolePositioniLattice(:,:)=0
   AvgElectronicDensity(:,:,:)=ZERO
   PRINT*, " ...done"
   ! Done
   !***********************************************


   !***********************************************
   !
   !------------------------------
   ! Initializing Conversion lattice-cartesian
   !------------------------------
   ! Calculation 
   PRINT*, "Calculating Transfer Matrix..."
   latticetoxyztransfermatrix(:,:)=ZERO 
   xyztolatticetransfermatrix(:,:)=ZERO
   mat3x3_aux(:,:)=ZERO
   latticetoxyztransfermatrix(:,:)= TRANSPOSE( cell(:,:))
   CALL invert3x3(dimensions,latticetoxyztransfermatrix(:,:),xyztolatticetransfermatrix(:,:))
   !  Test
   mat3x3_aux(:,:)= matmul(latticetoxyztransfermatrix(:,:),xyztolatticetransfermatrix(:,:)) 
   real_aux=ABS(mat3x3_aux(2,1))+ABS(mat3x3_aux(3,1))+ABS(mat3x3_aux(1,2))+ABS(mat3x3_aux(3,2))+ABS(mat3x3_aux(1,3))+ABS(mat3x3_aux(2,3))
   IF ( real_aux > EPS_m5 ) THEN
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Problem in the product of transfer matrices"
            PRINT*, "Expected   Wing Sum: ", ZERO
            PRINT*, "Calculated Wing Sum: ", real_aux
   ENDIF
   real_aux=ABS(mat3x3_aux(1,1))+ABS(mat3x3_aux(2,2))+ABS(mat3x3_aux(3,3))
   IF (  ( real_aux > (ONE*3.+ EPS_m5) ) .OR.  ( real_aux < (ONE*3.- EPS_m5) )) THEN
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Problem in the product of transfer matrices"
            PRINT*, "Expected   Trace: ", ONE*3.
            PRINT*, "Calculated Trace: ", real_aux
   ENDIF
   !  Printing Summary
   PRINT*, "Transfer Matrix.. Lattice to XYZ:"
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",latticetoxyztransfermatrix(1,1)  , " | " ,latticetoxyztransfermatrix(1,2)   , " | " ,latticetoxyztransfermatrix(1,3)  , " ) "
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",latticetoxyztransfermatrix(2,1)  , " | " ,latticetoxyztransfermatrix(2,2)   , " | " ,latticetoxyztransfermatrix(2,3)  , " ) "
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",latticetoxyztransfermatrix(3,1)  , " | " ,latticetoxyztransfermatrix(3,2)   , " | " ,latticetoxyztransfermatrix(3,3)  , " ) "  
   PRINT*, "Transfer Matrix.. XYZ to Lattice:"
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",xyztolatticetransfermatrix(1,1)  , " | " ,xyztolatticetransfermatrix(1,2)   , " | " ,xyztolatticetransfermatrix(1,3)  , " ) "
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",xyztolatticetransfermatrix(2,1)  , " | " ,xyztolatticetransfermatrix(2,2)   , " | " ,xyztolatticetransfermatrix(2,3)  , " ) "
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",xyztolatticetransfermatrix(3,1)  , " | " ,xyztolatticetransfermatrix(3,2)   , " | " ,xyztolatticetransfermatrix(3,3)  , " ) "  
   PRINT*, "Transfer Matrix.. Test Unity:"
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",mat3x3_aux(1,1)  , " | " ,mat3x3_aux(1,2)   , " | " ,mat3x3_aux(1,3)  , " ) "
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",mat3x3_aux(2,1)  , " | " ,mat3x3_aux(2,2)   , " | " ,mat3x3_aux(2,3)  , " ) "
   PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",mat3x3_aux(3,1)  , " | " ,mat3x3_aux(3,2)   , " | " ,mat3x3_aux(3,3)  , " ) "  
   ! Done
   !***********************************************


!***********************************************
!
!------------------------------
! Loop over hole position
!------------------------------
   ih=0
   PRINT*, "Loop over hole positions..."
   PRINT*, " ... Grid",  holegrid(1), holegrid(2), holegrid(3)
   ! Go through all hole locations 
   do iah=1, holegrid(1)
     do ibh=1, holegrid(2)
        do ich=1, holegrid(3)
            ih=ih+1
            PRINT*, "...Position...", ih, " of", (holegrid(1)*holegrid(2)*holegrid(3))
            !***********************************************
            ! Reading Electronic Density
            fname=' '
            WRITE(fname, '(A55)')  "                                                                                                    "
            WRITE(fname, '(A50,3i1)') TRIM(GenericNameInputFile), (iah-1), (ibh-1), (ich-1)
            PRINT*, "...Reading from file:", TRIM(fname)
            CALL GaussianFileread(trim(fname), iufile, dimensions, Nb_atoms,  grid, X0Cell, cell, atoms_type, atoms_position,   dV, VolCel, HolePositions(ih,:),  ElectronicDensity(:,:,:),norm(ih))
            ! Done
            !***********************************************
            AvgElectronicDensity(:,:,:) =  AvgElectronicDensity(:,:,:) +  ElectronicDensity(:,:,:)
            !***********************************************

            !***********************************************
            ! Calculating Hole Position
            vect3_aux(:)=(HolePositions(ih,:)-X0Cell(:))
            HolePositionLattice(ih,:)=ZERO
            HolePositionLattice(ih,:)=matmul(xyztolatticetransfermatrix(:,:),vect3_aux(:)) 
            holeposition_aux(:)= matmul(latticetoxyztransfermatrix(:,:),HolePositionLattice(ih,:))  + X0Cell(:)
            DO id=1, dimensions
                HolePositioniLattice(ih,id)=NINT( HolePositionLattice(ih,id) )
                vect3_aux(id)=REAL( HolePositioniLattice(ih,id) )
            ENDDO
            holeposition_aux2(:)=matmul(latticetoxyztransfermatrix(:,:),vect3_aux(:))  + X0Cell(:)
            ! Done
            !***********************************************

            !***********************************************
            ! Testing Hole Position
            IF ( SUM(HolePositions(ih,:)-holeposition_aux(:)) > EPS_m5 )  THEN
 	           PRINT*, "WARNING: Routine", TRIM(subroutinename)
 	           PRINT*, "WARNING: Problem in hole position"
 	           PRINT*, "Read position    [Bohr]   : ", HolePositions(ih,:)
 	           PRINT*, "Calculated position [Bohr]: ", holeposition_aux(:)
                   STOP
            ENDIF
            DO id=1, dimensions              
            	IF ( ABS( HolePositions(ih,id)- holeposition_aux2(id) )  > Maximum_norm )  THEN
 	           PRINT*, "WARNING: Routine", TRIM(subroutinename)
 	           PRINT*, "WARNING: Problem in hole position"
 	           PRINT'(A31,3F12.4)', "Read position    [Bohr]   :  ", HolePositions(ih,:)
 	           PRINT'(A31,3F12.4)', "Discretized position [Bohr]: ", holeposition_aux2(:)
      	           PRINT'(A31,F12.4)',  "Max Norm unit vector [Bohr]: ", Maximum_norm
                   STOP
                ENDIF
            ENDDO
            ! Done
            !***********************************************

            !***********************************************
            ! Printing Hole Position
            PRINT*, "...Expressing the hole coordinates in lattice unit"
            PRINT'(A38,3F12.4)', "...Hole coordinates [Bohr]           :", HolePositions(ih,:)
            PRINT'(A38,3F12.4)', "...Hole coordinates wrt Origin [Bohr]:", (HolePositions(ih,:)-X0Cell(:))
            PRINT'(A38,3F12.4)', "...Hole coordinates in lattice vector:", HolePositionLattice(ih,:)
            PRINT'(A38,3F12.4)', "...Error in the bases change:         ", (HolePositions(ih,:)-holeposition_aux(:))
            PRINT'(A38,3I7)',    "...Hole in the electron lattice      :", HolePositioniLattice(ih,:)
            PRINT'(A38,3F12.4)', "...Error in the discretization[Bohr] :", (HolePositions(ih,:)-holeposition_aux2(:))
            ! Done
            !***********************************************

 
	    !***********************************************
	    !
	    !------------------------------
	    ! Loop over electron position
	    !------------------------------
            PRINT*, "... Beginning of the electron position loop"
            DO ia=1,grid(1)
              DO ib=1,grid(2)
                 DO ic=1,grid(3)
		    !***********************************************
		    !
		    !------------------------------
		    ! Distance between electron and hole
		    !------------------------------
                    ! Calculating electron position in xyz coordinates
                    ElectronPosition(:)=ZERO
                    vect3_aux(1)=(real(ia)-0.5)
                    vect3_aux(2)=(real(ib)-0.5)
                    vect3_aux(3)=(real(ic)-0.5)
                    ElectronPosition(:)= matmul(latticetoxyztransfermatrix(:,:),vect3_aux(:))  + X0Cell(:)
                    distance=ZERO   
                    DO id=1,dimensions
                         distance=distance+ (( HolePositions(ih,id)- ElectronPosition(id) )**2)
                    ENDDO
                    distance=sqrt(distance)
                    ! Distance on the lattice 
                    ia_aux= ia-HolePositioniLattice(ih,1) + grid(1) ! ia-HolePositioniLattice(ih,1) runs from -grid(1)+1  to grid(1)-1
                    ib_aux= ib-HolePositioniLattice(ih,2) + grid(2) ! ia/b/c_aux from 1 to 2grid(1/2/3)-1 are the indices for cdf and 2 part corr
                    ic_aux= ic-HolePositioniLattice(ih,3) + grid(3) ! They indicate the distance vector and can represent negative distances
                    IF ( (ia_aux < 1).OR.(ia_aux > ((2*grid(1))-1)).OR.(ib_aux < 1).OR.(ib_aux > ((2*grid(2))-1)).OR.(ic_aux < 1).OR.(ic_aux > ((2*grid(3))-1))    )  THEN
                          PRINT*, "WARNING: value of the distance vector indices (The calculation will stop)"
                          PRINT*, "WARNING: Indices",ia_aux,"of", ((2*grid(1))-1),ib_aux,"of",((2*grid(2))-1),ic_aux ,"of",((2*grid(3))-1)
                          PRINT*, "WARNING: Elecron Grid", ia,"of",grid(1),  ib,"of",grid(2), ic,"of",grid(3) 
                          PRINT*, "WARNING: Hole Grid", HolePositioniLattice(ih,1), HolePositioniLattice(ih,2), HolePositioniLattice(ih,3)
                          STOP
                    ENDIF
		    !
		    !------------------------------
		    !  Distance between electron and hole
		    !------------------------------
	            ! Done
	            !***********************************************


		    !***********************************************
		    !
		    !------------------------------
		    ! Calculation of the 2 part corr and cdf
		    !------------------------------
                    ! Distance
                    DO id=1,Nb_Distances
                    	real_aux=(id-1)*DistanceStep                                                                        ! Find the position on the distance grid
                    	IF (distance.lt.real_aux+DistanceStep) then 
                       		cdf_distance(id)=cdf_distance(id)+ ElectronicDensity(ia,ib,ic)                              ! Store the component of the charge density in the cdf
                        	IF(distance.ge.real_aux) then
                          		twopartcorr_distance(id) = twopartcorr_distance(id) +  ElectronicDensity(ia,ib,ic)  ! Store the component of the charge density in the twoparticle correlation function
                         		counter_distance(id)=counter_distance(id)+1   
                       		ENDIF
                       ENDIF
                    ENDDO ! Nb_distances
                    ! Vector
                    twopartcorr(ia_aux,ib_aux,ic_aux)= twopartcorr(ia_aux,ib_aux,ic_aux)+ElectronicDensity(ia,ib,ic)
                    counter(ia_aux,ib_aux,ic_aux)= counter(ia_aux,ib_aux,ic_aux)+1                   
	            ! Done
	            !***********************************************
            	 ENDDO  ! ic
              ENDDO     ! ib
            ENDDO        ! ia
	    !
	    !------------------------------
	    ! Loop over electron position
	    !------------------------------
            ! Done
            !***********************************************
             PRINT*, "... End of the electron position loop"
    	   ENDDO              ! ich
   	   ENDDO                 ! ibh
	ENDDO                    ! iah
        PRINT*, "... End of the hole position loop"
        !
	!------------------------------
	! Loop over hole position
	!------------------------------
        ! Done
        !***********************************************





	!***********************************************
	!
	!------------------------------
	! Renormalization 
	!------------------------------
        ! Avg Density
        AvgElectronicDensity(:,:,:) =  AvgElectronicDensity(:,:,:)/(holegrid(1)*holegrid(2)*holegrid(3))
        ! Two Part Corr (d)
        PRINT*, "... Renormalizing twopartcorr_distance"
        DO id=1,Nb_Distances 
              ! Radius of the sphere
              distance= (id-1)*DistanceStep
              ! Difference in the Volumes of the sphere of radii distance and distance+DistanceStep 
              !UNCOMMENT NEXT LINES IF WANTING TO CORRECT FOR THE VOLUME OF THE SHELL
              !real_aux= ((distance)**3)*FOUR*PI/THREE
              !real_aux2= ((distance+DistanceStep)**3)*FOUR*PI/THREE 
              ! if the grid is fine enough, the difference in volume is dV * counter
   	      IF (counter_distance(id).eq.0) then
              	twopartcorr_distance(id)=ZERO
       	      ELSE
              	!twopartcorr_distance(id)=twopartcorr_distance(id)
        	twopartcorr_distance(id)=twopartcorr_distance(id)/(REAL(counter_distance(id)))
        	!twopartcorr_distance(id)=twopartcorr_distance(id) * (real_aux2-real_aux) / (dV*counter_distance(id))
       	      ENDIF
       	      cdf_distance_renormalized(id)=SUM(twopartcorr_distance(1:id))
        ENDDO !Nb_distances
        twopartcorr_distance(:)=twopartcorr_distance(:)/cdf_distance_renormalized(Nb_Distances)
        PRINT*, "... Done"







        PRINT*, "... Calculating reduced grid"
        !Finding the radius of the sphere containing 90 % of the electron
        DO id=1,Nb_Distances 
              IF(cdf_distance_renormalized(Nb_Distances+1-id) > 0.9) THEN
                  sphereradius = DistanceStep*(REAL(Nb_Distances+1-id))
                  reducedgridforplotting(1) =   NINT(2.0*sphereradius / vect_norm(1))
                  reducedgridforplotting(2) =   NINT(2.0*sphereradius / vect_norm(2))
                  reducedgridforplotting(3) =   NINT(2.0*sphereradius / vect_norm(3))
                  reducedgridforplottingstarting(1) = NINT(grid(1)-0.5*REAL(reducedgridforplotting(1)))
                  vect_aux(:)= REAL(reducedgridforplotting(:) +1 - (2*grid(:) ) )
		  reducedgridforplottingX0Cell(:) = matmul(latticetoxyztransfermatrix(:,:),  vect_aux(:))  + X0Cell(:)
                  vect_aux(:)= vect_aux(:) + REAL(reducedgridforplotting(:))
              ENDIF
        ENDDO !Nb_distances
        PRINT*, "... Done"
        PRINT*, "Reduced Grid starting at "

            PRINT'(A38,3F12.4)', "...Reduced Grid Origin  [Bohr]       :", reducedgridforplottingX0Cell(:)
            PRINT'(A38,3F12.4)', "...Reduced Grid Maximum [Bohr]       :", reducedgridforplottingX0Cell(:)

            PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",mat3x3_aux(1,1)  , " | " ,mat3x3_aux(1,2)   , " | " ,mat3x3_aux(1,3)  , " ) "
            PRINT'(A38,3F12.4)', "...Hole coordinates wrt Origin [Bohr]:", (HolePositions(ih,:)-X0Cell(:))
            PRINT'(A38,3F12.4)', "...Hole coordinates in lattice vector:", HolePositionLattice(ih,:)


		!***********************************************
		!
		!------------------------------
		! Output 1D data
		!------------------------------
		!PRINT*, "... Output 1-d files:"
		PRINT*, "...Writing 1D Correlation Function"
		outputtype="1D-distance"
		! WRITE(outputfname, '(A55)')  "                                                                                                    "
		WRITE(outputfname, '(A50,A1,A11,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".dat"
		PRINT*, "...Writing to file:", TRIM(outputfname)
		OPEN(iufile,file=trim(outputfname),form='formatted')
		REWIND(iufile)
		WRITE(iufile,*) "# Distance [Bohr], 2-part corr, cdf, NbPoints, cdf_from_num"
		DO id=1,Nb_Distances 
		      ! Radius of the sphere
		      distance= (id-1)*DistanceStep
		      WRITE(iufile,'(4ES13.5,I6)')  distance,  twopartcorr_distance(id), cdf_distance_renormalized(id), cdf_distance(id), counter_distance(id)
		ENDDO
		CLOSE(iufile)
	        PRINT*, "... Done"
	        ! Done
	        !***********************************************
        ! Two Part Corr (a,b,c) - recalculate with counter 0 corrections
        PRINT*, "... Renormalizing twopartcorr_abc"
        DO ic=1,((2*grid(3))-1)
		DO ib=1,((2*grid(2))-1)
			DO ia=1,((2*grid(1))-1)
				IF(counter(ia,ib,ic).eq.0) THEN
					twopartcorr(ia,ib,ic)=ZERO
				ELSE
					twopartcorr(ia,ib,ic)=twopartcorr(ia,ib,ic)/REAL(counter(ia,ib,ic))
				ENDIF
			ENDDO
		ENDDO
	ENDDO
        PRINT*, "... Done"
        !PRINT*, "... Recalculating cdf"
        !do ic=1,grid(3)
        !   Do ib=1,grid(2)!
        !      do ia=1,grid(1)
        !         cdf(ia,ib,ic)=SUM(twopartcorr(1:ia,1:ib,1:ic))
        !      end do  ! ic
        !   end do     ! ib!
	!end do        ! ia
        !PRINT*, "... Done"
        ! Two Part Corr (a,b,c) - Normalize here, we make the choice to NOT calculate the full cdf for efficiency purpose
	twopartcorr(:,:,:)=twopartcorr(:,:,:)/cdf_distance_renormalized(Nb_Distances)
	!------------------------------
	! Renormalization
	!------------------------------
        ! Done
        !***********************************************



	!***********************************************
	!
	!------------------------------
	! Calculation of Complementary Functions 
	!------------------------------
        ! Allocation
        ALLOCATE(TWODFunction_aux1(((2*grid(2))-1),((2*grid(3))-1)))
	ALLOCATE(TWODFunction_aux2(((2*grid(1))-1),((2*grid(3))-1)))
        ALLOCATE(TWODFunction_aux3(((2*grid(1))-1),((2*grid(2))-1)))
        ALLOCATE(ONEDFunction_aux1(((2*grid(1))-1)))
        ALLOCATE(ONEDFunction_aux2(((2*grid(2))-1)))
        ALLOCATE(ONEDFunction_aux3(((2*grid(3))-1)))
        ALLOCATE(vect1(2))
        ALLOCATE(vect2(2))
        Nb_atoms_aux = Nb_atoms  + (holegrid(1)*holegrid(2)*holegrid(3))
        ALLOCATE(atoms_position_aux1(Nb_atoms_aux,dimensions))
        ALLOCATE(atoms_position_aux2(Nb_atoms_aux,dimensions))
        ALLOCATE(atoms_type_aux(Nb_atoms))
        ! Initialize
	TWODFunction_aux1(:,:) = ZERO
	TWODFunction_aux2(:,:) = ZERO
	TWODFunction_aux3(:,:) = ZERO
	ONEDFunction_aux1(:) = ZERO
	ONEDFunction_aux2(:) = ZERO
	ONEDFunction_aux3(:) = ZERO
	!***********************************************
        ! Distance
        PRINT*, "... Calculating Average distance"
        AverageDistance=ZERO
        DO id=1,Nb_Distances 
              ! Radius of the sphere
                distance= (id-1)*DistanceStep
        	AverageDistance=AverageDistance+ ( twopartcorr_distance(id)*DistanceStep )
        	!twopartcorr_distance(id)=twopartcorr_distance(id) * (real_aux2-real_aux) / (dV*counter_distance(id))
        ENDDO !Nb_distances
        PRINT*, "... Done"
        PRINT*, "Avg Distance: ", AverageDistance, "[Bohr]"
        ! Done
	!***********************************************

	!***********************************************
        ! Dipole
	Dipole(:)=ZERO
        Quadrupole(:,:)=ZERO
        ! Quadrupole
        PRINT*, "... Calculating Dipole and quadrupole"
        DO ic=1,((2*grid(3))-1)
		DO ib=1,((2*grid(2))-1)
			DO ia=1,((2*grid(1))-1)
                            vect3_aux(:)=ZERO ! Vector r_h - r_e
		            vect3_aux2(1)=(real(ia-grid(1))) ! Vector r_h - r_e in lattice coordinates (running from -grid()+1 to grid()+1)
		            vect3_aux2(2)=(real(ib-grid(2)))
		            vect3_aux2(3)=(real(ic-grid(3)))
		            vect3_aux(:)= matmul(latticetoxyztransfermatrix(:,:),vect3_aux2(:)) ! Vector r_h - r_e
        		    Dipole(:)=Dipole(:)+ (vect3_aux(:)*twopartcorr(ia,ib,ic))
                            distance=sqrt((vect3_aux(1)**2) +(vect3_aux(2)**2)+(vect3_aux(3)**2))
                            Quadrupole(1,1)=Quadrupole(1,1)+(twopartcorr(ia,ib,ic)*( (3*(vect3_aux(1)**2)) - (distance**2) ))
                            Quadrupole(2,2)=Quadrupole(2,2)+(twopartcorr(ia,ib,ic)*( (3*(vect3_aux(2)**2)) - (distance**2) ))
                            Quadrupole(3,3)=Quadrupole(3,3)+(twopartcorr(ia,ib,ic)*( (3*(vect3_aux(3)**2)) - (distance**2) ))

                            Quadrupole(2,1)=Quadrupole(2,1)+(twopartcorr(ia,ib,ic)*(3*vect3_aux(2)*vect3_aux(1)))
                            Quadrupole(3,1)=Quadrupole(3,1)+(twopartcorr(ia,ib,ic)*(3*vect3_aux(3)*vect3_aux(1)))
                            Quadrupole(3,2)=Quadrupole(3,2)+(twopartcorr(ia,ib,ic)*(3*vect3_aux(3)*vect3_aux(2)))


			ENDDO
		ENDDO
       ENDDO
       Quadrupole(1,2)=Quadrupole(2,1)
       Quadrupole(1,3)=Quadrupole(3,1)
       Quadrupole(1,3)=Quadrupole(3,2)
       PRINT*, "... Done"
       PRINT*, "Dipole [eBohr]:"
       PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)',  " ( ",Dipole(1)  , " | " ,Dipole(2)   , " | " ,Dipole(3)  , " ) "
       PRINT*, "Quadrupole [eBohr^2]:"
       PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",Quadrupole(1,1)  , " | " ,Quadrupole(1,2)   , " | " ,Quadrupole(1,3)  , " ) "
       PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",Quadrupole(2,1)  , " | " ,Quadrupole(2,2)   , " | " ,Quadrupole(2,3)  , " ) "
       PRINT'(A3,F12.6,A3,F12.6,A3,F12.6,A3)', " ( ",Quadrupole(3,1)  , " | " ,Quadrupole(3,2)   , " | " ,Quadrupole(3,3)  , " ) "  
       ! Done
       !***********************************************


        ! Average 2D
		! Average over a
		DO ic=1,((2*grid(3))-1)
			   DO ib=1,((2*grid(2))-1)
				   DO  ia=1,((2*grid(1))-1)
				     TWODFunction_aux1(ib,ic) =  TWODFunction_aux1(ib,ic) + (twopartcorr(ia,ib,ic) / grid(1))
				   ENDDO
			   ENDDO
		ENDDO
		! Average over b
		DO ic=1,((2*grid(3))-1)
			   DO ia=1,((2*grid(1))-1)
				   DO  ib=1,((2*grid(2))-1)
				     TWODFunction_aux2(ia,ic) =  TWODFunction_aux2(ia,ic) +  (twopartcorr(ia,ib,ic) / grid(2))
				   ENDDO
			   ENDDO
		ENDDO
		! Average over c
		DO ia=1,((2*grid(1))-1)
			   DO  ib=1,((2*grid(2))-1)
				   DO ic=1,((2*grid(3))-1)
				     TWODFunction_aux3(ia,ib) =  TWODFunction_aux3(ia,ib) +  (twopartcorr(ia,ib,ic) / grid(3))
				   ENDDO
			   ENDDO
		ENDDO
        ! Average 1D
		! Average over a AND b
		DO ic=1,((2*grid(3))-1)
			   DO ib=1,((2*grid(2))-1)
				     ONEDFunction_aux3(ic) =  ONEDFunction_aux3(ic) + (TWODFunction_aux1(ib,ic)  / grid(2))
			   ENDDO
		ENDDO
		! Average over b AND c
		DO ia=1,((2*grid(1))-1)
		 	   DO ic=1,((2*grid(3))-1)
				     ONEDFunction_aux1(ia) = ONEDFunction_aux1(ia)+ (TWODFunction_aux2(ia,ic) / grid(3))
			   ENDDO
		ENDDO
		! Average over a AND c
		DO ib=1,((2*grid(2))-1)
			   DO ic=1,((2*grid(3))-1)
				     ONEDFunction_aux2(ib) =  ONEDFunction_aux2(ib) + (TWODFunction_aux1(ib,ic)  / grid(3))
			   ENDDO
		ENDDO

        ! Holes Positions
		atoms_position_aux1(:,:)=ZERO
		atoms_position_aux2(:,:)=ZERO
		atoms_type_aux(:)=0
		DO iatom=1, Nb_atoms
		       atoms_position_aux1(iatom,:)=atoms_position(iatom,:)
		       atoms_type_aux(iatom)=atoms_type(iatom)
		       vect3_aux(:)=(atoms_position_aux1(iatom,:)-X0Cell(:))
		       atoms_position_aux2(iatom,:)=matmul(xyztolatticetransfermatrix(:,:),vect3_aux(:))
		       atoms_position_aux2(iatom,:)=100.0*atoms_position_aux2(iatom,:)/grid(:)
		ENDDO
		DO ih=1,(holegrid(1)*holegrid(2)*holegrid(3))
		       atoms_position_aux1(Nb_atoms+ih,:)=HolePositions(ih,:)
		       atoms_type_aux(Nb_atoms+ih)=0
		       vect3_aux(:)=(atoms_position_aux1(Nb_atoms+ih,:)-X0Cell(:))
		       atoms_position_aux2(Nb_atoms+ih,:)=matmul(xyztolatticetransfermatrix(:,:),vect3_aux(:))
		       atoms_position_aux2(Nb_atoms+ih,:)=100.0*atoms_position_aux2(Nb_atoms+ih,:)/grid(:)
		ENDDO
		mat3x3_aux(:,:)=ZERO
		mat3x3_aux(1,1)=100.0*ONE/grid(1)
		mat3x3_aux(2,2)=100.0*ONE/grid(2)
		mat3x3_aux(3,3)=100.0*ONE/grid(3)
		vect3_aux(:)=ZERO

!***********************************************
!***********************************************
!        PRINT*, "... Output xyz files:"
!***********************************************



		   PRINT*, "... Output xyz file:"
		   outputtype="Coord"
		   WRITE(outputfname, '(A55)')  "                                                                                                    "
		   WRITE(outputfname, '(A50,A1,A5,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".xyz"
          	   CALL XYZFilewrite(outputfname, iufile, dimensions, (Nb_atoms+(holegrid(1)*holegrid(2)*holegrid(3))), atoms_type_aux, atoms_position_aux1)
		   ! Done
		   !***********************************************
 
	!------------------------------
	! Calculation of Complementary Functions
	!------------------------------
        ! Done
        !***********************************************

	!***********************************************
	!
	!------------------------------
	! Print Summary File
	!------------------------------
        ! Coordinates, Norms, distance, dipole, quadrupole
	PRINT*, "... Output Summary file:"
	outputtype="OUT"
	WRITE(outputfname, '(A55)')  "                                                                                                    "
	WRITE(outputfname, '(A50,A1,A3,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".txt"
	OPEN(iufile,file=trim(outputfname),form='formatted')
	REWIND(iufile)
	WRITE(iufile,*) "***************************************************"
	WRITE(iufile,*) " Copyright (C) 2012 The Molecular Foundry Berkeley"
	WRITE(iufile,*) " This file is distributed under the terms of the"
	WRITE(iufile,*) " GNU General Public License. See the file `License'"
	WRITE(iufile,*) " in the root directory of the present distribution,"
	WRITE(iufile,*) " or http://www.gnu.org/copyleft/gpl.txt ."
	WRITE(iufile,*) " Contributors  : Sahar Sharifzadeh, Pierre Darancet"
	WRITE(iufile,*) "**************************************************"
	WRITE(iufile,*) " This program calculates the electron-hole"
	WRITE(iufile,*) " 2-particle correlation function from the"
	WRITE(iufile,*) " output of Berkeley GW."	
	WRITE(iufile,*) "**************************************************"
	WRITE(iufile,*) " This general input file is: ",  TRIM(datafile)
	WRITE(iufile,*) "**************************************************"
	WRITE(iufile,*) " This general output file is: ",  TRIM(GenericNameOutputFile)
	WRITE(iufile,*) "Number of Hole positions:", (holegrid(1)*holegrid(2)*holegrid(3))
	WRITE(iufile,*) "# Coordinates [Bohr], Norm"
	DO ih=1,(holegrid(1)*holegrid(2)*holegrid(3))
	      WRITE(iufile,'(5ES13.5)')  HolePositions(ih,:), norm(ih)
	ENDDO
	WRITE(iufile,*) "Total Norm:", SUM(norm(1:(holegrid(1)*holegrid(2)*holegrid(3))))
       WRITE(iufile,*) "Dipole [eBohr]:"
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)')  " ( ",Dipole(1)  , " | " ,Dipole(2)   , " | " ,Dipole(3)  , " ) "
       PRINT*, "Quadrupole [eBohr^2]:"
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)') " ( ",Quadrupole(1,1)  , " | " ,Quadrupole(1,2)   , " | " ,Quadrupole(1,3)  , " ) "
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)') " ( ",Quadrupole(2,1)  , " | " ,Quadrupole(2,2)   , " | " ,Quadrupole(2,3)  , " ) "
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)') " ( ",Quadrupole(3,1)  , " | " ,Quadrupole(3,2)   , " | " ,Quadrupole(3,3)  , " ) "  
       CLOSE(iufile)

	!------------------------------
	! Print Summary File
	!------------------------------
        ! Done
        !***********************************************



	!***********************************************
	!
	!------------------------------
	! Print Complementary Functions 
	!------------------------------
        ! 
        ! Norms
		!
		!------------------------------
		! OutputNorm data
		!------------------------------
		PRINT*, "... Output Norm file:"
		outputtype="Norm"
		! WRITE(outputfname, '(A55)')  "                                                                                                    "
		WRITE(outputfname, '(A50,A1,A4,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".dat"
		PRINT*, "...Writing to file:", TRIM(outputfname)
		OPEN(iufile,file=trim(outputfname),form='formatted')
		REWIND(iufile)
		WRITE(iufile,*) "# Coordinates [Bohr], Norm"
		DO ih=1,(holegrid(1)*holegrid(2)*holegrid(3))
		      WRITE(iufile,'(5ES13.5)')  HolePositions(ih,:), norm(ih)
		ENDDO
		CLOSE(iufile)
	        PRINT*, "... Done"
        ! 1D

        PRINT*, "...Writing 1D avg Correlation Function"
        outputtype="1D-c"
       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A4,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".dat"
        PRINT*, "...Writing to file:", TRIM(outputfname)
        OPEN(iufile,file=trim(outputfname),form='formatted')
        REWIND(iufile)
        WRITE(iufile,*) "# Position on c, averaged 2part-corr"
        DO id=1,grid(3) 
              distance= REAL((id-1.0))/REAL(grid(3))
              WRITE(iufile,'(2ES13.5)')  distance,  ONEDFunction_aux3(id) 
        ENDDO
        CLOSE(iufile)
	!***********************************************
        outputtype="1D-b"
       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A4,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".dat"
        PRINT*, "...Writing to file:", TRIM(outputfname)
        OPEN(iufile,file=trim(outputfname),form='formatted')
        REWIND(iufile)
        WRITE(iufile,*) "# Position on b, averaged 2part-corr"
        DO id=1,grid(2) 
              distance= REAL((id-1.0))/REAL(grid(2))
              WRITE(iufile,'(2ES13.5)')  distance,  ONEDFunction_aux2(id) 
        ENDDO
        CLOSE(iufile)
	!***********************************************
        outputtype="1D-a"
       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A4,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".dat"
        PRINT*, "...Writing to file:", TRIM(outputfname)
        OPEN(iufile,file=trim(outputfname),form='formatted')
        REWIND(iufile)
        WRITE(iufile,*) "# Position on a, averaged 2part-corr"
        DO id=1,grid(1) 
              distance= REAL((id-1.0))/REAL(grid(1))
              WRITE(iufile,'(2ES13.5)')  distance,  ONEDFunction_aux1(id) 
        ENDDO
        CLOSE(iufile)
	!***********************************************
        PRINT*, "... Done"

        PRINT*, "... Output 2-d files:"
 
!        PRINT*, "......Writing 2D Correlation averaged over z"
!        outputtype="2D-xyc"
!        WRITE(outputfname, '(A50,A1,A5,A3)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".gp"        
!        PRINT*, "...Writing to file:", TRIM(outputfname)

!        PRINT*, "...... Done"
!        PRINT*, "......Writing 2D Correlation averaged over y"
!        outputtype="2D-xz"
!        WRITE(outputfname, '(A50,A1,A5,A3)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".gp"        
!        PRINT*, "...Writing to file:", TRIM(outputfname)

!        PRINT*, "...... Done"
!        PRINT*, "......Writing 2D Correlation averaged over x"
!        outputtype="2D-yz"
!        WRITE(outputfname, '(A50,A1,A5,A3)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".gp"        
!        PRINT*, "...Writing to file:", TRIM(outputfname)
!        PRINT*, "...... Done"
!***********************************************
        PRINT*, "......Writing 2D Correlation averaged over c"
        outputtype="2D-ab"
       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A5,A3)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".gp"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        vect1(1) = ONE /REAL(grid(1))
        vect1(2) = ZERO/REAL(grid(1))
        vect2(1) = ZERO/REAL(grid(2))
        vect2(2) = ONE /REAL(grid(2))
        CALL gnuplotwrite(iufile, TRIM(outputfname), 2, grid(1), grid(2), vect1(:), vect2(:), TWODFunction_aux3(:,:))
        PRINT*, "...... Done"
!***********************************************
        PRINT*, "......Writing 2D Correlation averaged over b"
        outputtype="2D-ac"
       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A5,A3)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".gp"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        vect1(1) = ONE /REAL(grid(1))
        vect1(2) = ZERO/REAL(grid(1))
        vect2(1) = ZERO/REAL(grid(3))
        vect2(2) = ONE /REAL(grid(3))
        CALL gnuplotwrite(iufile, TRIM(outputfname), 2, grid(1), grid(3), vect1(:), vect2(:), TWODFunction_aux2(:,:))
        PRINT*, "...... Done"
!***********************************************
        PRINT*, "......Writing 2D Correlation averaged over a"
        outputtype="2D-bc"
       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A5,A3)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".gp"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        vect1(1) = ONE /REAL(grid(2))
        vect1(2) = ZERO/REAL(grid(2))
        vect2(1) = ZERO/REAL(grid(3))
        vect2(2) = ONE /REAL(grid(3))
        CALL gnuplotwrite(iufile, TRIM(outputfname), 2, grid(2), grid(3), vect1(:), vect2(:), TWODFunction_aux1(:,:))
        PRINT*, "...... Done"
!***********************************************
        PRINT*, "... Done"
!***********************************************
 
!***********************************************
        PRINT*, "......Writing 3D Wavefunction on xyz"
        outputtype="3D_DEN"
        WRITE(outputfname, *) 
        !WRITE(outputfname, '(A100)')  "                                                                                                    "
        WRITE(outputfname,  '(A50,A1,A6,A5)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".cube"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        CALL GaussianFilewrite(TRIM(outputfname), iufile, dimensions, Nb_atoms_aux, grid, X0Cell, cell,atoms_type_aux,atoms_position_aux1, AvgElectronicDensity)
        PRINT*, "...... Done"

!***********************************************
! PRINT
!***********************************************

        grid_aux(1)=-1+(2*grid(1))
        grid_aux(2)=-1+(2*grid(2))
        grid_aux(3)=-1+(2*grid(3))


!***********************************************
        PRINT*, "......Writing 3D Correlation on abc"
        outputtype="3D-abc"
        !WRITE(outputfname, '(A100)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A6,A5)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".cube"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        CALL GaussianFilewrite(TRIM(outputfname), iufile, dimensions, Nb_atoms_aux, grid_aux, vect3_aux, mat3x3_aux,atoms_type_aux,atoms_position_aux2, twopartcorr)
        PRINT*, "...... Done"

!***********************************************


        PRINT*, "......Writing 3D Correlation on xyz"
        outputtype="3D-xyz"
        WRITE(outputfname, *) 
        !WRITE(outputfname, '(A100)')  "                                                                                                    "
        WRITE(outputfname,  '(A50,A1,A6,A5)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".cube"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        CALL GaussianFilewrite(TRIM(outputfname), iufile, dimensions, Nb_atoms_aux, reducedgridforplotting, reducedgridforplottingX0Cell, cell,atoms_type_aux,atoms_position_aux1, twopartcorr(reducedgridforplotting:,:,:))
        PRINT*, "...... Done"


!***********************************************
       PRINT*, "... Output 3-d files:"
!***********************************************
reducedgridforplotting
        PRINT*, "......Writing 3D Correlation on xyz"
        outputtype="3D-gen"
       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
        WRITE(outputfname, '(A50,A1,A6,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".txt"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        OPEN(iufile,file=trim(outputfname),form='formatted')
           do ia=1,reducedgridforplotting(1)
              do ib=1,reducedgridforplotting(2)
                 do ic=1,reducedgridforplotting(3)
		    ! distance between electron and hole
                    ElectronPosition(:)=ZERO
                    vect3_aux(1)=(real(ia))
                    vect3_aux(2)=(real(ib))
                    vect3_aux(3)=(real(ic))
                    ElectronPosition(:)= matmul(latticetoxyztransfermatrix(:,:),vect3_aux(:))  + reducedgridforplottingX0Cell(:)

                    WRITE(iufile,'(4ES14.7,I6)')  ElectronPosition(:), twopartcorr(ia+grid(1),ib+grid(2),ic+grid(3)) counter(ia+grid(1),ib+grid(2),ic+grid(3))
          	 ENDDO
       		ENDDO
	ENDDO
        CLOSE(iufile)


!***********************************************
!        PRINT*, "...... Done"
!***********************************************

!***********************************************
!        PRINT*, "...... Done"
!***********************************************
!        PRINT*, "... Done"


   PRINT*, "Deallocating Variables..."
   DEALLOCATE(HolePositions)
   DEALLOCATE( HolePositionLattice)
   DEALLOCATE(HolePositioniLattice)
   DEALLOCATE(norm)
   DEALLOCATE(twopartcorr)
   DEALLOCATE(twopartcorr_distance)
   DEALLOCATE(cdf)
   DEALLOCATE(cdf_distance)
   DEALLOCATE(cdf_distance_renormalized)
   DEALLOCATE(counter)
   DEALLOCATE(counter_distance)
   DEALLOCATE(ElectronicDensity)
   DEALLOCATE(AvgElectronicDensity)
   DEALLOCATE(TWODFunction_aux1)
   DEALLOCATE(TWODFunction_aux2)
   DEALLOCATE(TWODFunction_aux3)
   DEALLOCATE(ONEDFunction_aux1)
   DEALLOCATE(ONEDFunction_aux2)
   DEALLOCATE(ONEDFunction_aux3)
   DEALLOCATE(vect1)
   DEALLOCATE(vect2)
       DEALLOCATE(atoms_position_aux1)
       DEALLOCATE(atoms_position_aux2)
       DEALLOCATE(atoms_type_aux)


   PRINT*, " ...done"
   PRINT*, " End"
   END PROGRAM Calculate2PartCorr
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
      SUBROUTINE ReadInputFile(iufile, datafile,GenericNameInputFile,Nb_H_a,Nb_H_b,Nb_H_c,GenericNameOutputFile)
   !***********************************************
   ! For density calculations
   IMPLICIT NONE
! ...   Numerical constants
        
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0
   CHARACTER*9, PARAMETER :: subroutinename="ReadInput"
   Integer, parameter :: maxchar_in = 55
   Integer, parameter :: maxchar_out = 100
 

   INTEGER, INTENT(in)              :: iufile
   CHARACTER*12, INTENT(in)         :: datafile
   CHARACTER*(*), INTENT(out)        :: GenericNameInputFile
   INTEGER, INTENT(out)        :: Nb_H_a
   INTEGER, INTENT(out)        :: Nb_H_b
   INTEGER, INTENT(out)        :: Nb_H_c
   CHARACTER*(*), INTENT(out)        :: GenericNameOutputFile


   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) GenericNameInputFile
        PRINT*, "   ...READING: Generic Name for electronic density file: ", Trim(GenericNameInputFile)
        READ(iufile,*) Nb_H_a
        PRINT*, "   ...READING: Number of holes along a: ", Nb_H_a
        READ(iufile,*) Nb_H_b
        PRINT*, "   ...READING: Number of holes along b: ", Nb_H_b
        READ(iufile,*) Nb_H_c
        PRINT*, "   ...READING: Number of holes along c: ", Nb_H_c
        READ(iufile,*) GenericNameOutputFile
        PRINT*, "   ...READING: Generic Name for electronic density file: ", Trim(GenericNameOutputFile)

        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFile


      subroutine invert3x3(dimensions,A, Ainv)
IMPLICIT NONE
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0

   CHARACTER*9, PARAMETER :: subroutinename="Invert3x3"

      INTEGER, INTENT(in) ::dimensions
      REAL*8, INTENT(in)  :: A(dimensions,dimensions)              ! input
      REAL*8, INTENT(out) ::  Ainv(dimensions,dimensions)           ! output
      REAL*8 ::  Adet
      !
      IF (dimensions /= 3) THEN
         PRINT*, "Impossible to invert cell, dimensions /=3"
         STOP
      ENDIF

      Ainv(:,:) = 0d+0;
      Ainv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      Ainv(1,2) = A(3,2)*A(1,3) - A(1,2)*A(3,3)
      Ainv(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)

      Adet = Ainv(1,1)*A(1,1) + Ainv(1,2)*A(2,1) + Ainv(1,3)*A(3,1)
 
      IF (Adet==ZERO) THEN
         PRINT*, "Impossible to invert cell, det=0"
         STOP
      ENDIF
      Ainv(1,1) = Ainv(1,1)/Adet
      Ainv(1,2) = Ainv(1,2)/Adet
      Ainv(1,3) = Ainv(1,3)/Adet
      Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))/Adet
      Ainv(2,2) = (A(1,1)*A(3,3) - A(3,1)*A(1,3))/Adet
      Ainv(2,3) = (A(2,1)*A(1,3) - A(1,1)*A(2,3))/Adet
      Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))/Adet
      Ainv(3,2) = (A(3,1)*A(1,2) - A(1,1)*A(3,2))/Adet
      Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/Adet
      return     ! end of invert3x3 matrix code
      end subroutine invert3x3
!***********************************************
!***********************************************
!***********************************************
!***********************************************
      SUBROUTINE GaussianFilereadHeader1(FileInput, iufile, dimensions, Nb_atoms, grid, X0Cell, cell, VolCel, dV)
!   CALL  GaussianFilereadHeader1(fname, iufile, dimensions, Nb_atoms, grid, X0Cell, cell, VolCel, dV)

   !***********************************************
    IMPLICIT NONE
! ...   Numerical constants
        
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0

   CHARACTER*(*), INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile

   INTEGER,   INTENT(in)            :: dimensions
   INTEGER , INTENT(out)            :: grid(dimensions)
   INTEGER , INTENT(out)            :: Nb_atoms
   REAL*8,   INTENT(out)         :: X0Cell(dimensions)            ! position of the origin of the volumetric data. 
   REAL*8,   INTENT(out)         :: cell(dimensions,dimensions)  ! cell vectors
   REAL*8,   INTENT(out)         :: VolCel
   REAL*8,   INTENT(out)         :: dV
   CHARACTER*100  :: chr
  REAL*8         :: Crossaux(dimensions), aux(dimensions), real_io
        PRINT*, "...READING Incoming Cube File: Header", trim(FileInput)
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        READ(iufile,*) Nb_atoms, X0Cell(1),  X0Cell(2),  X0Cell(3) !    60  -40.000000  -30.000000  -30.000000
        PRINT'(A43,I5)', "    ...Number of atoms read              : ", Nb_atoms
        Nb_atoms=Nb_atoms-1
        PRINT'(A43,I5)', "    ...Number of actual atoms in the cell: ", Nb_atoms
        PRINT'(A25,3F12.6)', "    ...Origin[Bohr]   : ", X0Cell(:)
        READ(iufile,*) grid(1), cell(1,1), cell(1,2), cell(1,3)    !  300    0.267559    0.000000    0.000000
        READ(iufile,*) grid(2), cell(2,1), cell(2,2), cell(2,3)    !  300    0.000000    0.200669    0.000000
        READ(iufile,*) grid(3), cell(3,1), cell(3,2), cell(3,3)    !   50    0.000000    0.000000    1.632653
        CLOSE(iufile)
        PRINT*, "...Done"
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid    : ",  grid(1), " x " ,  grid(2), " x ",  grid(3)
        PRINT'(A24,F12.6,A3,F12.6,A3,F12.6)', "    ...A vector[Bohr]: ",  grid(1)*cell(1,1), " | " ,  grid(1)*cell(1,2), " | " ,  grid(1)*cell(1,3) 
        PRINT'(A24,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector[Bohr]: ",  grid(2)*cell(2,1), " | " ,  grid(2)*cell(2,2), " | " ,  grid(2)*cell(2,3) 
        PRINT'(A24,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector[Bohr]: ",  grid(3)*cell(3,1), " | " ,  grid(3)*cell(3,2), " | " ,  grid(3)*cell(3,3) 
        Crossaux(1) = (grid(1)*cell(1,2)*grid(2)*cell(2,3))  -(grid(2)*cell(2,2)* grid(1)*cell(1,3) )
        Crossaux(2) = (grid(1)*cell(1,3)*grid(2)*cell(2,1))  -(grid(2)*cell(2,3)*grid(1)*cell(1,1))
        Crossaux(3) = (grid(1)*cell(1,1)*grid(2)*cell(2,2))  - (grid(2)*cell(2,1)*grid(1)*cell(1,2))
        aux(1) = grid(3)*cell(3,1)
        aux(2) = grid(3)*cell(3,2)
        aux(3) = grid(3)*cell(3,3)
        VolCel = abs(dot_product(Crossaux(:), aux(:)))
        dV = VolCel / (grid(1)*grid(2)*grid(3))
        PRINT'(A36,F12.4)', "    ...Volume of the cell[Bohr^3]: ", VolCel
        PRINT'(A36,F12.6)', "    ...Voxel:            [Bohr^3] ", dV
END  SUBROUTINE GaussianFilereadHeader1

!***********************************************
      SUBROUTINE GaussianFilereadHeader2(FileInput, iufile, dimensions, Nb_atoms_check, grid_check, X0Cell_check, cell_check , atoms_type,atoms_position)
   !***********************************************
   IMPLICIT NONE

! ...   Numerical constants
        
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0
   CHARACTER*33, PARAMETER :: subroutinename="GaussianFilereadHeader2"

   CHARACTER*(*), INTENT(in)    :: FileInput
   INTEGER, INTENT(in)          :: iufile
   INTEGER,   INTENT(in)        :: dimensions
   INTEGER , INTENT(in)         :: grid_check(dimensions)
   INTEGER , INTENT(in)         :: Nb_atoms_check
   REAL*8,   INTENT(in)         :: X0Cell_check(dimensions)            ! position of the origin of the volumetric data. 
   REAL*8,   INTENT(in)         :: cell_check(dimensions,dimensions)  ! cell vectors
   !
   INTEGER,   INTENT(out)       ::  atoms_type(Nb_atoms_check)     !atoms_type( Nb_atoms )
   REAL*8,   INTENT(out)        ::  atoms_position(Nb_atoms_check,dimensions)


   REAL*8              :: barycenter(dimensions)
   REAL*8              :: CenterOfMass(dimensions)
   REAL*8              :: TotalMass

   INTEGER             :: grid_aux(dimensions)
   INTEGER             :: Nb_atoms_aux, idimension
   REAL*8              :: X0Cell_aux(dimensions)            ! position of the origin of the volumetric data. 
   REAL*8              :: cell_aux(dimensions,dimensions)  ! cell vectors
   INTEGER             :: iatom
  REAL*8         :: Crossaux(dimensions), aux(dimensions), real_io
   CHARACTER*100  :: chr
   Integer, parameter :: maxchar_in = 55
   Integer, parameter :: maxchar_out = 100
        PRINT*, "...READING Incoming Cube File: Header",trim(FileInput)
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        READ(iufile,*) Nb_atoms_aux, X0Cell_aux(1),  X0Cell_aux(2),  X0Cell_aux(3) !    60  -40.000000  -30.000000  -30.000000
        PRINT'(A25,3F12.6)', "    ...Origin [Bohr]   : ", X0Cell_aux(:)
        PRINT'(A43,I5)', "    ...Number of atoms read              : ", Nb_atoms_aux
        Nb_atoms_aux=Nb_atoms_aux-1
        PRINT'(A43,I5)', "    ...Number of actual atoms in the cell: ", Nb_atoms_aux
        READ(iufile,*) grid_aux(1), cell_aux(1,1), cell_aux(1,2), cell_aux(1,3)    !  300    0.267559    0.000000    0.000000
        READ(iufile,*) grid_aux(2), cell_aux(2,1), cell_aux(2,2), cell_aux(2,3)    !  300    0.000000    0.200669    0.000000
        READ(iufile,*) grid_aux(3), cell_aux(3,1), cell_aux(3,2), cell_aux(3,3)    !   50    0.000000    0.000000    1.632653
        PRINT*, "...Done"
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid    : ",  grid_aux(1), " x " ,  grid_aux(2), " x ",  grid_aux(3)
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid Ref: ",  grid_check(1), " x " ,  grid_check(2), " x ",  grid_check(3)
 

        IF ( Nb_atoms_aux /= Nb_atoms_check ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the number of actual atoms in the cell"
            PRINT*, "Reference: ", Nb_atoms_check
            PRINT*, "Read     : ", Nb_atoms_aux
            STOP
        ENDIF
        IF ( (ABS(X0Cell_aux(1)-X0Cell_check(1)) +  ABS(X0Cell_aux(2)-X0Cell_check(2))+ ABS(X0Cell_aux(3)-X0Cell_check(3))) > EPS_m5 ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the origin of the cell"
            PRINT*, "Reference: ", X0Cell_check(:)
            PRINT*, "Read     : ", X0Cell_aux(:)
            STOP
        ENDIF
       DO idimension = 1,dimensions
           IF (grid_check(idimension) /= grid_aux(idimension)) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the grid"
            PRINT*, "Reference: ", grid_check(:)
            PRINT*, "Read     : ", grid_aux(:)
            STOP
           ENDIF
           IF ( (ABS(cell_aux(idimension,1)- cell_check(idimension,1) ) + ABS(cell_aux(idimension,2)- cell_check(idimension,2))+ ABS(cell_aux(idimension,3)- cell_check(idimension,3))) > EPS_m5  ) THEN
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the lattice vector", idimension
            PRINT*, "Reference: ", cell_check(idimension,:)
            PRINT*, "Read     : ", cell_aux(idimension,:)
            STOP
           ENDIF
        ENDDO
        PRINT*, "    ...READING Atoms type and positions"
        barycenter(:)=ZERO
        CenterOfMass(:) = ZERO
        TotalMass= ZERO
    DO iatom = 1, Nb_atoms_aux
        atoms_type(iatom)=1111
        atoms_position(iatom,:)=10000*ONE 
        READ(iufile,*)   atoms_type(iatom), real_io, atoms_position(iatom,1), atoms_position(iatom,2), atoms_position(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146
        PRINT'(A18,I5,A6,F12.6,A1,F12.6,A1,F12.6,A2)', "        ...Atom: ",atoms_type(iatom)," at (",atoms_position(iatom,1),",",atoms_position(iatom,2),",",atoms_position(iatom,3)," )"
        barycenter(1)=        barycenter(1) + atoms_position(iatom,1)
        barycenter(2)=        barycenter(2) + atoms_position(iatom,2)
        barycenter(3)=        barycenter(3) + atoms_position(iatom,3)
        CenterOfMass(1) = CenterOfMass(1) + ( atoms_position(iatom,1)* atoms_type(iatom)*2 ) 
        CenterOfMass(2) = CenterOfMass(2) + ( atoms_position(iatom,2)* atoms_type(iatom)*2 ) 
        CenterOfMass(3) = CenterOfMass(3) + ( atoms_position(iatom,3)* atoms_type(iatom)*2 ) 
        TotalMass= TotalMass + atoms_type(iatom)*2
     ENDDO
        barycenter(:) = barycenter(:) / Nb_atoms_aux
        CenterOfMass(:) =CenterOfMass(:) / TotalMass
        PRINT'(A38,F12.6,A1,F12.6,A1,F12.6,A2)', "    ...Barycenter of the system at: (",barycenter(1),",",barycenter(2),",",barycenter(3), " )"
        PRINT'(A38,F12.6,A1,F12.6,A1,F12.6,A2)', "    ...Center of Mass[Bohr]       : (", CenterOfMass(1),",", CenterOfMass(2),",", CenterOfMass(3), " )"




END  SUBROUTINE GaussianFilereadHeader2
!***********************************************
      SUBROUTINE GaussianFileread(FileInput, iufile, dimensions, Nb_atoms_check, grid_check, X0Cell_check, cell_check, atoms_type_check, atoms_position_check, dV_check, VolCel_check, holeposition, wavefunction, norm)
!            CALL GaussianFileread(FileInput, iufile, dimensions, Nb_atoms,  grid, X0Cell, cell, atoms_type, atoms_position,  dV, VolCel, HolePositions(ih,:),  ElectronicDensity(:,:,:),norm)
   !***********************************************
    IMPLICIT NONE
! ...   Numerical constants
        
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0
        REAL*8, PARAMETER ::      EPS_m1  = 0.1d0
   CHARACTER*26, PARAMETER :: subroutinename="GaussianFileread"
   Integer, parameter :: maxchar_in = 55
   Integer, parameter :: maxchar_out = 100
   !
   CHARACTER*(*), INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: Nb_atoms_check
   INTEGER,   INTENT(in)            :: grid_check(dimensions)
   REAL*8,   INTENT(in)         :: X0Cell_check(dimensions)            ! position of the origin of the volumetric data. 
   REAL*8,   INTENT(in)         :: cell_check(dimensions,dimensions)  ! cell vectors
   INTEGER,   INTENT(in)        :: atoms_type_check(Nb_atoms_check)     !atoms_type( Nb_atoms )
   REAL*8,   INTENT(in)         ::  atoms_position_check(Nb_atoms_check,dimensions)
   REAL*8,   INTENT(in)         :: dV_check
   REAL*8,   INTENT(in)         :: VolCel_check

   REAL*8,   INTENT(out)         :: holeposition(dimensions)            
   REAL*8,   INTENT(out)         :: wavefunction(grid_check(1),grid_check(2),grid_check(3)) !wavefunction( grid(1),grid(2),grid(3) )
   REAL*8,   INTENT(out)         :: norm

! Local/Reading1
  CHARACTER*100  :: chr
  REAL*8      :: real_io
! Local/Loop
  INTEGER        :: iatom, ix, iy, iz, ia, ib, ic, idimension, imax
! Local/Calculating
  INTEGER     :: Nb_atoms_aux
  INTEGER     :: grid_aux(dimensions)
  REAL*8      :: X0Cell_aux(dimensions)            ! position of the origin of the volumetric data. 
  REAL*8      :: cell_aux(dimensions,dimensions)  ! cell vectors
  INTEGER        :: atoms_type_aux(Nb_atoms_check)     !atoms_type( Nb_atoms )
  REAL*8         ::  atoms_position_aux(Nb_atoms_check,dimensions)
  REAL*8         :: dV_aux
  REAL*8         :: VolCel_aux
  REAL*8         :: Calculated_Volume


  REAL*8         :: Crossaux(dimensions), aux(dimensions)

!***********************************************
!               MAIN BODY
!***********************************************
 ! READ
        PRINT*, "...READING Incoming Cube File", trim(FileInput)
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        READ(iufile,*) Nb_atoms_aux, X0Cell_aux(1),  X0Cell_aux(2),  X0Cell_aux(3) !    60  -40.000000  -30.000000  -30.000000
        PRINT'(A43,I5)', "    ...Number of atoms read              : ", Nb_atoms_aux
        Nb_atoms_aux=Nb_atoms_aux-1
        PRINT'(A43,I5)', "    ...Number of actual atoms in the cell: ", Nb_atoms_aux

        PRINT'(A43,I5)', "    ...Number of atoms reference: ", Nb_atoms_check

        PRINT'(A25,3F12.6)', "    ...Origin         : ", X0Cell_aux(:)
        READ(iufile,*) grid_aux(1), cell_aux(1,1), cell_aux(1,2), cell_aux(1,3)    !  300    0.267559    0.000000    0.000000
        READ(iufile,*) grid_aux(2), cell_aux(2,1), cell_aux(2,2), cell_aux(2,3)    !  300    0.000000    0.200669    0.000000
        READ(iufile,*) grid_aux(3), cell_aux(3,1), cell_aux(3,2), cell_aux(3,3)    !   50    0.000000    0.000000    1.632653
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid    : ",  grid_aux(1), " x " ,  grid_aux(2), " x ",  grid_aux(3)
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid Ref: ",  grid_check(1), " x " ,  grid_check(2), " x ",  grid_check(3)
        DO idimension = 1,dimensions
           IF (grid_check(idimension) /= grid_aux(idimension)) STOP
        ENDDO
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...A vector: ",  grid_aux(1)*cell_aux(1,1), " | " ,  grid_aux(1)*cell_aux(1,2), " | " ,  grid_aux(1)*cell_aux(1,3) 
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector: ",  grid_aux(2)*cell_aux(2,1), " | " ,  grid_aux(2)*cell_aux(2,2), " | " ,  grid_aux(2)*cell_aux(2,3) 
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector: ",  grid_aux(3)*cell_aux(3,1), " | " ,  grid_aux(3)*cell_aux(3,2), " | " ,  grid_aux(3)*cell_aux(3,3) 
        Crossaux(1) = (grid_aux(1)*cell_aux(1,2)*grid_aux(2)*cell_aux(2,3))  -(grid_aux(2)*cell_aux(2,2)*grid_aux(1)*cell_aux(1,3))
        Crossaux(2) = (grid_aux(1)*cell_aux(1,3)*grid_aux(2)*cell_aux(2,1))  -(grid_aux(2)*cell_aux(2,3)*grid_aux(1)*cell_aux(1,1))
        Crossaux(3) = (grid_aux(1)*cell_aux(1,1)*grid_aux(2)*cell_aux(2,2))  -(grid_aux(2)*cell_aux(2,1)*grid_aux(1)*cell_aux(1,2))
        aux(1) = grid_aux(3)*cell_aux(3,1)
        aux(2) = grid_aux(3)*cell_aux(3,2)
        aux(3) = grid_aux(3)*cell_aux(3,3)
        VolCel_aux = abs(dot_product(Crossaux(:), aux(:)))
        dV_aux = VolCel_aux / (grid_aux(1)*grid_aux(2)*grid_aux(3))
        PRINT'(A35,F12.4)', "    ...Volume of the cell[Bohr^3]: ", VolCel_aux
        PRINT'(A35,F12.6)', "    ...Voxel[Bohr^3]:              ", dV_aux

        IF ( Nb_atoms_aux /= Nb_atoms_check ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the number of actual atoms in the cell"
            PRINT*, "Reference: ", Nb_atoms_check
            PRINT*, "Read     : ", Nb_atoms_aux
            STOP
        ENDIF
        IF ( (ABS(X0Cell_aux(1)-X0Cell_check(1)) +  ABS(X0Cell_aux(2)-X0Cell_check(2))+ ABS(X0Cell_aux(3)-X0Cell_check(3))) > EPS_m5 ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the origin of the cell"
            PRINT*, "Reference: ", X0Cell_check(:)
            PRINT*, "Read     : ", X0Cell_aux(:)
            STOP
        ENDIF
       DO idimension = 1,dimensions
           IF (grid_check(idimension) /= grid_aux(idimension)) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the grid"
            PRINT*, "Reference: ", grid_check(:)
            PRINT*, "Read     : ", grid_aux(:)
            STOP
           ENDIF
           IF ( (ABS(cell_aux(idimension,1)- cell_check(idimension,1) ) + ABS(cell_aux(idimension,2)- cell_check(idimension,2))+ ABS(cell_aux(idimension,3)- cell_check(idimension,3))) > EPS_m5  ) THEN
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the lattice vector", idimension
            PRINT*, "Reference: ", cell_check(idimension,:)
            PRINT*, "Read     : ", cell_aux(idimension,:)
            STOP
           ENDIF
        ENDDO
        IF ( ABS(VolCel_aux-VolCel_check) > EPS_m5 ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the Volume of the cell"
            PRINT*, "Reference: ", VolCel_check
            PRINT*, "Read     : ", VolCel_aux
            STOP
        ENDIF

        PRINT*, "    ...READING Atoms type and positions"
     
    DO iatom = 1, Nb_atoms_aux
        READ(iufile,*)   atoms_type_aux(iatom), real_io, atoms_position_aux(iatom,1), atoms_position_aux(iatom,2), atoms_position_aux(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146
        IF (  atoms_type_check(iatom) /=  atoms_type_aux(iatom) ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the types of atom", iatom
            PRINT*, "Reference: ", atoms_type_check(iatom)
            PRINT*, "Read     : ", atoms_type_aux(iatom)
            STOP
        ENDIF
        IF (  (ABS(atoms_position_aux(iatom,1)- atoms_position_check(iatom,1) ) + ABS(atoms_position_aux(iatom,2)- atoms_position_check(iatom,2))+ ABS(atoms_position_aux(iatom,3)- atoms_position_check(iatom,3))) > EPS_m5   ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the positions of atom", iatom
            PRINT*, "Reference: ", atoms_position_check(iatom,:)
            PRINT*, "Read     : ", atoms_position_aux(iatom,:)
            STOP
        ENDIF
     ENDDO
     PRINT*, "    ...READING Hole position"
     READ(iufile,*)   iatom, real_io, holeposition(1), holeposition(2), holeposition(3)  ! 6    0.000000    2.924528    7.785081    1.274146
     PRINT'(A18,I5,A6,F12.6,A1,F12.6,A1,F12.6,A2)', "        ...Hole: ",iatom," at (",holeposition(1),",",holeposition(2),",",holeposition(3)," )"
       IF (   iatom /= 0  ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Problem in the specification of the hole"
            PRINT*, "Read     : ", iatom,  real_io, holeposition(:)
            STOP
        ENDIF

    Calculated_Volume=ZERO
    norm=ZERO
        PRINT*, "    ...READING Wavefunction"
    DO ix = 1,  grid_aux(1)
    DO iy = 1,  grid_aux(2)
    DO iz = 1,  grid_aux(3), 6 
        imax=MIN(5,(grid_aux(3)-iz))
!        PRINT*, ix, iy, iz, iz+imax
        READ(iufile,*)   wavefunction(ix,iy,iz:iz+imax) 
    ENDDO     
    ENDDO
    ENDDO
    DO ix = 1,  grid_aux(1)
    DO iy = 1,  grid_aux(2)
    DO iz = 1,  grid_aux(3) 
         norm=norm+ (dV_aux*wavefunction(ix,iy,iz))
         Calculated_Volume = Calculated_Volume+dV_aux
    ENDDO     
    ENDDO
    ENDDO
        PRINT'(A35,ES12.4)', "    ...Integrated Volume [Bohr^3]: ", Calculated_Volume
        PRINT'(A35,ES12.4)', "    ...Volume of the cell[Bohr^3]: ", VolCel_aux
        PRINT'(A35,ES12.6)', "    ...Norm:              ", norm
        PRINT'(A35,ES12.6)', "    ...Norm/Volume[Bohr^-3]      :", (norm/ Calculated_Volume)
        PRINT'(A35,ES12.6)', "    ...FFT*Norm/Volume[Bohr^-3]  :", (grid_aux(1)*grid_aux(2)*grid_aux(3)*norm/ Calculated_Volume)

        IF ( ABS(Calculated_Volume-VolCel_check) > EPS_m1 ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the Calculated Volume of the cell"
            PRINT*, "Reference : ", VolCel_check
            PRINT*, "Calculated: ", Calculated_Volume
            STOP
        ENDIF



        PRINT*, "            ...Done"
        PRINT*, "    ...CLOSING Input File"
        CLOSE(iufile)


END  SUBROUTINE GaussianFileread

!***********************************************
      SUBROUTINE gnuplotwrite(iufile, datafile, dimensions, nvect1, nvect2, vector1, vector2 , function_in)
   !***********************************************
   IMPLICIT NONE
! ...   Numerical constants
        
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0
   CHARACTER*12, PARAMETER :: subroutinename="gnuplotwrite"

   INTEGER, INTENT(in)              :: iufile
   CHARACTER*100, INTENT(in)        :: datafile
   CHARACTER*100  :: chr
   INTEGER,   INTENT(in)           :: dimensions
   INTEGER,     INTENT(in)         :: nvect1
   INTEGER,     INTENT(in)         :: nvect2
   REAL*8,   INTENT(in)         :: vector1(dimensions)
   REAL*8,   INTENT(in)         :: vector2(dimensions)
   REAL*8,   INTENT(in)         :: function_in(nvect1, nvect2) 
   REAL*8   :: X_value, Y_value
   INTEGER :: ia, ib     

       PRINT*, "    ...PRINTING Gnuplot Output File: ", TRIM(datafile)
       PRINT*, "       ...Mesh ", nvect1," x ", nvect2


        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        DO ia = 1,  nvect1
           DO ib = 1,  nvect2
                 X_value = (ia-1)* vector1(1)+(ib-1)* vector2(1)
                 Y_value = (ia-1)* vector1(2)+(ib-1)* vector2(2)
                 WRITE(iufile,'(3e15.5)') X_value, Y_value, function_in(ia,ib)  
           ENDDO
           WRITE(iufile,*)
        ENDDO


        CLOSE(iufile)


END  SUBROUTINE gnuplotwrite




!***********************************************
      SUBROUTINE GaussianFilewrite(FileInput, iufile, dimensions, Nb_atoms, grid, X0Cell, cell,atoms_type,atoms_position, wavefunction)
   !***********************************************

   IMPLICIT NONE
  !..   Numerical constants
        
        REAL*8, PARAMETER ::    ZERO = 0.0d0
        REAL*8, PARAMETER ::     ONE = 1.0d0
        REAL*8, PARAMETER ::      EPS_m5  = 0.00001d0
   CHARACTER*17, PARAMETER :: subroutinename="GaussianFilewrite"
  REAL*8, PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   REAL*8, PARAMETER ::  Ryd_to_eV = 13.605698066*ONE
!
 



   CHARACTER*100, INTENT(inout)        :: FileInput
   INTEGER, INTENT(in)              :: iufile

   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: Nb_atoms
   INTEGER,   INTENT(in)            :: grid(dimensions)
   REAL*8,   INTENT(in)         :: X0Cell(dimensions)            ! position of the origin of the volumetric data. 
   REAL*8,   INTENT(in)         :: cell(dimensions,dimensions)  ! cell vectors
   INTEGER,     INTENT(in)         :: atoms_type(Nb_atoms)     !atoms_type( Nb_atoms )
   REAL*8,   INTENT(in)         :: atoms_position(Nb_atoms,dimensions)
   REAL*8,   INTENT(in)         :: wavefunction(grid(1),grid(2),grid(3)) !wavefunction( grid(1),grid(2),grid(3) )

! Local/Reading1
  CHARACTER*100  :: chr
  REAL*8      :: real_io
! Local/Loop
  INTEGER        :: iatom, ix, iy, iz, ia, ib, ic, idimension, imax
! Local/Calculating
  REAL*8         :: dV,  R(dimensions)
  REAL*8         :: VolCel, distancetoat, distancemin
  REAL*8         :: Crossaux(dimensions), aux(dimensions)



       PRINT*, "    ...PRINTING Output File"
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        WRITE(iufile,*) TRIM(FileInput)                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        WRITE(iufile,*) TRIM(FileInput)                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        WRITE(iufile,'(I5,3F12.6)') Nb_atoms, X0Cell(1),  X0Cell(2),  X0Cell(3) !    60  -40.000000  -30.000000  -30.000000
        WRITE(iufile,'(I5,3F12.6)') grid(1), cell(1,1), cell(1,2), cell(1,3)    !  300    0.267559    0.000000    0.000000
        WRITE(iufile,'(I5,3F12.6)') grid(2), cell(2,1), cell(2,2), cell(2,3)    !  300    0.000000    0.200669    0.000000
        WRITE(iufile,'(I5,3F12.6)') grid(3), cell(3,1), cell(3,2), cell(3,3)    !   50    0.000000    0.000000    1.632653
        Crossaux(1) = (grid(1)*cell(1,2)*grid(2)*cell(2,3))  -(grid(2)*cell(2,2)* grid(1)*cell(1,3) )
        Crossaux(2) = (grid(1)*cell(1,3)*grid(2)*cell(2,1))  -(grid(2)*cell(2,3)*grid(1)*cell(1,1))
        Crossaux(3) = (grid(1)*cell(1,1)*grid(2)*cell(2,2))  - (grid(2)*cell(2,1)*grid(1)*cell(1,2))
        aux(1) = grid(3)*cell(3,1)
        aux(2) = grid(3)*cell(3,2)
        aux(3) = grid(3)*cell(3,3)
        VolCel = abs(dot_product(Crossaux(:), aux(:)))
        dV = VolCel / (grid(1)*grid(2)*grid(3))


    DO iatom = 1,Nb_atoms
        WRITE(iufile,'(I5,4F12.6)')   atoms_type(iatom), real(atoms_type(iatom)), atoms_position(iatom,1), atoms_position(iatom,2), atoms_position(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146
    ENDDO
    DO ia = 1,  grid(1)
    DO ib = 1,  grid(2)
    DO ic = 1,  grid(3), 6
        imax=MIN(5,(grid(3)-ic)) 
           WRITE(iufile,'(6 ES13.5)')   wavefunction(ia,ib,ic:ic+imax) 
    ENDDO     
    ENDDO
    ENDDO
        PRINT*, "            ...Done"
        PRINT*, "    ...CLOSING Output File"
        CLOSE(iufile)

END  SUBROUTINE GaussianFilewrite


!***********************************************
      SUBROUTINE XYZFilewrite(FileOutput, iufile, dimensions, Nb_atoms_check, atoms_type,atoms_position)
   !***********************************************
   IMPLICIT NONE

   CHARACTER*100, INTENT(in)        :: FileOutput
   INTEGER, INTENT(in)              :: iufile
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: Nb_atoms_check
   !
   INTEGER,   INTENT(in)           :: atoms_type(Nb_atoms_check)     !atoms_type( Nb_atoms )
   REAL*8,   INTENT(in)         :: atoms_position(Nb_atoms_check,dimensions)
! Local Variables
 CHARACTER*2        :: chr
   INTEGER :: iatom

        PRINT*, "...Writing XYZ File: "
        OPEN(iufile,file=trim(FileOutput),form='formatted')
        REWIND(iufile)
        WRITE(iufile,*) Nb_atoms_check               
        DO iatom = 1, Nb_atoms_check
           IF ( atoms_type(iatom) == 1 ) THEN 
               WRITE(chr,'(A2)') 'H '
           ELSE IF ( atoms_type(iatom) == 2) THEN 
               WRITE(chr,'(A2)') 'He'
           ELSE IF ( atoms_type(iatom) == 3) THEN 
               WRITE(chr,'(A2)') 'Li'
           ELSE IF ( atoms_type(iatom) == 4) THEN 
               WRITE(chr,'(A2)') 'Be'
           ELSE IF ( atoms_type(iatom) == 5) THEN 
               WRITE(chr,'(A2)') 'B '
           ELSE IF ( atoms_type(iatom) == 6) THEN 
               WRITE(chr,'(A2)') 'C '
           ELSE IF ( atoms_type(iatom) == 7) THEN 
               WRITE(chr,'(A2)') 'N '
           ELSE IF ( atoms_type(iatom) == 8) THEN 
               WRITE(chr,'(A2)') 'O '
           ELSE IF ( atoms_type(iatom) == 9) THEN 
               WRITE(chr,'(A2)') 'F '
           ELSE IF ( atoms_type(iatom) == 10) THEN 
               WRITE(chr,'(A2)') 'Ne'
           ELSE IF ( atoms_type(iatom) == 11) THEN 
               WRITE(chr,'(A2)') 'Na'
           ELSE IF ( atoms_type(iatom) == 12) THEN 
               WRITE(chr,'(A2)') 'Mg'
           ELSE IF ( atoms_type(iatom) == 13) THEN 
               WRITE(chr,'(A2)') 'Al'
           ELSE IF ( atoms_type(iatom) == 14) THEN 
               WRITE(chr,'(A2)') 'Si'
           ELSE IF ( atoms_type(iatom) == 15) THEN 
               WRITE(chr,'(A2)') 'P '
           ELSE IF ( atoms_type(iatom) == 16) THEN 
               WRITE(chr,'(A2)') 'S '
           ELSE IF ( atoms_type(iatom) == 17) THEN 
               WRITE(chr,'(A2)') 'Cl'
           ELSE IF ( atoms_type(iatom) == 18) THEN 
               WRITE(chr,'(A2)') 'Ar'
           ELSE IF ( atoms_type(iatom) == 79) THEN 
               WRITE(chr,'(A2)') 'Au'
           ELSE 
               PRINT*, "Warning, not identified"
               WRITE(chr,'(A2)') 'X '
           ENDIF
           WRITE(iufile,'(A2,A3,3F12.6)')   chr,"   ", atoms_position(iatom,1), atoms_position(iatom,2), atoms_position(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146 
        ENDDO
        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE XYZFilewrite



















!***********************************************
! PRINT
!***********************************************

!        grid_aux(1)=-1+(2*grid(1))
!        grid_aux(2)=-1+(2*grid(2))
!        grid_aux(3)=-1+(2*grid(3))
!
!        PRINT*, "......Writing 3D Correlation on xyz"
!        outputtype="3D-xyz"
!        WRITE(outputfname, *) 
!        !WRITE(outputfname, '(A100)')  "                                                                                                    "
!        WRITE(outputfname,  '(A50,A1,A6,A5)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".cube"        
!        PRINT*, "...Writing to file:", TRIM(outputfname)
!        CALL GaussianFilewrite(TRIM(outputfname), iufile, dimensions, Nb_atoms_aux, grid_aux, X0Cell, cell,atoms_type_aux,atoms_position_aux1, twopartcorr)
!        PRINT*, "...... Done"
!!***********************************************
!        PRINT*, "......Writing 3D Wavefunction on xyz"
!        outputtype="3D_DEN"
!        WRITE(outputfname, *) 
!        !WRITE(outputfname, '(A100)')  "                                                                                                    "
!        WRITE(outputfname,  '(A50,A1,A6,A5)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".cube"        
!        PRINT*, "...Writing to file:", TRIM(outputfname)
!        CALL GaussianFilewrite(TRIM(outputfname), iufile, dimensions, Nb_atoms_aux, grid, X0Cell, cell,atoms_type_aux,atoms_position_aux1, AvgElectronicDensity)
!        PRINT*, "...... Done"
!
!
!!***********************************************
!        PRINT*, "......Writing 3D Correlation on abc"
!        outputtype="3D-abc"
!        !WRITE(outputfname, '(A100)')  "                                                                                                    "
!        WRITE(outputfname, '(A50,A1,A6,A5)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".cube"        
!        PRINT*, "...Writing to file:", TRIM(outputfname)
!        CALL GaussianFilewrite(TRIM(outputfname), iufile, dimensions, Nb_atoms_aux, grid_aux, vect3_aux, mat3x3_aux,atoms_type_aux,atoms_position_aux2, twopartcorr)
!        PRINT*, "...... Done"
!!***********************************************
!       PRINT*, "... Output 3-d files:"
!!***********************************************
!        PRINT*, "......Writing 3D Correlation on xyz"
!        outputtype="3D-gen"
!       ! WRITE(outputfname, '(A55)')  "                                                                                                    "
!        WRITE(outputfname, '(A50,A1,A6,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".txt"        
!        PRINT*, "...Writing to file:", TRIM(outputfname)
!        OPEN(iufile,file=trim(outputfname),form='formatted')
!           do ia=(1-grid(1)),(grid(1)-1)
!              do ib=(1-grid(2)),(grid(2)-1)
!                 do ic=(1-grid(3)),(grid(3)-1)
!		    ! distance between electron and hole
!                    ElectronPosition(:)=ZERO
!                    vect3_aux(1)=(real(ia))
!                    vect3_aux(2)=(real(ib))
!                    vect3_aux(3)=(real(ic))
!                    ElectronPosition(:)= matmul(latticetoxyztransfermatrix(:,:),vect3_aux(:))  + X0Cell(:)
!                    WRITE(iufile,'(5ES13.5,I6)')  ElectronPosition(:), twopartcorr(ia+grid(1),ib+grid(2),ic+grid(3)), cdf(ia+grid(1),ib+grid(2),ic+grid(3)), counter(ia+grid(1),ib+grid(2),ic+grid(3))
!          	 ENDDO
!       		ENDDO
!	ENDDO
!        CLOSE(iufile)
!
!
!***********************************************
!        PRINT*, "...... Done"
!***********************************************




