!
! Copyright (C) 2012 The Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre Darancet
!***********************************************
!
! This program reads two cube files containing 3D functions
! Performs the two particle correlation functions
! Does the ratio 
! Calculates the axis, planar, distance average

   PROGRAM CalculateChargeTransfer
   !***********************************************
    ! USES SUBROUTINES:
    !invert3x3
    !_
    !GaussianFileread
    !ReadInput    !gnuplotwrite
    !GaussianFilereadHeader1
    !GaussianFilereadHeader2
    !GaussianFilewrite
    !GaussianFileread
  
   !***********************************************
   IMPLICIT NONE
   !***********************************************
   !** Parameters
   CHARACTER*5, PARAMETER :: VersionNumber="0.0.1"  
   Integer, parameter :: maxchar_in = 55
   Integer, parameter :: maxchar_out = 100
   INTEGER, PARAMETER :: dimensions = 3    
   INTEGER, PARAMETER :: iufile=100
   !REAL*8, PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   CHARACTER*15, PARAMETER :: datafile="INPUT_CTCALC.in"
   CHARACTER*23, PARAMETER :: subroutinename="CalculateChargeTransfer"

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

   REAL*8, ALLOCATABLE ::   twopartcorr_Volume(:,:,:)
   REAL*8, ALLOCATABLE ::   twopartcorr_distance(:),    twopartcorr_distance_Volume(:)


   !REAL*8, ALLOCATABLE ::   cdf(:,:,:)
   REAL*8 :: Normtwopartcorr
   REAL*8, ALLOCATABLE ::   cdf_distance_renormalized(:)
   INTEGER, ALLOCATABLE ::  counter(:,:,:)
   INTEGER, ALLOCATABLE ::  counter_distance(:)

   !** Local Variables
   ! Calculated Quantities
   !   
   INTEGER ::   Nb_Distances

   REAL*8  ::   DistanceStep
   REAL*8 ::   VectorPosition(dimensions)
   REAL*8 ::   ElectronHoleDistanceR(dimensions)
   REAL*8 ::   CutoffDistance
   REAL*8 ::   AverageDistance
   REAL*8 ::   INVOLUMEFRACTION
   REAL*8 ::   Dipole(dimensions)
   REAL*8 ::   Quadrupole(dimensions,dimensions)
   LOGICAL :: INTHEVOLUME
   REAL*8 ::   diagonal
   REAL*8 ::   distance
   REAL*8 ::   VolCel, dV
   REAL*8 ::   latticetoxyztransfermatrix(dimensions,dimensions)
   REAL*8 ::    xyztolatticetransfermatrix(dimensions,dimensions)
   ! Input Quantities
   !   
   REAL*8 :: norm(2)
   REAL*8, ALLOCATABLE ::   twopartcorr_INPUT(:,:,:)
   REAL*8 ::   cell(dimensions,dimensions)
   REAL*8 ::   X0Cell(dimensions)
   INTEGER ::   grid(dimensions)
   INTEGER ::   Nb_atoms

   REAL*8, ALLOCATABLE  ::   atoms_position(:,:)
   INTEGER, ALLOCATABLE ::  atoms_type(:)
   CHARACTER*(maxchar_in) ::fname=' ', NameInputFile=' '
   CHARACTER*(maxchar_out) ::GenericNameOutputFile=' ', outputfname=' ', outputtype=' '

   CHARACTER*(maxchar_out) :: Model 
   REAL*8 ::   axis_a, axis_b, axis_c


   ! Loop variables
   ! 
   INTEGER :: ia, ib, ic, id, id2, ia_aux, ib_aux, ic_aux, ih, iatom

   ! Auxiliary variables
   ! 
   REAL*8, ALLOCATABLE  :: vect1(:), vect2(:), TWODFunction_aux1(:,:), TWODFunction_aux2(:,:), TWODFunction_aux3(:,:),ONEDFunction_aux1(:), ONEDFunction_aux2(:), ONEDFunction_aux3(:)
   REAL*8 :: Maximum_norm, vect_norm(dimensions), mat3x3_aux(dimensions,dimensions), vect3_aux(dimensions), vect3_aux2(dimensions), MiddleOfTheCell(3), real_aux, real_aux2, ElectronHoleDistance
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
   PRINT*, " Contributors  : Pierre Darancet"
   PRINT*, "**************************************************"
   PRINT*, " This program reads two cube files containing 3D functions"
   PRINT*, " performs the two particle correlation functions"
   PRINT*, " calculates the ratio, and the axis-, planar- and distance-average"
   PRINT*, "**************************************************"
   PRINT*, " Version", TRIM(VersionNumber)
   PRINT*, " The general input file is: ",  TRIM(datafile)
   PRINT*, "**************************************************"
   PRINT*, " Beginning..."
!  !
   ! read input General Information
   !
   PRINT*, " Reading generic input file..."
   CALL ReadInputFile(iufile, datafile, NameInputFile,Model, axis_a, axis_b, axis_c ,GenericNameOutputFile)
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
   PRINT*, "... The program will be read from files:"
   PRINT*,  TRIM(NameInputFile), " for the Numerator"
   !  Electronic GRID
   PRINT*, " Reading first file for allocation..."
   fname=' '
   PRINT*, " File...", TRIM(NameInputFile)
   CALL  GaussianFilereadHeader1(trim(NameInputFile), iufile, dimensions, Nb_atoms, grid, X0Cell, cell, VolCel, dV)
   PRINT*, " ...done"
   !  Atomic GRID
   PRINT*, " Reading first file for atoms position..." 
   ALLOCATE(atoms_type(Nb_atoms))
   ALLOCATE(atoms_position(Nb_atoms,dimensions))
   CALL  GaussianFilereadHeader2(trim(NameInputFile), iufile, dimensions, Nb_atoms, grid, X0Cell, cell,atoms_type,atoms_position)
   !  Distance GRID
   PRINT*, " Calculating the distance grid..."
   diagonal=0.5*sqrt((abs(grid(1)*cell(1,1))  +  abs(grid(2)*cell(2,1)) +  abs(grid(3)*cell(3,1)))**2 +(abs(grid(1)*cell(1,2))  +  abs(grid(2)*cell(2,2)) +  abs(grid(3)*cell(3,2)))**2 +(abs(grid(1)*cell(1,3))  +  abs(grid(2)*cell(2,3)) +  abs(grid(3)*cell(3,3)))**2)
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
   ALLOCATE(counter(grid(1),grid(2),grid(3)))
   ALLOCATE(twopartcorr_Volume(grid(1),grid(2),grid(3)))
   ALLOCATE(twopartcorr_INPUT(grid(1),grid(2),grid(3)))
   ALLOCATE(twopartcorr_distance_Volume(Nb_Distances))
   ALLOCATE(twopartcorr_distance(Nb_Distances))
   ALLOCATE(cdf_distance_renormalized(Nb_Distances))
   ALLOCATE(counter_distance(Nb_Distances))
   PRINT*, " ...done"
   !  Printing Summary
   PRINT*, "... Largest distance in the supercell [Bohr]:", diagonal
   PRINT'(A35,3F12.4)', "... Norm of Lattice Vector [Bohr]: ", (vect_norm(id), id=1,dimensions)
   PRINT*, "... Maximal Norm of Lattice Vector [Bohr]:  ", Maximum_norm
   PRINT*, "... Step in distance [Bohr]:  ", DistanceStep
   PRINT*, "... Number of steps        :  ", Nb_Distances
   PRINT*, "... Largest distance on the grid [Bohr]:", DistanceStep*(Nb_Distances-1)
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

   twopartcorr_Volume(:,:,:)=ZERO
   twopartcorr_INPUT(:,:,:)=ZERO
   twopartcorr_distance_Volume(:)=ZERO
   twopartcorr_distance(:)=ZERO
   cdf_distance_renormalized(:)=ZERO
   
   counter(:,:,:)=0
   counter_distance(:)=0
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
    !***********************************************
    ! Reading Electronic Density
    PRINT*, "...Reading from file:", TRIM(NameInputFile)
    CALL GaussianFileread(trim(NameInputFile), iufile, dimensions, Nb_atoms,  grid, X0Cell, cell, atoms_type, atoms_position,   dV, VolCel,  twopartcorr_INPUT(:,:,:),norm(1))
    PRINT*, " ...done"
 
    !***********************************************

    !Calculating the middle of the cell for shifting purpose
    ! Coming out of the 2 part corr, the number of grid points should be odd, if not the code prints out an error message and leave
    IF ((MOD(grid(1),2)==0).OR.(MOD(grid(2),2)==0).OR.(MOD(grid(3),2)==0)) THEN
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Problem in the grid"
            PRINT*, "Grid expected to be odd values: ", grid(1), grid(2), grid(3)
            STOP
    ENDIF
    vect3_aux(1)=(real(grid(1))/TWO)
    vect3_aux(2)=(real(grid(2))/TWO)
    vect3_aux(3)=(real(grid(3))/TWO)
    MiddleOfTheCell(:)= matmul(latticetoxyztransfermatrix(:,:),vect3_aux(:)) +X0Cell(:) 

    !***********************************************
    !
    !------------------------------
    ! Loop over electron position
    !------------------------------
    PRINT*, "... Beginning of the position loop"
    DO ia=1,grid(1)
      DO ib=1,grid(2)
         DO ic=1,grid(3)
	    !***********************************************
	    !
	    !------------------------------
	    ! Distance between electron and hole
	    !------------------------------
            ! Calculating electron position in xyz coordinates
            VectorPosition(:)=ZERO
            vect3_aux(1)=(real(ia))
            vect3_aux(2)=(real(ib))
            vect3_aux(3)=(real(ic))
            ! The grid must be shifted so that the middle of the cell is considered to be (0.0,0.0,0.0)
            VectorPosition(:)= matmul(latticetoxyztransfermatrix(:,:),vect3_aux(:)) +X0Cell(:) - MiddleOfTheCell(:)
            distance=ZERO   
            DO id=1,dimensions
                 distance=distance+ (( VectorPosition(id) )**2)
            ENDDO
            distance=sqrt(distance)
            counter(ia,ib,ic)= counter(ia,ib,ic)+1                   
	    
	    !***********************************************
	    !
	    !------------------------------
	    ! Test if the position belongs to the volume 
	    !------------------------------
	    INTHEVOLUME=.FALSE.
	    IF ( TRIM(Model) == "Ellipsoid" ) THEN
		IF (  ((VectorPosition(1)/axis_a)**2 + (VectorPosition(2)/axis_b)**2 + (VectorPosition(3)/axis_c)**2) .LE. ONE ) INTHEVOLUME=.TRUE.
	    ELSE IF ( TRIM(Model) == "Cylinder" ) THEN
		IF ((  ( (VectorPosition(1)/axis_a)**2 + (VectorPosition(2)/axis_b)**2 ) .LE. ONE ) .AND. ( abs(VectorPosition(3)/axis_c) .LE. (ONE/TWO)) ) INTHEVOLUME=.TRUE.
	    ELSE

	    ENDIF 

		
	    !***********************************************
	    !
	    !------------------------------
	    ! Calculation of the 2 part corr ratio and cdf
	    !------------------------------
            IF ( INTHEVOLUME ) twopartcorr_Volume(ia,ib,ic)=twopartcorr_INPUT(ia,ib,ic)
            ! Distance
            DO id=1,Nb_Distances
            	real_aux=(id-1)*DistanceStep                                                                        ! Find the 
            	IF ((distance.lt.real_aux+DistanceStep).AND.(distance.ge.real_aux)) THEN
          		twopartcorr_distance(id) = twopartcorr_distance(id) +  twopartcorr_INPUT(ia,ib,ic)  ! Store the component 
         		counter_distance(id)=counter_distance(id)+1   
			IF ( INTHEVOLUME ) twopartcorr_distance_Volume(id)=twopartcorr_distance_Volume(id)+  twopartcorr_INPUT(ia,ib,ic) 
               ENDIF
            ENDDO ! Nb_distances
            ! Vector


!                    if ( (ia.eq.22).or.(ib.eq.22).or.(ic.eq.22) ) PRINT*, "counter....", counter(ia_aux,ib_aux,ic_aux) 
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


        PRINT*, "... Calculating the CDF"
        DO id=1,Nb_Distances 
   	      IF (counter_distance(id).eq.0) then
              	twopartcorr_distance(id)=ZERO
       	      ELSE
        	twopartcorr_distance(id)=twopartcorr_distance(id)/(REAL(counter_distance(id)))
       	      ENDIF
       	      cdf_distance_renormalized(id)=SUM(twopartcorr_distance(1:id))
        ENDDO !Nb_distances
        twopartcorr_distance(:)=twopartcorr_distance(:)/cdf_distance_renormalized(Nb_Distances)
        twopartcorr_distance_Volume(:)=twopartcorr_distance_Volume(:)/cdf_distance_renormalized(Nb_Distances)
        PRINT*, "... Done"
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
	      WRITE(iufile,'(3ES13.5,I6)')  distance,  twopartcorr_distance(id), cdf_distance_renormalized(id), counter_distance(id)
	ENDDO
	CLOSE(iufile)
        PRINT*, "... Done"

	PRINT*, "...Writing 1D Correlation Function"
	outputtype="1D-distance-involume"
	! WRITE(outputfname, '(A55)')  "                                                                                                    "
	WRITE(outputfname, '(A41,A1,A20,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".dat"
	PRINT*, "...Writing to file:", TRIM(outputfname)
	OPEN(iufile,file=trim(outputfname),form='formatted')
	REWIND(iufile)
	WRITE(iufile,*) "# Distance [Bohr], 2-part corr, cdf, NbPoints, cdf_from_num"
	DO id=1,Nb_Distances 
	      ! Radius of the sphere
	      distance= (id-1)*DistanceStep
	      WRITE(iufile,'(3ES13.5,I6)')  distance,  twopartcorr_distance_Volume(id), twopartcorr_distance(id), counter_distance(id)
	ENDDO
	CLOSE(iufile)
        PRINT*, "... Done"
        ! Done
        !***********************************************
        ! Two Part Corr (a,b,c) - recalculate with counter 0 corrections
        PRINT*, "... Renormalizing twopartcorr_abc"
        Normtwopartcorr=ZERO
        DO ic=1,grid(3)
		DO ib=1,grid(2)
		DO ia=1,grid(1)
                      Normtwopartcorr=Normtwopartcorr+twopartcorr_INPUT(ia,ib,ic)
		ENDDO
		ENDDO
	ENDDO
        PRINT*, "... Done"
        twopartcorr_Volume(:,:,:)=twopartcorr_Volume(:,:,:)/Normtwopartcorr
	twopartcorr_INPUT(:,:,:)=twopartcorr_INPUT(:,:,:)/Normtwopartcorr
        !***********************************************
        ! Two Part Corr (a,b,c) - recalculate with counter 0 corrections

        !***********************************************
        ! Two Part Corr Colume (a,b,c) - Calculate In volume fraction
        PRINT*, "... Calculaing In volume fraction"

	INVOLUMEFRACTION=ZERO
        DO ic=1,grid(3)
		DO ib=1,grid(2)
		DO ia=1,grid(1)
			INVOLUMEFRACTION=INVOLUMEFRACTION+twopartcorr_Volume(ia,ib,ic)
		ENDDO
		ENDDO
	ENDDO
        PRINT*, "... Done"
        !***********************************************
        PRINT*, "Volume Model (a,b,c)",  TRIM(Model), axis_a, axis_b, axis_c
        PRINT*, "FRACTION OF 2partCorr in volume:",  INVOLUMEFRACTION

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
        ALLOCATE(TWODFunction_aux1(grid(2),grid(3)))
	ALLOCATE(TWODFunction_aux2(grid(1),grid(3)))
        ALLOCATE(TWODFunction_aux3(grid(1),grid(2)))
        ALLOCATE(ONEDFunction_aux1(grid(1)))
        ALLOCATE(ONEDFunction_aux2(grid(2)))
        ALLOCATE(ONEDFunction_aux3(grid(3)))
        ALLOCATE(vect1(2))
        ALLOCATE(vect2(2))
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
              ! ss:
        	!AverageDistance=AverageDistance+ ( twopartcorr_distance(id)*distance*DistanceStep )
              ! PD
         	AverageDistance=AverageDistance+ ( twopartcorr_distance(id)*distance )
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
        do ic=1,grid(3)
           do ib=1,grid(2)
             do ia=1,grid(1)
                    vect3_aux(:)=ZERO 
	            vect3_aux2(1)=real(ia)
	            vect3_aux2(2)=real(ib)
	            vect3_aux2(3)=real(ic)
	            vect3_aux(:)= matmul(latticetoxyztransfermatrix(:,:),vect3_aux2(:))  +X0Cell(:) - MiddleOfTheCell(:)
		    Dipole(:)=Dipole(:)+ (vect3_aux(:)*twopartcorr_Volume(ia,ib,ic))
                    distance=sqrt((vect3_aux(1)**2) +(vect3_aux(2)**2)+(vect3_aux(3)**2))
                    Quadrupole(1,1)=Quadrupole(1,1)+(twopartcorr_Volume(ia,ib,ic)*( (3*(vect3_aux(1)**2)) - (distance**2) ))
                    Quadrupole(2,2)=Quadrupole(2,2)+(twopartcorr_Volume(ia,ib,ic)*( (3*(vect3_aux(2)**2)) - (distance**2) ))
                    Quadrupole(3,3)=Quadrupole(3,3)+(twopartcorr_Volume(ia,ib,ic)*( (3*(vect3_aux(3)**2)) - (distance**2) ))

                    Quadrupole(2,1)=Quadrupole(2,1)+(twopartcorr_Volume(ia,ib,ic)*(3*vect3_aux(2)*vect3_aux(1)))
                    Quadrupole(3,1)=Quadrupole(3,1)+(twopartcorr_Volume(ia,ib,ic)*(3*vect3_aux(3)*vect3_aux(1)))
                    Quadrupole(3,2)=Quadrupole(3,2)+(twopartcorr_Volume(ia,ib,ic)*(3*vect3_aux(3)*vect3_aux(2)))
	     ENDDO
	ENDDO
       ENDDO
       Quadrupole(1,2)=Quadrupole(2,1)
       Quadrupole(1,3)=Quadrupole(3,1)
       Quadrupole(2,3)=Quadrupole(3,2)
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
		DO ic=1,grid(3)
			   DO ib=1,grid(2)
				   DO  ia=1,grid(1)
				     TWODFunction_aux1(ib,ic) =  TWODFunction_aux1(ib,ic) + twopartcorr_Volume(ia,ib,ic)
				   ENDDO
			   ENDDO
		ENDDO
		! Average over b
		DO ic=1,grid(3)
			   DO ia=1,grid(1)
				   DO  ib=1,grid(2)
				     TWODFunction_aux2(ia,ic) =  TWODFunction_aux2(ia,ic) +  twopartcorr_Volume(ia,ib,ic)
				   ENDDO
			   ENDDO
		ENDDO
		! Average over c
		DO ia=1,grid(1)
			   DO  ib=1,grid(2)
				   DO ic=1,grid(3)
				     TWODFunction_aux3(ia,ib) =  TWODFunction_aux3(ia,ib) +  twopartcorr_Volume(ia,ib,ic)
				   ENDDO
			   ENDDO
		ENDDO
        ! Average 1D
		! Average over a AND b
		DO ic=1,grid(3)
			   DO ib=1,grid(2)
				     ONEDFunction_aux3(ic) =  ONEDFunction_aux3(ic) + TWODFunction_aux1(ib,ic)
			   ENDDO
		ENDDO
		! Average over b AND c
		DO ia=1,grid(1)
		 	   DO ic=1,grid(3)
				     ONEDFunction_aux1(ia) = ONEDFunction_aux1(ia)+ TWODFunction_aux2(ia,ic) 
			   ENDDO
		ENDDO
		! Average over a AND c
		DO ib=1,grid(2)
			   DO ic=1,grid(3)
				     ONEDFunction_aux2(ib) =  ONEDFunction_aux2(ib) + TWODFunction_aux1(ib,ic)
			   ENDDO
		ENDDO

 
!***********************************************
!***********************************************
!        PRINT*, "... Output xyz files:"
!***********************************************



		   PRINT*, "... Output xyz file:"
		   outputtype="Coord"
		   WRITE(outputfname, '(A55)')  "                                                                                                    "
		   WRITE(outputfname, '(A50,A1,A5,A4)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".xyz"
          	   CALL XYZFilewrite(outputfname, iufile, dimensions, Nb_atoms, atoms_type, atoms_position)
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
	WRITE(iufile,*) " Copyright (C) 2012-2013 The Molecular Foundry Berkeley"
	WRITE(iufile,*) " This file is distributed under the terms of the"
	WRITE(iufile,*) " GNU General Public License. See the file `License'"
	WRITE(iufile,*) " in the root directory of the present distribution,"
	WRITE(iufile,*) " or http://www.gnu.org/copyleft/gpl.txt ."
	WRITE(iufile,*) " Contributors  : Sahar Sharifzadeh, Pierre Darancet"
	WRITE(iufile,*) "**************************************************"
	WRITE(iufile,*) " Version Number: ", TRIM(VersionNumber)
	WRITE(iufile,*) "**************************************************"
	WRITE(iufile,*) " Input files: ",  TRIM(datafile), TRIM(NameInputFile)
	WRITE(iufile,*) "**************************************************"
	WRITE(iufile,*) " The general output file is: ",  TRIM(GenericNameOutputFile)
	WRITE(iufile,*) "Total Norm:", SUM(norm(1:2))
       WRITE(iufile,*) "Dipole [eBohr]:"
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)')  " ( ",Dipole(1)  , " | " ,Dipole(2)   , " | " ,Dipole(3)  , " ) "
       PRINT*, "Quadrupole [eBohr^2]:"
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)') " ( ",Quadrupole(1,1)  , " | " ,Quadrupole(1,2)   , " | " ,Quadrupole(1,3)  , " ) "
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)') " ( ",Quadrupole(2,1)  , " | " ,Quadrupole(2,2)   , " | " ,Quadrupole(2,3)  , " ) "
       WRITE(iufile,'(A3,F12.6,A3,F12.6,A3,F12.6,A3)') " ( ",Quadrupole(3,1)  , " | " ,Quadrupole(3,2)   , " | " ,Quadrupole(3,3)  , " ) "  
       WRITE(iufile,*) "Volume Model (a,b,c) [Bohr]",  TRIM(Model), " ( ", axis_a, " | ", axis_b, " | ", axis_c, " ) "
       WRITE(iufile,*) "FRACTION OF 2partCorr in volume:",  INVOLUMEFRACTION

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

        PRINT*, "... Output 2-d files:"
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
        PRINT*, "... Output 3-d files:"
!***********************************************

!***********************************************
! PRINT
!***********************************************

        PRINT*, "......Writing 3D Correlation on xyz"
        outputtype="3D-xyz"
        WRITE(outputfname, *) 
        !WRITE(outputfname, '(A100)')  "                                                                                                    "
        WRITE(outputfname,  '(A50,A1,A6,A5)') TRIM(GenericNameOutputFile),"_" ,  TRIM(outputtype), ".cube"        
        PRINT*, "...Writing to file:", TRIM(outputfname)
        CALL GaussianFilewrite(TRIM(outputfname), iufile, dimensions, Nb_atoms, grid, X0Cell, cell,atoms_type,atoms_position, twopartcorr_Volume)
        PRINT*, "...... Done"



   PRINT*, " ... twopartcorr"
   DEALLOCATE(twopartcorr_Volume)
   DEALLOCATE(twopartcorr_INPUT)

   PRINT*, " ... twopartcorr_distance"
   DEALLOCATE(twopartcorr_distance)
   DEALLOCATE(twopartcorr_distance_Volume)
   PRINT*, " ... cdf_distance_renormalized"
   DEALLOCATE(cdf_distance_renormalized)
   PRINT*, " ... counter"
   DEALLOCATE(counter)
   PRINT*, " ... counter_distance"
   DEALLOCATE(counter_distance)
   PRINT*, " ... TWODFunction_aux1"
   DEALLOCATE(TWODFunction_aux1)

   PRINT*, " ... TWODFunction_aux2"
   DEALLOCATE(TWODFunction_aux2)
   PRINT*, " ... TWODFunction_aux3"
   DEALLOCATE(TWODFunction_aux3)
   PRINT*, " ... ONEDFunction_aux1"
   DEALLOCATE(ONEDFunction_aux1)
   PRINT*, " ... ONEDFunction_aux2"
   DEALLOCATE(ONEDFunction_aux2)
   PRINT*, " ... ONEDFunction_aux3"
   DEALLOCATE(ONEDFunction_aux3)
   PRINT*, " ... vect1"
   DEALLOCATE(vect1)
   PRINT*, " ... vect2"
   DEALLOCATE(vect2)


   PRINT*, " ...done"
   PRINT*, " End"
   END PROGRAM CalculateChargeTransfer
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
      SUBROUTINE ReadInputFile(iufile, datafile,NameInputFile, Model, a_axis, b_axis, c_axis, GenericNameOutputFile)
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
   CHARACTER*15, INTENT(in)         :: datafile
   CHARACTER*(maxchar_in), INTENT(out)        :: NameInputFile
   CHARACTER*(maxchar_in), INTENT(out)        :: Model
   REAL*8, INTENT(out)            :: a_axis, b_axis, c_axis
   CHARACTER*(maxchar_out), INTENT(out)        :: GenericNameOutputFile


   CHARACTER*100  :: chr

        PRINT*, "...The Format of the input file is the following: "
        PRINT*, "   ...: Name for the incoming Cube file of 2 part corr"
        PRINT*, '   ...: Model used (Ellipsoid|Elliptic Cylinder):'
        PRINT*, '   ...: First  axis: a'
        PRINT*, '   ...: Second axis: b'
        PRINT*, '   ...: Thrid  axis: c'
        PRINT*, '   ...: Generic Name for the outcoming files:'
        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) NameInputFile
        PRINT*, "   ...READING: Name for the Cube file: ", Trim(NameInputFile)
        READ(iufile,*) Model
        PRINT*, "   ...READING: Model used (Ellipsoid|Cylinder): ", Trim(Model)
	IF ( ( TRIM(Model) == "Cylinder" ) .OR. ( TRIM(Model) == "Ellipsoid" ) ) THEN
	        PRINT*, "      ... Model supported! "
	ELSE 
            	PRINT*, "WARNING: Routine", TRIM(subroutinename)
	        PRINT*, "WARNING: Model Not supported"
        	PRINT*, "Model  : ", Model
		STOP
	ENDIF

        READ(iufile,*) a_axis
        PRINT*, "   ...READING: a axis length [Bohr]: ", a_axis
        READ(iufile,*) b_axis
        PRINT*, "   ...READING: b axis length [Bohr]: ", b_axis

	        READ(iufile,*) c_axis
        	PRINT*, "   ...READING: c axis length [Bohr]: ", c_axis

        READ(iufile,*) GenericNameOutputFile
        PRINT*, "   ...READING: Generic Name for the outputfiles: ", Trim(GenericNameOutputFile)

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
        PRINT*, "...READING Incoming Cube File: Header ", trim(FileInput)
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        READ(iufile,*) Nb_atoms, X0Cell(1),  X0Cell(2),  X0Cell(3) !    60  -40.000000  -30.000000  -30.000000
        PRINT'(A43,I5)', "    ...Number of atoms read              : ", Nb_atoms
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
        PRINT*, "...READING Incoming Cube File: Header ",trim(FileInput)
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        READ(iufile,*) Nb_atoms_aux, X0Cell_aux(1),  X0Cell_aux(2),  X0Cell_aux(3) !    60  -40.000000  -30.000000  -30.000000
        PRINT'(A25,3F12.6)', "    ...Origin [Bohr]   : ", X0Cell_aux(:)
        PRINT'(A43,I5)', "    ...Number of atoms read              : ", Nb_atoms_aux
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
      SUBROUTINE GaussianFileread(FileInput, iufile, dimensions, Nb_atoms_check, grid_check, X0Cell_check, cell_check, atoms_type_check, atoms_position_check, dV_check, VolCel_check, wavefunction, norm)

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
        ENDIF
        IF (  (ABS(atoms_position_aux(iatom,1)- atoms_position_check(iatom,1) ) + ABS(atoms_position_aux(iatom,2)- atoms_position_check(iatom,2))+ ABS(atoms_position_aux(iatom,3)- atoms_position_check(iatom,3))) > EPS_m5   ) THEN 
            PRINT*, "WARNING: Routine", TRIM(subroutinename)
            PRINT*, "WARNING: Difference in the positions of atom", iatom
            PRINT*, "Reference: ", atoms_position_check(iatom,:)
            PRINT*, "Read     : ", atoms_position_aux(iatom,:)
        ENDIF
     ENDDO
     PRINT*, "    ...READING Hole position"

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
        PRINT'(A35,ES13.6)', "    ...Norm:              ", norm
        PRINT'(A35,ES13.6)', "    ...Norm/Volume[Bohr^-3]      :", (norm/ Calculated_Volume)
        PRINT'(A35,ES13.6)', "    ...FFT*Norm/Volume[Bohr^-3]  :", (grid_aux(1)*grid_aux(2)*grid_aux(3)*norm/ Calculated_Volume)

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
 



   CHARACTER*(*), INTENT(inout)        :: FileInput
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
           ELSE IF ( atoms_type(iatom) == 0) THEN 
               WRITE(chr,'(A2)') 'X '
           ELSE 
               PRINT*, "Warning, not identified"
               WRITE(chr,'(A2)') 'X '
           ENDIF
           WRITE(iufile,'(A2,A3,3F12.6)')   chr,"   ", atoms_position(iatom,1), atoms_position(iatom,2), atoms_position(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146 
        ENDDO
        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE XYZFilewrite


