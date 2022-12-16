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
PROGRAM IMAGECHARGECORR

        IMPLICIT NONE

   CHARACTER*100, PARAMETER :: datafile="DATA.in"
   CHARACTER*20, PARAMETER :: subroutinename="CalculateImageCharge"

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


   INTEGER, PARAMETER :: FIRSTATOM=1
   INTEGER, PARAMETER :: Niteration=10000
   INTEGER, PARAMETER :: iufile=100
        REAL*8, PARAMETER :: Delta=0.04

         ! Input Parameters
         !Number of energies in the Siesta PDOS file
        CHARACTER*20   ::fname_in
        INTEGER :: BIAS=0
        INTEGER :: Nb_energy_in
        INTEGER :: norb
        INTEGER :: nat
        REAL*8  :: EDFTHOMO
        REAL*8  :: EDFTLUMO
        REAL*8  :: AuLeft
        REAL*8  :: AuRight
        REAL*8  :: Position_ImageCharge
        ! INPUT related VARIABLES        
        ! Useless things
        CHARACTER*100  :: chr
        CHARACTER*20   ::filename
        CHARACTER*1    :: chr11
        ! Important things
        ! Input energy grid
        REAL*8, ALLOCATABLE :: In_energygrid(:)
        ! Input PDOS
        REAL*8, ALLOCATABLE :: Input_orb_PDOS(:,:)
        ! Input Orbital center positions
        REAL*8, ALLOCATABLE :: Orb_center_coordinates(:,:)
        REAL*8, ALLOCATABLE :: Atoms_coordinates(:,:)
        INTEGER, ALLOCATABLE :: Orbtoatoms_coordinates(:)
        ! OUTPUT related VARIABLES       
   
         ! OUTPUT PDOS(z,Nb_energy_out)
        REAL*8, ALLOCATABLE :: PSEUDOCHARGEDENSITY(:,:) 
        REAL*8 :: IMAGECHARGE_HOMO, IMAGECHARGE_LUMO, DistanceSquared_xy, INTERPLANE_DISTANCE,PRODUCTDENSITY_HOMO,PRODUCTDENSITY_LUMO
        ! Test Variables
        ! LOOP Variables
        INTEGER :: ie, iat, iorb, k, iat2

        REAL*8 :: norm, minx, miny, minz, maxx, maxy,maxz, totminmax,tmp_real

        ! 
        ! End of Definitions 
        !

        ! ==================================================================================!
        ! ==================================================================================!
        ! ==================================================================================!

        ! 
        ! Main Body
        !

        open(iufile,file=trim(datafile),form='formatted')

        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) fname_in
        PRINT*, "   ...READING: FileName for PDOS: ", Trim(fname_in)
        READ(iufile,*) BIAS
        PRINT*, "   ...READING: BIAS: ", BIAS
        READ(iufile,*) Nb_energy_in
        PRINT*, "   ...READING: NB Energy: ",Nb_energy_in
        READ(iufile,*) norb
        PRINT*, "   ...READING: norb:",norb
        READ(iufile,*) nat
        PRINT*, "   ...READING: nat:", nat
        READ(iufile,*) EDFTHOMO
        PRINT*, "   ...READING: EDFTHOMO:",EDFTHOMO
        READ(iufile,*) EDFTLUMO
        PRINT*, "   ...READING: EDFTLUMO:",EDFTLUMO
        READ(iufile,*) AuLeft
        PRINT*, "   ...READING: Au Left[Ang]:" , AuLeft
        READ(iufile,*) AuRight 
        PRINT*, "   ...READING: AuRight[Ang]:",AuRight   
        READ(iufile,*) Position_ImageCharge
        PRINT*, "   ...READING: Position_ImageCharge[Ang]:", Position_ImageCharge
        CLOSE(iufile)
        PRINT*, "...Done"

      close(100)
        PRINT*, "Beginning..."

	IF (AuLeft>AuRight) THEN
		tmp_real=AuLeft
		AuLeft=AuRight
		AuRight=tmp_real
	ENDIF
        AuLeft=AuLeft*ANGSTROM_AU
        AuRight=AuRight*ANGSTROM_AU
        Position_ImageCharge=Position_ImageCharge*ANGSTROM_AU
        


        ! Allocation part       
        PRINT*, "Beginning Memory Allocations"
        ! Input grids
        ALLOCATE (In_energygrid(Nb_energy_in)  ) 
        ALLOCATE (Input_orb_PDOS(norb,Nb_energy_in))
        ! Output Grids
        ALLOCATE (PSEUDOCHARGEDENSITY(nat,2))
 
        ALLOCATE (Orb_center_coordinates(3,norb)) 
        ALLOCATE ( Orbtoatoms_coordinates(norb) )
        ALLOCATE ( Atoms_coordinates(3,nat) )
        PRINT*, "End of Allocations"
 


        ! Reading Input File

        PRINT*, "READING Incoming PDOS File"
      open(100,file=trim(fname_in),form='formatted')
      rewind(100)

        PRINT*, "opened"
      read(100,*)chr
      read(100,*)chr
      read(100,*)chr
      read(100,*)chr
      do ie = 1,Nb_energy_in
        read(100,*)In_energygrid(ie)
      end do
        PRINT*, "End of Energy grid Reading"

      read(100,*)chr
        PRINT*, "Orbital Loop"
      do iorb = 1,norb
       read(100,*)chr
       read(100,*)chr
 !       PRINT*, "Orbtoatom reading"

      read(100,"(a13,i3,a2)")chr,  Orbtoatoms_coordinates(iorb), chr11        
      print*,trim(chr),trim(chr11), Orbtoatoms_coordinates(iorb)
       read(100,*) chr
        PRINT*, "Orbtoatom read"
       read(100,"(a11,3f11.6,a2)")chr, Orb_center_coordinates(1,iorb), Orb_center_coordinates(2,iorb), &
     &  Orb_center_coordinates(3,iorb),chr11
       Atoms_coordinates(:,Orbtoatoms_coordinates(iorb))= Orb_center_coordinates(:,iorb)
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       read(100,*)chr
       do ie = 1, Nb_energy_in
           read(100,*)Input_orb_PDOS(iorb,ie)
       end do
       read(100,*)chr
       read(100,*)chr
      end do
      close(100)
        PRINT*, "END of READING"
    open(100,file="Positions_orbitals.dat",form='formatted')
          minx=Orb_center_coordinates(1,1)
          maxx=Orb_center_coordinates(1,1)
          miny=Orb_center_coordinates(2,1)
          maxy=Orb_center_coordinates(2,1)
          minz=Orb_center_coordinates(3,1)
          maxz=Orb_center_coordinates(3,1)
      do iorb = 2,norb
          IF (Orb_center_coordinates(1,iorb)< minx) minx= Orb_center_coordinates(1,iorb)
          IF (Orb_center_coordinates(2,iorb)< miny) miny= Orb_center_coordinates(2,iorb)
          IF (Orb_center_coordinates(3,iorb)< minz) minz= Orb_center_coordinates(3,iorb)

          IF (Orb_center_coordinates(1,iorb)> maxx)  maxx= Orb_center_coordinates(1,iorb)
          IF (Orb_center_coordinates(2,iorb)> maxy)  maxy= Orb_center_coordinates(2,iorb)
          IF (Orb_center_coordinates(3,iorb)> maxz)  maxz= Orb_center_coordinates(3,iorb)
      enddo 
    write(100,*)          minx, miny, minz, maxx, maxy, maxz
         totminmax=   (SUM(Orb_center_coordinates(1,:)) / norb        - ((minx + maxx)/2)) 
    write(100,*)  totminmax
               totminmax=   (SUM(Orb_center_coordinates(2,:)) / norb        - ((miny + maxy)/2)) 
    write(100,*)  totminmax
               totminmax=   (  SUM(Orb_center_coordinates(3,:)) / norb        - ((minz + maxz)/2) ) 
    write(100,*)  totminmax  
    do iorb = 1,norb
         write(100,*) (Orb_center_coordinates(1,iorb)-((maxx +minx)/2))  , (Orb_center_coordinates(2,iorb)-((maxy +miny)/2)), (Orb_center_coordinates(3,iorb)-((maxz +minz)/2))
      enddo
      close(100)
 


    ! OUTPUT PDOS calculation

    ! Orbital Loop
        PSEUDOCHARGEDENSITY(:,:)= 0.00000000000
	DO ie = 1, Nb_energy_in
		     IF ( ( in_energygrid(ie) <= EDFTHOMO + Delta ) .AND. ( in_energygrid(ie) >  EDFTHOMO - Delta ) ) THEN
                        DO iorb=1, norb
                           PSEUDOCHARGEDENSITY( Orbtoatoms_coordinates(iorb),1)= PSEUDOCHARGEDENSITY( Orbtoatoms_coordinates(iorb),1) + Input_orb_PDOS(iorb,ie)
                        ENDDO
                     ELSE IF ( ( in_energygrid(ie) <= EDFTLUMO + Delta ) .AND. ( in_energygrid(ie) >  EDFTLUMO - Delta ) ) THEN
                        DO iorb=1, norb
                           PSEUDOCHARGEDENSITY( Orbtoatoms_coordinates(iorb),2)= PSEUDOCHARGEDENSITY( Orbtoatoms_coordinates(iorb),2) + Input_orb_PDOS(iorb,ie)
                        ENDDO
                     ENDIF
        ENDDO ! ie



   !OUTPUT PART
    PRINT*, "PRINTING OUTPUT"
    
       OPEN ( 12, FILE="PDOS_perlevel.dat", FORM='formatted' )
       WRITE (12, '(a10)' ) "HOMO"
       norm=SUM( PSEUDOCHARGEDENSITY(:,1))
       PSEUDOCHARGEDENSITY(:,1)=PSEUDOCHARGEDENSITY(:,1)/norm
       WRITE (12, '(100(f15.9,a2))' ) (PSEUDOCHARGEDENSITY(iat,1)/norm , ", ", iat=1,nat)
       WRITE (12, '(2(e15.8))' )
       WRITE (12, '(2(f15.9))' ) norm, (SUM ( PSEUDOCHARGEDENSITY(:,1))/norm)
       WRITE (12, '(2(e15.8))' )
       WRITE (12, '(a10)' ) "LUMO"

       norm=SUM( PSEUDOCHARGEDENSITY(:,2))
       PSEUDOCHARGEDENSITY(:,2)=PSEUDOCHARGEDENSITY(:,2)/norm

       WRITE (12, '(100(f15.9,a2))' ) (PSEUDOCHARGEDENSITY(iat,2)/norm , ", ", iat=1,nat)
       WRITE (12, '(2(e15.8))' )
       WRITE (12, '(2(f15.9))' ) norm, (SUM ( PSEUDOCHARGEDENSITY(:,2))/norm)
       WRITE (12, '(2(e15.8))' )
      CLOSE(12)
    
       OPEN ( 12, FILE="PDOS_perlevelmathematica.dat", FORM='formatted' )
       WRITE (12, '(a4,i3,a1)' ) "num=",nat, ";"
       norm=SUM( PSEUDOCHARGEDENSITY(:,1))

       PSEUDOCHARGEDENSITY(:,1)=PSEUDOCHARGEDENSITY(:,1)/norm
       WRITE (12, '(a10)' ) "qHOMO = {"
       WRITE (12, '(100(f15.9,a2))' ) (PSEUDOCHARGEDENSITY(iat,1) , ", ", iat=1,nat-1), PSEUDOCHARGEDENSITY(nat,1), "};" 

       norm=SUM( PSEUDOCHARGEDENSITY(:,2))
       PSEUDOCHARGEDENSITY(:,2)=PSEUDOCHARGEDENSITY(:,2)/norm
       WRITE (12, '(a10)' ) "qLUMO = {"
       WRITE (12, '(100(f15.9,a2))' ) (PSEUDOCHARGEDENSITY(iat,2) , ", ", iat=1,nat-1), PSEUDOCHARGEDENSITY(nat,2), "};" 


       WRITE (12, '(a10)' ) "x = {"
       WRITE (12, '(100(f15.9,a2))' ) ((Atoms_coordinates(1,iat)*BOHR_RADIUS_ANGS), ", ", iat=1,nat-1), (Atoms_coordinates(1,nat)*BOHR_RADIUS_ANGS), "};" 
       WRITE (12, '(a10)' ) "y = {"
       WRITE (12, '(100(f15.9,a2))' ) ((Atoms_coordinates(2,iat)*BOHR_RADIUS_ANGS), ", ", iat=1,nat-1), (Atoms_coordinates(2,nat)*BOHR_RADIUS_ANGS), "};" 
       WRITE (12, '(a10)' ) "blist = {"
       WRITE (12, '(100(f15.9,a2))' ) ((Atoms_coordinates(3,iat)*BOHR_RADIUS_ANGS), ", ", iat=1,nat-1), (Atoms_coordinates(3,nat)*BOHR_RADIUS_ANGS), "};"
       WRITE (12, '(a5,f15.2,a1)' ) "AuL1=",AuLeft*BOHR_RADIUS_ANGS, ";"
       WRITE (12, '(a5,f15.2,a1)' ) "AuL2=",AuRight*BOHR_RADIUS_ANGS, ";"
       WRITE (12, '(a10,f15.2,a1)' ) "ImL1=AuL1+",Position_ImageCharge*BOHR_RADIUS_ANGS, ";"
       WRITE (12, '(a10,f15.2,a1)' ) "ImL2=AuL2-",Position_ImageCharge*BOHR_RADIUS_ANGS, ";"


      CLOSE(12)



	AuRight= AuRight - Position_ImageCharge ! ALREADY IN BOHR
	AuLeft= AuLeft + Position_ImageCharge ! ALREADY IN BOHR
	INTERPLANE_DISTANCE=ABS(AuRight-AuLeft)
	PRINT*, "INTERPLANE_DISTANCE [Ang]", INTERPLANE_DISTANCE*BOHR_RADIUS_ANGS
	PRINT*, "AuLeft",AuLeft
	Atoms_coordinates(3,:) = Atoms_coordinates(3,:) - AuLeft
	IMAGECHARGE_HOMO=ZERO
	IMAGECHARGE_LUMO=ZERO
 
  	DO     iat = 1, nat
		DO iat2=1, nat
			DistanceSquared_xy= (Atoms_coordinates(1,iat) - Atoms_coordinates(1,iat2))**2 + (Atoms_coordinates(2,iat) - Atoms_coordinates(2,iat2))**2 
			PRODUCTDENSITY_HOMO=-ONE*PSEUDOCHARGEDENSITY(iat2,1)*PSEUDOCHARGEDENSITY(iat,1) 
			PRODUCTDENSITY_LUMO=PSEUDOCHARGEDENSITY(iat2,2)*PSEUDOCHARGEDENSITY(iat,2)
			DO k=1,Niteration

				tmp_real=(ONE/ SQRT(DistanceSquared_xy + ( ( Atoms_coordinates(3,iat) - Atoms_coordinates(3,iat2) + (2.0*REAL(k)*INTERPLANE_DISTANCE))**2) ) ) - (ONE / SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2) - (2.0*REAL(k+1)*INTERPLANE_DISTANCE))**2) )) + (ONE / SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat2) - Atoms_coordinates(3,iat) + (2.0*REAL(k+1)*INTERPLANE_DISTANCE))**2) )) - (ONE / SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2) + (2.0*REAL(k)*INTERPLANE_DISTANCE))**2) ))
				IMAGECHARGE_HOMO=IMAGECHARGE_HOMO + (PRODUCTDENSITY_HOMO*tmp_real)
				IMAGECHARGE_LUMO=IMAGECHARGE_LUMO + (PRODUCTDENSITY_LUMO*tmp_real)

			ENDDO 
			tmp_real=  (-ONE/(SQRT(DistanceSquared_xy +  ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2) - (2.0*INTERPLANE_DISTANCE) )**2) ) ) )  + (ONE /(SQRT(DistanceSquared_xy + ((Atoms_coordinates(3,iat) - Atoms_coordinates(3,iat2) + (2.0*INTERPLANE_DISTANCE) )**2) )))- (ONE /(SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2)  )**2) ))) 
 
			IMAGECHARGE_HOMO=IMAGECHARGE_HOMO + (PRODUCTDENSITY_HOMO*tmp_real)
			IMAGECHARGE_LUMO=IMAGECHARGE_LUMO + (PRODUCTDENSITY_LUMO*tmp_real)


		ENDDO
	ENDDO 

	!The reason to multiply by Ry rather than Ha (a.u.) is because there is another factor of 2
	!due to image charge model. for simplicity, absorb that 1/2 to energy conversion factor.

	IMAGECHARGE_HOMO=IMAGECHARGE_HOMO* rytoev
	IMAGECHARGE_LUMO=IMAGECHARGE_LUMO* rytoev
	PRINT*, "Image Charge HOMO Iteration 10000", IMAGECHARGE_HOMO
	PRINT*, "Image Charge LUMO Iteration 10000", IMAGECHARGE_LUMO


	IMAGECHARGE_HOMO=ZERO
	IMAGECHARGE_LUMO=ZERO
 
	DO     iat = 1, nat
		DO iat2=1, nat
			DistanceSquared_xy= (Atoms_coordinates(1,iat) - Atoms_coordinates(1,iat2))**2 + (Atoms_coordinates(2,iat) - Atoms_coordinates(2,iat2))**2 
			PRODUCTDENSITY_HOMO=-ONE*PSEUDOCHARGEDENSITY(iat2,1)*PSEUDOCHARGEDENSITY(iat,1) 
			PRODUCTDENSITY_LUMO=PSEUDOCHARGEDENSITY(iat2,2)*PSEUDOCHARGEDENSITY(iat,2)
			DO k=1,(2*Niteration)

				tmp_real=(ONE/ SQRT(DistanceSquared_xy + ( ( Atoms_coordinates(3,iat) - Atoms_coordinates(3,iat2) + (2.0*REAL(k)*INTERPLANE_DISTANCE))**2) ) ) - (ONE / SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2) - (2.0*REAL(k+1)*INTERPLANE_DISTANCE))**2) )) + (ONE / SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat2) - Atoms_coordinates(3,iat) + (2.0*REAL(k+1)*INTERPLANE_DISTANCE))**2) )) - (ONE / SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2) + (2.0*REAL(k)*INTERPLANE_DISTANCE))**2) ))
				IMAGECHARGE_HOMO=IMAGECHARGE_HOMO + (PRODUCTDENSITY_HOMO*tmp_real)
				IMAGECHARGE_LUMO=IMAGECHARGE_LUMO + (PRODUCTDENSITY_LUMO*tmp_real)

			ENDDO 
			tmp_real=  (-ONE/(SQRT(DistanceSquared_xy +  ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2) - (2.0*INTERPLANE_DISTANCE) )**2) ) ) )  + (ONE /(SQRT(DistanceSquared_xy + ((Atoms_coordinates(3,iat) - Atoms_coordinates(3,iat2) + (2.0*INTERPLANE_DISTANCE) )**2) )))- (ONE /(SQRT( DistanceSquared_xy + ((Atoms_coordinates(3,iat) + Atoms_coordinates(3,iat2)  )**2) ))) 
 
			IMAGECHARGE_HOMO=IMAGECHARGE_HOMO + (PRODUCTDENSITY_HOMO*tmp_real)
			IMAGECHARGE_LUMO=IMAGECHARGE_LUMO + (PRODUCTDENSITY_LUMO*tmp_real)


		ENDDO
	ENDDO 


	!The reason to multiply by Ry rather than Ha (a.u.) is because there is another factor of 2
	!due to image charge model. for simplicity, absorb that 1/2 to energy conversion factor.

	IMAGECHARGE_HOMO=IMAGECHARGE_HOMO* rytoev
	IMAGECHARGE_LUMO=IMAGECHARGE_LUMO* rytoev
	PRINT*, "Image Charge HOMO Iteration 10000", IMAGECHARGE_HOMO
	PRINT*, "Image Charge LUMO Iteration 10000", IMAGECHARGE_LUMO







	!
	       WRITE(filename, '(I5,A12)') BIAS , ".IC.shift_mo"
       OPEN ( 12, FILE=TRIM(filename), FORM='formatted' )
       WRITE (12, '(A15)' ) "%block Shift_MO"
       WRITE (12, '(2I3,A50)' ) FIRSTATOM, (FIRSTATOM+nat-1),      "# (integers) molecule is from iatom1 to iatom2"
       WRITE (12, '(A70)' )  "  HOMOLEVEL  # (integer) index of homo level without counting spin degeneracy"
       WRITE (12, '(A16,F13.3,A5,F13.3,A18)' )  "TOTALHOMO #PBEH=", EDFTHOMO,  " ICH=", (IMAGECHARGE_HOMO), " [eV] Linear Term="
       WRITE (12, '(A16,F13.3,A5,F13.3,A18)' )  "TOTALLUMO #PBEL=", EDFTLUMO,  " ICL=", (IMAGECHARGE_LUMO), " [eV] Linear Term="
       WRITE (12, '(A18)' ) "%endblock Shift_MO"
       CLOSE(12)



     ! De-Allocation part       


        PRINT*, "Deallocating"
        ! Input grids
        DEALLOCATE (In_energygrid) 
        DEALLOCATE (Input_orb_PDOS)
        ! Output Grids
        DEALLOCATE (PSEUDOCHARGEDENSITY) 
        DEALLOCATE (Orb_center_coordinates) 
        DEALLOCATE (Atoms_coordinates)
        PRINT*, "End of Deallocations"

        PRINT*, "End of Program"
END PROGRAM IMAGECHARGECORR



