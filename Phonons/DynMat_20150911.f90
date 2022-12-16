module constants
IMPLICIT NONE
! Debug
LOGICAL, parameter :: debug = .TRUE.
! Kinds
integer, parameter :: dp = kind(1.0d0)
integer, parameter :: dpc = kind((1.0d0,1.0d0))
! Physical constants/conversion
Integer, parameter :: ndim = 3
Real(dp), parameter :: BOHR_ANG = 0.529177D0 
Real(dp), parameter :: BOHR_M = 0.529177D-10 
Real(dp), parameter :: Forces_AU_SI = 8.2387225D-8 ! Hartree/bohr to N 
Real(dp), parameter :: ANG_BOHR = 1.889727D0 
Real(dp), parameter :: THz_CM = 33.35641D0
Real(dp), parameter :: CM_meV = 0.12398D0
Real(dp), parameter :: THz_meV = 4.13567D0
REAL(dp), PARAMETER :: RYD_eV = 13.605826d0 ! eV 
REAL(dp), PARAMETER :: AU_THz  = 2.418D-5          ! THz
REAL(dp), PARAMETER :: AMU_KG = 1.660538921D-27 ! Atomic mass
REAL(dp), PARAMETER :: AMU_KGmodif= 1.660538921D-3 ! AMU to Thz
! Numbers
Real(dp), parameter :: PI = 3.14159265358979323846D0
REAL(dp), PARAMETER ::    ZERO = 0.0D0
REAL(dp), PARAMETER ::     ONE = 1.0D0
REAL(dp), PARAMETER ::     TWO = 2.0D0
! Tolerance
REAL(dp), parameter :: tolerance=0.00001D0
! Files assignments
Integer, parameter :: icor = 10
Integer, parameter :: icutm = 11
Integer, parameter :: icutp = 12
Integer, parameter :: istatat = 13
Integer, parameter :: iout = 30
Integer, parameter :: imodem = 37
Integer, parameter :: imodesl = 39
Integer, parameter :: ifixQm = 41 ! Conf with target amplitude
Integer, parameter :: ifixQp = 42 ! Conf with target amplitude
Integer, parameter :: ifixDm = 43 ! Conf with target mean displacement
Integer, parameter :: ifixDp = 44 ! Conf with target mean displacement.
Integer, parameter :: imod = 45
end module constants

module type_mod
use constants
PUBLIC :: atom, atomforce, aname
implicit none
! Types
TYPE :: atom
  CHARACTER(15) :: name    ! name of object
  REAL :: mass 
  REAL, Dimension(ndim) :: r ! coordinates
  !REAL, Dimension(ndim) :: f ! forces
  !INTEGER ::  status = 0 ! this variable used to store the indice of the type.  
  INTEGER :: type = 0    
END TYPE atom
!
TYPE :: atomforce
  REAL :: mass 
  REAL, Dimension(ndim) :: f ! forces
END TYPE atomforce
!
TYPE :: aname
 CHARACTER(15) :: name    ! name of atom
END TYPE aname

end module type_mod
!***********************************************
!
! Copyright (C) 2007-2008 UT Austin
!               2008-2012 Molecular Foundry
!               2012- Bowling Green State University
!               2015- Argonne Nationa l Lb 
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Alexey ZAYAK, Pierre DARANCET
!    Alexey ZAYAK: azayak@bgsu.edu
!    Pierre DARANCET: pdarancet@anl.gov   
!***********************************************
! Description
!program for creating and diagonalizing force constant matrix
!it needs the "forces_minus" and "forces_plus" file which 
!contain forces for all negative and positive
!displacements ..... 
!***********************************************
Program vwhole
use constants
use type_mod, only : atom, atomforce, aname
implicit none
! Atoms Descriptors
TYPE(atom), dimension(:), allocatable :: Atoms, StaticAt
TYPE(atomforce), dimension(:), allocatable :: ForceMatrix
TYPE(aname), dimension(:), allocatable :: AName ! Atomic Symbol
INTEGER, dimension(:), allocatable :: AType  ! Siesta type
Integer, dimension(:), allocatable ::  Natom ! Natom per type 
! Loop variables
Integer :: Write_counter, Stcounter
Integer :: Mode, idisp
Integer :: alpha, beta, ai, bi ! loop over Dm
Integer :: i, j ! loop over force constant
Integer :: x,y !loop over A
Integer :: t ! type
! I/O formatting
CHARACTER(30) :: STR_arg, F_name, M_name
CHARACTER(30) :: POS_name, filename
CHARACTER(30) :: File_n1, File_n2, File_n3
! Conversions Factors 
REAL(dp) :: A1,  B1,  Conv
REAL(dp) :: m1, m2, joint_mass
! Temporary calc. variables
REAL(dp) :: tmpreal
! Functions
!REAL(dp) :: whatismymass
! Reading variables
Integer :: readtmp
Integer :: inttmp
CHARACTER(30) :: chartmp
Real, Dimension(ndim) :: cor 
Integer :: N_arg
! LAPACK feeding variables
Integer :: LWork 
Integer, dimension(:), allocatable :: pivot
Real(dp), dimension(:), allocatable :: Work
Real(dp), dimension(:), allocatable :: c
Real(dp), dimension(:,:), allocatable ::  fB
Integer :: Info
! Output
Real(dp), dimension(:), allocatable :: displacement 
Real(dp), dimension(:), allocatable :: disp_spread, amp_spread 
Real(dp), dimension(:,:), allocatable ::  A ! the matrix
Real, dimension(:), allocatable :: mass
Real(dp), dimension(:), allocatable :: eigv

! Input: 
Integer :: N  !size of the dynamical matrix
Real(dp) :: Disp ! read in Ang, converted to Bohr
Integer :: NatTotal, StaticAtN, NATypes, Ndisp
CHARACTER(15) :: OutputMode ! Disp or Amp --define output type
Real(dp)  :: ScaleAmp  !read in Ang; Target Q in Bohr * sqrt(m)
Real(dp)  :: ScaleDisp   !read in Bohr; Target Disp in Bohr
! Matrix dimensions
Integer ::  Nat


! ==================================================================
! ==================================================================
! ==========================  START  ===============================
! ==================================================================
! ==================================================================
!
! reading arguments
!

open(UNIT=iout, file='out.dat',form='formatted', status='unknown')

IF ( debug ) THEN
      WRITE(*,*) "Reading arguments"
      WRITE(iout,*)  "Reading arguments"
END IF
N_arg = IARGC()          ! reads the number of arguments

if (N_arg .ne. 8) then
   write(*,*) '------------------------------------'
   write(*,*) 'This program needs 8 arguments'
   write(*,*) N_arg, 'arguments were entered'
   write(*,*) '1) Number of species'
   write(*,*) '2) Total number of atoms'
   write(*,*) '3) Number of displaced atoms'
   write(*,*) '4) Number of calculated configurations'
   write(*,*) '5) Dimension of Dyn Mat'
   write(*,*) '6) Input Displacement (Ang) '
   write(*,*) '7) Output Amplitude (Amp) OR Displacement (Disp)'
   write(*,*) '8) Output Amplitude/Displacement (au*Ang/Ang)'
   write(*,*) '------------------------------------'
   write(*,*) 'Input'
   write(*,*) '*.exe reads coordinates and forces'
   write(*,*) '- from cut_minus and cut_plus'
   write(*,*) '- from static.atoms'
   write(*,*) '- from siesta.cor'
   write(*,*) '------------------------------------'
   write(*,*) 'Have fun!'
   write(*,*) '------------------------------------'
   stop
end if


CALL GETARG(1, STR_arg) ! 1) Number of species 
 read(STR_arg,*) NATypes

CALL GETARG(2, STR_arg) ! Total number of atoms
 read(STR_arg,*) NatTotal

CALL GETARG(3, STR_arg) ! 3) Number of displaced atoms
 read(STR_arg,*) Nat

CALL GETARG(4, STR_arg) ! 4) Number of calculated configurations
 read(STR_arg,*) Ndisp 

CALL GETARG(5, STR_arg) ! 5) Dimension of Dyn Mat
 read(STR_arg,*) N
 LWork = (64+2)*N  ! Required by LAPACK

CALL GETARG(6, STR_arg) ! 6) Input Displacement (Ang) 
 read(STR_arg,*) Disp
 Disp = Disp * ANG_BOHR ! in Bohr

CALL GETARG(7, STR_arg) ! 7) Output mode: 
 read(STR_arg,*) OutputMode
 IF ( TRIM(OutputMode) == 'Disp' ) THEN
    IF ( debug ) THEN
        WRITE(*,*) "Reading arguments: Output mode Disp"
        WRITE(iout,*)  "Reading arguments: Output mode Disp"
    END IF
 ELSE IF ( TRIM(OutputMode) == 'Amp' ) THEN
    IF ( debug ) THEN
        WRITE(*,*) "Reading arguments: Output mode Amp"
        WRITE(iout,*)  "Reading arguments: Output mode Amp"
    END IF
 ELSE IF ( TRIM(OutputMode) == 'RescaledDisp' ) THEN 
    IF ( debug ) THEN
        WRITE(*,*) "Reading arguments: Output mode RescaledDisp"
        WRITE(iout,*) "Reading arguments: Output mode scaledDisp"
    END IF
 ELSE IF ( TRIM(OutputMode) == 'RescaledAmp' ) THEN
    IF ( debug ) THEN
        WRITE(*,*) "Reading arguments: Output mode scaledAmp"
        WRITE(iout,*) "Reading arguments: Output mode scaledAmp"
    END IF
 ELSE
     WRITE(*,*) 'Error: OutputMode:', TRIM(OutputMode)
     WRITE(*,*) 'Not understood'
     WRITE(iout,*) 'Error: OutputMode:', TRIM(OutputMode)
     WRITE(iout,*) 'Not understood'
     STOP
 END IF

CALL GETARG(8, STR_arg) ! 8)Amplitude / Displacement of Modes (au*Ang/Ang)
 read(STR_arg,*) ScaleDisp
 ScaleDisp = ScaleDisp * ANG_BOHR ! in Bohr
 ScaleAmp = ScaleDisp * ANG_BOHR ! in Bohr
IF ( debug ) THEN
    WRITE(*,*) "Done: Reading arguments"
    WRITE(iout,*)  "Done: Reading arguments"
END IF

StaticAtN = NatTotal - Nat

!
! DONE: reading arguments
!

!
! ALLOCATING
!

IF ( debug ) THEN
      WRITE(*,*) "Allocating arrays"
      WRITE(iout,*)  "Allocating arrays"
END IF
allocate(AName(NATypes+1)) ! array of names
allocate(AType(NATypes+1)) ! +1 is for static atoms
allocate(Natom(NATypes+1))
allocate(Atoms(NatTotal))  ! allocate array for atoms
allocate(ForceMatrix(Nat*(Ndisp+1)))  ! allocate array for atoms
!allocate(StaticAt(StaticAtN))
ALLOCATE(Work(LWork))
ALLOCATE(A(N,N))
ALLOCATE(fB(N,N))
ALLOCATE(c(N)) 
ALLOCATE(eigv(N))
ALLOCATE(displacement(N))
ALLOCATE(disp_spread(N))
ALLOCATE(amp_spread(N))
ALLOCATE(pivot(N))
allocate(mass(NATypes+1))
IF ( debug ) THEN
    WRITE(*,*) "DOne: Allocating arrays"
    WRITE(iout,*)  "Done: Allocating arrays"
END IF
!
! DONE: ALLOCATING
!
!
!============== READING SIESTA.COR  ===============
!

filename='siesta.cor' ! to be generalized
CALL ReadSiestaCor(icor,TRIM(filename),AName(1:NATypes),Natom(1:NATypes),mass(1:NATypes),NATypes)

!
AName(NATypes+1)='Au' !static atoms to be generalized
Natom(NATypes+1)=StaticAtN ! static atoms to be generalized
mass(NATypes+1)=1000.0 ! static atoms to be generalized

!
filename='cut_plus' ! to be generalized
filename2='static.atoms' ! to be generalized
CALL ReadAtoms(icutp,istatat,TRIM(filename),TRIM(filename2),Atoms(:),AName(:),Natom(:),mass(:),NATypes,Nat, StaticAtN)
!
!
filename='cut_plus' ! to be generalized
filename2='cut_minus' ! to be generalized
CALL ReadForces(icutm,icutp,TRIM(filename),TRIM(filename2),ForceMatrix,Atoms(1:Nat),Natom(1:NATypes),mass(1:NATypes),Disp,NATypes,Nat, Ndisp)

!---------------symmetrize/scale/prepare for diag
CALL WriteDynMat(iout,ForceMatrix,A,Disp,Nat, Ndisp,N)

! ------------------------------------

!Diagonalize this matrix
! to find all the eigenvalues and eigenvectors of the symmetric matrix
!'N':  only igenvalues, 'V' iegenvectors too
!'U':  Upper triangle of A is stored
! N : order of the matrix
! A : the matrix
! Work : workspace array
! LWork : size of Work
! Info : =0 if successful
IF ( debug ) THEN
   write(*,*) "Diagonalizing..."
   write(iout,*) "Diagonalizing..."
ENDIF
call DSYEV('V','U',N,A,N,c,Work,LWork,Info)
IF ( debug ) THEN
   write(*,*) " ... Done"
   write(iout,*) " ... Done"
ENDIF
! Diagonalization done

! Computing displacements and scaled displacements/amp

IF ( debug ) THEN
   write(*,*) " Computing freq and disp from eigen"
   write(iout,*) " Computing freq and disp from eigen"
ENDIF
write(iout,'(a70)') '#Mode Freq(cm-1) Freq(meV) u(A) Mass(amu), SprAmp (sqr(m)A), SprU(A)'
write(*,'(a70)') '#Mode Freq(cm-1) Freq(meV) u(A) Mass(amu), SprAmp (sqr(m)A), SprU(A)'
do Mode=1, N 
    if (c(Mode) .lt. ZERO ) then 
        eigv(Mode)=-sqrt(abs(c(Mode)))
    else
        eigv(Mode)=sqrt(c(Mode))
    endif
    ! Displacements  
    displacement(Mode)=ZERO
    ! Spread in displacement
    disp_spread(Mode)=ZERO
    ! Spread in Amplitude
    amp_spread(Mode)=ZERO
    do t=1, Nat
        ! Amplitude^2 on atom t 
        tmpreal=ZERO
        do j=1,ndim
            tmpreal=tmpreal+((A((j+(ndim*(t-1))),Mode))**2)
        end do
        ! Amplitude Spread = sum_t sqrt(Amp(t)^2)
        amp_spread(Mode)=amp_spread(Mode)+sqrt(tmpreal)
        ! Displacement Spread = sum_t sqrt((Amp(t)^2)/Mass(t))
        disp_spread(Mode)=disp_spread(Mode)+sqrt(tmpreal/Atoms(t)%mass)
        ! Displacement = sqrt( sum_t (Amp(t)^2)/Mass(t))
        ! sqrt(displacement) is taken at the end of the atom loop
        displacement(Mode)=displacement(Mode)+(tmpreal/Atoms(t)%mass)
    end do
    displacement(Mode)=BOHR_ANG*sqrt( displacement(Mode) )
    write(iout,'(7f15.8,i4)') ONE, &
           eigv(Mode)*THz_CM/TWO/PI,&
           eigv(Mode)*THz_meV/TWO/PI,& 
           displacement(Mode),&
           (BOHR_ANG/displacement(Mode))**2,&
           amp_spread(Mode), &
           disp_spread(Mode), &
           Mode
    write(*,'(7f15.8,i4)') ONE, &
           eigv(Mode)*THz_CM/TWO/PI,&
           eigv(Mode)*THz_meV/TWO/PI,& 
           displacement(Mode),&
           (BOHR_ANG/displacement(Mode))**2,&
           amp_spread(Mode), &
           disp_spread(Mode), &
           Mode
end do ! Mode loop
IF ( debug ) THEN
   write(*,*) " ... Done"
   write(iout,*) " ... Done"
ENDIF

!===================================
!
!   We now have everything we need, the rest is output
!
!===================================
IF ( TRIM(output_type) == 'nmd' ) THEN
    IF ( debug ) THEN
       write(*,*) " Writing mod.nmd"
       write(iout,*) " Writing mod.nmd"
    ENDIF
    filename="mod.nmd"
    CALL write_NMD_File(imod,TRIM(filename),Atoms(1:Nat),StaticAt,A,Nat,StaticAtN,ndim)
    IF ( debug ) THEN
        write(*,*) " Done writing mod.nmd"
        write(iout,*) " Done writing mod.nmd"
    ENDIF
ELSE 
    !===================================
    IF ( debug ) THEN
       write(*,*) " Writing Mode-specific output"
       write(iout,*) " Writing Mode-specific output"
    ENDIF
    !====================================
    
    do Mode=1, Nat*ndim, 1 ! Modeloop

write(M_name,'(i10)') Mode 
write(F_name,'(f12.2)') eigv(Mode)*THz_CM/TWO/PI
filename=TRIM(ADJUSTL(M_name))//'_'//TRIM(ADJUSTL(F_name))//'.'//TRIM(ADJUSTL(output_type))
!CALL Print_Conf_File
!open(UNIT=imodem,file=File_n1,form='formatted', status='unknown')
!          File_n2='froz_in.'// TRIM(ADJUSTL(M_name))
!          open(UNIT=38,file=File_n2,form='formatted',status='unknown')
!File_n3='sl_'// TRIM(ADJUSTL(M_name))// '_' //TRIM(ADJUSTL(F_name))// '.xyz'
!open(UNIT=imodesl,file=File_n3,form='formatted',status='unknown')
! Outpuit Files

IF ( TRIM(output_mode) == 'ScaleNormAmp' ) THEN
CALL writeoutmode(imodem,filename,TRIM(output_type),Atoms(1:Nat),StaticAt,A(:,Mode),Nat,StaticAtN,ndim)
ELSE IF  ( TRIM(output_mode) == 'ScaleNormDisp' ) THEN  

ELSE IF  ( TRIM(output_mode) == 'ScaleDisp' ) THEN  
ELSE 
ENDIF

POS_name = 'ConfFixAmp_m' // TRIM(ADJUSTL(M_name))
open (UNIT=ifixQm, FILE=POS_name, STATUS='unknown',form='formatted',action="write")
POS_name = 'ConfFixAmp_p' // TRIM(ADJUSTL(M_name))
open (UNIT=ifixQp, FILE=POS_name, STATUS='unknown',form='formatted',action="write")
POS_name = 'ConfFixDisp_m' // TRIM(ADJUSTL(M_name))
open (UNIT=ifixDm, FILE=POS_name, STATUS='unknown',form='formatted',action="write")
POS_name = 'ConfFixDisp_p' // TRIM(ADJUSTL(M_name))
open (UNIT=ifixDp, FILE=POS_name, STATUS='unknown',form='formatted',action="write")
write(imodem,*) Nat
write(imodem,*)
do i=1, Nat, 1 
    write(imodem,'(a2,6f10.4)') Atoms(i)%name, &
        Atoms(i)%r(:)*BOHR_ANG,&
        A((1+ndim*(i-1)):(ndim*i),Mode)*BOHR_ANG*&
        ScaleAmp/sqrt(Atoms(i)%mass)
end do
write(imodesl,*) (Nat+StaticAtN)
write(imodesl,*)
Write_counter = 0
do i=1, (Nat+StaticAtN), 1
  if (i .le. Nat) then
     Write_counter = Write_counter + 1 
     write(imodesl,'(a2,6f14.8)') Atoms(Write_counter)%name, &
        Atoms(Write_counter)%r(:)*BOHR_ANG,&
        A((1+ndim*(Write_counter-1)):(ndim*Write_counter),Mode)*&
        BOHR_ANG*ScaleAmp/sqrt(Atoms(Write_counter)%mass)
  else
     Write_counter = Write_counter + 1 ! to be generalized to non
     ! au leads
     write(imodesl,'(a2,6f14.8)') 'Au', &
         StaticAt(Write_counter-Nat)%r(:),&
         0.0,  0.0,  0.0           
  endif   
end do
! writing to Conf file_m
!---------------------------------------------------------------
  !Fix Q
  Write_counter = 0
  do i=1, (Nat+StaticAtN), 1

  if (i .le. Nat) then
    Write_counter = Write_counter + 1 
      write(ifixQm,'(3f14.8, i4,"  ",a2,"  ",i4)')  &
          (Atoms(Write_counter)%r(1)*BOHR_ANG-&
          ScaleAmp*A((1+ndim*(Write_counter-1)),Mode)*&
          BOHR_ANG/sqrt(Atoms(Write_counter)%mass)),&
          (Atoms(Write_counter)%r(2)*BOHR_ANG-&
          ScaleAmp*A((2+ndim*(Write_counter-1)),Mode)*&
          BOHR_ANG/sqrt(Atoms(Write_counter)%mass)),&
          (Atoms(Write_counter)%r(3)*BOHR_ANG-&
          ScaleAmp*A((3+ndim*(Write_counter-1)),Mode)*&
          BOHR_ANG/sqrt(Atoms(Write_counter)%mass)),&
          Atoms(Write_counter)%type,&
          Atoms(Write_counter)%name,&
          Write_counter
   else
      Write_counter = Write_counter + 1
      write(ifixQm,'(3f14.8, i4,"  ",a2,"  ",i4)') & 
          StaticAt(Write_counter-Nat)%r(1),&
          StaticAt(Write_counter-Nat)%r(2),&
          StaticAt(Write_counter-Nat)%r(3),&
          4,&
          "Au",&
          Write_counter
      ! to be generalized to non Au leads
   endif   
  end do
  !Fix Disp
  Write_counter = 0
  do i=1, (Nat+StaticAtN), 1
    if (i .le. Nat) then
      Write_counter = Write_counter + 1 
      write(ifixDm,'(3f14.8, i4,"  ",a2,"  ",i4)') &
        (Atoms(Write_counter)%r(1)*BOHR_ANG-&
        ScaleDisp*A((1+ndim*(Write_counter-1)),Mode)*&
        BOHR_ANG/(sqrt(Atoms(Write_counter)%mass)*&
        displacement(Mode))),&
        (Atoms(Write_counter)%r(2)*BOHR_ANG-&
        ScaleDisp*A((2+ndim*(Write_counter-1)),Mode)*&
        BOHR_ANG/(sqrt(Atoms(Write_counter)%mass)*&
        displacement(Mode))),&
        (Atoms(Write_counter)%r(3)*BOHR_ANG-&
        ScaleDisp*A((3+ndim*(Write_counter-1)),Mode)*&
        BOHR_ANG/(sqrt(Atoms(Write_counter)%mass)*&
        displacement(Mode))),&
        Atoms(Write_counter)%type,&
        Atoms(Write_counter)%name,&
        Write_counter
    else
      Write_counter = Write_counter + 1
      write(ifixDm,'(3f14.8, i4,"  ",a2,"  ",i4)')& 
         StaticAt(Write_counter-Nat)%r(1),&
         StaticAt(Write_counter-Nat)%r(2),&
         StaticAt(Write_counter-Nat)%r(3),&
         4,&
         "Au",&
         Write_counter
      ! Tp be generalized to non-Au leads
    endif   
  end do
! writing to Conf_p
!---------------------------------------------------------------
  !Fix Q
  Write_counter = 0
  do i=1, (Nat+StaticAtN), 1
    if (i .le. Nat) then
     Write_counter = Write_counter + 1 
       write(ifixQp,'(3f14.8, i4,"  ",a2,"  ",i4)') & 
           (Atoms(Write_counter)%r(1)*BOHR_ANG+&
           ScaleAmp*A((1+ndim*(Write_counter-1)),Mode)*&
           BOHR_ANG/sqrt(Atoms(Write_counter)%mass)),&
           (Atoms(Write_counter)%r(2)*BOHR_ANG+&
           ScaleAmp*A((2+ndim*(Write_counter-1)),Mode)*&
           BOHR_ANG/sqrt(Atoms(Write_counter)%mass)),&
           (Atoms(Write_counter)%r(3)*BOHR_ANG+&
           ScaleAmp*A((3+ndim*(Write_counter-1)),Mode)*&
           BOHR_ANG/sqrt(Atoms(Write_counter)%mass)),&
           Atoms(Write_counter)%type,&
           Atoms(Write_counter)%name,&
           Write_counter
    else
       Write_counter = Write_counter + 1
       write(ifixQp,'(3f14.8, i4,"  ",a2,"  ",i4)') & 
           StaticAt(Write_counter-Nat)%r(1),&
           StaticAt(Write_counter-Nat)%r(2),&
           StaticAt(Write_counter-Nat)%r(3),&
           4,&
           "Au",&
           Write_counter
       ! to be generalized to non-Au leads
    endif   
  end do
  !Fix Disp
  Write_counter = 0
  do i=1, (Nat+StaticAtN), 1
    if (i .le. Nat) then
      Write_counter = Write_counter + 1 
      write(ifixDp,'(3f14.8, i4,"  ",a2,"  ",i4)')& 
          (Atoms(Write_counter)%r(1)*BOHR_ANG+&
          ScaleDisp*A((1+ndim*(Write_counter-1)),Mode)*&
          BOHR_ANG/(sqrt(Atoms(Write_counter)%mass)*&
          displacement(Mode))),&
          (Atoms(Write_counter)%r(2)*BOHR_ANG+&
          ScaleDisp*A((2+ndim*(Write_counter-1)),Mode)*&
          BOHR_ANG/(sqrt(Atoms(Write_counter)%mass)* &
          displacement(Mode))),&
          (Atoms(Write_counter)%r(3)*BOHR_ANG+&
          ScaleDisp*A((3+ndim*(Write_counter-1)),Mode)*&
          BOHR_ANG/(sqrt(Atoms(Write_counter)%mass)*&
          displacement(Mode))),&
          Atoms(Write_counter)%type,&
          Atoms(Write_counter)%name,&
          Write_counter
    else
      Write_counter = Write_counter + 1
      write(ifixDp,'(3f14.8, i4,"  ",a2,"  ",i4)')&
          StaticAt(Write_counter-Nat)%r(1),&
          StaticAt(Write_counter-Nat)%r(2),&
          StaticAt(Write_counter-Nat)%r(3),&
          4,&
          "Au",&
          Write_counter
      ! to be generalized to non Au leads
    endif   
  end do
!--------------------------------------------------------------
        IF ( debug ) THEN
           write(*,*) " ... Done"
           write(iout,*) " ... Done"
        ENDIF
 
        write(imodem,*)'#frequency cm-1 ', eigv(Mode)*THz_CM/TWO/PI 
        close(ifixQp)
        close(ifixQm)
        close(ifixDp)
        close(ifixDm)
        close(imodem)
        close(imodesl)
!        close(38)
enddo
ENDIF
    !============== Deallocating ===============
deallocate(AName) ! array of names
deallocate(AType)
deallocate(Natom)
deallocate(Atoms) ! allocate array for atoms 
deallocate(StaticAt)
DEALLOCATE(Work)
DEALLOCATE(A)
DEALLOCATE(fB)
DEALLOCATE(c) 
DEALLOCATE(eigv)
DEALLOCATE(disp_spread)
DEALLOCATE(amp_spread)
DEALLOCATE(displacement)
DEALLOCATE(pivot)
deallocate(mass)
!==============/ Deallocating ===============

!..........................................

! 10     format(a2, 3f12.8)
! 11     format(i4, '  ', a2, 3f12.8, i3)
! 12     format(a16, a2, 3f12.8)
close(iout)
! 190    write(*,*) 'error when opening file'
!    STOP
WRITE(*,*) "END"

        CONTAINS
        
SUBROUTINE  ReadAtoms(icutp,istatat,filename,filename2,Atoms,AName,Natom,mass,NATypes,Nat, StaticAtN)
use constants, ONLY : debug, dp, ndim
use type_mod, ONLY  atom, aname 
IMPLICIT NONE
INTEGER, INTENT(in) :: icutp
INTEGER, INTENT(in) :: istatat
INTEGER, INTENT(in) :: Nat
INTEGER, INTENT(in) :: StaticAtN

CHARACTER(30), INTENT(in) :: filename
CHARACTER(30), INTENT(in) :: filename2
TYPE(aname), INTENT(in) :: AName(NATypes+1)
INTEGER, INTENT(in) :: Natom(NATypes+1)
REAL, INTENT(in) :: mass(NATypes+1)
TYPE(atom), INTENT(out) :: Atoms(Nat+StaticAtN)

INTEGER :: idisp, t, i
REAL(dp) :: cor(ndim), force(ndim)
Integer :: readtmp, Acounter
!
!---------------READING FORCE / Coordinates input files
!
IF ( debug ) THEN
    write(*,*) "Opening  cutp file"
ENDIF

open (UNIT=icutp,FILE=TRIM(filename),STATUS='OLD',action="read")
REWIND(icutp)
!----------------------------------------------------
! reading the zero displacement configuration ...
! the first entry
idisp=0
!..............Reading config number
read(icutp,*) readtmp
IF ( readtmp .ne. idisp  ) THEN
   PRINT*, 'ERROR!!!' 
   PRINT*, "Config p:", readtmp, "instead of", idisp
   STOP
ENDIF
!.............../Reading config number
Acounter = 0 ! Atom index
do t=1, NATypes, 1
  do i=1, Natom(t), 1
    Acounter = Acounter + 1
    read(icutp,*) cor(:),force(:)
    !Verifies that 0 config are the same for plus and minus
    !Atoms(Acounter)%status=t
    Atoms(Acounter)%r(:)=cor(:)*ANG_BOHR
    Atoms(Acounter)%name = AName(t)%name
    Atoms(Acounter)%mass = mass(t)
    ForceMatrix(Acounter)%mass = Atoms(Acounter)%mass
    Atoms(Acounter)%type = AType(t)
  end do !i
end do !t
REWIND(icutp)
CLOSE(icutp)
IF (Acounter .ne. Nat ) THEN
   PRINT*, 'ERROR!!!' 
   PRINT*, "Wrong Number of atoms", Acounter, "instead of", Nat
   STOP
ENDIF

IF ( debug ) THEN
    write(*,*) "Done reading cutp file"
ENDIF

IF ( debug ) THEN
    write(*,*) "Opening  static atom file"
ENDIF
!........................ static atoms
open (UNIT=istatat,FILE=TRIM(filename2),STATUS='OLD', action="read")
do t=Nat+1, Nat+StaticAtN
    read(istatat,*) cor(:), force(:)
    Atoms(t)%r(:)=cor(:)*ANG_BOHR
    Atoms(t)%name = AName(NATypes+1)%name
    Atoms(t)%mass = mass(NATypes+1)
    Atoms(t)%type = AType(NATypes+1)
end do !t
close(istatat)
! ......................../static atoms
END SUBROUTINE ReadAtoms

SUBROUTINE ReadForces(icutm,icutp,filem,filep,ForceMatrix,RefCoor,Natom,mass,Disp,NATypes,Nat, Ndisp)
use constants, ONLY : ZERO,ONE,TWO debug, dp, ndim, tolerance
use type_mod, ONLY:  atom, aname, atomforce
IMPLICIT NONE
INTEGER, INTENT(in) :: icutp
INTEGER, INTENT(in) :: icutm
INTEGER, INTENT(in) :: Nat
INTEGER, INTENT(in) :: NATypes
CHARACTER(30), INTENT(in) :: filem
CHARACTER(30), INTENT(in) :: filep
REAL, INTENT(in) :: mass(NATypes)
REAL(dp), INTENT(in) :: Disp
TYPE(atom), INTENT(in) :: RefCoor(Nat)
TYPE(atomforce), INTENT(out) :: ForceMatrix(Nat*(Ndisp+1))
INTEGER, INTENT(in) :: Natom(NATypes+1)
INTEGER, INTENT(in) :: Ndisp
INTEGER :: idisp, t, i, j, readtmp, Acounter, Configcounter
REAL(dp) :: checkdispMinus, checkdispPlus(ndim)
REAL(dp) :: ForcePlus(ndim), ForceMinus(ndim)
REAL(dp) :: cor(ndim)

IF ( debug ) THEN
    write(*,*) "Opening cutm and cutp files"
    write(iout,*) "Opening cutm and cutp files"
ENDIF

open (UNIT=icutm,FILE=TRIM(filem),STATUS='OLD',action="read")
open (UNIT=icutp,FILE=TRIM(filep),STATUS='OLD',action="read")
REWIND(icutm)
REWIND(icutp)
!----------------------------------------------------
! reading the zero displacement configuration ...
! the first entry

IF ( debug ) THEN
   write(*,*) " Reading Minus 0 displacement"
ENDIF
! Displacement number
idisp=0
!..............Reading config number
read(icutm,*) readtmp
IF ( readtmp .ne. idisp  ) THEN
   PRINT*, 'ERROR!!!' 
   PRINT*, "Config m:", readtmp, "config instead of",  idisp
   STOP
ENDIF
read(icutp,*) readtmp
IF ( readtmp .ne. idisp  ) THEN
   PRINT*, 'ERROR!!!' 
   PRINT*, "Config p:", readtmp, "config instead of",  idisp
   STOP
ENDIF
!.............../Reading config number
Acounter = 0 ! Atom index
checkdispMinus=ZERO
checkdispPlus=ZERO
do t=1, NATypes, 1
  do i=1, Natom(t), 1
    Acounter = Acounter + 1
    read(icutm,*) cor(:),ForceMinus(:)
    DO j=1,ndim
        checkdispMinus=checkdispMinus+abs(cor(j)-RefCoor(Acounter)%r(j))
    ENDDO
    read(icutp,*) cor(:),ForcePlus(:)
    DO j=1,ndim
        checkdispPlus=checkdispPlus+abs(cor(j)-RefCoor(Acounter)%r(j))
    ENDDO
    !Verifies that 0 config are the same for plus and minus
    IF ( checkdispMinus .gt. tolerance ) THEN
        PRINT*, 'ERROR!!!'
        PRINT*, 'Configs #0 are differing for atom:', Acounter
        PRINT('(a11,3f15.5)'), 'Ref Config', RefCoor(Acounter)%r(:)
        PRINT('(a11,3f15.5)'), 'M Config', cor(:)
        STOP
    ENDIF
    IF ( checkdispPlus .gt. tolerance ) THEN
        PRINT*, 'ERROR!!!'
        PRINT*, 'Configs #0 are differing for atom:', Acounter
        PRINT('(a11,3f15.5)'), 'Ref Config', RefCoor(Acounter)%r(:)
        PRINT('(a11,3f15.5)'), 'P Config', cor(:)
        STOP
    ENDIF
    !Atoms(Acounter)%status=t
    ForceMatrix(Acounter)%f(:)=(ForceMinus(:)-ForcePlus(:))/TWO
    ForceMatrix(Acounter)%mass = Atoms(Acounter)%mass
  end do !i
end do !t
! reading all other displacements ...        
do idisp=1, Ndisp, 1 ! over all displacements
    !
    read(icutm,*) readtmp
    IF ( readtmp .ne. idisp  ) THEN
       PRINT*, 'ERROR!!!' 
       PRINT*,"Config m:", readtmp,"Target read",idisp,"of", Ndisp
       STOP
    ENDIF
    read(icutp,*) readtmp
    IF ( readtmp .ne. idisp  ) THEN
       PRINT*, 'ERROR!!!' 
       PRINT*, "Config p:", readtmp,"Target read",idisp,"of",Ndisp
       STOP
    ENDIF
    IF ( debug  ) THEN
       write(*,*)  "reading over disp", idisp, NATypes, Ndisp
       write(iout,*)  "reading over disp", idisp, NATypes, Ndisp
    END IF
    ! read configuration
    checkdispPlus=ZERO
    checkdispMinus=ZERO
    Configcounter=0
    IF ( Acounter .ne. Nat*idisp ) THEN
       PRINT*, 'ERROR!!!' 
       PRINT*, "The atoms are not being read correctly"
       PRINT*, 'for displacement:', idisp
       PRINT*, "Expected: ",  Nat*idisp, "got:", Acounter
       STOP
    ENDIF
    do t=1, NATypes, 1 ! over all atoms types
        do i=1, Natom(t), 1 ! over the type
            Configcounter=Configcounter+1   
            Acounter = Acounter + 1
            read(icutm,*) cor(:), ForceMinus(:)
            ! Check if displacement makes sense
            DO j=1,ndim
               checkdispMinus=checkdispMinus+abs((ANG_BOHR*cor(j))-Atoms(Configcounter)%r(j))
            ENDDO
            read(icutp,*) cor(:), ForcePlus(:)
            ! Check if displacement makes sense
            DO j=1,ndim
               checkdispPlus=checkdispPlus+abs((ANG_BOHR*cor(j))-Atoms(Configcounter)%r(j))
            ENDDO
            !Atoms(Acounter)%r(:)=(ANG_BOHR*cor(:))-Atoms(Configcounter)%r(:)
            ForceMatrix(Acounter)%f(:)= (ForceMinus(:)-ForcePlus(:))/TWO
            !Atoms(Acounter)%name = Atoms(Configcounter)%name
            ForceMatrix(Acounter)%mass = Atoms(Configcounter)%mass
            !Atoms(Acounter)%type = Atoms(Configcounter)%type
        end do   !i atom per type  
    end do !t type
    IF ( ( checkdispMinus .gt. Disp+tolerance ) .OR. ( checkdispPlus .gt. Disp+tolerance  ) ) THEN
        PRINT*, 'ERROR!!!'
        PRINT*, 'Atom Displacements too large for config:', idisp
        PRINT*, 'Displacement P (Bohr)', checkdispPlus
        PRINT*, 'Displacement M (Bohr)', checkdispMinus
        PRINT*, 'Expected (Bohr):', Disp
        STOP
    ENDIF
end do ! idisp
IF ( debug ) THEN
    write(*,*) "Closing cutm and cutp files"
    write(iout,*) "Closing cutm and cutp files"
ENDIF
close(icutp)
close(icutm)
END SUBROUTINE ReadForces
!--------------- END READING FORCE / Coordinates input files
SUBROUTINE WriteDynMat(iout,ForceMatrix,DynMat,Disp,Nat, Ndisp,N)
use constants
use type_mod, ONLY : atomforce
IMPLICIT NONE
TYPE(atomforce), INTENT(in) :: ForceMatrix(Nat*(Ndisp+1))
REAL(dp), INTENT(out) :: DynMat(N,N)
INTEGER, INTENT(in) :: iout
INTEGER, INTENT(in) :: Ndisp
INTEGER, INTENT(in) :: Nat
INTEGER, INTENT(in) :: N
REAL(dp), INTENT(in) :: Disp
REAL(dp) :: joint_mass
INTEGER :: alpha, beta, i, j,x,y ai, bi 
REAL(dp) :: A1, B1, Conv

IF ( debug ) THEN
    write(*,*) "Writing force constant matrix"
ENDIF
! We want to write our force constants
! F/u ~ dF/du where u is our displacement in bohr
! Forces are read in eV/Ang --we convert them 
!     ang to bohr and eV to Ha (2* Rydbergs)
!  i.e conv_Siesta_AU = BOHR_ANG/(TWO*RYD_eV)
! au (Ha/Bohr) to Si 823.87225(14)×10−10 N
! since our objective is to get the force constants, i.e.
! F/disp with disp in bohrs, we will 
! write our conversion factor as Ha/(Bohr2) conversion factor
! Conv_Au_Si = (Forces_AU_SI/BOHR_M) 
 ! Forces in siesta are eV/Ang
 ! Atoms(i)%f(:)=823.87225*Atoms(i)%f(:)*BOHR_ANG/(TWO*RYD_eV) 
 ! atomic mass units to SI 1.6605402e-27
 ! Displacement is in Bohr 
 ! x 1.0e24 transformation to THz factor
 ! put here to get rid of the large powers
!Conv=1.0/(0.529e-2)
! next we multiply by 1/m in AMU 
! nevertheless, we want to obtain frequency Thz later on
! which induce a 10+24 factor
! we then use  AMU to KG * 10+24 to avoid large powers
! i.e. conv mass = AMU_KGmodif
!  Everything here is double precision: 
!Conv=(BOHR_ANG/TWO*RYD_eV)*(Forces_AU_SI/BOHR_M)*(ONE/AMU_KGmodif) 
Conv=(BOHR_ANG*Forces_AU_SI)/(TWO*RYD_eV*BOHR_M*AMU_KGmodif) 
! collecting force constant matrix ...
! Dimension is 3Nat*3Nat 
y=0
do alpha=1, Nat, 1
 do i=1, ndim, 1
   x=0
   y=y+1 ! increment for raws
   do beta=1, Nat, 1
     do j=1, ndim, 1
        ! Increments:
        x=x+1 ! increment for columns
        !ai = Nat+alpha+Nat*ndim*(beta-1)+Nat*(j-1)
        ai = alpha+ Nat*(j+ndim*(beta-1))
        A1=ForceMatrix(ai)%f(i)
        !bi = Nat+beta+Nat*ndim*(alpha-1)+Nat*(i-1)
        bi = beta+Nat*(i+ndim*(alpha-1))
        B1=ForceMatrix(bi)%f(j)
        IF (debug ) THEN
         write(iout,'(a6,i4,a5,i4)') 'alpha=',alpha,'beta=',beta 
         write(iout,'(a2,i1,a2,i1)') 'i=', i,'j=',j 
         write(iout,'(a2,i1,a2,i1)') 'x=', i,'y=',j 
         write(iout,'(a3,i5,a3,g15.5)') 'ai=', ai,'A1=',A1 
         write(iout,'(a3,i5,a3,g15.5)') 'bi=', ai,'B1=',A1 
        ENDIF
        joint_mass = sqrt(ForceMatrix(ai)%mass1*ForceMatrix(bi)%mass) !*1.66e-3
        !joint_mass = *1.66e-3
        !write(*,*) 'mass = ', joint_mass
        !Disp is in Bohr
        !Dm(alpha,i,beta,j) = 0.5*((A1+B1)/Disp)/joint_mass
        ! This formula used to read             
        !Dm(alpha,i,beta,j) = -0.5*((A2-A1)/(2.0*Disp)+ (B2-B1)/(2.0*Disp))/joint_mass
        ! since B2 and a2 were defined as -B1 and -B2...
        ! it was not necessary
        !A(x,y)=Dm(alpha,i,beta,j)
        A(x,y)= Conv*(A1+B1)/(TWO*Disp*joint_mass)
       end do!j 
   end do!beta
  end do!j  
end do!alpha   
! .............................
!--------------- DONE symmetrize/scale/prepare for diag

END SUBROUTINE WriteDynMat
!CALL write_NMD_File(imod,filename,Atoms,StaticAt,A,Nat,StaticAtN,ndim)
SUBROUTINE  write_NMD_File(ifile,fname,Atlist,StatAtlist,DynMat,Nat,NStAt,ndim)
use constants, ONLY : BOHR_ANG, debug
IMPLICIT NONE
INTEGER, INTENT(in) :: ifile
CHARACTER(30), INTENT(IN) :: fname
TYPE(atom), INTENT(IN) :: Atlist(Nat)
TYPE(atom), INTENT(IN) :: StatAtlist(NStAt)
REAL(dp), INTENT(IN) :: DynMat(Nat*ndim,Nat*ndim)
INTEGER, INTENT(in) :: Nat
INTEGER, INTENT(in) :: NStAt
INTEGER, INTENT(in) :: ndim
! Loop variables
Integer :: iMode,i, Write_Counter, StWrite_Counter

open(UNIT=ifile, file=fname,form='formatted',status='unknown')
!===================================
! writing NMD file
Write_counter = 0
StWrite_Counter = 0
write(ifile,'(a16)',advance="no") "atomnames  "
do i=1, (Nat+NStAt), 1
    if (i .le. Nat) then 
        Write_counter = Write_counter + 1
        write(ifile,'(a4)',advance="no") Atlist(Write_counter)%name
    else
        StWrite_counter = StWrite_counter + 1
        write(ifile,'(a4)',advance="no") StatAtlist(StWrite_counter)%name
    endif
end do
write(ifile,'(a1)') " "
Write_counter = 0
StWrite_counter = 0
write(ifile,'(a16)',advance="no") "coordinates  "
do i=1, (Nat+NStAt), 1
    if (i .le. Nat) then 
        Write_counter = Write_counter + 1
        write(ifile,'(3f14.8)',advance="no") Atoms(Write_counter)%r(:)*BOHR_ANG
    else
        StWrite_counter = StWrite_counter + 1
        write(ifile,'(3f14.8)',advance="no") StaticAt(StWrite_counter)%r(:)
    endif
end do
write(ifile,'(a1)') " "
do iMode=1, Nat*ndim, 1
    Write_counter = 0
    write(ifile,'(a8)',advance="no") "mode  "
    do i=1, (Nat+NStAt), 1
        if (i .le. Nat) then 
            Write_counter = Write_counter + 1
            write(ifile,'(3f14.8)',advance="no") DynMat((1+ndim*(Write_counter-1)):(ndim*Write_counter),Mode)*BOHR_ANG   !/sqrt(Atoms(Write_counter)%mass)
        else
            write(ifile,'(3f14.8)',advance="no") 0.0,  0.0,  0.0
        endif
    end do
end do
write(ifile,'(a1)') " "
!====================================
CLOSE(ifile)
END SUBROUTINE write_NMD_File

!====================================

SUBROUTINE   ReadSiestaCor(icor,fname,AName(:),Natom(:),mass(:),NATypes)
use constants, only : debug
use type_mod, only : atom, aname
IMPLICIT NONE
INTEGER, INTENT(in) :: icor
CHARACTER(30), INTENT(IN) :: fname
INTEGER, INTENT(in) :: NATypes
REAL, INTENT(out) :: mass(NATypes)
TYPE(aname), INTENT(out) :: AName(NATypes)
INTEGER, INTENT(out) :: Natom(NAtypes)
! Function
REAL(dp) :: whatismymass
! Loop variables
Integer :: t,Acounter

open (UNIT=icor,FILE=fname,STATUS='OLD',action="read")
IF ( debug ) THEN
    write(*,*) "Reading siesta.cor"
ENDIF

Acounter = 0
do t=1, NATypes, 1
  
  IF ( debug ) THEN
     WRITE(*,*) "Type=", t, 'Counter', Acounter
  END IF

  read(icor,*) AName(t), Natom(t)
  mass(t) = whatismymass(AName(t))
  write (*,*) AName(t),  Natom(t)
 ! PRINT*, AName(t),  Natom(t)
  do i=1, Natom(t), 1 
   Acounter = Acounter + 1
     !PRINT*, Acounter
     read(icor,*) cor(:), inttmp, chartmp, readtmp
  end do 
  AType(t)=inttmp ! take the type number of the last read atom
                 ! of the same type 
end do

if ( debug ) then 
    write(*,*) 'Number of atoms read', Acounter
endif
IF ( debug ) THEN
    write(*,*) "Closing siesta.cor"
    write(iout,*) "Closing siesta.cor"
ENDIF
CLOSE(icor)

END SUBROUTINE   ReadSiestaCor

End program vwhole


FUNCTION whatismymass(atomname) result(Mass)
IMPLICIT NONE
CHARACTER(15), INTENT(in) :: atomname    ! name of object 
integer, parameter :: dp = kind(1.0d0)
real(dp) :: Mass

IF ( TRIM(atomname) == 'H') THEN
    Mass=1.0079D0
    WRITE(*,*) 'Found: Type      Hydrogen'
    WRITE(*,*) 'Atomic #   1             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'He') THEN
    Mass=4.0026D0 
    WRITE(*,*) 'Found: Type        Helium'
    WRITE(*,*) 'Atomic #   2             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Li') THEN
    Mass=6.9410D0 
    WRITE(*,*) 'Found: Type       Lithium'
    WRITE(*,*) 'Atomic #   3             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Be') THEN
    Mass=9.0122D0  
    WRITE(*,*) 'Found: Type     Beryllium'
    WRITE(*,*) 'Atomic #   4             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'B') THEN
    Mass=10.8110D0 
    WRITE(*,*) 'Found: Type         Boron'
    WRITE(*,*) 'Atomic #   5             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'C') THEN
    Mass=12.0107D0 
    WRITE(*,*) 'Found: Type        Carbon'
    WRITE(*,*) 'Atomic #   6             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'N') THEN
    Mass=14.0067D0 
    WRITE(*,*) 'Found: Type      Nitrogen'
    WRITE(*,*) 'Atomic #   7             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'O') THEN
    Mass=15.9994D0 
    WRITE(*,*) 'Found: Type        Oxygen'
    WRITE(*,*) 'Atomic #   8             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'F') THEN
    Mass=18.9984D0 
    WRITE(*,*) 'Found: Type      Fluorine'
    WRITE(*,*) 'Atomic #   9             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ne') THEN
    Mass=20.1797D0 
    WRITE(*,*) 'Found: Type          Neon'
    WRITE(*,*) 'Atomic #  10             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Na') THEN
    Mass=22.9897D0 
    WRITE(*,*) 'Found: Type        Sodium'
    WRITE(*,*) 'Atomic #  11             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Mg') THEN
    Mass=24.3050D0 
    WRITE(*,*) 'Found: Type     Magnesium'
    WRITE(*,*) 'Atomic #  12             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Al') THEN
    Mass=26.9815D0 
    WRITE(*,*) 'Found: Type      Aluminum'
    WRITE(*,*) 'Atomic #  13             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Si') THEN
    Mass=28.0855D0 
    WRITE(*,*) 'Found: Type       Silicon'
    WRITE(*,*) 'Atomic #  14             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'P') THEN
    Mass=30.9738D0 
    WRITE(*,*) 'Found: Type    Phosphorus'
    WRITE(*,*) 'Atomic #  15             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'S') THEN
    Mass=32.0650D0 
    WRITE(*,*) 'Found: Type        Sulfur'
    WRITE(*,*) 'Atomic #  16             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Cl') THEN
    Mass=35.4530D0 
    WRITE(*,*) 'Found: Type      Chlorine'
    WRITE(*,*) 'Atomic #  17             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'K') THEN
    Mass=39.0983D0 
    WRITE(*,*) 'Found: Type     Potassium'
    WRITE(*,*) 'Atomic #  19             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ar') THEN
    Mass=39.9480D0 
    WRITE(*,*) 'Found: Type         Argon'
    WRITE(*,*) 'Atomic #  18             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ca') THEN
    Mass=40.0780D0 
    WRITE(*,*) 'Found: Type       Calcium'
    WRITE(*,*) 'Atomic #  20             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Sc') THEN
    Mass=44.9559D0 
    WRITE(*,*) 'Found: Type      Scandium'
    WRITE(*,*) 'Atomic #  21             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ti') THEN
    Mass=47.8670D0 
    WRITE(*,*) 'Found: Type      Titanium'
    WRITE(*,*) 'Atomic #  22             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'V') THEN
    Mass=50.9415D0 
    WRITE(*,*) 'Found: Type      Vanadium'
    WRITE(*,*) 'Atomic #  23             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Cr') THEN
    Mass=51.9961D0 
    WRITE(*,*) 'Found: Type      Chromium'
    WRITE(*,*) 'Atomic #  24             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Mn') THEN
    Mass=54.9380D0 
    WRITE(*,*) 'Found: Type     Manganese'
    WRITE(*,*) 'Atomic #  25             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Fe') THEN
    Mass=55.8450D0 
    WRITE(*,*) 'Found: Type          Iron'
    WRITE(*,*) 'Atomic #  26             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ni') THEN
    Mass=58.6934D0 
    WRITE(*,*) 'Found: Type        Nickel'
    WRITE(*,*) 'Atomic #  28             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Co') THEN
    Mass=58.9332D0 
    WRITE(*,*) 'Found: Type        Cobalt'
    WRITE(*,*) 'Atomic #  27             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Cu') THEN
    Mass=63.5460D0 
    WRITE(*,*) 'Found: Type        Copper'
    WRITE(*,*) 'Atomic #  29             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Zn') THEN
    Mass=65.3900D0 
    WRITE(*,*) 'Found: Type          Zinc'
    WRITE(*,*) 'Atomic #  30             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ga') THEN
    Mass=69.7230D0 
    WRITE(*,*) 'Found: Type       Gallium'
    WRITE(*,*) 'Atomic #  31             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ge') THEN
    Mass=72.6400D0 
    WRITE(*,*) 'Found: Type     Germanium'
    WRITE(*,*) 'Atomic #  32             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'As') THEN
    Mass=74.9216D0 
    WRITE(*,*) 'Found: Type       Arsenic'
    WRITE(*,*) 'Atomic #  33             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Se') THEN
    Mass=78.9600D0 
    WRITE(*,*) 'Found: Type      Selenium'
    WRITE(*,*) 'Atomic #  34             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Br') THEN
    Mass=79.9040D0 
    WRITE(*,*) 'Found: Type       Bromine'
    WRITE(*,*) 'Atomic #  35             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Kr') THEN
    Mass=83.8000D0 
    WRITE(*,*) 'Found: Type       Krypton'
    WRITE(*,*) 'Atomic #  36             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Rb') THEN
    Mass=85.4678D0 
    WRITE(*,*) 'Found: Type      Rubidium'
    WRITE(*,*) 'Atomic #  37             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Sr') THEN
    Mass=87.6200D0 
    WRITE(*,*) 'Found: Type     Strontium'
    WRITE(*,*) 'Atomic #  38             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Y') THEN
    Mass=88.9059D0 
    WRITE(*,*) 'Found: Type       Yttrium'
    WRITE(*,*) 'Atomic #  39             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Zr') THEN
    Mass=91.2240D0 
    WRITE(*,*) 'Found: Type     Zirconium'
    WRITE(*,*) 'Atomic #  40             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Nb') THEN
    Mass=92.9064D0
    WRITE(*,*) 'Found: Type       Niobium'
    WRITE(*,*) 'Atomic #  41             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Mo') THEN
    Mass=95.9400D0 
    WRITE(*,*) 'Found: Type    Molybdenum'
    WRITE(*,*) 'Atomic #  42             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Tc') THEN
    Mass=98.0000D0 
    WRITE(*,*) 'Found: Type    Technetium'
    WRITE(*,*) 'Atomic #  43             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ru') THEN
    Mass=101.0700D0
    WRITE(*,*) 'Found: Type     Ruthenium'
    WRITE(*,*) 'Atomic #  44             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Rh') THEN
    Mass=102.9055D0
    WRITE(*,*) 'Found: Type       Rhodium'
    WRITE(*,*) 'Atomic #  45             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Pd') THEN
    Mass=106.4200D0
    WRITE(*,*) 'Found: Type     Palladium'
    WRITE(*,*) 'Atomic #  46             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ag') THEN
    Mass=107.8682D0
    WRITE(*,*) 'Found: Type        Silver'
    WRITE(*,*) 'Atomic #  47             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Cd') THEN
    Mass=112.4110D0
    WRITE(*,*) 'Found: Type       Cadmium'
    WRITE(*,*) 'Atomic #  48             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'In') THEN
    Mass=114.8180D0
    WRITE(*,*) 'Found: Type        Indium'
    WRITE(*,*) 'Atomic #  49             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Sn') THEN
    Mass=118.7100D0
    WRITE(*,*) 'Found: Type           Tin'
    WRITE(*,*) 'Atomic #  50             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Sb') THEN
    Mass=121.7600D0
    WRITE(*,*) 'Found: Type      Antimony'
    WRITE(*,*) 'Atomic #  51             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'I') THEN
    Mass=126.9045D0
    WRITE(*,*) 'Found: Type        Iodine'
    WRITE(*,*) 'Atomic #  53             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Te') THEN
    Mass=127.6000D0
    WRITE(*,*) 'Found: Type     Tellurium'
    WRITE(*,*) 'Atomic #  52             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Xe') THEN
    Mass=131.2930D0
    WRITE(*,*) 'Found: Type         Xenon'
    WRITE(*,*) 'Atomic #  54             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Cs') THEN
    Mass=132.9055D0
    WRITE(*,*) 'Found: Type        Cesium'
    WRITE(*,*) 'Atomic #  55             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ba') THEN
    Mass=137.3270D0
    WRITE(*,*) 'Found: Type        Barium'
    WRITE(*,*) 'Atomic #  56             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'La') THEN
    Mass=138.9055D0
    WRITE(*,*) 'Found: Type     Lanthanum'
    WRITE(*,*) 'Atomic #  57             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ce') THEN
    Mass=140.1160D0
    WRITE(*,*) 'Found: Type        Cerium'
    WRITE(*,*) 'Atomic #  58             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Pr') THEN
    Mass=140.9077D0
    WRITE(*,*) 'Found: Type  Praseodymium'
    WRITE(*,*) 'Atomic #  59             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Nd') THEN
    Mass=144.2400D0
    WRITE(*,*) 'Found: Type     Neodymium'
    WRITE(*,*) 'Atomic #  60             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Pm') THEN
    Mass=145.0000D0
    WRITE(*,*) 'Found: Type    Promethium'
    WRITE(*,*) 'Atomic #  61             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Sm') THEN
    Mass=150.3600D0
    WRITE(*,*) 'Found: Type      Samarium'
    WRITE(*,*) 'Atomic #  62             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Eu') THEN
    Mass=151.9640D0
    WRITE(*,*) 'Found: Type      Europium'
    WRITE(*,*) 'Atomic #  63             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Gd') THEN
    Mass=157.2500D0
    WRITE(*,*) 'Found: Type    Gadolinium'
    WRITE(*,*) 'Atomic #  64             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Tb') THEN
    Mass=158.9253D0
    WRITE(*,*) 'Found: Type       Terbium'
    WRITE(*,*) 'Atomic #  65             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Dy') THEN
    Mass=162.5000D0
    WRITE(*,*) 'Found: Type    Dysprosium'
    WRITE(*,*) 'Atomic #  66             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ho') THEN
    Mass=164.9303D0
    WRITE(*,*) 'Found: Type       Holmium'
    WRITE(*,*) 'Atomic #  67             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Er') THEN
    Mass=167.2590D0
    WRITE(*,*) 'Found: Type        Erbium'
    WRITE(*,*) 'Atomic #  68             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Tm') THEN
    Mass=168.9342D0
    WRITE(*,*) 'Found: Type       Thulium'
    WRITE(*,*) 'Atomic #  69             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Yb') THEN
    Mass=173.0400D0
    WRITE(*,*) 'Found: Type     Ytterbium'
    WRITE(*,*) 'Atomic #  70             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Lu') THEN
    Mass=174.9670D0
    WRITE(*,*) 'Found: Type      Lutetium'
    WRITE(*,*) 'Atomic #  71             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Hf') THEN
    Mass=178.4900D0
    WRITE(*,*) 'Found: Type       Hafnium'
    WRITE(*,*) 'Atomic #  72             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ta') THEN
    Mass=180.9479D0
    WRITE(*,*) 'Found: Type      Tantalum'
    WRITE(*,*) 'Atomic #  73             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'W') THEN
    Mass=183.8400D0
    WRITE(*,*) 'Found: Type      Tungsten'
    WRITE(*,*) 'Atomic #  74             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Re') THEN
    Mass=186.2070D0
    WRITE(*,*) 'Found: Type       Rhenium'
    WRITE(*,*) 'Atomic #  75             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Os') THEN
    Mass=190.2300D0
    WRITE(*,*) 'Found: Type        Osmium'
    WRITE(*,*) 'Atomic #  76             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ir') THEN
    Mass=192.2170D0
    WRITE(*,*) 'Found: Type       Iridium'
    WRITE(*,*) 'Atomic #  77             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Pt') THEN
    Mass=195.0780D0
    WRITE(*,*) 'Found: Type      Platinum'
    WRITE(*,*) 'Atomic #  78             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Au') THEN
    Mass=196.9665D0
    WRITE(*,*) 'Found: Type          Gold'
    WRITE(*,*) 'Atomic #  79             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Hg') THEN
    Mass=200.5900D0
    WRITE(*,*) 'Found: Type       Mercury'
    WRITE(*,*) 'Atomic #  80             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Tl') THEN
    Mass=204.3833D0
    WRITE(*,*) 'Found: Type      Thallium'
    WRITE(*,*) 'Atomic #  81             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Pb') THEN
    Mass=207.2000D0
    WRITE(*,*) 'Found: Type          Lead'
    WRITE(*,*) 'Atomic #  82             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Bi') THEN
    Mass=208.9804D0
    WRITE(*,*) 'Found: Type       Bismuth'
    WRITE(*,*) 'Atomic #  83             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Po') THEN
    Mass=209.0000D0
    WRITE(*,*) 'Found: Type      Polonium'
    WRITE(*,*) 'Atomic #  84             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'At') THEN
    Mass=210.0000D0
    WRITE(*,*) 'Found: Type      Astatine'
    WRITE(*,*) 'Atomic #  85             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Rn') THEN
    Mass=222.0000D0
    WRITE(*,*) 'Found: Type         Radon'
    WRITE(*,*) 'Atomic #  86             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Fr') THEN
    Mass=223.0000D0
    WRITE(*,*) 'Found: Type      Francium'
    WRITE(*,*) 'Atomic #  87             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ra') THEN
    Mass=226.0000D0
    WRITE(*,*) 'Found: Type        Radium'
    WRITE(*,*) 'Atomic #  88             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Ac') THEN
    Mass=227.0000D0
    WRITE(*,*) 'Found: Type      Actinium'
    WRITE(*,*) 'Atomic #  89             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Pa') THEN
    Mass=231.0359D0
    WRITE(*,*) 'Found: Type  Protactinium'
    WRITE(*,*) 'Atomic #  91             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Th') THEN
    Mass=232.0381D0
    WRITE(*,*) 'Found: Type       Thorium'
    WRITE(*,*) 'Atomic #  90             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Np') THEN
    Mass=237.0000D0
    WRITE(*,*) 'Found: Type     Neptunium'
    WRITE(*,*) 'Atomic #  93             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'U') THEN
    Mass=238.0289D0
    WRITE(*,*) 'Found: Type       Uranium'
    WRITE(*,*) 'Atomic #  92             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Am') THEN
    Mass=243.0000D0
    WRITE(*,*) 'Found: Type     Americium'
    WRITE(*,*) 'Atomic #  95             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Pu') THEN
    Mass=244.0000D0
    WRITE(*,*) 'Found: Type     Plutonium'
    WRITE(*,*) 'Atomic #  94             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Cm') THEN
    Mass=247.0000D0
    WRITE(*,*) 'Found: Type        Curium'
    WRITE(*,*) 'Atomic #  96             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Bk') THEN
    Mass=247.0000D0
    WRITE(*,*) 'Found: Type     Berkelium'
    WRITE(*,*) 'Atomic #  97             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Cf') THEN
    Mass=251.0000D0
    WRITE(*,*) 'Found: Type   Californium'
    WRITE(*,*) 'Atomic #  98             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Es') THEN
    Mass=252.0000D0
    WRITE(*,*) 'Found: Type   Einsteinium'
    WRITE(*,*) 'Atomic #  99             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Fm') THEN
    Mass=257.0000D0
    WRITE(*,*) 'Found: Type       Fermium'
    WRITE(*,*) 'Atomic # 100             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Md') THEN
    Mass=258.0000D0
    WRITE(*,*) 'Found: Type   Mendelevium'
    WRITE(*,*) 'Atomic # 101             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'No') THEN
    Mass=259.0000D0
    WRITE(*,*) 'Found: Type      Nobelium'
    WRITE(*,*) 'Atomic # 102             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Rf') THEN
    Mass=261.0000D0
    WRITE(*,*) 'Found: Type Rutherfordium'
    WRITE(*,*) 'Atomic # 104             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Lr') THEN
    Mass=262.0000D0
    WRITE(*,*) 'Found: Type    Lawrencium'
    WRITE(*,*) 'Atomic # 103             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Db') THEN
    Mass=262.0000D0
    WRITE(*,*) 'Found: Type       Dubnium'
    WRITE(*,*) 'Atomic # 105             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Bh') THEN
    Mass=264.0000D0
    WRITE(*,*) 'Found: Type       Bohrium'
    WRITE(*,*) 'Atomic # 107             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Sg') THEN
    Mass=266.0000D0
    WRITE(*,*) 'Found: Type    Seaborgium'
    WRITE(*,*) 'Atomic # 106             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Mt') THEN
    Mass=268.0000D0
    WRITE(*,*) 'Found: Type    Meitnerium'
    WRITE(*,*) 'Atomic # 109             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Rg') THEN
    Mass=272.0000D0
    WRITE(*,*) 'Found: Type   Roentgenium'
    WRITE(*,*) 'Atomic # 111             Mass:',  Mass 
ELSE IF ( TRIM(atomname) == 'Hs') THEN
    Mass=277.0000D0
    WRITE(*,*) 'Found: Type       Hassium'
    WRITE(*,*) 'Atomic # 108             Mass:',  Mass 
ELSE 
  WRITE(*,*) 'ERROR MASS NOT FOUND for', atomname
  STOP
ENDIF
ENDFUNCTION

