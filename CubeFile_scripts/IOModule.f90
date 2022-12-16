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
   MODULE IO_Module
   !***********************************************
   IMPLICIT NONE
   PRIVATE 
   SAVE


   ! Public routines:
   PUBLIC                 ::  ReadInputFileOverlap
   PUBLIC                 ::  gnuplotwrite
   PUBLIC                 ::  datwrite
   PUBLIC                 ::  ReadInputFileDensity
   PUBLIC                 ::  ReadInputFilePotential
   PUBLIC                 ::  ReadInputFileWFinPotential
   PUBLIC                 ::  ReadInputFileDipole
   PUBLIC                 ::  ReadInputFilePlot
   PUBLIC                 ::  ReadInputFilePlotPotential
   PUBLIC                 ::  GaussianFilereadHeader
   PUBLIC                 ::  GaussianFilewrite
   PUBLIC                 ::  GaussianFileread
   PUBLIC                 ::  RhoFilereadHeader
   PUBLIC                 ::  RhoFileread
   PUBLIC                 ::  PotentialFilereadHeader
   PUBLIC                 ::  PotentialFileread
   PUBLIC                 ::  XYZFilereadHeader
   PUBLIC                 ::  XYZFileread

CONTAINS

!***********************************************
      SUBROUTINE ReadInputFileOverlap(iufile, datafile, FileInput, gold_position, Nb_distances, maxdistance, Reference_Atom,  FileOutputGaussian)

   !***********************************************
   ! For density calculations
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)               :: iufile
   CHARACTER*100, INTENT(in)         :: datafile
   CHARACTER*100, INTENT(out)        :: FileInput
   REAL(dbl), INTENT(out)            :: gold_position(3)
   INTEGER,   INTENT(out)            :: Nb_distances
   INTEGER,   INTENT(out)            :: Reference_Atom
   REAL(dbl), INTENT(out)            :: maxdistance
   CHARACTER*100, INTENT(out)        :: FileOutputGaussian

   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInput
        PRINT*, "...READING Incoming Data File: Input file                       :", Trim(FileInput)
        READ(iufile,*) gold_position(:)
        PRINT*, "...READING Incoming Data File: Gold Position                    :", gold_position(:)
        READ(iufile,*) Nb_distances  
        PRINT*, "...READING Incoming Data File:  Nb_distances for Overlap calc   :", Nb_distances  
        READ(iufile,*)  maxdistance 
        PRINT*, "...READING Incoming Data File:  maxdistance for Overlap calc    :", maxdistance
        READ(iufile,*)  Reference_Atom 
        PRINT*, "...READING Incoming Data File:  Reference_Atom for Overlap calc    :", Reference_Atom


        READ(iufile,*)  FileOutputGaussian 
        PRINT*, "...READING Incoming Data File:  Output file for DoubleCheck      :", Trim(FileOutputGaussian)
 
        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFileOverlap


!***********************************************
      SUBROUTINE gnuplotwrite(iufile, datafile, dimensions, nvect1, nvect2, vector1, vector2 , function_in)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)              :: iufile
   CHARACTER*100, INTENT(in)        :: datafile
   CHARACTER*100  :: chr
   INTEGER,   INTENT(in)           :: dimensions
   INTEGER,     INTENT(in)         :: nvect1
   INTEGER,     INTENT(in)         :: nvect2
   REAL(dbl),   INTENT(in)         :: vector1(dimensions)
   REAL(dbl),   INTENT(in)         :: vector2(dimensions)
   REAL(dbl),   INTENT(in)         :: function_in(nvect1, nvect2) 
   REAL(dbl)   :: X_value, Y_value
   INTEGER :: ia, ib     

       PRINT*, "    ...PRITING Gnuplot Output File: ", TRIM(datafile)
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
      SUBROUTINE datwrite(iufile, datafile, nvect, vect, function_in)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)              :: iufile
   CHARACTER*100, INTENT(in)        :: datafile
   CHARACTER*100  :: chr
   INTEGER,     INTENT(in)         :: nvect
   REAL(dbl),   INTENT(in)         :: vect
   REAL(dbl),   INTENT(in)         :: function_in(nvect) 
   REAL(dbl)   :: X_value
   INTEGER :: ia, ib     

       PRINT*, "    ...PRITING .dat Output File: ", TRIM(datafile)
        
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        DO ia = 1,  nvect
                 X_value = (ia-1)* vect
                 WRITE(iufile,*) X_value,  function_in(ia)  
        ENDDO


        CLOSE(iufile)


END  SUBROUTINE datwrite

!***********************************************
      SUBROUTINE ReadInputFileDensity(iufile, datafile, FileInput0bias, FileInputFinitebias, FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat, FileOutputDat0bias, FileOutputDatFinitebias)
   !***********************************************
   ! For density calculations
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)              :: iufile
   CHARACTER*100, INTENT(in)         :: datafile
   CHARACTER*100, INTENT(out)        :: FileInput0bias
   CHARACTER*100, INTENT(out)        :: FileInputFinitebias
   CHARACTER*100, INTENT(out)        :: FileInputxyz
   CHARACTER*100, INTENT(out)        :: FileOutputGaussian
   CHARACTER*100, INTENT(out)        :: FileOutputGnuplotx
   CHARACTER*100, INTENT(out)        :: FileOutputGnuploty
   CHARACTER*100, INTENT(out)        :: FileOutputDat
   CHARACTER*100, INTENT(out)        :: FileOutputDat0bias
   CHARACTER*100, INTENT(out)        :: FileOutputDatFinitebias

   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInput0bias
        PRINT*, "...READING Incoming Data File: Input file for 0 bias           :", Trim(FileInput0bias)
        READ(iufile,*) FileInputFinitebias
        PRINT*, "...READING Incoming Data File: Input file for Finite bias      :", Trim(FileInputFinitebias)
        READ(iufile,*) FileInputxyz
        PRINT*, "...READING Incoming Data File: Output file for xyz              :", Trim(FileInputxyz)
        READ(iufile,*)  FileOutputGaussian 
        PRINT*, "...READING Incoming Data File: Output file for Output Gaussian  :", Trim( FileOutputGaussian )
        READ(iufile,*)  FileOutputGnuplotx 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage x :", Trim(FileOutputGnuplotx)
        READ(iufile,*)  FileOutputGnuploty 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage y :", Trim(FileOutputGnuploty)
        READ(iufile,*)  FileOutputDat
        PRINT*, "...READING Incoming Data File: Output file for Average xy :", Trim(FileOutputDat)
        READ(iufile,*)  FileOutputDat0bias
        PRINT*, "...READING Incoming Data File: Output file for Average xy 0 bias      :", Trim(FileOutputDat0bias)
        READ(iufile,*)  FileOutputDatFinitebias
        PRINT*, "...READING Incoming Data File: Output file for Average xy Finite bias :", Trim(FileOutputDatFinitebias)

        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFileDensity
!***********************************************
      SUBROUTINE ReadInputFilePlot(iufile, datafile, FileInput, FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat)    
   !***********************************************
   ! For density calculations
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)               :: iufile
   CHARACTER*100, INTENT(in)         :: datafile
   CHARACTER*100, INTENT(out)        :: FileInput
   CHARACTER*100, INTENT(out)        :: FileInputxyz
   CHARACTER*100, INTENT(out)        :: FileOutputGaussian
   CHARACTER*100, INTENT(out)        :: FileOutputGnuplotx
   CHARACTER*100, INTENT(out)        :: FileOutputGnuploty
   CHARACTER*100, INTENT(out)        :: FileOutputDat

   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInput
        PRINT*, "...READING Incoming Data File: Input file                       :", Trim(FileInput)
        READ(iufile,*) FileInputxyz
        PRINT*, "...READING Incoming Data File: Output file for xyz              :", Trim(FileInputxyz)
        READ(iufile,*)  FileOutputGaussian 
        PRINT*, "...READING Incoming Data File: Output file for Output Gaussian  :", Trim( FileOutputGaussian )
        READ(iufile,*)  FileOutputGnuplotx 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage x :", Trim(FileOutputGnuplotx)
        READ(iufile,*)  FileOutputGnuploty 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage y :", Trim(FileOutputGnuploty)
        READ(iufile,*)  FileOutputDat
        PRINT*, "...READING Incoming Data File: Output file for Average xy :", Trim(FileOutputDat)

        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFilePlot
!***********************************************
      SUBROUTINE ReadInputFilePlotPotential(iufile, datafile, FileInput,FileInputGrid, FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat)    
   !***********************************************
   ! For density calculations
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)               :: iufile
   CHARACTER*100, INTENT(in)         :: datafile
   CHARACTER*100, INTENT(out)        :: FileInput
   CHARACTER*100, INTENT(out)        :: FileInputGrid
   CHARACTER*100, INTENT(out)        :: FileInputxyz
   CHARACTER*100, INTENT(out)        :: FileOutputGaussian
   CHARACTER*100, INTENT(out)        :: FileOutputGnuplotx
   CHARACTER*100, INTENT(out)        :: FileOutputGnuploty
   CHARACTER*100, INTENT(out)        :: FileOutputDat

   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInput
        PRINT*, "...READING Incoming Data File: Input file                       :", Trim(FileInput)
        READ(iufile,*) FileInputGrid
        PRINT*, "...READING Incoming Data File: Input file Grid                  :", Trim(FileInputGrid)
        READ(iufile,*) FileInputxyz
        PRINT*, "...READING Incoming Data File: Output file for xyz              :", Trim(FileInputxyz)
        READ(iufile,*)  FileOutputGaussian 
        PRINT*, "...READING Incoming Data File: Output file for Output Gaussian  :", Trim( FileOutputGaussian )
        READ(iufile,*)  FileOutputGnuplotx 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage x :", Trim(FileOutputGnuplotx)
        READ(iufile,*)  FileOutputGnuploty 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage y :", Trim(FileOutputGnuploty)
        READ(iufile,*)  FileOutputDat
        PRINT*, "...READING Incoming Data File: Output file for Average xy :", Trim(FileOutputDat)

        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFilePlotPotential
!***********************************************
      SUBROUTINE ReadInputFileDipole(iufile, datafile, FileInput, FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat)
   !***********************************************
   ! For density calculations
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)               :: iufile
   CHARACTER*100, INTENT(in)         :: datafile
   CHARACTER*100, INTENT(out)        :: FileInput
   CHARACTER*100, INTENT(out)        :: FileInputxyz
   CHARACTER*100, INTENT(out)        :: FileOutputGaussian
   CHARACTER*100, INTENT(out)        :: FileOutputGnuplotx
   CHARACTER*100, INTENT(out)        :: FileOutputGnuploty
   CHARACTER*100, INTENT(out)        :: FileOutputDat

   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File: ", TRIM(datafile)
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInput
        PRINT*, "...READING Incoming Data File: Input file                       :", Trim(FileInput)
        READ(iufile,*) FileInputxyz
        PRINT*, "...READING Incoming Data File: Output file for xyz              :", Trim(FileInputxyz)
        READ(iufile,*)  FileOutputGaussian 
        PRINT*, "...READING Incoming Data File: Output file for Output Gaussian  :", Trim( FileOutputGaussian )
        READ(iufile,*)  FileOutputGnuplotx 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage x :", Trim(FileOutputGnuplotx)
        READ(iufile,*)  FileOutputGnuploty 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage y :", Trim(FileOutputGnuploty)
        READ(iufile,*)  FileOutputDat
        PRINT*, "...READING Incoming Data File: Output file for Average xy :", Trim(FileOutputDat)

        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFileDipole
!***********************************************
      SUBROUTINE ReadInputFilePotential(iufile, datafile, FileInputGrid0bias, FileInputGridFinitebias,FileInput0bias, FileInputFinitebias,  FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat, FileOutputDat0bias, FileOutputDatFinitebias)
   !***********************************************

   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)               :: iufile
   CHARACTER*100, INTENT(in)         :: datafile
   CHARACTER*100, INTENT(out)        :: FileInputGrid0bias
   CHARACTER*100, INTENT(out)        :: FileInputGridFinitebias
   CHARACTER*100, INTENT(out)        :: FileInput0bias
   CHARACTER*100, INTENT(out)        :: FileInputFinitebias
   CHARACTER*100, INTENT(out)        :: FileInputxyz
   CHARACTER*100, INTENT(out)        :: FileOutputGaussian
   CHARACTER*100, INTENT(out)        :: FileOutputGnuplotx
   CHARACTER*100, INTENT(out)        :: FileOutputGnuploty
   CHARACTER*100, INTENT(out)        :: FileOutputDat
   CHARACTER*100, INTENT(out)        :: FileOutputDat0bias
   CHARACTER*100, INTENT(out)        :: FileOutputDatFinitebias

   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File:"
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInputGrid0bias
        PRINT*, "...READING Incoming Data File: Input file for Grid 0 bias      :", Trim(FileInputGrid0bias)
        READ(iufile,*) FileInputGridFinitebias
        PRINT*, "...READING Incoming Data File: Input file for Grid Finite bias :", Trim(FileInputGridFinitebias)
        READ(iufile,*) FileInput0bias
        PRINT*, "...READING Incoming Data File: Input file for 0 bias           :", Trim(FileInput0bias)
        READ(iufile,*) FileInputFinitebias
        PRINT*, "...READING Incoming Data File: Input file for Finite bias      :", Trim(FileInputFinitebias)
        READ(iufile,*) FileInputxyz
        PRINT*, "...READING Incoming Data File: Output file for xyz              :", Trim(FileInputxyz)
        READ(iufile,*)  FileOutputGaussian 
        PRINT*, "...READING Incoming Data File: Output file for Output Gaussian  :", Trim( FileOutputGaussian )
        READ(iufile,*)  FileOutputGnuplotx 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage x :", Trim(FileOutputGnuplotx)
        READ(iufile,*)  FileOutputGnuploty 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage y :", Trim(FileOutputGnuploty)
        READ(iufile,*)  FileOutputDat
        PRINT*, "...READING Incoming Data File: Output file for Average xy :", Trim(FileOutputDat)
        READ(iufile,*)  FileOutputDat0bias
        PRINT*, "...READING Incoming Data File: Output file for Average xy 0 bias      :", Trim(FileOutputDat0bias)
        READ(iufile,*)  FileOutputDatFinitebias
        PRINT*, "...READING Incoming Data File: Output file for Average xy Finite bias :", Trim(FileOutputDatFinitebias)

        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFilePotential
!***********************************************
      SUBROUTINE ReadInputFileWFinPotential(iufile, datafile, FileInputGrid0bias, FileInputGridFinitebias,FileInput0bias, FileInputFinitebias,  FileInputxyz ,  FileInputWF, FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat, FileOutputDat0bias, FileOutputDatFinitebias)
   !***********************************************

   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER, INTENT(in)               :: iufile
   CHARACTER*100, INTENT(in)         :: datafile
   CHARACTER*100, INTENT(out)        :: FileInputGrid0bias
   CHARACTER*100, INTENT(out)        :: FileInputGridFinitebias
   CHARACTER*100, INTENT(out)        :: FileInput0bias
   CHARACTER*100, INTENT(out)        :: FileInputFinitebias
   CHARACTER*100, INTENT(out)        :: FileInputxyz
   CHARACTER*100, INTENT(out)        :: FileInputWF
   CHARACTER*100, INTENT(out)        :: FileOutputGaussian
   CHARACTER*100, INTENT(out)        :: FileOutputGnuplotx
   CHARACTER*100, INTENT(out)        :: FileOutputGnuploty
   CHARACTER*100, INTENT(out)        :: FileOutputDat
   CHARACTER*100, INTENT(out)        :: FileOutputDat0bias
   CHARACTER*100, INTENT(out)        :: FileOutputDatFinitebias

   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Data File:"
        OPEN(iufile,file=trim(datafile),form='formatted')
        REWIND(iufile)
        READ(iufile,*) FileInputGrid0bias
        PRINT*, "...READING Incoming Data File: Input file for Grid 0 bias      :", Trim(FileInputGrid0bias)
        READ(iufile,*) FileInputGridFinitebias
        PRINT*, "...READING Incoming Data File: Input file for Grid Finite bias :", Trim(FileInputGridFinitebias)
        READ(iufile,*) FileInput0bias
        PRINT*, "...READING Incoming Data File: Input file for 0 bias           :", Trim(FileInput0bias)
        READ(iufile,*) FileInputFinitebias
        PRINT*, "...READING Incoming Data File: Input file for Finite bias      :", Trim(FileInputFinitebias)
        READ(iufile,*) FileInputxyz
        PRINT*, "...READING Incoming Data File: Output file for xyz              :", Trim(FileInputxyz)
        READ(iufile,*) FileInputWF
        PRINT*, "...READING Incoming Data File: Output file for WF               :", Trim(FileInputWF)
        READ(iufile,*)  FileOutputGaussian 
        PRINT*, "...READING Incoming Data File: Output file for Output Gaussian  :", Trim( FileOutputGaussian )
        READ(iufile,*)  FileOutputGnuplotx 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage x :", Trim(FileOutputGnuplotx)
        READ(iufile,*)  FileOutputGnuploty 
        PRINT*, "...READING Incoming Data File: Output file for GnuplotAverage y :", Trim(FileOutputGnuploty)
        READ(iufile,*)  FileOutputDat
        PRINT*, "...READING Incoming Data File: Output file for Average xy :", Trim(FileOutputDat)
        READ(iufile,*)  FileOutputDat0bias
        PRINT*, "...READING Incoming Data File: Output file for Average xy 0 bias      :", Trim(FileOutputDat0bias)
        READ(iufile,*)  FileOutputDatFinitebias
        PRINT*, "...READING Incoming Data File: Output file for Average xy Finite bias :", Trim(FileOutputDatFinitebias)

        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE ReadInputFileWFinPotential

!***********************************************
      SUBROUTINE GaussianFilewrite(FileInput, iufile,Cutoffdistance, dimensions, Nb_atoms, grid, X0Cell, cell,atoms_type,atoms_position, wavefunction)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok

   IMPLICIT NONE
   REAL(dbl), PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   REAL(dbl), PARAMETER ::  Ryd_to_eV = 13.605698066*ONE
   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   REAL(dbl), INTENT(in)               :: Cutoffdistance
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: Nb_atoms
   INTEGER,   INTENT(in)            :: grid(dimensions)
   REAL(dbl),   INTENT(in)         :: X0Cell(dimensions)            ! position of the origin of the volumetric data. 
   REAL(dbl),   INTENT(in)         :: cell(dimensions,dimensions)  ! cell vectors
   INTEGER,     INTENT(in)         :: atoms_type(Nb_atoms)     !atoms_type( Nb_atoms )
   REAL(dbl),   INTENT(in)         :: atoms_position(Nb_atoms,dimensions)
   REAL(dbl),   INTENT(in)         :: wavefunction(grid(1),grid(2),grid(3)) !wavefunction( grid(1),grid(2),grid(3) )

! Local/Reading1
  CHARACTER*100  :: chr
  REAL(dbl)      :: real_io
! Local/Loop
  INTEGER        :: iatom, ix, iy, iz, ia, ib, ic, idimension, imax
! Local/Calculating
  REAL(dbl)         :: dV, norm, R(dimensions)
  REAL(dbl)         :: VolCel, distancetoat, distancemin
  REAL(dbl)         :: Crossaux(dimensions), aux(dimensions)



       PRINT*, "    ...PRITING Output File"
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        WRITE(iufile,*) TRIM(FileInput)                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        WRITE(iufile,*) TRIM(FileInput)                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        WRITE(iufile,'(I5,3F12.6)') Nb_atoms, X0Cell(1)*Ang_to_Bohr,  X0Cell(2)*Ang_to_Bohr,  X0Cell(3)*Ang_to_Bohr !    60  -40.000000  -30.000000  -30.000000
        WRITE(iufile,'(I5,3F12.6)') grid(1), cell(1,1)*Ang_to_Bohr, cell(1,2)*Ang_to_Bohr, cell(1,3)*Ang_to_Bohr    !  300    0.267559    0.000000    0.000000
        WRITE(iufile,'(I5,3F12.6)') grid(2), cell(2,1)*Ang_to_Bohr, cell(2,2)*Ang_to_Bohr, cell(2,3)*Ang_to_Bohr    !  300    0.000000    0.200669    0.000000
        WRITE(iufile,'(I5,3F12.6)') grid(3), cell(3,1)*Ang_to_Bohr, cell(3,2)*Ang_to_Bohr, cell(3,3)*Ang_to_Bohr    !   50    0.000000    0.000000    1.632653
        Crossaux(1) = (grid(1)*cell(1,2)*grid(2)*cell(2,3))  -(grid(2)*cell(2,2)* grid(1)*cell(1,3) )
        Crossaux(2) = (grid(1)*cell(1,3)*grid(2)*cell(2,1))  -(grid(2)*cell(2,3)*grid(1)*cell(1,1))
        Crossaux(3) = (grid(1)*cell(1,1)*grid(2)*cell(2,2))  - (grid(2)*cell(2,1)*grid(1)*cell(1,2))
        aux(1) = grid(3)*cell(3,1)
        aux(2) = grid(3)*cell(3,2)
        aux(3) = grid(3)*cell(3,3)
        VolCel = dot_product(Crossaux(:), aux(:))
        dV = VolCel / (grid(1)*grid(2)*grid(3))



    DO iatom = 1,Nb_atoms
        WRITE(iufile,'(I5,4F12.6)')   atoms_type(iatom), real_io, atoms_position(iatom,1)*Ang_to_Bohr, atoms_position(iatom,2)*Ang_to_Bohr, atoms_position(iatom,3)*Ang_to_Bohr  ! 6    0.000000    2.924528    7.785081    1.274146
    ENDDO
    DO ia = 1,  grid(1)
    DO ib = 1,  grid(2)
    DO ic = 1,  grid(3), 6
        imax=MIN(5,(grid(3)-ic)) 
        distancetoat = ONE*VolCel 
        distancemin  = ONE*VolCel 
        DO iz = 1, imax+1        
            DO idimension=1,3
                 R(idimension) = (ia-1)*cell(1,idimension)+ (ib-1)*cell(2,idimension) + (ic-2+iz)*cell(3,idimension) +X0Cell(idimension)   
            ENDDO
            DO iatom = 1, Nb_atoms
                distancetoat = sqrt( ( R(1) - atoms_position(iatom,1) )**2 + ( R(2) -atoms_position(iatom,2))**2 + ( R(3) -atoms_position(iatom,3))**2 ) 
                IF (distancetoat < distancemin) distancemin = distancetoat
            ENDDO
        ENDDO
        IF ( distancemin <= Cutoffdistance ) THEN
           WRITE(iufile,'(6 ES13.5)')   wavefunction(ia,ib,ic:ic+imax) 
        ELSE
!           PRINT*, "Skipped", R(:), distancemin
           WRITE(iufile,'(6 ES13.5)')   (ZERO, iz=1,imax+1) 
        ENDIF 
    ENDDO     
    ENDDO
    ENDDO
        PRINT*, "            ...Done"
        PRINT*, "    ...CLOSING Output File"
        CLOSE(iufile)


END  SUBROUTINE GaussianFilewrite
!***********************************************
      SUBROUTINE GaussianFilereadHeader(FileInput, iufile,Cutoffdistance, dimensions, Nb_atoms, grid)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE
   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   REAL(dbl), INTENT(in)            :: Cutoffdistance
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER , INTENT(out)            :: grid(dimensions)
   INTEGER , INTENT(out)            :: Nb_atoms
   REAL(dbl)                        :: cell(dimensions,dimensions)  ! cell vectors
   REAL(dbl)                        :: X0Cell(dimensions)            ! position of the origin of the volumetric data.
   CHARACTER*100  :: chr

        PRINT*, "...READING Incoming Cube File: Header"
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        READ(iufile,*) Nb_atoms, X0Cell(1),  X0Cell(2),  X0Cell(3) !    60  -40.000000  -30.000000  -30.000000
        PRINT'(A25,I5)', "    ...Number of atoms: ", Nb_atoms
        PRINT'(A25,3F12.6)', "    ...Origin         : ", X0Cell(:)
        READ(iufile,*) grid(1), cell(1,1), cell(1,2), cell(1,3)    !  300    0.267559    0.000000    0.000000
        READ(iufile,*) grid(2), cell(2,1), cell(2,2), cell(2,3)    !  300    0.000000    0.200669    0.000000
        READ(iufile,*) grid(3), cell(3,1), cell(3,2), cell(3,3)    !   50    0.000000    0.000000    1.632653
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid    : ",  grid(1), " x " ,  grid(2), " x ",  grid(3)


        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE GaussianFilereadHeader
!***********************************************
      SUBROUTINE GaussianFileread(FileInput, iufile,Cutoffdistance, dimensions, Nb_atoms_check, grid_check, X0Cell, cell,atoms_type,atoms_position, wavefunction, norm)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE
   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   REAL(dbl), INTENT(in)               :: Cutoffdistance
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: Nb_atoms_check
   INTEGER,   INTENT(in)            :: grid_check(dimensions)
   !
   REAL(dbl),   INTENT(out)         :: X0Cell(dimensions)            ! position of the origin of the volumetric data. 
   REAL(dbl),   INTENT(out)         :: cell(dimensions,dimensions)  ! cell vectors
   INTEGER,   INTENT(out)           :: atoms_type(Nb_atoms_check)     !atoms_type( Nb_atoms )
   REAL(dbl),   INTENT(out)         ::  atoms_position(Nb_atoms_check,dimensions)
   REAL(dbl),   INTENT(out)         :: wavefunction(grid_check(1),grid_check(2),grid_check(3)) !wavefunction( grid(1),grid(2),grid(3) )
   REAL(dbl),   INTENT(out)         :: norm

 ! Input Check
 INTEGER               :: Nb_atoms  ! Number of atoms
 INTEGER               :: grid(dimensions)

! Local Variables

  REAL(dbl)                :: barycenter(dimensions)  ! barycenter
  REAL(dbl)                :: CenterOfMass(dimensions)  ! barycenter


! Local/Reading1
  CHARACTER*100  :: chr
  REAL(dbl)      :: real_io
! Local/Loop
  INTEGER        :: iatom, ix, iy, iz, ia, ib, ic, idimension, imax
! Local/Calculating
  REAL(dbl)         :: dV
  REAL(dbl)         :: VolCel, distancetoat, distancemin, Volume
  REAL(dbl)         :: Crossaux(dimensions), aux(dimensions)
  REAL(dbl)         :: TotalMass

!***********************************************
!               MAIN BODY
!***********************************************
 ! READ
        PRINT*, "...READING Incoming Cube File"
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        read(iufile,*)chr                                          !  JibinMolIsolatedconfig2_WF.WF103.cube                       
        READ(iufile,*) Nb_atoms, X0Cell(1),  X0Cell(2),  X0Cell(3) !    60  -40.000000  -30.000000  -30.000000
        PRINT'(A25,I5)', "    ...Number of atoms: ", Nb_atoms
        PRINT'(A25,I5)', "    ...Number of atoms refence: ", Nb_atoms_check
        IF (Nb_atoms /=  Nb_atoms_check) STOP
        PRINT'(A25,3F12.6)', "    ...Origin         : ", X0Cell(:)
        READ(iufile,*) grid(1), cell(1,1), cell(1,2), cell(1,3)    !  300    0.267559    0.000000    0.000000
        READ(iufile,*) grid(2), cell(2,1), cell(2,2), cell(2,3)    !  300    0.000000    0.200669    0.000000
        READ(iufile,*) grid(3), cell(3,1), cell(3,2), cell(3,3)    !   50    0.000000    0.000000    1.632653
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid    : ",  grid(1), " x " ,  grid(2), " x ",  grid(3)
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid Ref: ",  grid_check(1), " x " ,  grid_check(2), " x ",  grid_check(3)
        DO idimension = 1,dimensions
           IF (grid_check(idimension) /= grid(idimension)) STOP
        ENDDO
        PRINT'(A25,F12.6,A3,F12.6,A3,F12.6)', "    ...A vector [Bohr]: ",  grid(1)*cell(1,1), " | " ,  grid(1)*cell(1,2), " | " ,  grid(1)*cell(1,3) 
        PRINT'(A25,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector [Bohr]: ",  grid(2)*cell(2,1), " | " ,  grid(2)*cell(2,2), " | " ,  grid(2)*cell(2,3) 
        PRINT'(A25,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector [Bohr]: ",  grid(3)*cell(3,1), " | " ,  grid(3)*cell(3,2), " | " ,  grid(3)*cell(3,3) 
        Crossaux(1) = (grid(1)*cell(1,2)*grid(2)*cell(2,3))  -(grid(2)*cell(2,2)* grid(1)*cell(1,3) )
        Crossaux(2) = (grid(1)*cell(1,3)*grid(2)*cell(2,1))  -(grid(2)*cell(2,3)*grid(1)*cell(1,1))
        Crossaux(3) = (grid(1)*cell(1,1)*grid(2)*cell(2,2))  - (grid(2)*cell(2,1)*grid(1)*cell(1,2))
        aux(1) = grid(3)*cell(3,1)
        aux(2) = grid(3)*cell(3,2)
        aux(3) = grid(3)*cell(3,3)
        VolCel = dot_product(Crossaux(:), aux(:))
        dV = VolCel / (grid(1)*grid(2)*grid(3))
        PRINT'(A35,F12.4)', "    ...Volume of the cell [Bohr3]: ", VolCel
        PRINT'(A35,F12.6)', "    ...Voxel [Bohr3]             : ", dV
        PRINT*, "    ...READING Atoms type and positions"
        barycenter(:)=ZERO

        CenterOfMass(:) = ZERO
        TotalMass= ZERO

    DO iatom = 1, Nb_atoms
        atoms_type(iatom)=1111
        atoms_position(iatom,:)=10000*ONE 
        READ(iufile,*)   atoms_type(iatom), real_io, atoms_position(iatom,1), atoms_position(iatom,2), atoms_position(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146
        PRINT'(A18,I5,A7,I5,A6,F12.6,A1,F12.6,A1,F12.6,A2)', "        ...Atom:", iatom,  " Type: ",atoms_type(iatom)," at (",atoms_position(iatom,1),",",atoms_position(iatom,2),",",atoms_position(iatom,3)," )"
        barycenter(1)=        barycenter(1) + atoms_position(iatom,1)
        barycenter(2)=        barycenter(2) + atoms_position(iatom,2)
        barycenter(3)=        barycenter(3) + atoms_position(iatom,3)
        CenterOfMass(1) = CenterOfMass(1) + ( atoms_position(iatom,1)* atoms_type(iatom)*2 ) 
        CenterOfMass(2) = CenterOfMass(2) + ( atoms_position(iatom,2)* atoms_type(iatom)*2 ) 
        CenterOfMass(3) = CenterOfMass(3) + ( atoms_position(iatom,3)* atoms_type(iatom)*2 ) 
        TotalMass= TotalMass + atoms_type(iatom)*2
     ENDDO
        barycenter(:) = barycenter(:) / Nb_atoms
        CenterOfMass(:) =CenterOfMass(:) / TotalMass

        PRINT'(A45,F12.6,A1,F12.6,A1,F12.6,A2)', "    ...Barycenter of the system at [Bohr]: (",barycenter(1),",",barycenter(2),",",barycenter(3), " )"
        PRINT'(A45,F12.6,A1,F12.6,A1,F12.6,A2)', "    ...Center of Mass [Bohr]             : (", CenterOfMass(1),",", CenterOfMass(2),",", CenterOfMass(3), " )"
        PRINT*, "    ...READING Wavefunction"
    norm=ZERO
    Volume= ZERO
    DO ix = 1,  grid(1)
    DO iy = 1,  grid(2)
    DO iz = 1,  grid(3), 6 
        imax=MIN(5,(grid(3)-iz))
!        PRINT*, ix, iy, iz, iz+imax
        READ(iufile,*)   wavefunction(ix,iy,iz:iz+imax) 
    ENDDO     
    ENDDO
    ENDDO
    DO ix = 1,  grid(1)
    DO iy = 1,  grid(2)
    DO iz = 1,  grid(3) 
            norm = norm +  wavefunction(ix,iy,iz)*wavefunction(ix,iy,iz) * dV
            Volume=Volume+dV
    ENDDO     
    ENDDO
    ENDDO

    norm = sqrt(norm)
        PRINT*, "            ...Done"
        PRINT'(A35,F12.4)', "    ...Integrated Volume [Bohr3] : ", Volume
        PRINT'(A35,F12.4)', "    ...Volume of the cell [Bohr3]: ", VolCel
        PRINT*, "            ...Norm:  ", norm
        PRINT*, "    ...CLOSING Gaussian Input File"
        CLOSE(iufile)


END  SUBROUTINE GaussianFileread
!***********************************************
      SUBROUTINE RhoFilereadHeader(FileInput, iufile, Cutoffdistance, dimensions, grid)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok

   IMPLICIT NONE
   REAL(dbl), PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   REAL(dbl), PARAMETER ::  Ryd_to_eV = 13.605698066*ONE
   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   REAL(dbl), INTENT(in)            :: Cutoffdistance
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER , INTENT(out)            :: grid(dimensions)
!   INTEGER , INTENT(out)            :: Nb_atoms
   REAL(dbl)                        :: cell(dimensions,dimensions)  ! cell vectors
!   REAL(dbl)                        :: X0Cell(dimensions)            ! position of the origin of the volumetric data.
   CHARACTER*100  :: chr
   INTEGER :: idim1, idim2, ns
        PRINT*, "...READING Incoming .RHO File: Header ", TRIM(FileInput)
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        DO idim1=1, dimensions 
              READ(iufile,*) (cell(idim2,idim1),idim2=1,dimensions)
        ENDDO
        cell(:,:) = cell(:,:)/ Ang_to_Bohr
        READ(iufile,*) grid(1), grid(2),   grid(3),   ns
        PRINT*, "...Grid", grid(1), grid(2),   grid(3)
        CLOSE(iufile)
        PRINT*, "...Done"


END  SUBROUTINE RhoFilereadHeader
!***********************************************
      SUBROUTINE RhoFileread(FileInput, iu,Cutoffdistance, dimensions,  grid_check, X0Cell, cell, wavefunction, norm)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok

   IMPLICIT NONE
   REAL(dbl), PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   REAL(dbl), PARAMETER ::  Ryd_to_eV = 13.605698066*ONE
   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iu
   REAL(dbl), INTENT(in)            :: Cutoffdistance
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: grid_check(dimensions)
   !
   REAL(dbl),   INTENT(out)         :: X0Cell(dimensions)           ! position of the origin of the volumetric data. 
   REAL(dbl),   INTENT(out)         :: cell(dimensions,dimensions)  ! cell vectors
   REAL(dbl),   INTENT(out)         :: wavefunction(grid_check(1),grid_check(2),grid_check(3)) !wavefunction( grid(1),grid(2),grid(3) )
   REAL(dbl),   INTENT(out)         :: norm
 ! Input Check
 INTEGER               :: grid(dimensions)

! Local Variables
! Local/Reading1
  CHARACTER*100  :: chr
  REAL(dbl)      :: real_io
! Local/Loop
  INTEGER        :: iatom, ix, iy, iz, ia, ib, ic, idimension, imax, ns, idim1, idim2
! Local/Calculating
  REAL(dbl)         :: dV
  REAL(dbl)         :: VolCel
  REAL(dbl)         :: Crossaux(dimensions), aux(dimensions)


!***********************************************
!               MAIN BODY
!***********************************************

    X0Cell(:) = ZERO

        PRINT*, "...READING Incoming .RHO File: ", TRIM(FileInput)
        OPEN(iu,file=trim(FileInput),form='formatted')
        REWIND(iu)
        DO idim1=1, dimensions 
              READ(iu,*) (cell(idim2,idim1),idim2=1,dimensions)
        ENDDO
        cell(:,:) = cell(:,:)/ Ang_to_Bohr
        READ(iu,*) grid(1), grid(2),   grid(3),   ns
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid    : ",  grid(1), " x " ,  grid(2), " x ",  grid(3)
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid Ref: ",  grid_check(1), " x " ,  grid_check(2), " x ",  grid_check(3)
        DO idimension = 1,dimensions
           IF (grid_check(idimension) /= grid(idimension)) STOP
        ENDDO
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...A vector: ",  grid(1)*cell(1,1), " | " ,  grid(1)*cell(1,2), " | " ,  grid(1)*cell(1,3) 
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector: ",  grid(2)*cell(2,1), " | " ,  grid(2)*cell(2,2), " | " ,  grid(2)*cell(2,3) 
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector: ",  grid(3)*cell(3,1), " | " ,  grid(3)*cell(3,2), " | " ,  grid(3)*cell(3,3) 
!        Crossaux(1) = (grid(1)*cell(1,2)*grid(2)*cell(2,3))  -(grid(2)*cell(2,2)* grid(1)*cell(1,3) )
!        Crossaux(2) = (grid(1)*cell(1,3)*grid(2)*cell(2,1))  -(grid(2)*cell(2,3)*grid(1)*cell(1,1))
!        Crossaux(3) = (grid(1)*cell(1,1)*grid(2)*cell(2,2))  - (grid(2)*cell(2,1)*grid(1)*cell(1,2))
!        aux(1) = grid(3)*cell(3,1)
!        aux(2) = grid(3)*cell(3,2)
!        aux(3) = grid(3)*cell(3,3)
!        VolCel = dot_product(Crossaux(:), aux(:))
!        dV = VolCel / (grid(1)*grid(2)*grid(3))
!        PRINT'(A28,F12.4)', "    ...Volume of the cell: ", VolCel
!        PRINT'(A28,F12.6)', "    ...Voxel:              ", dV

        norm = ZERO
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...A vector: ",  cell(1,1), " | " , cell(1,2), " | " ,  cell(1,3) 
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector: ",  cell(2,1), " | " ,  cell(2,2), " | " ,  cell(2,3) 
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector: ",  cell(3,1), " | " ,  cell(3,2), " | " ,  cell(3,3) 
        Crossaux(1) = (cell(1,2)*cell(2,3))  -(cell(2,2)*cell(1,3) )
        Crossaux(2) = (cell(1,3)*cell(2,1))  -(cell(2,3)*cell(1,1))
        Crossaux(3) = (cell(1,1)*cell(2,2))  - (cell(2,1)*cell(1,2))
        aux(1) = cell(3,1)
        aux(2) = cell(3,2)
        aux(3) = cell(3,3)
        VolCel = dot_product(Crossaux(:), aux(:))
        dV = VolCel / (grid(1)*grid(2)*grid(3))
        PRINT'(A28,F12.4)', "    ...Volume of the cell: ", VolCel
        PRINT'(A28,F12.6)', "    ...Voxel:              ", dV
        PRINT*, "    ...READING Wavefunction"
        DO iz = 1,  grid(3)
   	 DO iy = 1,  grid(2)
          DO ix = 1,  grid(1)
            READ(iu,*)   wavefunction(ix,iy,iz) 
            norm = norm +  wavefunction(ix,iy,iz) * dV
          ENDDO     
         ENDDO
        ENDDO
        PRINT*, "            ...Done"
        PRINT*, "            ...Norm:  ", norm
        PRINT*, "    ...CLOSING Rho File"
        CLOSE(iu)
        PRINT*, "...Done"
   
        cell(1,:) = cell(1,:) / (grid(1))
        cell(2,:) = cell(2,:) / (grid(2))
        cell(3,:) = cell(3,:) / (grid(3))
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector: ",  grid(2)*cell(2,1), " | " ,  grid(2)*cell(2,2), " | " ,  grid(2)*cell(2,3) 
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector: ",  grid(3)*cell(3,1), " | " ,  grid(3)*cell(3,2), " | " ,  grid(3)*cell(3,3) 
       
 
END  SUBROUTINE RhoFileread

!***********************************************
      SUBROUTINE PotentialFilereadHeader(FileInput, iufile, Cutoffdistance, dimensions, grid, X0Cell, cell)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok

   IMPLICIT NONE
   REAL(dbl), PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   REAL(dbl), PARAMETER ::  Ryd_to_eV = 13.605698066*ONE
   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   REAL(dbl), INTENT(in)            :: Cutoffdistance
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER , INTENT(out)            :: grid(dimensions)
!   INTEGER , INTENT(out)            :: Nb_atoms
   REAL(dbl), INTENT(out)           :: cell(dimensions,dimensions)  ! cell vectors
   REAL(dbl), INTENT(out)           :: X0Cell(dimensions)            ! position of the origin of the volumetric data.
   CHARACTER*100  :: chr
   INTEGER :: idim1, idim2, ns
        PRINT*, "...READING Incoming .RHO File: Header  ", TRIM(FileInput)
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        DO idim1=1, dimensions 
              READ(iufile,*) (cell(idim2,idim1),idim2=1,dimensions)
        ENDDO
        cell(:,:) = cell(:,:)/ Ang_to_Bohr
        READ(iufile,*) grid(1), grid(2),   grid(3),   ns
        PRINT*, "...Grid", grid(1), grid(2),   grid(3)
        CLOSE(iufile)
        PRINT*, "...Done"

        X0Cell(:) = ZERO 
        cell(1,:) = cell(1,:) / (grid(1)-1)
        cell(2,:) = cell(2,:) / (grid(2)-1)
        cell(3,:) = cell(3,:) / (grid(3)-1)


END  SUBROUTINE PotentialFilereadHeader

!***********************************************
      SUBROUTINE PotentialFileread(FileInput, iu,Cutoffdistance, dimensions,  grid_check, wavefunction)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok

   IMPLICIT NONE
   REAL(dbl), PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   REAL(dbl), PARAMETER ::  Ryd_to_eV = 13.605698066*ONE
   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iu
   REAL(dbl), INTENT(in)            :: Cutoffdistance
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: grid_check(dimensions)
   !
   REAL(dbl),   INTENT(out)         :: wavefunction(grid_check(1),grid_check(2),grid_check(3)) !wavefunction( grid(1),grid(2),grid(3) )

 ! Input Check
 INTEGER               :: grid(dimensions)

! Local Variables
! Local/Reading1
  CHARACTER*100  :: chr
  CHARACTER*2    :: chr2
  REAL(dbl)      :: real_io
! Local/Loop
  INTEGER        :: iatom, ix, iy, iz, ia, ib, ic, idimension, imax, ns, idim1, idim2
! Local/Calculating
  REAL(dbl)         :: Crossaux(dimensions), aux(dimensions)


!***********************************************
!               MAIN BODY
!***********************************************



        PRINT*, "...READING Incoming Potential File:", TRIM(FileInput)
        OPEN(iu,file=trim(FileInput),form='formatted')
        REWIND(iu)
        READ(iu,*) ! Empty line
        READ(iu,*) chr2, grid(1), grid(2),   grid(3)
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid    : ",  grid(1), " x " ,  grid(2), " x ",  grid(3)
        PRINT'(A17,I5,A3,I5,A3,I5)', "    ...Grid Ref: ",  grid_check(1), " x " ,  grid_check(2), " x ",  grid_check(3)
        DO idimension = 1,dimensions
           IF (grid_check(idimension) /= grid(idimension)) STOP
        ENDDO
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...A vector: ",  grid(1)*cell(1,1), " | " ,  grid(1)*cell(1,2), " | " ,  grid(1)*cell(1,3) 
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector: ",  grid(2)*cell(2,1), " | " ,  grid(2)*cell(2,2), " | " ,  grid(2)*cell(2,3) 
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector: ",  grid(3)*cell(3,1), " | " ,  grid(3)*cell(3,2), " | " ,  grid(3)*cell(3,3) 
!        Crossaux(1) = (grid(1)*cell(1,2)*grid(2)*cell(2,3))  -(grid(2)*cell(2,2)* grid(1)*cell(1,3) )
!        Crossaux(2) = (grid(1)*cell(1,3)*grid(2)*cell(2,1))  -(grid(2)*cell(2,3)*grid(1)*cell(1,1))
!        Crossaux(3) = (grid(1)*cell(1,1)*grid(2)*cell(2,2))  - (grid(2)*cell(2,1)*grid(1)*cell(1,2))
!        aux(1) = grid(3)*cell(3,1)
!        aux(2) = grid(3)*cell(3,2)
!        aux(3) = grid(3)*cell(3,3)
!        VolCel = dot_product(Crossaux(:), aux(:))
!        dV = VolCel / (grid(1)*grid(2)*grid(3))
!        PRINT'(A28,F12.4)', "    ...Volume of the cell: ", VolCel
!        PRINT'(A28,F12.6)', "    ...Voxel:              ", dV

        PRINT*, "    ...READING Potential"
        DO iz = 1,  grid(3)
   	 DO iy = 1,  grid(2)
          DO ix = 1,  grid(1)
            READ(iu,*)   wavefunction(ix,iy,iz) 
            IF((wavefunction(ix,iy,iz)>1000000*ONE).OR.(wavefunction(ix,iy,iz)<-1000000*ONE)) PRINT*, "WARNING: bad value", ix,iy,iz
          ENDDO     
         ENDDO
        ENDDO
        PRINT*, "            ...Done"
        PRINT*, "    ...CLOSING Potential File", TRIM(FileInput)
        CLOSE(iu)
        PRINT*, "...Done"
        wavefunction(:,:,:) = wavefunction(:,:,:)*Ryd_to_eV
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector: ",  grid(2)*cell(2,1), " | " ,  grid(2)*cell(2,2), " | " ,  grid(2)*cell(2,3) 
!        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector: ",  grid(3)*cell(3,1), " | " ,  grid(3)*cell(3,2), " | " ,  grid(3)*cell(3,3) 
       
 
END  SUBROUTINE PotentialFileread

!***********************************************
      SUBROUTINE XYZFilereadHeader(FileInput, iufile, Nb_atoms)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   INTEGER,   INTENT(out)            :: Nb_atoms
   !
! Local Variables
 CHARACTER*100         :: chr


        PRINT*, "...READING Incoming XYZ File: "
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*) Nb_atoms               
        PRINT'(A25,I5)', "    ...Number of atoms: ", Nb_atoms
 
        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE XYZFilereadHeader

!***********************************************
      SUBROUTINE XYZFileread(FileInput, iufile, dimensions, Nb_atoms_check, atoms_type,atoms_position)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   CHARACTER*100, INTENT(in)        :: FileInput
   INTEGER, INTENT(in)              :: iufile
   INTEGER,   INTENT(in)            :: dimensions
   INTEGER,   INTENT(in)            :: Nb_atoms_check
   !
   INTEGER,   INTENT(out)           :: atoms_type(Nb_atoms_check)     !atoms_type( Nb_atoms )
   REAL(dbl),   INTENT(out)         :: atoms_position(Nb_atoms_check,dimensions)
! Local Variables
 INTEGER               :: Nb_atoms, iatom  ! Number of atoms
 CHARACTER*100         :: chr


        PRINT*, "...READING Incoming XYZ File: "
        OPEN(iufile,file=trim(FileInput),form='formatted')
        REWIND(iufile)
        read(iufile,*) Nb_atoms               
        PRINT'(A25,I5)', "    ...Number of atoms: ", Nb_atoms
        PRINT'(A25,I5)', "    ...Number of atoms refence: ", Nb_atoms_check
        IF (Nb_atoms /=  Nb_atoms_check) STOP
                                       
        read(iufile,*)                                                            
        DO iatom = 1, Nb_atoms
           READ(iufile,*)   chr, atoms_position(iatom,1), atoms_position(iatom,2), atoms_position(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146 

           IF ( TRIM(chr) == 'H' ) THEN 
               atoms_type(iatom) = 1
           ELSE IF ( TRIM(chr) == 'He' ) THEN 
               atoms_type(iatom) = 2
           ELSE IF ( TRIM(chr) == 'Li' ) THEN 
               atoms_type(iatom) = 3
           ELSE IF ( TRIM(chr) == 'Be' ) THEN 
               atoms_type(iatom) = 4
           ELSE IF ( TRIM(chr) == 'B' ) THEN 
               atoms_type(iatom) = 5
           ELSE IF ( TRIM(chr) == 'C' ) THEN 
               atoms_type(iatom) = 6
           ELSE IF ( TRIM(chr) == 'N' ) THEN 
               atoms_type(iatom) = 7
           ELSE IF ( TRIM(chr) == 'O' ) THEN 
               atoms_type(iatom) = 8
           ELSE IF ( TRIM(chr) == 'F' ) THEN 
               atoms_type(iatom) = 9
           ELSE IF ( TRIM(chr) == 'Ne' ) THEN 
               atoms_type(iatom) = 10
           ELSE IF ( TRIM(chr) == 'Na' ) THEN 
               atoms_type(iatom) = 11
           ELSE IF ( TRIM(chr) == 'Mg' ) THEN 
               atoms_type(iatom) = 12
           ELSE IF ( TRIM(chr) == 'Al' ) THEN 
               atoms_type(iatom) = 13
           ELSE IF ( TRIM(chr) == 'Si' ) THEN 
               atoms_type(iatom) = 14
           ELSE IF ( TRIM(chr) == 'P' ) THEN 
               atoms_type(iatom) = 15
           ELSE IF ( TRIM(chr) == 'S' ) THEN 
               atoms_type(iatom) = 16
           ELSE IF ( TRIM(chr) == 'Cl' ) THEN 
               atoms_type(iatom) = 17
           ELSE IF ( TRIM(chr) == 'Ar' ) THEN 
               atoms_type(iatom) = 18
           ELSE IF ( TRIM(chr) == 'Au' ) THEN 
               atoms_type(iatom) = 79
           ELSE 
               PRINT*, "Warning, not identified"
               atoms_type(iatom) = 10000
           ENDIF
           PRINT*, "Atom",   TRIM(chr), iatom
           PRINT*,  atoms_type(iatom), atoms_position(iatom,1), atoms_position(iatom,2), atoms_position(iatom,3)  ! 6    0.000000    2.924528    7.785081    1.274146     
        ENDDO
        CLOSE(iufile)
        PRINT*, "...Done"

END  SUBROUTINE XYZFileread

END MODULE IO_Module
