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
   PROGRAM PlotPotential
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   USE IO_module,                              ONLY : gnuplotwrite, datwrite, ReadInputFilePlotPotential, GaussianFilereadHeader, GaussianFilewrite, GaussianFileread, PotentialFilereadHeader,  RhoFileread, XYZFilereadHeader, XYZFileread
   IMPLICIT NONE
   INTEGER, PARAMETER :: dimensions = 3    
   INTEGER, PARAMETER :: iufile=100
   REAL*8, PARAMETER ::  Cutoffdistance=20004.5
   REAL*8, PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   CHARACTER*100, PARAMETER :: datafile="DATA.in"


   !
   ! Input variables
   !
   CHARACTER*100 ::  FileInput0bias
   CHARACTER*100 ::  FileInputGrid0bias
   CHARACTER*100 ::  FileInputxyz
   CHARACTER*100 ::  FileOutputGaussian0bias="rho0.cube"
   CHARACTER*100 ::  FileOutputGnuplotx
   CHARACTER*100 ::  FileOutputGnuploty
   CHARACTER*100 ::  FileOutputDat0bias

   REAL*8 ::   cell(dimensions,dimensions), norm
   REAL*8 ::   X0Cell(dimensions)
   INTEGER ::   mesh(dimensions)
   INTEGER ::   mesh_test(dimensions)
   INTEGER ::   Nb_atoms

   REAL*8, ALLOCATABLE ::   density_0bias(:,:,:)

   INTEGER, ALLOCATABLE ::   atoms_type(:)
   REAL*8, ALLOCATABLE ::   atoms_position(:,:)
   !
   ! Output variables
   !

   REAL*8, ALLOCATABLE ::   Integrated_density1(:,:)
   REAL*8, ALLOCATABLE ::   Integrated_density2(:,:)
   REAL*8, ALLOCATABLE ::   Integrated_density_xy(:) 
   ! 
   ! Local variables
   ! 
   INTEGER :: ia, ib, ic
   REAL*8, ALLOCATABLE  :: vect1(:), vect2(:)
!***********************************************
!
!------------------------------
! main body
!------------------------------
!  !
   ! read input General Information
   !      SUBROUTINE ReadInputFilePlot(iufile, datafile, FileInput, FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat)    
   CALL ReadInputFilePlotPotential(iufile, datafile, FileInput0bias,FileInputGrid0bias, FileInputxyz , FileOutputGaussian0bias, FileOutputGnuplotx , FileOutputGnuploty, FileOutputDat0bias)



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
   !CALL RhoFilereadHeader( FileInput0bias, iufile, Cutoffdistance, dimensions, mesh)

   CALL PotentialFilereadHeader( FileInputGrid0bias, iufile, Cutoffdistance, dimensions, mesh, X0Cell, cell)

   ALLOCATE( density_0bias(mesh(1),mesh(2),mesh(3)) )   
   density_0bias(:,:,:) = ZERO 
   CALL RhoFileread( FileInput0bias, iufile,Cutoffdistance, dimensions,  mesh, X0Cell, cell, density_0bias, norm)
   !CALL PotentialFileread( FileInput0bias, iufile,Cutoffdistance, dimensions,  mesh, density_0bias)
   density_0bias(:,:,:) =    density_0bias(:,:,:)
   CALL GaussianFilewrite(FileOutputGaussian0bias, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, density_0bias)


!  !
   !Difference  Average   on x and y lines
   !
   ALLOCATE( Integrated_density1(mesh(2),mesh(3)) )
   ALLOCATE( Integrated_density2(mesh(1),mesh(3)) )
   Integrated_density1(:,:) = ZERO
   Integrated_density2(:,:) = ZERO
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_density1(ib,ic) =  Integrated_density1(ib,ic) +  density_0bias(ia,ib,ic) / mesh(1)
		   ENDDO
	   ENDDO
   ENDDO

  DO ic=1,mesh(3)
       	  DO  ia=1,mesh(1)
	           DO ib=1,mesh(2)
		     Integrated_density2(ia,ic) =  Integrated_density2(ia,ic) +  density_0bias(ia,ib,ic) / mesh(2)
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
		     Integrated_density_xy(ic) =  Integrated_density_xy(ic) +  density_0bias(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_density_xy(ic) = Integrated_density_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO

!  !
   ! Print the 0bias file averaged in a .dat file
   !
   CALL datwrite(iufile, FileOutputDat0bias, mesh(3), cell(3,3), Integrated_density_xy(:))

!  !
   ! Deallocate
   !
   PRINT*, "Deallocate"
   DEALLOCATE(density_0bias)
   DEALLOCATE(atoms_type, atoms_position)
   DEALLOCATE(Integrated_density1,Integrated_density2, Integrated_density_xy)
   DEALLOCATE ( vect1,vect2 )



   PRINT*, "END"
END PROGRAM PlotPotential

