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
   PROGRAM PotentialDifference
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   USE IO_module,                              ONLY : gnuplotwrite, datwrite, ReadInputFilePotential, GaussianFilereadHeader, GaussianFilewrite, GaussianFileread, PotentialFilereadHeader,  PotentialFileread,  XYZFilereadHeader, XYZFileread
   IMPLICIT NONE
   INTEGER, PARAMETER :: dimensions = 3    
   CHARACTER*11, PARAMETER :: FileFormat = "unformatted"     
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
   CHARACTER*100 ::  FileOutputGaussian
   CHARACTER*100 ::  FileOutputGaussian0bias="rho0.cube"
   CHARACTER*100 ::  FileOutputGaussianFiniteBias="rhofinite.cube"
   CHARACTER*100 ::  FileOutputGnuplotx
   CHARACTER*100 ::  FileOutputGnuploty
   CHARACTER*100 ::  FileOutputDat
   CHARACTER*100 ::  FileOutputDat0bias
   CHARACTER*100 ::  FileOutputDatFinitebias
   CHARACTER(LEN=30)  :: Format
   REAL(dbl) ::   cell(dimensions,dimensions)
   REAL(dbl) ::   Position_aux1(dimensions)
   REAL(dbl) ::   Position_aux2(dimensions)

   REAL(dbl) ::   X0Cell(dimensions)
   INTEGER ::   mesh(dimensions)
   INTEGER ::   mesh_test(dimensions),vect_aux1(dimensions),vect_aux2(dimensions)
   INTEGER ::   Nb_atoms

   REAL(dbl), ALLOCATABLE ::   density_0bias(:,:,:)
   REAL(dbl), ALLOCATABLE ::   density_Finitebias(:,:,:)
   INTEGER, ALLOCATABLE ::   atoms_type(:)
   REAL(dbl), ALLOCATABLE ::   atoms_position(:,:)
   !
   ! Output variables
   !
   REAL(dbl), ALLOCATABLE ::   target_density(:,:,:)
   REAL(dbl), ALLOCATABLE ::   Integrated_density1(:,:)
   REAL(dbl), ALLOCATABLE ::   Integrated_density2(:,:)
   REAL(dbl), ALLOCATABLE ::   Integrated_density_xy(:) 
   INTEGER, ALLOCATABLE :: Atomsongrid(:,:)
   ! 
   ! Local variables
   ! 
   LOGICAL :: found
   INTEGER :: ia, ib, ic, iatom
   REAL(dbl), ALLOCATABLE  :: vect1(:), vect2(:)
!***********************************************
!
!------------------------------
! main body
!------------------------------
!  !
   ! read input General Information
   !
   CALL ReadInputFilePotential(iufile, datafile, FileInputGrid0bias, FileInputGridFinitebias, FileInput0bias, FileInputFinitebias, FileInputxyz , FileOutputGaussian , FileOutputGnuplotx , FileOutputGnuploty,FileOutputDat,FileOutputDat0bias,FileOutputDatFinitebias)



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
   CALL PotentialFilereadHeader( FileInputGrid0bias, iufile, Cutoffdistance, dimensions, mesh, X0Cell, cell,FileFormat)

   ALLOCATE( density_0bias(mesh(1),mesh(2),mesh(3)) )   
   density_0bias(:,:,:) = ZERO 
   CALL PotentialFileread( FileInput0bias, iufile,Cutoffdistance, dimensions,  mesh, density_0bias,FileFormat)
   density_0bias(:,:,:) =    density_0bias(:,:,:)
   CALL GaussianFilewrite(FileOutputGaussian0bias, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, density_0bias)

!  !
   ! read input file Finite bias
   !
   CALL PotentialFilereadHeader(FileInputGridFinitebias, iufile, Cutoffdistance, dimensions, mesh_test, X0Cell, cell,FileFormat)
   IF ( ( ABS(mesh_test(1)- mesh(1)) >  EPS_m5 ).OR.( ABS(mesh_test(2) - mesh(2)) >  EPS_m5 ).OR.( ABS(mesh_test(3)-mesh(3))>  EPS_m5 ) ) THEN
       PRINT*, "Discrepancies in the meshes" 
       PRINT*, mesh_test(1), mesh(1)
       PRINT*, mesh_test(2), mesh(2) 
       PRINT*, mesh_test(3), mesh(3) 
      STOP
   ENDIF

   ALLOCATE( density_Finitebias(mesh(1),mesh(2),mesh(3)) )   
   density_Finitebias(:,:,:) = ZERO 
   CALL PotentialFileread( FileInputFinitebias, iufile,Cutoffdistance, dimensions,  mesh,  density_Finitebias, FileFormat)

   density_Finitebias(:,:,:) =    density_Finitebias(:,:,:)
   CALL GaussianFilewrite(FileOutputGaussianFinitebias, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, density_Finitebias)


!  !
   ! Calculate the Difference
   !
   ALLOCATE( target_density(mesh(1),mesh(2),mesh(3)) )   
   target_density(:,:,:) = ZERO
   target_density(:,:,:) = density_Finitebias(:,:,:) - density_0bias(:,:,:)


!  !
   ! Output the values at the atomic positions
   !
   ALLOCATE (  Atomsongrid(Nb_atoms, dimensions))
   Atomsongrid(:,:)=1
   DO iatom=1,Nb_atoms
       found=.FALSE.
       icloop: DO ic=1,mesh(3)-1
          ibloop: DO ib=1,mesh(2)-1
            ialoop: DO ia=1,mesh(1)-1
                vect_aux1(1) = ia-1
                vect_aux1(2) = ib-1
                vect_aux1(3) = ic-1
                Position_aux1(:) = matmul( TRANSPOSE (cell(:,:)) , vect_aux1(:)  )    +  X0Cell(:)
                !PRINT*, Position_aux1(:) 
                vect_aux2(1) = ia
                vect_aux2(2) = ib
                vect_aux2(3) = ic
                Position_aux2(:) = matmul( TRANSPOSE (cell(:,:)) , vect_aux2(:)  )    +  X0Cell(:)
                !PRINT*,Position_aux2(:)
                IF ( ( Position_aux1(1) - atoms_position(iatom,1 ) <= ZERO  ) .AND. ( Position_aux2(1) - atoms_position(iatom,1 ) > ZERO  ))  THEN
                     IF ( ( Position_aux1(2) - atoms_position(iatom,2 ) <= ZERO  ) .AND. ( Position_aux2(2) - atoms_position(iatom,2 ) > ZERO  ))  THEN
                        IF ( ( Position_aux1(3) - atoms_position(iatom,3 ) <= ZERO  ) .AND. ( Position_aux2(3) - atoms_position(iatom,3 ) > ZERO ) )  THEN
                    found=.TRUE. 
                    Atomsongrid(iatom,:)=vect_aux2(:)
                    exit icloop
                ENDIF
                ENDIF
                ENDIF
            ENDDO ialoop
        ENDDO ibloop
       ENDDO icloop
       Format = "(3I6, 2F12.7)"
       PRINT*,"Atom:",iatom, "Type:",  atoms_type(iatom),"Actual position:", atoms_position(iatom,:)
       PRINT*," ... Found: ", found
       PRINT*," ... Approximated Position:", Atomsongrid(iatom,:)
       PRINT*, "... Delta:", target_density(Atomsongrid(iatom,1),Atomsongrid(iatom,2),Atomsongrid(iatom,3)) 
       PRINT*, "... Finite Bias:", density_Finitebias(Atomsongrid(iatom,1),Atomsongrid(iatom,2),Atomsongrid(iatom,3)) 
       PRINT*, "...     0  Bias:", density_0bias(Atomsongrid(iatom,1),Atomsongrid(iatom,2),Atomsongrid(iatom,3)) 
   
   ENDDO
   ! Atomsongrid(:,:)



   DEALLOCATE (Atomsongrid)
!  !
   ! Print the Difference in a Gaussian file
   !

   CALL GaussianFilewrite(FileOutputGaussian, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, target_density)

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
		     Integrated_density1(ib,ic) =  Integrated_density1(ib,ic) +  target_density(ia,ib,ic) / mesh(1)
		   ENDDO
	   ENDDO
   ENDDO

  DO ic=1,mesh(3)
       	  DO  ia=1,mesh(1)
	           DO ib=1,mesh(2)
		     Integrated_density2(ia,ic) =  Integrated_density2(ia,ic) +  target_density(ia,ib,ic) / mesh(2)
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
		     Integrated_density_xy(ic) =  Integrated_density_xy(ic) +  target_density(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_density_xy(ic) = Integrated_density_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO

!  !
   ! Print the Difference in a .dat file
   !
   CALL datwrite(iufile, FileOutputDat, mesh(3), cell(3,3), Integrated_density_xy(:))

   Integrated_density_xy(:) = ZERO 
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_density_xy(ic) =  Integrated_density_xy(ic) +   density_0bias(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_density_xy(ic) = Integrated_density_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO


!  !
   ! Print the 0bias file averaged in a .dat file
   !
   CALL datwrite(iufile, FileOutputDat0bias, mesh(3), cell(3,3), Integrated_density_xy(:))

   Integrated_density_xy(:) = ZERO 
   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1)
		     Integrated_density_xy(ic) =  Integrated_density_xy(ic) +   density_Finitebias(ia,ib,ic)
		   ENDDO
	   ENDDO
           Integrated_density_xy(ic) = Integrated_density_xy(ic)/ ( mesh(1)* mesh(2) )
   ENDDO

!  !
   ! Print the Finite Bias averaged in a .dat file
   !
   CALL datwrite(iufile, FileOutputDatFinitebias, mesh(3), cell(3,3), Integrated_density_xy(:))


!  !
   ! Deallocate
   !
   PRINT*, "Deallocate"
   DEALLOCATE(density_0bias, density_Finitebias)
   DEALLOCATE(atoms_type, atoms_position)
   DEALLOCATE(target_density, Integrated_density1,Integrated_density2, Integrated_density_xy)
   DEALLOCATE ( vect1,vect2 )



   PRINT*, "END"
END PROGRAM PotentialDifference

