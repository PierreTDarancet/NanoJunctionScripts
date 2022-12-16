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
   PROGRAM OverlapCalculations
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   USE IO_module,                              ONLY :  GaussianFilereadHeader, GaussianFilewrite, GaussianFileread, ReadInputFileOverlap

!!!!!!!!!!!!





   IMPLICIT NONE
   INTEGER, PARAMETER :: dimensions = 3    

   INTEGER, PARAMETER :: iufile=100
   REAL*8, PARAMETER ::  Cutoffdistance=20004.5

   REAL*8, PARAMETER ::  Ang_to_Bohr=1.889725989*ONE
   CHARACTER*100, PARAMETER :: datafile="DATA.in"


   !
   ! Input variables
   !
   CHARACTER*100 ::  FileInput
   CHARACTER*100 ::  FileOutputGaussian

   REAL*8 ::   cell(dimensions,dimensions), norm
   REAL*8 ::   X0Cell(dimensions), maxdistance
   REAL*8 ::   gold_position(dimensions)
   INTEGER ::   mesh(dimensions)
   INTEGER ::   mesh_test(dimensions)
   INTEGER ::   Nb_atoms
   INTEGER ::   Nb_distances
   INTEGER ::   Reference_atom

   REAL*8, ALLOCATABLE ::   density(:,:,:)
   REAL*8, ALLOCATABLE ::   distances(:)
   REAL*8, ALLOCATABLE ::   IntegratedVolumes(:)
   REAL*8, ALLOCATABLE ::   overlap(:)
   REAL*8, ALLOCATABLE ::   overlap_positiveValues(:)
   REAL*8, ALLOCATABLE ::   overlap_negativeValues(:)
   INTEGER, ALLOCATABLE ::   NumberOfPoints(:)


   INTEGER, ALLOCATABLE ::   atoms_type(:)
   REAL*8, ALLOCATABLE ::   atoms_position(:,:)
   !
   ! Output variables
   !
   ! 
   ! Local variables
   ! 
   INTEGER :: ia, ib, ic, id, idimension
   REAL*8, ALLOCATABLE  :: vect1(:), vect2(:)
   REAL*8 ::distancetoat, R(dimensions), Voxcell, Crossaux(dimensions), aux(dimensions), norm_check, volume_check
!***********************************************
!
!------------------------------
! main body
!------------------------------
!  !
   ! read input General Information
   !
   CALL ReadInputFileOverlap(iufile, datafile, FileInput, gold_position, Nb_distances, maxdistance,Reference_atom, FileOutputGaussian)


!  !
   ! read input file 
   !


   CALL GaussianFilereadHeader(FileInput, iufile, Cutoffdistance, dimensions, Nb_atoms, mesh)
   PRINT*, "Allocating..."
   ALLOCATE( atoms_type( Nb_atoms) )   
   ALLOCATE(atoms_position( Nb_atoms, dimensions) )   
   ALLOCATE( density(mesh(1),mesh(2),mesh(3)) )   
   ALLOCATE( distances(Nb_distances) )
   ALLOCATE( IntegratedVolumes(Nb_distances) )
   ALLOCATE( overlap(Nb_distances) )
   ALLOCATE( overlap_positiveValues(Nb_distances), overlap_negativeValues(Nb_distances)  )
   ALLOCATE( NumberOfPoints(Nb_distances) )

   PRINT*, "Initializing Distances..."
   DO ia = 1, Nb_distances
          distances(ia) =  ia*maxdistance/Nb_distances
   ENDDO
   overlap(:)   = ZERO
   overlap_positiveValues(:)   = ZERO
   overlap_negativeValues(:)   = ZERO
   IntegratedVolumes(:) = ZERO
   NumberOfPoints(:) = 0

   PRINT*, "Initializing Density..."
   density(:,:,:) = ZERO

   PRINT*, "Reading Density file..."
   CALL GaussianFileread( FileInput, iufile,Cutoffdistance, dimensions, Nb_atoms,  mesh, X0Cell, cell, atoms_type,atoms_position, density, norm)
   PRINT*, "Normalizing...", norm
   density(:,:,:) =    density(:,:,:)/norm

   PRINT*, "Writing Normalized Density..."
   CALL GaussianFilewrite(FileOutputGaussian, iufile,Cutoffdistance, dimensions, Nb_atoms, mesh, X0Cell, cell,atoms_type,atoms_position, density)


   PRINT*, "Calculating Overlap..."

        PRINT'(A45,F12.6,A1,F12.6,A1,F12.6,A2)', "    ...Center of the sphere at [Bohr]: (",(gold_position(1)+atoms_position(Reference_Atom,1)),",",(gold_position(2)+atoms_position(Reference_Atom,2)),",",(gold_position(3)+atoms_position(Reference_Atom,3)), " )"
        
        Crossaux(1) = (cell(1,2)*cell(2,3))  -(cell(2,2)*cell(1,3) )
        Crossaux(2) = (cell(1,3)*cell(2,1))  -(cell(2,3)*cell(1,1))
        Crossaux(3) = (cell(1,1)*cell(2,2))  - (cell(2,1)*cell(1,2))
        aux(1) = cell(3,1)
        aux(2) = cell(3,2)
        aux(3) = cell(3,3)
        Voxcell = dot_product(Crossaux(:), aux(:))
   PRINT*, "    ...Voxcell:",  Voxcell


    DO ia = 1,  mesh(1)
    DO ib = 1,  mesh(2)
    DO ic = 1,  mesh(3)

            DO idimension=1,3
                 R(idimension) = (ia-1)*cell(1,idimension)+ (ib-1)*cell(2,idimension) + (ic-1)*cell(3,idimension) +X0Cell(idimension)   
            ENDDO
            distancetoat = sqrt( ( R(1) - gold_position(1)-atoms_position(Reference_Atom,1) )**2 + ( R(2) - gold_position(2)-atoms_position(Reference_Atom,2) )**2 + ( R(3) - gold_position(3)-atoms_position(Reference_Atom,3))**2 ) 

            DO  id=1 ,  Nb_distances
                 IF (distancetoat < distances(id) ) THEN
                       overlap(id)   =    overlap(id) +    (density(ia,ib,ic)*density(ia,ib,ic) *Voxcell)
                       IntegratedVolumes(id) =  IntegratedVolumes(id) +   Voxcell 
                       NumberOfPoints(id) = NumberOfPoints(id)+1
                 ENDIF

            ENDDO
            IF (density(ia,ib,ic)>ZERO) THEN
		    DO  id=1 ,  Nb_distances
		         IF (distancetoat < distances(id) ) THEN
		              overlap_positiveValues(id)   =    overlap_positiveValues(id) +    (density(ia,ib,ic)* Voxcell)
		        ENDIF

		    ENDDO
            ELSE 
		    DO  id=1 ,  Nb_distances
		         IF (distancetoat < distances(id) ) THEN
		              overlap_negativeValues(id)   =    overlap_negativeValues(id) +    (density(ia,ib,ic) *Voxcell)
		        ENDIF

		    ENDDO
            ENDIF          
            norm_check=norm_check+(density(ia,ib,ic)*density(ia,ib,ic) )*Voxcell
            volume_check=volume_check+Voxcell
    ENDDO     
    ENDDO
    ENDDO

                       overlap(:)  = sqrt(overlap(:) )
                       overlap_positiveValues(:)  = overlap_positiveValues(:)
                       overlap_negativeValues(:)  = overlap_negativeValues(:) 
   PRINT*, "Calculating Overlap... done"
   PRINT*, "Norm Check", sqrt(norm_check), "Volume Check",volume_check


   PRINT*, "##############################################################################"
   PRINT*, "###                                RESULTS                                 ###"
   PRINT*, "##############################################################################"
   PRINT'(A34,I5,A6,I5)', " ### Sphere centered around atom #", Reference_Atom," Type ", atoms_type(Reference_Atom)
   PRINT'(A18,F12.6,A3,F12.6,A3,F12.6,A3)', " ### Coordinates (", atoms_position(Reference_Atom,1)," , ", atoms_position(Reference_Atom,2), " , ", atoms_position(Reference_Atom,3), " )"
   PRINT'(A18,F12.6,A3,F12.6,A3,F12.6,A3)', " ### Distance    (", gold_position(1)," , ",gold_position(2), " , ", gold_position(3), " )"
   PRINT*, "##############################################################################"
   PRINT*, " Dist [Bohr] | Overlap | Nb of Points | Eq. Radius[Bohr] | overlap+ | overlap-"
   DO  id=1 ,  Nb_distances   
               PRINT'(I5,A3,F12.6,A3,F12.6,A3,I7,A3,F12.6,A3,F12.6,A3,F12.6)', id," | ", distances(id)," | ", overlap(id)," | ", NumberOfPoints(id)," | ", ( ( IntegratedVolumes(id)*3.0/(PI*4.0) )**(1.0/3.0) )," | ", overlap_positiveValues(id)," | ", overlap_negativeValues(id)
   ENDDO
   PRINT*, "##############################################################################"

!  !
   ! Deallocate
   !
   PRINT*, "Deallocate"
   DEALLOCATE(density)
   DEALLOCATE(atoms_type, atoms_position)
   DEALLOCATE(IntegratedVolumes, overlap, distances )
   DEALLOCATE(overlap_positiveValues, overlap_negativeValues )




   PRINT*, "END"
END PROGRAM OverlapCalculations

