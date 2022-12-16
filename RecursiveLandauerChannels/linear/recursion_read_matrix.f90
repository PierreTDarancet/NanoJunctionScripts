!
!      Copyright (C) 2005 WanT Group
!
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!***************************************************************************
   SUBROUTINE recursion_read_matrix( file, a, dim1, dim2, coord, dim3, dim4)
   !***************************************************************************
  !
   ! First, the routine reads from stdin a namelist called MATRIX_DATA
   ! inside a NAME xml-iotk block, giving all the information related to the 
   ! required matrix A. Then the matrix is finally read from the wannier datafile.
   !
   USE KINDS
   USE parameters,       ONLY : nstrx
   USE files_module,     ONLY : file_open, file_close, file_delete
   USE io_module,        ONLY : aux_unit
   USE iotk_module
   USE timing_module
   USE parser_module

   IMPLICIT NONE

   ! 
   ! input variables
   !
   CHARACTER(*), INTENT(in)  :: file
   INTEGER,      INTENT(in)  :: dim1, dim2, dim3, dim4
   !
   !  output variables
   !
   COMPLEX(dbl), INTENT(out) :: a(dim1,dim2)
   REAL(dbl), INTENT(out) :: coord(dim3, dim4)
 
  !
   !  input/output variables
   !
   !
   ! local variables
   !
   INTEGER   :: i, j, ierr
   CHARACTER(23)             :: subname="recursion_read_matrix"
   INTEGER                   :: ldimwann

   INTEGER                   :: nrtot, ir
   INTEGER                   :: index, ivr_aux(3), nr_aux(3)

      !! for the moment, only the x direction of transport is implemented
   INTEGER                   :: transport_dir=1

   INTEGER,      ALLOCATABLE :: ivr(:,:)

   CHARACTER(nstrx)          :: attr
   LOGICAL                   :: found, check

   !
   ! end of declarations
   !
!
!----------------------------------------
! main Body
!----------------------------------------
!
   CALL timing( 'recursion_read_matrix', OPR='start' )

   !
   ! some checks
   !


   ! some checks

   IF (dim1 /= dim2) CALL errore(subname, 'bad dimensionality for dim1 or dim2', 1 )
   IF (dim1 /= dim4) CALL errore(subname, 'bad dimensionality for dim4', 2 )
   IF (dim3 /= 3) CALL errore(subname, 'bad dimensionality for dim3 : /= 3', 3 )

! reading form iotk-formatted .ham file produced by wannier
!

   !
   ! read basic quantities
   !



!
   CALL file_open( aux_unit, TRIM(file), PATH="/HAMILTONIAN/", &
                   ACTION="read", FORM="formatted" )
   !
   CALL iotk_scan_empty(aux_unit, "DATA", ATTR=attr, IERR=ierr)
      IF (ierr/=0) CALL errore(subname, 'searching DATA', ABS(ierr) )
   CALL iotk_scan_attr(attr,"dimwann",ldimwann, IERR=ierr)
      IF (ierr/=0) CALL errore(subname, 'searching dimwann', ABS(ierr) )
   CALL iotk_scan_attr(attr,"nr",nr_aux, IERR=ierr)
      IF (ierr/=0) CALL errore(subname, 'searching nr', ABS(ierr) )
   CALL iotk_scan_attr(attr,"nrtot",nrtot, IERR=ierr)
      IF (ierr/=0) CALL errore(subname, 'searching nrtot', ABS(ierr) )

   IF (ldimwann <=0 ) CALL errore(subname, 'invalid dimwann', 4)
   IF (nrtot <=0 ) CALL errore(subname, 'invalid nrtot', 5)
   IF ( dim1 /= ldimwann) CALL errore(subname, 'dimwan in input /= dimwan', 6 )


   !
   ALLOCATE( ivr(3,nrtot), STAT=ierr )
      IF (ierr/=0) CALL errore(subname, 'allocating ivr', ABS(ierr) )

   CALL iotk_scan_dat(aux_unit, "IVR", ivr, IERR=ierr)
      IF (ierr/=0) CALL errore(subname, 'searching indxws', ABS(ierr) )

   !
   ! get the desired R indexes
   !
   CALL iotk_scan_begin(aux_unit, "RHAM", IERR=ierr)
   IF (ierr/=0) CALL errore(subname, 'searching RHAM', ABS(ierr) )


      !
      ! set the 3D corresponding R vector
      !

      DO i=1,3
         ivr_aux(i) = 0
      ENDDO

      !
      ! search the 3D index corresponding to ivr_aux
      !
      found = .FALSE.
      DO ir = 1, nrtot
          ! 
          IF ( ALL( ivr(:,ir) == ivr_aux(:) ) )  THEN
               found = .TRUE.
               index = ir 
               EXIT 
          ENDIF
      ENDDO
      !
      IF ( .NOT. found ) CALL errore(subname, '3D R-vector not found', 4 )


      !
      ! read the 3D R matrix corresponding to index
      !
      CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index)), a, IERR=ierr)
      IF (ierr/=0) &
         CALL errore(subname, 'searching VR'//TRIM(iotk_index(index)), ABS(ierr) )


   CALL iotk_scan_end(aux_unit, "RHAM", IERR=ierr)
   IF (ierr/=0) CALL errore(subname, 'searching end of RHAM', ABS(ierr) )

   CALL iotk_scan_dat(aux_unit, "WANCENTER", coord, IERR=ierr)
      IF (ierr/=0) CALL errore(subname, 'searching coord', ABS(ierr) )



   CALL file_close( aux_unit, PATH="/HAMILTONIAN/", ACTION="read" )




!
! cleaning local workspace
!
   DEALLOCATE( ivr, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating ivr', ABS(ierr) )

   CALL timing( 'recursion_read_matrix', OPR='stop' )

END SUBROUTINE recursion_read_matrix
! 
