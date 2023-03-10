!
!      Copyright (C) 2005 WanT Group
!
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!***************************************************************************
   SUBROUTINE recursion_read_matrix( file, name, a, dim1, dim2, coord, dim3, dim4, dim5, dim6, nb_rep, size_C, size_ext)
   !***************************************************************************


  !
   ! First, the routine reads from stdin a namelist called MATRIX_DATA
   ! inside a NAME xml-iotk block, giving all the information related to the 
   ! required matrix A. Then the matrix is finally read from the wannier datafile.
   !
   USE KINDS
   USE constants,            ONLY : ZERO, CZERO, CONE
   USE parameters,       ONLY : nstrx
   USE files_module,     ONLY : file_open, file_close, file_delete
   USE io_module,        ONLY : aux_unit
   USE io_global_module,     ONLY : stdin
   USE iotk_module
   USE timing_module
   USE parser_module

   IMPLICIT NONE

   ! 
   ! input variables
   !
   CHARACTER(*), INTENT(in)  :: file, name
   INTEGER,      INTENT(in)  :: dim1, dim2, dim3, dim4, dim5, dim6, nb_rep
   REAL(dbl),    INTENT(in)  :: size_C, size_ext

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

   INTEGER                   :: nrtot, ir, iwan, jwan, iwan_shift, jwan_shift, irep
   INTEGER                   :: index, ivr_aux(3), nr_aux(3), index0, index1

      !! for the moment, only the x direction of transport is implemented
   INTEGER                   :: transport_dir=1

   INTEGER,      ALLOCATABLE :: ivr(:,:)
   COMPLEX(dbl), ALLOCATABLE :: lmatrix(:,:,:)
   REAL(dbl), ALLOCATABLE:: lcoord(:,:)
   CHARACTER(nstrx)          :: attr
   LOGICAL                   :: found, check
   INTEGER                   :: ncols, nrows
   INTEGER,      ALLOCATABLE :: icols(:), irows(:)
   CHARACTER(nstrx)          :: cols, rows

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
   ! initial values
   !

    a(:,:)=CZERO
    coord(:,:)=ZERO


   !
   ! some checks
   !

   IF (dim1 /= dim2) CALL errore(subname, 'bad dimensionality for dim1 or dim2', 1 )
   IF (dim1 /= dim4) CALL errore(subname, 'bad dimensionality for dim4', 2 )
   IF (dim3 /= 3) CALL errore(subname, 'bad dimensionality for dim3 : /= 3', 3 )
   IF (dim5 > dim1 ) CALL errore(subname, 'bad dimensionality for dim5 > dim1', 3 )
   IF (dim6 > dim1 ) CALL errore(subname, 'bad dimensionality for dim6 > dim1', 3 )


            SELECT CASE( TRIM(name) )
            CASE( "H_C" )

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
               !
               IF (ldimwann <=0 ) CALL errore(subname, 'invalid dimwann', 4)
               IF (nrtot <=0 ) CALL errore(subname, 'invalid nrtot', 5)
                !
                !
               IF ( dim5 /= ldimwann) CALL errore(subname, 'dimwan in input /= dimwan', 6 )
                !
               !
               ALLOCATE( lmatrix(ldimwann,ldimwann,1), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating lmatrix', ABS(ierr) )

               ALLOCATE( lcoord(dim3,ldimwann), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating lcoord', ABS(ierr) )

               ALLOCATE( ivr(3,nrtot), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating ivr', ABS(ierr) )
                 !
 
               ivr(:,:)=0
               lcoord(:,:)=ZERO
               lmatrix(:,:,:)=CZERO
              CALL iotk_scan_dat(aux_unit, "IVR", ivr, IERR=ierr)
                  IF (ierr/=0) CALL errore(subname, 'searching indxws', ABS(ierr) )
                 !
               !
               ! get the desired R indexes
               !
               CALL iotk_scan_begin(aux_unit, "RHAM", IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching RHAM', ABS(ierr) )
               !
                !
                  !
                  ! set the 3D corresponding R vector
                  !
                 !
                     ivr_aux(:) = 0
                !
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
                !
                  !
                  ! read the 3D R matrix corresponding to index
                  !
                  CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index)), lmatrix(:,:,1), IERR=ierr)
                  IF (ierr/=0) &
                     CALL errore(subname, 'searching VR'//TRIM(iotk_index(index)), ABS(ierr) )
                  !
                  !
                  DO jwan=1,dim5
                     DO iwan=1,dim5
                        a(iwan, jwan) = lmatrix( iwan, jwan, 1 )
                     ENDDO
                  ENDDO
                 !

                    !
                    !
               CALL iotk_scan_end(aux_unit, "RHAM", IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching end of RHAM', ABS(ierr) )
                  !
               CALL iotk_scan_dat(aux_unit, "WANCENTER", lcoord, IERR=ierr)
                  IF (ierr/=0) CALL errore(subname, 'searching coord', ABS(ierr) )

               DO iwan=1,dim5
                  coord(:, iwan) = lcoord( :, iwan )
               ENDDO

                 !
                 !
                 !
               CALL file_close( aux_unit, PATH="/HAMILTONIAN/", ACTION="read" )
                 !
                 !
                 !
                 !
            !
            ! cleaning local workspace
            !
                  !


            CASE( "H_R" )


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
               !
               IF (ldimwann <=0 ) CALL errore(subname, 'invalid dimwann', 7 )
               IF (nrtot <=0 ) CALL errore(subname, 'invalid nrtot', 8 )
                !
                !
               IF ( dim6 /= ldimwann) CALL errore(subname, 'dimwan in input /= dim_L', 9 )
                !
               !
               ALLOCATE( lmatrix(ldimwann,ldimwann,2), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating lmatrix', ABS(ierr) )

               ALLOCATE( lcoord(dim3,ldimwann), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating lcoord', ABS(ierr) )

               ALLOCATE( ivr(3,nrtot), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating ivr', ABS(ierr) )
                 !
 
               ivr(:,:)=0
               lcoord(:,:)=ZERO
               lmatrix(:,:,:)=CZERO
              CALL iotk_scan_dat(aux_unit, "IVR", ivr, IERR=ierr)
                  IF (ierr/=0) CALL errore(subname, 'searching indxws', ABS(ierr) )
                 !
               !
               ! get the desired R indexes
               !
               CALL iotk_scan_begin(aux_unit, "RHAM", IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching RHAM', ABS(ierr) )
               !

                  ivr_aux(:) = 0
                  !
                  ! search the 3D index corresponding to ivr_aux
                  !
                  found = .FALSE.
                  DO ir = 1, nrtot
                     ! 
                     IF ( ALL( ivr(:,ir) == ivr_aux(:) ) )  THEN
                           found = .TRUE.
                           index0 = ir 
                           EXIT 
                     ENDIF
                  ENDDO
                  IF ( .NOT. found ) CALL errore(subname, '3D R-vector not found', 5 )


                  ivr_aux(transport_dir) = 1
                  !
                  ! search the 3D index corresponding to ivr_aux
                  !
                  found = .FALSE.
                  DO ir = 1, nrtot
                     ! 
                     IF ( ALL( ivr(:,ir) == ivr_aux(:) ) )  THEN
                           found = .TRUE.
                           index1 = ir 
                           EXIT 
                     ENDIF
                  ENDDO
                  IF ( .NOT. found ) CALL errore(subname, '3D R-vector not found', 6 )
                  ! read the 3D R matrix corresponding to index
                  !
                  CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index0)), lmatrix(:,:,1), IERR=ierr)
                  IF (ierr/=0) &
                     CALL errore(subname, 'searching VR'//TRIM(iotk_index(index0)), ABS(ierr) )
                  CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index1)), lmatrix(:,:,2), IERR=ierr)
                  IF (ierr/=0) &
                     CALL errore(subname, 'searching VR'//TRIM(iotk_index(index1)), ABS(ierr) )
                              !


                  ! Initialize hamiltonian ON site
                  DO irep=1, nb_rep+1
                     DO jwan=1,ldimwann
                        jwan_shift= dim5 + jwan + (irep-1)*ldimwann
                        DO iwan=1,ldimwann
                           iwan_shift=dim5 + iwan + (irep-1)*ldimwann
                           a(iwan_shift, jwan_shift) = lmatrix( iwan, jwan, 1 )
                        ENDDO
                     ENDDO
                  ENDDO

                 ! hopping (only if nb rep =/0
                  IF (nb_rep > 0) THEN
                     DO irep=1, nb_rep
                        DO jwan=1,ldimwann
                           jwan_shift= dim5 + jwan + (irep)*ldimwann
                           DO iwan=1,ldimwann
                              iwan_shift=dim5 + iwan + (irep-1)*ldimwann
                              a(iwan_shift, jwan_shift) = lmatrix( iwan, jwan, 2 )
                              a(jwan_shift, iwan_shift) = CONJG (lmatrix( iwan, jwan, 2 ))
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF


               CALL iotk_scan_end(aux_unit, "RHAM", IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching end of RHAM', ABS(ierr) )
                  !
               CALL iotk_scan_dat(aux_unit, "WANCENTER", lcoord, IERR=ierr)
                  IF (ierr/=0) CALL errore(subname, 'searching coord', ABS(ierr) )

               DO irep=1, nb_rep+1
                  DO iwan=1,ldimwann
                     iwan_shift= dim5 + iwan + (irep-1)*ldimwann
                      ! seulement valide si trnaport selon l'axe x, pas dur ? g?n?ralizer mais flemme
                     coord(1, iwan_shift) = lcoord( 1, iwan ) + size_C + (irep-1)*size_ext
                     coord(2, iwan_shift) = lcoord( 2, iwan )
                     coord(3, iwan_shift) = lcoord( 3, iwan )
                  ENDDO
               ENDDO

                 !
                 !
                 !
               CALL file_close( aux_unit, PATH="/HAMILTONIAN/", ACTION="read" )
                 !



            CASE( "H_L" )

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
               !
               IF (ldimwann <=0 ) CALL errore(subname, 'invalid dimwann', 10)
               IF (nrtot <=0 ) CALL errore(subname, 'invalid nrtot', 11 )
                !
                !
               IF ( dim6 /= ldimwann) CALL errore(subname, 'dimwan in input /= dim_L', 12 )
                !
               !
               ALLOCATE( lmatrix(ldimwann,ldimwann,2), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating lmatrix', ABS(ierr) )

               ALLOCATE( lcoord(dim3,ldimwann), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating lcoord', ABS(ierr) )

               ALLOCATE( ivr(3,nrtot), STAT=ierr )
                  IF (ierr/=0) CALL errore(subname, 'allocating ivr', ABS(ierr) )
                 !
 
               ivr(:,:)=0
               lcoord(:,:)=ZERO
               lmatrix(:,:,:)=CZERO
              CALL iotk_scan_dat(aux_unit, "IVR", ivr, IERR=ierr)
                  IF (ierr/=0) CALL errore(subname, 'searching indxws', ABS(ierr) )
                 !
               !
               ! get the desired R indexes
               !
               CALL iotk_scan_begin(aux_unit, "RHAM", IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching RHAM', ABS(ierr) )
               !

                  ivr_aux(:) = 0
                  !
                  ! search the 3D index corresponding to ivr_aux
                  !
                  found = .FALSE.
                  DO ir = 1, nrtot
                     ! 
                     IF ( ALL( ivr(:,ir) == ivr_aux(:) ) )  THEN
                           found = .TRUE.
                           index0 = ir 
                           EXIT 
                     ENDIF
                  ENDDO
                  IF ( .NOT. found ) CALL errore(subname, '3D R-vector not found', 7 )


                  ivr_aux(transport_dir) = -1

                !

                  !
                  ! search the 3D index corresponding to ivr_aux
                  !
                  found = .FALSE.
                  DO ir = 1, nrtot
                     ! 
                     IF ( ALL( ivr(:,ir) == ivr_aux(:) ) )  THEN
                           found = .TRUE.
                           index1 = ir 
                           EXIT 
                     ENDIF
                  ENDDO
                  IF ( .NOT. found ) CALL errore(subname, '3D R-vector not found', 8 )


                  !
                  ! read the 3D R matrix corresponding to index
                  !
                  CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index1)), lmatrix(:,:,2), IERR=ierr)
                  IF (ierr/=0) &
                     CALL errore(subname, 'searching VR'//TRIM(iotk_index(index1)), ABS(ierr) )
                              !
                  ! read the 3D R matrix corresponding to index
                  !
                  CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index0)), lmatrix(:,:,1), IERR=ierr)
                  IF (ierr/=0) &
                     CALL errore(subname, 'searching VR'//TRIM(iotk_index(index0)), ABS(ierr) )
                  !
                  ! ON site
                  DO irep=1, nb_rep+1
                  DO jwan=1,ldimwann
                     jwan_shift= dim5 + jwan + (irep-1)*ldimwann
                     DO iwan=1,ldimwann
                        iwan_shift=dim5 + iwan + (irep-1)*ldimwann
                        a(iwan_shift, jwan_shift) = lmatrix( iwan, jwan, 1 )
                     ENDDO
                  ENDDO
                  ENDDO

                 ! hopping (only if nb rep =/0
                  IF (nb_rep > 0) THEN
                     DO irep=1, nb_rep
                        DO jwan=1,ldimwann
                           jwan_shift= dim5 + jwan + (irep)*ldimwann
                           DO iwan=1,ldimwann
                              iwan_shift=dim5 + iwan + (irep-1)*ldimwann
                              a(iwan_shift, jwan_shift) = lmatrix( iwan, jwan, 2 )
                              a(jwan_shift, iwan_shift) = CONJG (lmatrix( iwan, jwan, 2 ))
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF


               CALL iotk_scan_end(aux_unit, "RHAM", IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching end of RHAM', ABS(ierr) )
                  !
               CALL iotk_scan_dat(aux_unit, "WANCENTER", lcoord, IERR=ierr)
                  IF (ierr/=0) CALL errore(subname, 'searching coord', ABS(ierr) )

               DO irep=1, nb_rep+1
                  DO iwan=1,ldimwann
                     iwan_shift= dim5 + iwan + (irep-1)*ldimwann
                      ! seulement valide si trnaport selon l'axe x, pas dur ? g?n?ralizer mais flemme
                     coord(1, iwan_shift) = lcoord( 1, iwan ) - size_C - (irep - 1)*  size_ext
                     coord(2, iwan_shift) = lcoord( 2, iwan )
                     coord(3, iwan_shift) = lcoord( 3, iwan )
                  ENDDO
               ENDDO

                 !
                 !
                 !
               CALL file_close( aux_unit, PATH="/HAMILTONIAN/", ACTION="read" )
                 !



            CASE( "H_CR" )
                  !
                  !  initialize which WF are coupled in C and L
                  !  must be set in the input file until the form 
                  !  <H_CR rows="" cols="" \>
                     CALL iotk_scan_empty(stdin, TRIM(name), ATTR=attr, IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching for '//TRIM(name), ABS(ierr) )
                     !
                     CALL iotk_scan_attr(attr, 'cols', cols, FOUND=found, IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching for cols', ABS(ierr) )
                     IF( .NOT. found ) cols = 'all'
                     CALL change_case( cols, 'lower')
                     !
                     CALL iotk_scan_attr(attr, 'rows', rows, FOUND=found, IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching for rows', ABS(ierr) )
                     IF( .NOT. found ) rows = 'all'
                     CALL change_case( rows, 'lower')
                  !
                  !
                  !
                  ! parse the obtained data
                  !
                     !
                     ! deal with rows or cols = "all"
                     !
                     IF ( TRIM(rows) == "all" ) rows="1-"//TRIM( int2char(dim5))
                     IF ( TRIM(cols) == "all" ) cols="1-"//TRIM( int2char(dim6))
                     !
                     !
                     ! get the number of required rows and cols
                     CALL parser_replica( rows, nrows, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in rows string I',ABS(ierr))
                     CALL parser_replica( cols, ncols, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in cols string I',ABS(ierr))
                     ! 
                  !
                  !
                     ! Quelques doutes sur cette condition piquee dans WanT
                     ! Peut etre ? supprimer car nrow peut etre inferieur a dim 5
                     ! mais j ai pas tout compris de toute facon
                     !IF ( nrows /= dim5 ) CALL errore(subname,'invalid number of rows',3)
                     !IF ( ncols /= dim6 ) CALL errore(subname,'invalid number of cols',3)
                     !
                     ALLOCATE( irows(nrows), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating irows', ABS(ierr) )
                     ALLOCATE( icols(ncols), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating icols', ABS(ierr) )
                     !
                     ! get the actual indexes for rows and cols
                     CALL parser_replica( rows, nrows, irows, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in rows string II',ABS(ierr))
                     CALL parser_replica( cols, ncols, icols, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in cols string II',ABS(ierr))
                     !
                     ! simple check
                     DO i=1,nrows
                        IF ( irows(i) <=0 ) CALL errore(subname,'invalid irows(i)',i) 
                     ENDDO
                     DO i=1,ncols
                        IF ( icols(i) <=0 ) CALL errore(subname,'invalid icols(i)',i) 
                     ENDDO
                     !
                     ! find hamiltonian data in the .ham file
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
                  !
                     IF (ldimwann <=0 ) CALL errore(subname, 'invalid dimwann', 15)
                     IF (nrtot <=0 ) CALL errore(subname, 'invalid nrtot', 16 )
                     IF ( dim5 /= ldimwann) CALL errore(subname, 'dimwan in input /= dim_C', 17 )
                  !
                     !
                     !
                     DO i=1,ncols
                        IF ( icols(i) > dim6 ) CALL errore(subname, 'invalid icols(i)', i)
                     ENDDO
                     DO i=1,nrows
                        IF ( irows(i) > dim5 ) CALL errore(subname, 'invalid irows(i)', i)
                     ENDDO
                  !
                     !
                     ALLOCATE( ivr(3,nrtot), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating ivr', ABS(ierr) )
                     ALLOCATE( lmatrix(ldimwann,ldimwann,1), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating lmatrix', ABS(ierr) )
                  !
                     ivr(:,:)=0
                     lmatrix(:,:,:)=CZERO
                     !
                     CALL iotk_scan_dat(aux_unit, "IVR", ivr, IERR=ierr)
                        IF (ierr/=0) CALL errore(subname, 'searching indxws', ABS(ierr) )
                        !
                        !
                  !
                  !
                     !
                     ! get the desired R indexes
                     !
                     CALL iotk_scan_begin(aux_unit, "RHAM", IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching RHAM', ABS(ierr) )
                              !  set the 3D corresponding R vector
                              !
                       ivr_aux(:) = 0
                       ivr_aux(transport_dir) = 1
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
                        IF ( .NOT. found ) CALL errore(subname, '3D R-vector not found', 9)
                  !
                  !
                        !
                        ! read the 3D R matrix corresponding to index
                        !
                        CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index)), lmatrix(:,:,1), IERR=ierr)
                        IF (ierr/=0) &
                           CALL errore(subname, 'searching VR'//TRIM(iotk_index(index)), ABS(ierr) )
                  !
                        !
                        ! cut the total hamiltonian according to the required rows and cols
                        !
                        DO jwan=1,ncols
                           jwan_shift = icols(jwan) + ldimwann
                           DO iwan=1,nrows
                              !
                              a( irows(iwan) , jwan_shift) = lmatrix( irows(iwan), icols(jwan) , 1)
                              a(jwan_shift,  irows(iwan) ) = CONJG ( lmatrix( irows(iwan), icols(jwan) , 1 ) )
                           ENDDO
                        ENDDO
                  !
                  !
                     CALL iotk_scan_end(aux_unit, "RHAM", IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching end of RHAM', ABS(ierr) )
                  !
                     CALL file_close( aux_unit, PATH="/HAMILTONIAN/", ACTION="read" )
                  !
                  !
                  !
                  ! cleaning local workspace
                  !
                     DEALLOCATE( icols, irows, STAT=ierr)
                        IF (ierr/=0) CALL errore(subname, 'deallocating icols, irows', ABS(ierr) )
                  !


            CASE( "H_CL" )
                 !
                  !  initialize which WF are coupled in C and L
                  !  must be set in the input file until the form 
                  !  <H_CR rows="" cols="" \>
 
                     CALL iotk_scan_empty(stdin, TRIM(name), ATTR=attr, IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching for '//TRIM(name), ABS(ierr) )
                     !
                     CALL iotk_scan_attr(attr, 'cols', cols, FOUND=found, IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching for cols', ABS(ierr) )
                     IF( .NOT. found ) cols = 'all'
                     CALL change_case( cols, 'lower')
                     !
                     CALL iotk_scan_attr(attr, 'rows', rows, FOUND=found, IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching for rows', ABS(ierr) )
                     IF( .NOT. found ) rows = 'all'
                     CALL change_case( rows, 'lower')
                  !
                  !
                  !
                  ! parse the obtained data
                  !
                     !
                     ! deal with rows or cols = "all"
                     !
                     IF ( TRIM(rows) == "all" ) rows="1-"//TRIM( int2char(dim5))
                     IF ( TRIM(cols) == "all" ) cols="1-"//TRIM( int2char(dim6))
                     !
                     !
                     ! get the number of required rows and cols
                     CALL parser_replica( rows, nrows, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in rows string I',ABS(ierr))
                     CALL parser_replica( cols, ncols, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in cols string I',ABS(ierr))
                     ! 
                  !
                  !
                     ! Quelques doutes sur cette condition piquee dans WanT
                     ! Peut etre ? supprimer car nrow peut etre inferieur a dim 5
                     ! mais j ai pas tout compris de toute facon
                     !IF ( nrows /= dim5 ) CALL errore(subname,'invalid number of rows',3)
                     !IF ( ncols /= dim6 ) CALL errore(subname,'invalid number of cols',3)
                     !
                     ALLOCATE( irows(nrows), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating irows', ABS(ierr) )
                     ALLOCATE( icols(ncols), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating icols', ABS(ierr) )
                     !
                     ! get the actual indexes for rows and cols
                     CALL parser_replica( rows, nrows, irows, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in rows string II',ABS(ierr))
                     CALL parser_replica( cols, ncols, icols, IERR=ierr)
                     IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in cols string II',ABS(ierr))
                     !
                     ! simple check
                     DO i=1,nrows
                        IF ( irows(i) <=0 ) CALL errore(subname,'invalid irows(i)',i) 
                     ENDDO
                     DO i=1,ncols
                        IF ( icols(i) <=0 ) CALL errore(subname,'invalid icols(i)',i) 
                     ENDDO
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
                  !
                     IF (ldimwann <=0 ) CALL errore(subname, 'invalid dimwann', 18)
                     IF (nrtot <=0 ) CALL errore(subname, 'invalid nrtot', 19)
                     IF ( dim5 /= ldimwann) CALL errore(subname, 'dimwan in input /= dim_C', 20 )
                  !
                     !
                     !
                     DO i=1,ncols
                        IF ( icols(i) > dim5 ) CALL errore(subname, 'invalid icols(i)', i)
                     ENDDO
                     DO i=1,nrows
                        IF ( irows(i) > dim6 ) CALL errore(subname, 'invalid irows(i)', i)
                     ENDDO
                  !
                     !
                     ALLOCATE( ivr(3,nrtot), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating ivr', ABS(ierr) )
                     ALLOCATE( lmatrix(ldimwann,ldimwann,1), STAT=ierr )
                        IF (ierr/=0) CALL errore(subname, 'allocating lmatrix', ABS(ierr) )
                  !

                     ivr(:,:)=0
                     lmatrix(:,:,:)=CZERO
                  !
                     CALL iotk_scan_dat(aux_unit, "IVR", ivr, IERR=ierr)
                        IF (ierr/=0) CALL errore(subname, 'searching indxws', ABS(ierr) )
                  !
                  !
                  !
                     !
                     ! get the desired R indexes
                     !
                     CALL iotk_scan_begin(aux_unit, "RHAM", IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching RHAM', ABS(ierr) )
                              !  set the 3D corresponding R vector
                     !
                     ivr_aux(:) = 0
                     ivr_aux(transport_dir) = -1

                  !
                  !
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
                        IF ( .NOT. found ) CALL errore(subname, '3D R-vector not found', 9 )
                  !
                  !
                        !
                        ! read the 3D R matrix corresponding to index
                        !
                        CALL iotk_scan_dat( aux_unit, "VR"//TRIM(iotk_index(index)), lmatrix(:,:,1), IERR=ierr)
                        IF (ierr/=0) &
                           CALL errore(subname, 'searching VR'//TRIM(iotk_index(index)), ABS(ierr) )
                  !
                        !
                        ! cut the total hamiltonian according to the required rows and cols
                        !
                        DO jwan=1,ncols
                           !modif le 23042007
                           !jwan_shift = icols(jwan) + ldimwann
                           jwan_shift = jwan + ldimwann
                           ! euhhh ca semble bizarre quand meme
                           DO iwan=1,nrows
                              !
                              a( irows(iwan), jwan_shift) = lmatrix( irows(iwan), icols(jwan),1 )
                              a(jwan_shift,  irows(iwan) ) = CONJG ( lmatrix( irows(iwan), icols(jwan) , 1 ) )
                           ENDDO
                        ENDDO
                  !
                  !
                     CALL iotk_scan_end(aux_unit, "RHAM", IERR=ierr)
                     IF (ierr/=0) CALL errore(subname, 'searching end of RHAM', ABS(ierr) )
                  !
                     CALL file_close( aux_unit, PATH="/HAMILTONIAN/", ACTION="read" )
                  !
                  !
                  !
                  ! cleaning local workspace
                  !
                     DEALLOCATE( icols, irows, STAT=ierr)
                        IF (ierr/=0) CALL errore(subname, 'deallocating icols, irows', ABS(ierr) )

            CASE DEFAULT
                CALL errore(subname, 'invalid name = '//TRIM(name), ABS(ierr) )
            END SELECT




   IF (ALLOCATED(lcoord)) THEN
      DEALLOCATE( lcoord, STAT=ierr )
      IF (ierr/=0) CALL errore(subname, 'deallocating lcoord', ABS(ierr) )
   ENDIF
   DEALLOCATE( lmatrix, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating lmatrix', ABS(ierr) )

   DEALLOCATE( ivr, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating ivr', ABS(ierr) )




   CALL timing( 'recursion_read_matrix', OPR='stop' )

END SUBROUTINE recursion_read_matrix

