!
! Copyright (C) 2005 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=====================================================
   PROGRAM quick_plot
   !=====================================================
   !
   ! real space plot of the computed Wannier functions
   !
   USE kinds
   USE constants,          ONLY : ZERO, CZERO, ONE, TWO, CONE, CI, TPI, &
                                  bohr => bohr_radius_angs, EPS_m6, EPS_m4
   USE parameters,         ONLY : nstrx
   USE io_module,          ONLY : prefix, postfix, work_dir, stdin, stdout
   USE io_module,          ONLY : ioname, space_unit, wan_unit, dft_unit, &
                                  aux_unit, aux1_unit, aux2_unit, aux3_unit
   USE parser_module

   USE iotk_module
   USE files_module,       ONLY : file_open, file_close, file_delete

   IMPLICIT NONE

   !
   ! input variables
   !
   CHARACTER(nstrx) :: wann             ! contains the list of WF indexes to plot 
                                        ! in the fmt e.g. "1-3,4,7-9"
   REAL(dbl)        :: r1min, r1max     ! plot cell dim along a1 (cry units)
   REAL(dbl)        :: r2min, r2max     ! the same but for a2
   REAL(dbl)        :: r3min, r3max     ! the same but for a3
   REAL(dbl)        :: avec_x_1
   REAL(dbl)        :: avec_x_2
   REAL(dbl)        :: avec_x_3
   REAL(dbl)        :: avec_y_1
   REAL(dbl)        :: avec_y_2
   REAL(dbl)        :: avec_y_3
   REAL(dbl)        :: avec_z_1
   REAL(dbl)        :: avec_z_2
   REAL(dbl)        :: avec_z_3
   REAL(dbl)        :: alat
   REAL(dbl)        :: sigma
   CHARACTER( 20 )  :: datatype         ! ( "modulus" | "real" | "imaginary"  )
   CHARACTER( 20 )  :: output_fmt       ! ( "txt" | "plt" | "cube" | "xsf" )
   LOGICAL          :: assume_ncpp      ! If .TRUE. pp's are not read
   LOGICAL          :: locate_wf        ! move the centers of WF in a unit cell centered
                                        ! around the origin
   INTEGER          :: nb_x
   INTEGER          :: nb_y
   INTEGER          :: nb_z
   !
   ! input namelist
   !
   NAMELIST /INPUT/ prefix, postfix, work_dir, &
                    datatype, assume_ncpp, output_fmt, locate_wf, &
                    r1min, r1max, r2min, r2max, r3min, r3max, avec_x_1, avec_x_2, avec_x_3, avec_y_1, avec_y_2, avec_y_3, avec_z_1, &
                    avec_z_2, avec_z_3, alat, nb_x, nb_y, nb_z

   !
   ! local variables
   !

   REAL(dbl)    :: avecl(3,3), r0(3), r1(3), rmin(3), rmax(3)
   REAL(dbl)        :: avec(3,3)           ! dir lattice vects (Bohr)
   REAL(dbl),    ALLOCATABLE :: rwann_out(:,:,:)
   REAL(dbl),    ALLOCATABLE :: tautot(:,:)
   INTEGER   :: nat
   REAL(dbl),    ALLOCATABLE :: WF_like(:,:,:)
   REAL(dbl),    ALLOCATABLE :: wan_coord(:,:)
   CHARACTER(3), ALLOCATABLE :: symbtot(:)
   REAL(dbl),    ALLOCATABLE :: ix_2_x(:), iy_2_y(:), iz_2_z(:)
   CHARACTER( nstrx )  :: filename
   CHARACTER( 5 )      :: str, aux_fmt
   CHARACTER( 10 )      :: subname='quick_plot'
   LOGICAL             :: lfound
   LOGICAL             :: okp( 3 )
   !
   ! added for recursion
   !
   INTEGER :: dimwan, dim_recursion, dim_subspace, n_iter
   INTEGER :: iter, i_sub, i_recur, ierr, i_x, i_y, i_z
   COMPLEX(dbl), ALLOCATABLE :: state(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: coeff(:,:)
   CHARACTER(nstrx)          :: attr

   !
   ! end of declariations
   !

!
!------------------------------
! main body
!------------------------------
!
      !CALL startup(version_number,'plot')

     PRINT*, 'Point 0'
!
! ... Read INPUT namelist from stdin
!
      prefix                      = 'WanT'
      postfix                     = ' '
      work_dir                    = './'
      assume_ncpp                 = .FALSE.
      locate_wf                   = .TRUE.
      wann                        = ' '
      datatype                    = 'modulus'
      output_fmt                  = 'plt'
      r1min                       = -0.5
      r1max                       =  0.5
      r2min                       = -0.5
      r2max                       =  0.5
      r3min                       = -0.5
      r3max                       =  0.5
      sigma                       = 1.0
      avec_x_1                    = 1.0
      avec_x_2                    = ZERO
      avec_x_3                    = ZERO
      avec_y_1                    = ZERO
      avec_y_2                    = 1.0
      avec_y_3                    = ZERO
      avec_z_1                    = ZERO
      avec_z_2                    = ZERO
      avec_z_3                    = 1.0
      alat                        = 1.0
      nb_x                        = 100
      nb_y                        = 100
      nb_z                        = 100

!!!!!!!!!!!!!




      nat=1
      ALLOCATE (tautot( 3, nat), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating tautot ', ABS(ierr) )

      tautot(:,:) = ZERO
      ALLOCATE (symbtot(nat), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating symbtot ', ABS(ierr) )
      symbtot(:)='Au'

      READ(stdin, INPUT, IOSTAT=ierr)
      IF ( ierr /= 0 )  CALL errore(subname,'Unable to read namelist INPUT',ABS(ierr))

     PRINT*, 'Point 1'
      !
      ! Some checks
      !
      !IF ( LEN_TRIM( wann) == 0 ) CALL errore('plot', 'wann not supplied ', 1)
      !
      CALL change_case( datatype, 'lower')
      IF ( TRIM(datatype) /= "modulus" .AND. TRIM(datatype) /= "real" .AND. &
           TRIM(datatype) /= "imaginary"  ) &
           CALL errore(subname,'invalid DATATYPE = '//TRIM(datatype),2)
           !
           !
      CALL change_case(output_fmt,'lower')
      IF ( TRIM(output_fmt) /= "txt" .AND. TRIM(output_fmt) /= "plt" .AND. &
           TRIM(output_fmt) /= "cube" .AND. TRIM(output_fmt) /= "xsf" ) &
           CALL errore(subname, 'Invalid output_fmt = '//TRIM(output_fmt), 4)

      IF ( r1min > r1max ) CALL errore(subname, 'r1min > r1max',1)
      IF ( r2min > r2max ) CALL errore(subname, 'r2min > r2max',2)
      IF ( r3min > r3max ) CALL errore(subname, 'r3min > r3max',3)
      IF ( ABS(r1max -r1min) < EPS_m4 ) CALL errore(subname, 'r1 too small',1)
      IF ( ABS(r2max -r2min) < EPS_m4 ) CALL errore(subname, 'r2 too small',2)
      IF ( ABS(r3max -r3min) < EPS_m4 ) CALL errore(subname, 'r3 too small',3)
      rmin(1) = r1min 
      rmin(2) = r2min
      rmin(3) = r3min
      rmax(1) = r1max 
      rmax(2) = r2max
      rmax(3) = r3max
      avec(1,1) =  avec_x_1 * alat
      avec(2,1) =  avec_x_2 * alat
      avec(3,1) =  avec_x_3 * alat
      avec(1,2) =  avec_y_1 * alat
      avec(2,2) =  avec_y_2 * alat
      avec(3,2) =  avec_y_3 * alat
      avec(1,3) =  avec_z_1 * alat
      avec(2,3) =  avec_z_2 * alat
      avec(3,3) =  avec_z_3 * alat

      r0(1) = rmin(1)*avec(1,1)
      r0(2) = rmin(2)*avec(2,2)
      r0(3) = rmin(3)*avec(3,3)

     PRINT*, 'Point 2'
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! INFO dans le fichier .cen  ie coordonnees centres fonctions de wannier !!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL ioname('wannier_center',filename)
            CALL file_open(aux_unit,TRIM(filename),PATH="/CENTERS/",ACTION="read", &
                              FORM='formatted')

            CALL iotk_scan_empty(aux_unit, "DATADIM", ATTR=attr, IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching DATADIM', ABS(ierr) )
            CALL iotk_scan_attr(attr,"dimwann_total",dimwan, IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching dimwann', ABS(ierr) )
            CALL iotk_scan_attr(attr,"dim_rec",dim_recursion, IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching dim_recursion', ABS(ierr) )
            CALL iotk_scan_attr(attr,"dim_sub",dim_subspace, IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching dim_subspace', ABS(ierr) )
            !
            CALL iotk_scan_empty(aux_unit, "DATAFINAL", ATTR=attr, IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching DATAFINAL', ABS(ierr) )
            CALL iotk_scan_attr(attr,"final_iter",n_iter, IERR=ierr)
               IF (ierr/=0) CALL errore(subname, 'searching final_iter', ABS(ierr) )
            !

            ALLOCATE (  wan_coord(3, dim_recursion), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating wan_coord ',ABS(ierr) )

            CALL iotk_scan_dat(aux_unit,"RECCENTER", wan_coord, IERR=ierr) 
                  IF (ierr/=0) CALL errore(subname,'searching RECCENTER',ABS(ierr))

            !
       CALL file_close(aux_unit,PATH="/CENTERS/",ACTION="read")
       CALL ioname('wannier_center',filename,LPATH=.FALSE.)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     PRINT*, 'Point 3'
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ALLOCATE (  coeff(n_iter+1, dim_recursion), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating coeff ', ABS(ierr) )

       ALLOCATE (  state(n_iter+1,dim_subspace, dim_recursion), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating state ', ABS(ierr) )

       state(:,:,:) = CZERO

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! INFO dans le fichier .sta  ie etats de recursion dans la base de rec   !!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL ioname('states',filename)
      CALL file_open(aux2_unit,TRIM(filename),PATH="/STATES/",ACTION="read", &
                              FORM='formatted')

            !
            CALL iotk_scan_begin(aux2_unit,"STA")
            DO iter=1, n_iter+1
            !
               DO i_sub=1,dim_subspace
                  CALL iotk_scan_dat(aux2_unit,"ITER"//TRIM(iotk_index(iter-1))//TRIM(iotk_index(i_sub)), state(iter,i_sub,:))
               ENDDO
            !
            ENDDO
            !
            CALL iotk_scan_end(aux2_unit,"STA")
            !
       CALL file_close(aux2_unit,PATH="/STATES/",ACTION="read")
         !
      CALL ioname('states',filename,LPATH=.FALSE.)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     PRINT*, 'Point 4'
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      WRITE( stdout,"(/,' Recursion states on WF basis read from file : ',5x,a)") TRIM(filename)

      WRITE( stdout,"(/,2x,'Calculate coeff  :')")



      coeff(:,:) = CZERO
      DO iter=1, n_iter+1
      !
         DO i_recur=1, dim_recursion
            !
            DO i_sub=1, dim_subspace
               !
               coeff(iter,i_recur) = coeff(iter,i_recur) + (state(iter,i_sub,i_recur)*CONJG(state(iter,i_sub,i_recur)))
               ! A essayer pour voir le terme de phase
               !coeff(iter,i_recur) = coeff(iter,i_recur) + state(iter,i_sub,i_recur)
               !
            ENDDO
         coeff(iter,i_recur)= sqrt(coeff(iter,i_recur))
         ENDDO
      ENDDO

     PRINT*, 'Point 5'


       ALLOCATE (  WF_like( nb_x, nb_y, nb_z), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating WF_like ', ABS(ierr) )
       ALLOCATE (  ix_2_x(nb_x), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating ix_2_x ', ABS(ierr) )
       ALLOCATE (  iy_2_y(nb_y), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating iy_2_y ', ABS(ierr) )
       ALLOCATE (  iz_2_z(nb_z), STAT=ierr ) 
               IF( ierr /=0 ) CALL errore(subname, 'allocating iz_2_z ', ABS(ierr) )


             !!!!!!!!!!!!!!!!!!Bricolage a la pierre
             !!!!!!!!!!!!!!!!! en premiere approx ca marche deja pas dans le reseau hexagonal
      DO i_x=1, nb_x
        ix_2_x(i_x) = (( (i_x - 1)*avec(1,1) ) + rmin(1)*avec(1,1))/ nb_x
      ENDDO
      DO i_y=1, nb_y
        iy_2_y(i_y) = (( (i_y - 1)*avec(2,2) ) + rmin(2)*avec(2,2))/ nb_y
      ENDDO
      DO i_z=1, nb_z
        iz_2_z(i_z) = (( (i_z - 1)*avec(3,3) ) + rmin(3)*avec(3,3)) / nb_z
      ENDDO


     PRINT*, 'Point 6'

!       DO i_x=1, nb_x
!          DO i_y=1, nb_y
!             DO i_z=1, nb_z
!                DO i_recur=1,dim_recursion
!                   WF_like(i_recur, i_x, i_y, i_z) = (EXP(-((ix_2_x(i_x) - wan_coord(1, i_recur))**2 + &
!                         (iy_2_y(i_y)-wan_coord(2, i_recur))**2 + (iz_2_z(i_z) -wan_coord(3, i_recur))**2 )/(2*sigma)))/(SQRT(sigma))
! 
!                ENDDO
!             ENDDO
!          ENDDO
!       ENDDO
! 
!       WRITE( stdout,"(/,' Gaussian functions constructed : ')")


      !
      ! output and internal formats
      ! when .plt is required, a .cube file is written and then converted to .plt
      !
      SELECT CASE ( TRIM(output_fmt) )
      CASE ( "txt")
           aux_fmt = ".txt"
      CASE ( "xsf" )
           aux_fmt = ".xsf"
      CASE ( "plt", "cube" )
           aux_fmt = ".cube"
      CASE DEFAULT
           CALL errore('plot','invalid OUTPUT_FMT '//TRIM(output_fmt),4)
      END SELECT


     ALLOCATE ( rwann_out(nb_x,nb_y,nb_z), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating rwann_out ', ABS(ierr) )

 DO iter = 1, n_iter+1

     rwann_out(:,:,:) = ZERO


      SELECT CASE ( TRIM(datatype) )
      CASE( "modulus" )    
            str = "_WFM"
            DO i_recur = 1, dim_recursion
               !
               WF_like(:,:,:)=ZERO
               DO i_x=1, nb_x
                  DO i_y=1, nb_y
                     DO i_z=1, nb_z
                           WF_like(i_x, i_y, i_z) = (EXP(-((ix_2_x(i_x) - wan_coord(1, i_recur))**2 + &
                                 (iy_2_y(i_y)-wan_coord(2, i_recur))**2 + (iz_2_z(i_z) -wan_coord(3, i_recur))**2 )/(2*sigma)))/(SQRT(sigma))

                     ENDDO
                  ENDDO
               ENDDO

               rwann_out(:,:,:) = rwann_out(:,:,:) + coeff(iter, i_recur) * WF_like(:,:,:) *  WF_like(:,:,:)
            ENDDO 
      CASE( "real" )    
         str = "_WFR"
         !DO i_recur = 1, dim_recursion
               !
              ! rwann_out(:,:,:) = rwann_out(:,:,:) + REAL (coeff(iter, i_recur) *  WF_like(:,:,:))
         !ENDDO 
      CASE( "imaginary" )    
         str = "_WFI"
        ! rwann_out(:,:,:) = ONE* AIMAG(coeff(iter, i_recur)  ) * WF_like(i_recur,:,:,:)
       !  DO i_recur = 1, dim_recursion
               !
              ! rwann_out(:,:,:) = rwann_out(:,:,:) + AIMAG (coeff(iter, i_recur) *  WF_like(i_recur,:,:,:))
        ! ENDDO 
      CASE( "modulnocoeff" )    
            str = "_WFM"
        !    DO i_recur = 1, dim_recursion
               !
               !IF (ABS(coeff(iter, i_recur)) >= EPS_m6) THEN
               !    rwann_out(:,:,:) = rwann_out(:,:,:) + ONE * REAL( WF_like(i_recur,:,:,:) * WF_like(i_recur,:,:,:) )
               !ENDIF
        !    ENDDO 
      CASE DEFAULT
         CALL errore('plot','invalid DATATYPE '//TRIM(datatype),3)
      END SELECT 



      filename=TRIM(work_dir)//"/"//TRIM(prefix)//TRIM(postfix)
      IF ( iter <= 9 ) THEN
         filename=TRIM(filename)//TRIM(str)//"00"//TRIM(int2char(iter))
      ELSE IF ( iter <= 99 ) THEN
         filename=TRIM(filename)//TRIM(str)//"0"//TRIM(int2char(iter))
      ELSE IF ( iter <= 999 ) THEN
         filename=TRIM(filename)//TRIM(str)//TRIM(int2char(iter))
      ELSE
         CALL errore('plot','iter > 999', iter)
      ENDIF
      !
      !
      WRITE( stdout,"(2x,'writing WF for iter (',i4,') on file: ',a)") &
            iter, TRIM(filename)//TRIM(aux_fmt)
      OPEN ( aux_unit, FILE=TRIM(filename)//TRIM(aux_fmt), FORM='formatted', &
                     STATUS='unknown', IOSTAT=ierr )
      IF (ierr/=0) CALL errore('plot','opening file '//TRIM(filename)//TRIM(aux_fmt),1)


      SELECT CASE ( TRIM(output_fmt) )
      CASE( "cube", "plt" )

         ! 
         ! bohr
         avecl(:,1) = avec(:,1) / nb_x
         avecl(:,2) = avec(:,2) / nb_y
         avecl(:,3) = avec(:,3) / nb_z


!          WRITE(aux_unit, '( " WanT" )') 
!          WRITE(aux_unit, '( " plot output - cube format" )' ) 
!          WRITE(aux_unit, '(i4,3f12.6)' ) natot, r0(:) 
!          WRITE(aux_unit, '(i4,3f12.6)' ) (nb_x),  avec(:,1) 
!          WRITE(aux_unit, '(i4,3f12.6)' ) (nryh-nryl+1),  avecl(:,2) 
!          WRITE(aux_unit, '(i4,3f12.6)' ) (nrzh-nrzl+1),  avecl(:,3) 
! 
!          DO ia = 1, natot
!             CALL atomic_name2num( symbtot(ia), zatom )
!             WRITE(aux_unit, '(i4,4e13.5)' ) zatom, ONE, tautot( :, ia )
!          ENDDO
! 
!          DO nx = nrxl, nrxh
!          DO ny = nryl, nryh
!             WRITE( aux_unit, "(6e13.5)" ) rwann_out( nx, ny, : )
!          ENDDO
!          ENDDO

      CASE( "txt" )

!          WRITE(aux_unit, '( " 3 2" )') 
!          WRITE(aux_unit, '( 3i5 )' ) nrzh-nrzl+1, nryh-nryl+1, nrxh-nrxl+1
!          WRITE(aux_unit, '(6f10.4)' ) r0(3) * bohr, r1(3) * bohr,  &
!                                        r0(2) * bohr, r1(2) * bohr,  & 
!                                        r0(1) * bohr, r1(1) * bohr
! 
!          DO nz = nrzl, nrzh
!          DO ny = nryl, nryh
!          DO nx = nrxl, nrxh
!             WRITE( aux_unit, "(f20.10)" ) rwann_out( nx, ny, nz )
!          ENDDO
!          ENDDO
!          ENDDO
! 
      CASE( "xsf" )

         ! 
         ! bohr
         avecl(:,1) = avec(:,1)
         avecl(:,2) = avec(:,2)
         avecl(:,3) = avec(:,3)
         
         !
         ! tau is temporarily converted to bohr 
         ! avec and tau passed in bohr, but converted to Ang in the routine
         !
         !tautot = tautot * alat
         CALL xsf_struct ( avec, nat, tautot, symbtot, aux_unit )
         !tau = tau / alat
         !
         CALL xsf_datagrid_3d ( rwann_out(:,:,:), nb_x, nb_y, nb_z,    &
                                    r0, avecl(:,1), avecl(:,2), avecl(:,3), aux_unit )

      CASE DEFAULT
         CALL errore('plot','invalid OUTPUT_FMT '//TRIM(output_fmt),5)
      END SELECT
      !
      CLOSE(aux_unit)

 ENDDO


      DEALLOCATE( ix_2_x, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating ix_2_x', ABS(ierr) )
      DEALLOCATE( iy_2_y, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating iy_2_y', ABS(ierr) )
      DEALLOCATE( iz_2_z, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating iz_2_z', ABS(ierr) )
      DEALLOCATE (  coeff, STAT=ierr ) 
         IF( ierr /=0 ) CALL errore(subname, 'deallocating coeff ', ABS(ierr) )

      DEALLOCATE( wan_coord, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating wan_coord', ABS(ierr) )
      DEALLOCATE( WF_like , STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating WF_like', ABS(ierr) )
      DEALLOCATE (tautot, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating tautot', ABS(ierr) )
      DEALLOCATE (symbtot, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating symbtot', ABS(ierr) )

      DEALLOCATE( rwann_out, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating rwann_out', ABS(ierr) )
!       DEALLOCATE( rec2tot, STAT=ierr )
!          IF( ierr /=0 ) CALL errore('plot', 'deallocating rec2tot ', ABS(ierr) )
      DEALLOCATE( state, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating state ', ABS(ierr) )
      DEALLOCATE( coeff, STAT=ierr )
         IF( ierr /=0 ) CALL errore('plot', 'deallocating coeff', ABS(ierr) )


END PROGRAM quick_plot
