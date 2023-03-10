!
! Copyright (C) 2007 Institut N?el Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE convert_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE, EPS_m6, EPS_m4, EPS_m2
   !USE io_module, ONLY : title, prefix, postfix, work_dir
   USE io_global_module,     ONLY : stdin, stdout
   USE io_module,     ONLY :  ham_unit, mat_unit => aux_unit, &
                             sta_unit => aux2_unit, cen_unit => aux3_unit, ioname
   USE files_module, ONLY : file_open, file_close
   !USE parser_module
   USE iotk_module
   USE convert_control_variable_module,      ONLY : out_datafile_C, &
                                            out_datafile_L, &
                                            out_datafile_R, &
                                            in_datafile_L, &
                                            in_datafile_R, &
                                            dim_subspace, &
                                            n_iter_R,     &
                                            n_iter_L



   IMPLICIT NONE
   PRIVATE 
   SAVE

   ! variables for recursion
            ! For first iteration
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_L(:,:,:), B_matrix_L(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_R(:,:,:), B_matrix_R(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix(:,:,:), B_matrix(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  final_hamiltonian_C(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  final_hamiltonian_R(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  final_hamiltonian_L(:,:,:)

   INTEGER :: dimC
   INTEGER :: dimLR


   PUBLIC :: A_matrix_L, B_matrix_L
   PUBLIC :: A_matrix_R, B_matrix_R
   PUBLIC :: A_matrix, B_matrix
   PUBLIC :: final_hamiltonian_C
   PUBLIC :: final_hamiltonian_R
   PUBLIC :: final_hamiltonian_L

   ! general functions
   PUBLIC :: matrix_allocate
   PUBLIC :: matrix_deallocate
   PUBLIC :: print_want
   PUBLIC :: read_matrix
   PUBLIC :: build_hamiltonian
   PUBLIC :: build_matrix
   PUBLIC :: convert_summary_input


    !
   PUBLIC :: dimC
   PUBLIC :: dimLR

CONTAINS



!*******************************************************************
   SUBROUTINE matrix_allocate
   !*******************************************************************
      CHARACTER(15)      :: subname="matrix_allocate"
      INTEGER :: ierr


    ! Allocate hamiltonian

   dimLR = dim_subspace
   dimC = (n_iter_R+n_iter_L-3)*dim_subspace
   ALLOCATE ( final_hamiltonian_R(2,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating final_hamiltonian_R', ABS(ierr) )
   ALLOCATE ( final_hamiltonian_L(2,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating final_hamiltonian_L', ABS(ierr) )
   ALLOCATE ( final_hamiltonian_C(3,dimC,dimC), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating final_hamiltonian_C', ABS(ierr) )


    ! Allocate A-type matrices
   ALLOCATE ( A_matrix_L(n_iter_L,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_L', ABS(ierr) )
   ALLOCATE ( A_matrix_R(n_iter_R,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_R', ABS(ierr) )
   ALLOCATE ( A_matrix((n_iter_R+n_iter_L-1),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix', ABS(ierr) )
    ! Allocate B-type matrices
   ALLOCATE ( B_matrix_L(n_iter_L,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_L', ABS(ierr) )
   ALLOCATE ( B_matrix_R(n_iter_R,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_R', ABS(ierr) )
   ALLOCATE ( B_matrix((n_iter_R+n_iter_L),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix', ABS(ierr) )
   ! set alloc to true
    !alloc=.TRUE.

END SUBROUTINE matrix_allocate

!**********************************************************
   SUBROUTINE matrix_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(17)      :: subname="matrix_deallocate"
      INTEGER :: ierr

     IF ( ALLOCATED( final_hamiltonian_R ) ) THEN
           DEALLOCATE ( final_hamiltonian_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating final_hamiltonian_R', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( final_hamiltonian_L ) ) THEN
           DEALLOCATE ( final_hamiltonian_L, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating final_hamiltonian_L', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( final_hamiltonian_C ) ) THEN
           DEALLOCATE ( final_hamiltonian_C, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating final_hamiltonian_C', ABS(ierr) )
      ENDIF


     IF ( ALLOCATED( A_matrix_L  ) ) THEN
           DEALLOCATE ( A_matrix_L, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_L', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_L  ) ) THEN
           DEALLOCATE ( B_matrix_L, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_L', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( A_matrix_R  ) ) THEN
           DEALLOCATE ( A_matrix_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_R', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_R  ) ) THEN
           DEALLOCATE ( B_matrix_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_R', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( A_matrix  ) ) THEN
           DEALLOCATE ( A_matrix, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix  ) ) THEN
           DEALLOCATE ( B_matrix, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix', ABS(ierr) )
      ENDIF

   END SUBROUTINE matrix_deallocate


!*******************************************************************
   SUBROUTINE build_hamiltonian
   !*******************************************************************
       CHARACTER(17)      :: subname="build_hamiltonian"
       INTEGER :: icols, irows 
       INTEGER :: iter
       INTEGER :: icols_shift, irows_shift
       INTEGER :: range_iter
     ! Beginning

     ! set Hamiltonians to ZERO
     final_hamiltonian_L(:,:,:)=CZERO
     final_hamiltonian_C(:,:,:)=CZERO
     final_hamiltonian_R(:,:,:)=CZERO

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEFT PART
     ! generating the H_00L & H_01L
     ! ! Warning : only H_01L is read by WanT, we have to invert the B matrix another time
     ! ! should be improved 

     ! H00L
     final_hamiltonian_L(1,:,:)=A_matrix(1,:,:)
     !
     ! H01L
     ! ! Different for Right& Left ! !

     DO icols=1, dim_subspace
        DO irows=1, dim_subspace
           final_hamiltonian_L(2,irows,icols)=CONJG(B_matrix(1,icols,irows))
        ENDDO
     ENDDO
     !

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RIGHT PART
     ! generating the H_00R & H_01R
     !set iter as the last element of A & B Matrix
     iter = n_iter_L+ n_iter_R - 1
     ! H00R
     final_hamiltonian_R(1,:,:)=A_matrix(iter,:,:)
     !
     ! H01R
     ! ! Different for Right& Left ! !
     iter = n_iter_L+ n_iter_R 
     final_hamiltonian_R(2,:,:)=B_matrix(iter,:,:)
     !

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONDUCTOR PART
     ! generating the H_LC & H_CR
     !
     ! H_LC
     ! ! set iter as the second element of the B matrix
     iter = 2

     !Suppression du 21/03 Test sachant que ca marchait pas avant
         DO icols=1, dim_subspace
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ? rev?rifier pour le placement ligne colonne
            !!!!!!!!!!!!!!!!!!!!!!!!
     !       icols_shift=icols + (n_iter_R+n_iter_L-4)*dim_subspace
               icols_shift= icols

            DO irows=1, dim_subspace
   !            irows_shift=irows
               irows_shift= irows + (n_iter_R+n_iter_L-4)*dim_subspace

               final_hamiltonian_C(1,irows_shift,icols_shift)=CONJG(B_matrix(iter,icols,irows))
            ENDDO
         ENDDO
     !


     !
!      DO icols=1, dim_subspace
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! ? rev?rifier pour le placement ligne colonne
!         !!!!!!!!!!!!!!!!!!!!!!!!
!         icols_shift=icols
!         DO irows=1, dim_subspace
!            irows_shift=irows
!            final_hamiltonian_C(1,irows_shift,icols_shift)=B_matrix(iter,irows,icols)
!         ENDDO
!      ENDDO


     !
     ! H_CR
     ! ! set iter as the last but one element of the B matrix
     iter = n_iter_L+ n_iter_R - 1
     !
     !Suppression du 21/03 Test sachant que ca marchait pas avant

            DO icols=1, dim_subspace
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! ? rev?rifier pour le placement ligne colonne
               !!!!!!!!!!!!!!!!!!!!!!!!
               icols_shift=icols + (n_iter_R+n_iter_L-4)*dim_subspace
       !        icols_shift= icols
               DO irows=1, dim_subspace
                   irows_shift=irows
       ! irows_shift= irows + (n_iter_R+n_iter_L-4)*dim_subspace

                  final_hamiltonian_C(3,irows_shift,icols_shift)=B_matrix(iter,irows,icols)
               ENDDO
            ENDDO
     !

     !
!      DO icols=1, dim_subspace
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! ? rev?rifier pour le placement ligne colonne
!         !!!!!!!!!!!!!!!!!!!!!!!!
!         icols_shift= icols + (n_iter_R+n_iter_L-4)*dim_subspace
!         DO irows=1, dim_subspace
!            irows_shift= irows + (n_iter_R+n_iter_L-4)*dim_subspace
!            final_hamiltonian_C(3,irows_shift,icols_shift)=B_matrix(iter,irows,icols)
!         ENDDO
!      ENDDO

     ! generating the H_C
     !

     ! generating the off diagonal terms (B_matrix)
     ! last element for the B matrices
     range_iter= n_iter_L + n_iter_R - 2
     DO iter=3,range_iter
      DO icols=1, dim_subspace
         icols_shift=dim_subspace*(iter-2)+icols
         DO irows=1, dim_subspace
            irows_shift=dim_subspace*(iter-3)+irows
            final_hamiltonian_C(2,irows_shift,icols_shift)=B_matrix(iter,irows,icols)
            final_hamiltonian_C(2,icols_shift,irows_shift)=CONJG(B_matrix(iter,irows,icols))
            !
         ENDDO
      ENDDO
     ENDDO


     ! generating the diagonal terms (A_matrix)
     range_iter= n_iter_L + n_iter_R - 2
     DO iter=2,range_iter
      DO icols=1, dim_subspace
         icols_shift=dim_subspace*(iter-2)+icols
         DO irows=1, dim_subspace
            irows_shift=dim_subspace*(iter-2)+irows
            !
            final_hamiltonian_C(2,irows_shift,icols_shift)=A_matrix(iter,irows,icols)
            !
         ENDDO
      ENDDO
     ENDDO


  END SUBROUTINE build_hamiltonian


!*******************************************************************
   SUBROUTINE build_matrix
   !*******************************************************************

     ! This program builds the final ordered A_matrix & B_matrix 
     ! starting from the knowledge of the A/B_matrix_L/R read from recursion's files
       CHARACTER(12)      :: subname="build_matrix"
       LOGICAL   ::   test_equality
       INTEGER   ::   irows, icols
       INTEGER   ::   iter, shift_iter


     !Beginning
    WRITE(stdout,"(2x,70('='))" )
    WRITE(stdout,"(/,  7x,'Building ordered A & B matrices ')")
    WRITE(stdout,"(2x,70('='))" )
    WRITE(stdout,"()")
    WRITE(stdout,"()")

    !INitialize the final matrices : A_matrix & B_matrix
    A_matrix(:,:,:)=CZERO
    B_matrix(:,:,:)=CZERO
    !test equality for the first element of matrices R/L
    test_equality=.FALSE.
    DO icols=1, dim_subspace
      DO irows=1, dim_subspace
         IF (ABS(A_matrix_L(1,irows,icols) - A_matrix_R(1,irows,icols)) > EPS_m4 ) &
                test_equality=.TRUE.

!                   PRINT*, 'A_matrix_L(1,i,j)', irows, icols
!                   PRINT*, A_matrix_L(1,irows,icols)
! 
!                   PRINT*, 'A_matrix_R(1,i,j)', irows, icols
!                   PRINT*, A_matrix_R(1,irows,icols)

      ENDDO
    ENDDO


    IF (test_equality)  CALL errore(subname, 'initial subspace in L/R are different', 1 )

   ! initialize the first (Left) part of the A matrix
   ! initialize the first (Left) part of the B matrix=CONJG(TRANSPOSE(B_L))

    DO iter=1,n_iter_L
      ! corresponding iter in the R/L matrix 
      shift_iter = n_iter_L - iter + 1

      A_matrix(iter,:,:) = A_matrix_L(shift_iter,:,:)
      DO icols=1, dim_subspace
         DO irows=1, dim_subspace
            B_matrix(iter,irows,icols) = CONJG(B_matrix_L(shift_iter,icols,irows))
         ENDDO
      ENDDO

    ENDDO


    ! new test on the general A matrix
    iter=n_iter_L
    test_equality=.FALSE.
    DO icols=1, dim_subspace
      DO irows=1, dim_subspace
         IF (ABS(A_matrix(iter,irows,icols) - A_matrix_R(1,irows,icols)) > EPS_m4) &
                test_equality=.TRUE.

      ENDDO
    ENDDO
    IF (test_equality)  CALL errore(subname, 'initial subspace in L/R are different', 2 )



   ! initialize the second (Right) part of the A matrix
   ! initialize the second (Right) part of the B matrix=B_R
    shift_iter = n_iter_L 
    B_matrix(shift_iter+1,:,:) = B_matrix_R(1,:,:)
    DO iter=2,n_iter_R
      ! corresponding iter in the R/L matrix 
      shift_iter = n_iter_L + iter - 1 

      A_matrix(shift_iter,:,:) = A_matrix_R(iter,:,:)
      B_matrix(shift_iter+1,:,:) = B_matrix_R(iter,:,:)

    ENDDO
    WRITE(stdout,"()")
    WRITE(stdout,"()")
    WRITE(stdout,"(2x,70('='))" )
    WRITE( stdout,"(/, 7x,' Ordered A & B matrices built ')")
    WRITE(stdout,"(2x,70('='))" )

   END SUBROUTINE build_matrix



!**********************************************************
   SUBROUTINE print_want
   !**********************************************************
   IMPLICIT NONE
       CHARACTER(nstrx)   :: attr
       CHARACTER(nstrx)  :: name
       CHARACTER(17)      :: subname="HAMILTONIAN_WRITE"
       INTEGER            :: nrtot_C, nrtot_LR
       INTEGER            :: nr_C(3), nr_LR(3)
       INTEGER            :: ivr_C(3,3), ivr_LR(3,2)
       INTEGER            :: ik, ir, ierr




     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! initialization PART
       name='HAMILTONIAN'

       nrtot_C = 3
       nr_C(1) = 3
       nr_C(2) = 1
       nr_C(3) = 1 

       nrtot_LR = 2
       nr_LR(1) = 2
       nr_LR(2) = 1
       nr_LR(3) = 1 


       !ivr_C(3,nrtot)
       ivr_C(:,:)=0
       ivr_C(1,1)=-1
       ivr_C(1,2)=0
       ivr_C(1,3)=1

       !ivr_LR(3,)
       ivr_LR(:,:)=0
       ivr_LR(1,1)=0
       ivr_LR(1,2)=1


      CALL file_open(ham_unit,TRIM(out_datafile_C),PATH="/",ACTION="write", FORM='formatted')
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONDUCTOR PART
         !
         CALL iotk_write_begin(ham_unit,TRIM(name))
         CALL iotk_write_attr(attr,"dimwann",dimC,FIRST=.TRUE.) 
         CALL iotk_write_attr(attr,"nkpts",nrtot_C) 
         CALL iotk_write_attr(attr,"nk",nr_C) 
         CALL iotk_write_attr(attr,"nrtot",nrtot_C) 
         CALL iotk_write_attr(attr,"nr",nr_C) 
         CALL iotk_write_empty(ham_unit,"DATA",ATTR=attr)
               !
         CALL iotk_write_dat(ham_unit,"IVR", ivr_C, ATTR=attr, COLUMNS=3, IERR=ierr) 
               IF (ierr/=0) CALL errore(subname,'writing ivr',ABS(ierr))
         !
         CALL iotk_write_begin(ham_unit,"RHAM")
         DO ir = 1, nrtot_C
               CALL iotk_write_dat(ham_unit,"VR"//TRIM(iotk_index(ir)), final_hamiltonian_C(ir,:,:))
         ENDDO
         CALL iotk_write_end(ham_unit,"RHAM")
         CALL iotk_write_end(ham_unit,TRIM(name))
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL file_close(ham_unit,PATH="/",ACTION="write")
      WRITE( stdout,"(/,'  Hamiltonian of the Central region written on file : ',a)") TRIM(out_datafile_C)

      CALL file_open(ham_unit,TRIM(out_datafile_R),PATH="/",ACTION="write", FORM='formatted')
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RIGHT PART
       !
       CALL iotk_write_begin(ham_unit,TRIM(name))
       CALL iotk_write_attr(attr,"dimwann",dim_subspace,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"nkpts",nrtot_LR) 
       CALL iotk_write_attr(attr,"nk",nr_LR) 
       CALL iotk_write_attr(attr,"nrtot",nrtot_LR) 
       CALL iotk_write_attr(attr,"nr",nr_LR) 
       CALL iotk_write_empty(ham_unit,"DATA",ATTR=attr)
            !
       CALL iotk_write_dat(ham_unit,"IVR", ivr_LR, ATTR=attr, COLUMNS=3, IERR=ierr) 
            IF (ierr/=0) CALL errore(subname,'writing ivr',ABS(ierr))
       !
       CALL iotk_write_begin(ham_unit,"RHAM")
       DO ir = 1, nrtot_LR
             CALL iotk_write_dat(ham_unit,"VR"//TRIM(iotk_index(ir)), final_hamiltonian_R(ir,:,:))
       ENDDO
       CALL iotk_write_end(ham_unit,"RHAM")
       CALL iotk_write_end(ham_unit,TRIM(name))
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL file_close(ham_unit,PATH="/",ACTION="write")
      WRITE( stdout,"(/,'  Hamiltonian of the right region written on file : ',a)") TRIM(out_datafile_R)


      CALL file_open(ham_unit,TRIM(out_datafile_L),PATH="/",ACTION="write", FORM='formatted')
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEFT PART
       !
       CALL iotk_write_begin(ham_unit,TRIM(name))
       CALL iotk_write_attr(attr,"dimwann",dim_subspace,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"nkpts",nrtot_LR) 
       CALL iotk_write_attr(attr,"nk",nr_LR) 
       CALL iotk_write_attr(attr,"nrtot",nrtot_LR) 
       CALL iotk_write_attr(attr,"nr",nr_LR) 
       CALL iotk_write_empty(ham_unit,"DATA",ATTR=attr)
            !
       CALL iotk_write_dat(ham_unit,"IVR", ivr_LR, ATTR=attr, COLUMNS=3, IERR=ierr) 
            IF (ierr/=0) CALL errore(subname,'writing ivr',ABS(ierr))
       !
       CALL iotk_write_begin(ham_unit,"RHAM")
       DO ir = 1, nrtot_LR
             CALL iotk_write_dat(ham_unit,"VR"//TRIM(iotk_index(ir)), final_hamiltonian_L(ir,:,:))
       ENDDO
       CALL iotk_write_end(ham_unit,"RHAM")
       CALL iotk_write_end(ham_unit,TRIM(name))
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL file_close(ham_unit,PATH="/",ACTION="write")
      WRITE( stdout,"(/,'  Hamiltonian of the left region written on file : ',a)") TRIM(out_datafile_L)


   END SUBROUTINE print_want


!*******************************************************************
   SUBROUTINE read_matrix
   !*******************************************************************
       CHARACTER(11)      :: subname="read_matrix"
       CHARACTER(6) :: name="MATRIX"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: iter
       INTEGER            :: ierr
       INTEGER         :: final_iter, dim_sub
      CHARACTER( LEN=nstrx )  :: filename


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RIGHT PART
       CALL file_open(mat_unit,TRIM(in_datafile_R),PATH="/MATRIX/",ACTION="read", &
                              FORM='formatted')


       CALL iotk_scan_empty(mat_unit,"DATADIM",ATTR=attr)
       !CALL iotk_scan_attr(attr,"dimwan_total",dimwan_total,FIRST=.TRUE.) 
       !CALL iotk_scan_attr(attr,"dim_rec",dim_rec) 
       CALL iotk_scan_attr(attr,"dim_sub",dim_sub) 
       IF (dim_sub /= dim_subspace) CALL errore(subname, 'dim_subspace in input file and in .mat file are /=', ABS(dim_sub) )


       CALL iotk_scan_empty(mat_unit,"DATAFINAL",ATTR=attr)
       !CALL iotk_scan_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
       CALL iotk_scan_attr(attr,"final_iter",final_iter) 
       !CALL iotk_scan_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
       IF (final_iter+1 /= n_iter_R) CALL errore(subname, 'read n_iter_R and given n_iter_R are /=', ABS(final_iter+1) )
       !

       CALL iotk_scan_begin(mat_unit,"A_MAT")
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), A_matrix_R(1,:,:))
       DO iter=1, final_iter
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), A_matrix_R(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"A_MAT")
       !
       !
       CALL iotk_scan_begin(mat_unit,"B_MAT")
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), B_matrix_R(1,:,:))
       DO iter=1, final_iter
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), B_matrix_R(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"B_MAT")

      ! CALL iotk_scan_end(mat_unit,TRIM(name))

      CALL file_close(mat_unit,PATH="/MATRIX/",ACTION="read")

      !CALL ioname('matrix',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Matrices on recursion basis read from file : ',5x,a)") TRIM(in_datafile_R)
      ! End of right part

      ! reinitialize local variables
      final_iter = 0
      dim_sub = 0
      !

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEFT PART

       CALL file_open(mat_unit,TRIM(in_datafile_L),PATH="/MATRIX/",ACTION="read", &
                              FORM='formatted')


       CALL iotk_scan_empty(mat_unit,"DATADIM",ATTR=attr)
       !CALL iotk_scan_attr(attr,"dimwan_total",dimwan_total,FIRST=.TRUE.) 
       !CALL iotk_scan_attr(attr,"dim_rec",dim_rec) 
       CALL iotk_scan_attr(attr,"dim_sub",dim_sub) 
       IF (dim_sub /= dim_subspace) CALL errore(subname, 'dim_subspace in input file and in .mat file are /=', ABS(dim_sub) )


       CALL iotk_scan_empty(mat_unit,"DATAFINAL",ATTR=attr)
       !CALL iotk_scan_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
       CALL iotk_scan_attr(attr,"final_iter",final_iter) 
       !CALL iotk_scan_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
       IF (final_iter+1 /= n_iter_L) CALL errore(subname, 'read n_iter_L and given n_iter_L are /=', ABS(final_iter+1) )
       !

       CALL iotk_scan_begin(mat_unit,"A_MAT")
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), A_matrix_L(1,:,:))
       DO iter=1, final_iter
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), A_matrix_L(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"A_MAT")
       !
       !
       CALL iotk_scan_begin(mat_unit,"B_MAT")
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), B_matrix_L(1,:,:))
       DO iter=1, final_iter
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), B_matrix_L(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"B_MAT")

       !CALL iotk_scan_end(mat_unit,TRIM(name))

      CALL file_close(mat_unit,PATH="/MATRIX/",ACTION="read")

      !CALL ioname('matrix',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Matrices on recursion basis read from file : ',5x,a)") TRIM(in_datafile_L)



  END SUBROUTINE read_matrix

!**********************************************************
   SUBROUTINE convert_summary_input
   !**********************************************************
   ! 
   ! Print out all the informnatins obtained from the 
   ! input and initialization routines.
   !
  !
   ! input variables
   !
   INTEGER      :: unit
   !
   ! local variables
   !

!--------------------------------------------------------
   unit = stdout
   !
   ! <INPUT> section
   !
   WRITE(unit,"()")      
   WRITE(unit,"()")      
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',32x,'Main',32x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )
!       WRITE(unit,"(  7x,'     Calculation Title :',5x,a)") TRIM(title)
!       WRITE(unit,"(  7x,'                Prefix :',5x,a)") TRIM(prefix)
!       WRITE(unit,"(  7x,'               Postfix :',5x,a)") TRIM(postfix)
!          IF ( LEN_TRIM(work_dir) <= 65 ) THEN
!             WRITE(unit,"(  7x,'     Working directory :',5x,a)") TRIM(work_dir)
!          ELSE
!             WRITE(unit,"(  7x,'     Working directory :',5x,/,10x,a)") TRIM(work_dir)
!          ENDIF
   WRITE(unit,"()")
   WRITE(unit,"( /,2x,'<INPUT>')" )
   WRITE(unit,"(  7x,'iteration number for the left  part:',5x,i4)") n_iter_L
   WRITE(unit,"(  7x,'iteration number for the right part:',5x,i4)") n_iter_R

   WRITE(unit,"(  7x,'Hamiltonian data in the left  part read from file  :',5x,a)") TRIM(in_datafile_L)
   WRITE(unit,"(  7x,'Hamiltonian data in the right part read from file  :',5x,a)") TRIM(in_datafile_R)


   WRITE(unit,"(  7x,'WanT-like Hamiltonian data in the left       part written in file  :',5x,a)") TRIM(out_datafile_L)
   WRITE(unit,"(  7x,'WanT-like Hamiltonian data in the right      part written in file  :',5x,a)") TRIM(out_datafile_R)
   WRITE(unit,"(  7x,'WanT-like Hamiltonian data in the conductor  part written in file  :',5x,a)") TRIM(out_datafile_C)


   WRITE(unit,"(  7x,'Dimension of the subspace :',5x,i4)") dim_subspace


   WRITE(unit,"()")      
   WRITE(unit,"()")      

   END SUBROUTINE convert_summary_input





 END MODULE convert_module


