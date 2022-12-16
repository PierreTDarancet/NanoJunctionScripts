!
! Copyright (C) 2007 Institut Néel Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE in_matrix_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE, EPS_m6, EPS_m4, EPS_m2
   USE io_global_module,     ONLY : stdout
   USE io_module,     ONLY :  mat_unit => aux_unit, ioname
   USE files_module, ONLY : file_open, file_close
   !USE parser_module
   USE iotk_module
   USE T_control_module,      ONLY : in_datafile,       &
                                     dim_subspace,      &
                                     in_max_iter
 



   IMPLICIT NONE
   PRIVATE 
   SAVE

   ! variables for recursion
            ! For first iteration
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix(:,:,:), B_matrix(:,:,:)
    !
    !
   INTEGER :: in_dim
    !
   LOGICAL :: alloc

    !
   PUBLIC :: A_matrix, B_matrix
   PUBLIC :: in_dim
   PUBLIC :: alloc
   !PUBLIC :: total_iter
    !


   ! general functions
   PUBLIC :: in_matrix_allocate
   PUBLIC :: in_matrix_deallocate
   PUBLIC :: read_matrix
    !
    !



CONTAINS



!*******************************************************************
   SUBROUTINE in_matrix_allocate
   !*******************************************************************
      CHARACTER(18)      :: subname="in_matrix_allocate"
      INTEGER :: ierr


    ! Allocate hamiltonian

   !
   in_dim = in_max_iter*dim_subspace
   !
   !total_iter = in_max_iter_R + in_max_iter_L - 1

    ! Allocate A-type matrices
   ALLOCATE ( A_matrix(in_max_iter,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix', ABS(ierr) )
    ! Allocate B-type matrices
   ALLOCATE ( B_matrix(in_max_iter,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix', ABS(ierr) )
   ! set alloc to true
    alloc=.TRUE.

END SUBROUTINE in_matrix_allocate

!**********************************************************
   SUBROUTINE in_matrix_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(20)      :: subname="in_matrix_deallocate"
      INTEGER :: ierr


     IF ( ALLOCATED( A_matrix  ) ) THEN
           DEALLOCATE ( A_matrix, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix  ) ) THEN
           DEALLOCATE ( B_matrix, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix', ABS(ierr) )
      ENDIF

     alloc=.FALSE.

   END SUBROUTINE in_matrix_deallocate


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


       CALL file_open(mat_unit,TRIM(in_datafile),PATH="/MATRIX/",ACTION="read", &
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
       IF (final_iter+1 /= in_max_iter) CALL errore(subname, 'read n_iter and given n_iter are /=', ABS(final_iter+1) )
       !

       CALL iotk_scan_begin(mat_unit,"A_MAT")
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), A_matrix(1,:,:))
       DO iter=1, final_iter
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), A_matrix(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"A_MAT")
       !
       !
       CALL iotk_scan_begin(mat_unit,"B_MAT")
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), B_matrix(1,:,:))
       DO iter=1, final_iter
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), B_matrix(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"B_MAT")

      ! CALL iotk_scan_end(mat_unit,TRIM(name))

      CALL file_close(mat_unit,PATH="/MATRIX/",ACTION="read")

      !CALL ioname('matrix',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Matrices on recursion basis read from file : ',5x,a)") TRIM(in_datafile)
      ! End of right part

  END SUBROUTINE read_matrix


 END MODULE in_matrix_module


