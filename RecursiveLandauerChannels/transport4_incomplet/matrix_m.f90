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
   MODULE matrix_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE, EPS_m6, EPS_m4, EPS_m2
   !USE io_module, ONLY : title, prefix, postfix, work_dir
   USE io_global_module,     ONLY : stdout
   !USE parser_module
   USE iotk_module
   USE T_control_module,      ONLY : dim_subspace, &
                                              max_iter_R,   &
                                              max_iter_L
   USE in_matrix_module,              ONLY : A_matrix, B_matrix,   &
                                             total_iter




   IMPLICIT NONE
   PRIVATE 
   SAVE

   ! variables for recursion
            ! For first iteration
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_L(:,:,:), B_matrix_L(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_R(:,:,:), B_matrix_R(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_C(:,:,:), B_matrix_C(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_CR(:,:,:), B_matrix_CR(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_LC(:,:,:), B_matrix_LC(:,:,:)
   !

   LOGICAL :: alloc

   PUBLIC :: A_matrix_L, B_matrix_L
   PUBLIC :: A_matrix_R, B_matrix_R
   PUBLIC :: A_matrix_C, B_matrix_C
   PUBLIC :: A_matrix_LC, B_matrix_LC
   PUBLIC :: A_matrix_CR, B_matrix_CR
   PUBLIC :: alloc


   ! general functions
   PUBLIC :: matrix_allocate
   PUBLIC :: matrix_deallocate
   PUBLIC :: build_matrix


CONTAINS



!*******************************************************************
   SUBROUTINE matrix_allocate
   !*******************************************************************
      CHARACTER(15)      :: subname="matrix_allocate"
      INTEGER :: ierr, max_iter_C

   max_iter_C = total_iter -  max_iter_R -  max_iter_L

    ! Allocate A-type matrices
   ALLOCATE ( A_matrix_L(max_iter_L,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_L', ABS(ierr) )
   ALLOCATE ( A_matrix_R(max_iter_R,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_R', ABS(ierr) )
   ALLOCATE ( A_matrix_C(max_iter_C,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_C', ABS(ierr) )
   ALLOCATE ( A_matrix_CR(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_CR', ABS(ierr) )
   ALLOCATE ( A_matrix_LC(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_LC', ABS(ierr) )

    ! Allocate B-type matrices
   ALLOCATE ( B_matrix_L(max_iter_L,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_L', ABS(ierr) )
   ALLOCATE ( B_matrix_R(max_iter_R,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_R', ABS(ierr) )
   ALLOCATE ( B_matrix_C((max_iter_C - 1),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_C', ABS(ierr) )
   ALLOCATE ( B_matrix_CR(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_CR', ABS(ierr) )
   ALLOCATE ( B_matrix_LC(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_LC', ABS(ierr) )


   ! set alloc to true
    alloc=.TRUE.

END SUBROUTINE matrix_allocate

!**********************************************************
   SUBROUTINE matrix_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(17)      :: subname="matrix_deallocate"
      INTEGER :: ierr

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
     IF ( ALLOCATED( A_matrix_C  ) ) THEN
           DEALLOCATE ( A_matrix_C, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_C', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_C  ) ) THEN
           DEALLOCATE ( B_matrix_C, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_C', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( A_matrix_CR  ) ) THEN
           DEALLOCATE ( A_matrix_CR, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_CR', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_CR  ) ) THEN
           DEALLOCATE ( B_matrix_CR, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_CR', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( A_matrix_LC  ) ) THEN
           DEALLOCATE ( A_matrix_LC, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_LC', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_LC  ) ) THEN
           DEALLOCATE ( B_matrix_LC, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_LC', ABS(ierr) )
      ENDIF



     alloc = .FALSE.

   END SUBROUTINE matrix_deallocate

!*******************************************************************
   SUBROUTINE build_matrix
   !*******************************************************************

     ! This program builds the final ordered A_matrix & B_matrix 
     ! starting from the knowledge of the A/B_matrix_L/R read from recursion's files
       CHARACTER(13)      :: subname="build_matrix"
       INTEGER   ::   irows, icols
       INTEGER   ::   iter, shift_iter, max_iter_C


     !Beginning
    WRITE(stdout,"(2x,70('='))" )
    WRITE(stdout,"(/,  7x,'Building A & B matrices for recursion input')")
    WRITE(stdout,"(2x,70('='))" )
    WRITE(stdout,"()")
    WRITE(stdout,"()")

    !test equality for the first element of matrices R/L


    ! initialize the first (Left) part of the A matrix
    ! initialize the first (Left) part of the B matrix=CONJG(TRANSPOSE(B_L))
    A_matrix_R(:,:,:) = CZERO
    B_matrix_R(:,:,:) = CZERO
    !
    A_matrix_L(:,:,:) = CZERO
    B_matrix_L(:,:,:) = CZERO

    DO iter=1,max_iter_L
      ! corresponding iter in the R/L matrix 
      shift_iter = max_iter_L - iter + 1

      A_matrix_L(shift_iter,:,:) = A_matrix(iter,:,:) 
      !
      DO icols=1, dim_subspace
         !
         DO irows=1, dim_subspace
               !
               B_matrix_L(shift_iter,irows,icols) = CONJG(B_matrix(iter,icols,irows) )
               !
         ENDDO
         !
      ENDDO
      !
    ENDDO


 
   ! initialize the second (Right) part of the A matrix
   ! initialize the second (Right) part of the B matrix=B_R
    DO iter=1,max_iter_R
      ! corresponding iter in the R/L matrix 
      shift_iter = total_iter - max_iter_R + iter 

      A_matrix_R(iter,:,:) = A_matrix(shift_iter,:,:)
      DO icols=1, dim_subspace
          !
          DO irows=1, dim_subspace
              B_matrix_R(iter,irows,icols)= B_matrix(shift_iter+1,irows,icols)
          ENDDO
          !
      ENDDO

    ENDDO
    !
    A_matrix_C(:,:,:) = CZERO
    B_matrix_C(:,:,:) = CZERO
    !
    A_matrix_LC(:,:,:) = CZERO
    B_matrix_LC(:,:,:) = CZERO
    !
    A_matrix_CR(:,:,:) = CZERO
    B_matrix_CR(:,:,:) = CZERO
    !
    shift_iter = max_iter_L +1
    DO icols=1, dim_subspace
          !
          DO irows=1, dim_subspace
              B_matrix_LC(1,irows,icols) = B_matrix(shift_iter,irows,icols)
          ENDDO
          !
    ENDDO


    shift_iter = total_iter - max_iter_R +1
    DO icols=1, dim_subspace
          !
          DO irows=1, dim_subspace
              B_matrix_CR(1,irows,icols) = B_matrix(shift_iter,irows,icols)
          ENDDO
          !
    ENDDO


    max_iter_C = total_iter -  max_iter_R -  max_iter_L
    DO iter=1,max_iter_C
         !
         shift_iter =  max_iter_L + iter 
         !
         A_matrix_C(iter,:,:) = A_matrix(shift_iter,:,:)
         !
    ENDDO
    DO iter=1,max_iter_C-1
         shift_iter = max_iter_L + iter 
         DO icols=1, dim_subspace
            !
            DO irows=1, dim_subspace
               B_matrix_C(iter,irows,icols)= B_matrix(shift_iter+1,irows,icols)
            ENDDO
            !
         ENDDO
         !
    ENDDO



    WRITE(stdout,"()")
    WRITE(stdout,"()")
    WRITE(stdout,"(2x,70('='))" )
    WRITE( stdout,"(/, 7x,'A & B matrices for recursion built ')")
    WRITE(stdout,"(2x,70('='))" )


   END SUBROUTINE build_matrix


 END MODULE matrix_module


