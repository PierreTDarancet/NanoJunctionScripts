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
                                     max_iter_L,   &
                                     max_iter_CR,   &
                                     max_iter_CL,   &
                                     max_iter_CC

   USE in_matrix_module,              ONLY : A_matrix, B_matrix




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
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_CL(:,:,:), B_matrix_CL(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  B_matrix_L_CL(:,:,:), B_matrix_CR_R(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  B_matrix_CL_C(:,:,:), B_matrix_C_CR(:,:,:)
   !

   LOGICAL :: alloc

   PUBLIC :: A_matrix_L, B_matrix_L
   PUBLIC :: A_matrix_R, B_matrix_R
   PUBLIC :: A_matrix_C, B_matrix_C
   PUBLIC :: A_matrix_CL, B_matrix_CL
   PUBLIC :: A_matrix_CR, B_matrix_CR
   PUBLIC :: B_matrix_CR_R, B_matrix_L_CL
   PUBLIC :: B_matrix_C_CR, B_matrix_CL_C

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

!   max_iter_C = total_iter -  max_iter_R -  max_iter_L

    ! Allocate A-type matrices
   ALLOCATE ( A_matrix_L(max_iter_L,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_L', ABS(ierr) )
   ALLOCATE ( A_matrix_R(max_iter_R,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_R', ABS(ierr) )
   ALLOCATE ( A_matrix_C(max_iter_CC,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_C', ABS(ierr) )
   ALLOCATE ( A_matrix_CR(max_iter_CR,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_CR', ABS(ierr) )
   ALLOCATE ( A_matrix_CL(max_iter_CL,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_CL', ABS(ierr) )

    ! Allocate B-type matrices
   ALLOCATE ( B_matrix_L(max_iter_L,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_L', ABS(ierr) )
   ALLOCATE ( B_matrix_R(max_iter_R,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_R', ABS(ierr) )
   ALLOCATE ( B_matrix_C((max_iter_CC - 1),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_C', ABS(ierr) )
   ALLOCATE ( B_matrix_CR((max_iter_CR - 1),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_CR', ABS(ierr) )
   ALLOCATE ( B_matrix_CL((max_iter_CL - 1),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_CL', ABS(ierr) )
     !
   ALLOCATE ( B_matrix_CR_R(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_CR_R', ABS(ierr) )
   ALLOCATE ( B_matrix_L_CL(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_L_CL', ABS(ierr) )
   ALLOCATE ( B_matrix_C_CR(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_C_CR', ABS(ierr) )
   ALLOCATE ( B_matrix_CL_C(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_CL_C', ABS(ierr) )


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
     IF ( ALLOCATED( A_matrix_CL  ) ) THEN
           DEALLOCATE ( A_matrix_CL, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_CL', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_CL  ) ) THEN
           DEALLOCATE ( B_matrix_CL, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_CL', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( B_matrix_L_CL  ) ) THEN
           DEALLOCATE ( B_matrix_L_CL, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_L_CL', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_CL_C  ) ) THEN
           DEALLOCATE ( B_matrix_CL_C, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_CL_C', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_C_CR  ) ) THEN
           DEALLOCATE ( B_matrix_C_CR, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_C_CR', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_CR_R  ) ) THEN
           DEALLOCATE ( B_matrix_CR_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_CR_R', ABS(ierr) )
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
       INTEGER   ::   iter, shift_iter


     !Beginning
    WRITE(stdout,"(2x,70('='))" )
    WRITE(stdout,"(/,  7x,'Building A & B matrices for recursion input')")
    WRITE(stdout,"(2x,70('='))" )
    WRITE(stdout,"()")
    WRITE(stdout,"()")

!!!!!!             L !

    A_matrix_L(:,:,:) = CZERO

    DO iter=1,max_iter_L
      ! corresponding iter in the R/L matrix 
      shift_iter = max_iter_L - iter + 1

      A_matrix_L(iter,:,:) = A_matrix(shift_iter,:,:) 
      !
      !
    ENDDO

!!!!!!             CL !

    A_matrix_CL(:,:,:) = CZERO

    DO iter=1,max_iter_CL
      ! corresponding iter in the R/L matrix 
      shift_iter = max_iter_L + max_iter_CL - iter + 1

      A_matrix_CL(iter,:,:) = A_matrix(shift_iter,:,:) 
      !
    ENDDO

!!!!!!!!!!!      C !

    A_matrix_C(:,:,:) = CZERO

    DO iter=1,max_iter_CC
         !
         shift_iter =  max_iter_L + max_iter_CL + iter 
         !
         A_matrix_C(iter,:,:) = A_matrix(shift_iter,:,:)
         !
    ENDDO

!!!!!!!!!!!      CR !

    A_matrix_CR(:,:,:) = CZERO

    DO iter=1,max_iter_CR
      ! corresponding iter in the R/L matrix 
      shift_iter = max_iter_L + max_iter_CL + max_iter_CC + iter 

      A_matrix_CR(iter,:,:) = A_matrix(shift_iter,:,:)

    ENDDO

!!!!!!!!!!!      R !
    A_matrix_R(:,:,:) = CZERO

    DO iter=1,max_iter_R
      ! corresponding iter in the R/L matrix 
      shift_iter = max_iter_L + max_iter_CL + max_iter_CC + max_iter_CR + iter 

      A_matrix_R(iter,:,:) = A_matrix(shift_iter,:,:)

    ENDDO



!!!!!!!!!!!      L !

    B_matrix_L(:,:,:) = CZERO

    DO iter=1,max_iter_L
      shift_iter = max_iter_L - iter + 1

      !
      DO icols=1, dim_subspace
         !
         DO irows=1, dim_subspace
               !
               B_matrix_L(iter,irows,icols) = CONJG(B_matrix(shift_iter,icols,irows) )
               !
         ENDDO
         !
      ENDDO
      !
    ENDDO

!!!!!!!!!!!      L - CL!

    B_matrix_L_CL(:,:,:) = CZERO

    shift_iter = max_iter_L + 1

    !
    B_matrix_L_CL(1,:,:) = B_matrix(shift_iter,:,:)

!!!!!!             CL !

    B_matrix_CL(:,:,:) = CZERO

    DO iter=1,(max_iter_CL-1)
      shift_iter = max_iter_L + max_iter_CL - iter + 1

      DO icols=1, dim_subspace
         !
         DO irows=1, dim_subspace
                  !
                  B_matrix_CL(iter,irows,icols) = CONJG(B_matrix(shift_iter,icols,irows) )
                  !
         ENDDO
         !
      ENDDO
      !
    ENDDO

!!!!!!!!!!!       CL - C!

    B_matrix_CL_C(:,:,:) = CZERO

    shift_iter = max_iter_L + max_iter_CL + 1

    !
    B_matrix_CL_C(1,:,:) = B_matrix(shift_iter,:,:)

!!!!!!!!!!!        C!

    B_matrix_C(:,:,:) = CZERO

    DO iter=1,(max_iter_CC-1)
         shift_iter = max_iter_L + max_iter_CL + iter  + 1

         B_matrix_C(iter,:,:)= B_matrix(shift_iter,:,:)
         !
    ENDDO

!!!!!!!!!!!       C - CR!

    B_matrix_C_CR(:,:,:) = CZERO

    shift_iter = max_iter_L + max_iter_CL + max_iter_CC  + 1

    B_matrix_C_CR(1,:,:)= B_matrix(shift_iter,:,:)

!!!!!!!!!!!       CR!

    B_matrix_CR(:,:,:) = CZERO

    DO iter=1,(max_iter_CR-1)
         shift_iter = max_iter_L + max_iter_CL + max_iter_CC + iter  + 1

         B_matrix_CR(iter,:,:)= B_matrix(shift_iter,:,:)
         !
    ENDDO

!!!!!!!!!!!       CR - R!

    B_matrix_CR_R(:,:,:) = CZERO

    shift_iter = max_iter_L + max_iter_CL + max_iter_CC + max_iter_CR + 1

    B_matrix_CR_R(1,:,:)= B_matrix(shift_iter,:,:)

!!!!!!!!!!!       R!

    B_matrix_R(:,:,:) = CZERO

    DO iter=1,max_iter_R
         shift_iter = max_iter_L + max_iter_CL + max_iter_CC + max_iter_CR + iter  + 1

         B_matrix_R(iter,:,:)= B_matrix(shift_iter,:,:)
         !
    ENDDO



    WRITE(stdout,"()")
    WRITE(stdout,"()")
    WRITE(stdout,"(2x,70('='))" )
    WRITE( stdout,"(/, 7x,'A & B matrices for recursion built ')")
    WRITE(stdout,"(2x,70('='))" )


   END SUBROUTINE build_matrix


 END MODULE matrix_module


