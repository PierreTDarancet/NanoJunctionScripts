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
                                     max_iter_term,   &
                                     max_iter_renorm

   USE in_matrix_module,              ONLY : A_matrix, B_matrix




   IMPLICIT NONE
   PRIVATE 
   SAVE

   ! variables for recursion
            ! For first iteration
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_renorm(:,:,:), B_matrix_renorm(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_term(:,:,:), B_matrix_term(:,:,:)
   !
   COMPLEX(dbl), ALLOCATABLE ::  B_matrix_lead_bulk(:,:,:)
   !

   LOGICAL :: alloc

   PUBLIC :: A_matrix_renorm, B_matrix_renorm
   PUBLIC :: A_matrix_term, B_matrix_term
   PUBLIC :: B_matrix_lead_bulk

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
   ALLOCATE ( A_matrix_renorm(max_iter_renorm,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_renorm', ABS(ierr) )
   ALLOCATE ( A_matrix_term(max_iter_term,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_term', ABS(ierr) )

    ! Allocate B-type matrices
   ALLOCATE ( B_matrix_renorm((max_iter_renorm-1),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_renorm', ABS(ierr) )
   ALLOCATE ( B_matrix_term(max_iter_term,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_term', ABS(ierr) )
   ALLOCATE ( B_matrix_lead_bulk(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_lead_bulk', ABS(ierr) )

   ! set alloc to true
    alloc=.TRUE.

END SUBROUTINE matrix_allocate

!**********************************************************
   SUBROUTINE matrix_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(17)      :: subname="matrix_deallocate"
      INTEGER :: ierr

     IF ( ALLOCATED( A_matrix_renorm  ) ) THEN
           DEALLOCATE ( A_matrix_renorm, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_renorm', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( A_matrix_term  ) ) THEN
           DEALLOCATE ( A_matrix_term, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_term', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_renorm  ) ) THEN
           DEALLOCATE ( B_matrix_renorm, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_renorm', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_term  ) ) THEN
           DEALLOCATE ( B_matrix_term, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_term', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_lead_bulk  ) ) THEN
           DEALLOCATE ( B_matrix_lead_bulk, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_lead_bulk', ABS(ierr) )
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

!!!!!!            renorm !

    A_matrix_renorm(:,:,:) = CZERO

    DO iter=1,max_iter_renorm
      !
      A_matrix_renorm(iter,:,:) = A_matrix(iter,:,:) 
      !
    ENDDO

!!!!!!             term !

    A_matrix_term(:,:,:) = CZERO

    DO iter=1,max_iter_term
      ! corresponding iter in the R/L matrix 
      shift_iter = max_iter_renorm + iter

      A_matrix_term(iter,:,:) = A_matrix(shift_iter,:,:) 
      !
    ENDDO

!!!!!!!!!!!       renorm !

    B_matrix_renorm(:,:,:) = CZERO

    DO iter=1,(max_iter_renorm-1)
      !
      B_matrix_renorm(iter,:,:) = B_matrix(iter,:,:) 
      !
    ENDDO

!!!!!!             lead - bulk !

    B_matrix_lead_bulk(:,:,:) = CZERO
      !
      shift_iter = max_iter_renorm
      !
      B_matrix_lead_bulk(1,:,:) = B_matrix(shift_iter,:,:)
      !

!!!!!!!!!!!       term!

    B_matrix_term(:,:,:) = CZERO

    DO iter=1,max_iter_term
      !
      shift_iter = max_iter_renorm  + iter
      !
      B_matrix_term(iter,:,:)= B_matrix(shift_iter,:,:)
      !
    ENDDO

!!!!!!!!!!!

    WRITE(stdout,"()")
    WRITE(stdout,"()")
    WRITE(stdout,"(2x,70('='))" )
    WRITE( stdout,"(/, 7x,'A & B matrices for recursion built ')")
    WRITE(stdout,"(2x,70('='))" )


   END SUBROUTINE build_matrix


 END MODULE matrix_module


