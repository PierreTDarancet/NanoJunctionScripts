!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE distance_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE T_control_module,    ONLY : dim_subspace, &
                                            max_iter_R,   &
                                            max_iter_L
   USE matrix_module,                ONLY : A_matrix_L, A_matrix_R, &
                                            B_matrix_L, B_matrix_R, matrix_alloc => alloc
   USE kinds
   USE constants,            ONLY : CZERO, CONE, ZERO, EPS_m11


   IMPLICIT NONE
   PRIVATE 
   SAVE

! Variables
!
    INTEGER              :: nb_max_first

    INTEGER, ALLOCATABLE :: nb_first_L(:)
!
    INTEGER, ALLOCATABLE :: id_first_L(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_first_L(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_onsite_L(:)
!
    INTEGER, ALLOCATABLE :: nb_first_R(:)
!
    INTEGER, ALLOCATABLE :: id_first_R(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_first_R(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_onsite_R(:)
!

    LOGICAL :: distance_alloc

! Status

!
   PUBLIC :: nb_max_first
!
   PUBLIC :: nb_first_L
!
   PUBLIC :: id_first_L
!
   PUBLIC :: ene_first_L
!
   PUBLIC :: ene_onsite_L
!
   PUBLIC :: nb_first_R
!
   PUBLIC :: id_first_R
!
   PUBLIC :: ene_first_R
!
   PUBLIC :: ene_onsite_R
!
!
   PUBLIC :: distance_alloc



! functions


  PUBLIC :: distance_allocate
  PUBLIC :: distance_deallocate
  PUBLIC :: init_metric



CONTAINS


!*******************************************************************
   SUBROUTINE distance_allocate
   !*******************************************************************
      CHARACTER(17)      :: subname="distance_allocate"
      INTEGER :: ierr, dim_recursion_L, dim_recursion_R

   !
   ! allocate basic quantities
   !
   nb_max_first =  (dim_subspace *3) - 1
   dim_recursion_L = dim_subspace * max_iter_L
   dim_recursion_R = dim_subspace * max_iter_R
   IF  (distance_alloc) CALL errore(subname,'distance already allocated',1)

   !
   ALLOCATE(  nb_first_L(dim_recursion_L),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating nb_first_L',ABS(ierr))
   !
   ALLOCATE( id_first_L(nb_max_first, dim_recursion_L),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating id_first_L',ABS(ierr))
   !
   ALLOCATE( ene_first_L(nb_max_first, dim_recursion_L),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_first_L',ABS(ierr))
   !
   ALLOCATE( ene_onsite_L(dim_recursion_L),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_onsite_L',ABS(ierr))
   !
   !
   ALLOCATE(  nb_first_R(dim_recursion_R),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating nb_first_R',ABS(ierr))
   !
   ALLOCATE( id_first_R(nb_max_first, dim_recursion_R),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating id_first_R',ABS(ierr))
   !
   ALLOCATE( ene_first_R(nb_max_first, dim_recursion_R),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_first_R',ABS(ierr))
   !
   ALLOCATE( ene_onsite_R(dim_recursion_R),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_onsite_R',ABS(ierr))
   !


    distance_alloc=.TRUE.



END SUBROUTINE distance_allocate



!**********************************************************
   SUBROUTINE distance_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(19)      :: subname="distance_deallocate"
      INTEGER :: ierr

    IF (.NOT. distance_alloc) CALL errore(subname,'distance not allocated',1)

     IF ( ALLOCATED( nb_first_L ) ) THEN
           DEALLOCATE ( nb_first_L, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating nb_first_L', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( id_first_L ) ) THEN
           DEALLOCATE ( id_first_L , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating id_first_L', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_first_L ) ) THEN
           DEALLOCATE ( ene_first_L , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_first_L', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_onsite_L ) ) THEN
           DEALLOCATE ( ene_onsite_L, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_onsite_L', ABS(ierr) )
      ENDIF
     !
     IF ( ALLOCATED( nb_first_R ) ) THEN
           DEALLOCATE ( nb_first_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating nb_first_R', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( id_first_R ) ) THEN
           DEALLOCATE ( id_first_R , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating id_first_R', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_first_R ) ) THEN
           DEALLOCATE ( ene_first_R , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_first_R', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_onsite_R ) ) THEN
           DEALLOCATE ( ene_onsite_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_onsite_R', ABS(ierr) )
      ENDIF
     !

     distance_alloc=.FALSE.

   END SUBROUTINE distance_deallocate


!*******************************************************************
   SUBROUTINE init_metric()
   !*******************************************************************
      CHARACTER(11)      :: subname="init_metric"
      INTEGER :: ierr
      INTEGER :: i_rec, i_sub, iter, i_ind, i_ind_shift
      INTEGER :: dim_recursion_L, dim_recursion_R

    !
    dim_recursion_L = dim_subspace * max_iter_L
    !
    dim_recursion_R = dim_subspace * max_iter_R


    IF (.NOT. matrix_alloc) CALL errore(subname,'matrices not allocated',1)

    !
    IF (.NOT. distance_alloc) CALL errore(subname,'distance not allocated',1)
    !init
    nb_first_L(:)     =  0
    !
    id_first_L(:,:)   =  0
    !
    ene_first_L(:,:)  =  CZERO
    !
    i_sub = 1
    iter  =  1
      DO i_rec=1, dim_recursion_L
          !
          IF (i_sub == dim_subspace + 1) THEN
               i_sub = 1
               iter  = iter + 1
          ENDIF
          !

          IF ( iter ==  1 ) THEN
              nb_first_L(i_rec) = (dim_subspace * 2) - 1
              DO i_ind = 1, dim_subspace
                 i_ind_shift = i_ind + dim_subspace - 1
                 id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + dim_subspace + i_ind
                 !Attention
                 ene_first_L(i_ind_shift,i_rec) = CONJG(B_matrix_L(iter, i_ind, i_sub))
                 !ene_first_L(i_ind_shift,i_rec) = CONJG(B_matrix_L(iter, i_sub, i_ind))
                 !
              ENDDO
              DO i_ind = 1, dim_subspace
                 IF (i_ind < i_sub) THEN
                     i_ind_shift = i_ind 
                     id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_L(i_ind_shift,i_rec) = A_matrix_L(iter, i_sub, i_ind)
                 ELSE IF (i_ind > i_sub) THEN
                     i_ind_shift = i_ind - 1
                     id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_L(i_ind_shift,i_rec) = A_matrix_L(iter, i_sub, i_ind)
                 ELSE
                     ene_onsite_L(i_rec) =  A_matrix_L(iter, i_sub, i_sub)
                 ENDIF
              ENDDO
          ELSE IF ( iter == max_iter_L) THEN
              nb_first_L(i_rec) = (dim_subspace * 2) - 1
              DO i_ind = 1, dim_subspace
                 id_first_L(i_ind,i_rec) = i_rec - i_sub - dim_subspace + i_ind
                 ene_first_L(i_ind,i_rec) = B_matrix_L(iter-1, i_sub, i_ind)
              ENDDO
              DO i_ind = 1, dim_subspace
                 IF (i_ind < i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace
                     id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_L(i_ind_shift,i_rec) = A_matrix_L(iter, i_sub, i_ind)
                 ELSE IF (i_ind > i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace - 1
                     id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_L(i_ind_shift,i_rec) = A_matrix_L(iter, i_sub, i_ind)
                 ELSE
                     ene_onsite_L(i_rec) =  A_matrix_L(iter, i_sub, i_sub)
                 ENDIF
              ENDDO
          ELSE
              nb_first_L(i_rec) = (dim_subspace * 3) - 1
              DO i_ind = 1, dim_subspace
                 i_ind_shift = i_ind + 2*dim_subspace - 1
                 id_first_L(i_ind,i_rec) = i_rec - i_sub - dim_subspace + i_ind
                 ene_first_L(i_ind,i_rec) = B_matrix_L(iter-1, i_sub, i_ind)
                 id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + dim_subspace + i_ind
                 ene_first_L(i_ind_shift,i_rec) = CONJG(B_matrix_L(iter, i_ind, i_sub))
              ENDDO
              DO i_ind = 1, dim_subspace
                 IF (i_ind < i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace
                     id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_L(i_ind_shift,i_rec) = A_matrix_L(iter, i_sub, i_ind)
                 ELSE IF (i_ind > i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace - 1
                     id_first_L(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_L(i_ind_shift,i_rec) = A_matrix_L(iter, i_sub, i_ind)
                 ELSE
                     ene_onsite_L(i_rec) =  A_matrix_L(iter, i_sub, i_sub)
                 ENDIF
              ENDDO
          ENDIF
          !
          i_sub = i_sub + 1
          !
      ENDDO

    nb_first_R(:)     =  0
    !
    id_first_R(:,:)   =  0
    !
    ene_first_R(:,:)  =  CZERO
    !
    i_sub = 1
    iter  =  1
      DO i_rec=1, dim_recursion_R
          !
          IF (i_sub == dim_subspace + 1) THEN
               i_sub = 1
               iter  = iter + 1
          ENDIF
          !

          IF ( iter ==  1 ) THEN
              nb_first_R(i_rec) = (dim_subspace * 2) - 1
              DO i_ind = 1, dim_subspace
                 i_ind_shift = i_ind + dim_subspace - 1
                 id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + dim_subspace + i_ind
                 ene_first_R(i_ind_shift,i_rec) = CONJG(B_matrix_R(iter, i_ind, i_sub))
              ENDDO
              DO i_ind = 1, dim_subspace
                 IF (i_ind < i_sub) THEN
                     i_ind_shift = i_ind 
                     id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_R(i_ind_shift,i_rec) = A_matrix_R(iter, i_sub, i_ind)
                 ELSE IF (i_ind > i_sub) THEN
                     i_ind_shift = i_ind - 1
                     id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_R(i_ind_shift,i_rec) = A_matrix_R(iter, i_sub, i_ind)
                 ELSE
                     ene_onsite_R(i_rec) =  A_matrix_R(iter, i_sub, i_sub)
                 ENDIF
              ENDDO
          ELSE IF ( iter == max_iter_R) THEN
              nb_first_R(i_rec) = (dim_subspace * 2) - 1
              DO i_ind = 1, dim_subspace
                 id_first_R(i_ind,i_rec) = i_rec - i_sub - dim_subspace + i_ind
                 ene_first_R(i_ind,i_rec) = B_matrix_R(iter-1, i_sub, i_ind)
              ENDDO
              DO i_ind = 1, dim_subspace
                 IF (i_ind < i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace
                     id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_R(i_ind_shift,i_rec) = A_matrix_R(iter, i_sub, i_ind)
                 ELSE IF (i_ind > i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace - 1
                     id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_R(i_ind_shift,i_rec) = A_matrix_R(iter, i_sub, i_ind)
                 ELSE
                     ene_onsite_R(i_rec) =  A_matrix_R(iter, i_sub, i_sub)
                 ENDIF
              ENDDO
          ELSE
              nb_first_R(i_rec) = (dim_subspace * 3) - 1
              DO i_ind = 1, dim_subspace
                 i_ind_shift = i_ind + 2*dim_subspace - 1
                 id_first_R(i_ind,i_rec) = i_rec - i_sub - dim_subspace + i_ind
                 ene_first_R(i_ind,i_rec) = B_matrix_R(iter-1, i_sub, i_ind)
                 id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + dim_subspace + i_ind
                 ene_first_R(i_ind_shift,i_rec) = CONJG(B_matrix_R(iter, i_ind, i_sub))
              ENDDO
              DO i_ind = 1, dim_subspace
                 IF (i_ind < i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace
                     id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_R(i_ind_shift,i_rec) = A_matrix_R(iter, i_sub, i_ind)
                 ELSE IF (i_ind > i_sub) THEN
                     i_ind_shift = i_ind + dim_subspace - 1
                     id_first_R(i_ind_shift,i_rec) = i_rec - i_sub + i_ind
                     ene_first_R(i_ind_shift,i_rec) = A_matrix_R(iter, i_sub, i_ind)
                 ELSE
                     ene_onsite_R(i_rec) =  A_matrix_R(iter, i_sub, i_sub)
                 ENDIF
              ENDDO
          ENDIF
          !
          i_sub = i_sub + 1
          !
      ENDDO



END SUBROUTINE init_metric


END MODULE distance_module


