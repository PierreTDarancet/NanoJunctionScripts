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
   SUBROUTINE  green_tridiag( g_tilde_11, g_tilde_1N, g_tilde_N1, ene, A_matrix, B_matrix, sgm, max_iter, dim_sub) 
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE, EPS_m6, EPS_m4, EPS_m2, CI
   USE util_module,          ONLY : mat_mul, mat_sv
   USE timing_module, ONLY : timing
   IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! DEBUG : somme uniquement la partie  diagonale de A et la partie reelle de B 
!!!!!!!!!! valable un iquement dans le cas considere et pas general pour un sou
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Le but est d'obtenir l'element de matrice P1 G_tilde (z) PN
! les 1 et N sont censés correspondre aux indices des matrices A et B
! Sgm est la self energy des contacts L et R et est donc projetee sur PN
! ATTENTION les g_diag et g_off_diag NE SONT dans l'ensemble PAS les elements de matrice
! de la fonction Green TOTALE MAIS des elements de matrices de la fonction de Green 
! d'une restriction du Hamiltonien


   !
   ! Input variables
   !
   INTEGER,      INTENT(in)    :: dim_sub
   INTEGER,      INTENT(in)    :: max_iter
   COMPLEX(dbl), INTENT(in)    :: A_matrix(max_iter,dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(in)    :: B_matrix((max_iter-1),dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(in)    :: sgm(dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(in)    :: ene
   !
   ! Output variables
   !
   COMPLEX(dbl), INTENT(out)   :: g_tilde_1N(dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(out)   :: g_tilde_N1(dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(out)   :: g_tilde_11(dim_sub,dim_sub)

   !
   ! local variables
   !
   INTEGER                   :: i_sub, iter, shift_iter, j_sub
   COMPLEX(dbl)              :: i_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: h_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: g_diag(max_iter,dim_sub, dim_sub)
   COMPLEX(dbl)              :: g_offdiag((max_iter-1),dim_sub, dim_sub)
   COMPLEX(dbl)              :: g_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: sgm_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: bconj_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: delta
   CHARACTER(13)             :: subname='green_tridiag'

   ! Ajout debug
   REAL(dbl)                 :: A_matrix_debug(max_iter,dim_sub)
   COMPLEX(dbl)              :: B_matrix_debug((max_iter-1),dim_sub,dim_sub)
   COMPLEX(dbl)              :: h_work_herm(dim_sub,dim_sub)
   COMPLEX(dbl)              :: h_work_antiherm(dim_sub,dim_sub)
!
!----------------------------------------
! main Body
!----------------------------------------
!
   CALL timing('green_tridiag',OPR='start')

   ! Init
   g_work(:,:) = CZERO
   g_diag(:,:,:) = CZERO
   h_work(:,:) = CZERO
   i_work(:,:) = CZERO
   ! diag element

PRINT*, "B_matrix"
   DO iter = 1, (max_iter-1)
      !
      DO i_sub = 1, dim_sub
         DO j_sub = 1, dim_sub
             !
             PRINT*,REAL(B_matrix(iter,i_sub,j_sub)), AIMAG(B_matrix(iter,i_sub,j_sub))
         ENDDO
      ENDDO
   ENDDO

  PRINT*, "/B_matrix"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUG
   A_matrix_debug(:,:)= ZERO
   B_matrix_debug(:,:,:) = CZERO

   DO iter = 1, max_iter
      !
      DO i_sub = 1, dim_sub
         !
         A_matrix_debug(iter,i_sub) = REAL (A_matrix(iter,i_sub,i_sub))
         !
         IF ( ABS( A_matrix_debug(iter,i_sub)) < EPS_m6 )  THEN
             A_matrix_debug(iter,i_sub) = ZERO
         ENDIF
         !
      ENDDO
   ENDDO

   DO iter = 1, (max_iter-1)
      !
      DO j_sub = 1, dim_sub
         DO i_sub = 1, dim_sub
            !
            B_matrix_debug(iter,i_sub,j_sub) = REAL(B_matrix(iter,i_sub,j_sub))
            !
            IF ( ABS( B_matrix_debug(iter,i_sub, j_sub)) < EPS_m6 )  THEN
               B_matrix_debug(iter,i_sub,j_sub) = CZERO
            ENDIF
   
         ENDDO
      ENDDO
   ENDDO
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/DEBUG

      ! First element
      g_diag(max_iter,:,:) = CZERO
      DO i_sub = 1, dim_sub
         !
         g_diag(max_iter,i_sub,i_sub)= CONE
         i_work(i_sub,i_sub) = ene - A_matrix_debug(max_iter,i_sub)
         !
      ENDDO

      h_work(:,:) = i_work(:,:)  - sgm(:,:)


      CALL mat_sv(dim_sub, dim_sub, h_work, g_diag(max_iter,:,:))
      !
!       PRINT*, "Sgm Bulk"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(sgm(i_sub,j_sub)), AIMAG(sgm(i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/Sgm Bulk"
      h_work_herm(:,:) = 0.5000 * (h_work(:,:) + CONJG(TRANSPOSE(h_work)))

      h_work_antiherm(:,:) = 0.5000 * (h_work(:,:) - CONJG(TRANSPOSE(h_work)))


      PRINT*, "Sgm herm"
      DO i_sub = 1, dim_sub
         !
         DO j_sub = 1, dim_sub
            !
            PRINT*, REAL(h_work_herm(i_sub,j_sub)), AIMAG(h_work_herm(i_sub,j_sub))
            !
         ENDDO
      ENDDO
      PRINT*, "/Sgm herm"
      PRINT*, "Sgm antiherm"
      DO i_sub = 1, dim_sub
         !
         DO j_sub = 1, dim_sub
            !
            PRINT*, REAL(h_work_antiherm(i_sub,j_sub)), AIMAG(h_work_antiherm(i_sub,j_sub))
            !
         ENDDO
      ENDDO
      PRINT*, "/Sgm antiherm"
! 
!       PRINT*, "h_work"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(h_work(i_sub,j_sub)), AIMAG(h_work(i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/h_work"
! 
!       PRINT*, "i_work"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(i_work(i_sub,j_sub)), AIMAG(i_work(i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/i_work"
! 
!       PRINT*, "g_diag22"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(g_diag(max_iter,i_sub,j_sub)), AIMAG(g_diag(max_iter,i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/g_diag22"




      DO iter = 1, (max_iter-1)
         shift_iter = max_iter - iter
         !
         h_work(:,:) = CZERO
         i_work(:,:) = CZERO
         sgm_work(:,:) = CZERO
         g_work(:,:) = CZERO
         bconj_work(:,:) = CZERO
         g_diag(shift_iter,:,:) = CZERO
         !
         bconj_work(:,:) = TRANSPOSE( B_matrix_debug(shift_iter,:,:) )
         CALL mat_mul(g_work, bconj_work, 'N', g_diag((shift_iter+1),:,:), 'N', dim_sub, dim_sub, dim_sub)
         CALL mat_mul(sgm_work, g_work, 'N', B_matrix_debug(shift_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
         !
         !
 
         DO i_sub = 1, dim_sub
            !
            g_diag(shift_iter,i_sub,i_sub)= CONE
            i_work(i_sub,i_sub) = ene - A_matrix_debug(shift_iter,i_sub)
            !
         ENDDO

         !
         h_work(:,:) = i_work(:,:)  - sgm_work(:,:)
         !
!          DO i_sub = 1, dim_sub
!             !
!             DO j_sub = 1, dim_sub
!                !
!                PRINT*, REAL(sgm_work(i_sub,j_sub)), AIMAG(sgm_work(i_sub,j_sub))
!                !
!             ENDDO
!          ENDDO
!      PRINT*, "Sgm"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(sgm_work(i_sub,j_sub)), AIMAG(sgm_work(i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/Sgm Bulk"
! 
!       PRINT*, "h_work"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(h_work(i_sub,j_sub)), AIMAG(h_work(i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/h_work"
! 
!       PRINT*, "i_work"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(i_work(i_sub,j_sub)), AIMAG(i_work(i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/i_work"
! 
      h_work_herm(:,:) = 0.5000 * (h_work(:,:) + CONJG(TRANSPOSE(h_work)))

      h_work_antiherm(:,:) = 0.5000 * (h_work(:,:) - CONJG(TRANSPOSE(h_work)))


      PRINT*, "Sgm herm"
      DO i_sub = 1, dim_sub
         !
         DO j_sub = 1, dim_sub
            !
            PRINT*, REAL(h_work_herm(i_sub,j_sub)), AIMAG(h_work_herm(i_sub,j_sub))
            !
         ENDDO
      ENDDO
      PRINT*, "/Sgm herm"
      PRINT*, "Sgm antiherm"
      DO i_sub = 1, dim_sub
         !
         DO j_sub = 1, dim_sub
            !
            PRINT*, REAL(h_work_antiherm(i_sub,j_sub)), AIMAG(h_work_antiherm(i_sub,j_sub))
            !
         ENDDO
      ENDDO
      PRINT*, "/Sgm antiherm"

         !
         CALL mat_sv(dim_sub, dim_sub, h_work, g_diag(shift_iter,:,:))
         !



!       PRINT*, "g_diag"
!       DO i_sub = 1, dim_sub
!          !
!          DO j_sub = 1, dim_sub
!             !
!             PRINT*, REAL(g_diag(shift_iter,i_sub,j_sub)), AIMAG(g_diag(shift_iter,i_sub,j_sub))
!             !
!          ENDDO
!       ENDDO
!       PRINT*, "/g_diag"
       ENDDO

   ! off diag element

      ! First element
      g_offdiag((max_iter-1),:,:) =  CZERO
      g_work(:,:) = CZERO
      bconj_work(:,:) = TRANSPOSE( B_matrix_debug((max_iter-1),:,:) )
      !
      CALL mat_mul(g_work, bconj_work, 'N', g_diag(max_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
      CALL mat_mul(g_offdiag((max_iter-1),:,:), g_diag((max_iter-1),:,:) , 'N', g_work, 'N', dim_sub, dim_sub, dim_sub)
      !
      !
      DO iter = 1, (max_iter-2)
         shift_iter = max_iter - 1 - iter
         g_offdiag(shift_iter,:,:) =  CZERO
         g_work(:,:) = CZERO
         bconj_work(:,:) = TRANSPOSE( B_matrix_debug(shift_iter,:,:) )
        CALL mat_mul(g_work, bconj_work, 'N', g_offdiag((shift_iter+1),:,:), 'N', dim_sub, dim_sub, dim_sub)
        CALL mat_mul(g_offdiag(shift_iter,:,:), g_diag(shift_iter,:,:) , 'N', g_work, 'N', dim_sub, dim_sub, dim_sub)
      ENDDO


   ! Final step
   !
   g_tilde_1N(:,:) = g_offdiag(1,:,:)
   !

   ! off diag element

      ! First element
      g_offdiag(:,:,:) =  CZERO
      g_work(:,:) = CZERO
      !bconj_work(:,:) = CONJG( TRANSPOSE( B_matrix((max_iter-1),:,:) ) )
      !
      CALL mat_mul(g_work, g_diag(max_iter,:,:), 'N', B_matrix_debug((max_iter-1),:,:), 'N', dim_sub, dim_sub, dim_sub)
      CALL mat_mul(g_offdiag((max_iter-1),:,:), g_work , 'N', g_diag((max_iter-1),:,:), 'N', dim_sub, dim_sub, dim_sub)
      !
      !
      DO iter = 1, (max_iter-2)
         shift_iter = max_iter - 1 - iter
         g_offdiag(shift_iter,:,:) =  CZERO
         g_work(:,:) = CZERO
         !bconj_work(:,:) = CONJG( TRANSPOSE( B_matrix(shift_iter,:,:) ) )
        CALL mat_mul(g_work, g_offdiag((shift_iter+1),:,:), 'N', B_matrix_debug(shift_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
        CALL mat_mul(g_offdiag(shift_iter,:,:), g_work  , 'N', g_diag(shift_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
      ENDDO


   ! Final step
   !
   g_tilde_N1(:,:) = g_offdiag(1,:,:)
   !
   g_tilde_11(:,:) = g_diag(1,:,:)
   !
   CALL timing('green_tridiag',OPR='stop')
   !

END SUBROUTINE green_tridiag
