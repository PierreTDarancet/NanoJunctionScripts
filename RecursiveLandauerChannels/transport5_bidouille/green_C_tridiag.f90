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
   SUBROUTINE  green_CC_tridiag( g_CC_1N, g_CC_N1, ene, A_matrix_in, B_matrix, sgm_L, sgm_R, max_iter, dim_sub) 
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE, EPS_m6, EPS_m4, EPS_m2, CI
   USE util_module,          ONLY : mat_mul, mat_sv
   USE timing_module, ONLY : timing
   IMPLICIT NONE

  
! Le but est d'obtenir l'element de matrice P1 G_CC (z) PN
! les 1 et N sont censés correspondre aux indices des matrices A et B
! Sgm_L/R est la self energy des contacts L et R et est donc projetee sur P1/PN
! ATTENTION les g_diag et g_off_diag NE SONT dans l'ensemble PAS les elements de matrice
! de la fonction Green TOTALE MAIS des elements de matrices de la fonction de Green 
! d'une restriction du Hamiltonien


   !
   ! Input variables
   !
   INTEGER,      INTENT(in)    :: dim_sub
   INTEGER,      INTENT(in)    :: max_iter
   COMPLEX(dbl), INTENT(in)    :: A_matrix_in(max_iter,dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(in)    :: B_matrix((max_iter-1),dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(in)    :: sgm_L(dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(in)    :: sgm_R(dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(in)    :: ene
   !
   ! Output variables
   !
   COMPLEX(dbl), INTENT(out)   :: g_CC_1N(dim_sub,dim_sub)
   COMPLEX(dbl), INTENT(out)   :: g_CC_N1(dim_sub,dim_sub)

   !
   ! local variables
   !
   INTEGER                   :: i_sub, iter, shift_iter
   COMPLEX(dbl)              :: i_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: h_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: g_diag(max_iter,dim_sub, dim_sub)
   COMPLEX(dbl)              :: g_offdiag((max_iter-1),dim_sub, dim_sub)
   COMPLEX(dbl)              :: g_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: sgm_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: bconj_work(dim_sub,dim_sub)
   COMPLEX(dbl)              :: delta
   CHARACTER(16)             :: subname='green_CC_tridiag'
   COMPLEX(dbl)              :: A_matrix(max_iter,dim_sub,dim_sub)
!
!----------------------------------------
! main Body
!----------------------------------------
!
   CALL timing('green_tridiag',OPR='start')
   A_matrix(:,:,:)= A_matrix_in(:,:,:)
   ! Init
   g_work(:,:) = CZERO
   g_diag(:,:,:) = CZERO
   h_work(:,:) = CZERO
   i_work(:,:) = CZERO
   sgm_work(:,:) = CZERO
   ! diag element

      ! Integration des self energies
      sgm_work(:,:) = A_matrix(1,:,:) + sgm_L(:,:)
      A_matrix(1,:,:) = sgm_work(:,:)
      sgm_work(:,:) = A_matrix(max_iter,:,:) + sgm_R(:,:)
      A_matrix(max_iter,:,:) = sgm_work(:,:)
      ! First element
      g_diag(max_iter,:,:) = CZERO
      DO i_sub = 1, dim_sub
         !
         g_diag(max_iter,i_sub,i_sub)= CONE
         i_work(i_sub,i_sub) = ene
         !
      ENDDO

      h_work(:,:) = i_work(:,:) - A_matrix(max_iter,:,:)

      CALL mat_sv(dim_sub, dim_sub, h_work, g_diag(max_iter,:,:))
      !

      DO iter = 1, (max_iter-1)
         shift_iter = max_iter - iter
         !
         bconj_work(:,:) = CONJG( TRANSPOSE( B_matrix(shift_iter,:,:) ) )
         CALL mat_mul(g_work, bconj_work, 'N', g_diag((shift_iter+1),:,:), 'N', dim_sub, dim_sub, dim_sub)
         CALL mat_mul(sgm_work, g_work, 'N', B_matrix(shift_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
         !
         h_work(:,:) = CZERO
         h_work(:,:) = i_work(:,:) - A_matrix(shift_iter,:,:) - sgm_work(:,:)
         !
         g_diag(shift_iter,:,:) = CZERO
         DO i_sub = 1, dim_sub
             !
             g_diag(shift_iter,i_sub,i_sub)= CONE
             !
         ENDDO
         !
         CALL mat_sv(dim_sub, dim_sub, h_work, g_diag(shift_iter,:,:))
         !
      ENDDO

   ! off diag element

      ! First element
      g_offdiag((max_iter-1),:,:) =  CZERO
      g_work(:,:) = CZERO
      bconj_work(:,:) = CONJG( TRANSPOSE( B_matrix((max_iter-1),:,:) ) )
      CALL mat_mul(g_work, bconj_work, 'N', g_diag(max_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
      CALL mat_mul(g_offdiag((max_iter-1),:,:), g_diag((max_iter-1),:,:) , 'N', g_work, 'N', dim_sub, dim_sub, dim_sub)
      !
      !
      DO iter = 1, (max_iter-2)
         shift_iter = max_iter - 1 - iter
         g_offdiag(shift_iter,:,:) =  CZERO
         g_work(:,:) = CZERO
         bconj_work(:,:) = CONJG( TRANSPOSE( B_matrix(shift_iter,:,:) ) )
        CALL mat_mul(g_work, bconj_work, 'N', g_offdiag((shift_iter+1),:,:), 'N', dim_sub, dim_sub, dim_sub)
        CALL mat_mul(g_offdiag(shift_iter,:,:), g_diag(shift_iter,:,:) , 'N', g_work, 'N', dim_sub, dim_sub, dim_sub)
      ENDDO


   ! Final step
   !
   g_CC_1N(:,:) = g_offdiag(1,:,:)
   !

   ! off diag element

      ! First element
      g_offdiag(:,:,:) =  CZERO
      g_work(:,:) = CZERO
      

      CALL mat_mul(g_work, g_diag(max_iter,:,:), 'N', B_matrix((max_iter-1),:,:), 'N', dim_sub, dim_sub, dim_sub)
      CALL mat_mul(g_offdiag((max_iter-1),:,:), g_work , 'N', g_diag((max_iter-1),:,:), 'N', dim_sub, dim_sub, dim_sub)
      !
      !
      DO iter = 1, (max_iter-2)
         shift_iter = max_iter - 1 - iter
         g_offdiag(shift_iter,:,:) =  CZERO
         g_work(:,:) = CZERO
         !bconj_work(:,:) = CONJG( TRANSPOSE( B_matrix(shift_iter,:,:) ) )
        CALL mat_mul(g_work, g_offdiag((shift_iter+1),:,:), 'N', B_matrix(shift_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
        CALL mat_mul(g_offdiag(shift_iter,:,:), g_work , 'N', g_diag(shift_iter,:,:), 'N', dim_sub, dim_sub, dim_sub)
      ENDDO


   g_CC_N1(:,:) = g_offdiag(1,:,:)


   CALL timing('green_tridiag',OPR='stop')
   !

END SUBROUTINE green_CC_tridiag
