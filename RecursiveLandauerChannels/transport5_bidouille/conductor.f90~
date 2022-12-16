!
!      Copyright (C) 2004 WanT Group
!      Copyright (C) 1999 Marco Buongiorno Nardelli
!      Copyright (C) 2007 Pierre Darancet Institut Néel
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!***********************************************
   PROGRAM conductor5
   !***********************************************
   USE constants,            ONLY : PI, ZERO, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1
   USE distance_module,      ONLY : nb_first_L, id_first_L, ene_first_L, ene_onsite_L, &
                                    nb_first_R, id_first_R, ene_first_R, ene_onsite_R, &
                                    distance_allocate, init_metric
   USE files_module,         ONLY : file_open, file_close
   USE in_matrix_module,     ONLY : in_matrix_allocate, read_matrix, in_build_matrix
   USE io_module,            ONLY : stdout, stdin, sgm_unit => aux_unit,   &
                                    dos_unit => aux1_unit, cond_unit => aux2_unit
   USE iotk_module
   USE kinds,                ONLY : dbl
   USE matrix_module,        ONLY : matrix_allocate, build_matrix, &
                                    A_matrix_C, B_matrix_C, A_matrix_CR, B_matrix_CR, &
                                    A_matrix_CL, B_matrix_CL, B_matrix_L_CL, B_matrix_CR_R, &
                                    B_matrix_CL_C, B_matrix_C_CR
   USE parameters,           ONLY : nstrx 
   USE T_control_module,     ONLY : nprint, &
                                    max_iter_L, max_iter_R, max_iter_CL, max_iter_CR, max_iter_CC, dim_subspace, &
                                    cut_chain
   USE T_egrid_module,       ONLY : egrid_init, ne, egrid, delta, delta_lead
   USE T_input_module,       ONLY : input_manager
   USE timing_module,        ONLY : timing, timing_overview, global_list, timing_upto_now
   USE util_module,          ONLY : mat_mul, mat_hdiag
   USE version_module,       ONLY : version_number

   IMPLICIT NONE

   !
   ! local variables
   !
   COMPLEX(dbl)     :: ene, ene_lead
   CHARACTER(nstrx) :: filename
   INTEGER          :: i, ie, ierr, ncount, i_L, i_R, iter, i_rows, i_cols, i_sub
   REAL(dbl)        :: real_iter
   !
   COMPLEX(dbl),    ALLOCATABLE :: work(:,:)
   COMPLEX(dbl),    ALLOCATABLE :: work2(:,:)
   COMPLEX(dbl),    ALLOCATABLE :: g_aux(:,:)
   COMPLEX(dbl),    ALLOCATABLE :: B_aux(:,:)
   !




   REAL(dbl),    ALLOCATABLE :: conduct(:,:)

   !
   !
   ! Variables associees au calcul de G_L et G_R
   !
   REAL(dbl),    ALLOCATABLE :: rec_var_L(:,:,:), rec_var_R(:,:,:)
   COMPLEX(dbl), ALLOCATABLE :: rec_chain_L(:,:,:,:), rec_chain_R(:,:,:,:)
   INTEGER,      ALLOCATABLE :: cut_iter_R(:,:), cut_iter_L(:,:)
   !

   !
   ! Fonctions de Green
   !
   COMPLEX(dbl), ALLOCATABLE :: P1_g_tilde_CR_P1(:,:), P1_g_tilde_CR_PN(:,:), PN_g_tilde_CR_P1(:,:)
   COMPLEX(dbl), ALLOCATABLE :: P1_g_tilde_CL_P1(:,:), P1_g_tilde_CL_PN(:,:), PN_g_tilde_CL_P1(:,:)
   COMPLEX(dbl), ALLOCATABLE :: gL(:,:)
   COMPLEX(dbl), ALLOCATABLE :: gR(:,:)
   COMPLEX(dbl), ALLOCATABLE :: P1_g_CC_PN(:,:), PN_g_CC_P1(:,:)


   !
   ! Variables associees au calcul des Self energies
   !
   COMPLEX(dbl), ALLOCATABLE :: sgm_L(:,:), sgm_R(:,:)
   COMPLEX(dbl), ALLOCATABLE :: sgm_CL(:,:), sgm_CR(:,:)
   !

   !
   ! Taux d'injection
   !
   COMPLEX(dbl), ALLOCATABLE :: gamma_L(:,:), gamma_R(:,:)
   COMPLEX(dbl), ALLOCATABLE :: gamma_tilde_L(:,:), gamma_tilde_R(:,:)
   !



   !
   ! Variables associees aux fichiers de sortie
   !
   REAL(dbl),    ALLOCATABLE :: sigma_R(:,:), sigma_L(:,:)
   REAL(dbl),    ALLOCATABLE :: gamma_L_sum(:,:), gamma_R_sum(:,:)
   REAL(dbl),    ALLOCATABLE :: gamma_tilde_L_sum(:,:), gamma_tilde_R_sum(:,:)
   REAL(dbl),    ALLOCATABLE :: dos_CL(:), dos_CR(:), dos_L(:,:), dos_R(:,:)
   !REAL(dbl),    ALLOCATABLE :: dos_CC(:)
   REAL(dbl), ALLOCATABLE :: sigma_CR(:,:), sigma_CL(:,:)
   !


   !!!!!!!!! DEBUG 
   COMPLEX(dbl), ALLOCATABLE :: test_mat(:,:), z(:,:)
   REAL(dbl), ALLOCATABLE :: eig(:)

!
!------------------------------
! main body
!------------------------------
!
   CALL startup(version_number,'conductor')
   !
   ! read input file
   !
   CALL input_manager()
   !
   ! init
   !
   CALL in_matrix_allocate()
   CALL read_matrix()
   CALL in_build_matrix()
   CALL matrix_allocate()
   CALL build_matrix()
   CALL distance_allocate()
   CALL init_metric()
   CALL egrid_init()

   CALL summary( stdout )
   !
   ! local variable allocations
   !

!       ALLOCATE ( dos_CC(dim_subspace,ne), STAT=ierr )
!          IF( ierr /=0 ) CALL errore('conductor','allocating dos_CC', ABS(ierr) )
   ALLOCATE ( dos_CL(ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_CL', ABS(ierr) )
   ALLOCATE ( dos_CR(ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_CR', ABS(ierr) )
   ALLOCATE ( dos_R(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_R', ABS(ierr) )
   ALLOCATE ( dos_L(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_L', ABS(ierr) )
      !
   ALLOCATE ( sigma_R(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sigma_R', ABS(ierr) )
   ALLOCATE ( sigma_L(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sigma_L', ABS(ierr) )
   ALLOCATE ( sigma_CR(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sigma_CR', ABS(ierr) )
   ALLOCATE ( sigma_CL(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sigma_CL', ABS(ierr) )
   ALLOCATE ( gamma_R_sum(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_R_sum', ABS(ierr) )
   ALLOCATE ( gamma_L_sum(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_L_sum', ABS(ierr) )
   ALLOCATE ( gamma_tilde_R_sum(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_tilde_R_sum', ABS(ierr) )
   ALLOCATE ( gamma_tilde_L_sum(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_tilde_L_sum', ABS(ierr) )
      !
   ALLOCATE ( conduct(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating conduct', ABS(ierr) )
   ALLOCATE ( work(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating work', ABS(ierr) )
   ALLOCATE ( work2(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating work2', ABS(ierr) )
   ALLOCATE ( g_aux(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating g_aux', ABS(ierr) )
   ALLOCATE ( B_aux(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating B_aux', ABS(ierr) )
      !
   ALLOCATE ( rec_chain_L(2, max_iter_L, dim_subspace, dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating rec_chain_L', ABS(ierr) )
   ALLOCATE ( rec_chain_R(2, max_iter_R, dim_subspace, dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating rec_chain_R', ABS(ierr) )
   ALLOCATE ( rec_var_L(max_iter_L, dim_subspace, dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating rec_var_L', ABS(ierr) )
   ALLOCATE ( rec_var_R(max_iter_R, dim_subspace, dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating rec_var_R', ABS(ierr) )
   ALLOCATE ( cut_iter_L(dim_subspace, dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating cut_iter_L', ABS(ierr) )
   ALLOCATE ( cut_iter_R(dim_subspace, dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating cut_iter_R', ABS(ierr) )
      !
   ALLOCATE ( sgm_L(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sgm_L', ABS(ierr) )
   ALLOCATE ( sgm_R(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sgm_R', ABS(ierr) )
   ALLOCATE ( sgm_CL(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sgm_CL', ABS(ierr) )
   ALLOCATE ( sgm_CR(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sgm_CR', ABS(ierr) )
      !
   ALLOCATE ( gamma_R(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_R', ABS(ierr) )
   ALLOCATE ( gamma_L(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_L', ABS(ierr) )
   ALLOCATE ( gamma_tilde_R(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_tilde_R', ABS(ierr) )
   ALLOCATE ( gamma_tilde_L(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_tilde_L', ABS(ierr) )
      !
   ALLOCATE ( gL(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gL', ABS(ierr) )
   ALLOCATE ( gR(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gR', ABS(ierr) )
   ALLOCATE ( P1_g_tilde_CR_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_tilde_CR_P1', ABS(ierr) )
   ALLOCATE ( P1_g_tilde_CR_PN(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_tilde_CR_PN', ABS(ierr) )
   ALLOCATE ( PN_g_tilde_CR_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating PN_g_tilde_CR_P1', ABS(ierr) )
   ALLOCATE ( P1_g_tilde_CL_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_tilde_CL_P1', ABS(ierr) )
   ALLOCATE ( P1_g_tilde_CL_PN(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_tilde_CL_PN', ABS(ierr) )
   ALLOCATE ( PN_g_tilde_CL_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating PN_g_tilde_CL_P1', ABS(ierr) )
   ALLOCATE ( P1_g_CC_PN(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_CC_PN', ABS(ierr) )
   ALLOCATE ( PN_g_CC_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating PN_g_CC_P1', ABS(ierr) )




    !!!!!!!!!!!!! DEBUG
   ALLOCATE ( test_mat(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating testmat', ABS(ierr) )
   test_mat(:,:) = CZERO
   ALLOCATE ( eig(dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating testmat', ABS(ierr) )
   ALLOCATE ( z(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating z', ABS(ierr) )
   z(:,:) = CZERO
    !!!!!!!!!!!!! DEBUG



!          ! Calcul des chaines de recursion
!          WRITE(stdout,"()")
!          WRITE(stdout,"()")
!          WRITE(stdout, "('Rec on h00_R...')")
!          CALL scalar_recursion(rec_chain_R(:,:,:,:), nb_first_R(:), id_first_R(:,:), ene_first_R(:,:), ene_onsite_R(:), (max_iter_R*dim_subspace), max_iter_R, dim_subspace, rec_var_R(:,:,:), cut_chain, cut_iter_R(:,:) )
!          WRITE(stdout,"()")
!          !
!          filename = 'rec_var_R.dat'
!          OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
!          DO i_cols=1,dim_subspace 
!             WRITE( dos_unit, "('icols = ',i5,' ')") i_cols
!             DO i_rows = 1, dim_subspace
!                   WRITE( dos_unit, "('irows = ',i5,' ')") i_rows
!                   DO iter = 1, max_iter_R
!                      real_iter = REAL(iter, dbl)
!                      WRITE ( dos_unit, '(6(f15.9))' ) real_iter, rec_var_R(iter, i_rows, i_cols), rec_chain_R(1,iter,i_rows, i_cols), rec_chain_R(2,iter,i_rows, i_cols)
!                   ENDDO
!             ENDDO
!          ENDDO
!          CLOSE( dos_unit )
!          !
!          !
!          WRITE(stdout, "('Rec on  h00_L...')")
!          CALL scalar_recursion(rec_chain_L(:,:,:,:), nb_first_L(:), id_first_L(:,:), ene_first_L(:,:), ene_onsite_L(:), (max_iter_L*dim_subspace), max_iter_L, dim_subspace, rec_var_L(:,:,:), cut_chain, cut_iter_L(:,:) )
!          WRITE(stdout,"()")
!          !
!          filename = 'rec_var_L.dat'
!          OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
!          DO i_cols=1,dim_subspace
!             WRITE( dos_unit, "('icols = ',i5,' ')") i_cols
!             DO i_rows = 1, dim_subspace
!                   WRITE( dos_unit, "('irows = ',i5,' ')") i_rows
!                   DO iter = 1, max_iter_L
!                      real_iter = REAL(iter, dbl)
!                      WRITE ( dos_unit, '(6(f15.9))' ) real_iter, rec_var_L(iter, i_rows, i_cols), rec_chain_L(1,iter,i_rows, i_cols), rec_chain_L(2,iter,i_rows, i_cols)
!                   ENDDO
!             ENDDO
!          ENDDO
!          CLOSE( dos_unit )
   !

   energy_loop: &
   DO ie = 1, ne
      ncount = ie
      !
      ! grids and misc
      !
      ene =  egrid(ie)  + delta * CI
      ene_lead = egrid(ie)  + delta_lead * CI
           IF ( MOD( ncount, nprint) == 0 .OR. ncount == 1 ) THEN
!                WRITE( stdout,"(2x,'Energy step = ',i5) ") ncount
                WRITE(stdout,"(2x, 'Computing E( ',i5,' ) = ', f9.5, ' eV' )") ncount, egrid(ie)
                CALL timing_upto_now(stdout)
           ENDIF

      conduct(:,ie) = ZERO
      ! 
      ! construct leads self-energies 
      ! 

      CALL green_tout_neuf( rec_chain_R, ene_lead, gR, dim_subspace, max_iter_R, cut_chain, cut_iter_R)
      !
      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_CR_R(1,:,:)) )
      CALL mat_mul(work, B_aux, 'N', gR, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(sgm_R, work, 'N', B_matrix_CR_R(1,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)

      ! 
      CALL green_tout_neuf( rec_chain_L, ene_lead, gL, dim_subspace, max_iter_L, cut_chain, cut_iter_L)
      !
      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_L_CL(1,:,:)) )
      CALL mat_mul(work, B_matrix_L_CL(1,:,:), 'N', gL, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(sgm_L, work, 'N', B_aux, 'N', dim_subspace, dim_subspace, dim_subspace) 

      !
      ! gamma_L and gamma_R
      !
      gamma_L(:,:) = CI * (  sgm_L(:,:) - CONJG( TRANSPOSE(sgm_L(:,:)) )   )
      gamma_R(:,:) = CI * (  sgm_R(:,:) - CONJG( TRANSPOSE(sgm_R(:,:)) )   )

      !!!debug
         !
          !
         eig(:) =  CZERO
         !
      CALL mat_hdiag(z(:,:) , eig(:), gamma_L(:,:), dim_subspace) 
      test_mat(:,ie)=eig(:)
      !!!debug


      ! Calcul des fonctions G_tilde_CR et G_tilde CL
      CALL green_tridiag(P1_g_tilde_CR_P1(:,:), P1_g_tilde_CR_PN(:,:), PN_g_tilde_CR_P1(:,:), ene, A_matrix_CR(:,:,:), B_matrix_CR(:,:,:), sgm_R(:,:), max_iter_CR, dim_subspace) 
      CALL green_tridiag(P1_g_tilde_CL_P1(:,:), P1_g_tilde_CL_PN(:,:), PN_g_tilde_CL_P1(:,:), ene, A_matrix_CL(:,:,:), B_matrix_CL(:,:,:), sgm_L(:,:), max_iter_CL, dim_subspace) 

      ! Calcul des Sigma CR et des Sigma_CL
      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_C_CR(1,:,:)) )
      CALL mat_mul(work, B_aux, 'N', P1_g_tilde_CR_P1, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(sgm_CR, work, 'N', B_matrix_C_CR(1,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)
      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_CL_C(1,:,:)) )
      CALL mat_mul(work, B_matrix_CL_C(1,:,:), 'N', P1_g_tilde_CL_P1, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(sgm_CL, work, 'N', B_aux(:,:), 'N', dim_subspace, dim_subspace, dim_subspace)

      !
      ! gamma_tilde_L and gamma_tilde_R
      !
      ! gamma_tilde_R = H_C_CR * G^r_CR (z) * Gamma_R * G^a_CR (z) * H_CR_C
      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_C_CR(1,:,:)) )
      CALL mat_mul(work, B_aux, 'N', P1_g_tilde_CR_PN, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(work2, work, 'N', gamma_R, 'N', dim_subspace, dim_subspace, dim_subspace)
      g_aux(:,:)=  CONJG( TRANSPOSE(P1_g_tilde_CR_PN(:,:)))
      CALL mat_mul(work, work2, 'N', g_aux, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(gamma_tilde_R, work, 'N',B_matrix_C_CR(1,:,:) , 'N', dim_subspace, dim_subspace, dim_subspace)
      ! gamma_tilde_L = H_C_CL *P1 G^a_CL (z) * Gamma_L * G^r_CL (z) * H_CL_C
      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_CL_C(1,:,:)) )
      g_aux(:,:) =  CONJG( TRANSPOSE(PN_g_tilde_CL_P1(:,:)))
      CALL mat_mul(work, B_matrix_CL_C(1,:,:), 'N', g_aux, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(work2, work, 'N', gamma_L, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(work, work2, 'N', PN_g_tilde_CL_P1, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(gamma_tilde_L, work, 'N',B_aux , 'N', dim_subspace, dim_subspace, dim_subspace)
      !

      !
      ! Calcul de G_CC
      !
      CALL green_CC_tridiag( P1_g_CC_PN(:,:), PN_g_CC_P1(:,:), ene, A_matrix_C(:,:,:), B_matrix_C(:,:,:), sgm_CL(:,:), sgm_CR(:,:), max_iter_CC, dim_subspace) 


      !
      CALL transmittance_min(dim_subspace, gamma_tilde_L, gamma_tilde_R, P1_g_CC_PN(:,:), conduct(:,ie) )
      !


      ! Fichiers de sortie

      !
      ! Compute density of states for the conductor layer
      !
      !
      dos_L(:,ie)=ZERO
      DO i = 1, dim_subspace
         dos_L(i,ie) =  -  AIMAG( gL(i,i) ) / PI
      ENDDO
      !
      dos_R(:,ie)=ZERO
      DO i = 1, dim_subspace
         dos_R(i,ie) =  -  AIMAG ( gR(i,i) ) / PI
      ENDDO
      !
      dos_CL(ie)=ZERO
      DO i = 1, dim_subspace
         dos_CL(ie) = dos_CL(ie) -  AIMAG( P1_g_tilde_CL_P1(i,i) ) / PI
      ENDDO
      !
      dos_CR(ie)=ZERO
      DO i = 1, dim_subspace
         dos_CR(ie) = dos_CR(ie) -  AIMAG( P1_g_tilde_CR_P1(i,i) ) / PI
      ENDDO

      !
      sigma_L(:,ie) = ZERO
      DO i = 1, dim_subspace
         sigma_L(i,ie) =  -  REAL(SUM( sgm_L(:,i) ))
      ENDDO
      !
      sigma_R(:,ie) = ZERO
      DO i = 1, dim_subspace
         sigma_R(i,ie) =  - REAL(SUM( sgm_R(:,i) ))
      ENDDO
      !
      sigma_CL(:,ie) = ZERO
      DO i = 1, dim_subspace
         sigma_CL(i,ie) = -  REAL(SUM( sgm_CL(:,i) ))
      ENDDO
      !
      sigma_CR(:,ie) = ZERO
      DO i = 1, dim_subspace
         sigma_CR(i,ie) =  - REAL(SUM( sgm_CR(:,i) ))
      ENDDO


      !
      gamma_R_sum(:,ie)=ZERO
      DO i = 1, dim_subspace
         gamma_R_sum(i,ie) = gamma_R(i,i)
      ENDDO
      !
      gamma_L_sum(:,ie)=ZERO
      DO i = 1, dim_subspace
         gamma_L_sum(i,ie) =  gamma_L(i,i)
      ENDDO
      !
      gamma_tilde_R_sum(:,ie)=ZERO
      DO i = 1, dim_subspace
         gamma_tilde_R_sum(i,ie) = gamma_tilde_R(i,i)
      ENDDO
      !
      gamma_tilde_L_sum(:,ie)=ZERO
      DO i = 1, dim_subspace
         gamma_tilde_L_sum(i,ie) = gamma_tilde_L(i,i)
      ENDDO
      !



   ENDDO energy_loop

   !
   ! close sgm file
   !


!
! ... write DOS and CONDUCT data on files
!

   filename = 'cond.dat'
   OPEN ( cond_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( cond_unit, '(30(f15.9))' ) egrid(ie), SUM( conduct(:,ie) ), (conduct(i_sub,ie), i_sub=1,dim_subspace)
   ENDDO
   CLOSE( cond_unit )
   filename = 'sigma_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( sigma_L(:,ie) ), (sigma_L(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'sigma_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( sigma_R(:,ie) ), (sigma_R(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'sigma_CL.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( sigma_CL(:,ie) ), (sigma_CL(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'sigma_CR.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( sigma_CR(:,ie) ), (sigma_CR(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'gamma_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( gamma_R_sum(:,ie) ), (gamma_R_sum(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'gamma_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( gamma_L_sum(:,ie) ), (gamma_L_sum(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'gamma_tilde_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( gamma_tilde_R_sum(:,ie) ), (gamma_tilde_R_sum(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'gamma_tilde_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( gamma_tilde_L_sum(:,ie) ), (gamma_tilde_L_sum(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( dos_L(:,ie) ), ( dos_L(i_L,ie), i_L=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( dos_R(:,ie) ), ( dos_R(i_R,ie), i_R=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_CL.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), dos_CL(ie)
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_CR.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), dos_CR(ie) 
   ENDDO
   CLOSE( dos_unit )
   filename = 'cut_iter.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO i_L = 1, dim_subspace
      WRITE ( dos_unit, '(15i4)' ) ( cut_iter_L(i_L,i_R),  cut_iter_R(i_L,i_R), i_R=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
!    filename = 'test.dat'
!    OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
!    DO ie = 1, ne
!        WRITE ( dos_unit, '(12(f15.9))' ) egrid(ie), (test_mat(i_sub,ie), i_sub=1,dim_subspace )
!    ENDDO
!    CLOSE( dos_unit )



!
! ...  Finalize timing
!
   CALL timing('conductor',OPR='stop')
   CALL timing_overview(stdout,LIST=global_list,MAIN_NAME='conductor')


!
!...  free memory
!

     !!!!!!!!!!!!!  DEBUG !!!!!!!!!!!!!!!
   DEALLOCATE ( test_mat, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating test_mat', ABS(ierr) )
   DEALLOCATE ( eig, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating eig', ABS(ierr) )

   DEALLOCATE ( z, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating z', ABS(ierr) )
     !!!!!!!!!!!!!  DEBUG !!!!!!!!!!!!!!!

   DEALLOCATE ( gamma_L_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_L_sum', ABS(ierr) )
   DEALLOCATE ( gamma_R_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_R_sum', ABS(ierr) )
   DEALLOCATE ( gamma_tilde_L_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_tilde_L_sum', ABS(ierr) )
   DEALLOCATE ( gamma_tilde_R_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_tilde_R_sum', ABS(ierr) )
   DEALLOCATE ( sigma_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sigma_L', ABS(ierr) )
   DEALLOCATE ( sigma_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sigma_R', ABS(ierr) )
   DEALLOCATE ( sigma_CL, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sigma_CL', ABS(ierr) )
   DEALLOCATE ( sigma_CR, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sigma_CR', ABS(ierr) )
   DEALLOCATE ( dos_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_L', ABS(ierr) )
   DEALLOCATE ( dos_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_R', ABS(ierr) )
   DEALLOCATE ( dos_CL, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_CL', ABS(ierr) )
   DEALLOCATE ( dos_CR, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_CR', ABS(ierr) )
        !
   DEALLOCATE ( conduct, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating conduct', ABS(ierr) )
   DEALLOCATE ( work, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating work', ABS(ierr) )
   DEALLOCATE ( work2, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating work2', ABS(ierr) )
   DEALLOCATE ( g_aux, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating g_aux', ABS(ierr) )
   DEALLOCATE ( B_aux, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating B_aux', ABS(ierr) )
        !
   DEALLOCATE ( rec_chain_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating rec_chain_L', ABS(ierr) )
   DEALLOCATE ( rec_chain_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating rec_chain_R', ABS(ierr) )
   DEALLOCATE ( rec_var_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating rec_var_L', ABS(ierr) )
   DEALLOCATE ( rec_var_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating rec_var_R', ABS(ierr) )
   DEALLOCATE ( cut_iter_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating cut_iter_L', ABS(ierr) )
   DEALLOCATE ( cut_iter_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating cut_iter_R', ABS(ierr) )
        !
   DEALLOCATE ( sgm_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sgm_L', ABS(ierr) )
   DEALLOCATE ( sgm_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sgm_R', ABS(ierr) )
   DEALLOCATE ( sgm_CL, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sgm_CL', ABS(ierr) )
   DEALLOCATE ( sgm_CR, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sgm_CR', ABS(ierr) )
   DEALLOCATE ( gamma_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_L', ABS(ierr) )
   DEALLOCATE ( gamma_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_R', ABS(ierr) )
   DEALLOCATE ( gamma_tilde_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_tilde_L', ABS(ierr) )
   DEALLOCATE ( gamma_tilde_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_tilde_R', ABS(ierr) )
        !
   DEALLOCATE ( gR, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gR', ABS(ierr) )
   DEALLOCATE ( gL, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gL', ABS(ierr) )
   DEALLOCATE ( P1_g_tilde_CR_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_tilde_CR_P1', ABS(ierr) )
   DEALLOCATE ( P1_g_tilde_CR_PN, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_tilde_CR_PN', ABS(ierr) )
   DEALLOCATE ( PN_g_tilde_CR_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating PN_g_tilde_CR_P1', ABS(ierr) )
   DEALLOCATE ( P1_g_tilde_CL_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_tilde_CL_P1', ABS(ierr) )
   DEALLOCATE ( P1_g_tilde_CL_PN, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_tilde_CL_PN', ABS(ierr) )
   DEALLOCATE ( PN_g_tilde_CL_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating PN_g_tilde_CL_P1', ABS(ierr) )
   DEALLOCATE ( P1_g_CC_PN, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_CC_PN', ABS(ierr) )
   DEALLOCATE ( PN_g_CC_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating PN_g_CC_P1', ABS(ierr) )

   CALL cleanup()

END PROGRAM conductor5
  
