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
   PROGRAM gamma_c
   !***********************************************
   USE constants,            ONLY : PI, ZERO, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1
   USE files_module,         ONLY : file_open, file_close
   USE in_matrix_module,     ONLY : in_matrix_allocate, read_matrix
   USE io_module,            ONLY : stdout, stdin, sgm_unit => aux_unit,   &
                                    dos_unit => aux1_unit, cond_unit => aux2_unit
   USE iotk_module
   USE kinds,                ONLY : dbl
   USE matrix_module,        ONLY : matrix_allocate, build_matrix, &
                                    A_matrix_term, B_matrix_term, B_matrix_renorm, A_matrix_renorm, &
                                    B_matrix_lead_bulk
   USE parameters,           ONLY : nstrx 
   USE T_control_module,     ONLY : nprint, a_analytique, b_analytique, &
                                    max_iter_renorm, max_iter_term, dim_subspace, &
                                    print_gamma, print_gamma_tilde, &
                                    method_sigma
   USE T_egrid_module,       ONLY : egrid_init, ne, egrid, delta_lead
   USE T_input_module,       ONLY : input_manager
   USE timing_module,        ONLY : timing, timing_overview, global_list, timing_upto_now
   USE util_module,          ONLY : mat_mul, mat_hdiag
   USE version_module,       ONLY : version_number
   USE T_gamma_module,       ONLY : gamma_allocate, gamma_print, gamma_tilde_print, &
                                    sigma_print, sigma_tilde_print, &
                                    sigma_s, sigma_tilde_s, gamma_s, gamma_tilde_s


   IMPLICIT NONE

   !
   ! local variables
   !
   COMPLEX(dbl)     :: ene_lead
   CHARACTER(nstrx) :: filename
   INTEGER          :: i, ie, ierr, ncount, iter, i_rows, i_cols, i_sub, j_sub, i_rows2, i_cols2
   REAL(dbl)        :: real_iter
   !
   COMPLEX(dbl),    ALLOCATABLE :: work(:,:)
   COMPLEX(dbl),    ALLOCATABLE :: g_aux(:,:)
   COMPLEX(dbl),    ALLOCATABLE :: B_aux(:,:)
   !
   INTEGER          :: dim_term
   INTEGER          :: niterx
   !
   ! Fonctions de Green
   !
   COMPLEX(dbl), ALLOCATABLE :: P1_g_tilde_P1(:,:), P1_g_tilde_PN(:,:), PN_g_tilde_P1(:,:)
   COMPLEX(dbl), ALLOCATABLE :: g_term(:,:)
   !
   ! Variables associees au calcul des Self energies
   !
   COMPLEX(dbl), ALLOCATABLE :: sgm_term(:,:)
   COMPLEX(dbl), ALLOCATABLE :: sgm_renorm(:,:)
   !
   ! Taux d'injection
   !
   COMPLEX(dbl), ALLOCATABLE :: gamma_aux(:,:)
   COMPLEX(dbl), ALLOCATABLE :: gamma_tilde_aux(:,:)
   !
   !
   ! Variables associees aux fichiers de sortie
   !
   REAL(dbl),    ALLOCATABLE :: sigma_test(:,:)
   REAL(dbl),    ALLOCATABLE :: gamma_sum(:,:)
   REAL(dbl),    ALLOCATABLE :: gamma_tilde_sum(:,:)
   REAL(dbl),    ALLOCATABLE :: dos_term(:,:), dos_renorm(:)
!    REAL(dbl),    ALLOCATABLE :: sigma_CR(:,:)
   !
   COMPLEX(dbl), ALLOCATABLE :: tot(:,:), h00_term(:,:)
   COMPLEX(dbl), ALLOCATABLE :: tott(:,:), aux00(:,:), aux01(:,:), s00(:,:)
   COMPLEX(dbl), ALLOCATABLE :: g_term_aux(:,:)


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
   CALL matrix_allocate()
   CALL build_matrix()
   CALL egrid_init()

   CALL summary( stdout )
   !
   ! local variable allocations
   !

   ALLOCATE ( dos_renorm(ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_renorm', ABS(ierr) )
   ALLOCATE ( dos_term(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_term', ABS(ierr) )
      !
   ALLOCATE ( sigma_test(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sigma_test', ABS(ierr) )
   ALLOCATE ( gamma_sum(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_sum', ABS(ierr) )
   ALLOCATE ( gamma_tilde_sum(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_tilde_sum', ABS(ierr) )
      !
   ALLOCATE ( work(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating work', ABS(ierr) )
   ALLOCATE ( g_aux(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating g_aux', ABS(ierr) )
   ALLOCATE ( B_aux(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating B_aux', ABS(ierr) )
      !
      !
   ALLOCATE ( sgm_term(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sgm_term', ABS(ierr) )
   ALLOCATE ( sgm_renorm(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sgm_renorm', ABS(ierr) )
      !
   ALLOCATE ( gamma_aux(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma', ABS(ierr) )
   ALLOCATE ( gamma_tilde_aux(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_tilde', ABS(ierr) )
      !
   ALLOCATE ( g_term(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating g_term', ABS(ierr) )
   ALLOCATE ( P1_g_tilde_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_tilde_P1', ABS(ierr) )
   ALLOCATE ( P1_g_tilde_PN(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_tilde_PN', ABS(ierr) )
   ALLOCATE ( PN_g_tilde_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating PN_g_tilde_P1', ABS(ierr) )
       !
       CALL gamma_allocate()
       !
   IF ( TRIM(method_sigma) == 'WanT' ) THEN 
         !
         dim_term = dim_subspace * max_iter_term
         !
         ALLOCATE ( h00_term(dim_term,dim_term), STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','allocating h_term', ABS(ierr) )
         ALLOCATE ( tot(dim_term,dim_term), STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','allocating tot', ABS(ierr) )
         ALLOCATE ( tott(dim_term,dim_term), STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','allocating tott', ABS(ierr) )
         ALLOCATE ( g_term_aux(dim_term,dim_term), STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','allocating g_term_aux', ABS(ierr) )
         ALLOCATE ( aux00(dim_term,dim_term), STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','allocating aux00', ABS(ierr) )
         ALLOCATE ( aux01(dim_term,dim_term), STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','allocating aux01', ABS(ierr) )
         ALLOCATE ( s00(dim_term,dim_term), STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','allocating s00', ABS(ierr) )


         niterx=200
         !
         DO i=1, max_iter_term
            !
            DO i_sub=1, dim_subspace
               !
               i_rows = (dim_subspace * (i-1)) + i_sub
               !
               DO j_sub=1, dim_subspace
                  !
                  i_cols = (dim_subspace * (i-1)) + j_sub
                  !
                  h00_term(i_rows,i_cols) = A_matrix_term(i,i_sub,j_sub)
                  !
               ENDDO
               !
            ENDDO
            !
         ENDDO
         !
         IF (max_iter_term >1) THEN
            DO i=1, (max_iter_term-1)
               !
               DO i_sub=1, dim_subspace
                  !
                  i_rows = (dim_subspace * (i-1)) + i_sub
                  !
                  DO j_sub=1, dim_subspace
                     !
                     i_cols = (dim_subspace * (i)) + j_sub
                     !
                     h00_term(i_rows,i_cols) = CONJG (B_matrix_term(i,j_sub,i_sub))
                     h00_term(i_cols,i_rows) = B_matrix_term(i,j_sub,i_sub)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
         ENDIF
         !
         !
         i= max_iter_term
            !
            DO i_sub=1, dim_subspace
               !
               i_rows = (dim_subspace * (i-1)) + i_sub
               !
               DO j_sub=1, dim_subspace
                  !
                  i_cols = j_sub
                  !
                  aux01(i_rows,i_cols) = CONJG (B_matrix_term(i,j_sub,i_sub))
                  !
               ENDDO
               !
            ENDDO
            !
         !
         s00(:,:) = CZERO
         !
         DO i=1, dim_term
            !
            s00(i,i) = CONE
            !
         ENDDO
   ENDIF
   !
   !PRINT*, 'h00_term'
   !PRINT*, h00_term(:,:)
   !PRINT*, 'aux01'
   !PRINT*, aux01(:,:)
   !



   energy_loop: &
   DO ie = 1, ne
      ncount = ie
      !
      ! grids and misc
      !
      ene_lead = egrid(ie)  + delta_lead * CI
           IF ( MOD( ncount, nprint) == 0 .OR. ncount == 1 ) THEN
                WRITE(stdout,"(2x, 'Computing E( ',i5,' ) = ', f9.5, ' eV' )") ncount, egrid(ie)
                CALL timing_upto_now(stdout)
           ENDIF

      ! 
      ! construct leads self-energies 
      ! 

      IF ( TRIM(method_sigma) == 'analytique' ) THEN
         !
         CALL green_analytique(ene_lead, g_term, dim_subspace, a_analytique, b_analytique)
         !
      ELSE IF ( TRIM(method_sigma) == 'WanT' ) THEN
         !
         aux00(:,:)  = h00_term(:,:) - ene_lead * s00(:,:)
         !
         CALL transfer( dim_term, niterx, tot, tott, aux00, aux01 )
         CALL green_want( dim_term, tot, tott, aux00, aux01, ene_lead, g_term_aux, 1 )
         !
         g_term(1:dim_subspace, 1:dim_subspace) =  g_term_aux(1:dim_subspace, 1:dim_subspace)
         !
      ELSE
         !
         CALL errore('conductor','incorrect method_sigma',1 )
         !
      ENDIF
      !
      !PRINT*, B_matrix_lead_bulk(1,:,:)

      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_lead_bulk(1,:,:)) )
      CALL mat_mul(work, B_aux, 'N', g_term, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(sgm_term, work, 'N', B_matrix_lead_bulk(1,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)
      !  gamma
      gamma_aux(:,:) = CI * (  sgm_term(:,:) - CONJG( TRANSPOSE( sgm_term(:,:) ) )   )
      CALL green_tridiag( P1_g_tilde_P1(:,:), P1_g_tilde_PN(:,:), PN_g_tilde_P1(:,:), ene_lead, A_matrix_renorm(:,:,:), B_matrix_renorm(:,:,:), sgm_term(:,:), max_iter_renorm, dim_subspace) 
      ! Calcul des Sigma CR et des Sigma_CL
      !B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_C_lead(1,:,:)) )
      !CALL mat_mul(work, B_aux, 'N', P1_g_tilde_P1, 'N', dim_subspace, dim_subspace, dim_subspace)
      !CALL mat_mul(sgm_renorm, work, 'N', B_matrix_C_lead(1,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)
      sgm_renorm(:,:) = P1_g_tilde_P1(:,:)
      ! gamma_tilde_
      ! gamma_tilde_R = H_C_CR * G^r_CR (z) * Gamma_R * G^a_CR (z) * H_CR_C
      !B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_C_lead(1,:,:)) )
      !CALL mat_mul(work, B_aux, 'N', P1_g_tilde_PN, 'N', dim_subspace, dim_subspace, dim_subspace)
      !CALL mat_mul(gamma_tilde, work, 'N',B_matrix_C_lead(1,:,:) , 'N', dim_subspace, dim_subspace, dim_subspace)
      ! gamma_tilde_L = H_C_CL *P1 G^a_CL (z) * Gamma_L * G^r_CL (z) * H_CL_C
      CALL mat_mul(work, P1_g_tilde_PN, 'N', gamma_aux, 'N', dim_subspace, dim_subspace, dim_subspace)
      g_aux(:,:)=  CONJG( TRANSPOSE(P1_g_tilde_PN(:,:)))
      CALL mat_mul(gamma_tilde_aux, work, 'N', g_aux, 'N', dim_subspace, dim_subspace, dim_subspace)


      ! Fichiers de sortie

      IF ( print_gamma_tilde ) THEN 
            !
            sigma_tilde_s(ie,:,:) =  sgm_renorm(:,:)
            !
            gamma_tilde_s(ie,:,:) = gamma_tilde_aux(:,:)
            !
      ENDIF

      IF ( print_gamma ) THEN 
            !
            sigma_s(ie,:,:) =  sgm_term(:,:)
            !
            gamma_s(ie,:,:) = gamma_aux(:,:)
            !
      ENDIF
      !
      !
      dos_term(:,ie)=ZERO
      DO i = 1, dim_subspace
         dos_term(i,ie) =  -  AIMAG ( g_term(i,i) ) / PI
      ENDDO
      !
      !
      dos_renorm(ie)=ZERO
      DO i = 1, dim_subspace
         dos_renorm(ie) = dos_renorm(ie) -  AIMAG( P1_g_tilde_P1(i,i) ) / PI
      ENDDO
      !
      !
      sigma_test(:,ie) = ZERO
      DO i = 1, dim_subspace
         sigma_test(i,ie) =  - REAL(SUM( sgm_term(:,i) ))
      ENDDO
      !
      !
      !sigma_renorm(:,ie) = ZERO
      !DO i = 1, dim_subspace
      !   sigma_renorm(i,ie) =  - REAL(SUM( sgm_renorm(:,i) ))
      !ENDDO
      !
      !
      gamma_sum(:,ie)=ZERO
      DO i = 1, dim_subspace
         gamma_sum(i,ie) = gamma_aux(i,i)
      ENDDO
      !
      gamma_tilde_sum(:,ie)=ZERO
      DO i = 1, dim_subspace
         gamma_tilde_sum(i,ie) = gamma_tilde_aux(i,i)
      ENDDO
      !
      !
   ENDDO energy_loop

   IF ( print_gamma ) THEN 
         !
         !
         filename = 'gamma.dat'
         OPEN ( cond_unit, FILE=TRIM(filename), FORM='formatted' )
         CALL gamma_print(cond_unit,TRIM(filename))
         CLOSE( cond_unit )
         !
         !
         filename = 'sigma.dat'
         OPEN ( cond_unit, FILE=TRIM(filename), FORM='formatted' )
         CALL sigma_print(cond_unit,TRIM(filename))
         CLOSE( cond_unit )
         !
   ENDIF   

   IF ( print_gamma_tilde ) THEN 
         !
         !
         filename = 'gamma_tilde.dat'
         OPEN ( cond_unit, FILE=TRIM(filename), FORM='formatted' )
         CALL gamma_tilde_print(cond_unit,TRIM(filename))
         CLOSE( cond_unit )
         !
         !
         filename = 'sigma_tilde.dat'
         OPEN ( cond_unit, FILE=TRIM(filename), FORM='formatted' )
         CALL sigma_tilde_print(cond_unit,TRIM(filename))
         CLOSE( cond_unit )
         !
   ENDIF   


   filename = 'sigma_test.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( sigma_test(:,ie) ), (sigma_test(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
!    filename = 'sigma_renorm.dat'
!    OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
!    DO ie = 1, ne
!        WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( sigma_renorm(:,ie) ), (sigma_renorm(i_sub,ie), i_sub=1,dim_subspace )
!    ENDDO
!    CLOSE( dos_unit )
   filename = 'gamma_test.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( gamma_sum(:,ie) ), (gamma_sum(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'gamma_tilde_test.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( gamma_tilde_sum(:,ie) ), (gamma_tilde_sum(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_term.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(10(f15.9))' ) egrid(ie), SUM( dos_term(:,ie) ), ( dos_term(i_sub,ie), i_sub=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_renorm.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), dos_renorm(ie) 
   ENDDO
   CLOSE( dos_unit )



!
! ...  Finalize timing
!
   CALL timing('conductor',OPR='stop')
   CALL timing_overview(stdout,LIST=global_list,MAIN_NAME='conductor')


!
!...  free memory
!
   IF (TRIM(method_sigma) == 'WanT')   THEN
        !PRINT*, '1'
        DEALLOCATE ( tot, STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','deallocating tot', ABS(ierr) )
        !PRINT*, '2'
        DEALLOCATE ( tott, STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','deallocating tott', ABS(ierr) )
        !PRINT*, '3'
        DEALLOCATE ( g_term_aux, STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','deallocating g_term_aux', ABS(ierr) )
        !PRINT*, '4'
        DEALLOCATE ( h00_term, STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','deallocating h_term', ABS(ierr) )
        !PRINT*, '5'
        DEALLOCATE ( aux00, STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','deallocating aux00', ABS(ierr) )
        !PRINT*, '6'
        DEALLOCATE ( aux01, STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','deallocating aux01', ABS(ierr) )
        !PRINT*, '7'
        DEALLOCATE ( s00, STAT=ierr )
            IF( ierr /=0 ) CALL errore('conductor','deallocating s00', ABS(ierr) )
   ENDIF
   !PRINT*, '1'
   DEALLOCATE ( dos_term, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_term', ABS(ierr) )
   !PRINT*, '2'
   DEALLOCATE ( dos_renorm, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_renorm', ABS(ierr) )
        !
   !PRINT*, '3'
   DEALLOCATE ( work, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating work', ABS(ierr) )
   !PRINT*, '4'
   DEALLOCATE ( g_aux, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating g_aux', ABS(ierr) )
   !PRINT*, '5'
   DEALLOCATE ( B_aux, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating B_aux', ABS(ierr) )
        !
        !
   !PRINT*, '6'
   DEALLOCATE ( sigma_test, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sigma', ABS(ierr) )
!    DEALLOCATE ( sigma_renorm, STAT=ierr )
!         IF( ierr /=0 ) CALL errore('conductor','deallocating sigma_renorm', ABS(ierr) )
   !PRINT*, '7'
   DEALLOCATE ( sgm_term, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sgm_term', ABS(ierr) )
   !PRINT*, '8'
   DEALLOCATE ( sgm_renorm, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sgm_renorm', ABS(ierr) )
   !PRINT*, '9'
   DEALLOCATE ( gamma_aux, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma', ABS(ierr) )
   !PRINT*, '1'
   DEALLOCATE ( gamma_tilde_aux, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_tilde', ABS(ierr) )
   !PRINT*, '2'
   DEALLOCATE ( gamma_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_sum', ABS(ierr) )
   !PRINT*, '3'
   DEALLOCATE ( gamma_tilde_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_tilde_sum', ABS(ierr) )
   !PRINT*, '4'
        !
   DEALLOCATE ( g_term, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating g', ABS(ierr) )
   !PRINT*, '5'
   DEALLOCATE ( P1_g_tilde_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_tilde_P1', ABS(ierr) )
   !PRINT*, '6'
   DEALLOCATE ( P1_g_tilde_PN, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_tilde_PN', ABS(ierr) )
   !PRINT*, '7'
   DEALLOCATE ( PN_g_tilde_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating PN_g_tilde_P1', ABS(ierr) )
   !PRINT*, '8'
   CALL cleanup()

END PROGRAM gamma_c
  
