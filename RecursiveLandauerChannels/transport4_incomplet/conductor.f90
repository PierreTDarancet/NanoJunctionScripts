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
   PROGRAM conductor4
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
   USE matrix_module,        ONLY : matrix_allocate, build_matrix
   USE parameters,           ONLY : nstrx 
   USE T_control_module,     ONLY : conduct_formula, nprint, &
                                    max_iter_L, max_iter_R, dim_subspace, &
                                    cut_chain
   USE T_egrid_module,       ONLY : egrid_init, ne, egrid, delta, delta_lead
   USE T_hamiltonian_module, ONLY : dimL, dimR, dimC,  &
                                    h00_C, h_LC, h_CR, s00_C, &
                                    hamiltonian_init, hamiltonian_allocate
   USE T_input_module,       ONLY : input_manager
   USE T_workspace_module,   ONLY : aux00_C, aux_LC, aux_CL, aux_CR, aux_RC,    &
                                    gR, gL, gC, gamma_R, gamma_L, sgm_L, sgm_R, &
                                    workspace_allocate
   USE timing_module,        ONLY : timing, timing_overview, global_list, timing_upto_now
   USE util_module,          ONLY : mat_mul, mat_sv, mat_hdiag
   USE version_module,       ONLY : version_number

   IMPLICIT NONE

   !
   ! local variables
   !
   COMPLEX(dbl)     :: ene, ene_lead
   CHARACTER(nstrx) :: filename
   INTEGER          :: i, ie, ik, ierr, ncount, i_L, i_R, iter, i_rows, i_cols
   REAL(dbl)        :: real_iter
   !
   INTEGER,      ALLOCATABLE :: cut_iter_R(:,:), cut_iter_L(:,:)
   REAL(dbl),    ALLOCATABLE :: dos(:,:), conduct(:,:), sigma_R(:,:), sigma_L(:,:), gamma_L_sum(:,:), gamma_R_sum(:,:)
   REAL(dbl),    ALLOCATABLE :: dos_L(:,:), dos_R(:,:)
   REAL(dbl),    ALLOCATABLE :: rec_var_L(:,:,:), rec_var_R(:,:,:)
   REAL(dbl),    ALLOCATABLE :: cond_aux(:)
   COMPLEX(dbl), ALLOCATABLE :: work(:,:), rec_chain_L(:,:,:,:), rec_chain_R(:,:,:,:)
   !

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
   ! Set up the layer hamiltonians
   CALL hamiltonian_allocate()
   CALL hamiltonian_init()
   ! write input data on output file
   CALL summary( stdout )
   !
   ! local variable allocations
   !
   CALL workspace_allocate()


   ALLOCATE ( dos(dimC,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos', ABS(ierr) )
   ALLOCATE ( dos_L(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_L', ABS(ierr) )
   ALLOCATE ( dos_R(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating dos_R', ABS(ierr) )
   ALLOCATE ( sigma_R(dimC,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sigma_R', ABS(ierr) )
   ALLOCATE ( sigma_L(dimC,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating sigma_L', ABS(ierr) )
   ALLOCATE ( conduct(dimC,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating conduct', ABS(ierr) )
   ALLOCATE ( cond_aux(dimC), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating cond_aux', ABS(ierr) )
   ALLOCATE ( work(dimC,dimC), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating work', ABS(ierr) )
   ALLOCATE ( gamma_R_sum(dimC,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_R_sum', ABS(ierr) )
   ALLOCATE ( gamma_L_sum(dimC,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating gamma_L_sum', ABS(ierr) )
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

   WRITE(stdout,"()")


   WRITE(stdout,"()")
   WRITE(stdout, "('Rec on h00_R...')")
   CALL scalar_recursion(rec_chain_R(:,:,:,:), nb_first_R(:), id_first_R(:,:), ene_first_R(:,:), ene_onsite_R(:), dimR, max_iter_R, dim_subspace, rec_var_R(:,:,:), cut_chain, cut_iter_R(:,:) )
   WRITE(stdout,"()")
   !
   filename = 'rec_var_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO i_cols=1,dim_subspace 
        !
        WRITE( dos_unit, "('icols = ',i5,' ')") i_cols
        !
        DO i_rows = 1, dim_subspace
            !
            WRITE( dos_unit, "('irows = ',i5,' ')") i_rows
            !
            DO iter = 1, max_iter_R
               !
               real_iter = REAL(iter, dbl)
               !
               WRITE ( dos_unit, '(6(f15.9))' ) real_iter, rec_var_R(iter, i_rows, i_cols), rec_chain_R(1,iter,i_rows, i_cols), rec_chain_R(2,iter,i_rows, i_cols)
               !
            ENDDO
       ENDDO


   ENDDO
   CLOSE( dos_unit )
   !
   !
   WRITE(stdout, "('Rec on  h00_L...')")
   CALL scalar_recursion(rec_chain_L(:,:,:,:), nb_first_L(:), id_first_L(:,:), ene_first_L(:,:), ene_onsite_L(:), dimL, max_iter_L, dim_subspace, rec_var_L(:,:,:), cut_chain, cut_iter_L(:,:) )
   WRITE(stdout,"()")
   !
   filename = 'rec_var_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )

   DO i_cols=1,dim_subspace
        !
        WRITE( dos_unit, "('icols = ',i5,' ')") i_cols
        !
        DO i_rows = 1, dim_subspace
            !
            WRITE( dos_unit, "('irows = ',i5,' ')") i_rows
            !
            DO iter = 1, max_iter_L
               !
               real_iter = REAL(iter, dbl)
               !
               WRITE ( dos_unit, '(6(f15.9))' ) real_iter, rec_var_L(iter, i_rows, i_cols), rec_chain_L(1,iter,i_rows, i_cols), rec_chain_L(2,iter,i_rows, i_cols)
               !
            ENDDO
       ENDDO
   ENDDO
   CLOSE( dos_unit )
   !

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

      dos(:,ie) = ZERO
      conduct(:,ie) = ZERO
      sigma_L(:,ie) = ZERO
      sigma_R(:,ie) = ZERO

      !
      aux00_C(:,:)  = h00_C(:,:)  -ene * s00_C(:,:)
      !
      aux_LC(:,:) = h_LC(:,:) 
      aux_CR(:,:) = h_CR(:,:) 
      !
      aux_CL(:,:) = CONJG( TRANSPOSE( h_LC(:,:) ) )
      aux_RC(:,:) = CONJG( TRANSPOSE( h_CR(:,:) ) )


      ! 
      ! construct leads self-energies 
      ! 
      ! ene

      CALL green_tout_neuf( rec_chain_R, ene_lead, gR, dim_subspace, max_iter_R, cut_chain, cut_iter_R)
      !
      CALL mat_mul(work, aux_CR, 'N', gR, 'N', dimC, dim_subspace, dim_subspace)
      CALL mat_mul(sgm_R, work, 'N', aux_RC, 'N', dimC, dimC, dim_subspace)

      ! ene
      CALL green_tout_neuf( rec_chain_L, ene_lead, gL, dim_subspace, max_iter_L, cut_chain, cut_iter_L)
      !
      CALL mat_mul(work, aux_CL, 'N', gL, 'N', dimC, dim_subspace, dim_subspace)
      CALL mat_mul(sgm_L, work, 'N', aux_LC, 'N', dimC, dimC, dim_subspace) 

      !
      ! gamma_L and gamma_R
      !
      gamma_L(:,:) = CI * (  sgm_L(:,:) - CONJG( TRANSPOSE(sgm_L(:,:)) )   )
      gamma_R(:,:) = CI * (  sgm_R(:,:) - CONJG( TRANSPOSE(sgm_R(:,:)) )   )

      !
      ! Construct the conductor green's function
      ! gC = work^-1  (retarded)
      !
      work(1:dimC,1:dimC) = -aux00_C(:,:) -sgm_L(:,:) -sgm_R(:,:) 

      gC(:,:) = CZERO
      DO i = 1, dimC
         gC(i,i)= CONE
      ENDDO

      CALL mat_sv(dimC, dimC, work, gC)

      !
      ! Compute density of states for the conductor layer
      !
      DO i = 1, dimC
         dos(i,ie) = dos(i,ie) -  AIMAG( gC(i,i) ) / PI
      ENDDO
      !
      DO i = 1, dim_subspace
         dos_L(i,ie) = dos_L(i,ie) -  AIMAG( gL(i,i) ) / PI
      ENDDO
      !
      DO i = 1, dimC
         sigma_L(i,ie) = sigma_L(i,ie) -  REAL(SUM( sgm_L(:,i) ))
      ENDDO
      !
      DO i = 1, dim_subspace
         dos_R(i,ie) = dos_R(i,ie) -  AIMAG ( gR(i,i) ) / PI
      ENDDO
      !
      DO i = 1, dimC
         sigma_R(i,ie) = sigma_R(i,ie) - REAL(SUM( sgm_R(:,i) ))
      ENDDO
      !
      DO i = 1, dimC
         gamma_R_sum(i,ie) = gamma_R_sum(i,ie) + gamma_R(i,i)
      ENDDO
      !
      DO i = 1, dimC
         gamma_L_sum(i,ie) = gamma_L_sum(i,ie) +  gamma_L(i,i)
      ENDDO
      !

      !
      CALL transmittance(dimC, gamma_L, gamma_R, gC,  TRIM(conduct_formula), cond_aux )
      conduct(:,ie) = conduct(:,ie) + cond_aux(:)


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
       WRITE ( cond_unit, '(2(f15.9))' ) egrid(ie), SUM( conduct(:,ie) )
   ENDDO
   CLOSE( cond_unit )

   filename = 'dos.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( dos(:,ie) )
   ENDDO
   CLOSE( dos_unit )
   filename = 'sigma_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( sigma_L(:,ie) )
   ENDDO
   CLOSE( dos_unit )
   filename = 'sigma_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( sigma_R(:,ie) )
   ENDDO
   CLOSE( dos_unit )
   filename = 'gamma_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( gamma_R_sum(:,ie) )
   ENDDO
   CLOSE( dos_unit )
   filename = 'gamma_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( gamma_L_sum(:,ie) )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_L.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( dos_L(:,ie) )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_R.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( dos_R(:,ie) )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_L_decomp.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(7(f15.9))' ) egrid(ie), ( dos_L(i_L,ie), i_L=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'dos_R_decomp.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( dos_unit, '(7(f15.9))' ) egrid(ie), ( dos_R(i_R,ie), i_R=1,dim_subspace )
   ENDDO
   CLOSE( dos_unit )
   filename = 'cut_iter.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   DO i_L = 1, dim_subspace
      WRITE ( dos_unit, '(15i4)' ) ( cut_iter_L(i_L,i_R),  cut_iter_R(i_L,i_R), i_R=1,dim_subspace )
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
   DEALLOCATE ( gamma_L_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_L_sum', ABS(ierr) )
   DEALLOCATE ( gamma_R_sum, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating gamma_R_sum', ABS(ierr) )
   DEALLOCATE ( dos_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_L', ABS(ierr) )
   DEALLOCATE ( dos_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos_R', ABS(ierr) )
   DEALLOCATE ( dos, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos', ABS(ierr) )
   DEALLOCATE ( sigma_L, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sigma', ABS(ierr) )
   DEALLOCATE ( sigma_R, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating sigma', ABS(ierr) )
   DEALLOCATE ( conduct, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating conduct', ABS(ierr) )
   DEALLOCATE ( cond_aux, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating cond_aux', ABS(ierr) )
   DEALLOCATE ( work, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating work', ABS(ierr) )
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

   CALL cleanup()

END PROGRAM conductor4
  
