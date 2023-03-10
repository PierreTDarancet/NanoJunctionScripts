!
!      Copyright (C) 2004 WanT Group
!      Copyright (C) 1999 Marco Buongiorno Nardelli
!      Copyright (C) 2007-8 Pierre Darancet Institut N?el
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!***********************************************
   PROGRAM conductor6_imp_gamma
   !***********************************************
   USE constants,            ONLY : PI, ZERO, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1
   USE files_module,         ONLY : file_open, file_close
   USE in_matrix_module,     ONLY : in_matrix_allocate, read_matrix,  A_matrix_C, B_matrix_C
   USE io_module,            ONLY : stdout, stdin, sgm_unit => aux_unit,   &
                                    dos_unit => aux1_unit, cond_unit => aux2_unit
   USE iotk_module
   USE kinds,                ONLY : dbl
   USE parameters,           ONLY : nstrx 
   USE T_control_module,     ONLY : nprint,  dim_subspace, in_max_iter_C, &
                                    imp_gamma, print_gamma
   USE T_egrid_module,       ONLY : egrid_init, ne, egrid, delta
   USE T_input_module,       ONLY : input_manager
   USE timing_module,        ONLY : timing, timing_overview, global_list, timing_upto_now
   USE version_module,       ONLY : version_number
   USE T_gamma_module,       ONLY : gamma_allocate, gamma_L_read, gamma_R_read, &
                                    sigma_L_read, sigma_R_read, &
                                    sigma_L, sigma_R, gamma_L_s, gamma_R_s


   IMPLICIT NONE

   !
   ! local variables
   !
   COMPLEX(dbl)     :: ene
   CHARACTER(nstrx) :: filename
   INTEGER          :: i, ie, ierr, ncount, i_sub
   !
   REAL(dbl),    ALLOCATABLE :: conduct(:,:)
   !
   ! Fonctions de Green
   !
   COMPLEX(dbl), ALLOCATABLE :: P1_g_CC_PN(:,:), PN_g_CC_P1(:,:)

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
   CALL egrid_init()

   CALL summary( stdout )
   !
   ! local variable allocations
   !

!       ALLOCATE ( dos_CC(dim_subspace,ne), STAT=ierr )
!          IF( ierr /=0 ) CALL errore('conductor','allocating dos_CC', ABS(ierr) )
      !
      !
   ALLOCATE ( conduct(dim_subspace,ne), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating conduct', ABS(ierr) )
      !
   ALLOCATE ( P1_g_CC_PN(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating P1_g_CC_PN', ABS(ierr) )
   ALLOCATE ( PN_g_CC_P1(dim_subspace,dim_subspace), STAT=ierr )
      IF( ierr /=0 ) CALL errore('conductor','allocating PN_g_CC_P1', ABS(ierr) )

   !
   IF ( imp_gamma ) THEN 
       !
       CALL gamma_allocate()
       !
       PRINT*, '1'
       filename = 'gamma_L.dat'
       CALL gamma_L_read(cond_unit,TRIM(filename))
       !
       PRINT*, '2'
       filename = 'gamma_R.dat'
       CALL gamma_R_read(cond_unit,TRIM(filename))
       !
       PRINT*, '3'
       filename = 'sigma_L.dat'
       CALL sigma_L_read(cond_unit,TRIM(filename))
       !
       PRINT*, '4'
       filename = 'sigma_R.dat'
       CALL sigma_R_read(cond_unit,TRIM(filename))
       !
   ENDIF   

   energy_loop: &
   DO ie = 1, ne
      ncount = ie
      !
      ! grids and misc
      !
      ene =  egrid(ie)  + delta * CI
           IF ( MOD( ncount, nprint) == 0 .OR. ncount == 1 ) THEN
                WRITE(stdout,"(2x, 'Computing E( ',i5,' ) = ', f9.5, ' eV' )") ncount, egrid(ie)
                CALL timing_upto_now(stdout)
           ENDIF

      conduct(:,ie) = ZERO
      ! 
      ! construct leads self-energies 
      ! 

      !
      ! Calcul de G_CC
      !
      CALL green_CC_tridiag( P1_g_CC_PN(:,:), PN_g_CC_P1(:,:), ene, A_matrix_C(:,:,:), B_matrix_C(:,:,:), sigma_L(ie,:,:), sigma_R(ie,:,:), in_max_iter_C, dim_subspace) 

      !
      CALL transmittance_min(dim_subspace, gamma_L_s(ie,:,:), gamma_R_s(ie,:,:), P1_g_CC_PN(:,:), conduct(:,ie) )
      !

   ENDDO energy_loop


!
! ... write DOS and CONDUCT data on files
!

   filename = 'cond.dat'
   OPEN ( cond_unit, FILE=TRIM(filename), FORM='formatted' )
   DO ie = 1, ne
       WRITE ( cond_unit, '(30(f15.9))' ) egrid(ie), SUM( conduct(:,ie) ), (conduct(i_sub,ie), i_sub=1,dim_subspace)
   ENDDO
   CLOSE( cond_unit )


!
! ...  Finalize timing
!
   CALL timing('conductor',OPR='stop')
   CALL timing_overview(stdout,LIST=global_list,MAIN_NAME='conductor')


!
!...  free memory
!

   DEALLOCATE ( conduct, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating conduct', ABS(ierr) )
        !
   DEALLOCATE ( P1_g_CC_PN, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating P1_g_CC_PN', ABS(ierr) )
   DEALLOCATE ( PN_g_CC_P1, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating PN_g_CC_P1', ABS(ierr) )

   CALL cleanup()

END PROGRAM conductor6_imp_gamma
  
