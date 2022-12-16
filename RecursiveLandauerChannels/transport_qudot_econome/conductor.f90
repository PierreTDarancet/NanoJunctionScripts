!
!      Copyright (C) 2004 WanT Group
!      Copyright (C) 1999 Marco Buongiorno Nardelli
!      Copyright (C) 2007-8 Pierre Darancet Institut Néel
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!***********************************************
   PROGRAM conductor_qudot_imp_gamma
   !***********************************************
   USE constants,            ONLY : PI, ZERO, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1
   USE files_module,         ONLY : file_open, file_close
   USE in_matrix_module,     ONLY : in_matrix_allocate, read_matrix,  &
                                    m_calcul_gamma, m_calcul_sigma, m_calcul_g_CC, &
                                    m_calcul_transmittance
   USE io_module,            ONLY : stdout, stdin, sgm_unit => aux_unit,   &
                                    dos_unit => aux1_unit, cond_unit => aux2_unit
   USE iotk_module
   USE kinds,                ONLY : dbl
   USE parameters,           ONLY : nstrx 
   USE T_control_module,     ONLY : nprint,  dim_subspace, in_max_iter_C, &
                                    imp_gamma, print_gamma, dimC, datafile_C_form, &
                                    in_datafile_gamma_L, in_datafile_gamma_R, &
                                    in_datafile_sigma_L, in_datafile_sigma_R
   USE T_egrid_module,       ONLY : egrid_init, ne, egrid, delta

   USE in_hamiltonian_module, ONLY : hamiltonian_allocate, read_hamiltonian, &
                                    h_calcul_transmittance
   USE T_input_module,       ONLY : input_manager
   USE timing_module,        ONLY : timing, timing_overview, global_list, timing_upto_now
   USE version_module,       ONLY : version_number
   USE T_gamma_module,       ONLY : gamma_allocate, gamma_L_read, gamma_R_read, &
                                    sigma_L_read, sigma_R_read, &
                                    gL => sigma_L, gR => sigma_R, gamma_L_s, gamma_R_s


   IMPLICIT NONE

   !
   ! local variables
   !
   COMPLEX(dbl)     :: ene
   CHARACTER(nstrx) :: filename
   INTEGER          :: i, ie, ierr, ncount, i_sub
   !
   REAL(dbl),    ALLOCATABLE :: conduct(:,:)
   REAL(dbl),    ALLOCATABLE :: dos(:,:)
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
   IF  (TRIM(datafile_C_form) =='matrix') THEN 
      !
      ALLOCATE ( conduct(ne,dim_subspace), STAT=ierr )
         IF( ierr /=0 ) CALL errore('conductor','allocating conduct', ABS(ierr) )
      ALLOCATE ( dos(ne,dim_subspace), STAT=ierr )
         IF( ierr /=0 ) CALL errore('conductor','allocating dos', ABS(ierr) )

      CALL in_matrix_allocate()
      CALL read_matrix()
      !
   ELSE IF (TRIM(datafile_C_form) =='hamiltonian') THEN 
      !
      ALLOCATE ( conduct(ne,dimC), STAT=ierr )
         IF( ierr /=0 ) CALL errore('conductor','allocating conduct', ABS(ierr) )
      ALLOCATE ( dos(ne,dimC), STAT=ierr )
         IF( ierr /=0 ) CALL errore('conductor','allocating dos', ABS(ierr) )

      CALL hamiltonian_allocate()
      CALL read_hamiltonian()
      !
   ELSE
      !
      CALL errore('conductor','problem in datafile_C_form', 10 )
      !
   ENDIF
   !
   conduct(:,:) = ZERO
   dos(:,:) = ZERO
   !
   CALL egrid_init()

   CALL summary( stdout )
   !
   ! local variable allocations
   !

   IF ( imp_gamma ) THEN 
       !
       CALL gamma_allocate()
       !
       PRINT*, '1'
       filename = TRIM(in_datafile_gamma_L)
       CALL gamma_L_read(cond_unit,TRIM(filename))
       !
       PRINT*, '2'
       filename =  TRIM(in_datafile_gamma_R)
       CALL gamma_R_read(cond_unit,TRIM(filename))
       !
       PRINT*, '3'
       filename =  TRIM(in_datafile_sigma_L)
       CALL sigma_L_read(cond_unit,TRIM(filename))
       !
       PRINT*, '4'
       filename =  TRIM(in_datafile_sigma_R)
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

      ! 
      ! construct leads self-energies 
      ! 
      IF (TRIM(datafile_C_form) =='matrix') THEN 
         !
         CALL m_calcul_gamma( gamma_L_s(ie,:,:) , gamma_R_s(ie,:,:) , ie )
         CALL m_calcul_sigma( gL(ie,:,:) , gR(ie,:,:), ie)
         CALL m_calcul_g_CC( dos(ie,:) , ene, ie)
         CALL m_calcul_transmittance( conduct(ie,:) , ie)
         !
      ELSE IF (TRIM(datafile_C_form) =='hamiltonian') THEN 
         !
!          CALL h_calcul_gamma( gamma_L_s(ie,:,:) , gamma_R_s(ie,:,:) , ie )
!          CALL h_calcul_sigma( gL(ie,:,:) , gR(ie,:,:), ie)
!          CALL h_calcul_g_CC( dos(ie,:) , ene, ie)
         CALL h_calcul_transmittance( gamma_L_s(ie,:,:) , gamma_R_s(ie,:,:) , gL(ie,:,:) , gR(ie,:,:), dos(ie,:) , ene, conduct(ie,:) )
         !
      ELSE
         CALL errore('conductor','problem in datafile_C_form', 10 )
      ENDIF

   ENDDO energy_loop
!
! ... write DOS and CONDUCT data on files
!

   filename = 'cond.dat'
   OPEN ( cond_unit, FILE=TRIM(filename), FORM='formatted' )
   !
   IF (TRIM(datafile_C_form) =='matrix') THEN 
      DO ie = 1, ne
         WRITE ( cond_unit, '(30(f15.9))' ) egrid(ie), SUM( conduct(ie,:) ), (conduct(ie,i_sub), i_sub=1,dim_subspace)
      ENDDO
   ELSE IF (TRIM(datafile_C_form) =='hamiltonian') THEN 
      DO ie = 1, ne
         WRITE ( cond_unit, '(30(f15.9))' ) egrid(ie), SUM( conduct(ie,:) )
      ENDDO
   ELSE
      CALL errore('conductor','problem in datafile_C_form', 10 )
   ENDIF
   !
   CLOSE( cond_unit )
   !
   filename = 'dos.dat'
   OPEN ( dos_unit, FILE=TRIM(filename), FORM='formatted' )
   !
   IF (TRIM(datafile_C_form) =='matrix') THEN 
      DO ie = 1, ne
         WRITE ( dos_unit, '(30(f15.9))' ) egrid(ie), SUM( dos(ie,:) ), (dos(ie,i_sub), i_sub=1,dim_subspace)
      ENDDO
   ELSE IF (TRIM(datafile_C_form) =='hamiltonian') THEN 
      DO ie = 1, ne
         WRITE ( dos_unit, '(2(f15.9))' ) egrid(ie), SUM( dos(ie,:) )
      ENDDO
   ELSE
      CALL errore('conductor','problem in datafile_C_form', 10 )
   ENDIF
   !
   CLOSE( dos_unit )


!! ...  Finalize timing
!
   CALL timing('conductor',OPR='stop')
   CALL timing_overview(stdout,LIST=global_list,MAIN_NAME='conductor')
!
!...  free memory
!
   DEALLOCATE ( conduct, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating conduct', ABS(ierr) )
        !
   DEALLOCATE ( dos, STAT=ierr )
        IF( ierr /=0 ) CALL errore('conductor','deallocating dos', ABS(ierr) )
        !
   !
   CALL cleanup()

END PROGRAM conductor_qudot_imp_gamma
  
