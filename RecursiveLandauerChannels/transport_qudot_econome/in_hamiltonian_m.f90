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
   MODULE in_hamiltonian_module
!*********************************************
   USE parameters,            ONLY : nstrx
   USE kinds
   USE constants,             ONLY : ZERO, CZERO, CONE, ONE, EPS_m6, EPS_m4, EPS_m2
   USE io_global_module,      ONLY : stdout
   USE io_module,             ONLY :  mat_unit => aux_unit, ioname
   USE files_module,          ONLY : file_open, file_close
   USE iotk_module
   USE T_control_module,      ONLY : in_datafile_C,     &
                                     dim_subspace,      &
                                     dimC
   USE T_egrid_module,        ONLY : ne, egrid
   USE util_module,           ONLY : mat_mul
!   USE timing_module,         ONLY : timing

   IMPLICIT NONE
   PRIVATE 
   SAVE

   COMPLEX(dbl), ALLOCATABLE ::  h_CC(:,:)

!   COMPLEX(dbl), ALLOCATABLE ::  g_CC(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  h_LC(:,:), h_CR(:,:)

!   COMPLEX(dbl), ALLOCATABLE ::  gamma_tilde_L(:,:,:), gamma_tilde_R(:,:,:)

!   COMPLEX(dbl), ALLOCATABLE ::  sigma_tilde_L(:,:,:), sigma_tilde_R(:,:,:)

   !
   LOGICAL :: alloc
   !

   !
   PUBLIC :: h_CC, h_LC, h_CR
!   PUBLIC :: g_CC
   !
!   PUBLIC :: gamma_tilde_L, gamma_tilde_R
!   PUBLIC :: sigma_tilde_L, sigma_tilde_R

   PUBLIC :: alloc
    !
   ! general functions
   PUBLIC :: hamiltonian_allocate
   PUBLIC :: hamiltonian_deallocate
   PUBLIC :: read_hamiltonian
!    PUBLIC :: h_calcul_gamma
!    PUBLIC :: h_calcul_sigma
!    PUBLIC :: h_calcul_g_CC
   PUBLIC :: h_calcul_transmittance
    !
    !



CONTAINS



!*******************************************************************
   SUBROUTINE hamiltonian_allocate
   !*******************************************************************
      CHARACTER(20)      :: subname="hamiltonian_allocate"
      INTEGER :: ierr


    ! Allocate A-type matrices
   ALLOCATE ( h_CC(dimC,dimC), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating in_h_CC', ABS(ierr) )
   ALLOCATE ( h_LC(dim_subspace,dimC), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating h_LC', ABS(ierr) )
   ALLOCATE ( h_CR(dimC,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating h_CR', ABS(ierr) )

!    ALLOCATE ( g_CC(ne,dimC,dimC), STAT=ierr )
!         IF( ierr /=0 ) CALL errore(subname, 'allocating in_h_CC', ABS(ierr) )


!    ALLOCATE ( gamma_tilde_R(ne,dimC,dimC), STAT=ierr )
!         IF( ierr /=0 ) CALL errore(subname, 'allocating gamma_tilde_R', ABS(ierr) )
!    ALLOCATE ( gamma_tilde_L(ne,dimC,dimC), STAT=ierr )
!         IF( ierr /=0 ) CALL errore(subname, 'allocating gamma_tilde_L', ABS(ierr) )
! 
!    ALLOCATE ( sigma_tilde_R(ne,dimC,dimC), STAT=ierr )
!         IF( ierr /=0 ) CALL errore(subname, 'allocating sigma_tilde_R', ABS(ierr) )
!    ALLOCATE (  sigma_tilde_L(ne,dimC,dimC), STAT=ierr )
!         IF( ierr /=0 ) CALL errore(subname, 'allocating sigma_tilde_L', ABS(ierr) )

   h_CC(:,:) = CZERO
   h_LC(:,:) = CZERO
   h_CR(:,:) = CZERO

!    g_CC(:,:,:) = CZERO

!    gamma_tilde_R(:,:,:) = CZERO
!    gamma_tilde_L(:,:,:) = CZERO
!    sigma_tilde_R(:,:,:) = CZERO
!    sigma_tilde_L(:,:,:) = CZERO

   ! set alloc to true
    alloc=.TRUE.

END SUBROUTINE hamiltonian_allocate

!**********************************************************
   SUBROUTINE hamiltonian_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(22)      :: subname="hamiltonian_deallocate"
      INTEGER :: ierr


     IF ( ALLOCATED( h_CC  ) ) THEN
           DEALLOCATE ( h_CC, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating h_CC', ABS(ierr) )
      ENDIF
     !
     IF ( ALLOCATED( h_CR  ) ) THEN
           DEALLOCATE ( h_CR, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating h_CR', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( h_LC  ) ) THEN
           DEALLOCATE ( h_LC, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating h_LC', ABS(ierr) )
      ENDIF

!      IF ( ALLOCATED( g_CC  ) ) THEN
!            DEALLOCATE ( g_CC, STAT=ierr )
!            IF( ierr /=0 ) CALL errore(subname, 'deallocating g_CC', ABS(ierr) )
!       ENDIF
!      !
!      IF ( ALLOCATED( gamma_tilde_R ) ) THEN
!            DEALLOCATE ( gamma_tilde_R, STAT=ierr )
!            IF( ierr /=0 ) CALL errore(subname, 'deallocating gamma_tilde_R', ABS(ierr) )
!       ENDIF
! 
!      IF ( ALLOCATED(gamma_tilde_L  ) ) THEN
!            DEALLOCATE (gamma_tilde_L , STAT=ierr )
!            IF( ierr /=0 ) CALL errore(subname, 'deallocating gamma_tilde_L', ABS(ierr) )
!       ENDIF
! 
!      IF ( ALLOCATED( sigma_tilde_R ) ) THEN
!            DEALLOCATE ( sigma_tilde_R, STAT=ierr )
!            IF( ierr /=0 ) CALL errore(subname, 'deallocating sigma_tilde_R', ABS(ierr) )
!       ENDIF
! 
!      IF ( ALLOCATED(sigma_tilde_L  ) ) THEN
!            DEALLOCATE (sigma_tilde_L , STAT=ierr )
!            IF( ierr /=0 ) CALL errore(subname, 'deallocating sigma_tilde_L', ABS(ierr) )
!       ENDIF



     alloc=.FALSE.

   END SUBROUTINE hamiltonian_deallocate


!*******************************************************************
   SUBROUTINE read_hamiltonian
   !*******************************************************************
       CHARACTER(11)      :: subname="read_matrix"
       CHARACTER(6) :: name="MATRIX"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: iter
       INTEGER            :: ierr
       INTEGER         :: final_iter, dim_sub
      CHARACTER( LEN=nstrx )  :: filename

       CALL file_open(mat_unit,TRIM(in_datafile_C),PATH="/MATRIX/",ACTION="read", &
                              FORM='formatted')


 !      CALL iotk_scan_empty(mat_unit,"DATADIM",ATTR=attr)
       !CALL iotk_scan_attr(attr,"dimwan_total",dimwan_total,FIRST=.TRUE.) 
       !CALL iotk_scan_attr(attr,"dim_rec",dim_rec) 
!        CALL iotk_scan_attr(attr,"dim_sub",dim_sub) 
!        IF (dim_sub /= dim_subspace) CALL errore(subname, 'dim_subspace in input file and in .mat file are /=', ABS(dim_sub) )
! 
! 
!        CALL iotk_scan_empty(mat_unit,"DATAFINAL",ATTR=attr)
!        !CALL iotk_scan_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
!        CALL iotk_scan_attr(attr,"final_iter",final_iter) 
!        !CALL iotk_scan_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
!        IF (final_iter+1 /= in_max_iter_C) CALL errore(subname, 'read n_iter_C and given n_iter_C are /=', ABS(final_iter+1) )
!        !
! 
        CALL iotk_scan_begin(mat_unit,"H_MAT")
        CALL iotk_scan_dat(mat_unit,"H_LC", h_LC(:,:))
        CALL iotk_scan_dat(mat_unit,"H_CC", h_CC(:,:))
        CALL iotk_scan_dat(mat_unit,"H_CR", h_CR(:,:))
        CALL iotk_scan_end(mat_unit,"H_MAT")
!        !
!        !
! 
        ! CALL iotk_scan_end(mat_unit,TRIM(name))

      CALL file_close(mat_unit,PATH="/MATRIX/",ACTION="read")

      !CALL ioname('matrix',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'Hamiltonian Matrices read from file : ',5x,a)") TRIM(in_datafile_C)
      ! End of right part

  END SUBROUTINE read_hamiltonian


   !*******************************************************************
   SUBROUTINE h_calcul_transmittance(gamma_L_bare, gamma_R_bare, gL, gR, dos, ene, conduct)
   !*******************************************************************
   !
   ! input/output variables
   !
   REAL(dbl),    INTENT(out) :: dos(dimC)
   REAL(dbl),    INTENT(out) ::  conduct(dimC)

   COMPLEX(dbl), INTENT(in) :: gamma_L_bare(dim_subspace,dim_subspace)
   COMPLEX(dbl), INTENT(in) :: gamma_R_bare(dim_subspace,dim_subspace)

   COMPLEX(dbl), INTENT(in)  :: gL(dim_subspace,dim_subspace)
   COMPLEX(dbl), INTENT(in)  :: gR(dim_subspace,dim_subspace)

   COMPLEX(dbl), INTENT(in) :: ene
   !
   ! local variables
   !

   COMPLEX(dbl) :: h_aux1(dimC,dim_subspace), work1(dimC,dim_subspace), h_aux2(dim_subspace,dimC)


   COMPLEX(dbl) :: tmp(dimC,dimC), tmp1(dimC,dimC)
   COMPLEX(dbl) :: ga_aux(dimC,dimC)
   REAL(dbl)    :: dos_aux(dimC)

   COMPLEX(dbl) :: gamma_tilde_L(dimC,dimC)
   COMPLEX(dbl) :: gamma_tilde_R(dimC,dimC)

   COMPLEX(dbl) :: sigma_tilde_L(dimC,dimC)
   COMPLEX(dbl) :: sigma_tilde_R(dimC,dimC)

   COMPLEX(dbl) :: g_CC(dimC,dimC)
   INTEGER :: i, j, ierr
   !
   ! end of declarations
   !
     !!!!!!!!!!!!!!!!!!!!GAMMA
      gamma_tilde_L(:,:) = CZERO
      gamma_tilde_R(:,:) = CZERO

      !PRINT*, 'Calcul Gamma 1'

      h_aux1(:,:) =  CONJG( TRANSPOSE(h_LC(:,:)) )
      CALL mat_mul(work1, h_aux1 , 'N', gamma_L_bare(:,:), 'N', dimC, dim_subspace, dim_subspace)
      CALL mat_mul(gamma_tilde_L(:,:), work1, 'N', h_LC(:,:), 'N', dimC, dimC, dim_subspace)

      work1(:,:) = CZERO

      !PRINT*, 'Calcul Gamma 2'
      h_aux2(:,:) =  CONJG( TRANSPOSE(h_CR(:,:)) )
      CALL mat_mul(work1, h_CR(:,:)  , 'N', gamma_R_bare(:,:), 'N', dimC, dim_subspace, dim_subspace)
      CALL mat_mul(gamma_tilde_R(:,:), work1 , 'N', h_aux2 , 'N', dimC, dimC, dim_subspace)

     !!!!!!!!!!!!!!!!!!!!!!!SIGMA

      sigma_tilde_L(:,:) = CZERO
      sigma_tilde_R(:,:) = CZERO
      !PRINT*, 'Calcul Sigma 1'

      work1(:,:) = CZERO
      h_aux1(:,:) =  CONJG( TRANSPOSE(h_LC(:,:)) )
      CALL mat_mul(work1, h_aux1, 'N', gL, 'N', dimC, dim_subspace, dim_subspace)
      CALL mat_mul(sigma_tilde_L(:,:), work1, 'N', h_LC(:,:), 'N', dimC, dimC, dim_subspace)
      !

      work1(:,:) = CZERO
      !PRINT*, 'Calcul Sigma 2'
      h_aux2(:,:) =  CONJG( TRANSPOSE(h_CR(:,:)) )
      CALL mat_mul(work1, h_CR, 'N', gR, 'N', dimC, dim_subspace, dim_subspace)
      CALL mat_mul(sigma_tilde_R(:,:), work1, 'N', h_aux2(:,:), 'N', dimC, dimC, dim_subspace)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!GCC

     g_CC(:,:) = CZERO
     dos_aux(:) = ZERO
     CALL green_CC_inv( g_CC(:,:), dos_aux(:), ene, h_CC(:,:), sigma_tilde_L(:,:), sigma_tilde_R(:,:), dimC) 

     dos(:) = dos_aux(:)

!
!------------------------------
! main body
!------------------------------
!
   !CALL timing('transmittance_min', OPR='start')
   ! 
   ! adding the identity matrix
   ! 
   ! gL * gr -> tmp1
   !
   CALL mat_mul(tmp1, gamma_tilde_L(:,:), 'N', g_CC(:,:), 'N', dimC, dimC, dimC)
   !
   ! gL * gr * gR -> tmp
   !
   CALL mat_mul(tmp, tmp1, 'N', gamma_tilde_R(:,:), 'N', dimC, dimC, dimC)
   !
   ! gL * gr * gR * ga -> tmp
   !
   ga_aux(:,:) =  CONJG( TRANSPOSE( g_CC(:,:) ))
   CALL mat_mul(tmp1, tmp, 'N', ga_aux, 'N', dimC, dimC, dimC)
   !
   !
   DO i=1,dimC
      conduct(i) = REAL( tmp1(i,i) )
   ENDDO

END SUBROUTINE h_calcul_transmittance


 END MODULE in_hamiltonian_module


