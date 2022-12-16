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
   MODULE in_matrix_module
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
                                     in_max_iter_C
   USE T_egrid_module,        ONLY : ne, egrid
   USE util_module,           ONLY : mat_mul
!   USE timing_module,         ONLY : timing

   IMPLICIT NONE
   PRIVATE 
   SAVE

   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_C(:,:,:), B_matrix_C(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  B_matrix_L_C(:,:,:), B_matrix_C_R(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  gamma_tilde_L(:,:,:), gamma_tilde_R(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  sigma_tilde_L(:,:,:), sigma_tilde_R(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  P1_g_CC_PN(:,:,:), PN_g_CC_P1(:,:,:)
   !
   LOGICAL :: alloc
   !

   !
   PUBLIC :: A_matrix_C, B_matrix_C, B_matrix_L_C, B_matrix_C_R
   PUBLIC :: P1_g_CC_PN, PN_g_CC_P1
   PUBLIC :: gamma_tilde_L, gamma_tilde_R
   PUBLIC :: sigma_tilde_L, sigma_tilde_R

   PUBLIC :: alloc
    !
   ! general functions
   PUBLIC :: in_matrix_allocate
   PUBLIC :: in_matrix_deallocate
   PUBLIC :: read_matrix
   PUBLIC :: m_calcul_gamma
   PUBLIC :: m_calcul_sigma
   PUBLIC :: m_calcul_g_CC
   PUBLIC :: m_calcul_transmittance
    !
    !



CONTAINS



!*******************************************************************
   SUBROUTINE in_matrix_allocate
   !*******************************************************************
      CHARACTER(18)      :: subname="in_matrix_allocate"
      INTEGER :: ierr


    ! Allocate A-type matrices
   ALLOCATE ( A_matrix_C(in_max_iter_C,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating in_A_matrix_C', ABS(ierr) )
   ALLOCATE ( B_matrix_C((in_max_iter_C-1),dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating in_B_matrix_C', ABS(ierr) )
   ALLOCATE ( B_matrix_L_C(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating in_B_matrix_L_C', ABS(ierr) )
   ALLOCATE ( B_matrix_C_R(1,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating in_B_matrix_C_R', ABS(ierr) )

   ALLOCATE ( P1_g_CC_PN(ne,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating P1_g_CC_PN', ABS(ierr) )
   ALLOCATE ( PN_g_CC_P1(ne,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating PN_g_CC_P1', ABS(ierr) )

   ALLOCATE ( gamma_tilde_R(ne,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating gamma_tilde_R', ABS(ierr) )
   ALLOCATE ( gamma_tilde_L(ne,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating gamma_tilde_L', ABS(ierr) )
   ALLOCATE ( sigma_tilde_R(ne,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating sigma_tilde_R', ABS(ierr) )
   ALLOCATE ( sigma_tilde_L(ne,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating sigma_tilde_L', ABS(ierr) )


   ! set alloc to true
    alloc=.TRUE.




END SUBROUTINE in_matrix_allocate

!**********************************************************
   SUBROUTINE in_matrix_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(20)      :: subname="in_matrix_deallocate"
      INTEGER :: ierr


     IF ( ALLOCATED( A_matrix_C  ) ) THEN
           DEALLOCATE ( A_matrix_C, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_C', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_C  ) ) THEN
           DEALLOCATE ( B_matrix_C, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_C', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_L_C  ) ) THEN
           DEALLOCATE ( B_matrix_L_C, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_L_C', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_C_R  ) ) THEN
           DEALLOCATE ( B_matrix_C_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_C_R', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( P1_g_CC_PN  ) ) THEN
           DEALLOCATE ( P1_g_CC_PN, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating P1_g_CC_PN', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( PN_g_CC_P1  ) ) THEN
           DEALLOCATE ( PN_g_CC_P1, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating PN_g_CC_P1', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( gamma_tilde_R ) ) THEN
           DEALLOCATE ( gamma_tilde_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating gamma_tilde_R', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED(gamma_tilde_L  ) ) THEN
           DEALLOCATE (gamma_tilde_L , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating gamma_tilde_L', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( sigma_tilde_R ) ) THEN
           DEALLOCATE ( sigma_tilde_R, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating sigma_tilde_R', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED(sigma_tilde_L  ) ) THEN
           DEALLOCATE (sigma_tilde_L , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating sigma_tilde_L', ABS(ierr) )
      ENDIF




     alloc=.FALSE.

   END SUBROUTINE in_matrix_deallocate


!*******************************************************************
   SUBROUTINE read_matrix
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


       CALL iotk_scan_empty(mat_unit,"DATADIM",ATTR=attr)
       !CALL iotk_scan_attr(attr,"dimwan_total",dimwan_total,FIRST=.TRUE.) 
       !CALL iotk_scan_attr(attr,"dim_rec",dim_rec) 
       CALL iotk_scan_attr(attr,"dim_sub",dim_sub) 
       IF (dim_sub /= dim_subspace) CALL errore(subname, 'dim_subspace in input file and in .mat file are /=', ABS(dim_sub) )


       CALL iotk_scan_empty(mat_unit,"DATAFINAL",ATTR=attr)
       !CALL iotk_scan_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
       CALL iotk_scan_attr(attr,"final_iter",final_iter) 
       !CALL iotk_scan_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
       IF (final_iter+1 /= in_max_iter_C) CALL errore(subname, 'read n_iter_C and given n_iter_C are /=', ABS(final_iter+1) )
       !

       CALL iotk_scan_begin(mat_unit,"A_MAT")
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), A_matrix_C(1,:,:))
       DO iter=1, final_iter
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), A_matrix_C(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"A_MAT")
       !
       !
       CALL iotk_scan_begin(mat_unit,"B_MAT")
       !!! !!!!!!!! ATTENTION FINAL ITER non lue pour B_matrix_C
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), B_matrix_C(1,:,:))
       DO iter=1, (final_iter-1)
        CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), B_matrix_C(iter+1,:,:))
       ENDDO
       CALL iotk_scan_end(mat_unit,"B_MAT")

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       !!!!!!!! ATTENTION FINAL ITER = B_matrix_LC
       !!!!!!!! ATTENTION FINAL ITER +1  = B_matrix_CR
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(final_iter)), B_matrix_L_C(1,:,:))
       CALL iotk_scan_dat(mat_unit,"ITER"//TRIM(iotk_index(final_iter+1)), B_matrix_C_R(1,:,:))



      ! CALL iotk_scan_end(mat_unit,TRIM(name))

      CALL file_close(mat_unit,PATH="/MATRIX/",ACTION="read")

      !CALL ioname('matrix',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Matrices on recursion basis read from file : ',5x,a)") TRIM(in_datafile_C)
      ! End of right part

  END SUBROUTINE read_matrix

!*******************************************************************
   SUBROUTINE m_calcul_gamma( gamma_L_bare, gamma_R_bare, iene  )
   !*******************************************************************
      COMPLEX(dbl), INTENT(in) :: gamma_L_bare(dim_subspace,dim_subspace)
      COMPLEX(dbl), INTENT(in) :: gamma_R_bare(dim_subspace,dim_subspace)
      INTEGER, INTENT(in)      :: iene
      CHARACTER(14)      :: subname="m_calcul_gamma"
      INTEGER :: ierr
      COMPLEX(dbl), ALLOCATABLE :: B_aux(:,:), work(:,:)

      ALLOCATE ( B_aux(dim_subspace,dim_subspace), STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'allocating B_aux', ABS(ierr) )
      ALLOCATE ( work(dim_subspace,dim_subspace), STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'allocating work', ABS(ierr) )

      gamma_tilde_R(iene,:,:) = CZERO
      gamma_tilde_L(iene,:,:) = CZERO

      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_L_C(1,:,:)) )
      CALL mat_mul(work, B_matrix_L_C(1,:,:), 'N', gamma_L_bare(:,:), 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(gamma_tilde_L(iene,:,:), work, 'N', B_aux, 'N', dim_subspace, dim_subspace, dim_subspace)

      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_C_R(1,:,:)) )
      CALL mat_mul(work, B_aux , 'N', gamma_R_bare(:,:), 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(gamma_tilde_R(iene,:,:), work, 'N', B_matrix_C_R(1,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)

      DEALLOCATE ( B_aux, STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'deallocating B_aux', ABS(ierr) )
      DEALLOCATE ( work, STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'deallocating work', ABS(ierr) )


END SUBROUTINE m_calcul_gamma

!*******************************************************************
   SUBROUTINE m_calcul_sigma( gL, gR, iene  )
   !*******************************************************************
      COMPLEX(dbl), INTENT(in) :: gL(dim_subspace,dim_subspace)
      COMPLEX(dbl), INTENT(in) :: gR(dim_subspace,dim_subspace)
      INTEGER, INTENT(in)      :: iene
      CHARACTER(14)      :: subname="m_calcul_gamma"
      INTEGER :: ierr
      COMPLEX(dbl), ALLOCATABLE :: B_aux(:,:), work(:,:)

      ALLOCATE ( B_aux(dim_subspace,dim_subspace), STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'allocating B_aux', ABS(ierr) )
      ALLOCATE ( work(dim_subspace,dim_subspace), STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'allocating work', ABS(ierr) )

      sigma_tilde_R(iene,:,:) = CZERO
      sigma_tilde_L(iene,:,:) = CZERO

      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_L_C(1,:,:)) )
      CALL mat_mul(work, B_matrix_L_C(1,:,:), 'N', gL, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(sigma_tilde_L(iene,:,:), work, 'N', B_aux, 'N', dim_subspace, dim_subspace, dim_subspace)
      !
      B_aux(:,:) =  CONJG( TRANSPOSE(B_matrix_C_R(1,:,:)) )
      CALL mat_mul(work, B_aux, 'N', gR, 'N', dim_subspace, dim_subspace, dim_subspace)
      CALL mat_mul(sigma_tilde_R(iene,:,:), work, 'N', B_matrix_C_R(1,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)

      DEALLOCATE ( B_aux, STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'deallocating B_aux', ABS(ierr) )
      DEALLOCATE ( work, STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname,'deallocating work', ABS(ierr) )


END SUBROUTINE m_calcul_sigma

!*******************************************************************
   SUBROUTINE m_calcul_g_CC( dos, ene, iene)
   !*******************************************************************
      COMPLEX(dbl), INTENT(in) :: ene
      INTEGER, INTENT(in)      :: iene
      REAL(dbl), INTENT(out) :: dos(dim_subspace)
      CHARACTER(13)      :: subname="m_calcul_g_CC"
      COMPLEX(dbl) :: dos_aux(dim_subspace)
      INTEGER :: ierr

     CALL green_CC_tridiag( P1_g_CC_PN(iene,:,:), PN_g_CC_P1(iene, :,:), dos_aux(:), ene, A_matrix_C(:,:,:), B_matrix_C(:,:,:), sigma_tilde_L(iene,:,:), sigma_tilde_R(iene,:,:), in_max_iter_C, dim_subspace) 

     dos(:) = dos_aux(:)

END SUBROUTINE m_calcul_g_CC

   !*******************************************************************
   SUBROUTINE m_calcul_transmittance(conduct, iene)
   !*******************************************************************
   !
   ! input/output variables
   !
   REAL(dbl),    INTENT(out)::  conduct(dim_subspace)
   INTEGER, INTENT(in)      :: iene

   !
   ! local variables
   !
   COMPLEX(dbl) :: tmp(dim_subspace,dim_subspace), tmp1(dim_subspace,dim_subspace)
   COMPLEX(dbl) :: ga_aux(dim_subspace,dim_subspace)
   INTEGER :: i, j, ierr
   !
   ! end of declarations
   !

!
!------------------------------
! main body
!------------------------------
!
!   CALL timing('transmittance_min', OPR='start')
   ! 
   ! adding the identity matrix
   ! 
   ! gL * gr -> tmp1
   !
   CALL mat_mul(tmp1, gamma_tilde_L(iene,:,:), 'N', P1_g_CC_PN(iene,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)
   !
   ! gL * gr * gR -> tmp
   !
   CALL mat_mul(tmp, tmp1, 'N', gamma_tilde_R(iene,:,:), 'N', dim_subspace, dim_subspace, dim_subspace)
   !
   ! gL * gr * gR * ga -> tmp
   !
   ga_aux(:,:) =  CONJG( TRANSPOSE( P1_g_CC_PN(iene,:,:) ))
   CALL mat_mul(tmp1, tmp, 'N', ga_aux, 'N', dim_subspace, dim_subspace, dim_subspace)
   !
   !
   DO i=1,dim_subspace
      conduct(i) = REAL( tmp1(i,i) )
   ENDDO

END SUBROUTINE m_calcul_transmittance


 END MODULE in_matrix_module


