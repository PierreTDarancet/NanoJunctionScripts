!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE recursion_function_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE dim_variable_module,        ONLY : dim_subspace, dim_recursion
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE,  EPS_m4, EPS_m7, EPS_m10, EPS_m14
   USE iotk_module
   USE util_module,          ONLY : mat_hdiag, mat_mul, mat_sv


   IMPLICIT NONE
   PRIVATE 
   SAVE

   ! general functions

    PUBLIC :: calculate_variation
    PUBLIC :: calculate_recursion_state
    PUBLIC :: calculate_A_matrix
    PUBLIC :: normalize_recursion_state
    PUBLIC :: normalize_recursion_state_rpp
    PUBLIC :: calculate_total_variation
    PUBLIC :: calculate_recursion_state_diag_A



CONTAINS

 !*******************************************************************
    SUBROUTINE calculate_variation(A_matrix_new, A_matrix, variation, debug_m)
    !*******************************************************************
   !  Input variables
        COMPLEX(dbl), INTENT(in) :: A_matrix_new(dim_subspace,dim_subspace)
        COMPLEX(dbl), INTENT(in) :: A_matrix(dim_subspace,dim_subspace)
        LOGICAL,      INTENT(in) :: debug_m
   ! output variables
       REAL(dbl), INTENT(out) :: variation(dim_subspace)
    ! local variables
    CHARACTER(19) :: subname="calculate_variation"
    INTEGER       :: ierr, i_sub
    REAL(dbl) :: eig_new(dim_subspace), eig(dim_subspace)
    COMPLEX(dbl) :: z(dim_subspace,dim_subspace)

      ! init
     eig_new(:)=ZERO
     eig(:)=ZERO
      !

    ! CALL timing('mat_hdiag', OPR='start')
   IF (dim_subspace /= 1) THEN
      !
     CALL mat_hdiag(z(:,:) , eig_new(:), A_matrix_new(:,:), dim_subspace )
      !

     CALL mat_hdiag(z(:,:) , eig(:), A_matrix(:,:), dim_subspace )

   ELSE IF (dim_subspace == 1) THEN 
     eig_new(1)=A_matrix_new(1,1)
     eig(1)=A_matrix(1,1)
   ELSE
      CALL errore(subname,'problem in dim_subspace',ABS(dim_subspace))

   ENDIF
     !CALL timing('mat_hdiag', OPR='stop')

     DO i_sub=1, dim_subspace
        !
         variation(i_sub)= ABS(eig_new(i_sub) - eig(i_sub))
        !
     ENDDO



    END SUBROUTINE calculate_variation


 !*******************************************************************
    SUBROUTINE calculate_total_variation(A_matrix_new, A_matrix, variation, debug_m)
    !*******************************************************************
   !  Input variables
        COMPLEX(dbl), INTENT(in) :: A_matrix_new(dim_subspace,dim_subspace)
        COMPLEX(dbl), INTENT(in) :: A_matrix(dim_subspace,dim_subspace)
        LOGICAL,      INTENT(in) :: debug_m
   ! output variables
       REAL(dbl), INTENT(out) :: variation(dim_subspace)
    ! local variables
    CHARACTER(65) :: subname="calculate_total_variation"
    INTEGER       :: ierr, i_sub, j_sub
    REAL(dbl) :: row_new(dim_subspace), row(dim_subspace)
    COMPLEX(dbl) :: z(dim_subspace,dim_subspace)

      ! init
     row_new(:)=ZERO
     row(:)=ZERO
      !

    ! CALL timing('mat_hdiag', OPR='start')
   IF (dim_subspace /= 1) THEN
      !
      DO j_sub=1, dim_subspace
         row_new(i_sub) = row_new(i_sub) + A_matrix_new(i_sub,j_sub)
         !
         row(i_sub) = row(i_sub) + A_matrix(i_sub,j_sub)
         !
      ENDDO
   ELSE IF (dim_subspace == 1) THEN 
     row_new(1)=A_matrix_new(1,1)
     row(1)=A_matrix(1,1)
   ELSE
      CALL errore(subname,'problem in dim_subspace',ABS(dim_subspace))

   ENDIF
     !CALL timing('mat_hdiag', OPR='stop')

     DO i_sub=1, dim_subspace
        !
         variation(i_sub)= ABS(row_new(i_sub) - row(i_sub))
        !
     ENDDO



    END SUBROUTINE calculate_total_variation


 !*******************************************************************
    SUBROUTINE calculate_recursion_state(Hamiltonian, state, A_matrix, B_matrix, state_old, state_new, debug_m )
    !*******************************************************************
  ! input variables
      COMPLEX(dbl), INTENT(in) :: Hamiltonian(dim_recursion,dim_recursion)
      COMPLEX(dbl), INTENT(in) :: state(dim_subspace, dim_recursion)
      COMPLEX(dbl), INTENT(in) :: A_matrix(dim_subspace,dim_subspace)
      COMPLEX(dbl), INTENT(in) :: B_matrix(dim_subspace,dim_subspace)
      COMPLEX(dbl), INTENT(in) :: state_old(dim_subspace, dim_recursion)
      LOGICAL,      INTENT(in) :: debug_m
   ! output variables
       !
      COMPLEX(dbl), INTENT(out) :: state_new(dim_subspace, dim_recursion)
    !
   ! local variables
    CHARACTER(25) :: subname="calculate_recursion_state"
    INTEGER       :: ierr, i_recur, i_sub, j_sub
    INTEGER       :: dim_vect=1
    COMPLEX(dbl)  :: state_new_aux(dim_subspace, dim_recursion)
    COMPLEX(dbl)  :: state_new_aux2(dim_subspace, dim_recursion)
    REAL(dbl)   :: test_orthogo_aux, test_orthogo

       ! peut ?tre optimis? : d?j? calcul? au step pr?c?dent (calcul A_matrix)
  DO i_sub=1, dim_subspace

     !CALL mat_mul(state_new_aux(i_sub,:) , Hamiltonian, 'N', state(i_sub,:) , 'N', dim_recursion, dim_vect, dim_recursion)
      !    CALL mat_mul(A, aux_CL, 'N', gL, 'N', dimC, dimL, dimL)
     state_new_aux(i_sub,:) = matmul(Hamiltonian, state(i_sub,:))
  ENDDO
! 

  DO i_sub=1, dim_subspace
        state_new_aux2(i_sub, : )  = state_new_aux(i_sub, :) 
        DO i_recur=1, dim_recursion
         DO  j_sub=1, dim_subspace
            state_new_aux2(i_sub, i_recur )  = state_new_aux2(i_sub,i_recur) - A_matrix(j_sub,i_sub) *state(j_sub,i_recur)
            !state_new_aux2(i_sub, i_recur )  = state_new_aux(i_sub, i_recur) - A_matrix(j_sub,i_sub) *state(i_sub,i_recur)

         ENDDO
        ENDDO
   ENDDO


  DO i_sub=1, dim_subspace

        state_new(i_sub, : )  = state_new_aux2(i_sub, :) 
        DO i_recur=1, dim_recursion
         DO  j_sub=1, dim_subspace
            !state_new_aux2(i_sub, i_recur )  = state_new_aux(i_sub, i_recur) - A_matrix(j_sub,i_sub) *state(i_sub,i_recur)
            state_new(i_sub, i_recur)  = state_new(i_sub,i_recur)    -  CONJG( B_matrix(i_sub,j_sub)) * state_old(j_sub, i_recur)
         ENDDO
        ENDDO
   ENDDO

    test_orthogo_aux=ZERO
    test_orthogo=ZERO
    ! test the new states
    DO i_sub=1, dim_subspace
         DO j_sub=1, dim_subspace
            test_orthogo_aux = dot_product(state_new(i_sub,:) , state_old(j_sub,:))
            test_orthogo = test_orthogo + REAL(ABS(test_orthogo_aux), dbl)
         ENDDO
    ENDDO

    IF ( (test_orthogo >= EPS_m4).OR.( test_orthogo <= - EPS_m4) ) THEN
        IF (debug_m) THEN
            CALL errore(subname,'new vectors not orthogonal to n-1 states',-100.00)
        ELSE
            CALL errore(subname,'new vectors not orthogonal to n-1 states',1)
        ENDIF
    ENDIF

    test_orthogo_aux=ZERO
    test_orthogo=ZERO
    DO i_sub=1, dim_subspace
         DO j_sub=1, dim_subspace
            test_orthogo_aux = dot_product(state_new(i_sub,:) , state(j_sub,:))
            test_orthogo = test_orthogo + REAL(ABS(test_orthogo_aux), dbl)
         ENDDO
    ENDDO

    IF ( (test_orthogo >= EPS_m4).OR.( test_orthogo <= - EPS_m4) ) THEN
        IF (debug_m) THEN
            CALL errore(subname,'new vectors not orthogonal to n states',-100.00)
        ELSE
            CALL errore(subname,'new vectors not orthogonal to n states',2)

        ENDIF
    ENDIF



    END SUBROUTINE  calculate_recursion_state
 !*******************************************************************
    SUBROUTINE normalize_recursion_state( state, B_matrix, test_norm, debug_m)
    !*******************************************************************
        LOGICAL,      INTENT(in) :: debug_m 
   !input/output variables
      COMPLEX(dbl), INTENT(inout) :: state(dim_subspace, dim_recursion)
   ! output variables
      COMPLEX(dbl), INTENT(out) :: B_matrix(dim_subspace,dim_subspace)
      REAL(dbl), INTENT(out)  :: test_norm
   ! local variables
    CHARACTER(25) :: subname="normalize_recursion_state"
    COMPLEX(dbl) :: aux_state(dim_subspace, dim_recursion)
    COMPLEX(dbl) :: S_matrix(dim_subspace,dim_subspace)
    COMPLEX(dbl) :: cu(dim_subspace,dim_subspace)
    COMPLEX(dbl) :: norm
    REAL(dbl)   :: eig(dim_subspace)
    INTEGER       :: ierr, j_sub, i_sub, i_recur
    INTEGER       :: dim_vect=1
    REAL(dbl)   :: test_norm_aux, test_orthogo_aux, test_orthogo
    !
    !

   ! Initialize local variables
   B_matrix(:,:)= CZERO
   !
   aux_state(:,:)=state(:,:)
   !

   ! if dim subspace > 1 orthogonalization procedure else only normalization
   IF (dim_subspace /= 1 ) THEN
         !calculate S matrix
         ! scalar < state(i,:)|state(j,:)>
      DO j_sub=1, dim_subspace
         DO i_sub=1, dim_subspace
            !
            S_matrix(i_sub,j_sub) = (dot_product( state(i_sub,:) , state(j_sub,:) ))
            !
         ENDDO
      ENDDO

         ! diagonalize S matrix in order to find 2 orthogonal vectors

         CALL mat_hdiag(cu(:,:) , eig(:), S_matrix(:,:), dim_subspace )
         !
         ! calculate the  eigenvectors
         state(:,:)=CZERO
         DO i_sub=1, dim_subspace
            !
            DO j_sub=1, dim_subspace
                  ! State(i) = Sigma_j cu_transposed(i,j) * Aux_State(j)
                  state(i_sub,:)= state(i_sub,:) + ((CONJG(cu(j_sub, i_sub))) * aux_state(j_sub, :))
                  !
            ENDDO
            !
      !
            ! normalization
            norm = sqrt(dot_product( state(i_sub,:) , state(i_sub,:) ))
            state(i_sub, :)=   state(i_sub, :) / norm
      !
         ENDDO
      !    !
            !
      ELSE IF (dim_subspace == 1 ) THEN
            ! normalization
         DO i_sub=1, dim_subspace
            norm = sqrt(dot_product( state(i_sub,:) , state(i_sub,:) ))
            state(i_sub, :)=   state(i_sub, :) / norm
            !
         ENDDO

      ELSE
          CALL errore(subname,'problem in dim_subspace',ABS(dim_subspace))

      ENDIF
      !
      ! calculate B_matrix
      !
      DO i_sub=1, dim_subspace
         DO j_sub=1, dim_subspace
         !
         B_matrix(i_sub,j_sub) = (dot_product( state(i_sub,:) , aux_state(j_sub,:) ))
         !
         ENDDO
      ENDDO
      !

     ! numerical cut-off in order to avoid  diagonalization procedure's noise
      DO i_sub=1, dim_subspace
         DO i_recur=1, dim_recursion
            IF (ABS(state(i_sub, i_recur)) < EPS_m10 ) THEN
               state(i_sub, i_recur) = CZERO
            ENDIF
         ENDDO 
      ENDDO
      ! test normalisation :
   !
      test_norm = ZERO
   
      DO i_sub=1, dim_subspace
         !
         !CALL mat_mul( test_norm_aux , state(i_sub,:), 'C', state(i_sub,:) , 'N', dim_vect, dim_vect, dim_recursion)
         test_norm_aux = dot_product(state(i_sub,:) , state(i_sub,:))
         !
         test_norm = test_norm + REAL(ABS(test_norm_aux)/dim_subspace, dbl)
      ENDDO
      IF ( (test_norm >= ONE + EPS_m4).OR.( test_norm <= ONE - EPS_m4) ) THEN
         IF (debug_m) THEN
            CALL errore(subname,'problem in normalization ',-100)
         ELSE
            CALL errore(subname,'problem in normalization ',1)
         ENDIF
      ENDIF

     test_orthogo_aux = ZERO
     test_orthogo = ZERO

     DO i_sub=1, dim_subspace
         DO j_sub=1, dim_subspace
            test_orthogo_aux = dot_product(state(i_sub,:) , state(j_sub,:))
            test_orthogo = test_orthogo + REAL(ABS(test_orthogo_aux)/dim_subspace, dbl)
         ENDDO
    ENDDO


      IF ( (test_orthogo >= ONE + EPS_m4).OR.( test_orthogo <= ONE - EPS_m4) ) THEN
         IF (debug_m) THEN
            CALL errore(subname,'vector not orthogonal',-100)
         ELSE
            CALL errore(subname,'vector not orthogonal',2)
         ENDIF
      ENDIF

    END SUBROUTINE  normalize_recursion_state

 !*******************************************************************
    SUBROUTINE normalize_recursion_state_rpp( state, B_matrix, test_norm, debug_m)
    !*******************************************************************
      LOGICAL,      INTENT(in) :: debug_m
    !input/output variables
      COMPLEX(dbl), INTENT(inout) :: state(dim_subspace, dim_recursion)
   ! output variables
      COMPLEX(dbl), INTENT(out) :: B_matrix(dim_subspace,dim_subspace)
      REAL(dbl), INTENT(out)  :: test_norm
   ! local variables
    CHARACTER(29) :: subname="normalize_recursion_state_rpp"
    COMPLEX(dbl) :: aux_state(dim_subspace, dim_recursion)
    COMPLEX(dbl) :: S_matrix(dim_subspace,dim_subspace)
    COMPLEX(dbl) :: cu(dim_subspace,dim_subspace)
    COMPLEX(dbl) :: norm
    REAL(dbl)   :: eig(dim_subspace)
    REAL(dbl)    :: real_part
    INTEGER       :: ierr, j_sub, i_sub, i_recur
    INTEGER       :: dim_vect=1
    REAL(dbl)   :: test_norm_aux, test_orthogo_aux, test_orthogo
    !
    !

   ! Initialize local variables
   B_matrix(:,:)= CZERO
   !
   aux_state(:,:)=state(:,:)
   !

   ! if dim subspace > 1 orthogonalization procedure else only normalization
   IF (dim_subspace /= 1 ) THEN
         !calculate S matrix
         ! scalar < state(i,:)|state(j,:)>
      DO j_sub=1, dim_subspace
         DO i_sub=1, dim_subspace
            !
            S_matrix(i_sub,j_sub) = (dot_product( state(i_sub,:) , state(j_sub,:) ))
            !
         ENDDO
      ENDDO
    PRINT*, S_matrix
         ! diagonalize S matrix in order to find 2 orthogonal vectors

         CALL mat_hdiag(cu(:,:) , eig(:), S_matrix(:,:), dim_subspace )
         !
         ! calculate the  eigenvectors
         state(:,:)=CZERO
         DO i_sub=1, dim_subspace
            !
            DO j_sub=1, dim_subspace
                  ! State(i) = Sigma_j cu_transposed(i,j) * Aux_State(j)
                  state(i_sub,:)= state(i_sub,:) + ((CONJG(cu(j_sub, i_sub))) * aux_state(j_sub, :))
                  !
            ENDDO
            !
      !
            ! normalization
            norm = sqrt(dot_product( state(i_sub,:) , state(i_sub,:) ))
            state(i_sub, :)=   state(i_sub, :) / norm
      !
         ENDDO
      !    !
            !
      ELSE IF (dim_subspace == 1 ) THEN
            ! normalization
         DO i_sub=1, dim_subspace
            norm = sqrt(dot_product( state(i_sub,:) , state(i_sub,:) ))
            state(i_sub, :)=   state(i_sub, :) / norm
            !
         ENDDO

      ELSE
          CALL errore(subname,'problem in dim_subspace',ABS(dim_subspace))

      ENDIF


      ! calculate & change real part if necessary to avoid a phase  factor
      DO i_sub=1, dim_subspace
        real_part = ZERO

        DO i_recur=1, dim_recursion
            real_part =  real_part +  REAL(state(i_sub, i_recur), dbl)
        ENDDO

        IF  ( real_part < ZERO ) THEN
            state(i_sub, :) =   - state(i_sub, :)
        ENDIF

      ENDDO
      !

      !
      ! calculate B_matrix
      !
      DO i_sub=1, dim_subspace
         DO j_sub=1, dim_subspace
         !
         B_matrix(i_sub,j_sub) = (dot_product( state(i_sub,:) , aux_state(j_sub,:) ))
         !
         ENDDO
      ENDDO
      !

     ! numerical cut-off in order to avoid  diagonalization procedure's noise
      DO i_sub=1, dim_subspace
         DO i_recur=1, dim_recursion
            IF (ABS(state(i_sub, i_recur)) < EPS_m14 ) THEN
               state(i_sub, i_recur) = CZERO
            ENDIF
         ENDDO 
      ENDDO
      ! test normalisation :
   !
      test_norm = ZERO
   
      DO i_sub=1, dim_subspace
         !
         !CALL mat_mul( test_norm_aux , state(i_sub,:), 'C', state(i_sub,:) , 'N', dim_vect, dim_vect, dim_recursion)
         test_norm_aux = dot_product(state(i_sub,:) , state(i_sub,:))
         !
         test_norm = test_norm + REAL(ABS(test_norm_aux)/dim_subspace, dbl)
      ENDDO
   


      IF ( (test_norm >= ONE + EPS_m4).OR.( test_norm <= ONE - EPS_m4) ) THEN
         IF (debug_m) THEN
            CALL errore(subname,'problem in normalization ',-100)
         ELSE
            CALL errore(subname,'problem in normalization ',1)
         ENDIF
      ENDIF

     test_orthogo_aux = ZERO
     test_orthogo = ZERO

     DO i_sub=1, dim_subspace
         DO j_sub=1, dim_subspace
            test_orthogo_aux = dot_product(state(i_sub,:) , state(j_sub,:))
            test_orthogo = test_orthogo + REAL(ABS(test_orthogo_aux)/dim_subspace, dbl)
         ENDDO
    ENDDO


      IF ( (test_orthogo >= ONE + EPS_m4).OR.( test_orthogo <= ONE - EPS_m4) ) THEN
         IF (debug_m) THEN
            CALL errore(subname,'vectors not orthogonal',-100)
         ELSE
            CALL errore(subname,'vectors not orthogonal',2)
         ENDIF
      ENDIF


    END SUBROUTINE  normalize_recursion_state_rpp


 !*******************************************************************
    SUBROUTINE calculate_A_matrix( Hamiltonian , state, A_matrix, debug_m)
    !*******************************************************************
   !
      COMPLEX(dbl), INTENT(in) :: Hamiltonian(dim_recursion,dim_recursion)
      COMPLEX(dbl), INTENT(in) :: state(dim_subspace, dim_recursion)
      LOGICAL,      INTENT(in) :: debug_m
   ! output variables
      COMPLEX(dbl), INTENT(out) :: A_matrix(dim_subspace,dim_subspace)
!      COMPLEX(dbl), INTENT(out) :: state_new(dim_subspace, dim_recursion)
   ! local variables
      COMPLEX(dbl) :: state_new(dim_subspace, dim_recursion)
    CHARACTER(18) :: subname="calculate_A_matrix"
    INTEGER       :: ierr, i_sub, j_sub, i_recur
    INTEGER       :: dim_vect=1

  DO i_sub=1, dim_subspace
      !!
      !CALL mat_mul(state_new(i_sub,:),Hamiltonian, 'N', state(i_sub,:) , 'N', dim_recursion, dim_vect, dim_recursion)
      !    CALL mat_mul(A, aux_CL, 'N', gL, 'N', dimC, dimL, dimL)
     state_new(i_sub,:) = matmul(Hamiltonian , state(i_sub,:))
  ENDDO

    A_matrix(:,:) = CZERO

    ! ajouter hermitianit?
  DO j_sub=1, dim_subspace
   DO i_sub=1, dim_subspace
      !
      ! A_matrix( i_sub , j_sub ) =  A_matrix(i_sub , j_sub) + state_new(i_sub , i_recur) * state*(j_sub , i_recur)
      !CALL mat_mul(A_matrix(i_sub,j_sub), state(i_sub,:), 'C', state_new(j_sub,:) , 'N', dim_vect, dim_vect, dim_recursion)
      A_matrix(i_sub,j_sub) = dot_product(state(i_sub,:), state_new(j_sub,:))
      !
   ENDDO
  ENDDO

 END SUBROUTINE calculate_A_matrix




 !*******************************************************************
    SUBROUTINE calculate_recursion_state_diag_A(Hamiltonian, state, A_matrix, B_matrix_old, state_old, state_new, B_matrix, test_norm, debug_m )
    !*******************************************************************

  ! input variables
      COMPLEX(dbl), INTENT(in) :: Hamiltonian(dim_recursion,dim_recursion)
      COMPLEX(dbl), INTENT(in) :: state(dim_subspace, dim_recursion)
      COMPLEX(dbl), INTENT(in) :: A_matrix(dim_subspace,dim_subspace)
      COMPLEX(dbl), INTENT(in) :: B_matrix_old(dim_subspace,dim_subspace)
      COMPLEX(dbl), INTENT(in) :: state_old(dim_subspace, dim_recursion)
      LOGICAL,      INTENT(in) :: debug_m
    !input/output variables

   ! output variables
     COMPLEX(dbl), INTENT(out) :: B_matrix(dim_subspace,dim_subspace)
     COMPLEX(dbl), INTENT(out) :: state_new(dim_subspace, dim_recursion)
      REAL(dbl), INTENT(out)  :: test_norm
   ! local variables
    CHARACTER(32) :: subname="normalize_recursion_state_diag_A"
    COMPLEX(dbl) :: S_matrix(dim_subspace,dim_subspace)
    COMPLEX(dbl) :: cu(dim_subspace,dim_subspace)
    COMPLEX(dbl) :: state_new_aux(dim_subspace, dim_recursion)
    COMPLEX(dbl) :: state_new_aux2(dim_subspace, dim_recursion)
    COMPLEX(dbl) :: state_new_aux3(dim_subspace, dim_recursion)
    COMPLEX(dbl) :: A_matrix_aux(dim_subspace,dim_subspace)
    COMPLEX(dbl) :: norm
    REAL(dbl)   :: eig(dim_subspace)
    REAL(dbl)   :: test_norm_aux
    INTEGER       :: ierr, j_sub, i_sub, i_recur
    INTEGER       :: dim_vect=1



   ! Calculate recursion states
      !
      DO i_sub=1, dim_subspace
      !
         state_new_aux(i_sub,:) = matmul(Hamiltonian, state(i_sub,:))
      ENDDO
      ! 
      !
      DO i_sub=1, dim_subspace
            state_new_aux2(i_sub, : )  = state_new_aux(i_sub, :) 
            DO i_recur=1, dim_recursion
               DO  j_sub=1, dim_subspace
                  state_new_aux2(i_sub, i_recur )  = state_new_aux2(i_sub,i_recur) - A_matrix(j_sub,i_sub) *state(j_sub,i_recur)
               ENDDO
            ENDDO
         ENDDO
      !
      !
      DO i_sub=1, dim_subspace
            state_new(i_sub, : )  = state_new_aux2(i_sub, :) 
            DO i_recur=1, dim_recursion
               DO  j_sub=1, dim_subspace
                  state_new(i_sub, i_recur)  = state_new(i_sub,i_recur)    -  CONJG( B_matrix(i_sub,j_sub)) * state_old(j_sub, i_recur)
               ENDDO
            ENDDO
         ENDDO


   ! Calculate_S_matrix

   state_new_aux(:,:)=state_new(:,:)
   !

   ! if dim subspace > 1 orthogonalization procedure else only normalization
   IF (dim_subspace /= 1 ) THEN
         !calculate S matrix
         ! scalar < state_new(i,:)|state_new(j,:)>
      DO j_sub=1, dim_subspace
         DO i_sub=1, dim_subspace
            !
            S_matrix(i_sub,j_sub) = (dot_product( state_new(i_sub,:) , state_new(j_sub,:) ))
            !
         ENDDO
      ENDDO

         ! diagonalize S matrix in order to find 2 orthogonal vectors

         CALL mat_hdiag(cu(:,:) , eig(:), S_matrix(:,:), dim_subspace )
         !
         ! calculate the  eigenvectors
   ! calculate new recursion state

         state_new_aux2(:,:)=CZERO
         DO i_sub=1, dim_subspace
            !
            DO j_sub=1, dim_subspace
                  ! State(i) = Sigma_j cu_transposed(i,j) * Aux_State(j)
                  state_new_aux2(i_sub,:)= state_new_aux2(i_sub,:) + ((CONJG(cu(j_sub, i_sub))) * state_new_aux(j_sub, :))
                  !
            ENDDO
            !
      !
            ! normalization
            norm = sqrt(dot_product( state_new_aux2(i_sub,:) , state_new_aux2(i_sub,:) ))
            state_new_aux2(i_sub, :)=   state_new_aux2(i_sub, :) / norm
      !
         ENDDO
      !    !
   state_new_aux3(:,:)=CZERO

   ! Calculate A_matrix_aux
         DO i_sub=1, dim_subspace
               !!
               !CALL mat_mul(state_new(i_sub,:),Hamiltonian, 'N', state(i_sub,:) , 'N', dim_recursion, dim_vect, dim_recursion)
               !    CALL mat_mul(A, aux_CL, 'N', gL, 'N', dimC, dimL, dimL)
            state_new_aux3(i_sub,:) = matmul(Hamiltonian , state_new_aux2(i_sub,:))
         ENDDO
         !
            A_matrix_aux(:,:) = CZERO
         !
            ! ajouter hermitianit?
         DO j_sub=1, dim_subspace
            DO i_sub=1, dim_subspace
               !
               A_matrix_aux(i_sub,j_sub) = dot_product(state_new_aux2(i_sub,:), state_new_aux3(j_sub,:))
               !
            ENDDO
         ENDDO
 

   ! diagonalize A matrix

         CALL mat_hdiag(cu(:,:) , eig(:), A_matrix_aux(:,:), dim_subspace )
         !
      ! calculate final recursion state
         state_new(:,:)=CZERO
         DO i_sub=1, dim_subspace
            !
            DO j_sub=1, dim_subspace
                  ! State(i) = Sigma_j cu_transposed(i,j) * Aux_State(j)
                  state_new(i_sub,:)= state_new(i_sub,:) + ((CONJG(cu(j_sub, i_sub))) * state_new_aux2(j_sub, :))
                  !
            ENDDO
            !
        ENDDO



            !
      ELSE IF (dim_subspace == 1 ) THEN
            ! normalization
         DO i_sub=1, dim_subspace
            norm = sqrt(dot_product( state_new(i_sub,:) , state_new(i_sub,:) ))
            state_new(i_sub, :)=   state_new(i_sub, :) / norm
            !
         ENDDO

      ELSE
          CALL errore(subname,'problem in dim_subspace',ABS(dim_subspace))

      ENDIF



      !
      ! calculate B_matrix
      !
      DO i_sub=1, dim_subspace
         DO j_sub=1, dim_subspace
         !
         B_matrix(i_sub,j_sub) = (dot_product( state_new(i_sub,:) , state_new_aux(j_sub,:) ))
         !
         ENDDO
      ENDDO
      !

     ! numerical cut-off in order to avoid  diagonalization procedure's noise
      DO i_sub=1, dim_subspace
         DO i_recur=1, dim_recursion
            IF (ABS(state_new(i_sub, i_recur)) < EPS_m10 ) THEN
               state_new(i_sub, i_recur) = CZERO
            ENDIF
         ENDDO 
      ENDDO
      ! test normalisation :
   !
      test_norm = ZERO
      !
      DO i_sub=1, dim_subspace
         !
         !CALL mat_mul( test_norm_aux , state(i_sub,:), 'C', state(i_sub,:) , 'N', dim_vect, dim_vect, dim_recursion)
         test_norm_aux = dot_product(state_new(i_sub,:) , state_new(i_sub,:))
         !
         test_norm = test_norm + REAL(ABS(test_norm_aux)/dim_subspace, dbl)
      ENDDO
      !
      IF ( (test_norm >= ONE + EPS_m4).OR.( test_norm <= ONE - EPS_m4) ) &
            CALL errore(subname,'problem in normalization norm=',ABS(test_norm))




    END SUBROUTINE  calculate_recursion_state_diag_A





END MODULE recursion_function_module


