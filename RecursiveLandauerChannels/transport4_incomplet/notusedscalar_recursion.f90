!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*******************************************************************
   SUBROUTINE scalar_recursion(rec_chain, nb_first, id_first, ene_first, ene_onsite, dim_rec, max_iter, dim_subspace, var_sum)
   !*******************************************************************

   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, EPS_m4, CI, ONE
   USE timing_module



   IMPLICIT NONE
   !
   ! Input variables
   !

   INTEGER, INTENT(in)      :: max_iter
   INTEGER, INTENT(in)      :: dim_subspace
   INTEGER, INTENT(in)      :: dim_rec
   !COMPLEX(dbl), INTENT(in) :: rec_hamiltonian(dim_rec,dim_rec)
   COMPLEX(dbl), INTENT(in) :: nb_first(dim_rec)
   COMPLEX(dbl), INTENT(in) :: id_first(((dim_subspace*3)-1),dim_rec)
   COMPLEX(dbl), INTENT(in) :: ene_onsite(dim_rec)
   COMPLEX(dbl), INTENT(in) :: ene_first(((dim_subspace*3)-1),dim_rec)

   !
   ! Output variables
   REAL(dbl), INTENT(out)    :: var_sum(max_iter,dim_subspace,dim_subspace)
   COMPLEX(dbl), INTENT(out) :: rec_chain(2,max_iter,dim_subspace,dim_subspace)
   !

   !
   ! local variables
   !
   CHARACTER(16) :: subname="scalar_recursion"
   INTEGER       :: ierr, i_sub, j_sub, k_sub
   COMPLEX(dbl)  :: trial_state(dim_rec)
   INTEGER       :: nb_max_first 

   !
   ! end of declarations
   !

!
!----------------------------------------
! main Body
!----------------------------------------
!
   CALL timing( 'scalar_recursion', OPR='start' )

   nb_max_first = (dim_subspace*3)-1
   PRINT*, 'Entering scalar recursion part'
   !
   ! allocations
   !
!    ALLOCATE ( trial_state(dim_recursion), STAT=ierr )
!         IF( ierr /=0 ) CALL errore(subname, 'allocating recursion_state', ABS(ierr) )


    DO i_sub=1, dim_subspace
       !

       DO j_sub=1, dim_subspace
           PRINT*, 'Recursion for i_sub='
           PRINT*, i_sub
           PRINT*, 'Recursion for j_sub='
           PRINT*, j_sub
           !
           trial_state(:)=CZERO
           !

           IF  ( i_sub == j_sub ) THEN
               !
               trial_state(i_sub) = CONE
               !
           ELSE IF ( i_sub > j_sub ) THEN
               !
               trial_state(i_sub) = CONE
               trial_state(j_sub) = CONE
               trial_state(:)=trial_state(:) / SQRT(2*ONE)
               !
           ELSE IF ( i_sub < j_sub ) THEN
               !
               trial_state(i_sub) = CONE
               trial_state(j_sub) = CI
               trial_state(:)=trial_state(:) / SQRT(2*ONE)
               !
           ELSE
               !
               CALL errore(subname, 'problem with i_sub and j_sub', 1 )
               !
           ENDIF


           DO k_sub=1,dim_subspace
               PRINT*, 'Trial state for i_rec ='
               PRINT*,  k_sub, trial_state(k_sub)
           ENDDO
           !
           CALL recursion_loop(trial_state(:), nb_first(:), id_first(:,:), ene_first(:,:), ene_onsite(:), dim_rec, nb_max_first, max_iter, rec_chain(:,:,i_sub,j_sub), var_sum(:,i_sub,j_sub) )
           !

       ENDDO
       !
    ENDDO
    !



   !
   ! local cleaning
   !
!      IF ( ALLOCATED( trial_state  ) ) THEN
!            DEALLOCATE ( trial_state, STAT=ierr )
!            IF( ierr /=0 ) CALL errore(subname, 'deallocating trial_state', ABS(ierr) )
!       ENDIF

   CALL timing( 'scalar_recursion', OPR='stop' )


END SUBROUTINE scalar_recursion
