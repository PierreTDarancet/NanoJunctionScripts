!
! Copyright (C) 2004 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! <INFO>
!*********************************************
   MODULE summary_module
!*********************************************
    USE parameters, ONLY : nstrx
    USE kinds, ONLY : dbl
    USE constants,  ONLY : ZERO, PI, TPI, BOHR => bohr_radius_angs
    USE parser_module, ONLY: log2char
    USE io_module, ONLY : title, prefix, postfix, work_dir

    USE control_variable_module, ONLY : max_iter, min_iter, conv_criteria, &
                                        calculation_mode, variation_mode, &
                                        print_hamiltonian, print_center, print_state, &
                                        print_variation, print_A_sum, print_B_sum, print_matrix, &
                                        print_overlap, debug_mode, print_distance, print_orbital, &
                                        use_second, use_third
    USE dim_variable_module,  ONLY :  n_orb, dim_subspace, dim_recursion, limit_0, &
                                      limit_1, limit_2, limit_3, nb_max_first, nb_max_second, &
                                      nb_max_third
    USE subspace_variable_module, ONLY : x_limit, perform_recursion_x_greater, &
                                         perform_recursion_x_lesser
    USE orbital_module,          ONLY : orb, orb_recursion

   USE identity_module,      ONLY : orbitale


   IMPLICIT NONE
   PRIVATE

! 
! Print out all the informnatins obtained from the 
! input and initialization routines.
!
! output layout :
!
!
! contains:
! SUBROUTINE  summary(unit[,linput][,llattice][,latoms][,lpseudo][,lkpoints][,leig])
! </INFO>
!

    PUBLIC :: summary_input
    PUBLIC :: summary_recursion_input
    PUBLIC :: summary_recursion
!       PUBLIC :: summary_recursion
!       PUBLIC :: summary_output
             
!
! end delcarations
!

   CONTAINS

!
! subroutines
!
   !

!**********************************************************
   SUBROUTINE summary_input(unit)
   !**********************************************************
   ! 
   ! Print out all the informnatins obtained from the 
   ! input and initialization routines.
   !

   !
   ! input variables
   !
   INTEGER,   INTENT(in)         :: unit
   !
   ! local variables
   !

!--------------------------------------------------------

   !
   ! <INPUT> section
   !
   WRITE(unit,"()")      
   WRITE(unit,"()")      
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',32x,'Main',32x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )
   WRITE(unit,"(  7x,'     Calculation Title :',5x,a)") TRIM(title)
   WRITE(unit,"(  7x,'                Prefix :',5x,a)") TRIM(prefix)
   WRITE(unit,"(  7x,'               Postfix :',5x,a)") TRIM(postfix)
   IF ( LEN_TRIM(work_dir) <= 65 ) THEN
      WRITE(unit,"(  7x,'     Working directory :',5x,a)") TRIM(work_dir)
   ELSE
      WRITE(unit,"(  7x,'     Working directory :',5x,/,10x,a)") TRIM(work_dir)
   ENDIF
   WRITE(unit,"()")
   WRITE(unit,"( /,2x,'<INPUT>')" )
   WRITE(unit,"(  7x,'DEBUG MODE              :',5x,a)") log2char(debug_mode)

   WRITE(unit,"( /,2x,'<INPUT for hamiltonian part>')" )
   WRITE(unit,"(  7x,'Number of orbital :',5x,i4)") n_orb
   WRITE(unit,"(  7x,'Use second        :',5x,a)") log2char(use_second)
   WRITE(unit,"(  7x,'Use third         :',5x,a)") log2char(use_third)
   WRITE(unit,"(  7x,'Max number of first  neighbors :',5x,i4)") nb_max_first
   WRITE(unit,"(  7x,'Max number of second neighbors :',5x,i4)") nb_max_second
   WRITE(unit,"(  7x,'Max number of third  neighbors :',5x,i4)") nb_max_third

   WRITE(unit,"( /,2x,'<INPUT for recursion part>')" )
   WRITE(unit,"(  7x,'X limit                            :',5x,f10.5)")  x_limit
   WRITE(unit,"(  7x,'Dimension of initial subspace      :',5x,i4)") dim_subspace
   WRITE(unit,"(  7x,'Dimension of recursion Hamiltonian :',5x,i4)") dim_recursion
   WRITE(unit,"(  7x,'Recursion for x >                  :',5x,a)") log2char(perform_recursion_x_greater)
   WRITE(unit,"(  7x,'Recursion for x <                  :',5x,a)") log2char(perform_recursion_x_lesser)
   WRITE(unit,"(  7x,'Max iteration number:',5x,i4)") max_iter
   WRITE(unit,"(  7x,'Min iteration number:',5x,i4)") min_iter
   WRITE(unit,"(  7x,'Convergence criteria:',5x,f10.5)") conv_criteria
   WRITE(unit,"(  7x,'Method for variation  :',5x,a)") TRIM(variation_mode)
   WRITE(unit,"(  7x,'Method for recursion states construction  :',5x,a)") TRIM(calculation_mode)
!   WRITE(unit,"( /,2x,'<INPUT for central region >')" )
!   WRITE(unit,"(  7x,'Hamiltonian data in the central region read from file  :',5x,a)") TRIM(datafile_H)
!   WRITE(unit,"(  7x,'Size of the central part:',5x,f10.5)") cell_size_C

!   WRITE(unit,"(  7x,'Number of replica in the right part :',5x,i4)") nb_replica_R
!   WRITE(unit,"(  7x,'Dimension of input Hamiltonian in the right part :',5x,i4)") dim_R
!   WRITE(unit,"(  7x,'Size of the right part  :',5x,f10.5)") cell_size_R
!   WRITE(unit,"(  7x,'Number of replica in the left part  :',5x,i4)") nb_replica_L
!   WRITE(unit,"(  7x,'Dimension of input Hamiltonian in the left part  :',5x,i4)") dim_L
!   WRITE(unit,"(  7x,'Size of the left part   :',5x,f10.5)") cell_size_L
!   WRITE(unit,"(  7x,'Total Dimension of input Hamiltonian :',5x,i4)") dimwan_total
   WRITE(unit,"()")      
   WRITE(unit,"()")
   WRITE(unit,"( /,2x,'<OUTPUT FILES>')" )
   WRITE(unit,"(  7x,'Print output file for matrix           :',5x,a)") log2char(print_matrix)
   WRITE(unit,"(  7x,'Print output file for hamiltonian      :',5x,a)") log2char(print_hamiltonian)
   WRITE(unit,"(  7x,'Print output file for centers          :',5x,a)") log2char(print_center)
   WRITE(unit,"(  7x,'Print output file for states           :',5x,a)") log2char(print_state)
   WRITE(unit,"(  7x,'Print output file for variation        :',5x,a)") log2char(print_variation)
   WRITE(unit,"(  7x,'Print output file for A matrix summary :',5x,a)") log2char(print_A_sum)
   WRITE(unit,"(  7x,'Print output file for B matrix summary :',5x,a)") log2char(print_B_sum)
   WRITE(unit,"(  7x,'Print output file for overlap matrix   :',5x,a)") log2char(print_overlap)
   WRITE(unit,"(  7x,'Print output file for distance matrix  :',5x,a)") log2char(print_distance)
   WRITE(unit,"(  7x,'Print output file for orbital          :',5x,a)") log2char(print_orbital)
   WRITE(unit,"()")      
   WRITE(unit,"()")
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',24x,'Enter recursion part',24x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )

   END SUBROUTINE summary_input


!**********************************************************
   SUBROUTINE summary_recursion_input(unit,A_mat, B_mat, state, wan_num, subrecur2tot, norm, alloc)
   !**********************************************************
   !
   ! Print out all the informnatins obtained from the 
   ! input and initialization routines.
   !

   !
   ! input variables
   !
   INTEGER,   INTENT(in)         :: wan_num(dim_subspace)

   INTEGER,   INTENT(in)         :: unit
   LOGICAL, INTENT(in) :: alloc
   COMPLEX(dbl), INTENT(in)  :: A_mat(dim_subspace,dim_subspace), B_mat(dim_subspace,dim_subspace)
   COMPLEX(dbl), INTENT(in)  :: state(dim_subspace,dim_recursion)
   !
   INTEGER, INTENT(in) :: subrecur2tot(dim_recursion)
   REAL(dbl), INTENT(in) :: norm


   ! local variables
   !
   INTEGER   :: i_sub, j_sub, k_sub
   INTEGER   :: i_recur, j_recur
   INTEGER   :: i_wan, j_wan
   INTEGER   :: iter=0
!--------------------------------------------------------
                                 

   !
   !
   WRITE(unit,"()")      
   WRITE(unit,"()")      
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',23x,'Recursion Input summary',22x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )
   WRITE(unit,"()")
   WRITE(unit,"(  7x,'Recursion allocate   :',5x,a)") log2char(alloc)
   WRITE(unit,"()")
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',28x,'Coordinates',29x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )
   WRITE(unit,"( /,2x,'<Coordinates of WF in INPUT file>')" )
   WRITE(unit,"()")
   DO i_wan=1,n_orb
       WRITE( unit, "(7x, 'Coordinates of WF n° ',i4,' =  ( ',3f9.5,' ) ')") i_wan, orb(i_wan)%coord(:)
   ENDDO
   WRITE(unit,"()")
   WRITE(unit,"()")
   WRITE(unit, "(7x, 'WF used in Recursion :  ( ',300i4,' ) ')") (subrecur2tot(i_recur), i_recur=1,dim_recursion )
   WRITE(unit,"()")
!    WRITE( unit,"( /,2x,'<Coordinates of WF used in recursion>')" )
!    DO i_recur=1,dim_recursion
!        WRITE( unit, "(7x, 'Coordinates of WF n° ',i4,' =  ( ',3f9.5,' ) ')") i_recur, orb_recursion(i_recur)%coord(:)
!    ENDDO
!    WRITE(unit,"()")
   WRITE(unit,"()")
   WRITE( unit, "(7x, 'Initial Subspace WF :   ( ',300i4,' ) ')") ( wan_num(i_sub), i_sub=1,dim_subspace )
   WRITE(unit,"()")
   WRITE( unit,"( /,2x,'<Coordinates of WF in initial subspace>')" )
   DO i_sub=1,dim_subspace
       WRITE( unit, "(7x, 'Coordinates of subspace WF n° ',i4,' =  ( ',3f9.5,' ) ')") i_sub, orb(wan_num(i_sub))%coord(:)
   ENDDO
   WRITE(unit,"()")
   WRITE(unit,"()")
!    WRITE(unit,"(2x,70('='))" )
!    WRITE(unit,"(2x,'=',26x,'Hamiltonian DATA',26x,'=')" )
!    WRITE(unit,"(2x,70('='),/)" )
!    WRITE( unit,"( /,2x,'<Read Hamiltonian>')" )
!    DO i_wan=1,dimwan_total
!        WRITE( unit, " (7x, 'WF ',i4,': H_tot = ( ',300e9.5,' )   ')") i_wan, &
!                                                      ( totham(i_wan,j_wan), j_wan=1,dimwan_total )
!    ENDDO
!    WRITE( unit,"( /,2x,'</Read Hamiltonian>')" )
!    WRITE(unit,"()")      
   WRITE(unit,"()")      
!    WRITE( unit,"( /,2x,'<Recursion Hamiltonian>')" )
!    DO i_recur=1,dim_recursion
!        WRITE( unit, " (7x, 'WF ',i4,': H_rec = ( ',300e9.5,' )   ')") i_recur, &
!                                                      ( recham(i_recur,j_recur), j_recur=1,dim_recursion )
!    ENDDO
!    WRITE( unit,"( /,2x,'</Recursion Hamiltonian>')" )
! 
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',27x,'Initial states',27x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )
   DO i_recur=1,dim_recursion
        WRITE( unit, " (7x, 'on WF (',i4,') :  ( ',15f9.5,' ) ( ',15f9.5,')   ')") i_recur, &
                                 (REAL(state(i_sub, i_recur)), AIMAG(state(i_sub, i_recur)) , i_sub=1,dim_subspace)
   ENDDO
   WRITE(unit,"()")      
!   WRITE( unit, " (7x, 'Total Hamiltonian Hermitianity = ',f9.5,'  ')")  herm_tot
!   WRITE(unit,"()")      
!   WRITE(unit,"()")      
!   WRITE( unit, " (7x, 'Rec   Hamiltonian Hermitianity = ',f9.5,'  ')")  herm_rec
!   WRITE(unit,"()")      
 

  WRITE( unit, "(7x, 'Iteration = (',i4,')  ')") iter
   WRITE( unit, " (7x, 'Cols : A matrix  ',15i15,' ')") ( j_sub, j_sub=1,dim_subspace )
   DO i_sub=1,dim_subspace
       WRITE( unit, " (7x, 'Line ',i4,': A = ( ',15f9.5,' ) ( ',15f9.5,' )  ')") i_sub, &
                                                     ( REAL(A_mat(i_sub,j_sub)), AIMAG(A_mat(i_sub,j_sub)), j_sub=1,dim_subspace )
   ENDDO

   WRITE( unit, " (7x, 'Cols : B matrix ',15i15,'  ')")  ( k_sub, k_sub=1,dim_subspace )
   DO i_sub=1,dim_subspace
       WRITE( unit, " (7x, 'Line ',i4,': B = ( ',15f9.5,' ) ( ',15f9.5,' )  ')") i_sub, &
                                                     ( REAL(B_mat(i_sub,k_sub)), AIMAG(B_mat(i_sub,k_sub)), k_sub=1,dim_subspace )
   ENDDO
   WRITE(unit,"()")      
   WRITE( unit, " (7x, 'State Norm = ',f9.5,'  ')")  norm
   WRITE(unit,"()")      



   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',24x,'Enter recursion loop',24x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )

   END SUBROUTINE summary_recursion_input

   
   !**********************************************************
      SUBROUTINE summary_recursion(unit, iter, A_mat, B_mat, variation, var, norm)
      !**********************************************************

   !

   ! 
   ! Print out all the informnatins obtained from the 
   ! recursion routines.
   !

   !
   ! input variables
   !
   INTEGER,   INTENT(in)         :: unit, iter
   REAL(dbl), INTENT(in)         :: variation(dim_subspace), var, norm
   COMPLEX(dbl), INTENT(in)  :: A_mat(dim_subspace,dim_subspace), B_mat(dim_subspace,dim_subspace)

   !
   ! local variables
   !
   INTEGER :: i_sub, j_sub, k_sub
     !

   WRITE( unit, "(7x, 'Iteration = (',i4,')  ')") iter

   WRITE( unit, " (7x, 'Cols : A matrix ',15i15,' ')") ( j_sub, j_sub=1,dim_subspace )
   DO i_sub=1,dim_subspace
       WRITE( unit, " (7x, 'Line ',i4,': A = ( ',15f9.5,' ) ( ',15f9.5,' )  ')") i_sub, &
                                                     ( REAL(A_mat(i_sub,j_sub)), AIMAG(A_mat(i_sub,j_sub)), j_sub=1,dim_subspace )
   ENDDO

   WRITE( unit, " (7x, 'Cols : B matrix ',15i15,' ')")  ( k_sub, k_sub=1,dim_subspace )
   DO i_sub=1,dim_subspace
       WRITE( unit, " (7x, 'Line ',i4,': B = ( ',15f9.5,' ) ( ',15f9.5,' )  ')") i_sub, &
                                                     ( REAL(B_mat(i_sub,k_sub)), AIMAG(B_mat(i_sub,k_sub)), k_sub=1,dim_subspace )
   ENDDO
   WRITE(unit,"()")      
   WRITE( unit, " (7x, 'State Norm = ',f15.9,'  ')")  norm
   WRITE(unit,"()")      
          ! variation print
     WRITE( unit, "(  2x,'Variation decomposition : ')")
   DO i_sub=1,dim_subspace
     WRITE( unit, "(  4x,'Variation on ',i4,'    =   ', f13.6 ) " ) i_sub, variation(i_sub)
   ENDDO

     WRITE( unit, "(  4x,'Total Variation      =   ', f13.6 ) " ) var
     WRITE( unit, "()")
   IF (.NOT.debug_mode) THEN
         IF ((var <= conv_criteria) .AND. (iter > 1) .AND. (iter > min_iter) ) THEN
            WRITE(unit,"(2x,70('='))" )
            WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
            WRITE( unit, "( 2x,'Total Variation < Convergence criterion')" )
            WRITE( unit, "()")
            WRITE( unit, "( 2x,'Take decision to leave') " )
            WRITE( unit, "()")
            WRITE(unit,"(2x,70('='))" )
         ELSE IF ((var <= conv_criteria) .AND. ((iter <= min_iter) .OR. (iter==1))) THEN
            WRITE(unit,"(2x,70('='))" )
            WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
            WRITE( unit, "( 2x,'Min number of iterations not yet reached')" )
            WRITE( unit, "()")
            WRITE( unit, "( 2x,'Take decision to continue') " )
            WRITE( unit, "()")
            WRITE(unit,"(2x,70('='))" )
         ELSE IF (iter==max_iter) THEN
            WRITE(unit,"(2x,70('='))" )
            WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
            WRITE( unit, "( 2x,'Max number of iterations reached')" )
            WRITE( unit, "()")
            WRITE( unit, "( 2x,'Take decision to leave') " )
            WRITE( unit, "()")
            WRITE(unit,"(2x,70('='))" )
         ELSE 
            WRITE(unit,"(2x,70('='))" )
            WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
            WRITE( unit, "( 2x,'Total Variation > Convergence criterion')" )
            WRITE( unit, "()")
            WRITE( unit, "( 2x,'Take decision to continue') " )
            WRITE( unit, "()")
            WRITE(unit,"(2x,70('='))" )
         ENDIF
    ELSE 
         IF (iter==max_iter) THEN
            WRITE(unit,"(2x,70('='))" )
            WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
            WRITE( unit, "( 2x,'Max number of iterations reached')" )
            WRITE( unit, "()")
            WRITE( unit, "( 2x,'Take decision to leave') " )
            WRITE( unit, "()")
            WRITE(unit,"(2x,70('='))" )
         ELSE
            WRITE(unit,"(2x,70('='))" )
            WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
            WRITE( unit, "( 2x,'DEBUG MODE')" )
            WRITE( unit, "()")
            WRITE( unit, "( 2x,'Take decision to continue') " )
            WRITE( unit, "()")
            WRITE(unit,"(2x,70('='))" )
        ENDIF
    ENDIF
!          ELSE IF  THEN
!             WRITE(unit,"(2x,70('='))" )
!             WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
!             WRITE( unit, "( 2x,'End of the First iteration')" )
!             WRITE( unit, "()")
!             WRITE( unit, "( 2x,'Take decision to continue') " )
!             WRITE( unit, "()")
!             WRITE(unit,"(2x,70('='))" )

   END SUBROUTINE summary_recursion

END MODULE summary_module

