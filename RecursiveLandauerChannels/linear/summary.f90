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
    USE control_variable_module, ONLY : datafile_H, max_iter, conv_criteria
    USE dim_variable_module,  ONLY :  dimwan, dim_subspace, dim_recursion 
    USE subspace_variable_module, ONLY : x_limit, perform_recursion_x_greater, &
                                           perform_recursion_x_lesser



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
   WRITE(unit,"(  7x,'Hamiltonian data read from file  :',5x,a)") TRIM(datafile_H)
   WRITE(unit,"(  7x,'X limit :',5x,f10.5)")  x_limit
   WRITE(unit,"(  7x,'Dimension of initial subspace:',5x,i4)") dim_subspace
   WRITE(unit,"(  7x,'Recursion x >       :',5x,a)") log2char(perform_recursion_x_greater)
   WRITE(unit,"(  7x,'Recursion x <       :',5x,a)") log2char(perform_recursion_x_lesser)
   WRITE(unit,"(  7x,'Max iteration number:',5x,i4)") max_iter
   WRITE(unit,"(  7x,'Convergence criteria:',5x,f10.5)") conv_criteria
   WRITE(unit,"(  7x,'Dimension of input Hamiltonian:',5x,i4)") dimwan
   WRITE(unit,"(  7x,'Dimension of recursion Hamiltonian:',5x,i4)") dim_recursion

   WRITE(unit,"()")      
   WRITE(unit,"()")      
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',24x,'Enter recursion part',24x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )

   END SUBROUTINE summary_input


!**********************************************************
   SUBROUTINE summary_recursion_input(unit,A_mat, B_mat, state, recham, totham, wancoord, wan_num, wancoordrec, subrecur2tot, norm, alloc)
   !**********************************************************
   !
   ! Print out all the informnatins obtained from the 
   ! input and initialization routines.
   !

   !
   ! input variables
   !
   INTEGER,   INTENT(in)         :: unit
   LOGICAL, INTENT(in) :: alloc
   COMPLEX(dbl), INTENT(in)  :: A_mat(dim_subspace,dim_subspace), B_mat(dim_subspace,dim_subspace)
   COMPLEX(dbl), INTENT(in)  :: state(dim_subspace,dim_recursion)
   COMPLEX(dbl), INTENT(in)  :: recham(dim_recursion,dim_recursion), totham(dimwan,dimwan)
   !
   REAL(dbl), INTENT(in) :: wancoord(3,dimwan)
   INTEGER, INTENT(in) :: wan_num(dim_subspace)
   REAL(dbl), INTENT(in) :: wancoordrec(3,dim_recursion)
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
   DO i_wan=1,dimwan
       WRITE( unit, "(7x, 'Coordinates of WF n° ',i4,' =  ( ',3f9.5,' ) ')") i_wan, wancoord(:,i_wan)
   ENDDO
   WRITE(unit,"()")
   WRITE(unit,"()")
   WRITE(unit, "(7x, 'WF used in Recursion :  ( ',300i4,' ) ')") (subrecur2tot(i_recur), i_recur=1,dim_recursion )
   WRITE(unit,"()")
   WRITE( unit,"( /,2x,'<Coordinates of WF used in recursion>')" )
   DO i_recur=1,dim_recursion
       WRITE( unit, "(7x, 'Coordinates of WF n° ',i4,' =  ( ',3f9.5,' ) ')") i_recur, wancoordrec(:,i_recur)
   ENDDO
   WRITE(unit,"()")
   WRITE(unit,"()")
   WRITE( unit, "(7x, 'Initial Subspace WF :   ( ',300i4,' ) ')") ( wan_num(i_sub), i_sub=1,dim_subspace )
   WRITE(unit,"()")
   WRITE( unit,"( /,2x,'<Coordinates of WF in initial subspace>')" )
   DO i_sub=1,dim_subspace
       WRITE( unit, "(7x, 'Coordinates of subspace WF n° ',i4,' =  ( ',3f9.5,' ) ')") i_sub, wancoord(:,wan_num(i_sub))
   ENDDO
   WRITE(unit,"()")
   WRITE(unit,"()")
!    WRITE(unit,"(2x,70('='))" )
!    WRITE(unit,"(2x,'=',26x,'Hamiltonian DATA',26x,'=')" )
!    WRITE(unit,"(2x,70('='),/)" )
!    WRITE( unit,"( /,2x,'<Read Hamiltonian>')" )
!    DO i_wan=1,dimwan
!        WRITE( unit, " (7x, 'WF ',i4,': H_tot = ( ',300e9.5,' )   ')") i_wan, &
!                                                      ( totham(i_wan,j_wan), j_wan=1,dimwan )
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
   WRITE(unit,"(2x,'=',28x,'Initial state',27x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )
   DO i_recur=1,dim_recursion
        WRITE( unit, " (7x, 'on WF (',i4,') :  ( ',15f9.5,' ) ( ',15f9.5,')   ')") i_recur, &
                                 (REAL(state(i_sub, i_recur)), AIMAG(state(i_sub, i_recur)) , i_sub=1,dim_subspace)
   ENDDO
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

   ! !**********************************************************
   !    SUBROUTINE summary_output(unit)
   !    !**********************************************************
   !    END SUBROUTINE summary_output

!   WRITE(unit,"()")
!   WRITE( unit,"(7x,'Print info each ', i3,' energy step' )" ) nprint
!   WRITE(unit,"()")

      !    WRITE( unit,"( 2x,'</INPUT>',/)" )
      ! 
      !    WRITE( unit,"( /,2x,'<ENERGY_GRID>')" )
      !    WRITE(unit,"(  7x,'Dimension   :',5x,i6)") ne
      !    WRITE(unit,"(  7x,'Min Energy  :',5x,f10.5)") emin
      !    WRITE(unit,"(  7x,'Max Energy  :',5x,f10.5)") emax
      !    WRITE(unit,"(  7x,'Energy Step :',5x,f10.5)") de
      !    WRITE(unit,"(  7x,'Delta       :',5x,f10.5)") delta
      !    WRITE( unit,"( 2x,'</ENERGY_GRID>',/)" )
      ! 
      !    IF ( kpoints_alloc ) THEN
      !        !
      !        WRITE( unit, " ( /,2x,'<K-POINTS>')" )
      !        WRITE( unit, "(7x, 'nkpts_par = ',i4 ) " ) nkpts_par
      !        WRITE( unit, "(7x, 'nrtot_par = ',i4 ) " ) nrtot_par
      !        !
      !        ! 3D vectors initialization
      !        !
      !        nk_par3D(:) = 1
      !        nr_par3D(:) = 1
      !        vkpt_par3D(:) = ZERO
      !        vr_par3D(:) = ZERO
      !        !
      !        SELECT CASE (transport_dir)
      ! 
      !        CASE (1)
      !            !
      !            ! transport direction along the x axis
      !            !
      !            nk_par3D(2) = nk_par(1) 
      !            nk_par3D(3) = nk_par(2) 
      !            !
      !            WRITE( unit, "(7x, 'Parallel kpoints grid:      nk = (',3i3,' )') " ) nk_par3D(:) 
      !            !
      !            DO ik=1,nkpts_par
      !               vkpt_par3D(2) = vkpt_par(1,ik)
      !               vkpt_par3D(3) = vkpt_par(2,ik)
      !               WRITE( unit, " (7x, 'k point ', i4, ':   ( ',3f9.5,' ),   weight = ', f8.4 )") &
      !               ik, ( vkpt_par3D(i), i=1,3 ), wk_par(ik)
      !            ENDDO
      !            !
      !            nr_par3D(2) = nr_par(1) 
      !            nr_par3D(3) = nr_par(2) 
      !            !
      !            WRITE( unit, "(/,7x, 'Parallel R vector grid:      nr = (',3i3,' )') " ) nr_par3D(:) 
      !            !
      !            DO ik=1,nrtot_par
      !               vr_par3D(2) = vr_par(1,ik)
      !               vr_par3D(3) = vr_par(2,ik)
      !               WRITE( unit, " (7x, 'R vector', i4, ':   ( ',3f9.5,' ),   weight = ', f8.4 )") &
      !               ik, ( vr_par3D(i), i=1,3 ), wr_par(ik)
      !            ENDDO
      !     
      !        CASE (2)
      !            !
      !            ! transport direction along the y axis
      !            !
      !            nk_par3D(1) = nk_par(1) 
      !            nk_par3D(3) = nk_par(2) 
      !            WRITE( unit, "(7x, 'Parallel kpoints grid:      nk = (',3i3,' )') " ) nk_par3D(:) 
      !            DO ik=1,nkpts_par
      !               vkpt_par3D(1) = vkpt_par(1,ik)
      !               vkpt_par3D(3) = vkpt_par(2,ik)
      !               WRITE( unit, " (7x, 'k point ', i4, ':   ( ',3f9.5,' ),   weight = ', f8.4 )") &
      !               ik, ( vkpt_par3D(i), i=1,3 ), wk_par(ik)
      !            ENDDO
      !            !
      !            nr_par3D(1) = nr_par(1) 
      !            nr_par3D(3) = nr_par(2) 
      !            WRITE( unit, "(/, 7x, 'Parallel R vector grid:      nr = (',3i3,' )') " ) nr_par3D(:) 
      !            DO ik=1,nrtot_par
      !               vr_par3D(1) = vr_par(1,ik)
      !               vr_par3D(3) = vr_par(2,ik)
      !               WRITE( unit, " (7x, 'R vector', i4, ':   ( ',3f9.5,' ),   weight = ', f8.4 )") &
      !               ik, ( vr_par3D(i), i=1,3 ), wr_par(ik)
      !            ENDDO
      !     
      !        CASE (3)
      !            !
      !            ! transport direction along the z axis
      !            !
      !            nk_par3D(1) = nk_par(1) 
      !            nk_par3D(2) = nk_par(2) 
      !            WRITE( unit, "(7x, 'Parallel kpoints grid:      nk = (',3i3,' )') " ) nk_par3D(:) 
      !            DO ik=1,nkpts_par
      !               vkpt_par3D(1) = vkpt_par(1,ik)
      !               vkpt_par3D(2) = vkpt_par(2,ik)
      !               WRITE( unit, " (7x, 'k point ', i4, ':   ( ',3f9.5,' ),   weight = ', f8.4 )") &
      !               ik, ( vkpt_par3D(i), i=1,3 ), wk_par(ik)
      !            ENDDO
      !            WRITE(unit,"()")
      !            nr_par3D(1) = nr_par(1) 
      !            nr_par3D(2) = nr_par(2) 
      !            WRITE( unit, "(7x, 'Parallel R vector grid:      nr = (',3i3,' )') " ) nr_par3D(:) 
      !            DO ik=1,nrtot_par
      !               vr_par3D(1) = vr_par(1,ik)
      !               vr_par3D(2) = vr_par(2,ik)
      !               WRITE( unit, " (7x, 'R vector', i4, ':   ( ',3f9.5,' ),   weight = ', f8.4 )") &
      !               ik, ( vr_par3D(i), i=1,3 ), wr_par(ik)
      !            ENDDO
      !        END SELECT
      !    ENDIF
      !    WRITE( unit, " ( 2x,'</K-POINTS>',/)" )

   
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

   IF ((var <= conv_criteria) .AND. (iter>1) ) THEN
   WRITE(unit,"(2x,70('='))" )
   WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
     WRITE( unit, "( 2x,'Total Variation < Convergence criterion')" )
     WRITE( unit, "()")
     WRITE( unit, "( 2x,'Take decision to leave') " )
     WRITE( unit, "()")
   WRITE(unit,"(2x,70('='))" )
   ELSE IF (iter==1) THEN
   WRITE(unit,"(2x,70('='))" )
   WRITE( unit, "(7x, 'At Iteration  (',i4,')  ')") iter
     WRITE( unit, "( 2x,'End of the First iteration')" )
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

   END SUBROUTINE summary_recursion

END MODULE summary_module

