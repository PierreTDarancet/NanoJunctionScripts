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
    USE control_variable_module, ONLY : datafile_H, use_second, &
                                        use_third, print_hamiltonian, &
                                        print_distance, print_orbital
    USE dim_variable_module,  ONLY : cell_size,  n_orb, &
                                     limit_0, limit_1, &
                                     limit_2, limit_3
    USE orbital_module,  ONLY : orb, orbital_alloc


!

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
   INTEGER :: i_orb
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
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"( /,2x,'<INPUT>')" )
   WRITE(unit,"(  7x,'Dimension  of  Hamiltonian     :',5x,i4)") n_orb
   WRITE(unit,"(  7x,'Use second neightbour coupling :',5x,a)") log2char(use_second)
   WRITE(unit,"(  7x,'Use third  neightbour coupling :',5x,a)") log2char(use_third)
   WRITE(unit,"(  7x,'Size of the elementary cell    :',5x,f10.5)") cell_size
   WRITE(unit,"(  7x,'Limit for coupling                   :',5x,f10.5)") limit_0
   WRITE(unit,"(  7x,'Limit for first  neightbour coupling :',5x,f10.5)") limit_1
   WRITE(unit,"(  7x,'Limit for second neightbour coupling :',5x,f10.5)") limit_2
   WRITE(unit,"(  7x,'Limit for third  neightbour coupling :',5x,f10.5)") limit_3
   WRITE(unit,"( /,2x,'<OUTPUT FILES>')" )
   WRITE(unit,"(  7x,'Hamiltonian data of the central part written in file :',5x,a)") TRIM(datafile_H)
   WRITE(unit,"(  7x,'Print hamiltonian :',5x,a)") log2char(print_hamiltonian)
   WRITE(unit,"(  7x,'Print orbital     :',5x,a)") log2char(print_orbital)
   WRITE(unit,"(  7x,'Print distance    :',5x,a)") log2char(print_distance)

   WRITE(unit,"( /,2x,'</INPUT>')" )
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"()")      
   WRITE(unit,"()")
   WRITE(unit,"()")
   WRITE(unit,"(  7x,'Orbital allocate   :',5x,a)") log2char(orbital_alloc)
   WRITE(unit,"(2x,70('='))" )
   WRITE(unit,"(2x,'=',28x,'Coordinates',29x,'=')" )
   WRITE(unit,"(2x,70('='),/)" )
   WRITE(unit,"( /,2x,'<Coordinates of WF in INPUT file>')" )
   WRITE(unit,"()")
   DO i_orb=1, n_orb
       WRITE( unit, "(7x, 'WF n° ',i4,' Type ',a,' Coordinates =  ( ',3f9.5,' ) ')") i_orb, orb(i_orb)%nature, orb(i_orb)%coord(:)
   ENDDO
   WRITE(unit,"()")
   WRITE(unit,"()")

   END SUBROUTINE summary_input



END MODULE summary_module


