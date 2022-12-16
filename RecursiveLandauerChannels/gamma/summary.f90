!
! Copyright (C) 2005 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!**********************************************************
   SUBROUTINE summary(unit)
   !**********************************************************
   ! 
   ! Print out all the informnatins obtained from the 
   ! input and initialization routines.
   !
   USE kinds, ONLY : dbl
   USE constants,  ONLY : ZERO
   USE parser_module, ONLY: log2char
   USE T_control_module, ONLY :  in_datafile, &
                                 nprint, dim_subspace, max_iter_renorm, &
                                 max_iter_term,  in_max_iter, &
                                 method_sigma, print_gamma, print_gamma_tilde, a_analytique, b_analytique
   USE T_egrid_module, ONLY : ne, emin, emax,  de, delta_lead
   IMPLICIT NONE

   !
   ! input variables
   !
   INTEGER,   INTENT(in)         :: unit

   !
   ! local variables
   !
   INTEGER      :: i, ik

!--------------------------------------------------------

   !
   ! <INPUT> section
   !
   WRITE(unit,"()")      
   WRITE( unit,"( /,2x,'<INPUT>')" )
   WRITE(unit,"()")
   WRITE( unit,"(7x,'Print info each ', i3,' energy step' )" ) nprint
   WRITE(unit,"()")
   WRITE( unit,"( /,2x,'<INPUT FILES>')" )
   WRITE(unit,"(  7x,'Subspace dimensionnality:',5x,i4)") dim_subspace
   WRITE(unit,"(  7x,'Recursion data read from file  :',5x,a)") TRIM(in_datafile)
   WRITE(unit,"(  7x,'iteration number for the right part:',5x,i4)") in_max_iter
   WRITE( unit,"( 2x,'</INPUT FILES>',/)" )
   WRITE( unit,"( /,2x,'<Calculation>')" )
   WRITE(unit,"(  7x,'Max iteration number for the term  part   :',5x,i4)") max_iter_term
   WRITE(unit,"(  7x,'Max iteration number for the renorm  part :',5x,i4)") max_iter_renorm
   WRITE(unit,"(  7x,'Method for Sigma                          :',5x,a)") TRIM(method_sigma)
   WRITE(unit,"(  7x,'On site value :',5x,f10.5)") a_analytique
   WRITE(unit,"(  7x,'Hopping value :',5x,f10.5)") b_analytique
   WRITE( unit,"( 2x,'</Calculation>',/)" )
   WRITE(unit,"(  7x,'Print gamma        :',5x,a)") TRIM(log2char(print_gamma))
   WRITE(unit,"(  7x,'Print gamma_tilde  :',5x,a)") TRIM(log2char(print_gamma_tilde))
   WRITE( unit,"( /,2x,'</BULK PART>')" )
   WRITE( unit,"( 2x,'</INPUT>',/)" )
   WRITE( unit,"( /,2x,'<ENERGY_GRID>')" )
   WRITE(unit,"(  7x,'Dimension   :',5x,i6)") ne
   WRITE(unit,"(  7x,'Min Energy  :',5x,f10.5)") emin
   WRITE(unit,"(  7x,'Max Energy  :',5x,f10.5)") emax
   WRITE(unit,"(  7x,'Energy Step :',5x,f10.5)") de
   WRITE(unit,"(  7x,'Delta_lead  :',5x,f10.5)") delta_lead
   WRITE( unit,"( 2x,'</ENERGY_GRID>',/)" )

END SUBROUTINE summary

