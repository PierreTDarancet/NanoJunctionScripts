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
   USE T_control_module, ONLY :  in_datafile_L,  in_datafile_R, &
                                 nprint, dim_subspace, max_iter_R, &
                                 max_iter_L,  in_max_iter_R, in_max_iter_L, &
                                 conv_criterion, numerical_cut_off, &
                                 max_iter_CR, max_iter_CL, cut_chain, max_iter_CC
   USE T_egrid_module, ONLY : ne, emin, emax, delta, de, delta_lead
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
   WRITE(unit,"(  7x,'Right part data read from file  :',5x,a)") TRIM(in_datafile_R)
   WRITE(unit,"(  7x,'Left  part data read from file  :',5x,a)") TRIM(in_datafile_L)
   WRITE(unit,"(  7x,'iteration number for the right part:',5x,i4)") in_max_iter_R
   WRITE(unit,"(  7x,'iteration number for the left  part:',5x,i4)") in_max_iter_L
   WRITE( unit,"( 2x,'</INPUT FILES>',/)" )
   WRITE( unit,"( /,2x,'<RECURSION PART>')" )
   WRITE(unit,"(  7x,'Max iteration number for the left  part           :',5x,i4)") max_iter_L
   WRITE(unit,"(  7x,'Max iteration number for the conductor-left  part :',5x,i4)") max_iter_CL
   WRITE(unit,"(  7x,'Max iteration number for the conductor part       :',5x,i4)") max_iter_CC
   WRITE(unit,"(  7x,'Max iteration number for the conductor-right part :',5x,i4)") max_iter_CR
   WRITE(unit,"(  7x,'Max iteration number for the right part           :',5x,i4)") max_iter_R
   WRITE( unit,"( 2x,'</RECURSION PART>',/)" )
   WRITE(unit,"(  7x,'Convergency criterion  :',5x,f10.5)") conv_criterion
   WRITE(unit,"(  7x,'Numerical cut off  :',5x,f10.5)") numerical_cut_off
   WRITE(unit,"(  7x,'CUT CHAIN  :',5x,a)") TRIM(log2char(cut_chain))
   WRITE( unit,"( /,2x,'</BULK PART>')" )
   WRITE( unit,"( 2x,'</INPUT>',/)" )
   WRITE( unit,"( /,2x,'<ENERGY_GRID>')" )
   WRITE(unit,"(  7x,'Dimension   :',5x,i6)") ne
   WRITE(unit,"(  7x,'Min Energy  :',5x,f10.5)") emin
   WRITE(unit,"(  7x,'Max Energy  :',5x,f10.5)") emax
   WRITE(unit,"(  7x,'Energy Step :',5x,f10.5)") de
   WRITE(unit,"(  7x,'Delta       :',5x,f10.5)") delta
   WRITE(unit,"(  7x,'Delta_lead  :',5x,f10.5)") delta_lead
   WRITE( unit,"( 2x,'</ENERGY_GRID>',/)" )

END SUBROUTINE summary

