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
   USE T_control_module, ONLY :  in_datafile_C, nprint, dim_subspace,  &
                                 in_max_iter_C, print_gamma, imp_gamma
   USE T_egrid_module, ONLY : ne, emin, emax, delta, de
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
   WRITE(unit,"(  7x,'Hamiltonian data read from file :',5x,a)") TRIM(in_datafile_C)
   WRITE(unit,"(  7x,'iteration number for the hamiltonian part:',5x,i4)") in_max_iter_C
   WRITE( unit,"( 2x,'</INPUT FILES>',/)" )
   WRITE(unit,"(  7x,'Import Gamma  :',5x,a)") TRIM(log2char(imp_gamma))

   WRITE( unit,"( /,2x,'<ENERGY_GRID>')" )
   WRITE(unit,"(  7x,'Dimension   :',5x,i6)") ne
   WRITE(unit,"(  7x,'Min Energy  :',5x,f10.5)") emin
   WRITE(unit,"(  7x,'Max Energy  :',5x,f10.5)") emax
   WRITE(unit,"(  7x,'Energy Step :',5x,f10.5)") de
   WRITE(unit,"(  7x,'Delta       :',5x,f10.5)") delta
   WRITE( unit,"( 2x,'</ENERGY_GRID>',/)" )
   WRITE( unit,"( 2x,'</INPUT>',/)" )

END SUBROUTINE summary

