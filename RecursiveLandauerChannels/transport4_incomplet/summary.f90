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
   USE T_control_module, ONLY :  calculation_type, conduct_formula, &
                                 in_datafile_L,  in_datafile_R, &
                                 transport_dir,  &
                                 nprint, dim_subspace, max_iter_R, &
                                 max_iter_L,  in_max_iter_R, in_max_iter_L, &
                                 debug_mode, conv_criterion, numerical_cut_off, &
                                 max_iter_CR, max_iter_LC, cut_chain
   USE T_hamiltonian_module, ONLY :     dimC, &
                                        dimR, &
                                        dimL, &
                                        dim_CC
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
   WRITE(unit,"(  7x,'DEBUG MODE  :',5x,a)") TRIM(log2char(debug_mode))
   WRITE(unit,"(  7x,'Calculation Type    :',5x,a)") TRIM(calculation_type)
   WRITE(unit,"(  7x,'Conductance Formula :',5x,a)") TRIM(conduct_formula)
   WRITE(unit,"(  7x,'Transport Direction :',5x,i2)") transport_dir
   !WRITE(unit,"(  7x,'Max iteration number:',5x,i4)") niterx
   WRITE(unit,"()")
   WRITE( unit,"(7x,'Print info each ', i3,' energy step' )" ) nprint
   WRITE(unit,"()")
   WRITE( unit,"( /,2x,'<INPUT FILES>')" )
   WRITE(unit,"(  7x,'Right part data read from file  :',5x,a)") TRIM(in_datafile_R)
   WRITE(unit,"(  7x,'Left  part data read from file  :',5x,a)") TRIM(in_datafile_L)
   WRITE(unit,"(  7x,'Subspace dimensionnality:',5x,i4)") dim_subspace
   WRITE(unit,"(  7x,'iteration number for the right part:',5x,i4)") in_max_iter_R
   WRITE(unit,"(  7x,'iteration number for the left  part:',5x,i4)") in_max_iter_L
   WRITE( unit,"( 2x,'</INPUT FILES>',/)" )
   WRITE( unit,"( /,2x,'<RECURSION PART>')" )
   WRITE(unit,"(  7x,'Max iteration number for the right part:',5x,i4)") max_iter_R
   WRITE(unit,"(  7x,'Max iteration number for the left  part:',5x,i4)") max_iter_L
   WRITE(unit,"(  7x,'Convergency criterion  :',5x,f10.5)") conv_criterion
   WRITE(unit,"(  7x,'Numerical cut off  :',5x,f10.5)") numerical_cut_off
   WRITE(unit,"(  7x,'CUT CHAIN  :',5x,a)") TRIM(log2char(cut_chain))
   WRITE( unit,"( 2x,'</RECURSION PART>',/)" )
   IF (TRIM(conduct_formula) == 'mayou') THEN
      WRITE( unit,"( /,2x,'<MAYOU FORMULA>')" )
      WRITE(unit,"(  7x,'Max iteration number for the C right part:',5x,i4)") max_iter_CR
      WRITE(unit,"(  7x,'Max iteration number for the C left  part:',5x,i4)") max_iter_LC
      WRITE(unit,"(  7x,'Dimension of the conductor part:',5x,i6)") dim_CC
      WRITE( unit,"( 2x,'</MAYOU FORMULA>',/)" )
   ENDIF
   WRITE( unit,"( 2x,'</INPUT>',/)" )

   WRITE( unit,"( /,2x,'<ENERGY_GRID>')" )
   WRITE(unit,"(  7x,'Dimension   :',5x,i6)") ne
   WRITE(unit,"(  7x,'Min Energy  :',5x,f10.5)") emin
   WRITE(unit,"(  7x,'Max Energy  :',5x,f10.5)") emax
   WRITE(unit,"(  7x,'Energy Step :',5x,f10.5)") de
   WRITE(unit,"(  7x,'Delta       :',5x,f10.5)") delta
   WRITE(unit,"(  7x,'Delta_lead  :',5x,f10.5)") delta_lead
   WRITE( unit,"( 2x,'</ENERGY_GRID>',/)" )

   WRITE( unit,"( /,2x,'<OUT DIMENSION>')" )
   WRITE(unit,"(  7x,'Dimension of the conductor part:',5x,i6)") dimC
   WRITE(unit,"(  7x,'Dimension of the right     part:',5x,i6)") dimR
   WRITE(unit,"(  7x,'Dimension of the left      part:',5x,i6)") dimL
   WRITE( unit,"( 2x,'</OUT DIMENSION>',/)" )

END SUBROUTINE summary

