!
! Copyright (C) 2009 Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!***********************************************
   MODULE input_hamiltonian
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI

   IMPLICIT NONE

   COMPLEX(dbl), PARAMETER      ::   HOMO_left    = CMPLX( -1.385, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   LUMO_left    = CMPLX(  5.03 , ZERO, dbl)
   COMPLEX(dbl), PARAMETER      ::   HOMO_middle  = CMPLX( -2.572, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   LUMO_middle  = CMPLX(  4.186, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   HOMO_right   = CMPLX(  1.77 , ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   LUMO_right   = CMPLX( -4.522, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   H_lL         = CMPLX(  0.5, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   H_lR         = CMPLX(  0.1, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   H_mL         = CMPLX(  0.3, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   H_mR         = CMPLX(  0.3, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   H_rL         = CMPLX(  0.1, ZERO, dbl) 
   COMPLEX(dbl), PARAMETER      ::   H_rR         = CMPLX(  0.5, ZERO, dbl) 
   !
   !
   PUBLIC      ::   HOMO_left   
   PUBLIC      ::   LUMO_left   
   PUBLIC      ::   HOMO_middle 
   PUBLIC      ::   LUMO_middle 
   PUBLIC      ::   HOMO_right  
   PUBLIC      ::   LUMO_right  
   PUBLIC      ::   H_lL        
   PUBLIC      ::   H_lR        
   PUBLIC      ::   H_rL        
   PUBLIC      ::   H_rR        


  END  MODULE input_hamiltonian
