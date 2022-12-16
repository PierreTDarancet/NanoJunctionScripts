!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE identity_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE
   IMPLICIT NONE
   PRIVATE
   SAVE
   !
   TYPE orbitale
        CHARACTER(2)     :: nature      ! ("s","pz","px",py...)
        REAL(dbl)        :: coord(3)    ! coordinates
        LOGICAL          :: init
   END TYPE orbitale
   !
   TYPE  orbital_def
         CHARACTER(2)  ::   name
         COMPLEX(dbl) ::   onsite
         COMPLEX(dbl) ::      s(3)
         COMPLEX(dbl) ::      px(3)
         COMPLEX(dbl) ::      py(3)
         COMPLEX(dbl) ::      pz(3)
         COMPLEX(dbl) ::      d(3)
         COMPLEX(dbl) ::      f(3)
   END TYPE orbital_def
   !


 PUBLIC :: orbitale
 PUBLIC :: orbital_def


END MODULE identity_module