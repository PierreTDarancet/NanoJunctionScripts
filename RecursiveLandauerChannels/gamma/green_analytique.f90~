!
! Copyright (C) 2006 LEPES-CNRS
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!***********************************************************************************
   SUBROUTINE green_analytique( ene, g,  dim_mat,a_analytique, b_analytique )
   !***********************************************************************************
   !
   !  Construct green's functions
   !     
   !
   USE kinds
   USE constants, ONLY : CZERO, CONE, CI, ZERO
   USE timing_module, ONLY : timing
   IMPLICIT NONE

   !
   ! I/O variables
   !
   INTEGER,      INTENT(in)    :: dim_mat
   COMPLEX(dbl), INTENT(in)    :: ene
   REAL(dbl), INTENT(in)    :: a_analytique
   REAL(dbl), INTENT(in)    :: b_analytique

   COMPLEX(dbl), INTENT(out)   :: g(dim_mat,dim_mat)

   !
   ! local variables
   !
   INTEGER                   :: i_mat, j_mat, ierr, iter, shift_iter
   REAL(dbl)                 :: spectra_min, spectra_max
   COMPLEX(dbl)              :: work(dim_mat)
   COMPLEX(dbl)              :: delta
   CHARACTER(15)             :: subname='green_tout_neuf'

!
!----------------------------------------
! main Body
!----------------------------------------
!
   CALL timing('green_analytique',OPR='start')
!  local variables
   work(:) = CZERO
   g(:,:)=CZERO
!  output variable
    


    DO i_mat=1, dim_mat
      !
      spectra_min= a_analytique - 2 * ABS (b_analytique)
      spectra_max= a_analytique + 2 * ABS (b_analytique)
      !
      IF ( ( REAL(ene) <= spectra_max ) .AND. ( REAL (ene) >=  spectra_min )) THEN
            !
            ! b² - 4ac
            !
            delta = (a_analytique - ene)**2 - (4 * (b_analytique)**2)
            IF (REAL(delta) > ZERO )  CALL errore(subname, 'delta > 0', INT(ABS(delta+1)) )
            !
            work(i_mat)= (  -(a_analytique - ene) - CI* SQRT ( -delta ) ) / (2* (b_analytique**2))
            !
            work(i_mat)= CONE / (ene - a_analytique - ((b_analytique **2) * work(i_mat)))
            !
      ELSE 
            !
            work(i_mat) = CZERO
            !
      ENDIF
    ENDDO
    ! 
!     PRINT*, a_analytique
!     PRINT*, b_analytique
!     PRINT*, delta
    !
    g(:,:) = CZERO
    !
    DO i_mat=1, dim_mat
        !
        g(i_mat,i_mat) = work(i_mat)

        !
    ENDDO
    !

!


   CALL timing('green_analytique',OPR='stop')
    !

END SUBROUTINE green_analytique

