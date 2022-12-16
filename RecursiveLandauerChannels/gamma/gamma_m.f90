!
! Copyright (C) 2005 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*********************************************
   MODULE T_gamma_module
!*********************************************
   USE kinds,           ONLY : dbl
   USE parameters,      ONLY : nstrx
   USE constants,            ONLY : PI, ZERO, CZERO, CONE, CI, EPS_m5
   USE T_control_module,     ONLY : dim_subspace
   USE T_egrid_module,       ONLY : ne, emin, emax
   USE iotk_module

   IMPLICIT NONE
   PRIVATE 
   SAVE
!
! Contains transport hamiltonian data
! 
    ! 
    !
    COMPLEX(dbl), ALLOCATABLE :: sigma_s(:,:,:), sigma_tilde_s(:,:,:)
    !
    COMPLEX(dbl), ALLOCATABLE :: gamma_s(:,:,:), gamma_tilde_s(:,:,:)
    !

    !
    LOGICAL :: alloc = .FALSE.


!
! end delcarations
!

   PUBLIC :: sigma_s, sigma_tilde_s, gamma_tilde_s, gamma_s
   !
   PUBLIC :: alloc
   !
   PUBLIC :: gamma_allocate
   PUBLIC :: gamma_deallocate
   PUBLIC :: gamma_print
   PUBLIC :: gamma_tilde_print
   PUBLIC :: sigma_print
   PUBLIC :: sigma_tilde_print


CONTAINS

!
! subroutines
!

!**********************************************************
   SUBROUTINE gamma_allocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(14)      :: subname="gamma_allocate"
      INTEGER  :: ierr

      IF ( alloc )       CALL errore(subname,'already allocated', 1 )

      !
      !
      ALLOCATE ( sigma_s(ne,dim_subspace,dim_subspace), STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, ' allocating sigma_s', 1 )
      ALLOCATE ( sigma_tilde_s(ne,dim_subspace,dim_subspace), STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, ' allocating sigma_tilde_s', 1 )

      ALLOCATE ( gamma_tilde_s(ne,dim_subspace,dim_subspace), STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, ' allocating gamma_tilde_s', 1 )
      ALLOCATE ( gamma_s(ne,dim_subspace,dim_subspace), STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, ' allocating gamma_s', 1 )

      !
      !
      sigma_tilde_s(:,:,:) = CZERO
      sigma_s(:,:,:) = CZERO
      !
      gamma_tilde_s(:,:,:) = CZERO
      gamma_s(:,:,:) = CZERO
      !

      alloc = .TRUE.
      !
   END SUBROUTINE gamma_allocate


!**********************************************************
   SUBROUTINE gamma_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(22)      :: subname="hamiltonian_deallocate"
      INTEGER :: ierr

      IF ( .NOT. alloc ) RETURN

      DEALLOCATE ( sigma_s, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating sigma_s', 1 )
      DEALLOCATE ( sigma_tilde_s, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating sigma_tilde_s', 1 )
          !
      DEALLOCATE ( gamma_s, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating gamma_s', 1 )
      DEALLOCATE ( gamma_tilde_s, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating gamma_tilde_s', 1 )
          !

      alloc = .FALSE.   

   END SUBROUTINE gamma_deallocate


!**********************************************************
   SUBROUTINE gamma_print(unit, name)
   !**********************************************************
   IMPLICIT NONE
       INTEGER,         INTENT(in) :: unit
       CHARACTER(*),    INTENT(in) :: name
       CHARACTER(nstrx)   :: attr
       CHARACTER(11)      :: subname="gamma_print"
       INTEGER            :: ierr, ien
 
       CALL iotk_write_begin(unit,TRIM(name))

       CALL iotk_write_begin(unit,"DATA")
       CALL iotk_write_attr(attr,"nomega",ne) 
       CALL iotk_write_attr(attr,"emin",emin) 
       CALL iotk_write_attr(attr,"emax",emax) 
       CALL iotk_write_empty(unit,"DATA",ATTR=attr)

       CALL iotk_write_begin(unit,"gamma_term")

!        CALL iotk_write_begin(unit,"sigma_rp")
!           DO ien = 1, ne
!              CALL iotk_write_dat(unit,"sigma_rp"//TRIM(iotk_index(ien)), sigma_rp(ien,:,:))
!           ENDDO
!        CALL iotk_write_end(unit,"sigma_rp")
!        CALL iotk_write_begin(unit,"sigma_ip")
!            DO ien = 1, ne
!              CALL iotk_write_dat(unit,"sigma_ip"//TRIM(iotk_index(ien)), sigma_ip(ien,:,:))
!           ENDDO
!        CALL iotk_write_end(unit,"sigma_ip")

       DO ien = 1, ne
             !
             CALL iotk_write_dat(unit,"gamma"//TRIM(iotk_index(ien)), gamma_s(ien,:,:))
       ENDDO


       CALL iotk_write_end(unit,"gamma_term")
       CALL iotk_write_end(unit,TRIM(name))

   END SUBROUTINE gamma_print

!**********************************************************
   SUBROUTINE gamma_tilde_print(unit, name)
   !**********************************************************
   IMPLICIT NONE
       INTEGER,         INTENT(in) :: unit
       CHARACTER(*),    INTENT(in) :: name
       CHARACTER(nstrx)   :: attr
       CHARACTER(17)      :: subname="gamma_tilde_print"
       INTEGER            :: ierr, ien
 
       CALL iotk_write_begin(unit,TRIM(name))

       CALL iotk_write_begin(unit,"DATA")
       CALL iotk_write_attr(attr,"nomega",ne) 
       CALL iotk_write_attr(attr,"emin",emin) 
       CALL iotk_write_attr(attr,"emax",emax) 
       CALL iotk_write_empty(unit,"DATA",ATTR=attr)

       CALL iotk_write_begin(unit,"gamma_renorm")

       DO ien = 1, ne
             !
             CALL iotk_write_dat(unit,"gamma"//TRIM(iotk_index(ien)), gamma_tilde_s(ien,:,:))
       ENDDO

       CALL iotk_write_end(unit,"gamma_renorm")
       CALL iotk_write_end(unit,TRIM(name))

   END SUBROUTINE gamma_tilde_print


!**********************************************************
   SUBROUTINE sigma_print(unit, name)
   !**********************************************************
   IMPLICIT NONE
       INTEGER,         INTENT(in) :: unit
       CHARACTER(*),    INTENT(in) :: name
       CHARACTER(nstrx)   :: attr
       CHARACTER(13)      :: subname="sigma_print"
       INTEGER            :: ierr, ien
 
       CALL iotk_write_begin(unit,TRIM(name))

       CALL iotk_write_begin(unit,"DATA")
       CALL iotk_write_attr(attr,"nomega",ne) 
       CALL iotk_write_attr(attr,"emin",emin) 
       CALL iotk_write_attr(attr,"emax",emax) 
       CALL iotk_write_empty(unit,"DATA",ATTR=attr)

       CALL iotk_write_begin(unit,"sigma")

       DO ien = 1, ne
             !
             CALL iotk_write_dat(unit,"sigma"//TRIM(iotk_index(ien)), sigma_s(ien,:,:))
       ENDDO


       CALL iotk_write_end(unit,"sigma")
       CALL iotk_write_end(unit,TRIM(name))

   END SUBROUTINE sigma_print

!**********************************************************
   SUBROUTINE sigma_tilde_print(unit, name)
   !**********************************************************
   IMPLICIT NONE
       INTEGER,         INTENT(in) :: unit
       CHARACTER(*),    INTENT(in) :: name
       CHARACTER(nstrx)   :: attr
       CHARACTER(17)      :: subname="sigma_tilde_print"
       INTEGER            :: ierr, ien
 
       CALL iotk_write_begin(unit,TRIM(name))

       CALL iotk_write_begin(unit,"DATA")
       CALL iotk_write_attr(attr,"nomega",ne) 
       CALL iotk_write_attr(attr,"emin",emin) 
       CALL iotk_write_attr(attr,"emax",emax) 
       CALL iotk_write_empty(unit,"DATA",ATTR=attr)

       CALL iotk_write_begin(unit,"sigma")

       DO ien = 1, ne
             !
             CALL iotk_write_dat(unit,"sigma"//TRIM(iotk_index(ien)), sigma_tilde_s(ien,:,:))
       ENDDO

       CALL iotk_write_end(unit,"sigma_R")
       CALL iotk_write_end(unit,TRIM(name))

   END SUBROUTINE sigma_tilde_print

END MODULE T_gamma_module

