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
   MODULE orbital_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : CZERO, CONE, ZERO
   USE identity_module,      ONLY : orbitale
   USE dim_variable_module,  ONLY : n_orb, dim_recursion
   USE iotk_module
   USE files_module, ONLY : file_open, file_close
   USE io_module,            ONLY : orb_unit => aux2_unit, &
                                    ioname
   USE io_global_module,     ONLY : stdin, stdout


   IMPLICIT NONE
   PRIVATE 
   SAVE

! Variables
  
    TYPE(orbitale), ALLOCATABLE :: orb(:)
    TYPE(orbitale), ALLOCATABLE :: orb_recursion(:)
    LOGICAL :: orbital_alloc

! Status

   PUBLIC :: orb
!
   PUBLIC :: orb_recursion
!
   PUBLIC :: orbital_alloc
!
! functions

  PUBLIC :: orbital_allocate
  PUBLIC :: orbital_deallocate
  PUBLIC :: orbital_init
  PUBLIC :: print_orbital


CONTAINS


!*******************************************************************
   SUBROUTINE orbital_allocate
   !*******************************************************************
      CHARACTER(16)      :: subname="orbital_allocate"
      INTEGER :: ierr


   IF (orbital_alloc)   CALL errore(subname, 'orb  already allocated', 1 )
   !
   ! allocate basic quantities
   !

   ALLOCATE( orb(n_orb),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating orb',ABS(ierr))
   !
   ALLOCATE( orb_recursion(dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating orb_recursion',ABS(ierr))


   orbital_alloc=.TRUE.


END SUBROUTINE orbital_allocate


!**********************************************************
   SUBROUTINE orbital_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(18)      :: subname="orbital_deallocate"
      INTEGER :: ierr

     IF (.NOT. orbital_alloc) CALL errore(subname, 'orb  already deallocated', 1 )


     IF ( ALLOCATED( orb  ) ) THEN
           DEALLOCATE ( orb, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating orb', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( orb_recursion  ) ) THEN
           DEALLOCATE ( orb_recursion, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating orb_recursion', ABS(ierr) )
      ENDIF

     orbital_alloc =.FALSE.

   END SUBROUTINE orbital_deallocate

!*******************************************************************
   SUBROUTINE orbital_init
   !*******************************************************************
      CHARACTER(12)      :: subname="orbital_init"
      INTEGER :: ierr, i_orb, i_rec


   IF (.NOT. orbital_alloc)   CALL errore(subname, 'orb  not allocated', 1 )
   !
   ! allocate basic quantities
   !


   DO i_orb=1, n_orb
       orb(i_orb)%nature=' '
       orb(i_orb)%coord(:)=ZERO
       orb(i_orb)%init = .FALSE.

   ENDDO

   DO i_rec=1, dim_recursion
       orb_recursion(i_rec)%nature=' '
       orb_recursion(i_rec)%coord(:)=ZERO
       orb_recursion(i_rec)%init = .FALSE.

   ENDDO


END SUBROUTINE orbital_init

!*******************************************************************
   SUBROUTINE print_orbital
   !*******************************************************************
       CHARACTER(13)      :: subname="print_orbital"
       CHARACTER(7) :: name="ORBITAL"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: iter
       INTEGER            :: ierr, i_orb, i_rec
      CHARACTER( LEN=nstrx )  :: filename
      INTEGER            :: nrtot_C
      INTEGER            :: nr_C(3)
      INTEGER            :: ivr_C(3,3)



         nrtot_C = 3
         nr_C(1) = 3
         nr_C(2) = 1
         nr_C(3) = 1 
   !
         !ivr_C(3,nrtot)
         ivr_C(:,:)=0
         ivr_C(1,1)=-1
         ivr_C(1,2)=0
         ivr_C(1,3)=1
 
       CALL ioname('orbital',filename)

       CALL file_open(orb_unit,TRIM(filename),PATH="/",ACTION="write", &
                              FORM='formatted')

         !
       CALL iotk_write_begin(orb_unit,TRIM(name))
       CALL iotk_write_attr(attr,"dimwann",n_orb,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"nkpts",nrtot_C) 
       CALL iotk_write_attr(attr,"nk",nr_C) 
       CALL iotk_write_attr(attr,"nrtot",nrtot_C) 
       CALL iotk_write_attr(attr,"nr",nr_C) 
       CALL iotk_write_empty(orb_unit,"DATA",ATTR=attr)
        !
       CALL iotk_write_dat(orb_unit,"IVR", ivr_C, ATTR=attr, COLUMNS=3, IERR=ierr) 
               IF (ierr/=0) CALL errore(subname,'writing ivr',ABS(ierr))
       !
       CALL iotk_write_begin(orb_unit,"ORB")

       DO i_orb=1,n_orb
            CALL iotk_write_dat(orb_unit,"ORB"//TRIM(iotk_index(i_orb)),orb(i_orb)%nature)
            CALL iotk_write_dat(orb_unit,"ORB"//TRIM(iotk_index(i_orb)),orb(i_orb)%coord(:))
       ENDDO


       CALL iotk_write_end(orb_unit,"ORB")
         !

       CALL iotk_write_begin(orb_unit,"ORB_REC")

       DO i_rec=1,dim_recursion
            CALL iotk_write_dat(orb_unit,"ORB_REC"//TRIM(iotk_index(i_rec)),orb_recursion(i_rec)%nature)
            CALL iotk_write_dat(orb_unit,"ORB_REC"//TRIM(iotk_index(i_rec)),orb_recursion(i_rec)%coord(:))
       ENDDO


       CALL iotk_write_end(orb_unit,"ORB_REC")
       !

       CALL iotk_write_end(orb_unit,TRIM(name))

      CALL file_close(orb_unit,PATH="/",ACTION="write")

      CALL ioname('orbital',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Orbital written on file : ',5x,a)") TRIM(filename)


END SUBROUTINE print_orbital



END MODULE orbital_module


