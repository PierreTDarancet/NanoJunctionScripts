!
! Copyright (C) 2005 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
 

!**********************************************************
   PROGRAM util_input
   !**********************************************************
   USE constants, ONLY: ZERO, CZERO
   USE kinds,                ONLY : dbl
   USE parameters,           ONLY : nstrx
   USE iotk_module
   USE files_module,	    ONLY : file_open, file_close, file_delete
   USE kinds, ONLY : dbl
   USE io_module
   USE parser_module
   IMPLICIT NONE
  ! ... 

!
! Contains correlation self-energy data
! 
    !
    !
    COMPLEX(dbl), ALLOCATABLE :: Hamiltonian(:,:,:)
    REAL(dbl), ALLOCATABLE :: coord(:,:)
    COMPLEX(dbl), ALLOCATABLE :: onsite(:,:), hopping(:,:)
    INTEGER, ALLOCATABLE :: ivr(:,:)
    INTEGER      ::  nrtot, dimwann, dim_subspace
      CHARACTER( LEN=nstrx )  :: filename
    CHARACTER(nstrx)   :: attr, subname
    INTEGER            :: ir, ierr, unit, nr_aux(3), maxiwan, iwan, shift, shift1, shift2

    filename="1D"
    subname="util"
        nrtot=3
        nr_aux(1)=3
        nr_aux(2)=1
        nr_aux(3)=1
        unit=19
        dim_subspace=2
        dimwann=20
 
   ALLOCATE( ivr(3,nrtot), STAT=ierr )
      IF (ierr/=0) CALL errore(subname, 'allocating ivr', ABS(ierr) )
   ALLOCATE( coord(3,dimwann), STAT=ierr )
      IF (ierr/=0) CALL errore(subname, 'allocating coord', ABS(ierr) )
   ALLOCATE( onsite(dim_subspace,dim_subspace), STAT=ierr )
      IF (ierr/=0) CALL errore(subname, 'allocating onsite', ABS(ierr) )
   ALLOCATE( hopping(dim_subspace,dim_subspace), STAT=ierr )
      IF (ierr/=0) CALL errore(subname, 'allocating hopping', ABS(ierr) )

       ivr(:,:)=0
       ivr(1,1)=-1
       ivr(1,2)=0
       ivr(1,3)=1

       !ivr(1,2)=0

       onsite(1,1)= 1.0
        onsite(1,2)= 2.0
!        onsite(1,3)= 3.0
        onsite(2,2)= 10.0
        onsite(2,1)= 2.0
!        onsite(2,3)= 3.0
!        onsite(3,1)= 3.0
!        onsite(3,2)= 3.0
!        onsite(3,3)= 5.0

       hopping(1,1)= CMPLX(3.0, 1.0)
        hopping(1,2)= 1.0
!        hopping(1,3)= 2.0
! 
        hopping(2,1)= 1.0
        hopping(2,2)=CMPLX(4.0, 1.0)
!        hopping(2,3)= 2.0
! 
!        hopping(3,1)= 2.0
!        hopping(3,2)= 2.0
!        hopping(3,3)= 3.0


   ALLOCATE( Hamiltonian(dimwann, dimwann, nrtot), STAT=ierr )
      IF (ierr/=0) CALL errore(subname, 'allocating Hamiltonian', ABS(ierr) )
PRINT*, 'AVANT INIT0'


  Hamiltonian(:,:,:) = CZERO
         !!!!!!!!!!!!!!!!! ivr = 0 0 0 (ir=2)
  ir=2
  Hamiltonian(1,1,ir)= onsite(1,1)
  Hamiltonian(1,2,ir)= onsite(1,2)
!   Hamiltonian(1,3,ir)= onsite(1,3)
   Hamiltonian(2,2,ir)= onsite(2,2)
   Hamiltonian(2,1,ir)= onsite(2,1)
!   Hamiltonian(2,3,ir)= onsite(2,3)
!   Hamiltonian(3,1,ir)= onsite(3,1)
!   Hamiltonian(3,2,ir)= onsite(3,2)
!   Hamiltonian(3,3,ir)= onsite(3,3)
  Hamiltonian(1,1+dim_subspace,ir)= hopping(1,1)
   Hamiltonian(1,2+dim_subspace,ir)= hopping(1,2)
!   Hamiltonian(1,3+dim_subspace,ir)= hopping(1,3)
   Hamiltonian(2,2+dim_subspace,ir)= hopping(2,2)
   Hamiltonian(2,1+dim_subspace,ir)= hopping(2,1)
!   Hamiltonian(2,3+dim_subspace,ir)= hopping(2,3)
! 
! 
!   Hamiltonian(3,1+dim_subspace,ir)= hopping(3,1)
!   Hamiltonian(3,2+dim_subspace,ir)= hopping(3,2)
!   Hamiltonian(3,3+dim_subspace,ir)= hopping(3,3)
!PRINT*, 'AVANT INIT'
  maxiwan=INT(dimwann/dim_subspace) - 1
  DO iwan=1, maxiwan
     shift=iwan*dim_subspace
     shift1=(iwan+1)*dim_subspace
     shift2=(iwan-1)*dim_subspace

             Hamiltonian(1+shift,1+shift,ir) = onsite(1,1)+iwan*0.000001
             Hamiltonian(2+shift,2+shift,ir) = onsite(2,2)+iwan*0.000002

     !Hamiltonian(1+shift,1+shift,ir) = onsite(1,1)
     !Hamiltonian(2+shift,2+shift,ir) = onsite(2,2)
     Hamiltonian(1+shift,2+shift,ir) = onsite(1,2)
     Hamiltonian(2+shift,1+shift,ir) = onsite(2,1)
      !      Hamiltonian(1+shift,3+shift,ir) = onsite(1,3)
      !      Hamiltonian(3+shift,1+shift,ir) = onsite(3,1)
      !      Hamiltonian(3+shift,3+shift,ir) = onsite(3,3)
      !      Hamiltonian(2+shift,3+shift,ir) = onsite(2,3)
      !      Hamiltonian(3+shift,2+shift,ir) = onsite(3,2)
     IF ( iwan <= (maxiwan -1) ) THEN

            Hamiltonian(1+shift,1+shift1,ir)= hopping(1,1)
            Hamiltonian(1+shift,2+shift1,ir)= hopping(1,2)
            Hamiltonian(2+shift,2+shift1,ir)= hopping(2,2)
            Hamiltonian(2+shift,1+shift1,ir)= hopping(2,1)
            !      Hamiltonian(1+shift,3+shift1,ir)= hopping(1,3)
            !      Hamiltonian(3+shift,1+shift1,ir)= hopping(3,1)
            !      Hamiltonian(2+shift,3+shift1,ir)= hopping(2,3)
            !      Hamiltonian(3+shift,2+shift1,ir)= hopping(3,2)
            !      Hamiltonian(3+shift,3+shift1,ir)= hopping(3,3)
    ENDIF

         Hamiltonian(1+shift,1+shift2,ir)= CONJG (hopping(1,1))
         Hamiltonian(1+shift,2+shift2,ir)= CONJG (hopping(2,1))
         !      Hamiltonian(1+shift,3+shift2,ir)= CONJG (hopping(3,1))
         Hamiltonian(2+shift,2+shift2,ir)= CONJG (hopping(2,2))
         Hamiltonian(2+shift,1+shift2,ir)= CONJG (hopping(1,2))
         !      Hamiltonian(2+shift,3+shift2,ir)= CONJG (hopping(3,2))
         !      Hamiltonian(3+shift,1+shift2,ir)= CONJG (hopping(1,3))
         !      Hamiltonian(3+shift,2+shift2,ir)= CONJG (hopping(2,3))
         !      Hamiltonian(3+shift,3+shift2,ir)= CONJG (hopping(3,3))


  ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ir=1
  iwan=maxiwan
     shift=iwan*dim_subspace

     Hamiltonian(1,1+shift,ir) = CONJG (hopping(1,1))
     Hamiltonian(2,2+shift,ir) = CONJG (hopping(2,2))
     Hamiltonian(1,2+shift,ir) = CONJG (hopping(2,1))
     Hamiltonian(2,1+shift,ir) = CONJG (hopping(1,2))
  !shiter colonne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ir=3
!   Hamiltonian(1+shift,1,ir)= onsite(1,1)
!   Hamiltonian(1+shift,2,ir)= onsite(1,2)
!   Hamiltonian(1,3,ir)= onsite(1,3)
!    Hamiltonian(2+shift,2,ir)= onsite(2,2)
!    Hamiltonian(2+shift,1,ir)= onsite(2,1)
!   Hamiltonian(2,3,ir)= onsite(2,3)
!   Hamiltonian(3,1,ir)= onsite(3,1)
!   Hamiltonian(3,2,ir)= onsite(3,2)
!   Hamiltonian(3,3,ir)= onsite(3,3)
  Hamiltonian(1+shift,1,ir)= hopping(1,1)
   Hamiltonian(1+shift,2,ir)= hopping(1,2)
!   Hamiltonian(1,3+dim_subspace,ir)= hopping(1,3)
   Hamiltonian(2+shift,2,ir)= hopping(2,2)
   Hamiltonian(2+shift,1,ir)= hopping(2,1)
!   Hamiltonian(2,3+dim_subspace,ir)= hopping(2,3)
!   !shiter ligne
! 
!   Hamiltonian(3,1+dim_subspace,ir)= hopping(3,1)
!   Hamiltonian(3,2+dim_subspace,ir)= hopping(3,2)
!   Hamiltonian(3,3+dim_subspace,ir)= hopping(3,3)
!PRINT*, 'AVANT INIT'



  coord(:,:)=ZERO
  DO iwan=1, dimwann
     coord(1,iwan)=iwan*1.0001
     coord(3,iwan)=1.0000000
  ENDDO


   CALL file_open( unit, TRIM(filename), PATH="/", &
                   ACTION="write", FORM="formatted" )


   CALL iotk_write_begin(unit,"HAMILTONIAN")


   !
   CALL iotk_write_attr(attr,"dimwann",dimwann)
   CALL iotk_write_attr(attr,"nr",nr_aux)
   CALL iotk_write_attr(attr,"nrtot",nrtot)
   CALL iotk_write_empty(unit,"DATA",ATTR=attr, IERR=ierr)
      IF (ierr/=0) CALL errore(subname, 'writting DATA', ABS(ierr) )


   !
PRINT*, 'AVANT IVR'
   CALL iotk_write_attr(attr,"units","crystal",FIRST=.TRUE.)
   CALL iotk_write_dat(unit,"IVR", ivr, ATTR=attr, COLUMNS=3) 


   !
   ! get the desired R indexes
   !
PRINT*, 'AVANT HAM'

   CALL iotk_write_begin(unit,"RHAM")
   DO ir=1,nrtot
      CALL iotk_write_dat(unit,"VR"//TRIM(iotk_index(ir)),  Hamiltonian(:,:,ir))
   ENDDO
   CALL iotk_write_end(unit,"RHAM")
   CALL iotk_write_dat(unit,"WANCENTER", coord, ATTR=attr, COLUMNS=3, IERR=ierr) 
      IF (ierr/=0) CALL errore(subname, 'writting coord', ABS(ierr) )

   CALL iotk_write_end(unit,"HAMILTONIAN")
   CALL file_close( unit, PATH="/", ACTION="write" )

!
! cleaning local workspace
!
   DEALLOCATE( ivr, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating ivr', ABS(ierr) )
   DEALLOCATE( coord, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating coord', ABS(ierr) )
   DEALLOCATE(Hamiltonian, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating Hamiltonian', ABS(ierr) )
   DEALLOCATE(hopping, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating Hopping', ABS(ierr) )
   DEALLOCATE(onsite, STAT=ierr)
      IF (ierr/=0) CALL errore(subname, 'deallocating onsite', ABS(ierr) )


   END PROGRAM util_input


