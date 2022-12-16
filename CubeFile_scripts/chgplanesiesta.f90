! This program is to read the charge or potential data from RHO or VT files
!    and calculate the average values in planes. For SIESTA
! INPUT: RHO or VT      -- read charge or local potential
!        planechg.in    -- which file to read, the orientation of the plane
! OUTPUT: chgplane.dat  -- the average values of the data in each plane
! Shenyuan 02-12-2009 at LBL

PROGRAM avechgsiesta
IMPLICIT NONE
REAL(8) :: alatt,vec(3,3),vol,tt(3),lvec,area
REAL, ALLOCATABLE :: chg(:,:,:,:),chgplane(:,:),chgave(:,:,:),chgsum(:)
	! chgave(iz,iave,is),posave(iz,iave,is)
REAL(8), ALLOCATABLE :: lave(:),pos(:),posave(:,:)
REAL :: ryd
INTEGER :: i,j,k,plane,scale,average,lsub,ladd,iposave,iave
INTEGER :: ngrid(3),atom(10),natom,np1,np2,nng,ns,is,lxy,lryd, check
INTEGER, ALLOCATABLE :: lnn(:)
CHARACTER(30) :: Fin,Fop,Fout,A
CHARACTER(60) :: atoms

! ****************************   INITIAL   ******************************

alatt=0.52917720859d0	! Bohr radius, convert the position to Angstrom
ryd=13.60569193		! 1 Ryd = 13.61 eV
Fin='planechg.in'
Fout='planechg.dat'
OPEN(10,FILE=Fin)

READ(10,'(A)')Fop		! which file to be read: CHG or LOCPOT
Fop=TRIM(Fop)
READ(10,'(A)')Fout		! which file to be written: 
Fout=TRIM(Fout)
OPEN(11,FILE=Fout)
READ(10,*)plane		! the orientation of the plane, 1-x,2-y,3-z
READ(10,*)scale		! scale the data by 1/vol or not, 1-yes, other-no
READ(10,*)lxy		! average in xy plane: 1-yes,other-no
READ(10,*)lryd		! multiply data by ryd: 1-yes,other-no
READ(10,*)average		! whether to average over z: 0-no, n-# of average
IF(average > 0) THEN
	ALLOCATE(lave(average))
	ALLOCATE(lnn(average))
	READ(10,*)lave(1:average)
END IF 

OPEN(unit=12,FILE=Fop,status='old',form='formatted')
PRINT*, "Opened"
READ(12,*) vec
PRINT*, "vec", vec 
READ(12,*) ngrid,ns
PRINT*, " ngrid,ns", ngrid,ns

ALLOCATE(chg(ngrid(1),ngrid(2),ngrid(3),ns))
ALLOCATE(chgplane(ngrid(plane),ns))
ALLOCATE(pos(ngrid(plane)))
ALLOCATE(chgsum(ns))
IF(average > 0) THEN
	ALLOCATE(chgave(ngrid(plane),average,ns))
	ALLOCATE(posave(ngrid(plane),average))
END IF
 check = 0 
DO is=1,ns
  PRINT*, "is = ", is
	DO k=1,ngrid(3)
              ! PRINT*, "ik = ", k
		DO j=1,ngrid(2)
                     !   PRINT*, "j = ", j
                        DO i=1,ngrid(1)
			READ(12,*) chg(i,j,k,is)

                        ENDDO
		END DO
	END DO
END DO
PRINT*, check
IF(lryd == 1)chg=chg*ryd	! convert Ryd to eV
!WRITE(11,*)SUM(chg)

	np1=mod(plane+1,3)		! determine the vectors in plane
	IF(np1 == 0) np1=3
	np2=mod(plane+2,3)
	IF(np2 == 0) np2=3

        PRINT*, "Average over planes ", np1, np2
	lvec=alatt*SQRT(vec(plane,1)**2+vec(plane,2)**2+vec(plane,3)**2)
        PRINT*, "lvec = ", lvec
	nng=ngrid(np1)*ngrid(np2)
        PRINT*, "nng = ", nng
 
                                        check = 0 
DO is=1,ns
	DO i=1,ngrid(plane)		! calculate the sum and average in plane
		pos(i)=lvec*DBLE(i-1)/DBLE(ngrid(plane))
		chgplane(i,is)=0.0d0
		do j=1,ngrid(np1)
			do k=1,ngrid(np2)
				IF(plane == 1) THEN
					chgplane(i,is)=chgplane(i,is)+chg(i,j,k,is)
				ELSE IF(plane == 2) THEN
					chgplane(i,is)=chgplane(i,is)+chg(k,i,j,is)
				ELSE
					chgplane(i,is)=chgplane(i,is)+chg(j,k,i,is)
                                      check=check+1
				END IF
			END DO
		END DO
	END DO
END DO

PRINT*, check
!WRITE(11,*)SUM(chgplane)
IF(lxy == 1)chgplane=chgplane/REAL(nng)
!WRITE(11,*)SUM(chgplane)
IF(scale == 1) THEN
        PRINT*, "Code is scaling by volume"
	CALL crossprod(vec(np1,:),vec(np2,:),tt)
	CALL dotprod(tt,tt,area)
	area=SQRT(area)/DBLE(nng)/alatt
	chgplane=chgplane*area
END IF

IF(average > 0)THEN		! calculate the average value within lave
	DO iave=1,average
                PRINT*, "The code is averaging along the main axis"
		lnn(iave)=NINT(lave(iave)/lvec*DBLE(ngrid(plane)))
		chgsum=0.0d0
		DO i=1,lnn(iave)
			DO is=1,ns
				chgsum(is)=chgsum(is)+chgplane(i,is)
			END DO
		END DO
		DO i=1,ngrid(plane)
			lsub=i
			ladd=i+lnn(iave)
			iposave=i+(lnn(iave)-1)/2
			IF(ladd > ngrid(plane)) ladd=ladd-ngrid(plane)
			IF(iposave > ngrid(plane)) iposave=iposave-ngrid(plane)
			IF(mod(lnn(iave),2) == 0) THEN
				IF(iposave == ngrid(plane)) THEN
					posave(iposave,iave)=pos(iposave)+lvec/DBLE(ngrid(plane))/2.0d0
				ELSE
					posave(iposave,iave)=(pos(iposave)+pos(iposave+1))/2.0d0
				END IF
			ELSE
				posave(iposave,iave)=pos(iposave)
			END IF
			DO is=1,ns
				chgave(iposave,iave,is)=chgsum(is)/REAL(lnn(iave))
				chgsum(is)=(chgplane(ladd,is)-chgplane(lsub,is))+chgsum(is)
			END DO
		END DO
	END DO
END IF

IF(average == 0) THEN
	!WRITE(11,*)ngrid(np1),ngrid(np2),ngrid(plane),ns
	DO is=1,ns
		WRITE(11,200)(i,pos(i),chgplane(i,is),i=1,ngrid(plane))
	END DO
ELSE
	WRITE(11,*)ngrid(np1),ngrid(np2),ngrid(plane),ns,(lnn(iave),iave=1,average)
	DO is =1,ns
		DO i=1,ngrid(plane)
			WRITE(11,100)i,pos(i),chgplane(i,is),((posave(i,iave),chgave(i,iave,is)),iave=1,average)
		END DO
		WRITE(11,'(/)')
	END DO
END IF

DEALLOCATE(chg)
DEALLOCATE(chgplane)
DEALLOCATE(pos)
IF(average > 0) THEN
	DEALLOCATE(chgave)
	DEALLOCATE(posave)
END IF

100	FORMAT(I5,F11.6,E20.11,6(F11.6,E20.11))
200	FORMAT(I5,F11.6,E20.11)

CLOSE(10)
CLOSE(11)
CLOSE(12)
STOP
END PROGRAM avechgsiesta

!  ****  SUBROUTINE dotproduce  ****
SUBROUTINE dotprod(a,b,d)
IMPLICIT NONE
DOUBLE PRECISION :: a(3),b(3),d
d=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
RETURN
END SUBROUTINE dotprod

!  ****  SUBROUTINE crossprod  ****
SUBROUTINE crossprod(a,b,c)
IMPLICIT NONE
DOUBLE PRECISION :: a(3),b(3),c(3)
c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)
RETURN
END SUBROUTINE crossprod
