PROGRAM dipoleFromCube
  IMPLICIT none

!indexes begin at 1 because this is fortran
INTEGER :: atoms,i,ix,iy,iz,gridsize=0,lineaxis=0,planesize=0,negCount=0,j=0
REAL*8 :: rVecCent(3)
INTEGER :: axisindex(3)
INTEGER :: voxelcount(3)
REAL*8:: xcoord,ycoord
REAL*8:: tempVector(3)
REAL*8:: rVector(3)
REAL*8:: dipole(3)
REAL*8:: dipoleCenterOfCharge(3)
REAL*8:: ionDipoleZero(3),ionDipoleCenterOfCharge(3)
REAL*8:: totalIonCharge
REAL*8:: voxelvector(3,3)
REAL*8:: totalcharge=0
REAL*8:: centerOfCharge(3)
REAL*8:: totalMass=0
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: grid
REAL*8, ALLOCATABLE :: atomData(:,:)
REAL*8:: total=0,vVoxel=0,nFac=0
INTEGER :: isChargeDensity=0
CHARACTER *256 :: BUFFER,psiIn
CHARACTER *256 :: OutLine

CALL GETARG(1, BUFFER)
READ (BUFFER,'(a256)') psiIn

!skip the first two line

!write(6,*) 'Tring to read ', psiIn

open(100, FILE=psiIn)

read(100,*)
read(100,*)
read(100,*) , atoms

read(100,*)  (voxelcount(i),voxelvector(i,1),voxelvector(i,2),&
     &voxelvector(i,3), i=1, 3)
do i=1,3
write(6,*) voxelcount(i),voxelvector(i,1),voxelvector(i,2),&
     &voxelvector(i,3)
end do

centerOfCharge=0

allocate(atomData(atoms,6))

DO i=1,atoms,1
   read(100,*) (atomData(i,j),j=1,5)
   atomData(i,6) = atomData(i,2)*2 !roughly the mass

   if(atomData(i,2).gt.1.and.atomData(i,2).lt.11) then !fix because we are using pseudopotentials, very dirty
      atomData(i,2)=atomData(i,2)-2
   else if(atomData(i,2).gt.10.and.atomData(i,2).lt.19) then
      atomData(i,2)=atomData(i,2)-10
   endif

   DO j = 1,3
      centerOfCharge(j) = centerOfCharge(j) + atomData(i,j+2) * atomData(i,2)
   END DO
   totalMass= totalMass + atomData(i,6)

   totalIonCharge = totalIonCharge + atomData(i,2)   
END DO
centerOfCharge = centerOfCharge/totalIonCharge 


gridsize=voxelcount(1)*voxelcount(2)*voxelcount(3)
allocate(grid(voxelcount(1),voxelcount(2),voxelcount(3)))

!But the gaussian file is always x,y,z inner loop 
!so we may as well read the whole thing in
read(100,*) (((grid(ix,iy,iz),iz = 1,voxelcount(3)),&
      &iy=1,voxelcount(2)),ix=1,voxelcount(1))

!now check the norm
!grid=ABS(grid)

total=sum(grid)
nFac=1/total !assume single occupany of each state, norm is actually -1

write(6,*) 'Grid Total :', total

!getting the volume per voxel don't assume anything other than the voxels
!are parallel pipeds
tempVector(1)=voxelvector(2,2)*voxelvector(3,3)&
     &-voxelvector(2,3)*voxelvector(3,2)
tempVector(2)=voxelvector(2,3)*voxelvector(3,1)&
     &-voxelvector(2,1)*voxelvector(3,3)
tempVector(3)=voxelvector(2,1)*voxelvector(3,2)&
     &-voxelvector(2,2)*voxelvector(3,1)
DO i=1,3,2
   vVoxel = vVoxel + voxelvector(1,i)*tempVector(i)
END DO

write(6,*) 'Volume per Voxel :', vVoxel
write(6,*) 'Total Charge Density:', vVoxel * total

rVecCent(:) = voxelVector(1,:)*(voxelcount(:)) / 2 &
     + voxelVector(2,:)*(voxelcount(:)) / 2 &
     + voxelVector(3,:)*(voxelcount(:)) / 2

write(6,*) 'Center of Grid: ', rVecCent

do iz = 1,voxelcount(3)
   do iy = 1,voxelcount(2)
      do ix = 1,voxelcount(1)
         rVector(1) = voxelvector(1,1)*(ix-1) + &
              voxelvector(2,1)*(iy-1) + &
              voxelvector(3,1)*(iz-1) 
         rVector(2) = voxelvector(1,2)*(ix-1) + &
              voxelvector(2,2)*(iy-1) + &
              voxelvector(3,2)*(iz-1) 
         rVector(3) = voxelvector(1,3)*(ix-1) + &              
              voxelvector(2,3)*(iy-1) + &
              voxelvector(3,3)*(iz-1) 
         dipole(:) = dipole(:) + -ABS(grid(ix,iy,iz)) * vVoxel * rVector(:) 
         totalcharge = totalcharge + -ABS(grid(ix,iy,iz)) * vVoxel
      end do
   end do
end do

if(ABS(totalcharge - totalIonCharge).lt.0.001) then
   isChargeDensity = 1
end if


write(6,*) ' total charge: ', totalcharge
write(6,*) ' total ion charge: ', totalIonCharge
write(outLine,*) ' dipole: ', dipole
write(6,*) outLine
!write(6,*) ' dipole total: ', dipole*totalIonCharge

dipoleCenterOfCharge(:) = dipole(:) - totalcharge * centerOfCharge(:)

write(6,*) ' Center of Charge: ', centerOfCharge
!write(6,*) ' dipole Center of Grid: ', dipoleCenterOfGrid
!write(6,*) ' dipole Center of Grid total: ', dipoleCenterOfGrid*totalIonCharge

write(outLine,*) ' dpCC: ', dipoleCenterOfCharge
write(6,*) outLine
DO i=1,atoms,1
   do j = 1,3 
      ionDipoleZero(j) = ionDipoleZero(j) + atomData(i,2)*(atomData(i,j+2))
   end do   
END DO

ionDipoleCenterOfCharge(:) = ionDipoleZero(:) - totalIonCharge * centerOfCharge(:)

write(6,*) 'Total Ion Charge:', totalIonCharge
write(6,*) 'Ions Dipole Center Of Charge:', ionDipoleCenterOfCharge
!write(6,*) 'Ions Dipole Center of Grid:', ionDipoleCenterOfGrid
!write(6,*) 'Ion Dipole Corner of Box:', ionDipoleZero

if(isChargeDensity.eq.1) then
   write(6,*) 'Corner Total Dipole:', ionDipoleZero + dipole
   write(6,*) 'Center of Charge Total Dipole:', ionDipoleCenterOfCharge + dipoleCenterOfCharge
end if


!output for plotting

iy = 36
iz = 36

open(101,FILE="plotLineOut")
do i = 1,voxelcount(1)
  xcoord = voxelvector(1,1)*(i-.5) + &
              voxelvector(2,1)*(iy-.5) + &
              voxelvector(3,1)*(iz-.5)
  write(101,*) xcoord , grid(i,iy,iz)
end do
close(101)

open(101,FILE="plotPlaneOut")
do i = 1,voxelcount(1)
   do j = 1,voxelcount(2)
      xcoord = voxelvector(1,1)*(i-.5) + &
           voxelvector(2,1)*(j-.5) + &
           voxelvector(3,1)*(iz-.5)
      ycoord = voxelvector(1,2)*(i-.5) + &
           voxelvector(2,2)*(j-.5) + &
           voxelvector(3,2)*(iz-.5)
      write(101,*) xcoord, ycoord, grid(i,j,iz)
   end do
end do

END PROGRAM dipoleFromCube


