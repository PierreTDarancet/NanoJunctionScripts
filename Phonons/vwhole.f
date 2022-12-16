
      module constants
       integer, parameter :: dp = kind(1.0d0)
       integer, parameter :: dpc = kind((1.0d0,1.0d0))
      end module constants



!program for creating and diagonalizing force constant matrix
!it needs the "forces_minus" and "forces_plus" file which 
!contain forces for all negative and positive
!displacements ..... 

      Program vwhole

      use constants

      implicit none


       TYPE :: atom
        CHARACTER(15) :: name    ! name of object
        REAL :: mass 
        REAL, Dimension(3) :: r ! coordinates
        REAL, Dimension(3) :: f
        INTEGER ::  status = 0      ! 0 if inside, -1 if outside
       END TYPE atom


       TYPE :: name
        CHARACTER(15) :: name    ! name of atom
       END TYPE name


       TYPE(name), dimension(:), allocatable :: AName 
       REAL(dp), dimension(:,:,:,:), allocatable :: Dm
       Real, dimension(:), allocatable :: mass


       Integer :: NATypes, Ntotal, Acounter, PhCount, Ndisp
       Integer :: MAcounter, PAcounter, Write_counter
       Integer :: DNumber, Disp_order, zz, tg, tz, ff, FCount
       CHARACTER(25) :: STR_arg, File_n1, File_n2, F_name, M_name
       CHARACTER(25) :: POS_name
       CHARACTER(25) :: STR_temp
       Real, Dimension(3) :: cor 
       Integer, dimension(:), allocatable ::  Natom 
       TYPE(atom), dimension(:), allocatable :: Atoms, Sorted 
       TYPE(atom), dimension(:), allocatable :: MAtoms, PAtoms
       TYPE(atom), dimension(:), allocatable :: TempAtoms1, TempAtoms2 

       


       Real, Dimension(3) :: RX = (/ -1.0, 1.0, 1.0 /)
       Real, Dimension(3) :: RY = (/ 1.0, -1.0, 1.0 /)
       Real, Dimension(3) :: Atemp0, Atemp1
       Real, Dimension(3) :: Ashift, Myshift
       Real :: shift
       Real(dp), Dimension(3) :: LatPar


       Integer :: hh, Zcount, vv, Tcount, Status_sum
       Integer :: moved_atom, AddCount, Nentries, NewCount
       Integer :: i, j, k, t, l, d, h, SCount, Nat, Mode, N_new_disp
       Integer :: N_arg, EScount
       REAL(dp) :: A1, A2, B1, B2, Conv
       REAL :: m1, m2, joint_mass
       Integer :: alpha, beta, ai, bi


        Real(dp), parameter :: Disp = 0.01889726 ! in a.u.!to be changed
        Real(dp), parameter :: tol = 0.001
        Real(dp), parameter :: ScaleAmp = 0.2
        Real(dp) :: Moved
      
        Integer :: x,y
        Integer, parameter :: N = 9 !size of the dynamical matrix!to be changed
        Integer, parameter :: LWork = (64+2)*N
        Integer, parameter :: Nmax = 10
        Integer, parameter :: LDA = N
        Real(dp), dimension(LWork) :: Work
	Real(dp), dimension(N,N) ::  A,fB ! the matrix
        Real(dp), dimension(N) :: b, c, eigv
        Integer, dimension(N) :: pivot
	integer :: kk, Info, jj, found




      
        


!
! reading arguments
!

       N_arg = IARGC()          ! reads the number of arguments

       if (N_arg .ne. 3) then
          write(*,*) 'this program needs 3 arguments'
          write(*,*) N_arg, 'arguments were entered'
          write(*,*) 'enter number of species and'
          write(*,*) 'total number of atoms'
          write(*,*) 'number of calculated configurations'
          write(*,*) '------------------------------------'
          write(*,*)
          write(*,*) 'Usage'
          write(*,*) 'unfold.exe reads coordinates and forces'
          write(*,*) 'from forces'
          write(*,*) 'in a.u.'
          write(*,*) 'generates a complete forces list'
          write(*,*) 'mode number'
          write(*,*)
          write(*,*) 'have fun!'
          stop
       end if


       CALL GETARG(1, STR_arg) ! returns the argument 1
!       write(*,*) STR_arg
        read(STR_arg,*) NATypes
!       write(*,*) NATypes, ' types of atoms'
       CALL GETARG(2, STR_arg) ! returns the argument 2
!       write(*,*) STR_arg
        read(STR_arg,*) Nat
!       write(*,*) Nat, ' atoms'
       CALL GETARG(3, STR_arg) ! returns the argument 3
        read(STR_arg,*) Nentries


    

   

       Ndisp = Nentries ! all_atoms*3*2

       Ntotal = Nat*Nat*3+Nat ! total number of forces
       ! 30 entries
       ! i.e. Ndisp * number of atoms (5)
       ! plus one more entry, which is
       ! the initial configuration (5 atoms )

       allocate(AName(NATypes)) ! array of names
       allocate(Natom(NATypes))
       allocate(Atoms(Ntotal))  ! allocate array for atoms 
       allocate(MAtoms(Ntotal))  ! allocate array for atoms 
       allocate(PAtoms(Ntotal))  ! allocate array for atoms 

       allocate(Sorted(Ntotal))  ! allocate array for atoms 
       allocate(mass(NATypes))
       allocate(Dm(Nat,3,Nat,3))
       allocate(TempAtoms1(Nat))
       allocate(TempAtoms2(Nat))





!        write(*,*) 'enter names,  masses and # of each species:'
!        do i=1, NATypes, 1
!           read(5,*)  AName(i)%name, mass(i), Natom(i)
!        enddo   



!        AName(1)%name = 'P'
!        AName(2)%name = 'Si'
!        AName(3)%name = 'H'
!        mass(1) = 30.97
!        mass(2) = 28.09
!        mass(3) = 1.008
!        Natom(1) = 1
!        Natom(2) = 28
!        Natom(3) = 36


!        AName(1)%name = 'Si'
!        AName(2)%name = 'H'
!        mass(1) = 28.09
!        mass(2) = 1.008
!        Natom(1) = 87 
!        Natom(2) = 76




!        Natom(1) = 87 
!        Natom(2) = 76

!        AName(1)%name = 'Si'
!        mass(1) = 28.09
!        Natom(1) = 6


        AName(1)%name = 'C'
        AName(2)%name = 'O'
        mass(1) = 14.006 
        mass(2) = 15.9994
        Natom(1) = 1 
        Natom(2) = 2


!      LatPar(1) = 10.2165184020996094
!      LatPar(2) = 11.7970190048217773
!      LatPar(3) = 25.0000000000000000


        open (UNIT=11,FILE='forces_minus',STATUS='OLD',action="read",
     1         ERR=190)
        open (UNIT=12,FILE='forces_plus',STATUS='OLD',action="read",
     1         ERR=190)


        open(UNIT=30, file='out.dat',form='formatted',status='unknown')






!----------------------------------------------------
! reading the zero displacement configuration ...
! the first entry
          MAcounter = 0
!..............
!          read(11,*)
          read(11,*)
          do t=1, NATypes, 1
            do i=1, Natom(t), 1
              MAcounter = MAcounter + 1
              read(11,*) MAtoms(MAcounter)%r(:),
     1                   MAtoms(MAcounter)%f(:)
              MAtoms(MAcounter)%name = AName(t)%name
              MAtoms(MAcounter)%mass = mass(t)
            end do !i
          end do !t
!          read(11,*)
!          read(11,*)

!........................
          PAcounter = 0
!          read(12,*)
          read(12,*)
          do t=1, NATypes, 1
            do i=1, Natom(t), 1
              PAcounter = PAcounter + 1
              read(12,*) PAtoms(PAcounter)%r(:),
     1                   PAtoms(PAcounter)%f(:)
              PAtoms(PAcounter)%name = AName(t)%name
              PAtoms(PAcounter)%mass = mass(t)
            end do !i
          end do !t
!          read(12,*)
!          read(12,*)
! ..................................


          Acounter = 0
          do t=1, NATypes, 1
            do i=1, Natom(t), 1
              Acounter = Acounter + 1
              Atoms(Acounter)%status=t
              Atoms(Acounter)%r(:)=PAtoms(Acounter)%r(:)*1.889726
              Atoms(Acounter)%f(:)= 
     1            (MAtoms(Acounter)%f(:)-PAtoms(Acounter)%f(:))/2.0
              Atoms(Acounter)%name = PAtoms(Acounter)%name
              Atoms(Acounter)%mass = mass(t)
            end do !i
          end do !t



!          write(*,*)'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'



!------------------------------------------------------------------
          AddCount = 0 ! counter of additional configurations
          NewCount = 0


! reading all other displacements ...        
        do d=1, Ndisp, 1 ! over all displacements
!          read(11,*)
          read(11,*)

!          read(12,*)
          read(12,*)

          write(*,*)
          SCount = 0 ! counter of atoms in each configuration


!          write(*,*) 'reading configuration', d
! read configuration
          do t=1, NATypes, 1 ! over all atoms types
            do i=1, Natom(t), 1 ! over the type
              Acounter = Acounter + 1
              SCount = SCount + 1

              read(11,*) MAtoms(Acounter)%r(:),
     1                   MAtoms(Acounter)%f(:)
              MAtoms(Acounter)%name = AName(t)%name
              MAtoms(Acounter)%mass = mass(t)
              MAtoms(Acounter)%status = t

              read(12,*) PAtoms(Acounter)%r(:),
     1                   PAtoms(Acounter)%f(:)
              PAtoms(Acounter)%name = AName(t)%name
              PAtoms(Acounter)%mass = mass(t)
              PAtoms(Acounter)%status = t


              Atoms(Acounter)%status=t
              Atoms(Acounter)%r(:)=PAtoms(Acounter)%r(:)
              Atoms(Acounter)%f(:)= 
     1            (MAtoms(Acounter)%f(:)-PAtoms(Acounter)%f(:))/2.0
              Atoms(Acounter)%name = PAtoms(Acounter)%name
              Atoms(Acounter)%mass = PAtoms(Acounter)%mass
              

              write(*,'(a2,3g12.5,3g15.7)') Atoms(Acounter)%name, 
     1            Atoms(Acounter)%r(:), Atoms(Acounter)%f(:)

              ! subtract the initial configuration, so that
              ! we get displacements
              Ashift(:) =  Atoms(Acounter)%r(:)-Atoms(SCount)%r(:)
              shift = sqrt(Ashift(1)*Ashift(1)+
     1                    Ashift(2)*Ashift(2)+
     2                    Ashift(3)*Ashift(3))
              

              ! determine which atom was shifted
              if (shift .ge. 0.009) then
                 ! this is our atom
                 Myshift = Ashift
                 !write(*,*) ' => ',Ashift(:)
                 moved_atom = SCount
              end if   
              
              
            end do               !i  
          end do !t
!          read(11,*)
!          read(11,*)

!          read(12,*)
!          read(12,*)

!---------------------------------------------------------


          end do ! d










       

       do i=1, Ntotal, 1
        Atoms(i)%f(:)=Atoms(i)%f(:)*4.119/1.889726/13.36054 ! Ry/a.u. to Newtons  
!        write(30,'(a2,i4,3g15.5,3g15.5)') Atoms(i)%name, 
!     1        Atoms(i)%status,Atoms(i)%r(:), Atoms(i)%f(:)
!        write(30,*) 'i = ',i
       end do













       


!       m1=28.086*1.66e-3 ! 1.66e-27 - atomic mass to kg
                         ! x 1.0e24 transformation to THz factor
                         ! put here to get rid of the large powers
       Conv=1.0/(5.29e-3)

! collecting force constant matrix ...
       do alpha=1, Nat, 1
        do i=1, 3, 1
          do beta=1, Nat, 1
            do j=1, 3, 1 
            
                    


               ai = Nat+alpha+Nat*3*(beta-1)+Nat*(j-1)
            !   write(*,*) 'ai = ',ai
            A1=Conv*Atoms(ai)%f(i)
            A2=-A1

               bi = Nat+beta+Nat*3*(alpha-1)+Nat*(i-1)
            !   write(*,*) 'bi = ',bi
            B1=Conv*Atoms(bi)%f(j)
            B2=-B1
            
!            write(*,*) 'A1 = ', A1, '     B1 = ',B1, '  Disp = ',Disp

            m1 = Atoms(ai)%mass
            m2 = Atoms(bi)%mass
            joint_mass = sqrt(m1*m2)*1.66e-3
            !joint_mass = 29.08*1.66e-3
!          write(*,*) 'mass = ', joint_mass
               Dm(alpha,i,beta,j) = -0.5*((A2-A1)/(2.0*Disp)+
     1         (B2-B1)/(2.0*Disp))/joint_mass



              end do!j 
          end do!beta
         end do!j  
       end do!alpha   
! .............................


!       do i=1, Nat*3, 1 
!              read(40,*) A(1:18,i)
!       enddo


       y=0
! writing the force matrix out 
       do alpha=1, Nat, 1
            do i=1, 3, 1
               x=0
               y=y+1 ! increment for raws
          do beta=1, Nat, 1
              do j=1, 3, 1 
               x=x+1 ! increment for columns
               A(x,y)=Dm(alpha,i,beta,j)
!              write(*,'(2f10.4)') Dm(alpha,i,beta,j), Dm(beta,j,alpha,i) 

              end do!j
            end do!beta 
          end do!i  
       end do!alpha   


!       do i=1, Nat*3, 1 
!              read(40,*) fB(1:18,i)
!       enddo

!        do i=1, N
!           write(*,'(18f12.4)') A(1:18,i)
!        end do
!        write(*,*)
!        do i=1, N
!           write(*,'(18f12.4)') fB(1:18,i)
!           write(*,*)
!        end do

!        A(:,4:18)=fB(:,4:18)


!        do i=1, N
!           write(*,'(18f12.4)') A(1:18,i)
!        end do
!        write(*,*)
!        do i=1, N
!           write(*,'(18f12.4)') fB(1:18,i)
!           write(*,*)
!        end do



! ------------------------------------

!Diagonalize this matrix
! to find all the eigenvalues and eigenvectors of the symmetric matrix
!'N':  only igenvalues, 'V' iegenvectors too
!'U':  Upper triangle of A is stored
! N : order of the matrix
! A : the matrix
! Work : workspace array
! LWork : size of Work
! Info : =0 if successful

	call DSYEV('V','U',N,A,N,c,Work,LWork,Info)


	do i=(1+6), N ! skip first six modes
           if (c(i) .lt. 0.0 ) then 
             eigv(i)=-sqrt(abs(c(i)))
            else
             eigv(i)=sqrt(c(i))
            endif
	   write(30,'(3f15.5)') 1.0, eigv(i)*33.3566/2.0/3.1415
	end do







        do Mode=7, Nat*3, 1
          write(M_name,'(i10)') Mode 
          write(F_name,'(f12.2)') eigv(Mode)*33.3566/2.0/3.1415
          write(*,*) F_name
          File_n1='m_'// TRIM(ADJUSTL(M_name))// '_' //
     1        TRIM(ADJUSTL(F_name))// '.xyz'
          open(UNIT=37,file=File_n1,form='formatted',status='unknown')
!          File_n2='froz_in.'// TRIM(ADJUSTL(M_name))
!          open(UNIT=38,file=File_n2,form='formatted',status='unknown')
! POSCARs
        POS_name = 'Conf_m' // TRIM(ADJUSTL(M_name))
        open (UNIT=41, FILE=POS_name, STATUS='unknown',
     1         form='formatted',action="write",ERR=190)
        POS_name = 'Conf_p' // TRIM(ADJUSTL(M_name))
        open (UNIT=42, FILE=POS_name, STATUS='unknown',
     1         form='formatted',action="write",ERR=190)



          
! writing mimus configuration
!---------------------------------------------------------------

        Write_counter = 0
        do i=1, (Nat), 1

          Write_counter = Write_counter + 1 
            write(41,'(3f14.8, i4,"  ",a2,"  ",i4)')  
     1           (Atoms(Write_counter)%r(1)*0.52917-4.44-
     2           A((1+3*(Write_counter-1)),Mode)*ScaleAmp*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(2)*0.52917-5.485-
     2           A((2+3*(Write_counter-1)),Mode)*ScaleAmp*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(3)*0.52917-7.09-
     2           A((3+3*(Write_counter-1)),Mode)*ScaleAmp*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           Atoms(Write_counter)%status,
     2           Atoms(Write_counter)%name,
     3           Write_counter
        end do
!--------------------------------------------------------------
! writing plus configuration
!---------------------------------------------------------------

        Write_counter = 0
        do i=1, (Nat), 1

          Write_counter = Write_counter + 1 
            write(42,'(3f14.8, i4,"  ",a2,"  ",i4)')  
     1           (Atoms(Write_counter)%r(1)*0.52917-4.44+
     2           A((1+3*(Write_counter-1)),Mode)*ScaleAmp*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(2)*0.52917-5.485+
     2           A((2+3*(Write_counter-1)),Mode)*ScaleAmp*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(3)*0.52917-7.09+
     2           A((3+3*(Write_counter-1)),Mode)*ScaleAmp*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           Atoms(Write_counter)%status,
     2           Atoms(Write_counter)%name,
     3           Write_counter
        end do
!--------------------------------------------------------------

        Moved=0.0
           
        write(37,*) Nat
        write(37,*)
        do i=1, Nat, 1
            write(37,'(a2,6f10.4)') Atoms(i)%name, 
     1           Atoms(i)%r(:)*0.52917,
     2           A((1+3*(i-1)):(3*i),Mode)*0.52917/sqrt(Atoms(i)%mass)

            Moved = Moved + (A((1+3*(i-1)),Mode)*0.52917/
     1       sqrt(Atoms(i)%mass))**2.0 +
     2        (A((2+3*(i-1)),Mode)*0.52917/
     3       sqrt(Atoms(i)%mass))**2.0 +
     4       (A((3+3*(i-1)),Mode)*0.52917/
     5       sqrt(Atoms(i)%mass))**2.0



        end do
         write(30,'(F14.7)') sqrt(Moved)






 
        write(37,*)'#Energy ', eigv(Mode)*33.3566/2.0/3.1415
        close(42)
        close(41)
        close(37)
        enddo




!-----------------------------------------
        goto 200

 190    write(*,*) 'error when opening file'

 200    continue





        deallocate(TempAtoms2)
        deallocate(TempAtoms1)
        deallocate(Dm)
        deallocate(mass)
        deallocate(Sorted)
        deallocate(Atoms)
        deallocate(MAtoms)
        deallocate(PAtoms)
        deallocate(AName)
        deallocate(Natom)


!..........................................

 10     format(a2, 3f10.5)
 11     format(i4, '  ', a2, 3f10.5, i3)
 12     format(a16, a2, 3f10.5)



        close(30)
        close(12)
        close(11)



        stop

        End program vwhole
