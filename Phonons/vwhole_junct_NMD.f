       
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
        INTEGER ::  status = 0  
        INTEGER :: type = 0    
       END TYPE atom


       TYPE :: name
        CHARACTER(15) :: name    ! name of atom
       END TYPE name


       TYPE(name), dimension(:), allocatable :: AName 
       REAL(dp), dimension(:,:,:,:), allocatable :: Dm
       Real, dimension(:), allocatable :: mass
       INTEGER, dimension(:), allocatable :: AType 


       Integer :: NatTotal, StaticAtN
       Integer :: NATypes, Ntotal, Acounter, PhCount, Ndisp
       Integer :: MAcounter, PAcounter, Write_counter, Stcounter
       Integer :: DNumber, Disp_order, zz, tg, tz, ff, FCount

       CHARACTER(25) :: STR_arg, File_n1, File_n2, F_name, M_name
       CHARACTER(25) :: POS_name
       CHARACTER(25) :: File_n3
       CHARACTER(25) :: STR_temp
       Real, Dimension(3) :: cor 
       Integer, dimension(:), allocatable ::  Natom 
       TYPE(atom), dimension(:), allocatable :: Atoms, StaticAt
       TYPE(atom), dimension(:), allocatable :: MAtoms, PAtoms
       Real(dp), Dimension(3) :: LatPar

       


       Real, Dimension(3) :: Atemp0, Atemp1
       Real, Dimension(3) :: Ashift, Myshift
       Real :: shift


       Integer :: hh, Zcount, vv, Tcount, Status_sum
       Integer :: moved_atom, AddCount, Nentries, NewCount
       Integer :: i, j, k, t, l, d, h, SCount, Nat, Mode, N_new_disp
       Integer :: N_arg, EScount
       REAL(dp) :: A1, A2, B1, B2, Conv
       REAL :: m1, m2, joint_mass
       Integer :: alpha, beta, ai, bi


        Real(dp), parameter :: Disp =  0.02834589   !0.01889726  ! 0.03779452 ! in a.u.
        Real(dp), parameter :: tol = 0.001
        Real(dp), parameter :: ScaleAmp = 1.0
      
        Integer :: x,y
        Integer, parameter :: N = 90 !size of the dynamical matrix
        Integer, parameter :: LWork = (64+2)*N
        Integer, parameter :: Nmax = 10
        Integer, parameter :: LDA = N
        Real(dp), dimension(LWork) :: Work
	Real(dp), dimension(N,N) ::  A,fB ! the matrix
        Real(dp), dimension(N) :: b, c, eigv
        Integer, dimension(N) :: pivot
        Integer :: kk, Info, jj, found
!
! reading arguments
!

       N_arg = IARGC()          ! reads the number of arguments

       if (N_arg .ne. 4) then
          write(*,*) 'this program needs two arguments'
          write(*,*) N_arg, 'arguments were entered'
          write(*,*) 'enter number of species and'
          write(*,*) 'total number of atoms'
          write(*,*) 'total number of displaced atoms'
          write(*,*) 'number of calculated configurations'
          write(*,*) '------------------------------------'
          write(*,*)
          write(*,*) 'Usage'
          write(*,*) '*.exe reads coordinates and forces'
          write(*,*) 'from forces_minus and forces_plus'
          write(*,*) 'reads static.atoms'
          write(*,*)
          write(*,*) 'have fun!'
          stop
       end if


       CALL GETARG(1, STR_arg) ! returns the argument 1
        read(STR_arg,*) NATypes

       CALL GETARG(2, STR_arg) ! returns the argument 2
        read(STR_arg,*) NatTotal

       CALL GETARG(3, STR_arg) ! returns the argument 2
        read(STR_arg,*) Nat

       CALL GETARG(4, STR_arg) ! returns the argument 3
        read(STR_arg,*) Nentries


    

   

       Ndisp = Nentries ! all_atoms*3*2
       StaticAtN = NatTotal - Nat

       Ntotal = Nat*Nat*3+Nat ! total number of forces
       ! 30 entries
       ! i.e. Ndisp * number of atoms (5)
       ! plus one more entry, which is
       ! the initial configuration (5 atoms )

       allocate(AName(NATypes)) ! array of names
       allocate(AType(NATypes))
       allocate(Natom(NATypes))
       allocate(Atoms(Ntotal))  ! allocate array for atoms 
       allocate(MAtoms(Ntotal))  ! allocate array for atoms 
       allocate(PAtoms(Ntotal))  ! allocate array for atoms 
       allocate(StaticAt(StaticAtN))

       allocate(mass(NATypes))
       allocate(Dm(Nat,3,Nat,3))
!        AName(1)%name = 'P'
!        mass(1) = 30.97
!        AName(1)%name = 'Si'
!        mass(1) = 28.09
!        AName(1)%name = 'Au'
!        mass(1) = 196.97
!        AName(2)%name = 'S'
!        mass(2) = 32.066
!        AName(3)%name = 'C'
!        mass(3) = 12.011
!        AName(4)%name = 'H'
!        mass(4) = 1.008
!        Natom(1) = 4
!        Natom(2) = 1
!        Natom(3) = 6
!        Natom(4) = 5
!      LatPar(1) = 10.2165184020996094
!      LatPar(2) = 11.7970190048217773
!      LatPar(3) = 25.0000000000000000
        AName(1)%name = 'Au'
        AName(2)%name = 'H'
        AName(3)%name = 'S'
        AName(4)%name = 'C'
        AName(5)%name = 'Si'
        mass(1) = 196.97
        mass(2) = 1.008
        mass(3) = 32.066
        mass(4) = 12.011
        mass(5) = 28.09
        AType(1) = 1 
        AType(2) = 2
        AType(3) = 3 
        AType(4) = 4
        AType(5) = 5
        Natom(1) = 38
        Natom(2) = 16
        Natom(3) = 2
        Natom(4) = 6
        Natom(5) = 2



        open (UNIT=11,FILE='cut_minus',STATUS='OLD',action="read",
     1         ERR=190)
        open (UNIT=12,FILE='cut_plus',STATUS='OLD',action="read",
     1         ERR=190)


        open (UNIT=13,FILE='static.atoms',STATUS='OLD',action="read",
     1         ERR=190)



        open(UNIT=30, file='out.dat',form='formatted',status='unknown')
        open(UNIT=45, file='mod.nmd',form='formatted',status='unknown')











!----------------------------------------------------
! reading the zero displacement configuration ...
! the first entry
          MAcounter = 0
!..............
          do t=1, NATypes, 1
            do i=1, Natom(t), 1
              MAcounter = MAcounter + 1
              read(11,*) MAtoms(MAcounter)%r(:),
     1                   MAtoms(MAcounter)%f(:)
              MAtoms(MAcounter)%name = AName(t)%name
              MAtoms(MAcounter)%mass = mass(t)
              MAtoms(MAcounter)%type = AType(t)
            end do !i
          end do !t
          read(11,*)


!........................
          PAcounter = 0
          do t=1, NATypes, 1
            do i=1, Natom(t), 1
              PAcounter = PAcounter + 1
              read(12,*) PAtoms(PAcounter)%r(:),
     1                   PAtoms(PAcounter)%f(:)
              PAtoms(PAcounter)%name = AName(t)%name
              PAtoms(PAcounter)%mass = mass(t)
              PAtoms(PAcounter)%type = AType(t)
            end do !i
          end do !t
          read(12,*)
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
              Atoms(Acounter)%type = AType(t)
            end do !i
          end do !t





!........................
          Stcounter = 0
          do t=1, StaticAtN, 1
              Stcounter = Stcounter + 1
              read(13,*) StaticAt(Stcounter)%r(:),
     1                   StaticAt(Stcounter)%f(:)
          end do !t
! ..................................




          


!------------------------------------------------------------------
          AddCount = 0 ! counter of additional configurations
          NewCount = 0


! reading all other displacements ...        
        do d=1, Ndisp, 1 ! over all displacements

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
              MAtoms(Acounter)%type = AType(t)
              MAtoms(Acounter)%status = t

              read(12,*) PAtoms(Acounter)%r(:),
     1                   PAtoms(Acounter)%f(:)
              PAtoms(Acounter)%name = AName(t)%name
              PAtoms(Acounter)%mass = mass(t)
              PAtoms(Acounter)%type = AType(t)
              PAtoms(Acounter)%status = t


              Atoms(Acounter)%status=t
              Atoms(Acounter)%r(:)=PAtoms(Acounter)%r(:)
              Atoms(Acounter)%f(:)= 
     1            (MAtoms(Acounter)%f(:)-PAtoms(Acounter)%f(:))/2.0
              Atoms(Acounter)%name = PAtoms(Acounter)%name
              Atoms(Acounter)%mass = PAtoms(Acounter)%mass
              Atoms(Acounter)%type = PAtoms(Acounter)%type
              

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
          read(11,*)
          read(12,*)

!---------------------------------------------------------


          end do ! d




          





       

       do i=1, Ntotal, 1
        Atoms(i)%f(:)=Atoms(i)%f(:)*4.119/1.889726/13.36054 
! Ry/a.u. to Newtons   

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


!===================================
! writing NMD file
        Write_counter = 0
        write(45,'(a16)',advance="no") "atomnames  "
        do i=1, (Nat+StaticAtN), 1
        if (i .le. Nat) then 
          Write_counter = Write_counter + 1
          write(45,'(a4)',advance="no")
     1       Atoms(Write_counter)%name
         else
            write(45,'(a4)',advance="no") "Au  "
         endif
        end do
        write(45,'(a1)') " "
!====================================

!===================================
! writing NMD file
        Write_counter = 0
        write(45,'(a16)',advance="no") "coordinates  "
        do i=1, (Nat+StaticAtN), 1
        if (i .le. Nat) then 
          Write_counter = Write_counter + 1
          write(45,'(3f14.8)',advance="no")
     1       Atoms(Write_counter)%r(:)*0.52917
         else
            write(45,'(3f14.8)',advance="no")
     2       StaticAt(i-Write_counter)%r(:)
         endif
        end do
        write(45,'(a1)') " "
!====================================





        do Mode=7, Nat*3, 1
          write(M_name,'(i10)') Mode 
          write(F_name,'(f12.2)') eigv(Mode)*33.3566/2.0/3.1415
          write(*,*) F_name
          File_n1='m_'// TRIM(ADJUSTL(M_name))// '_' //
     1        TRIM(ADJUSTL(F_name))// '.xyz'
          open(UNIT=37,file=File_n1,form='formatted',status='unknown')
          File_n2='froz_in.'// TRIM(ADJUSTL(M_name))
!          open(UNIT=38,file=File_n2,form='formatted',status='unknown')
          File_n3='sl_'// TRIM(ADJUSTL(M_name))// '_' //
     1        TRIM(ADJUSTL(F_name))// '.xyz'
          open(UNIT=39,file=File_n3,form='formatted',status='unknown')
! POSCARs
        POS_name = 'Conf_m' // TRIM(ADJUSTL(M_name))
        open (UNIT=41, FILE=POS_name, STATUS='unknown',
     1         form='formatted',action="write",ERR=190)
        POS_name = 'Conf_p' // TRIM(ADJUSTL(M_name))
        open (UNIT=42, FILE=POS_name, STATUS='unknown',
     1         form='formatted',action="write",ERR=190)

          

           
        write(37,*) Nat
        write(37,*)
        do i=1, Nat, 1
            write(37,'(a2,6f10.4)') Atoms(i)%name, 
     1           Atoms(i)%r(:)*0.52917,
     2           A((1+3*(i-1)):(3*i),Mode)*0.52917*ScaleAmp
        end do




        write(39,*) (Nat+StaticAtN)
        write(39,*)
        Write_counter = 0
        do i=1, (Nat+StaticAtN), 1

        if (i .le. Nat) then
          Write_counter = Write_counter + 1 
            write(39,'(a2,6f14.8)') Atoms(Write_counter)%name, 
     1           Atoms(Write_counter)%r(:)*0.52917,
     2           A((1+3*(Write_counter-1)):(3*Write_counter),Mode)*
     3           0.52917*ScaleAmp
         else
            Write_counter = Write_counter + 1
            write(39,'(a2,6f14.8)') 'Au', 
     1           StaticAt(Write_counter-Nat)%r(:),
     2           0.0,  0.0,  0.0           
         endif   
        end do


!===================================
! writing NMD file
        Write_counter = 0
        write(45,'(a8)',advance="no") "mode  "
        do i=1, (Nat+StaticAtN), 1
        if (i .le. Nat) then 
          Write_counter = Write_counter + 1
          write(45,'(3f14.8)',advance="no")
     1           A((1+3*(Write_counter-1)):(3*Write_counter),Mode)*
     2           0.52917*ScaleAmp   !/sqrt(Atoms(Write_counter)%mass)
         else
            write(45,'(3f14.8)',advance="no")
     2           0.0,  0.0,  0.0
         endif
        end do
        write(45,'(a1)') " "
!====================================



! writing to POSCAR_m
!---------------------------------------------------------------

        Write_counter = 0
        do i=1, (Nat+StaticAtN), 1

        if (i .le. Nat) then
          Write_counter = Write_counter + 1 
            write(41,'(3f14.8, i4,"  ",a2,"  ",i4)')  
     1           (Atoms(Write_counter)%r(1)*0.52917-
     2           A((1+3*(Write_counter-1)),Mode)*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(2)*0.52917-
     2           A((2+3*(Write_counter-1)),Mode)*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(3)*0.52917-
     2           A((3+3*(Write_counter-1)),Mode)*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           Atoms(Write_counter)%type,
     2           Atoms(Write_counter)%name,
     3           Write_counter
         else
            Write_counter = Write_counter + 1
            write(41,'(3f14.8, i4,"  ",a2,"  ",i4)')  
     1           StaticAt(Write_counter-Nat)%r(1),
     2           StaticAt(Write_counter-Nat)%r(2),
     3           StaticAt(Write_counter-Nat)%r(3),
     1           4,
     2           "Au",
     3           Write_counter
         endif   
        end do
!--------------------------------------------------------------
! writing to POSCAR_p
!---------------------------------------------------------------

        Write_counter = 0
        do i=1, (Nat+StaticAtN), 1

        if (i .le. Nat) then 
          Write_counter = Write_counter + 1 
            write(42,'(3f14.8, i4,"  ",a2,"  ",i4)')  
     1           (Atoms(Write_counter)%r(1)*0.52917+
     2           A((1+3*(Write_counter-1)),Mode)*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(2)*0.52917+
     2           A((2+3*(Write_counter-1)),Mode)*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           (Atoms(Write_counter)%r(3)*0.52917+
     2           A((3+3*(Write_counter-1)),Mode)*
     3           0.52917/sqrt(Atoms(Write_counter)%mass)),
     1           Atoms(Write_counter)%type,
     2           Atoms(Write_counter)%name,
     3           Write_counter
         else
            Write_counter = Write_counter + 1
            write(42,'(3f14.8, i4,"  ",a2,"  ",i4)')  
     1           StaticAt(Write_counter-Nat)%r(1),
     2           StaticAt(Write_counter-Nat)%r(2),
     3           StaticAt(Write_counter-Nat)%r(3),
     1           4,
     2           "Au",
     3           Write_counter
         endif   
        end do
!--------------------------------------------------------------





 
        write(37,*)'#frequency cm-1 ', eigv(Mode)*33.3566/2.0/3.1415
        close(42)
        close(41)
        close(37)
!        close(38)
        enddo




!-----------------------------------------
        goto 200

 190    write(*,*) 'error when opening file'

 200    continue





        deallocate(Dm)
        deallocate(mass)
        deallocate(StaticAt)
        deallocate(Atoms)
        deallocate(MAtoms)
        deallocate(PAtoms)
        deallocate(AName)
        deallocate(AType)
        deallocate(Natom)


!..........................................

 10     format(a2, 3f12.8)
 11     format(i4, '  ', a2, 3f12.8, i3)
 12     format(a16, a2, 3f12.8)


        close(45)
        close(30)
        close(13)
        close(12)
        close(11)



        stop

        End program vwhole
