
      module constants
       integer, parameter :: dp = kind(1.0d0)
       integer, parameter :: dpc = kind((1.0d0,1.0d0))
      end module constants





      Program select_junct

      use constants

      implicit none


       TYPE :: atom
        REAL, Dimension(3) :: r ! coordinates
        REAL, Dimension(3) :: f
       END TYPE atom


       TYPE :: name
        CHARACTER(15) :: name    ! name of atom
       END TYPE name


       TYPE(name), dimension(:), allocatable :: AName 
       REAL(dp), dimension(:,:,:,:), allocatable :: Dm
       Real, dimension(:), allocatable :: mass


       Integer :: NATypes, Ntotal, Acounter, PhCount, Ndisp
       CHARACTER(25) :: STR_arg, File_n1, File_n2, F_name, M_name
       CHARACTER(25) :: STR_temp
       TYPE(atom), dimension(:), allocatable :: Atoms


       







       Integer :: i, j, k, t, l, d, h, SCount, Nat
       Integer :: N_arg, Nentries
       Integer :: x,y
       Integer :: kk, Info, jj, found




      
        


!
! reading arguments
!

       N_arg = IARGC()          ! reads the number of arguments

       if (N_arg .ne. 2) then
          write(*,*) 'this program needs two arguments'
          write(*,*) N_arg, 'arguments were entered'
          write(*,*) '--'
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
        read(STR_arg,*) Nat
       CALL GETARG(2, STR_arg) ! returns the argument 3
        read(STR_arg,*) Nentries


    

   

       Ndisp = Nentries ! 

       allocate(Atoms(Nat))  ! allocate array for atoms 





        open (UNIT=11,FILE='forces_plus',STATUS='OLD',action="read",
     1         ERR=190)


        open(UNIT=30, file='static.atoms',form='formatted',
     1       status='unknown')






!----------------------------------------------------
! reading the zero displacement configuration ...
! the first entry
          Acounter = 0

!..............
          read(11,*)
          do t=1, Nat, 1

              Acounter = Acounter + 1
              read(11,*) Atoms(Acounter)%r(:),
     1                   Atoms(Acounter)%f(:)

        l=Acounter      
        ! selecting atoms which are moved
        if (l .le. 30) then 
              write(*,'(3g12.5,3g15.7)') 
     1            Atoms(l)%r(:), Atoms(l)%f(:)


        else ! storing atoms which are not moved
              write(30,'(3g12.5,3g15.7)') 
     1            Atoms(l)%r(:), 0.0,  0.0, 0.0
        endif   

          end do !t
          write(*,*) '--'
! ..................................







!------------------------------------------------------------------







! reading all other displacements ...        
        do d=1, Ndisp, 1 ! over all displacements
          read(11,*)


          Acounter = 0
          do t=1, Nat, 1

              Acounter = Acounter + 1
              read(11,*)  Atoms(Acounter)%r(:), 
     1                   Atoms(Acounter)%f(:)

        l=Acounter      
        ! selecting atoms which are moved
        if (l .le. 30) then 
              write(*,'(3g12.5,3g15.7)') 
     1            Atoms(l)%r(:), Atoms(l)%f(:)

        endif   

          end do !t


          write(*,*) '--'


!---------------------------------------------------------


          end do ! d





!-----------------------------------------
        goto 200

 190    write(*,*) 'error when opening file'

 200    continue





        deallocate(Atoms)



!..........................................

 10     format(a2, 3f12.8)
 11     format(i4, '  ', a2, 3f12.8, i3)
 12     format(a16, a2, 3f12.8)



        close(30)
        close(12)
        close(11)



        stop

        End program select_junct
