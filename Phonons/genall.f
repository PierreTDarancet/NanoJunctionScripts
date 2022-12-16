
      module constants
       integer, parameter :: dp = kind(1.0d0)
       integer, parameter :: dpc = kind((1.0d0,1.0d0))
      end module constants



      Program gencont
      use constants
      implicit none


       TYPE :: atom
        CHARACTER(15) :: name    ! name of object
        REAL(dp), Dimension(3) :: r ! coordinates
        INTEGER ::  type_index = 0 
        INTEGER ::  total_index = 0      ! 0 if inside, -1 if outside
       END TYPE atom


       TYPE :: name
        CHARACTER(15) :: name    ! name of atom
       END TYPE name


       TYPE(name), dimension(:), allocatable :: AName 


       Integer :: NATypes, Ntotal, Acounter, PhCount, WrCount
       CHARACTER(15) :: STR_arg, CONT_name, CONT_num
       Integer, dimension(:), allocatable ::  Natom 
       TYPE(atom), dimension(:), allocatable :: Atoms 
       !CHARACTER(100) :: STR_debug

       Real(dp), Dimension(3) :: Shift
       REAL(dp), Dimension(3) :: temp1
       Integer :: i, j, k, t, l
       Integer :: N_arg, Nmove
       Real(dp) :: temp0
       Real(dp), Dimension(3) :: LatPar

       Real(dp) :: Disp  ! in a.u.
      


!
! reading arguments
!

       N_arg = IARGC()          ! reads the number of arguments

       if (N_arg .ne. 2) then
          write(*,*) 'this program  needs four arguments'
          write(*,*) N_arg, 'arguments were entered'
          write(*,*) 'enter number of species and'
          write(*,*) 'total number of atoms'
          write(*,*) ' '
          write(*,*) ' '
          write(*,*) '------------------------------------'
          write(*,*)
          write(*,*) 'Usage'
          write(*,*) 'gencont.exe 3  20 '
          write(*,*) ''
          write(*,*) 'displacement is in Angstroms'
          write(*,*)
          write(*,*) 'have fun!'
          stop
       end if


       CALL GETARG(1, STR_arg) ! returns the argument 1
!    write(*,*) STR_arg
       read(STR_arg,*) NATypes
!    write(*,*) NATypes, ' types of atoms'
       CALL GETARG(2, STR_arg) ! returns the argument 2
!    write(*,*) STR_arg
       read(STR_arg,*) Ntotal
!    write(*,*) Ntotal, ' atoms'
!       CALL GETARG(3, STR_arg) ! returns the argument 2
!    write(*,*) STR_arg
!       read(STR_arg,*) Disp



       Disp = 0.015

       allocate(AName(NATypes)) ! array of names
       allocate(Natom(NATypes))
       allocate(Atoms(Ntotal))  ! allocate array for atoms 


        open (UNIT=11,FILE='siesta.cor',STATUS='OLD',action="read",
     1         ERR=190)

        PRINT*, "Reading siesta.cor"

        Acounter = 0

        do t=1, NATypes, 1
          !PRINT*, "Type=", t
          read(11,*) AName(t), Natom(t)
           write (*,*) AName(t),  Natom(t)
         ! PRINT*, AName(t),  Natom(t)
          do i=1, Natom(t), 1 
           Acounter = Acounter + 1
             !PRINT*, Acounter
             if (Acounter .gt. Ntotal) then 
               write(*,*) 'ERROR!'
               write(*,*) 'Number of atoms is wrong'
               write(*,*) 'exiting ...'
               goto 200
             endif
             !read(11,*) STR_debug
             !PRINT*, STR_debug
             read(11,*) Atoms(Acounter)%r(:), 
     1           Atoms(Acounter)%type_index,
     2            Atoms(Acounter)%name,
     3           Atoms(Acounter)%total_index
          end do  
        enddo


        write (*,*) Ntotal
        write (*,*) 'Coordinates'

        Nmove = Ntotal 
       

        PhCount = 0 ! counting distorted structures
        WrCount = 0 ! counting written structures      
        
        write(CONT_num,'(i10)') WrCount
        write(*,*) 'CONT_num = ', ADJUSTL(CONT_num)
        CONT_name = 'Plus_' // ADJUSTL(CONT_num)
        open (UNIT=40, FILE=CONT_name, STATUS='unknown',
     1         form='formatted',action="write",ERR=190)



      do i=1, Ntotal, 1

       do j=1, 3, 1
         temp1(j)=Atoms(i)%r(j)
       end do

       write(40,13)  temp1(:),  Atoms(i)%type_index,
     1               Atoms(i)%name,
     2               Atoms(i)%total_index
      end do

        close(40)


        do l=1, Ntotal, 1 ! displace all atoms

         if (l .le. Nmove) then   ! move only Nmove first atoms
  
          do k=1, 3, 1      ! in three directions
             PhCount = PhCount + 1

              temp0 = Atoms(l)%r(k)
             

              Atoms(l)%r(k) = Atoms(l)%r(k) + Disp


              WrCount = WrCount + 1





        write(CONT_num,'(i10)') WrCount
        write(*,*) 'CONT_num = ', ADJUSTL(CONT_num)
        CONT_name = 'Plus_' // ADJUSTL(CONT_num)
        open (UNIT=40, FILE=CONT_name, STATUS='unknown',
     1         form='formatted',action="write",ERR=190)



      do i=1, Ntotal, 1

       do j=1, 3, 1
         temp1(j)=Atoms(i)%r(j)
       end do

       write(40,13)  temp1(:),  Atoms(i)%type_index,
     1               Atoms(i)%name,
     2               Atoms(i)%total_index
      end do

        close(40)



             Atoms(l)%r(k) = Atoms(l)%r(k) - 2.0*Disp


        write(CONT_num,'(i10)') WrCount
        write(*,*) 'CONT_num = ', ADJUSTL(CONT_num)
        CONT_name = 'Minus_' // ADJUSTL(CONT_num)
        open (UNIT=41, FILE=CONT_name, STATUS='unknown',
     1         form='formatted',action="write",ERR=190)



      do i=1, Ntotal, 1

       do j=1, 3, 1
         temp1(j)=Atoms(i)%r(j)
       end do

       write(41,13)  temp1(:),  Atoms(i)%type_index,
     1               Atoms(i)%name,
     2               Atoms(i)%total_index
      end do

        close(41)




    

                   write(*,*) WrCount,' structures are written'  

                   ! restore the coordinates
                   Atoms(l)%r(k) = temp0
       
          end do ! k=1 ->3 

          endif ! l <= Nmove

        end do ! l=1 -> all atoms  







!-----------------------------------------
        goto 200

 190    write(*,*) 'error when opening file'

 200    continue

        deallocate(Atoms)
        deallocate(AName)
        deallocate(Natom)


!..........................................

 10     format(a2, 3f12.6)
 11     format(i4, '  ', a2, 3f12.6, i3)
 12     format(a16, a2, 3f12.6)
 13     format(3f12.6, i4 ,' ', a2, i4)




        close(11)



        stop

        End program gencont
