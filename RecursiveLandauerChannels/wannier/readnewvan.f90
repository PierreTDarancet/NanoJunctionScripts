!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine readnewvan (is, iunps)
  !---------------------------------------------------------------------
  !
  !     This routine reads the quantities which defines the Vanderbilt
  !     pseudopotential from the file produced by the atomic program.
  !     It is compatible only with the ld1 atomic code
  !
  USE kinds, only: dbl
  USE parameters,   ONLY: nchix, lmaxx, nbrx, ndmx, npsx, lqmax
  USE constants,    ONLY: fpi
  USE atom_module,  ONLY: zmesh, mesh, xmin, dx, r, rab, chi, oc, nchi, &
       lchi, rho_at, rho_atc, nlcc
  USE pseud_module, ONLY: zp, lmax, lloc
  USE uspp_param,   ONLY: vloc_at, dion, betar, qqq, qfunc, nqlc, &
       rinner, nbeta, kkbeta, lll, psd, tvanp
  USE funct_module
  !
  implicit none
  !    First the arguments passed to the subroutine

  integer :: is, iunps
  ! The number of the pseudopotential
  ! the unit with the pseudopotential

  integer :: nb, mb, n, ir, pseudotype, ios, nwfs, ndum, l, ikk
  ! counters on beta functions
  ! counter on mesh points
  ! counters on mesh points
  ! the type of pseudopotential
  ! I/O control
  ! the number of pseudo wavefunctions
  ! dummy integer variable
  ! counter on angular momentum
  ! the kkbeta for each beta

  real(kind=dbl) :: x, etotps, rdum
  ! auxiliary variable
  ! total energy of the pseudoatom
  ! dummy real variable

  logical :: rel
  ! if true the atomic calculation is relativistic

  character (len=75) :: titleps
  ! the title of the pseudo
  character (len=2) :: adum

  if (is.lt.0.or.is.gt.npsx) call errore ('readnewvan', 'Wrong is number', 1)
  read (iunps, '(a75)', err = 100, iostat = ios) titleps

  psd (is) = titleps (7:8)
  read (iunps, '(i5)', err = 100, iostat = ios) pseudotype
  if (pseudotype.eq.3) then
     tvanp (is) = .true.
  else
     tvanp (is) = .false.
  endif
  read (iunps, '(2l5)', err = 100, iostat = ios) rel, nlcc (is)
  read (iunps, '(4i5)', err = 100, iostat = ios) iexch, icorr, igcx, &
       igcc

  dft = '?'

  read (iunps, '(2e17.11,i5)') zp (is) , etotps, lmax (is)

  read (iunps, '(4e17.11,i5)', err = 100, iostat = ios) xmin (is) , &
       rdum, zmesh (is) , dx (is) , mesh (is)


  if (mesh (is) > ndmx) call errore ('readnewvan', 'mesh is too big', 1)

  read (iunps, '(2i5)', err = 100, iostat = ios) nwfs, nbeta (is)
  if (nbeta (is) .gt.nbrx) call errore ('readnewvan', 'nbeta is too large', 1)

  if (nwfs.gt.nchix) call errore ('readnewvan', 'nwfs is too large', 1)
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (rdum, nb = 1, nwfs)

  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (rdum, nb = 1, nwfs)
  do nb = 1, nwfs
     read (iunps, '(a2,2i3,f6.2)', err = 100, iostat = ios) &
          adum, ndum, lchi (nb, is) , oc (nb, is)
     lll (nb, is) = lchi (nb, is)
     !
     ! workaround to distinguish bound states from unbound states
     !
     if (oc (nb, is) <= 0.d0) oc (nb, is) = -1.0
  enddo
  kkbeta (is) = 0
  do nb = 1, nbeta (is)
     read (iunps, '(i6)', err = 100, iostat = ios) ikk
     kkbeta (is) = max (kkbeta (is), ikk)
     read (iunps, '(1p4e19.11)', err = 100, iostat = ios) &
          (betar (ir, nb, is) , ir = 1, ikk)
     do ir = ikk + 1, mesh (is)
        betar (ir, nb, is) = 0.d0
     enddo
     do mb = 1, nb
        read (iunps, '(1p4e19.11)', err = 100, iostat = ios) dion (nb, mb, is)
        dion (mb, nb, is) = dion (nb, mb, is)
        if (pseudotype.eq.3) then
           read (iunps, '(1p4e19.11)', err = 100, iostat = ios) qqq (nb, mb, is)
           qqq (mb, nb, is) = qqq (nb, mb, is)
           read (iunps, '(1p4e19.11)', err = 100, iostat = ios) &
                (qfunc (n, nb, mb, is) , n = 1, mesh (is) )
           do n = 1, mesh (is)
              qfunc (n, mb, nb, is) = qfunc (n, nb, mb, is)
           enddo
        else
           qqq (nb, mb, is) = 0.d0
           qqq (mb, nb, is) = 0.d0
           do n = 0, mesh (is)
              qfunc (n, nb, mb, is) = 0.d0
              qfunc (n, mb, nb, is) = 0.d0
           enddo
        endif
     enddo
  enddo
  !
  !   reads the local potential
  !
  lloc (is) = 0
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) rdum, &
       (vloc_at (ir, is) , ir = 1, mesh (is) )
  !
  !     reads the atomic charge
  !
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (rho_at (ir, &
       is) , ir = 1, mesh (is) )

  !
  !  if present reads the core charge
  !
  if (nlcc (is) ) then
     read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (rho_atc ( &
          ir, is) , ir = 1, mesh (is) )
  endif
  !
  !    read the pseudo wavefunctions of the atom
  !
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) ( (chi (ir, &
       nb, is) , ir = 1, mesh (is) ) , nb = 1, nwfs)

100 call errore ('readnewvan', 'Reading pseudo file', abs (ios) )
  !
  !    set several variables for compatibility with the rest of the
  !    code
  !
  nchi (is) = nwfs

  nqlc (is) = 2 * lmax (is) + 1

  if (nqlc (is) .gt.lqmax.or.nqlc (is) .lt.0) call errore (' readnewvan ',&
       'Wrong  nqlc', nqlc (is) )
  do l = 1, nqlc (is)
     rinner (l, is) = 0.d0
  enddo
  !
  !    compute the radial mesh
  !
  do ir = 1, mesh (is)
     x = xmin (is) + real (ir - 1, dbl) * dx (is)
     r (ir, is) = exp (x) / zmesh (is)
     rab (ir, is) = dx (is) * r (ir, is)
  enddo
  !
  !    For compatibility with rho_atc in the non-US case
  !
  if (nlcc (is) ) then
     do ir = 1, mesh (is)
        rho_atc (ir, is) = rho_atc (ir, is) / fpi / r (ir, is) **2
     enddo

  endif

  close (iunps)
  return
end subroutine readnewvan

