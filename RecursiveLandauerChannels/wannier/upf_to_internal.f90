!
! Copyright (C) 2004 PWSCF-CP-FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module is USEd, for the time being, as an interface
! between the UPF pseudo type and the pseudo variables internal representation

!=----------------------------------------------------------------------------=!
  MODULE upf_to_internal
!=----------------------------------------------------------------------------=!

  IMPLICIT NONE
  SAVE

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

!
!---------------------------------------------------------------------
subroutine set_pseudo_upf (is, upf)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal PWscf variables
  !
  ! PWSCF modules
  !
  USE atom_module,  ONLY: mesh, msh, r, rab, &
       chi, oc, nchi, lchi, jchi, rho_at, rho_atc, nlcc
  USE pseud_module, ONLY: lloc, lmax, zp
  USE uspp_param, ONLY: vloc_at, dion, betar, qqq, qfcoef, qfunc, nqf, nqlc, &
       rinner, nbeta, kkbeta, lll, jjj, psd, tvanp
  USE funct_module, ONLY: dft, which_dft
  !
  USE ions_module, ONLY: zv
  USE spin_orb_module, ONLY: lspinorb
  USE pseudo_types_module
  !
  implicit none
  !
  real(kind=dbl), parameter :: rcut = 10.0_dbl
  integer :: is, ir
  !
  !     Local variables
  !
  integer :: nb
  TYPE (pseudo_upf) :: upf
  !
  !
  zp(is)  = upf%zp
  psd (is)= upf%psd
  tvanp(is)=upf%tvanp
  nlcc(is) = upf%nlcc
  dft = upf%dft
  call which_dft (upf%dft)
  mesh(is) = upf%mesh
  !
  nchi(is) = upf%nwfc
  lchi(1:upf%nwfc, is) = upf%lchi(1:upf%nwfc)
  oc(1:upf%nwfc, is) = upf%oc(1:upf%nwfc)
  chi(1:upf%mesh, 1:upf%nwfc, is) = upf%chi(1:upf%mesh, 1:upf%nwfc)
  !
  nbeta(is)= upf%nbeta
  kkbeta(is)=0
  do nb=1,upf%nbeta
     kkbeta(is)=max(upf%kkbeta(nb),kkbeta(is))
  end do
  betar(1:upf%mesh, 1:upf%nbeta, is) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  dion(1:upf%nbeta, 1:upf%nbeta, is) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !
  lmax(is) = upf%lmax
  nqlc(is) = upf%nqlc
  nqf (is) = upf%nqf
  lll(1:upf%nbeta,is) = upf%lll(1:upf%nbeta)
  rinner(1:upf%nqlc,is) = upf%rinner(1:upf%nqlc)
  qqq(1:upf%nbeta,1:upf%nbeta,is) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
  qfunc (1:upf%mesh, 1:upf%nbeta, 1:upf%nbeta, is) = &
       upf%qfunc(1:upf%mesh,1:upf%nbeta,1:upf%nbeta)
  qfcoef(1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta, is ) = &
       upf%qfcoef( 1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta )
  !
  r  (1:upf%mesh, is) = upf%r  (1:upf%mesh)
  rab(1:upf%mesh, is) = upf%rab(1:upf%mesh)

  if (lspinorb.and..not.upf%has_so) &
     call errore('upf_to_internal','At least one non s.o. pseudo',-1)
   
  lspinorb=lspinorb.and.upf%has_so
  if (upf%has_so) then
     jchi(1:upf%nwfc, is) = upf%jchi(1:upf%nwfc)
     jjj(1:upf%nbeta, is) = upf%jjj(1:upf%nbeta)
  else
     jchi(1:upf%nwfc, is) = 0.d0
     jjj(1:upf%nbeta, is) = 0.d0
  endif
  !
  if ( upf%nlcc) then
     rho_atc(1:upf%mesh, is) = upf%rho_atc(1:upf%mesh)
  else
     rho_atc(:,is) = 0.d0
  end if
  rho_at (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)
  !!! TEMP
  lloc(is) = 0
  !!!
  vloc_at(1:upf%mesh,is) = upf%vloc(1:upf%mesh)

  do ir = 1, mesh (is)
    if (r (ir, is) .gt.rcut) then
        msh (is) = ir
        goto 5
    endif
  enddo
  msh (is) = mesh (is)
  !
  ! force msh to be odd for simpson integration
  !
5 msh (is) = 2 * ( (msh (is) + 1) / 2) - 1

  zv(is) = zp(is)  !!! maybe not needed: it is done in setup

end subroutine set_pseudo_upf

!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
