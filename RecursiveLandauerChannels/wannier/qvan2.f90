!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine qvan2 (ngy, ih, jh, np, qmod, qg, ylmk0)
  !-----------------------------------------------------------------------
  !
  !    This routine computes the fourier transform of the Q function assum
  !    that the radial fourier trasform is already computed and stored
  !    in qrad.
  !
  !    The formula implemented here is
  !
  !     q(g,l,k) = sum_lm (-i)^l ap(lm,l,k) yr_lm(g^) qrad(g,l,l,k)
  !
  !
#include "machine.h"
  USE kinds, ONLY: dbl
  USE constants, ONLY : ZERO, CZERO, ONE, TWO, THREE, CI, EPS_m6
  USE us_module, ONLY: dq, qrad
  USE uspp_param, ONLY: lmaxq, nbrx
  USE uspp, ONLY: nlx, lpl, lpx, ap, indv, nhtolm
  IMPLICIT NONE

  integer :: ngy, ih, jh, np
  ! input: the number of G vectors to compute
  ! input: the first index of Q
  ! input: the second index of Q
  ! input: the number of the pseudopotential

  real(kind=dbl) :: ylmk0 (ngy, lmaxq * lmaxq), qmod (ngy)
  ! the spherical harmonics
  ! input: moduli of the q+g vectors
  complex(kind=dbl) :: qg (ngy)
  ! output: the fourier transform of interest
  !
  !     here the local variables
  !

  complex(kind=dbl) :: sig
  ! (-i)^L

  integer :: nb, mb, nmb, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3
  ! the atomic index corresponding to ih
  ! the atomic index corresponding to jh
  ! combined index (nb,mb)
  ! the lm corresponding to ih
  ! the lm corresponding to jh
  ! counter on g vectors
  ! the actual LM
  ! the angular momentum L
  ! the possible LM's compatible with ih,j
  ! counters for interpolation table

  real(kind=dbl) :: sixth, dqi, qm, px, ux, vx, wx, uvx, pwx, work
  ! 1 divided by six
  ! 1 divided dq
  ! qmod/dq
  ! measures for interpolation table
  ! auxiliary variables for intepolation
  ! auxiliary variable
  !
  LOGICAL :: ltest
  !
  !     compute the indices which correspond to ih,jh
  !
  sixth = ONE / ( TWO * THREE )
  dqi = 1 / dq
  nb = indv (ih, np)
  mb = indv (jh, np)
  if (nb.ge.mb) then
     nmb = nb * (nb - 1) / 2 + mb
  else
     nmb = mb * (mb - 1) / 2 + nb
  endif
  ivl = nhtolm(ih, np)
  jvl = nhtolm(jh, np)
  if (nb.gt.nbrx) call errore (' qvan2 ', ' nb.gt.nbrx ', nb)
  if (mb.gt.nbrx) call errore (' qvan2 ', ' mb.gt.nbrx ', mb)
  if (ivl.gt.nlx) call errore (' qvan2 ', ' ivl.gt.nlx  ', ivl)
  if (jvl.gt.nlx) call errore (' qvan2 ', ' jvl.gt.nlx  ', jvl)
  qg(:) = CZERO
  !
  !    and make the sum over the non zero LM
  !
  do lm = 1, lpx (ivl, jvl)
     lp = lpl (ivl, jvl, lm)
     !
     !     extraction of angular momentum l from lp:
     !
     if (lp.eq.1) then
        l = 1
     elseif ( (lp.ge.2) .and. (lp.le.4) ) then
        l = 2
     elseif ( (lp.ge.5) .and. (lp.le.9) ) then
        l = 3
     elseif ( (lp.ge.10) .and. (lp.le.16) ) then
        l = 4
     elseif ( (lp.ge.17) .and. (lp.le.25) ) then
        l = 5
     elseif ( (lp.ge.26) .and. (lp.le.36) ) then
        l = 6
     elseif ( (lp.ge.37) .and. (lp.le.49) ) then
        l = 7
     else
        call errore (' qvan ', ' lp > 49 ', lp)
     endif
     sig = (-CI) ** (l - 1)
     sig = sig * ap (lp, ivl, jvl)
     do ig = 1, ngy
        !
        ! calculate quantites depending on the module of G only when needed
        !
        IF ( ig > 1 ) ltest = ABS( qmod(ig) - qmod(ig-1) ) > EPS_m6
        !
        IF ( ig == 1 .OR. ltest ) THEN
           qm = qmod (ig) * dqi
           px = qm - int (qm)
           ux = ONE - px
           vx = TWO - px
           wx = THREE - px
           i0 = INT (qm) + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           uvx = ux * vx * sixth
           pwx = px * wx * 0.5_dbl
           work = qrad (i0, nmb, l, np) * uvx * wx + &
                  qrad (i1, nmb, l, np) * pwx * vx - &
                  qrad (i2, nmb, l, np) * pwx * ux + &
                  qrad (i3, nmb, l, np) * px * uvx
        endif
        qg (ig) = qg (ig) + sig * ylmk0 (ig, lp) * work
     enddo
  enddo

  return
end subroutine qvan2

