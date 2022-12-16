!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE uspp_param
  !
  ! ... Ultrasoft and Norm-Conserving pseudopotential parameters
  !  
  USE kinds,      ONLY : dbl
  USE parameters, ONLY : lqmax, nbrx, npsx, nqfx, ndmx
  !
  SAVE
  !
  CHARACTER(LEN=2 ) ::  psd(npsx)   ! name of the pseudopotential

  REAL(dbl) :: &
       dion(nbrx,nbrx,npsx),       &! D_{mu,nu} parameters (in the atomic case)
       betar(ndmx,nbrx,npsx),      &! radial beta_{mu} functions
       jjj(nbrx,npsx),             &! total angular momentum of the beta function
       qqq(nbrx,nbrx,npsx),        &! q_{mu,nu} parameters (in the atomic case)
       qfunc(ndmx,nbrx,nbrx,npsx), &! Q_{mu,nu}(|r|) function for |r|> r_L
       qfcoef(nqfx,lqmax,nbrx,nbrx,npsx), &! coefficients for Q for |r|<r_L
       vloc_at(ndmx,npsx),                &! local potential
       rinner(lqmax,npsx)                  ! values of r_L
  INTEGER :: &
       nbeta(npsx),          &! number of beta functions
       nh(npsx),             &! number of beta functions per atomic type
       nhm,                  &! max number of different beta functions per atom
       kkbeta(npsx),         &! point where the beta are zero
       nqf(npsx),            &! number of coefficients for Q
       nqlc(npsx),           &! number of angular momenta in Q
       ifqopt(npsx),         &! level of q optimization
       lll(nbrx,npsx),       &! angular momentum of the beta function
       iver(3,npsx)           ! version of the atomic code
  INTEGER :: &
       lmaxkb,               &! max angular momentum
       lmaxq                  ! max angular momentum + 1 for Q functions
  LOGICAL :: &
       tvanp(npsx),          &! if .TRUE. the atom is of Vanderbilt type
       newpseudo(npsx)        ! if .TRUE. multiple projectors are allowed
END MODULE uspp_param
!
MODULE uspp
  !
  ! Ultrasoft PPs:
  ! - Clebsch-Gordan coefficients "ap", auxiliary variables "lpx", "lpl"
  ! - beta and q functions of the solid
  !
  USE kinds, ONLY: dbl
  USE parameters, ONLY: lmaxx, lqmax
  IMPLICIT NONE
  PRIVATE
  SAVE
  PUBLIC :: nlx, lpx, lpl, ap, aainit, indv, nhtol, nhtolm, nkb, nkbus, &
       vkb, vkb_ik, dvan, deeq, qq, qb, nhtoj, becsum, uspp_deallocate
  PUBLIC :: qq_so, dvan_so, deeq_nc 
  

  INTEGER, PARAMETER :: &
       nlx  = (lmaxx+1)**2, &! maximum number of combined angular momentum
       mx   = 2*lqmax-1      ! maximum magnetic angular momentum of Q
  INTEGER ::             &! for each pair of combined momenta lm(1),lm(2): 
       lpx(nlx,nlx),     &! maximum combined angular momentum LM
       lpl(nlx,nlx,mx)    ! list of combined angular momenta  LM
  REAL(dbl) :: &
       ap(lqmax*lqmax,nlx,nlx)
  ! Clebsch-Gordan coefficients for spherical harmonics
  !
  INTEGER :: nkb,        &! total number of beta functions, with struct.fact.
             nkbus        ! as above, for US-PP only
  !
  INTEGER, ALLOCATABLE ::&
       indv(:,:),        &! indes linking  atomic beta's to beta's in the solid
       nhtol(:,:),       &! correspondence n <-> angular momentum l
       nhtolm(:,:)        ! correspondence n <-> combined lm index for (l,m)
  !
  INTEGER :: vkb_ik       ! indicate for which ik vkb has been calculated
  COMPLEX(dbl), ALLOCATABLE   :: &
       vkb(:,:)                ! all beta functions in reciprocal space
                               ! indeces: plw, betaf
  REAL(dbl), ALLOCATABLE :: &
       dvan(:,:,:),           &! the D functions of the solid
       deeq(:,:,:,:),         &! the integral of V_eff and Q_{nm} 
       becsum(:,:,:),         &! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
       qq(:,:,:),             &! the q functions in the solid
       nhtoj(:,:)              ! correspondence n <-> total angular momentum
  COMPLEX(dbl), ALLOCATABLE :: & 
       qb(:,:,:,:)             ! the b FT of the Q(r) for each kpt (i,j,ia,ib)
  !
  COMPLEX(dbl), ALLOCATABLE :: & ! variables for spin-orbit/noncolinear case:
       qq_so(:,:,:,:),           &! Q_{nm}
       dvan_so(:,:,:,:),         &! D_{nm}
       deeq_nc(:,:,:,:)           ! \int V_{eff}(r) Q_{nm}(r) dr 
  !
  ! spin-orbit coupling: qq and dvan are complex, qq has additional spin index
  ! noncolinear magnetism: deeq is complex (even in absence of spin-orbit)
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  subroutine aainit(lli)
    !-----------------------------------------------------------------------
    !
    ! this routine computes the coefficients of the expansion of the product
    ! of two real spherical harmonics into real spherical harmonics.
    !
    !     Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
    !
    ! On output:
    ! ap     the expansion coefficients
    ! lpx    for each input limi,ljmj is the number of LM in the sum
    ! lpl    for each input limi,ljmj points to the allowed LM
    !
    ! The indices limi,ljmj and LM assume the order for real spherical
    ! harmonics given in routine ylmr2
    !
    implicit none
    !
    ! input: the maximum li considered
    !  
    integer :: lli
    !
    ! local variables
    !
    integer :: llx, l, li, lj
    real(dbl) , allocatable :: r(:,:), rr(:), ylm(:,:), mly(:,:)
    ! an array of random vectors: r(3,llx)
    ! the norm of r: rr(llx)
    ! the real spherical harmonics for array r: ylm(llx,llx)
    ! the inverse of ylm considered as a matrix: mly(llx,llx)
    real(dbl) :: dum
    !
    ! ANDREA: mortal error is given
    if (lli < 0) call errore('aainit','lli not allowed',-lli)

    if (lli*lli > nlx) call errore('aainit','nlx is too small ',lli*lli)

    llx = (2*lli-1)**2
    if (2*lli-1 > lqmax) &
         call errore('aainit','ap leading dimension is too small',llx)

    allocate (r( 3, llx ))    
    allocate (rr( llx ))    
    allocate (ylm( llx, llx ))    
    allocate (mly( llx, llx ))    

    r(:,:)   = 0.d0
    ylm(:,:) = 0.d0
    mly(:,:) = 0.d0
    ap(:,:,:)= 0.d0

    ! - generate an array of random vectors (uniform deviate on unitary sphere)

    call gen_rndm_r(llx,r,rr)

    ! - generate the real spherical harmonics for the array: ylm(ir,lm)

    call ylmr2(llx,llx,r,rr,ylm)

    !-  store the inverse of ylm(ir,lm) in mly(lm,ir)

    call invmat(llx, ylm, mly, dum)

    !-  for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
    do li = 1, lli*lli
       do lj = 1, lli*lli
          lpx(li,lj)=0
          do l = 1, llx
             ap(l,li,lj) = compute_ap(l,li,lj,llx,ylm,mly)
             if (abs(ap(l,li,lj)) > 1.d-3) then
                lpx(li,lj) = lpx(li,lj) + 1
                if (lpx(li,lj) > mx) &
                     call errore('aainit','mx dimension too small', lpx(li,lj))
                lpl(li,lj,lpx(li,lj)) = l
             end if
          end do
       end do
    end do
    
    deallocate(mly)
    deallocate(ylm)
    deallocate(rr)
    deallocate(r)
    
    return
  end subroutine aainit
  !
  !-----------------------------------------------------------------------
  subroutine gen_rndm_r(llx,r,rr)
    !-----------------------------------------------------------------------
    ! - generate an array of random vectors (uniform deviate on unitary sphere)
    !
    USE constants, ONLY: tpi
    implicit none
    !
    ! first the I/O variables
    !
    integer :: llx         ! input: the dimension of r and rr
    
    real(dbl) :: &
         r(3,llx),  &! output: an array of random vectors
         rr(llx)    ! output: the norm of r
    !
    ! here the local variables
    !
    integer :: ir
    real(dbl) :: costheta, sintheta, phi
    real(dbl), external :: rndm
    
    do ir = 1, llx
       costheta = 2.d0 * rndm() - 1.d0
       sintheta = sqrt ( 1.d0 - costheta*costheta)
       phi = tpi * rndm()
       r (1,ir) = sintheta * cos(phi)
       r (2,ir) = sintheta * sin(phi)
       r (3,ir) = costheta
       rr(ir)   = 1.d0
    end do
    
    return
  end subroutine gen_rndm_r
  !
  !-----------------------------------------------------------------------
  function compute_ap(l,li,lj,llx,ylm,mly)
    !-----------------------------------------------------------------------
    !-  given an l and a li,lj pair compute ap(l,li,lj)
    implicit none
    !
    ! first the I/O variables
    !
    integer :: &
         llx,         &! the dimension of ylm and mly
         l,li,lj       ! the arguments of the array ap
    
    real(dbl) :: &
         compute_ap,  &! this function
         ylm(llx,llx),&! the real spherical harmonics for array r
         mly(llx,llx)  ! the inverse of ylm considered as a matrix
    !
    ! here the local variables
    !
    integer :: ir
    
    compute_ap = 0.d0
    do ir = 1,llx
       compute_ap = compute_ap + mly(l,ir)*ylm(ir,li)*ylm(ir,lj)
    end do
    
    return
  end function compute_ap

  SUBROUTINE uspp_deallocate()
    IF( ALLOCATED( nhtol ) )    DEALLOCATE( nhtol )
    IF( ALLOCATED( indv ) )     DEALLOCATE( indv )
    IF( ALLOCATED( nhtolm ) )   DEALLOCATE( nhtolm )
    IF( ALLOCATED( nhtoj ) )    DEALLOCATE( nhtoj )
    IF( ALLOCATED( vkb ) )      DEALLOCATE( vkb )
    IF( ALLOCATED( becsum ) )   DEALLOCATE( becsum )
    IF( ALLOCATED( qq ) )       DEALLOCATE( qq )
    IF( ALLOCATED( qb ) )       DEALLOCATE( qb )
    IF( ALLOCATED( dvan ) )     DEALLOCATE( dvan )
    IF( ALLOCATED( deeq ) )     DEALLOCATE( deeq )
    IF( ALLOCATED( qq_so ) )    DEALLOCATE( qq_so )
    IF( ALLOCATED( dvan_so ) )  DEALLOCATE( dvan_so )
    IF( ALLOCATED( deeq_nc ) )  DEALLOCATE( deeq_nc )
  END SUBROUTINE uspp_deallocate

end module uspp

