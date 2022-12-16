! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
c $Id: constr.f,v 1.6 2003/06/23 09:46:16 ordejon Exp $

      subroutine constr( cell, na, isa, amass, xa, stress, fa, ntcon )
c *****************************************************************
c User-written routine to implement specific geometric constraints,
c by orthogonalizing the forces and stress to undesired changes.
c Arguments:
c real*8  cell(3,3)    : input lattice vectors (Bohr)
c integer na           : input number of atoms
c integer isa(na)      : input species indexes
c real*8  amass(na)    : input atomic masses
c real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
c real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
c real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
c integer ntcon        : total number of positions constr. imposed
c *****************************************************************
      implicit         none
      integer          na, isa(na), ntcon, i
      double precision amass(na), cell(3,3), fa(3,na),
     .                 stress(3,3), xa(3,na), fz1, fz2
      fz1 = 0.
      fz2 = 0.
      do i=1,16
         fz1 = fz1 + fa(3,na-128+2*i-1)
         fz2 = fz2 + fa(3,na-128+2*i)
      enddo
      fz1=fz1/16
      fz2=fz2/16

      do i=1,48
         fa(2,na-96+2*i-1)=0.0
         fa(2,na-96+2*i)=0.0
         fa(2,na-96+2*i-1)=0.0
         fa(2,na-96+2*i)=0.0
         fa(3,na-96+2*i-1)=fz1
         fa(3,na-96+2*i)=fz2
      enddo
      
c Write here your problem-specific code.

      end

