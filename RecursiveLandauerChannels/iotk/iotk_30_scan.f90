# 1 "iotk_scan.spp"
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 28 "iotk_scan.spp"
#include "iotk_auxmacros.h"
# 30 "iotk_scan.spp"

! SISTEMARE RAW

# 35 "iotk_scan.spp"

# 37 "iotk_scan.spp"
recursive subroutine iotk_scan_begin_x(unit,name,attr,dummy,found,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_scan_interf
  use iotk_files_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                        intent(in)  :: unit
  character(len=*),               intent(in)  :: name
#ifdef __WORKAROUND6
  character(len=*),     optional              :: attr
#else
  character(len=*),     optional, intent(out) :: attr
#endif
  type(iotk_dummytype), optional              :: dummy
  logical,              optional, intent(out) :: found
  integer,              optional, intent(out) :: ierr
  character(iotk_namlenx) :: namel
  character(iotk_attlenx) :: attrl
  character(iotk_vallenx) :: link
  logical :: link_binary,link_raw
  integer :: link_unit
  logical :: binary
  integer :: ierrl,iostat
  logical :: link_found,foundl
  type(iotk_unit), pointer :: this_unit
  integer :: lunit
  character(iotk_fillenx) :: oldfile
  ierrl = 0
  if(present(attr)) attr(1:1)=iotk_eos
  foundl = .false.
  call iotk_strcpy(namel,iotk_strtrim(name),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 73 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this_unit)
  if(associated(this_unit)) then
    if(this_unit%raw) then
      foundl = .true.
      goto 1
    end if
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 87 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  call iotk_scan(lunit, 1,1,namel,attrl,binary,foundl,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 92 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  if(.not.foundl)  then
    call iotk_scan(lunit,-1,1,namel,attrl,binary,foundl,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 98 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 98 "iotk_scan.spp"
call iotk_error_msg(ierrl,'')
# 98 "iotk_scan.spp"
call iotk_error_write(ierrl,"namel",namel)
      goto 1
    end if
    if(.not.foundl) goto 1
  end if
  call iotk_scan_attr(attrl,"iotk_link",link,found=link_found,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 105 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  link_binary=.false.
  if(link_found) then
    call iotk_scan_attr(attrl,"iotk_raw",link_raw,default=.false.,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 112 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
      goto 1
    end if
    if(link_raw) then
      call iotk_scan_attr(attrl,"iotk_binary",link_binary,default=.false.,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 118 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
        goto 1
      end if
    end if
    call iotk_free_unit(link_unit,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 124 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
      goto 1
    end if
    inquire(unit=lunit,name=oldfile,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 129 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
      goto 1
    end if
    call iotk_open_read(link_unit,file=iotk_complete_filepath(link,oldfile),attr=attrl, &
      binary=link_binary,raw=link_raw,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 135 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
      goto 1
    end if
    call iotk_unit_parent(parent=lunit,son=link_unit,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 140 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
      goto 1
    end if
  end if
  if(present(attr)) call iotk_strcpy(attr,attrl,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 146 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_begin",__FILE__,__LINE__)
# 153 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 153 "iotk_scan.spp"
call iotk_error_msg(ierrl,'Tag not found')
# 153 "iotk_scan.spp"
call iotk_error_write(ierrl,"namel",namel)
    ierrl = - ierrl
  end if
  if(ierrl==0 .and. foundl .and. associated(this_unit)) then
    this_unit%level = this_unit%level + 1
!write(0,*) "LEVEL=",this_unit%level,"incrementato"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_begin_x

# 168 "iotk_scan.spp"
recursive subroutine iotk_scan_end_x(unit,name,dummy,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_files_interf
  use iotk_scan_interf
  use iotk_misc_interf
  use iotk_str_interf
  use iotk_unit_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
  character(iotk_namlenx) :: namel
  logical :: binary,foundl,raw
  character(iotk_attlenx) :: attrl
  integer :: ierrl
  integer :: lunit
  type(iotk_unit), pointer :: this_unit
  ierrl = 0
  raw = .false.
  call iotk_strcpy(namel,iotk_strtrim(name),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_end",__FILE__,__LINE__)
# 191 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this_unit)
  if(associated(this_unit)) then
    if(associated(this_unit%parent) .and. this_unit%level == 0) then
      this_unit => this_unit%parent
      call iotk_close_read(lunit,ierr=ierrl)
      if(ierrl/=0) goto 1
      lunit = iotk_phys_unit(unit)
      call iotk_unit_get(lunit,pointer=this_unit)
    end if
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) goto 1
  call iotk_scan(lunit,1,2,namel,attrl,binary,foundl,ierrl)
  if(ierrl/=0 .or. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_end",__FILE__,__LINE__)
# 209 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 209 "iotk_scan.spp"
call iotk_error_msg(ierrl,'foundl')
    goto 1
  end if
  if(iotk_strlen(attrl)/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_end",__FILE__,__LINE__)
# 213 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 213 "iotk_scan.spp"
call iotk_error_msg(ierrl,'An end tag should not contain attributes')
# 213 "iotk_scan.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 213 "iotk_scan.spp"
call iotk_error_write(ierrl,"attr",attrl)
    goto 1
  end if
  if(associated(this_unit)) this_unit%level = this_unit%level - 1
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_end_x

# 226 "iotk_scan.spp"
subroutine iotk_scan_pi_x(unit,name,attr,dummy,found,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_misc_interf
  use iotk_unit_interf
  use iotk_str_interf
  implicit none
  integer,                        intent(in)  :: unit
  character(*),                   intent(in)  :: name
#ifdef __WORKAROUND6
  character(len=*),     optional              :: attr
#else
  character(len=*),     optional, intent(out) :: attr
#endif
  type(iotk_dummytype), optional              :: dummy
  logical,              optional, intent(out) :: found
  integer,              optional, intent(out) :: ierr
  character(iotk_namlenx) :: namel
  character(iotk_attlenx) :: attrl
  type(iotk_unit), pointer :: this
  logical :: binary,foundl
  integer :: ierrl,lunit
  ierrl = 0
  if(present(attr)) attr(1:1)=iotk_eos
  call iotk_strcpy(namel,iotk_strtrim(name),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_pi",__FILE__,__LINE__)
# 253 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  if(associated(this)) then
    if(this%raw) then
      foundl=.true.
      goto 1
    end if
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) goto 1
  call iotk_scan(lunit,1,5,namel,attrl,binary,foundl,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_pi",__FILE__,__LINE__)
# 268 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  if(.not.foundl)  then
    call iotk_scan(lunit,-1,5,namel,attrl,binary,foundl,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_pi",__FILE__,__LINE__)
# 274 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 274 "iotk_scan.spp"
call iotk_error_msg(ierrl,'')
# 274 "iotk_scan.spp"
call iotk_error_write(ierrl,"namel",namel)
      goto 1
    end if
    if(.not.foundl) goto 1
  end if
  if(present(attr)) call iotk_strcpy(attr,attrl,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_pi",__FILE__,__LINE__)
# 281 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
  end if
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_pi",__FILE__,__LINE__)
# 287 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 287 "iotk_scan.spp"
call iotk_error_msg(ierrl,'Tag not found')
# 287 "iotk_scan.spp"
call iotk_error_write(ierrl,"namel",namel)
    ierrl = - ierrl
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_pi_x

# 298 "iotk_scan.spp"
subroutine iotk_scan_empty_x(unit,name,attr,dummy,found,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_misc_interf
  use iotk_str_interf
  use iotk_unit_interf
  implicit none
  integer,                        intent(in)  :: unit
  character(len=*),               intent(in)  :: name
#ifdef __WORKAROUND6
  character(len=*),     optional              :: attr
#else
  character(len=*),     optional, intent(out) :: attr
#endif
  type(iotk_dummytype), optional              :: dummy
  logical,              optional, intent(out) :: found
  integer,              optional, intent(out) :: ierr
  character(iotk_namlenx) :: namel
  character(iotk_attlenx) :: attrl
  type(iotk_unit),pointer :: this
  logical :: binary,foundl
  integer :: ierrl,lunit
  ierrl = 0
  if(present(attr)) attr(1:1)=iotk_eos
  call iotk_strcpy(namel,iotk_strtrim(name),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_empty",__FILE__,__LINE__)
# 325 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  if(associated(this)) then
    if(this%raw) then
      foundl=.true.
      goto 1
    end if
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_empty",__FILE__,__LINE__)
# 338 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  call iotk_scan(lunit,1,3,namel,attrl,binary,foundl,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_empty",__FILE__,__LINE__)
# 343 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
    goto 1
  end if
  if(.not.foundl)  then
    call iotk_scan(lunit,-1,3,namel,attrl,binary,foundl,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_empty",__FILE__,__LINE__)
# 349 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 349 "iotk_scan.spp"
call iotk_error_msg(ierrl,'')
# 349 "iotk_scan.spp"
call iotk_error_write(ierrl,"namel",(namel(1:iotk_strlen(namel))))
      goto 1
    end if
    if(.not.foundl) goto 1
  end if
  if(present(attr)) call iotk_strcpy(attr,attrl,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_empty",__FILE__,__LINE__)
# 356 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
  end if
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_empty",__FILE__,__LINE__)
# 362 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 362 "iotk_scan.spp"
call iotk_error_msg(ierrl,'Tag not found')
# 362 "iotk_scan.spp"
call iotk_error_write(ierrl,"namel",namel)
    ierrl = - ierrl
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_empty_x

# 373 "iotk_scan.spp"
subroutine iotk_scan_tag_x(unit,direction,control,tag,binary,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  integer,                 intent(in)  :: unit
  integer,                 intent(in)  :: direction
  integer,                 intent(out) :: control
  character(iotk_taglenx), intent(out) :: tag
  logical,                 intent(in)  :: binary
  integer,                 intent(out) :: ierr

  integer(iotk_header_kind) :: header
  integer :: taglen,pos,pos1,res,length,iostat
  character(2) :: begin,end
  character(iotk_linlenx) :: line
  character(4) :: predelim
  character(3) :: postdelim
  logical :: found
  ierr = 0
  iostat = 0
  tag  = " "
  if(binary) then
    found = .false.
    do
      if(direction<0) then
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 403 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 403 "iotk_scan.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
      read(unit,iostat=iostat) header
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 409 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 409 "iotk_scan.spp"
call iotk_error_msg(ierr,' ')
# 409 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
      control = modulo(header,iotk_ncontrol+1)
      if(control/=0 .and. control/=128) then
        found = .true.
        taglen  = modulo(header/(iotk_ncontrol+1),iotk_taglenx+1)
        read(unit,iostat=iostat) header,tag(1:taglen)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 418 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 418 "iotk_scan.spp"
call iotk_error_msg(ierr,' ')
# 418 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
!!!!!!
! AGGIUNGO QUESTO PER TOGLIERE I DELIMITATORI DI TAG
! IN MODO CHE SI POSSANO OPZIONALMENTE METTERE NEI FILE BINARI
        select case(control)
        case(1)
          predelim="<" ; postdelim=">"
        case(2)
          predelim="</" ; postdelim=">"
        case(3)
          predelim="<" ; postdelim="/>"
        case(4)
          predelim="<!--" ; postdelim="-->"
        case(5)
          predelim="<?" ; postdelim="?>"
        end select
        pos = index(tag(1:taglen),trim(predelim))
        if(pos/=0) pos=pos+len(trim(predelim))
        if(pos==0) pos=1
        pos1= index(tag(1:taglen),trim(postdelim),back=.true.)
        if(pos1/=0) pos1=pos1-1
        if(pos1==0) pos1=taglen
        tag(1:1+pos1-pos) = tag(pos:pos1)
        taglen=1+pos1-pos
!!!!!!
        if(taglen<len(tag)) tag(taglen+1:taglen+1)=iotk_eos
      end if
      if(direction<0) then
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 450 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 450 "iotk_scan.spp"
call iotk_error_msg(ierr,' ')
# 450 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
      if(found) exit
    end do
    if(direction<0) then
      backspace(unit,iostat=iostat)
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 459 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 459 "iotk_scan.spp"
call iotk_error_msg(ierr,' ')
# 459 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  else
! RISISTEMARE IN MODO CHE SI POSSA AVERE NELLA TAG ANCHE < e >
    if(direction>=0) then
      taglen = 0
      do
        call iotk_getline(unit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 470 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
!        pos = scan(line(1:length),"<")
        pos = iotk_strscan(line,"<")
        if(pos/=0) exit
      end do
      do
!        pos1 = scan(line(pos+1:length),">") + pos
        pos1 = iotk_strscan(line(pos+1:),">") + pos
        if(pos1/=pos) exit
        if(taglen+length-pos+1>len(tag)) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 482 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 482 "iotk_scan.spp"
call iotk_error_msg(ierr,'Tag too long')
          return
        end if
        tag(taglen+1:taglen+1) = " "
        tag(taglen+2:taglen+length-pos+1) = line(pos+1:length)
        taglen = taglen+length-pos+1
        pos = 0
        call iotk_getline(unit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 491 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")    
          return
        end if
      end do
      if(taglen+pos1-pos>len(tag)) then
        call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 496 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 496 "iotk_scan.spp"
call iotk_error_msg(ierr,'Tag too long')
        return
      end if
      tag(taglen+1:taglen+1) = " "
      tag(taglen+2:taglen+pos1-pos) = line(pos+1:pos1-1)
      taglen =taglen+pos1-pos
      res = len_trim(line(1:length))-pos1 ! LA LUNGHEZZA E' TRIMMATA. IN QUESTO MODO SI VA A CAPO
                                          ! SE CI SONO SOLO SPAZI
      if(res>0) then
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 507 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 507 "iotk_scan.spp"
call iotk_error_msg(ierr,' ')
# 507 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
        call iotk_getline(unit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 512 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 517 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 517 "iotk_scan.spp"
call iotk_error_msg(ierr,' ')
# 517 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
        res = length-res
        read(unit,"(a)",iostat=iostat,advance='no') line(1:res)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 523 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 523 "iotk_scan.spp"
call iotk_error_msg(ierr,'length')
# 523 "iotk_scan.spp"
call iotk_error_write(ierr,"res",res)
# 523 "iotk_scan.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
!      pos = verify(tag," ")
!      pos1 = len_trim(tag(1:taglen))
!      pos1 = taglen
       pos = 2
       pos1=taglen
    else
      call iotk_getline(unit,line,length,ierr)
      if(ierr/=0) then
        call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 535 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
        return
      end if
      res = length
!write(0,*) ">>>",res
      do
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 543 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
        call iotk_getline(unit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 548 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
!write(0,*) ">>>%",length,res
        pos = length - res
        pos = scan(line(1:pos),">",back=.true.)
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 556 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
        if(pos/=0) exit
        res = 0
      end do
      taglen=len(tag)+1
      do
        pos1 = scan(line(1:pos-1),"<",back=.true.)
        res = taglen
        if(pos1>0) exit
!CHECK
        tag(res-1:res-1) = " "
        tag(res-pos:res-2) = line(1:pos-1)
        taglen=taglen-pos
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 573 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
        call iotk_getline(unit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 578 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 583 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
        pos = length+1
      end do
!CHECK
      tag(res-1:res-1) = " "
      tag(res-pos+pos1:res-2) = line(pos1+1:pos-1)
      tag(1:len(tag)-res+pos-pos1+1)  =tag(res-pos+pos1:len(tag))
!write(0,*) "%%%%"//tag(1:len(tag)-res+pos-pos1+1)//"%%%%"
      read(unit,"(a)",iostat=iostat,advance="no") line(1:pos1-1)
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_tag",__FILE__,__LINE__)
# 595 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
        return
      end if
!      pos1 = len_trim(tag(1:len(tag)-res+pos-pos1+1))
      pos1 = len(tag)-res+pos-pos1+1-1
!      pos = verify(tag," ")
      pos = 1
    end if
    tag(pos1+1:pos1+1) = iotk_eos
!    write(0,*) "**",direction,"%"//(tag(1:iotk_strlen(tag)))//"%",pos,pos1
! UNA VOLTA RISISTEMATO SOPRA, FARE CONTROLLI PIU' STRINGENTI QUI
    if(tag(pos:pos)=="/" .and. tag(pos1:pos1)/="/") then
      control = 2
      tag = tag(pos+1:pos1)//iotk_eos
    else if(tag(pos:pos)/="/" .and. tag(pos1:pos1)=="/") then
      control = 3
      tag = tag(pos:pos1-1)//iotk_eos
    else if(tag(pos:pos)=="?" .and. tag(pos1:pos1)=="?") then
      control = 5
      tag = tag(pos+1:pos1-1)//iotk_eos
    else if(tag(pos:pos+2)=="!--" .and. tag(pos1-1:pos1)=="--") then
      control = 4
      tag = tag(pos+3:pos1-2)//iotk_eos
    else
      control = 1
      tag = tag(pos:pos1)//iotk_eos
    end if
!    write(0,*) "**",control,"%"//(tag(1:iotk_strlen(tag)))//"%"
  end if
end subroutine iotk_scan_tag_x

# 627 "iotk_scan.spp"
subroutine iotk_scan_x(unit,direction,control,name,attr,binary,found,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  integer,                 intent(in)  :: unit
  integer,                 intent(in)  :: direction
  integer,                 intent(in)  :: control
  character(iotk_namlenx), intent(in)  :: name
  character(iotk_attlenx), intent(out) :: attr
  logical,                 intent(in)  :: binary
  logical,                 intent(out) :: found
  integer,                 intent(out) :: ierr

  character(iotk_taglenx) :: tag
  character(iotk_namlenx) :: r_name
  integer :: level,r_control,pos,pos1
  logical :: lall,match

  found=.false.
  ierr = 0
  if(control==2 .and. direction<0) then
    call iotk_error_issue(ierr,"iotk_scan",__FILE__,__LINE__)
# 651 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
    return
  end if
  level = 0
  ierr = 0
  do
    lall=.false.
    if(direction>=0 .and. level==0) lall=.true.
    if(direction<0  .and. level==0 .and. control/=1) lall=.true.
    if(direction<0  .and. level==1 .and. control==1) lall=.true.
    call iotk_scan_tag(unit,direction,r_control,tag,binary,ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_scan",__FILE__,__LINE__)
# 663 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
      return
    end if
    if(r_control==4) cycle
    if(lall .or. r_control==5) then
      call iotk_tag_parse(tag,r_name,attr,ierr)
      if(ierr/=0) then
        call iotk_error_issue(ierr,"iotk_scan",__FILE__,__LINE__)
# 670 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
# 670 "iotk_scan.spp"
call iotk_error_msg(ierr,'direction')
# 670 "iotk_scan.spp"
call iotk_error_write(ierr,"control",control)
        return
      end if
    end if
    match = lall .and. r_control==control .and. iotk_strcomp(r_name,iotk_strtrim(name))
    if(r_control==5) then
      if(r_name=="iotk") then
        call iotk_check_iotk_attr(unit,attr,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan",__FILE__,__LINE__)
# 679 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
          return
        end if
      end if
    end if
    select case(direction)
    case(0:)
      select case(r_control)
      case(1)
        if(level==0 .and. match) exit
        level = level + 1
      case(2)
        if(level==0 .and. match) exit
        if(level==0) then
          call iotk_scan_tag(unit,-1,r_control,tag,binary,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan",__FILE__,__LINE__)
# 695 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
            return
          end if
          return
        end if
        level = level - 1
      case(3)
        if(level==0 .and. match) exit
      case(5)
        if(level==0 .and. match) exit
      end select
    case(:-1)
      select case(r_control)
      case(2)
        level = level + 1
      case(1)
        if(level==1 .and. match) exit
        if(level==0) then
          call iotk_scan_tag(unit,+1,r_control,tag,binary,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan",__FILE__,__LINE__)
# 715 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
            return
          end if
          return
        end if
        level = level - 1
      case(3)
        if(level==0 .and. match) exit
      case(5)
        if(level==0 .and. match) exit
      end select
    end select
  end do
  if(direction<0) then
    call iotk_scan_tag(unit,+1,r_control,tag,binary,ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_scan",__FILE__,__LINE__)
# 731 "iotk_scan.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.15 ")
    end if
  end if
  found=.true.
end subroutine iotk_scan_x

# 738 "iotk_scan.spp"
subroutine iotk_getline_x(unit,line,length,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer,           intent(in)  :: unit
#ifdef __WORKAROUND6
  character(len=*)               :: line
#else
  character(len=*),  intent(out) :: line
#endif
  integer, optional, intent(out) :: length
  integer, optional, intent(out) :: ierr
  integer :: iostat
#if defined __IOTK_WORKAROUND1
  character(len=iotk_linlenx) :: buffer
#else
  character(len=iotk_getline_buffer) :: buffer
#endif
  integer :: pos,buflen,ierrl,pos1
  logical :: eor
  pos = 0
  ierrl=0
#ifdef __IOTK_WORKAROUND1
! Prima soluzione: Lettura advancing
  read(unit,"(a)",iostat=iostat) buffer
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_getline",__FILE__,__LINE__)
# 765 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 765 "iotk_scan.spp"
call iotk_error_msg(ierrl,'')
# 765 "iotk_scan.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 2
  end if
  buflen = len_trim(buffer)
  line(1:buflen) = buffer(1:buflen)
  line(buflen+1:buflen+1) = iotk_eos
  if(present(length)) length = buflen
#else
  do
    eor = .true.
    read(unit,"(a)",iostat=iostat,eor=1,size=buflen,advance="no") buffer
3   continue
    eor = .false.
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_getline",__FILE__,__LINE__)
# 779 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
      call iotk_error_write(ierrl,"iostat",iostat)
      goto 2
    end if
1   continue
    if(buflen==0) exit
    pos1 = min(pos+buflen,len(line))
    line(pos+1:pos1) = buffer(1:pos1-pos)
    pos = pos1
    if(eor .or. pos>=len(line)) exit
  end do
  if(pos<len(line)) line(pos+1:pos+1) = iotk_eos
  if(present(length)) length = pos
  if(pos>=len(line)) then
    call iotk_error_issue(ierrl,"iotk_getline",__FILE__,__LINE__)
# 793 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 793 "iotk_scan.spp"
call iotk_error_msg(ierrl,'Line too long')
    read(unit,*,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_getline",__FILE__,__LINE__)
# 796 "iotk_scan.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.15 ")
# 796 "iotk_scan.spp"
call iotk_error_msg(ierrl,'iostat')
      goto 2
    end if
  end if
#endif
2 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
  
end subroutine iotk_getline_x
