! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module ascii_util_module
USE nrtype
implicit none
private
public::file_open
public::split_line
public::get_vlines
contains


 ! *********************************************************************************************************
 ! public subroutine file_open: open file
 ! *********************************************************************************************************
 subroutine file_open(infile,unt,err,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)              :: infile      ! filename
 integer(i4b),intent(out)             :: unt         ! file unit
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 logical(lgt)                         :: xist        ! .TRUE. if the file exists
 logical(lgt)                         :: xopn        ! .TRUE. if the file is already open
 ! initialize errors
 err=0; message="f-file_open/"
 ! check if the file exists
 inquire(file=trim(infile),exist=xist) ! Check for existence of file
 if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
 end if
 ! check if the file is already open
 inquire(file=trim(infile),opened=xopn) ! Check if the file is open
 if(xopn)then
  message=trim(message)//"FileAlreadyOpen['"//trim(infile)//"']"
  err=20; return
 end if
 ! open file
 open(newunit=unt,file=trim(infile),status="old",action="read",iostat=err)
 if(err/=0)then
   message=trim(message)//"OpenError['"//trim(infile)//"']"
   err=20; return
 end if
 end subroutine file_open


 ! *********************************************************************************************************
 ! public subroutine split_line: split a line of characters into an vector of "words"
 ! *********************************************************************************************************
 subroutine split_line(inline,words,err,message)
 ! do not know how many "words", so use linked lists
 implicit none
 ! declare dummy arguments
 character(*),intent(in)              :: inline     ! line of characters
 character(*),intent(out),allocatable :: words(:) ! vector of "words"
 integer(i4b),intent(out)             :: err      ! error code
 character(*),intent(out)             :: message  ! error message
 ! declare local variables
 integer(i4b),parameter  :: cLen=8192
 character(len=cLen)     :: temp                  ! temporary line of characters
 integer(i4b)            :: iword                 ! loop through words
 integer(i4b),parameter  :: maxWords=1000         ! maximum number of words in a line
 integer(i4b)            :: i1                    ! index at the start of a given word
 character(len=256)      :: cword                 ! the current word
 integer(i4b)            :: nWords                ! number of words in the character string
 ! define pointers for linked list
 type node
  character(len=256)     :: chardat
  integer(i4b)           :: ix
  type(node),pointer     :: next=>null()
 end type node
 type(node),pointer      :: list=>null()
 type(node),pointer      :: current=>null()
 type(node),pointer      :: previous=>null()
 ! start procedure here
 err=0; message='split_line/'
 temp=inline  ! initialize string of characters
 i1=1         ! initialize the index at the start of the first word
 ! ***** loop through the character string
 do iword=1,maxWords
  ! extract a given "word"
  temp=adjustl(temp(i1:len_trim(temp))); if(len_trim(temp)==0) exit
  read(temp,*) cword
  i1  =len_trim(cword)+1
  ! add the variable to the linked list
  if(iword==1)then
   allocate(list); list=node(cword,iword,null())
   current=>list
  else
   allocate(current%next); current%next=node(cword,iword,null())
   current=>current%next
  end if
  ! check that the line has fewer words than maxWords
  if (iword==maxWords)then; err=20; message=trim(message)//"exceedMaxWords [line = "//trim(inline)//"]"; return; end if
 end do
 ! ***** allocate space for the list of words
 nWords = current%ix
 allocate(words(nWords),stat=err)
 if(err/=0)then; err=30; message=trim(message)//"problemAllocateWords"; return; end if
 ! ***** save the list in a vector, and deallocate space as we go...
 current=>list
 do while(associated(current))
  words(current%ix) = current%chardat
  previous=>current; current=>current%next
  deallocate(previous)
 end do
 end subroutine split_line


 ! *********************************************************************************************************
 ! public subroutine get_vlines: get valid lines of data from file and store as a vector of charater strings
 ! *********************************************************************************************************
 subroutine get_vlines(unt,vlines,err,message)
 ! do not know how many valid lines, so use linked lists
 implicit none
 ! declare dummy arguments
 integer(i4b),intent(in)              :: unt         ! file unit
 character(*),intent(out),allocatable :: vlines(:)   ! vector of character strings
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 integer(i4b)            :: iline                    ! loop through lines in the file
 integer(i4b),parameter  :: maxLines=1000000         ! maximum number of valid lines in a file
 character(len=2048)     :: temp                     ! character data or a given line
 integer(i4b)            :: icount                   ! counter for the valid lines
 integer(i4b)            :: iend                     ! index to indicate end of the file
 ! define pointers for linked list
 type node
  character(len=2048)    :: chardat
  integer(i4b)           :: ix
  type(node),pointer     :: next=>null()
 end type node
 type(node),pointer      :: list=>null()
 type(node),pointer      :: current=>null()
 type(node),pointer      :: previous=>null()
 ! start procedure here
 err=0; message='get_vlines/'
 ! ***** get the valid lines of data from the file and store in linked lists *****
 icount=0  ! initialize the counter for the valid lines
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)exit    ! read line of data
  if (temp(1:1)=='!' .or. temp == '')cycle            ! skip comment and empty lines
  icount = icount+1
  ! add the variable to the linked list
  if(.not.associated(list))then
   allocate(list,previous,current); list=node(temp,icount,null())
   current=>list
  else
   allocate(current%next)
   current%next=node(temp,icount,null())
   current=>current%next
  end if
  if (iline==maxLines)then; err=20; message=trim(message)//"exceedMaxLines"; return; end if
 end do  ! looping through the lines in the file (exit clause above will kick in)
 ! ***** allocate space for the valid lines *****
 allocate(vlines(icount),stat=err)
 if(err/=0)then; err=30; message=trim(message)//"problemAllocateVlines"; return; end if
 ! ***** save the list in a vector, and deallocate space as we go... *****
 current=>list
 do while(associated(current))
  vlines(current%ix) = current%chardat
  previous=>current; current=>current%next
  deallocate(previous)
 end do
 if(associated(list)) nullify(list)
 end subroutine get_vlines


end module ascii_util_module
