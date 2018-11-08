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

module checkStruc_module
USE nrtype
USE globalData,only:integerMissing
implicit none
private
public::checkStruc
contains


 ! ************************************************************************************************
 ! public subroutine checkStruc: check data structures 
 ! ************************************************************************************************
 subroutine checkStruc(err,message)
 ! ascii utilities
 USE ascii_util_module,only:split_line
 ! summary of data structures
 USE globalData,only:structInfo
 ! metadata structures
 USE globalData,only:time_meta,forc_meta,attr_meta,type_meta        ! metadata structures
 USE globalData,only:prog_meta,diag_meta,flux_meta,deriv_meta       ! metadata structures
 USE globalData,only:mpar_meta,indx_meta                            ! metadata structures
 USE globalData,only:bpar_meta,bvar_meta                            ! metadata structures
 ! named variables defining strructure elements
 USE var_lookup,only:iLookTIME,iLookFORCE,iLookATTR,iLookTYPE       ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookPROG,iLookDIAG,iLookFLUX,iLookDERIV       ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookPARAM,iLookINDEX                          ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookBPAR,iLookBVAR                            ! named variables showing the elements of each data structure
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iStruct     ! index of data structure
 integer(i4b),parameter               :: nStruct=size(structInfo)  ! number of data structures
 character(len=8192)                  :: longString  ! string containing the indices defined in the structure constructor
 character(len=32),allocatable        :: words(:)    ! vector of words extracted from the long string
 integer(i4b)                         :: ix          ! index of the variable in the data structure
 integer(i4b)                         :: ixTest      ! test the structure constructor = (1,2,3,...,nVar)
 character(len=256)                   :: cmessage    ! error message of downwind routine
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="checkStruc/"

 ! -----
 ! * check that the structure constructors are correct...
 ! ------------------------------------------------------

 ! loop through data structures 
 do iStruct=1,nStruct
  ! convert the lookup structures to a character string
  ! expect the lookup structures to be a vector (1,2,3,...,n)
  select case(trim(structInfo(iStruct)%structName))
   case('time');  write(longString,*) iLookTIME
   case('forc');  write(longString,*) iLookFORCE
   case('attr');  write(longString,*) iLookATTR
   case('type');  write(longString,*) iLookTYPE
   case('mpar');  write(longString,*) iLookPARAM
   case('bpar');  write(longString,*) iLookBPAR
   case('bvar');  write(longString,*) iLookBVAR
   case('indx');  write(longString,*) iLookINDEX
   case('prog');  write(longString,*) iLookPROG
   case('diag');  write(longString,*) iLookDIAG
   case('flux');  write(longString,*) iLookFLUX
   case('deriv'); write(longString,*) iLookDERIV
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  ! check that the length of the lookup structure matches the number of variables in the data structure
  call split_line(longString,words,err,cmessage) ! convert the long character string to a vector of "words"
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  if(size(words)/=structInfo(iStruct)%nVar)then; err=20; message=trim(message)//'unexpected number of elements'; return; end if
  ! check that the elements in the lookup structure are sequential integers (1,2,3,...,n)
  do ix=1,structInfo(iStruct)%nVar
   read(words(ix),*) ixTest  ! convert character to integer; store in ixTest
   if(ixTest/=ix)then ! expect that the ix-th word is equal to ix
    write(message,'(a,i0,a)')trim(message)//'problem with structure constructor iLook'//trim(structInfo(iStruct)%lookName)//' [element=',ix,']'
    err=20; return
   end if
  end do
 end do  ! looping through data structures

 ! -----
 ! * check that the metadata is fully populated...
 ! -----------------------------------------------

 ! loop through data structures
 do iStruct=1,nStruct
  ! check that the metadata is fully populated 
  select case(trim(structInfo(iStruct)%structName))
   case('time');  call checkPopulated(iStruct,time_meta,err,cmessage)
   case('forc');  call checkPopulated(iStruct,forc_meta,err,cmessage) 
   case('attr');  call checkPopulated(iStruct,attr_meta,err,cmessage) 
   case('type');  call checkPopulated(iStruct,type_meta,err,cmessage) 
   case('mpar');  call checkPopulated(iStruct,mpar_meta,err,cmessage) 
   case('bpar');  call checkPopulated(iStruct,bpar_meta,err,cmessage) 
   case('bvar');  call checkPopulated(iStruct,bvar_meta,err,cmessage) 
   case('indx');  call checkPopulated(iStruct,indx_meta,err,cmessage) 
   case('prog');  call checkPopulated(iStruct,prog_meta,err,cmessage) 
   case('diag');  call checkPopulated(iStruct,diag_meta,err,cmessage) 
   case('flux');  call checkPopulated(iStruct,flux_meta,err,cmessage) 
   case('deriv'); call checkPopulated(iStruct,deriv_meta,err,cmessage) 
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors) 
 end do  ! looping through data structures


 contains

  ! ************************************************************************************************
  ! internal subroutine checkPopulated: check that the metadata is fully populated...
  ! ************************************************************************************************
  subroutine checkPopulated(iStruct,metadata,err,message)
  ! access the data type for the metadata structures
  USE data_types,only:var_info 
  ! get index from character string
  USE get_ixname_module,only: get_ixUnknown! variable lookup structure
  implicit none
  ! dummy variables
  integer(i4b),intent(in)   :: iStruct     ! index of data structure
  type(var_info)            :: metadata(:) ! metadata structure 
  integer(i4b),intent(out)  :: err         ! error code
  character(*),intent(out)  :: message     ! error message
  ! local variables
  integer(i4b)              :: iVar        ! index of variable within a data structure
  integer(i4b)              :: jVar        ! index of variable within a data structure (returned from the variable name)
  character(LEN=100)        :: typeName    ! name of variable type to be returned by get_ixUnknown
  character(len=256)        :: cmessage    ! error message of downwind routine
  ! initialize error control
  err=0; message='checkPopulated/'
 
  ! loop through variables
  do iVar=1,size(metadata)

   ! check that this variable is populated 
   if (trim(metadata(iVar)%varname)=='empty') then
    write(message,'(a,i0,a)') trim(message)//trim(structInfo(iStruct)%structName)//'_meta structure is not populated for named variable # ',iVar, ' in structure iLook'//trim(structInfo(iStruct)%lookName)
    err=20; return
   end if

   ! look for the populated variable
   call get_ixUnknown(trim(metadata(iVar)%varname),typeName,jVar,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors) 

   ! check that the variable was found at all
   if (jVar==integerMissing) then
    message = trim(message)//'cannot find variable '//trim(metadata(iVar)%varname)//' in structure '//trim(structInfo(iStruct)%structName)//'_meta; you need to add variable to get_ix'//trim(structInfo(iStruct)%structName)
    err=20; return
   end if
   
   ! check that the variable was found in the correct structure
   if (trim(structInfo(iStruct)%structName)/=typeName) then
    message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' from structure '//trim(structInfo(iStruct)%structName)//'_meta is in structure '//trim(typeName)//'_meta'
    err=20; return
   end if

   ! check that the variable index is correct
   ! This can occur because (1) the code in popMetadat is corrupt (e.g., mis-match in look-up variable); or (2) var_lookup is corrupt.
   if (jVar/=iVar) then
    write(message,'(a,i0,a,i0,a)') trim(message)//'variable '//trim(metadata(iVar)%varname)//' has index ', iVar, ' (expect index ', jVar, '); problem possible in popMetadat, get_ix'//trim(structInfo(iStruct)%structName)//', or var_lookup'
    err=20; return
   end if

  end do  ! looping through variables in structure iStruct

  end subroutine checkPopulated

 end subroutine checkStruc

end module checkStruc_module
