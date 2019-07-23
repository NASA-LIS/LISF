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

module childStruc_module
USE nrtype
USE globalData,only:integerMissing ! missing value
USE nr_utility_module,only:arth    ! get a sequence of numbers

implicit none
private
public::childStruc
contains


 ! ************************************************************************************************
 ! public subroutine childStruc: create a child data structure
 ! ************************************************************************************************
 subroutine childStruc(metaParent,mask,                        & ! input
                       metaChild,parent2child_map,err,message)   ! output
 USE data_types,only:var_info               ! data type for the metadata structure
 USE data_types,only:extended_info          ! data type for the extended metadata structure
 implicit none
 ! input variables
 type(var_info),intent(in)                   :: metaParent(:)       ! parent metadata structure
 logical(lgt),intent(in)                     :: mask(:)             ! variables desired
 ! output variables
 type(extended_info),allocatable,intent(out) :: metaChild(:)        ! child metadata structure
 integer(i4b),allocatable,intent(out)        :: parent2child_map(:) ! index of the child variable
 integer(i4b),intent(out)                    :: err                 ! error code
 character(*),intent(out)                    :: message             ! error message
 ! local variables
 integer(i4b)                                :: nParent             ! number of elements in the parent data structure
 integer(i4b)                                :: nChild              ! number of elements in the child data structure
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="childStruc/"

 ! check the size of the input structures
 nParent = size(metaParent)
 if(size(mask)/=nParent)then
  message=trim(message)//'size of mask vector does not match the size of the parent structure'
  err=20; return
 end if

 ! allocate space for the child metadata structure
 nChild = count(mask)
 if(allocated(metaChild)) deallocate(metaChild)
 allocate(metaChild(nChild),stat=err)
 if(err/=0)then
  message=trim(message)//'problem allocating space for the child metadata structure'
  err=20; return
 end if

 ! define mapping with the parent data structure
 metaChild(:)%ixParent = pack(arth(1,1,nParent), mask)

 ! copy across the metadata from the parent structure
 metaChild(:)%var_info = metaParent(metaChild(:)%ixParent)

 ! allows to map from the parent to the child - must carry this around outside
 if(allocated(parent2child_map)) then; err=20; message=trim(message)//'child map already allocated'; return; end if; 
 allocate(parent2child_map(nParent))
 parent2child_map(:) = integerMissing
 if(nChild>0) parent2child_map(metaChild(:)%ixParent) = arth(1,1,nChild)

 end subroutine childStruc

end module childStruc_module
