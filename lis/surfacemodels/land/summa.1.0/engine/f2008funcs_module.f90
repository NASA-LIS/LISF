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

module f2008funcs_module
USE nrtype
implicit none
private
public::cloneStruc
public::findIndex

! define generic interface
interface cloneStruc
 module procedure cloneStruc_rv, cloneStruc_iv
end interface cloneStruc

contains

 ! ************************************************************************************************
 ! public function findIndex: find the first index within a vector
 ! ************************************************************************************************
 function findIndex(vector,desiredValue,missingValue)
 ! finds the first index within a vector
 !  -- if the index does not exist, returns zero
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc
 implicit none
 ! dummy variables
 integer(i4b),intent(in)            :: vector(:)    ! vector to search
 integer(i4b),intent(in)            :: desiredValue ! desired value in the vector
 integer(i4b),intent(in),optional   :: missingValue ! desired missing value if desiredValue is not found
 integer(i4b)                       :: findIndex    ! first index of the desired value in the vector
 ! local variables
 integer(i4b),dimension(1)          :: vecIndex     ! first index of the desired value in the vector (vec of length=1) 

 ! check if the value exisits
 if(any(vector==desiredValue))then

  ! get the index: merge provides a vector with 1s where mask is true and 0s otherwise, so maxloc(merge) is the first index of value=1
  ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc
  vecIndex=maxloc( merge(1, 0, vector==desiredValue) )

 ! value does not exist
 else
  if(present(missingValue))then
   vecIndex=missingValue
  else
   vecIndex=0
  endif
 endif

 ! return function value (extract into a scalar)
 findIndex=vecIndex(1)

 end function findIndex

 ! ************************************************************************************************
 ! public subroutine cloneStruc_rv: clone a data structure (real vector)
 ! ************************************************************************************************
 subroutine cloneStruc_rv(dataVec,lowerBound,source,mold,err,message)
 implicit none
 ! input-output: data vector for allocation/population
 real(dp),intent(inout),allocatable     :: dataVec(:)            ! data vector
 ! input
 integer(i4b),intent(in)                :: lowerBound            ! lower bound
 real(dp),intent(in),optional           :: source(lowerBound:)   ! dataVec = shape of source + elements of source
 real(dp),intent(in),optional           :: mold(lowerBound:)     ! dataVec = shape of mold
 ! error control
 integer(i4b),intent(out)               :: err                   ! error code
 character(*),intent(out)               :: message               ! error message
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b),dimension(1)              :: upperBound            ! upper bound of the data vector (array)
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="cloneStruc_rv/"

 ! check that source and mold are present
 if(.not.present(source) .and. .not.present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (neither given)'
  err=20; return
 end if

 ! check that source and mold are not both present
 if(present(source) .and. present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (both given)'
  err=20; return
 end if

 ! get the upper bounds of the source or the mold vector
 if(present(source))then; upperBound=ubound(source); end if
 if(present(mold))  then; upperBound=ubound(mold);   end if

 ! reallocate spcae
 if(allocated(dataVec)) deallocate(dataVec)
 allocate(dataVec(lowerBound:upperBound(1)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the data vector'; return; end if

 ! copy data
 if(present(source)) dataVec(lowerBound:upperBound(1)) = source(lowerBound:upperBound(1))

 end subroutine cloneStruc_rv

 ! ************************************************************************************************
 ! public subroutine cloneStruc_iv: clone a data structure (integer vector)
 ! ************************************************************************************************
 subroutine cloneStruc_iv(dataVec,lowerBound,source,mold,err,message)
 implicit none
 ! input-output: data vector for allocation/population
 integer(i4b),intent(inout),allocatable :: dataVec(:)            ! data vector
 ! input
 integer(i4b),intent(in)                :: lowerBound            ! lower bound
 integer(i4b),intent(in),optional       :: source(lowerBound:)   ! dataVec = shape of source + elements of source
 integer(i4b),intent(in),optional       :: mold(lowerBound:)     ! dataVec = shape of mold
 ! error control
 integer(i4b),intent(out)               :: err                   ! error code
 character(*),intent(out)               :: message               ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b),dimension(1)              :: upperBound            ! upper bound of the data vector (array)
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="cloneStruc_iv/"

 ! check that source and mold are present
 if(.not.present(source) .and. .not.present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (neither given)'
  err=20; return
 end if

 ! check that source and mold are not both present
 if(present(source) .and. present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (both given)'
  err=20; return
 end if

 ! get the upper bounds of the source or the mold vector
 if(present(source))then; upperBound=ubound(source); end if
 if(present(mold))  then; upperBound=ubound(mold);   end if

 ! reallocate spcae
 if(allocated(dataVec)) deallocate(dataVec)
 allocate(dataVec(lowerBound:upperBound(1)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the data vector'; return; end if

 ! copy data
 if(present(source)) dataVec(lowerBound:upperBound(1)) = source(lowerBound:upperBound(1))

 end subroutine cloneStruc_iv

end module f2008funcs_module
