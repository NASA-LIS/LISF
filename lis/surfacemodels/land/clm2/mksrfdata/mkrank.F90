!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "preproc.h"
#include "LIS_misc.h"

subroutine mkrank (n, a, miss, iv, num)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! return indices of largest [num] values in array [a] 
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
!
! $Id: mkrank.F90,v 1.5 2004/05/07 22:18:36 jim Exp $ 
!
!-----------------------------------------------------------------------

  use LIS_precisionMod
  implicit none

! ------------------------ input variables ------------------------
  integer , intent(in) :: n        !array length
  real(r8), intent(in) :: a(0:n)   !array to be ranked
  integer , intent(in) :: miss     !missing data value
  integer , intent(in) :: num      !number of largest values requested
! -----------------------------------------------------------------

! ------------------------ output variables -----------------------
  integer iv(num)      !index to [num] largest values in array [a]
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  real(r8) a_max       !maximum value in array
  integer i            !array index
  real(r8) delmax      !tolerance for finding if larger value
  integer m            !do loop index
  integer k            !do loop index
  logical exclude      !true if data value has already been chosen
! -----------------------------------------------------------------

  delmax = 1.e-06

! -----------------------------------------------------------------
! Find index of largest non-zero number
! -----------------------------------------------------------------

  iv(1) = miss
  a_max = -9999.

  do i = 0, n
     if (a(i)>0. .and. (a(i)-a_max)>delmax) then
        a_max = a(i)
        iv(1)  = i
     end if
  end do

! iv(1) = miss indicates no values > 0. this is an error

  if (iv(1) == miss) then
     write (6,*) 'MKRANK error: iv(1) = missing'
     call endrun
  end if

! -----------------------------------------------------------------
! Find indices of the next [num]-1 largest non-zero number.
! iv(m) = miss if there are no more values > 0
! -----------------------------------------------------------------

  do m = 2, num
     iv(m) = miss
     a_max = -9999.
     do i = 0, n

! exclude if data value has already been chosen

        exclude = .false.
        do k = 1, m-1
           if (i == iv(k)) exclude = .true.
        end do

! if not already chosen, see if it is the largest of the remaining values

        if (.not. exclude) then
           if (a(i)>0. .and. (a(i)-a_max)>delmax) then
              a_max = a(i)
              iv(m)  = i
           end if
        end if
     end do
  end do

  return
end subroutine mkrank


