!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_numericalMethodsMod

!BOP
!
! !MODULE: LDT_numericalMethodsMod
!
! !DESCRIPTION:
!  The routines in this module contains a few borrowed routines from
!  the Fortran Numerical Recipes and other numerical methods books.
!
! !REVISION HISTORY:
  implicit none

  PRIVATE

  public :: LDT_quicksort_1arr
  public :: LDT_rand_func

!EOP
  contains

!BOP
!
! !ROUTINE: LDT_quicksort_1arr
! \label{LDT_quicksort_1arr}
!
! !DESCRIPTION:
! Quick sort routine from:
!
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
!
! Obtained Quicksort routine from the site::
!  http://users.bigpond.net.au/amiller/qsort.f90
!
! !REVISION HISTORY:
!  08 Aug 2012: KR Arsenault : Specification in LDT.
!
! !INTERFACE:
  RECURSIVE SUBROUTINE LDT_quicksort_1arr( array_len, list )
! EOP

  implicit none
!- Inputs::
   integer, intent (in)                         :: array_len
   real, dimension (array_len), intent(inout)   :: list

!- Start recursive quicksort routine::
   call quicksort_1(1, size(list))

! --------
   contains
! --------
    recursive subroutine quicksort_1(left_end, right_end)

     integer, intent(in):: left_end, right_end
!   Local variables
     integer            :: i, j
     real               :: reference, temp
     integer, parameter :: max_simple_sort_size = 7

    IF (right_end < left_end + max_simple_sort_size) THEN
     ! Use interchange sort for small lists
       call interchange_sort(left_end, right_end)
    ELSE 
     ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1
       DO
       ! Scan list from left end until element >= reference is found
         DO;  i = i + 1
           IF (list(i) >= reference) EXIT
         END DO
       ! Scan list from right end until element <= reference is found
         DO;  j = j - 1
           IF (list(j) <= reference) EXIT
         END DO
         IF (i < j) THEN
         ! Swap two out-of-order elements::
           temp = list(i)
           list(i) = list(j)
           list(j) = temp
         ELSE IF (i == j) THEN
           i = i + 1
           EXIT
         ELSE
           EXIT
         END IF
       END DO
       IF (left_end < j)  call quicksort_1(left_end, j)
       IF (i < right_end) call quicksort_1(i, right_end)
    END IF
  end subroutine quicksort_1

  subroutine interchange_sort(left_end, right_end)
   INTEGER, INTENT(IN) :: left_end, right_end
!  Local variables
   INTEGER             :: i, j
   REAL                :: temp

   DO i = left_end, right_end - 1
      DO j = i+1, right_end
        IF (list(i) > list(j)) THEN
           temp = list(i)
           list(i) = list(j)
           list(j) = temp
        END IF
      END DO
   END DO
  end subroutine interchange_sort

 end subroutine LDT_quicksort_1arr

! =========================================================

 subroutine LDT_rand_func(idum,rand)

!  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
!  any negative value to initialize or reinitialize the sequence.
!  This function is taken from W.H. Press', "Numerical Recipes" p. 199.

    implicit none
    integer :: seed
    integer :: count, count_rate, count_max
    integer :: idum
    real   :: rand

    call system_clock(count,count_rate,count_max)
    seed = -count-count_max
    call random_seed(seed)
    call random_number(rand)
#if 0 
    real, parameter :: mbig=4000000.,mseed=1618033.,mz=0.
    real, parameter :: fac=1./4000000.

!  According to Knuth, any large mbig, and any smaller (but still large)
!  mseed can be substituted for the above values.
    integer         ::  ma(55)
    integer         ::  iff
    integer         ::  idum
    integer         ::  i,j,k,ii,jj
    real            ::  mk
    integer         ::  inext
    integer         ::  inextp
    real            ::  mj
    real            ::  rand

    iff = 0
    ma = 0

    if (idum.lt.0 .or. iff.eq.0) then
       iff=1
       mj=mseed-float(abs(idum))
       mj=mod(mj,mbig)
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
       enddo
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.mz) ma(i)=ma(i)+mbig
          enddo
       enddo
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56) inext=1
    inextp=inextp+1
    inextp=inextp+1
    if(inextp.eq.56) inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.mz) mj=mj+mbig
    ma(inext)=mj
    rand=mj*fac
    return
#endif

 end subroutine LDT_rand_func

end module LDT_numericalMethodsMod
