!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
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
  public :: LDT_gasdev

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

! Other subroutines ...

!BOP
! 
! !ROUTINE: LDT_gasdev
! \label{LDT_gasdev}
! 
! !DESCRIPTION:
! 
! Returns a normally distributed deviate with zero mean and unit
! variance, using ran2(idum) as the source of uniform deviates. 
! This routine is adapted from the Numerical Recipies for Fortran 
!
! !REVISION HISTORY: 
!  27 Feb 2005: Sujay Kumar : Specification in LIS.
!
! !INTERFACE:
    FUNCTION LDT_gasdev(idum)
!EOP
      INTEGER idum
      REAL LDT_gasdev
      INTEGER iset
      REAL fac,gset,rsq,v1,v2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1        v1=2.*ran2(idum)-1.
         v2=2.*ran2(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         LDT_gasdev=v2*fac
         iset=1
      else
         LDT_gasdev=gset
         iset=0
      endif
      return
    END FUNCTION LDT_gasdev


!BOP
! 
! !ROUTINE: ran2
! \label{ran2}
!
! !DESCRIPTION: 
! Long period (>2!1e18) random number generator of L'Ecuyer with
! Bays-Durham shuffle and added safeguards. Returns a uniform
! random deviate between 0.0 and 1.0 (exclusive of the endpoint
! values). Call with ``idum'' a negative integer to initialize;
! thereafter, do not alter ``idum'' between successive deviates
! in a sequence. RNMX should approximate the largest floating 
! value that is less than 1. This function is adapted from 
! Numerical Recipies for Fortran. 
! 
! !REVISION HISTORY: 
!  27 Feb 2005: Sujay Kumar : Specification in LIS.
!
! !INTERFACE:
    FUNCTION ran2(idum)
!EOP

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
           NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
         end do
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
    END function ran2
   

end module LDT_numericalMethodsMod
