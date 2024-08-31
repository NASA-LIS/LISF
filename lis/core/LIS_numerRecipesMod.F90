!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_numerRecipesMod

!BOP
!
! !MODULE: LIS_numerRecipesMod
!
! !DESCRIPTION:
!  The code in this file contains a few borrowed routines from 
!  the Fortran Numerical Recipes book
!   
! !REVISION HISTORY:
  implicit none
  
  PRIVATE

  public :: LIS_gasdev
  public :: LIS_rand_func
  public :: LIS_indexArrayReal
  public :: LIS_reverse
  public :: LIS_quicksort
  
  !EOP

  integer, parameter :: max_simple_sort_size = 7
  
  interface LIS_quicksort
     module procedure LIS_quicksort_1int
     module procedure LIS_quicksort_1real
     module procedure LIS_quicksort_matrix_int
     module procedure LIS_quicksort_matrix_real
     module procedure LIS_quicksort_2VEC_int
     module procedure LIS_quicksort_2VEC_real     
  end interface
    
  contains

!BOP
! 
! !ROUTINE: LIS_gasdev
! \label{LIS_gasdev}
! 
! !DESCRIPTION:
! 
! Returns a normally distributed deviate with zero mean and unit
! variance, using ran2(idum) as the source of uniform deviates. 
! This routine is adapted from the Numerical Recipies for Fortran 
!
! !REVISION HISTORY: 
!  27 Feb 2005: Sujay Kumar : Specification in LIS.
!  25 Jun 2021: Hiroko Kato Beaudoing and Sarith Mahanama added:
!               LIS_quicksort (from LIS_numericalMethodsMod
!               LIS_indexArrayReal, LIS_reverse    
!      
! !INTERFACE:
    FUNCTION LIS_gasdev(idum)
!EOP
      INTEGER idum
      REAL LIS_gasdev
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
         LIS_gasdev=v2*fac
         iset=1
      else
         LIS_gasdev=gset
         iset=0
      endif
      return
    END FUNCTION LIS_gasdev


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
   

  subroutine LIS_rand_func(idum,rand)

!  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
!  any negative value to initialize or reinitialize the sequence.
!  This function is taken from W.H. Press', ``Numerical Recipes'' p. 199.

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
    if(inextp.eq.56) inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.mz) mj=mj+mbig
    ma(inext)=mj
    rand=mj*fac
    return
#endif
  end subroutine LIS_rand_func

  SUBROUTINE LIS_quicksort_2VEC_int(A,B)
    implicit none
    integer, dimension(:), intent (inout) :: A
    integer, dimension(:), intent (inout) :: B
    integer, allocatable, dimension(:,:)  :: C
    integer                  :: i

    i = size (B)
    allocate (C(1:i,1))
    c (:,1) = B
    call LIS_quicksort_matrix_int (A,C)
    B(:) = C(:,1)
    deallocate (c)

  END SUBROUTINE LIS_quicksort_2VEC_int 

  ! -------------------------------------------------

  SUBROUTINE LIS_quicksort_2VEC_real (A,B)

    implicit none
    integer, dimension(:), intent (inout) :: A
    real,    dimension(:), intent (inout) :: B
    real,    allocatable, dimension(:,:)  :: C
    integer                  :: i

    i = size (B)
    allocate (C(1:i,1))
    c (:,1) = B
    call LIS_quicksort_matrix_real (A,C)
    B(:) = C(:,1)
    deallocate (c)

  END SUBROUTINE LIS_quicksort_2VEC_real

  ! ------------------------------------------------------------------
 
  !BOP
  !
  ! !ROUTINE: LIS_quicksort_1arr
  ! \label{LIS_quicksort_1arr}
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
  
  RECURSIVE SUBROUTINE LIS_quicksort_1int(list)
    ! EOP
    
    implicit none
    !- Inputs::
    
    integer, dimension (:), intent(inout)   :: list

    !- Start recursive quicksort routine::
    call quicksort_1(1, size(list))
   
  contains    
    ! --------
    recursive subroutine quicksort_1(left_end, right_end)

     integer, intent(in):: left_end, right_end
!   Local variables
     integer            :: i, j
     real               :: reference, temp

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

 end subroutine LIS_quicksort_1int

!#####################################################################

  RECURSIVE SUBROUTINE LIS_quicksort_matrix_int (list, matrix)
! EOP

  implicit none
!- Inputs::
   integer, dimension (:), intent(inout)   :: list
   integer, dimension (:,:), intent (inout):: matrix
   integer                                 :: NY

   NY = SIZE (matrix,2)
!- Start recursive quicksort routine::
   call quicksort_2D(1, size(list))

! --------
   contains
! --------
    recursive subroutine quicksort_2D(left_end, right_end)

     integer, intent(in):: left_end, right_end
!   Local variables
     integer            :: i, j
     real               :: reference, temp, temp2D(NY)

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
           temp2D(:) = matrix(i,:)

           list(i) = list(j)
           matrix(i,:) = matrix(j,:)

           list(j) = temp
           matrix(j,:) = temp2D(:)

         ELSE IF (i == j) THEN
           i = i + 1
           EXIT
         ELSE
           EXIT
         END IF
       END DO
       IF (left_end < j)  call quicksort_2D(left_end, j)
       IF (i < right_end) call quicksort_2D(i, right_end)
    END IF
  end subroutine quicksort_2D

  subroutine interchange_sort(left_end, right_end)
   INTEGER, INTENT(IN) :: left_end, right_end
!  Local variables
   INTEGER             :: i, j
   REAL                :: temp, temp2D(NY)

   DO i = left_end, right_end - 1
      DO j = i+1, right_end
        IF (list(i) > list(j)) THEN
           temp = list(i)
           temp2D(:) = matrix(i,:)

           list(i) = list(j)
           matrix(i,:) = matrix(j,:)

           list(j) = temp
           matrix(j,:) = temp2D(:)
        END IF
      END DO
   END DO
  end subroutine interchange_sort

end subroutine LIS_quicksort_matrix_int

!#####################################################################

  RECURSIVE SUBROUTINE LIS_quicksort_matrix_real (list, matrix)
! EOP

  implicit none
!- Inputs::
   integer, dimension (:), intent(inout)   :: list
   real,  dimension (:,:), intent (inout)  :: matrix
   integer                                 :: NY

   NY = SIZE (matrix,2)
!- Start recursive quicksort routine::
   call quicksort_2D(1, size(list))

! --------
   contains
! --------
    recursive subroutine quicksort_2D(left_end, right_end)

     integer, intent(in):: left_end, right_end
!   Local variables
     integer            :: i, j
     real               :: reference, temp, temp2D(NY)

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
           temp2D(:) = matrix(i,:)

           list(i) = list(j)
           matrix(i,:) = matrix(j,:)

           list(j) = temp
           matrix(j,:) = temp2D(:)

         ELSE IF (i == j) THEN
           i = i + 1
           EXIT
         ELSE
           EXIT
         END IF
       END DO
       IF (left_end < j)  call quicksort_2D(left_end, j)
       IF (i < right_end) call quicksort_2D(i, right_end)
    END IF
  end subroutine quicksort_2D

  subroutine interchange_sort(left_end, right_end)
   INTEGER, INTENT(IN) :: left_end, right_end
!  Local variables
   INTEGER             :: i, j
   REAL                :: temp, temp2D(NY)

   DO i = left_end, right_end - 1
      DO j = i+1, right_end
        IF (list(i) > list(j)) THEN
           temp = list(i)
           temp2D(:) = matrix(i,:)

           list(i) = list(j)
           matrix(i,:) = matrix(j,:)

           list(j) = temp
           matrix(j,:) = temp2D(:)
        END IF
      END DO
   END DO
  end subroutine interchange_sort

end subroutine LIS_quicksort_matrix_real

!#####################################################################

  RECURSIVE SUBROUTINE LIS_quicksort_1real(list)
! EOP

  implicit none
!- Inputs::
   real, dimension (:), intent(inout)   :: list

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

 end subroutine LIS_quicksort_1real
  
 ! =========================================================
! LIS_indexArrayReal sorts Array in in ascending order and returns Index containing coreresponding array indices adapted from NR indexx_

 subroutine LIS_indexArrayReal(n,Array,Index)
   implicit none
   integer, intent(in)  :: n
   real   , intent(in)  :: Array(n)
   integer, intent(out) :: Index(n)
   integer, parameter   :: nn=15, nstack=50
   integer              :: k,i,j,indext,jstack,l,r
   integer              :: istack(nstack)
   real                 :: a
   do j = 1,n
      Index(j) = j
   end do
   jstack=0
   l=1
   r=n
   do
      if (r-l < nn) then
         do j=l+1,r
            indext=Index(j)
            a=Array(indext)
            do i=j-1,l,-1
               if (Array(Index(i)) <= a) exit
               Index(i+1)=Index(i)
            end do
            Index(i+1)=indext
         end do
         if (jstack == 0) return
         r=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+r)/2
         call swap(Index(k),Index(l+1))
         call exchangeIndex(Index(l),Index(r))
         call exchangeIndex(Index(l+1),Index(r))
         call exchangeIndex(Index(l),Index(l+1))
         i=l+1
         j=r
         indext=Index(l+1)
         a=Array(indext)
         do
            do
               i=i+1
               if (Array(Index(i)) >= a) exit
            end do
            do
               j=j-1
               if (Array(Index(j)) <= a) exit
            end do
            if (j < i) exit
            call swap(Index(i),Index(j))
         end do
         Index(l+1)=Index(j)
         Index(j)=indext
         jstack=jstack+2
         if (jstack > nstack) then
            write(*,*) 'NSTACK too small in indexArrayReal()'   ! xxx
            error stop
         end if
         if (r-i+1 >= j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         end if
      end if
   end do
 contains
   subroutine exchangeIndex(i,j)
     integer, intent(inout) :: i,j
     integer                :: swp
     if (Array(j) < Array(i)) then
        swp=i
        i=j
        j=swp
     end if
   end subroutine exchangeIndex
   pure elemental subroutine swap(a,b)
     implicit none
     integer, intent(inout) :: a,b
     integer :: dum
     dum=a
     a=b
     b=dum
   end subroutine swap
 end subroutine LIS_indexArrayReal

 ! reverse

 subroutine LIS_reverse (a)

   implicit none

   integer, dimension (:), intent (inout) :: a
   INTEGER            :: Head, Tail, i, n, temp           

   n = size (a)

   Head = 1                             ! start with the beginning
   Tail = n                             ! start with the end
   DO                                   ! for each pair...
      IF (Head >= Tail)  EXIT           !    if Head crosses Tail, exit
      Temp    = a(Head)                 !    otherwise, swap them
      a(Head) = a(Tail)
      a(Tail) = Temp
      Head    = Head + 1                !    move forward
      Tail    = Tail - 1                !    move backward
   END DO                               ! loop bac
 END subroutine LIS_reverse
 
  end module LIS_numerRecipesMod
