!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------

! --------------------------------------------------------------------
! MODULE  Sorting:
!    This module can sort a set of numbers.  The method used is 
! usually referred to as "selection" method.
! --------------------------------------------------------------------

module LVT_SortingMod

  IMPLICIT  NONE
  
  PRIVATE   
  
  PUBLIC :: LVT_sort
  PUBLIC :: LVT_quicksort

CONTAINS

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  LVT_sort(x, Size)
      IMPLICIT  NONE
      REAL,    DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location
      real :: temp

      DO i = 1, Size-1             ! except for the last
         Location = FindMinimum(x, i, Size)  ! find min from this to last
!         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
         temp = x(i)
         x(i) = x(location)
         x(location) = temp
      END DO
    END SUBROUTINE  LVT_sort

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      REAL,    DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      REAL                               :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)          ! assume the first is the min
      Location = Start             ! record its position
      DO i = Start+1, End          ! start with next elements
         IF (x(i) < Minimum) THEN  !   if x(i) less than the min?
            Minimum  = x(i)        !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      REAL,    INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

   recursive subroutine LVT_quicksort(a, size, first, last)
     implicit none
     integer    :: size
     real       :: a(size)
     integer    :: first, last
     real       :: x, t
     
     integer i, j
     
     x = a( (first+last) / 2 )
     i = first
     j = last
     do
        do while (a(i) < x)
           i=i+1
        end do
        do while (x < a(j))
           j=j-1
        end do
        if (i >= j) exit
        t = a(i);  a(i) = a(j);  a(j) = t
        i=i+1
        j=j-1
     end do
     if (first < i-1) call LVT_quicksort(a, size, first, i-1)
     if (j+1 < last)  call LVT_quicksort(a, size, j+1, last)
   end subroutine LVT_quicksort

 end module LVT_SortingMod
