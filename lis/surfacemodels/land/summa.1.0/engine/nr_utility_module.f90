module nr_utility_module
USE nrtype
! contains functions that should really be part of the fortran standard, but are not
implicit none
INTERFACE arth
! MODULE PROCEDURE arth_r,arth_d, arth_i
!Hack by SVK
 MODULE PROCEDURE arth_d, arth_i
END INTERFACE
! (everything private unless otherwise specifed)
private
! build vectors of regularly spaced numbers
public::arth
contains

 ! *************************************************************************************************
 ! * the arth function, used to build a vector of regularly spaced numbers 
 ! *************************************************************************************************
#if 0 
 FUNCTION arth_r(first,increment,n)
 implicit none
 REAL(SP), INTENT(IN) :: first,increment
 INTEGER(I4B), INTENT(IN) :: n
 REAL(SP), DIMENSION(n) :: arth_r
 INTEGER(I4B) :: k
 arth_r(1)=first
 if(n>1)then
  do k=2,n
   arth_r(k) = arth_r(k-1) + increment
  end do
 end if
 END FUNCTION arth_r
#endif
 ! ------------------------------------------------------------------------------------------------
 FUNCTION arth_d(first,increment,n)
 implicit none
 REAL(DP), INTENT(IN) :: first,increment
 INTEGER(I4B), INTENT(IN) :: n
 REAL(DP), DIMENSION(n) :: arth_d
 INTEGER(I4B) :: k
 arth_d(1)=first
 if(n>1)then
  do k=2,n
   arth_d(k) = arth_d(k-1) + increment
  end do
 end if
 END FUNCTION arth_d
 ! ------------------------------------------------------------------------------------------------
 FUNCTION arth_i(first,increment,n)
 implicit none
 INTEGER(I4B), INTENT(IN) :: first,increment,n
 INTEGER(I4B), DIMENSION(n) :: arth_i
 INTEGER(I4B) :: k
 arth_i(1)=first
 if(n>1)then
  do k=2,n
   arth_i(k) = arth_i(k-1) + increment
  end do
 end if
 END FUNCTION arth_i

end module nr_utility_module
