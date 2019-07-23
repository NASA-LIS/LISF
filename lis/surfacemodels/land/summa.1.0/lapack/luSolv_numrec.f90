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

module luSolv_module
USE nrtype
implicit none
private
public::ludcmp,lubksb
contains


 ! *********************************************************************************************************************
 ! public subroutine ludcmp: Replaces an n-by-n matrix, a, with the LU decomposition of a row-wise permutation of itself
 ! *********************************************************************************************************************
 SUBROUTINE ludcmp(a,indx,d,err,message)
 IMPLICIT NONE
 ! input/output
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
 INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
 REAL(DP), INTENT(OUT) :: d
 integer(i4b),intent(out)           :: err
 character(*),intent(out)           :: message
 ! local variables
 REAL(DP), DIMENSION(size(a,1)) :: vv
 REAL(DP), PARAMETER :: TINY=1.0e-20_dp
 INTEGER(I4B) :: j,n,imax
 ! initialize error control
 err=0; message='ludcmp/'
 n = size(indx)
 if(size(a,1) /= n .or. size(a,2) /= n)then
  message=trim(message)//'mismatch in size of matrices'
  err=20; return
 endif
 d=1.0_dp
 vv=maxval(abs(a),dim=2)
 if (any(vv == 0.0))then
  message=trim(message)//'singular matrix'
  err=20; return
 endif
 vv=1.0_dp/vv
 do j=1,n
  imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
  if (j /= imax) then
   call swap(a(imax,:),a(j,:))
   d=-d
   vv(imax)=vv(j)
  end if
  indx(j)=imax
  if (a(j,j) == 0.0) a(j,j)=TINY
  a(j+1:n,j)=a(j+1:n,j)/a(j,j)
  a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
 end do
 END SUBROUTINE ludcmp


 ! *********************************************************************************************************************
 ! public subroutine lubksb: Solves the set of n linear equations Ax = b (used with ludcmp to do this)
 ! *********************************************************************************************************************
 SUBROUTINE lubksb(a,indx,b,err,message)
 IMPLICIT NONE
 ! input/output
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
 INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
 integer(i4b),intent(out)           :: err
 character(*),intent(out)           :: message
 ! local
 INTEGER(I4B) :: i,n,ii,ll
 REAL(DP) :: summ
 ! initialize error control
 err=0; message='lubksb/'
 n = size(indx)
 if(size(a,1) /= n .or. size(a,2) /= n)then
  message=trim(message)//'mismatch in size of matrices'
  err=20; return
 endif
 ii=0
 do i=1,n
  ll=indx(i)
  summ=b(ll)
  b(ll)=b(i)
  if (ii /= 0) then
   summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
  else if (summ /= 0.0) then
   ii=i
  end if
  b(i)=summ
 end do
 do i=n,1,-1
  b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
 end do
 END SUBROUTINE lubksb


 ! *********************************************************************************************************************
 ! private subroutine swap: swap two reals
 ! *********************************************************************************************************************
 SUBROUTINE swap(a,b)
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
 REAL(DP), DIMENSION(SIZE(a)) :: dum
 dum=a
 a=b
 b=dum
 END SUBROUTINE swap


 ! *********************************************************************************************************************
 ! private subroutine outerprod: outer product of vectors a and b
 ! *********************************************************************************************************************
 FUNCTION outerprod(a,b)
 REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
 REAL(DP), DIMENSION(size(a),size(b)) :: outerprod
 outerprod = spread(a,dim=2,ncopies=size(b)) * &
  spread(b,dim=1,ncopies=size(a))
 END FUNCTION outerprod


 ! *********************************************************************************************************************
 ! private subroutine imaxloc: returns the index of the maxloc on an array
 ! *********************************************************************************************************************
 FUNCTION imaxloc(arr)
 REAL(DP), DIMENSION(:), INTENT(IN) :: arr
 INTEGER(I4B) :: imaxloc
 INTEGER(I4B), DIMENSION(1) :: imax
 imax=maxloc(arr(:))
 imaxloc=imax(1)
 END FUNCTION imaxloc


end module luSolv_module
