!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! !ROUTINE: enks_ludcmp
!  \label{enks_ludcmp}
! 
! !DESCRIPTION:
!  LU decomposition after Numerical Recipes 
!  
!  Given a matrix a[1..n][1..n], this routine replaces it by the
!  LU decomposition of a rowwise permutation of itself. a and n are
!  inputs. a is output, arranged as in equation (2.3.14) (see book).
!  indx[1..n] is an output vector that records the row permutation
!  effected by the partial pivoting; This routine is used in combination
!  with lubksb to solve linear equations or invert a matrix.
!  
! !REVISION HISTORY:
!  edited 20 Apr 01, reichle
!
!  - eliminated d (for computation of determinant)
!  - eliminated np (can use automatic arrays nowadays...)
!
!  edited 03 Jun 02, reichle
!  
!  - eliminated parameter NMAX, using dynamic allocation
!
! !INTERFACE: 
SUBROUTINE enks_ludcmp(a,n,indx)
!EOP
  implicit none

  INTEGER n,indx(n)
  REAL a(n,n),TINY
  PARAMETER (TINY=1.0e-20)
  INTEGER i,imax,j,k
  REAL aamax,dum,sum1,vv(n)
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     enddo
     if (aamax.eq.0.) print*,'singular matrix in enks_ludcmp'
     vv(i)=1./aamax
  enddo
  do j=1,n
     do i=1,j-1
        sum1=a(i,j)
        do k=1,i-1
           sum1=sum1-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum1
     enddo
     aamax=0.
     do i=j,n
        sum1=a(i,j)
        do k=1,j-1
           sum1=sum1-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum1
        dum=vv(i)*abs(sum1)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     enddo
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        enddo
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(a(j,j).eq.0.)a(j,j)=TINY
     if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        enddo
     endif
  enddo
  return
END SUBROUTINE enks_ludcmp

  
!BOP
! !ROUTINE: enks_lubksb
! \label{enks_lubksb}
! 
! !DESCRIPTION:
!
!  LU backsubstitution after Numerical Recipes 
! 
!  Solves the set of n linear equations A X = B. Here a[1..n][1..n] is
!  input, not as the matrix A but rather as its LU decomposition, 
!  determined by the routine enks\_ludcmp. indx[1..n] is input as the permutation
!  vector returned by enks\_ludcmp. b[1..n] is input as the right-hand side
!  vector B, and returns the solution vector X. a,n, and indx ar not 
!  modified by this routine and can be left in place for successive calls
!  with different right-hand sides b. This routine takes into account
!  the possibility that b will begin with many zero elements, so it is
!  efficient for use in matrix inversion.
!  
! !REVISION HISTORY:
!  edited 20 Apr 01, reichle
!  
! !INTERFACE: 
SUBROUTINE enks_lubksb(a,n,indx,b)
!EOP  
  implicit none

  INTEGER n,indx(n)
  REAL a(n,n),b(n)
  INTEGER i,ii,j,ll
  REAL sum1
  ii=0
  do i=1,n
     ll=indx(i)
     sum1=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum1=sum1-a(i,j)*b(j)
        enddo
     else if (sum1.ne.0.) then
        ii=i
     endif
     b(i)=sum1
  enddo
  do i=n,1,-1
     sum1=b(i)
     do j=i+1,n
        sum1=sum1-a(i,j)*b(j)
     enddo
     b(i)=sum1/a(i,i)
  enddo
  return
END SUBROUTINE enks_lubksb


