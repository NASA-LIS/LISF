!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
! by Matthew Garcia, GEST Research Associate
!    Hydrological Sciences Branch (Code 614.3)
!    NASA-GSFC 
!
! Calculation of multiquadric-biharmonic interpolation weights
! -- Cartesian distance calculation           -- function cartdist
! -- Establish station and domain matrices    -- subroutine mqbmatrix
! -- matrix inversion for MQ-B multiplication -- subroutine matinvert
!
! Based on original programs for MQ-B calculation by Siryanone (1988) and 
!   Haider (1994)
!
! References
!   Hardy, R.L., 1990:  Theory and applications of the multiquadric-biharmonic 
!      method.  Computers Math. Applic., v. 19, pp. 163-208.
!   Haider, S.K., 1994:  Spatial storm characteristics and basin response.  
!      M.S. Thesis, Dept. of Hydrology and Water Resources, University of 
!      Arizona.  261 pp.
!   Supachai, S., 1988:  Comparative studies of kriging, multiquadric-
!      biharmonic, and other methods for solving mineral resource problems.
!      Ph.D. Dissertation, Department of Earth Sciences, Iowa State University.
!      355 pp.
!
! Version History
!   05 May 2005  Matthew Garcia  Initial specification
!
module mqb_module
!
  implicit none
!
contains
!
!
real(4) function cartdist2(x1,y1,x2,y2)
!
! argument variables
  real(4), intent(IN) :: x1,y1,x2,y2
!
! local variables
  real(4) :: xdist,ydist,dsq
!
  xdist = x1 - x2
  ydist = y1 - y2
  dsq = (xdist * xdist) + (ydist * ydist)
  cartdist2 = sqrt(dsq)
end function cartdist2
!
!
real(4) function mqweight(x1,y1,z1,x2,y2,z2,param)
!
! argument variables
  real(4), intent(IN) :: x1,y1,z1,x2,y2,z2,param
!
! local variables
  real(4) :: xdist,ydist,zdist,wgtsq
!
  xdist = x1 - x2
  ydist = y1 - y2
  zdist = z1 - z2
  wgtsq = (xdist * xdist) + (ydist * ydist) + (zdist * zdist) + param
  mqweight = sqrt(wgtsq)
end function mqweight
!
!
subroutine mqbmatrix(nn,npts,lx,ry,inc,sdata,stnmtx,wgtarr,locarr)
!
! argument variables
  integer, intent(IN) :: nn,npts
  real(4), intent(IN) :: lx,ry,inc
  real(4), intent(INOUT) :: sdata(:,:)
  real(4), intent(INOUT) :: stnmtx(:,:)
  real(4), intent(INOUT) :: wgtarr(:,:,:)
  real(4), intent(INOUT) :: locarr(:,:)
!
! local variables
  integer :: s,i,j
  real :: det,prodsum,prodavg
  real :: xloc,yloc
  real, allocatable :: a(:,:),inv_a(:,:),prod(:,:)
!
  print *,'MSG: mqb_module -- establishing MQ-B station matrix'
  allocate(a(nn+1,nn+1))
  a = 0.0
  do i = 1,nn
    do j = 1,nn
      a(i,j) = cartdist2(sdata(2,i),sdata(3,i),sdata(2,j),sdata(3,j))
    end do
  end do
  do i = 1,nn
    a(i,nn+1) = 1.0
    a(nn+1,i) = 1.0
  end do
  a(nn+1,nn+1) = 0.0
  allocate(inv_a(nn+1,nn+1))  
  call gjinvert(a,inv_a,nn+1)
  print *,'MSG: mqb_module -- inverse MQ-B station matrix obtained'
  stnmtx = inv_a
!
!  allocate(prod(nn+1,nn+1))  
!  prod = matmul(a,inv_a)
!  prodsum = sum(prod)
!  prodavg = prodsum / (nn + 1)
!  print *,'DBG: mqb_module -- product element average =',prodavg
!  deallocate(prod)
!
  print *,'MSG: mqb_module -- establishing MQ-B domain matrix'
  wgtarr = 0.0
  do i = 1,nn
    do j = 1,npts
      wgtarr(1,j,i) = cartdist2(sdata(2,i),sdata(3,i),locarr(j,1),locarr(j,2))
    end do
  end do
!
  deallocate(a)
  deallocate(inv_a)
  return
!
end subroutine mqbmatrix
!
!
! Gauss-Jordan Method for Matrix Inversion
!
! Inversion of n-square matrix
!
! Loosely based on original code by David Goodrich, USDA ARS, for use with 
!   MQ-B program by Haider (1994)
!
! References
!   Haider, S.K., 1994:  Spatial storm characteristics and basin response.  
!      M.S. Thesis, Dept. of Hydrology and Water Resources, University of 
!      Arizona.  261 pp.
!
! Version History
!   05 May 2005  Matthew Garcia  Initial specification
!
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine gjinvert(a_orig,inv_a,n)
!
! argument variables
  integer, intent(IN) :: n
  real, intent(IN) :: a_orig(n,n)
  real, intent(OUT) :: inv_a(n,n)
!
! local variables
  real :: big,dum,pivinv
  real :: a(n,n),b(n,1)
  integer :: irow,icol
  integer :: i,j,k,l,ll
  integer :: m = 0
  integer :: ipiv(n),indxr(n),indxc(n)
!
  a = a_orig
  ipiv = 0
  do i = 1,n
    big = 0
    do j = 1,n
      if (ipiv(j).ne.1) then
        do k = 1,n
          if (ipiv(k).eq.0) then
            if (abs(a(j,k)).ge.big) then
              big = abs(a(j,k))
              irow = j
              icol = k
            endif
          else if (ipiv(k).gt.1) then
            print *,'ERR: mqb_module -- error in matrix inversion'
            stop
          endif
        end do
      endif
    end do
!
    ipiv(icol) = ipiv(icol)+1
    if (irow.ne.icol) then
      do l = 1,n
        dum = a(irow,l)
        a(irow,l) = a(icol,l)
        a(icol,l) = dum
      end do
      do l = 1,m
        dum = b(irow,l)
        b(irow,l) = b(icol,l)
        b(icol,l) = dum
      end do  
    endif
    indxr(i)=irow
    indxc(i)=icol
    if (a(icol,icol).eq.0) then
      print *,'ERR: mqb_module -- singular matrix encountered'
      stop
    endif
!
    pivinv = 1 / a(icol,icol)
    a(icol,icol) = 1
    do l = 1,n
      a(icol,l) = a(icol,l) * pivinv
    end do
    do l = 1,m
      b(icol,l) = b(icol,l) * pivinv
    end do
!
    do ll = 1,n
      if (ll.ne.icol) then
        dum = a(ll,icol)
        a(ll,icol) = 0
        do l = 1,n
          a(ll,l) = a(ll,l) - a(icol,l) * dum
        end do
        do l = 1,m
          b(ll,l) = b(ll,l) - b(icol,l) * dum
        end do
      endif
    end do
  end do
!
  do j = 1,n
    l = n + 1 - j
    if (indxr(l).ne.indxc(l)) then
      do k = 1,n
        dum = a(k,indxr(l))
        a(k,indxr(l)) = a(k,indxc(l))
        a(k,indxc(l)) = dum
      end do
    endif
  end do
!
  inv_a = a
  return
end subroutine gjinvert
!
!
end module mqb_module
