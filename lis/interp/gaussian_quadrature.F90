!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine gaqd(nlat,edges)

  integer, intent(in) :: nlat
  real, intent(out) :: edges(nlat+2)
  real :: lam(nlat), e(nlat), a
  integer j,i

! initialize diagonal and subdiagonal elements
  e = 0.
  lam = 0.

! subdiagonal elements 
  do j = 2,nlat
    e(j) = sqrt(((j-1.)**2.)/((2.*j-1.)*(2.*j-3.)))
  enddo
! eigenvalue decomposition
  call tqli(lam,e,nlat)

! sort in descending order
  do j=2,nlat
    a = lam(j)
    do i=j-1,1,-1
      if (lam(i).le.a) then
        goto 3
      endif
      lam(i+1) = lam(i)
   end do
   i=0
 3 lam(i+1) = a
  end do

! radial points converted to degrees latitude
  do j = 1,nlat/2+1
   edges(j+1) =        -(acos(lam(j))*180.0/3.1415927 - 90.0)
   edges(nlat+2 - j) = (acos(lam(j))*180.0/3.1415927 - 90.0)
  enddo
  edges(1) = -90.
  edges(nlat+2) = 90.
! round to nearest 1e-5 degree
  do j = 1,nlat
    edges(j) = nint(edges(j)*10000.)/10000.
  enddo

! reverse direction of array
!  tmp = edges
!  do j = 1,nlat
!   edges
!  enddo


  return
end subroutine gaqd

subroutine tqli(d,e,n)
  integer n
  real d(n),e(n)
  integer i,iter,k,l,m
  real b,c,dd,f,g,p,r,s,pythag,sgng

  do i=2,n
    e(i-1) = e(i)
  enddo
  e(n) = 0

  do l=1,n
    iter = 0
 1  do m=l,n-1
      dd = abs(d(m))+abs(d(m+1))
      if (abs(e(m))+dd.eq.dd) goto 2
    enddo
    m=n
 2  if (m.ne.l) then
      iter = iter + 1
      g = (d(l+1)-d(l))/(2.*e(l))
      r = sqrt(g**2+1.)
      sgng = sign(1.,g)
      if (sgng.eq.0.) then 
        sgng = 1.
      endif
      g = d(m)-d(l)+e(l)/(g+r*sgng)
      s = 1.
      c = 1.
      p = 0.

      do i=m-1,l,-1
        f = s*e(i)
        b = c*e(i)
        r = sqrt(f**2+g**2)
        e(i+1) = r

        if (r.eq.0.) then
          d(i+1) = d(i+1) - p
          e(m) = 0.
          goto 1
        endif
        s = f/r
        c = g/r
        g = d(i+1)-p
        r = (d(i)-g)*s+2.*c*b
        d(i+1) = g+(s*r)
        p = s*r
        g = c*r-b
      enddo
    
      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.
      goto 1
    endif
  enddo

  return
end subroutine


