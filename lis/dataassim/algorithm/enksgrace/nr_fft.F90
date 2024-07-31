!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! adapted from f77 Numerical Recipes
!
! reichle, 8 Jun 2005
!
! -----------------------------------------------------------------

module nr_fft
  
  implicit none

contains

  SUBROUTINE fourn(nsize,data,nn,ndim,isign)
    
    ! Replaces data by its ndim-dimensional discrete Fourier transform, 
    ! if isign is input as 1. nn(1:ndim) is an integer array containing the 
    ! lengths of each dimension (number of complex values), which MUST all 
    ! be powers of 2. data is a real array of length twice the
    ! product of these lengths, in which the data are stored as in a 
    ! multidimensional complex FORTRAN array. If isign is input as -1, 
    ! data is replaced by its inverse transform times
    ! the product of the lengths of all dimensions.
    
    integer nsize
    INTEGER isign,ndim,nn(ndim)
    REAL data(nsize)
    INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1, &
         k2,n,nprev,nrem,ntot
    REAL tempi,tempr
    DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
    ntot=1
    do idim=1,ndim
       ntot=ntot*nn(idim)
    end do
    nprev=1
    do idim=1,ndim
       n=nn(idim)
       nrem=ntot/(n*nprev)
       ip1=2*nprev
       ip2=ip1*n
       ip3=ip2*nrem
       i2rev=1
       do i2=1,ip2,ip1
          if(i2.lt.i2rev)then
             do i1=i2,i2+ip1-2,2
                do i3=i1,ip3,ip2
                   i3rev=i2rev+i3-i2
                   tempr=data(i3)
                   tempi=data(i3+1)
                   data(i3)=data(i3rev)
                   data(i3+1)=data(i3rev+1)
                   data(i3rev)=tempr
                   data(i3rev+1)=tempi
                end do
             end do
          end if
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
             i2rev=i2rev-ibit
             ibit=ibit/2
             goto 1
          end if
          i2rev=i2rev+ibit
       end do
       ifp1=ip1
2      if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,ip1
             do i1=i3,i3+ip1-2,2
                do i2=i1,ip3,ifp2
                   k1=i2
                   k2=k1+ifp1
                   tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                   tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                   data(k2)=data(k1)-tempr
                   data(k2+1)=data(k1+1)-tempi
                   data(k1)=data(k1)+tempr
                   data(k1+1)=data(k1+1)+tempi
                end do
             end do
             wtemp=wr
             wr=wr*wpr-wi*wpi+wr
             wi=wi*wpr+wtemp*wpi+wi
          end do
          ifp1=ifp2
          goto 2
       endif
       nprev=n*nprev
    end do
    return
  END SUBROUTINE fourn
  
end module nr_fft


! ========================== EOF ================================


