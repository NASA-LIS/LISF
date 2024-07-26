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
! !ROUTINE: clsmf25_relayer_snow
! \label{clsmf25_relayer_snow}
!
! !REVISION HISTORY:
!  10 Sep 2015 : Adapted from the work of Manuela Girotto (GMAO)
!
! !INTERFACE:
subroutine clsmf25_relayer_snow(t, htsnn, wesn, sndz)
! !USES:
  use LIS_coreMod
  use clsmf25_constants 
!EOP
  
      implicit none
      integer :: t
!      integer, parameter ::  N_snow = 3
      real, intent(inout) :: htsnn(N_snow),wesn(N_snow),sndz(N_snow)
      
      real, dimension(size(sndz),2)   :: ds
      real, dimension(size(sndz)+1)   :: sdold, sdnew
      real, dimension(size(sndz)+1,2) :: h, s
      
      integer :: i, j
      
!      real, parameter :: dz1max = 0.08   ! [m]
!      real, parameter :: wemin  = 13.0   ! [kg m-2]
      real, parameter :: small  = 1.e-20 
      real :: areasc,dz
      
!**** Initialize some variables.
      
      h  = 0.
      s  = 0.
      ds = 0.
      dz = 0.
      
      areasc = min(sum(wesn)/wemin,1.)
      
!**** Compute specific heat & water contents of old layers.

      do i=1,N_snow
         if (sndz(i) > 0.) then
            h(i,1) = htsnn(i)/sndz(i)
            h(i,2) =  wesn(i)/sndz(i)
         endif
      enddo
      
!**** Obtain old & new layer thicknesses & boundaries.
      
      sdold = 0.
      sdnew = 0.
      
      do i=N_snow,1,-1
         sdold(i) = sdold(i+1) + sndz(i)
      enddo
      
      sndz = sdold(1)/float(N_snow)
!      if(t.eq.1401.and.LIS_localPet.eq.140) then 
!         print*, 'rel1 ',sndz
!      endif
      if(sndz(1) > dz1max) then
         sndz(2:) = (sdold(1)-dz1max)/float(N_snow-1)
         sndz(1)  = dz1max
      endif
!      if(t.eq.1401.and.LIS_localPet.eq.140) then 
!         print*, 'rel2 ',sndz
!      endif
      
      do i=N_snow,1,-1
         sdnew(i) = sdnew(i+1) + sndz(i)
      enddo
      
!**** Since the snow boundary has moved, redistribute heat  
!     contents & water equivalents of old to new snow layers.
                                                              
      do i=1,N_snow
         
         j = i
         dz=sdnew(i+1)-sdold(i+1)
         if(dz < 0.) j = i + 1
         s(i+1,:) = h(j,:)*dz
         ds(i,:)  = s(i,:) - s(i+1,:)
      enddo

      htsnn = htsnn + ds(:,1)
      wesn  = wesn  + ds(:,2) 
      
!      if(sum(wesn) < wemin) sndz = sndz /(areasc + small)

      return

    end subroutine clsmf25_relayer_snow


#if 0 
  implicit none
!  integer, parameter ::  N_snow = 3
  real, intent(inout) :: htsnn(N_snow),wesn(N_snow),sndz(N_snow)
  
  real, dimension(size(sndz),2)   :: ds
  real, dimension(size(sndz)+1)   :: sdold, sdnew
  real, dimension(size(sndz)+1,2) :: h, s
  
  integer :: i, j
  
  !      real, parameter :: dz1max = 0.08   ! [m]
  !      real, parameter :: wemin  = 13.0   ! [kg m-2]
  real, parameter :: small  = 1.e-20 
  real :: areasc,dz
  
  !**** Initialize some variables.
  
  h  = 0.
  s  = 0.
  ds = 0.
  dz = 0.
  
  areasc = min(sum(wesn)/wemin,1.)
  
  !**** Compute specific heat & water contents of old layers.
  
  do i=1,N_snow
     if (sndz(i) > 0.) then
        h(i,1) = htsnn(i)/sndz(i)
        h(i,2) =  wesn(i)/sndz(i)
     endif
  enddo
  
  !**** Obtain old & new layer thicknesses & boundaries.
  
  sdold = 0.
  sdnew = 0.
  
  do i=N_snow,1,-1
     sdold(i) = sdold(i+1) + sndz(i)
  enddo
  
  sndz = sdold(1)/float(N_snow)
  if(sndz(1) > dz1max) then
     sndz(2:) = (sdold(1)-dz1max)/float(N_snow-1)
     sndz(1)  = dz1max
  endif
  
  do i=N_snow,1,-1
     sdnew(i) = sdnew(i+1) + sndz(i)
  enddo
  
  !**** Since the snow boundary has moved, redistribute heat  
  !     contents & water equivalents of old to new snow layers.
  
  do i=1,N_snow
     
     j = i
     dz=sdnew(i+1)-sdold(i+1)
     if(dz < 0.) j = i + 1
     s(i+1,:) = h(j,:)*dz
     ds(i,:)  = s(i,:) - s(i+1,:)
  enddo
  
  htsnn = htsnn + ds(:,1)
  wesn  = wesn  + ds(:,2) 
  
!! MG September 15th 
!! Commenting the following out because of the new adjustment of SD [using WEMIN conditions]
!!      if(sum(wesn) < wemin) sndz = sndz /(areasc + small)
  return
  
end subroutine clsmf25_relayer_snow

#endif
