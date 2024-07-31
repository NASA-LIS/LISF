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
! !ROUTINE: mapsib2umd
!  \label{mapsib2umd}
!
! !REVISION HISTORY:
!  13 Apr 2001: Urszula Jambor; Initial code, based on calc_albedo.f scheme
!  03 Feb 2002: Jon Gottschalck; Added sections to use Koster tilespace 
!               (may be better to put the tile space distinction in 
!               umd_sibalb.F90 at a later time)
! !INTERFACE:
 subroutine mapsib2umd ()
! !USES: 
  use sibalb_module    ! SiB-based coefficients for albedo calculation.

  implicit none
!
! DESCRIPTION:
!  This routine assigns then maps 7 SiB vegetation types to 13 UMD 
!  vegetation types, defining 4 coefficient arrays.
!
!  **NOTE: current form of this routine only has mapping based on North 
!  America, and not also for the globe.**
!
!EOP

!===  Begin variable definitions 

  integer :: i, j

  real :: ALVDRS = 0.100 ! vis direct solar rad soil albedo
  real :: ALIDRS = 0.200 ! IR, direct solar rad. soil albedo
  real :: ALVDRD = 0.300 ! vis direct solar rad desert albedo
  real :: ALIDRD = 0.350 ! IR, direct solar rad.desert albedo
  real :: ALVDRI = 0.700 ! vis direct solar rad ice albedo
  real :: ALIDRI = 0.700 ! IR, direct solar rad. ice albedo

!=== End of variable definitions

!=== Begin assigning coefficients ========================================

!=== (Data statements for ALVDR described in full; data statements for
!===  other constants follow same framework.)

!===   BROADLEAF EVERGREEN (ITYP=1); GREEN=0.33; LAI: .5-7
  sib%ALVDRold(1:14,1,1) = (/ 0.0808, 0.0796, 0.0792, 0.0790,     &
                              (0.0789, I=1,10) /)

!===   BROADLEAF EVERGREEN (ITYP=1); GREEN=0.67; LAI: .5-7
  sib%ALVDRold(1:14,2,1) = (/ 0.0788, 0.0775, 0.0771, 0.0769,     &
                              (0.0768, I=1,10) /)

!===   BROADLEAF DECIDUOUS (ITYP=2); GREEN=0.33; LAI: .5-7
  sib%ALVDRold(1:14,1,2) = (/ 0.0803, 0.0790, 0.0785, 0.0784,     &
                              (0.0783, I=1,3),(0.0782, I=1,7) /)

!===   BROADLEAF DECIDUOUS (ITYP=2); GREEN=0.67; LAI: .5-7
  sib%ALVDRold(1:14,2,2) = (/ 0.0782, 0.0770, 0.0765, 0.0763,     &
                              (0.0762, I=1,10) /)

!===   NEEDLELEAF (ITYP=3); GREEN=0.33; LAI=.5-7
  sib%ALVDRold(1:14,1,3) = (/ 0.0758, 0.0746, 0.0742, 0.0740,     &
                              (0.0739, I=1,10) /)

!===   NEEDLELEAF (ITYP=3); GREEN=0.67; LAI=.5-7
  sib%ALVDRold(1:14,2,3) = (/ 0.0683, 0.0672, 0.0667, 0.0665,     &
                              0.0665, (0.0664, I=1,9) /)

!===   GROUNDCOVER (ITYP=4); GREEN=0.33; LAI=.5-7
  sib%ALVDRold(1:14,1,4) = (/ 0.2436, 0.2470, 0.2486, 0.2494,     &
                              0.2498, 0.2500, 0.2501, 0.2501,     &
                              (0.2502, I=1,6) /)

!===   GROUNDCOVER (ITYP=4); GREEN=0.67; LAI=.5-7
  sib%ALVDRold(1:14,2,4) = (/ (0.1637, I=1,14) /)

!===   BROADLEAF SHRUBS (ITYP=5); GREEN=0.33,LAI=.5-7
  sib%ALVDRold(1:14,1,5) = (/ 0.0807, 0.0798, 0.0794, 0.0792,     &
                              0.0792, (0.0791, I=1,9) /)

!===   BROADLEAF SHRUBS (ITYP=5); GREEN=0.67,LAI=.5-7
  sib%ALVDRold(1:14,2,5) = (/ 0.0787, 0.0777, 0.0772, 0.0771,     &
                              (0.0770, I=1,10) /)

!===   DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.33,LAI=.5-7
  sib%ALVDRold(1:14,1,6) = (/ 0.0802, 0.0791, 0.0787, 0.0786,     &
                              (0.0785, I=1,10) /)

!===   DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.67,LAI=.5-7
  sib%ALVDRold(1:14,2,6) = (/ 0.0781, 0.0771, 0.0767, 0.0765,     &
                              0.0765, (0.0764, I=1,9) /)

!===   BARE SOIL
  sib%ALVDRold(1:14,1,7) = (/ (ALVDRS, I=1,14) /)
  sib%ALVDRold(1:14,2,7) = (/ (ALVDRS, I=1,14) /)

!===   DESERT
  sib%ALVDRold(1:14,1,8) = (/ (ALVDRD, I=1,14) /)
  sib%ALVDRold(1:14,2,8) = (/ (ALVDRD, I=1,14) /)

!===   ICE
  sib%ALVDRold(1:14,1,9) = (/ (ALVDRI, I=1,14) /)
  sib%ALVDRold(1:14,2,9) = (/ (ALVDRI, I=1,14) /)


!  IF (LDAS%KOSTER .EQ. 1) THEN

!   do j=1,2
!    do i=1,14
!     sib%ALVDR(i,j,1)=sib%ALVDRold(i,j,1)
!     sib%ALVDR(i,j,2)=sib%ALVDRold(i,j,2)
!     sib%ALVDR(i,j,3)=sib%ALVDRold(i,j,3)
!     sib%ALVDR(i,j,4)=sib%ALVDRold(i,j,4)
!     sib%ALVDR(i,j,5)=sib%ALVDRold(i,j,5)
!     sib%ALVDR(i,j,6)=sib%ALVDRold(i,j,6)
!     sib%ALVDR(i,j,7)=sib%ALVDRold(i,j,7)
!     sib%ALVDR(i,j,8)=sib%ALVDRold(i,j,8)
!     sib%ALVDR(i,j,9)=sib%ALVDRold(i,j,9)
!    enddo
!   enddo

! ELSE

  do j=1,2
     do i=1,14
        sib%ALVDR(i,j,1)=sib%ALVDRold(i,j,3)
        
        sib%ALVDR(i,j,2)=sib%ALVDRold(i,j,1)
        
        sib%ALVDR(i,j,3)=sib%ALVDRold(i,j,3)
        
        sib%ALVDR(i,j,4)=sib%ALVDRold(i,j,2)
        
        sib%ALVDR(i,j,5)=(.5*sib%ALVDRold(i,j,2) +    &
                          .5*sib%ALVDRold(i,j,3))
        
        sib%ALVDR(i,j,6)=(.5702*sib%ALVDRold(i,j,3) + &
                          .2200*sib%ALVDRold(i,j,4) + &
                          .2098*sib%ALVDRold(i,j,2))
        
        sib%ALVDR(i,j,7)=(.6372*sib%ALVDRold(i,j,4) + &
                          .1955*sib%ALVDRold(i,j,3) + &
                          .1564*sib%ALVDRold(i,j,2) + &
                          .0109*sib%ALVDRold(i,j,5))
           
        sib%ALVDR(i,j,8)=(.4403*sib%ALVDRold(i,j,8) + &
                          .4365*sib%ALVDRold(i,j,4) + &
                          .0743*sib%ALVDRold(i,j,6) + &
                          .0489*sib%ALVDRold(i,j,5))
        
        sib%ALVDR(i,j,9)=(.8506*sib%ALVDRold(i,j,8) + &
                          .0950*sib%ALVDRold(i,j,5) + &
                          .0399*sib%ALVDRold(i,j,4) + &
                          .0145*sib%ALVDRold(i,j,6))
        
        sib%ALVDR(i,j,10)=sib%ALVDRold(i,j,4)
        
        sib%ALVDR(i,j,11)=sib%ALVDRold(i,j,4)
        
        sib%ALVDR(i,j,12)=sib%ALVDRold(i,j,7)
        
        sib%ALVDR(i,j,13)=(.7114*sib%ALVDRold(i,j,4) + &
                           .1055*sib%ALVDRold(i,j,2) + &
                           .0723*sib%ALVDRold(i,j,3) + &
                           .0526*sib%ALVDRold(i,j,8) + &
                           .0178*sib%ALVDRold(i,j,6) + &
                           .0077*sib%ALVDRold(i,j,5) + &
                           .0327*sib%ALVDRold(i,j,7))
     end do !i
  end do !j
  
!  ENDIF

!****
!**** -------------------------------------------------

  sib%BTVDRold(1:14,1,1) = (/ 0.0153, 0.0372, 0.0506, 0.0587, 0.0630, &
                              0.0652, 0.0663, 0.0668, 0.0671, 0.0672, &
                              (0.0673, I=1,4) /)
  sib%BTVDRold(1:14,2,1) = (/ 0.0135, 0.0354, 0.0487, 0.0568, 0.0611, &
                              0.0633, 0.0644, 0.0650, 0.0652, 0.0654, &
                              0.0654, (0.0655, I=1,3) /)
  
  sib%BTVDRold(1:14,1,2) = (/ 0.0148, 0.0357, 0.0462, 0.0524, 0.0554, &
                              0.0569, 0.0576, 0.0579, 0.0580, 0.0581, &
                              0.0581, (0.0582, I=1,3) /)
  sib%BTVDRold(1:14,2,2) = (/ 0.0131, 0.0342, 0.0446, 0.0508, 0.0539, &
                              0.0554, 0.0560, 0.0564, 0.0565,         &
                              (0.0566, I=1,5) /)
  
  sib%BTVDRold(1:14,1,3) = (/ 0.0108, 0.0334, 0.0478, 0.0571, 0.0624, &
                              0.0652, 0.0666, 0.0673, 0.0677, 0.0679, &
                              (0.0680, I=1,4) /)
  sib%BTVDRold(1:14,2,3) = (/ 0.0034, 0.0272, 0.0408, 0.0501, 0.0554, &
                              0.0582, 0.0597, 0.0604, 0.0608, 0.0610, &
                              (0.0611, I=1,4) /)

  sib%BTVDRold(1:14,1,4) = (/ 0.2050, 0.2524, 0.2799, 0.2947, 0.3022, &
                              0.3059, 0.3076, 0.3085, 0.3088, 0.3090, &
                              (0.3091, I=1,4) /)
  sib%BTVDRold(1:14,2,4) = (/ 0.1084, 0.1404, 0.1617, 0.1754, 0.1837, &
                              0.1887, 0.1915, 0.1931, 0.1940, 0.1946, &
                              0.1948, 0.1950, 0.1951, 0.1951 /)
  
  sib%BTVDRold(1:14,1,5) = (/ 0.0203, 0.0406, 0.0548, 0.0632, 0.0679, &
                              0.0703, 0.0716, 0.0722, 0.0726, 0.0727, &
                              0.0728, 0.0728, 0.0728, 0.0729 /)
  sib%BTVDRold(1:14,2,5) = (/ 0.0184, 0.0385, 0.0526, 0.0611, 0.0658, &
                              0.0683, 0.0696, 0.0702, 0.0705, 0.0707, &
                              (0.0708, I=1,4) /)
  
  sib%BTVDRold(1:14,1,6) = (/ 0.0199, 0.0388, 0.0494, 0.0554, 0.0584, &
                              0.0599, 0.0606, 0.0609, 0.0611,         &
                              (0.0612, I=1,5) /)
  sib%BTVDRold(1:14,2,6) = (/ 0.0181, 0.0371, 0.0476, 0.0537, 0.0568, &
                              0.0583, 0.0590, 0.0593, 0.0595, 0.0595, &
                              (0.0596, I=1,4) /)
  
  sib%BTVDRold(1:14,1,7) = (/ (0., I=1,14) /)
  sib%BTVDRold(1:14,2,7) = (/ (0., I=1,14) /)
  
  sib%BTVDRold(1:14,1,8) = (/ (0., I=1,14) /)
  sib%BTVDRold(1:14,2,8) = (/ (0., I=1,14) /)
  
  sib%BTVDRold(1:14,1,9) = (/ (0., I=1,14) /)
  sib%BTVDRold(1:14,2,9) = (/ (0., I=1,14) /)
  
!  IF (LDAS%KOSTER .EQ. 1) THEN

!   do j=1,2
!    do i=1,14
!     sib%BTVDR(i,j,1)=sib%BTVDRold(i,j,1)
!     sib%BTVDR(i,j,2)=sib%BTVDRold(i,j,2)
!     sib%BTVDR(i,j,3)=sib%BTVDRold(i,j,3)
!     sib%BTVDR(i,j,4)=sib%BTVDRold(i,j,4)
!     sib%BTVDR(i,j,5)=sib%BTVDRold(i,j,5)
!     sib%BTVDR(i,j,6)=sib%BTVDRold(i,j,6)
!     sib%BTVDR(i,j,7)=sib%BTVDRold(i,j,7)
!     sib%BTVDR(i,j,8)=sib%BTVDRold(i,j,8)
!     sib%BTVDR(i,j,9)=sib%BTVDRold(i,j,9)
!    enddo
!   enddo

! ELSE

  do j=1,2
     do i=1,14
        sib%BTVDR(i,j,1)=sib%BTVDRold(i,j,3)
        
        sib%BTVDR(i,j,2)=sib%BTVDRold(i,j,1)
        
        sib%BTVDR(i,j,3)=sib%BTVDRold(i,j,3)
        
        sib%BTVDR(i,j,4)=sib%BTVDRold(i,j,2)
        
        sib%BTVDR(i,j,5)=(.5*sib%BTVDRold(i,j,2) +    &
                          .5*sib%BTVDRold(i,j,3))
        
        sib%BTVDR(i,j,6)=(.5702*sib%BTVDRold(i,j,3) + &
                          .2200*sib%BTVDRold(i,j,4) + &
                          .2098*sib%BTVDRold(i,j,2))
        
        sib%BTVDR(i,j,7)=(.6372*sib%BTVDRold(i,j,4) + &
                          .1955*sib%BTVDRold(i,j,3) + &
                          .1564*sib%BTVDRold(i,j,2) + &
                          .0109*sib%BTVDRold(i,j,5))
        
        sib%BTVDR(i,j,8)=(.4403*sib%BTVDRold(i,j,8) + &
                          .4365*sib%BTVDRold(i,j,4) + &
                          .0743*sib%BTVDRold(i,j,6) + &
                          .0489*sib%BTVDRold(i,j,5))
        
        sib%BTVDR(i,j,9)=(.8506*sib%BTVDRold(i,j,8) + &
                          .0950*sib%BTVDRold(i,j,5) + &
                          .0399*sib%BTVDRold(i,j,4) + &
                          .0145*sib%BTVDRold(i,j,6))
        
        sib%BTVDR(i,j,10)=sib%BTVDRold(i,j,4)
        
        sib%BTVDR(i,j,11)=sib%BTVDRold(i,j,4)
        
        sib%BTVDR(i,j,12)=sib%BTVDRold(i,j,7)
        
        sib%BTVDR(i,j,13)=(.7114*sib%BTVDRold(i,j,4) + & 
                           .1055*sib%BTVDRold(i,j,2) + &
                           .0723*sib%BTVDRold(i,j,3) + &
                           .0526*sib%BTVDRold(i,j,8) + &
                           .0178*sib%BTVDRold(i,j,6) + &
                           .0077*sib%BTVDRold(i,j,5) + &
                           .0327*sib%BTVDRold(i,j,7))
     end do !i
  end do !j
  
!  ENDIF

!****
!**** -----------------------------------------------------------

  sib%GMVDRold(1:14,1,1) = (/ 0.0814, 0.1361, 0.2078, 0.2650, 0.2986, &
                              0.3169, 0.3265, 0.3313, 0.3337, 0.3348, &
                              0.3354, 0.3357, 0.3358, 0.3358 /)
  sib%GMVDRold(1:14,2,1) = (/ 0.0760, 0.1336, 0.2034, 0.2622, 0.2969, &
                              0.3159, 0.3259, 0.3309, 0.3333, 0.3346, &
                              0.3352, 0.3354, 0.3356, 0.3356 /) 
  
  sib%GMVDRold(1:14,1,2) = (/ 0.0834, 0.1252, 0.1558, 0.1927, 0.2131, &
                              0.2237, 0.2290, 0.2315, 0.2327, 0.2332, &
                              0.2335, 0.2336, 0.2336, 0.2337 /)
  sib%GMVDRold(1:14,2,2) = (/ 0.0789, 0.1235, 0.1531, 0.1912, 0.2122, &
                              0.2232, 0.2286, 0.2312, 0.2324, 0.2330, &
                              0.2333, 0.2334, 0.2335, 0.2335 /) 
  
  sib%GMVDRold(1:14,1,3) = (/ 0.0647, 0.1342, 0.2215, 0.2968, 0.3432, &
                              0.3696, 0.3838, 0.3912, 0.3950, 0.3968, &
                              0.3978, 0.3982, 0.3984, 0.3985 /)
  sib%GMVDRold(1:14,2,3) = (/ 0.0258, 0.1227, 0.1999, 0.2825, 0.3339, &
                              0.3634, 0.3794, 0.3877, 0.3919, 0.3940, &
                              0.3950, 0.3956, 0.3958, 0.3959 /)
  
  sib%GMVDRold(1:14,1,4) = (/ 0.3371, 0.5762, 0.7159, 0.7927, 0.8324, &
                              0.8526, 0.8624, 0.8671, 0.8693, 0.8704, &
                              0.8709, 0.8710, 0.8712, 0.8712 /)
  sib%GMVDRold(1:14,2,4) = (/ 0.2634, 0.4375, 0.5532, 0.6291, 0.6763, &
                              0.7048, 0.7213, 0.7310, 0.7363, 0.7395, &
                              0.7411, 0.7420, 0.7426, 0.7428 /)
  
  sib%GMVDRold(1:14,1,5) = (/ 0.0971, 0.1544, 0.2511, 0.3157, 0.3548, &
                              0.3768, 0.3886, 0.3948, 0.3978, 0.3994, &
                              0.4001, 0.4006, 0.4007, 0.4008 /)
  sib%GMVDRold(1:14,2,5) = (/ 0.0924, 0.1470, 0.2458, 0.3123, 0.3527, &
                              0.3756, 0.3877, 0.3942, 0.3974, 0.3990, &
                              0.3998, 0.4002, 0.4004, 0.4005 /)
  
  sib%GMVDRold(1:14,1,6) = (/ 0.0970, 0.1355, 0.1841, 0.2230, 0.2447, &
                              0.2561, 0.2617, 0.2645, 0.2658, 0.2664, &
                              0.2667, (0.2669, I=1,3) /)
  sib%GMVDRold(1:14,2,6) = (/ 0.0934, 0.1337, 0.1812, 0.2213, 0.2437, &
                              0.2554, 0.2613, 0.2642, 0.2656, 0.2662, &
                              0.2665, 0.2667, 0.2667, 0.2668 /)
  
  sib%GMVDRold(1:14,1,7) = (/ (1., I=1,14) /)
  sib%GMVDRold(1:14,2,7) = (/ (1., I=1,14) /)
  
  sib%GMVDRold(1:14,1,8) = (/ (1., I=1,14) /)
  sib%GMVDRold(1:14,2,8) = (/ (1., I=1,14) /)
  
  sib%GMVDRold(1:14,1,9) = (/ (1., I=1,14) /)
  sib%GMVDRold(1:14,2,9) = (/ (1., I=1,14) /)

!  IF (LDAS%KOSTER .EQ. 1) THEN

!   do j=1,2
!    do i=1,14
!     sib%GMVDR(i,j,1)=sib%GMVDRold(i,j,1)
!     sib%GMVDR(i,j,2)=sib%GMVDRold(i,j,2)
!     sib%GMVDR(i,j,3)=sib%GMVDRold(i,j,3)
!     sib%GMVDR(i,j,4)=sib%GMVDRold(i,j,4)
!     sib%GMVDR(i,j,5)=sib%GMVDRold(i,j,5)
!     sib%GMVDR(i,j,6)=sib%GMVDRold(i,j,6)
!     sib%GMVDR(i,j,7)=sib%GMVDRold(i,j,7)
!     sib%GMVDR(i,j,8)=sib%GMVDRold(i,j,8)
!     sib%GMVDR(i,j,9)=sib%GMVDRold(i,j,9)
!    enddo
!   enddo

! ELSE
    
  do j=1,2
     do i=1,14
        sib%GMVDR(i,j,1)=sib%GMVDRold(i,j,3)
        
        sib%GMVDR(i,j,2)=sib%GMVDRold(i,j,1)
        
        sib%GMVDR(i,j,3)=sib%GMVDRold(i,j,3)
        
        sib%GMVDR(i,j,4)=sib%GMVDRold(i,j,2)
        
        sib%GMVDR(i,j,5)=(.5*sib%GMVDRold(i,j,2) +  &
                          .5*sib%GMVDRold(i,j,3))
        
        sib%GMVDR(i,j,6)=(.5702*sib%GMVDRold(i,j,3) +  &
                          .2200*sib%GMVDRold(i,j,4) +  &
                          .2098*sib%GMVDRold(i,j,2))
        
        sib%GMVDR(i,j,7)=(.6372*sib%GMVDRold(i,j,4) +  &
                          .1955*sib%GMVDRold(i,j,3) +  &
                          .1564*sib%GMVDRold(i,j,2) +  &
                          .0109*sib%GMVDRold(i,j,5))
        
        sib%GMVDR(i,j,8)=(.4403*sib%GMVDRold(i,j,8) +  &
                          .4365*sib%GMVDRold(i,j,4) +  &
                          .0743*sib%GMVDRold(i,j,6) +  &
                          .0489*sib%GMVDRold(i,j,5))
        
        sib%GMVDR(i,j,9)=(.8506*sib%GMVDRold(i,j,8) +  &
                          .0950*sib%GMVDRold(i,j,5) +  &
                          .0399*sib%GMVDRold(i,j,4) +  &
                          .0145*sib%GMVDRold(i,j,6))
        
        sib%GMVDR(i,j,10)=sib%GMVDRold(i,j,4)
        
        sib%GMVDR(i,j,11)=sib%GMVDRold(i,j,4)
        
        sib%GMVDR(i,j,12)=sib%GMVDRold(i,j,7)
        
        sib%GMVDR(i,j,13)=(.7114*sib%GMVDRold(i,j,4) +  &
                           .1055*sib%GMVDRold(i,j,2) +  &
                           .0723*sib%GMVDRold(i,j,3) +  &
                           .0526*sib%GMVDRold(i,j,8) +  &
                           .0178*sib%GMVDRold(i,j,6) +  &
                           .0077*sib%GMVDRold(i,j,5) +  &
                           .0327*sib%GMVDRold(i,j,7))
     end do !i
  end do !j
 
!  ENDIF

!****
!****  -----------------------------------------------------------

  sib%ALIDRold(1:14,1,1) = (/ 0.2867, 0.2840, 0.2828, 0.2822, 0.2819, &
                              0.2818, 0.2817, 0.2817, (0.2816, I=1,6) /)
  sib%ALIDRold(1:14,2,1) = (/ 0.3564, 0.3573, 0.3577, 0.3580, 0.3581, &
                              0.3581, (0.3582, I=1,8) /)

  sib%ALIDRold(1:14,1,2) = (/ 0.2848, 0.2819, 0.2804, 0.2798, 0.2795, &
                              0.2793, 0.2793, (0.2792, I=1,7) /)
  sib%ALIDRold(1:14,2,2) = (/ 0.3544, 0.3550, 0.3553, 0.3555, 0.3555, &
                              (0.3556, I=1,9) /)

  sib%ALIDRold(1:14,1,3) = (/ 0.2350, 0.2311, 0.2293, 0.2285, 0.2281, &
                              0.2280, (0.2279, I=1,8) /)
  sib%ALIDRold(1:14,2,3) = (/ 0.2474, 0.2436, 0.2418, 0.2410, 0.2406, &
                              0.2405, (0.2404, I=1,3),(0.2403, I=1,5) /)

  sib%ALIDRold(1:14,1,4) = (/ 0.5816, 0.6157, 0.6391, 0.6556, 0.6673, &
                              0.6758, 0.6820, 0.6866, 0.6899, 0.6924, &
                              0.6943, 0.6956, 0.6966, 0.6974 /)
  sib%ALIDRold(1:14,2,4) = (/ 0.5489, 0.5770, 0.5955, 0.6079, 0.6163, &
                              0.6221, 0.6261, 0.6288, 0.6308, 0.6321, &
                              0.6330, 0.6337, 0.6341, 0.6344 /)

  sib%ALIDRold(1:14,1,5) = (/ 0.2845, 0.2837, 0.2832, 0.2831, 0.2830, &
                              (0.2829, I=1,9) /)
  sib%ALIDRold(1:14,2,5) = (/ 0.3532, 0.3562, 0.3578, 0.3586, 0.3590, &
                              0.3592, 0.3594, 0.3594, 0.3594,         &
                              (0.3595, I=1,5) /)

  sib%ALIDRold(1:14,1,6) = (/ 0.2825, 0.2812, 0.2806, 0.2803, 0.2802, &
                              (0.2801, I=1,9) /)
  sib%ALIDRold(1:14,2,6) = (/ 0.3512, 0.3538, 0.3552, 0.3559, 0.3562, &
                              0.3564, 0.3565, 0.3565, (0.3566, I=1,6) /)

  sib%ALIDRold(1:14,1,7) = (/ (ALIDRS, I=1,14) /)
  sib%ALIDRold(1:14,2,7) = (/ (ALIDRS, I=1,14) /)
  
  sib%ALIDRold(1:14,1,8) = (/ (ALIDRD, I=1,14) /)
  sib%ALIDRold(1:14,2,8) = (/ (ALIDRD, I=1,14) /)

  sib%ALIDRold(1:14,1,9) = (/ (ALIDRI, I=1,14) /)
  sib%ALIDRold(1:14,2,9) = (/ (ALIDRI, I=1,14) /)

!  IF (LDAS%KOSTER .EQ. 1) THEN

!   do j=1,2
!    do i=1,14
!     sib%ALIDR(i,j,1)=sib%ALIDRold(i,j,1)
!     sib%ALIDR(i,j,2)=sib%ALIDRold(i,j,2)
!     sib%ALIDR(i,j,3)=sib%ALIDRold(i,j,3)
!     sib%ALIDR(i,j,4)=sib%ALIDRold(i,j,4)
!     sib%ALIDR(i,j,5)=sib%ALIDRold(i,j,5)
!     sib%ALIDR(i,j,6)=sib%ALIDRold(i,j,6)
!     sib%ALIDR(i,j,7)=sib%ALIDRold(i,j,7)
!     sib%ALIDR(i,j,8)=sib%ALIDRold(i,j,8)
!     sib%ALIDR(i,j,9)=sib%ALIDRold(i,j,9)
!    enddo
!   enddo

! ELSE

  do j=1,2
     do i=1,14
        sib%ALIDR(i,j,1)=sib%ALIDRold(i,j,3)
        
        sib%ALIDR(i,j,2)=sib%ALIDRold(i,j,1)
        
        sib%ALIDR(i,j,3)=sib%ALIDRold(i,j,3)
        
        sib%ALIDR(i,j,4)=sib%ALIDRold(i,j,2)
        
        sib%ALIDR(i,j,5)=(.5*sib%ALIDRold(i,j,2) +    &
                          .5*sib%ALIDRold(i,j,3))
        
        sib%ALIDR(i,j,6)=(.5702*sib%ALIDRold(i,j,3) + &
                          .2200*sib%ALIDRold(i,j,4) + &
                          .2098*sib%ALIDRold(i,j,2))
        
        sib%ALIDR(i,j,7)=(.6372*sib%ALIDRold(i,j,4) + &
                          .1955*sib%ALIDRold(i,j,3) + &
                          .1564*sib%ALIDRold(i,j,2) + &
                          .0109*sib%ALIDRold(i,j,5))
        
        sib%ALIDR(i,j,8)=(.4403*sib%ALIDRold(i,j,8) + &
                          .4365*sib%ALIDRold(i,j,4) + &
                          .0743*sib%ALIDRold(i,j,6) + &
                          .0489*sib%ALIDRold(i,j,5))
        
        sib%ALIDR(i,j,9)=(.8506*sib%ALIDRold(i,j,8) + &
                          .0950*sib%ALIDRold(i,j,5) + &
                          .0399*sib%ALIDRold(i,j,4) + &
                          .0145*sib%ALIDRold(i,j,6))
        
        sib%ALIDR(i,j,10)=sib%ALIDRold(i,j,4)
        
        sib%ALIDR(i,j,11)=sib%ALIDRold(i,j,4)
        
        sib%ALIDR(i,j,12)=sib%ALIDRold(i,j,7)
        
        sib%ALIDR(i,j,13)=(.7114*sib%ALIDRold(i,j,4) + &
                           .1055*sib%ALIDRold(i,j,2) + &
                           .0723*sib%ALIDRold(i,j,3) + &
                           .0526*sib%ALIDRold(i,j,8) + &
                           .0178*sib%ALIDRold(i,j,6) + &
                           .0077*sib%ALIDRold(i,j,5) + &
                           .0327*sib%ALIDRold(i,j,7))
     end do !i
  end do !j

!  ENDIF

!****
!**** -----------------------------------------------------------

  sib%BTIDRold(1:14,1,1) = (/ 0.1291, 0.1707, 0.1969, 0.2125, 0.2216, &
                              0.2267, 0.2295, 0.2311, 0.2319, 0.2323, &
                              0.2326, 0.2327, 0.2327, 0.2328 /)
  sib%BTIDRold(1:14,2,1) = (/ 0.1939, 0.2357, 0.2598, 0.2735, 0.2810, &
                              0.2851, 0.2874, 0.2885, 0.2892, 0.2895, &
                              0.2897, (0.2898, I=1,3) /)
  
  sib%BTIDRold(1:14,1,2) = (/ 0.1217, 0.1522, 0.1713, 0.1820, 0.1879, &
                              0.1910, 0.1926, 0.1935, 0.1939, 0.1942, &
                              0.1943, 0.1943, 0.1944, 0.1944 /)
  sib%BTIDRold(1:14,2,2) = (/ 0.1781, 0.2067, 0.2221, 0.2301, 0.2342, &
                              0.2363, 0.2374, 0.2379, 0.2382, 0.2383, &
                              0.2384, 0.2384, 0.2385, 0.2385 /)
  
  sib%BTIDRold(1:14,1,3) = (/ 0.0846, 0.1299, 0.1614, 0.1814, 0.1935, &
                              0.2004, 0.2043, 0.2064, 0.2076, 0.2082, &
                              0.2085, 0.2087, 0.2087, 0.2088 /)
  sib%BTIDRold(1:14,2,3) = (/ 0.0950, 0.1410, 0.1722, 0.1921, 0.2042, &
                              0.2111, 0.2151, 0.2172, 0.2184, 0.2191, &
                              0.2194, 0.2196, 0.2197, 0.2197 /)
  
  sib%BTIDRold(1:14,1,4) = (/ 0.5256, 0.7444, 0.9908, 1.2700, 1.5680, &
                              1.8505, 2.0767, 2.2211, 2.2808, 2.2774, &
                              2.2362, 2.1779, 2.1160, 2.0564 /)
  sib%BTIDRold(1:14,2,4) = (/ 0.4843, 0.6714, 0.8577, 1.0335, 1.1812, &
                              1.2858, 1.3458, 1.3688, 1.3685, 1.3546, &
                              1.3360, 1.3168, 1.2989, 1.2838 /)
  
  sib%BTIDRold(1:14,1,5) = (/ 0.1498, 0.1930, 0.2201, 0.2364, 0.2460, &
                              0.2514, 0.2544, 0.2560, 0.2569, 0.2574, &
                              0.2577, 0.2578, 0.2579, 0.2579 /)
  sib%BTIDRold(1:14,2,5) = (/ 0.2184, 0.2656, 0.2927, 0.3078, 0.3159, &
                              0.3202, 0.3224, 0.3235, 0.3241, 0.3244, &
                              0.3245, (0.3246, I=1,3) /)
  
  sib%BTIDRold(1:14,1,6) = (/ 0.1369, 0.1681, 0.1860, 0.1958, 0.2010, &
                              0.2038, 0.2053, 0.2060, 0.2064, 0.2066, &
                              0.2067, (0.2068, I=1,3) /)
  sib%BTIDRold(1:14,2,6) = (/ 0.1969, 0.2268, 0.2416, 0.2488, 0.2521, &
                              0.2537, 0.2544, 0.2547, 0.2548,         &
                              (0.2549, I=1,5) /)
  
  sib%BTIDRold(1:14,1,7) = (/ (0., I=1,14) /)
  sib%BTIDRold(1:14,2,7) = (/ (0., I=1,14) /)
  
  sib%BTIDRold(1:14,1,8) = (/ (0., I=1,14) /)
  sib%BTIDRold(1:14,2,8) = (/ (0., I=1,14) /)
  
  sib%BTIDRold(1:14,1,9) = (/ (0., I=1,14) /)
  sib%BTIDRold(1:14,2,9) = (/ (0., I=1,14) /)
  
!  IF (LDAS%KOSTER .EQ. 1) THEN

!   do j=1,2
!    do i=1,14
!     sib%BTIDR(i,j,1)=sib%BTIDRold(i,j,1)
!     sib%BTIDR(i,j,2)=sib%BTIDRold(i,j,2)
!     sib%BTIDR(i,j,3)=sib%BTIDRold(i,j,3)
!     sib%BTIDR(i,j,4)=sib%BTIDRold(i,j,4)
!     sib%BTIDR(i,j,5)=sib%BTIDRold(i,j,5)
!     sib%BTIDR(i,j,6)=sib%BTIDRold(i,j,6)
!     sib%BTIDR(i,j,7)=sib%BTIDRold(i,j,7)
!     sib%BTIDR(i,j,8)=sib%BTIDRold(i,j,8)
!     sib%BTIDR(i,j,9)=sib%BTIDRold(i,j,9)
!    enddo
!   enddo

! ELSE

  do j=1,2
     do i=1,14
        sib%BTIDR(i,j,1)=sib%BTIDRold(i,j,3)
        
        sib%BTIDR(i,j,2)=sib%BTIDRold(i,j,1)
        
        sib%BTIDR(i,j,3)=sib%BTIDRold(i,j,3)
        
        sib%BTIDR(i,j,4)=sib%BTIDRold(i,j,2)
        
        sib%BTIDR(i,j,5)=(.5*sib%BTIDRold(i,j,2) +    &
                          .5*sib%BTIDRold(i,j,3))
        
        sib%BTIDR(i,j,6)=(.5702*sib%BTIDRold(i,j,3) + &
                          .2200*sib%BTIDRold(i,j,4) + &
                          .2098*sib%BTIDRold(i,j,2))
        
        sib%BTIDR(i,j,7)=(.6372*sib%BTIDRold(i,j,4) + &
                          .1955*sib%BTIDRold(i,j,3) + &
                          .1564*sib%BTIDRold(i,j,2) + &
                          .0109*sib%BTIDRold(i,j,5))
        
        sib%BTIDR(i,j,8)=(.4403*sib%BTIDRold(i,j,8) + &
                          .4365*sib%BTIDRold(i,j,4) + &
                          .0743*sib%BTIDRold(i,j,6) + &
                          .0489*sib%BTIDRold(i,j,5))
        
        sib%BTIDR(i,j,9)=(.8506*sib%BTIDRold(i,j,8) + &
                          .0950*sib%BTIDRold(i,j,5) + &
                          .0399*sib%BTIDRold(i,j,4) + &
                          .0145*sib%BTIDRold(i,j,6))
        
        sib%BTIDR(i,j,10)=sib%BTIDRold(i,j,4)
        
        sib%BTIDR(i,j,11)=sib%BTIDRold(i,j,4)
        
        sib%BTIDR(i,j,12)=sib%BTIDRold(i,j,7)
        
        sib%BTIDR(i,j,13)=(.7114*sib%BTIDRold(i,j,4) + &
                           .1055*sib%BTIDRold(i,j,2) + &
                           .0723*sib%BTIDRold(i,j,3) + &
                           .0526*sib%BTIDRold(i,j,8) + &
                           .0178*sib%BTIDRold(i,j,6) + &
                           .0077*sib%BTIDRold(i,j,5) + &
                           .0327*sib%BTIDRold(i,j,7))
     end do !i
  end do !j
  
!  ENDIF

!****
!**** --------------------------------------------------------------

  sib%GMIDRold(1:14,1,1) = (/ 0.1582, 0.2581, 0.3227, 0.3635, 0.3882, &
                              0.4026, 0.4108, 0.4154, 0.4179, 0.4193, &
                              0.4200, 0.4204, 0.4206, 0.4207 /)
  sib%GMIDRold(1:14,2,1) = (/ 0.1934, 0.3141, 0.3818, 0.4200, 0.4415, &
                              0.4533, 0.4598, 0.4633, 0.4651, 0.4662, &
                              0.4667, 0.4671, 0.4672, 0.4672 /)
  
  sib%GMIDRold(1:14,1,2) = (/ 0.1347, 0.1871, 0.2277, 0.2515, 0.2651, &
                              0.2727, 0.2768, 0.2790, 0.2801, 0.2808, &
                              0.2811, 0.2812, 0.2813, 0.2814 /)
  sib%GMIDRold(1:14,2,2) = (/ 0.1440, 0.2217, 0.2629, 0.2839, 0.2947, &
                              0.3003, 0.3031, 0.3046, 0.3054, 0.3058, &
                              0.3060, 0.3061, 0.3061, 0.3062 /)
  
  sib%GMIDRold(1:14,1,3) = (/ 0.1372, 0.2368, 0.3235, 0.3839, 0.4229, &
                              0.4465, 0.4602, 0.4679, 0.4722, 0.4745, &
                              0.4758, 0.4764, 0.4768, 0.4770 /)
  sib%GMIDRold(1:14,2,3) = (/ 0.1435, 0.2524, 0.3370, 0.3955, 0.4332, &
                              0.4563, 0.4697, 0.4773, 0.4815, 0.4839, &
                              0.4851, 0.4858, 0.4861, 0.4863 /)
  
  sib%GMIDRold(1:14,1,4) = (/ 0.4298, 0.9651, 1.6189, 2.4084, 3.2992, &
                              4.1928, 4.9611, 5.5095, 5.8085, 5.9069, &
                              5.8726, 5.7674, 5.6346, 5.4944 /)
  sib%GMIDRold(1:14,2,4) = (/ 0.4167, 0.8974, 1.4160, 1.9414, 2.4147, &
                              2.7803, 3.0202, 3.1468, 3.1954, 3.1932, &
                              3.1676, 3.1328, 3.0958, 3.0625 /)
  
  sib%GMIDRold(1:14,1,5) = (/ 0.1959, 0.3203, 0.3985, 0.4472, 0.4766, &
                              0.4937, 0.5034, 0.5088, 0.5117, 0.5134, &
                              0.5143, 0.5147, 0.5150, 0.5152 /)
  sib%GMIDRold(1:14,2,5) = (/ 0.2328, 0.3859, 0.4734, 0.5227, 0.5498, &
                              0.5644, 0.5720, 0.5761, 0.5781, 0.5792, &
                              0.5797, 0.5800, 0.5802, 0.5802 /)
  
  sib%GMIDRold(1:14,1,6) = (/ 0.1447, 0.2244, 0.2698, 0.2953, 0.3094, &
                              0.3170, 0.3211, 0.3233, 0.3244, 0.3250, &
                              0.3253, 0.3255, 0.3256, 0.3256 /)
  sib%GMIDRold(1:14,2,6) = (/ 0.1643, 0.2624, 0.3110, 0.3347, 0.3461, &
                              0.3517, 0.3543, 0.3556, 0.3562, 0.3564, &
                              0.3565, 0.3566, 0.3566, 0.3566 /)
  
  sib%GMIDRold(1:14,1,7) = (/ (1., I=1,14) /)
  sib%GMIDRold(1:14,2,7) = (/ (1., I=1,14) /)
  
  sib%GMIDRold(1:14,1,8) = (/ (1., I=1,14) /)
  sib%GMIDRold(1:14,2,8) = (/ (1., I=1,14) /)
  
  sib%GMIDRold(1:14,1,9) = (/ (1., I=1,14) /)
  sib%GMIDRold(1:14,2,9) = (/ (1., I=1,14) /)


!  IF (LDAS%KOSTER .EQ. 1) THEN

!   do j=1,2
!    do i=1,14
!     sib%GMIDR(i,j,1)=sib%GMIDRold(i,j,1)
!     sib%GMIDR(i,j,2)=sib%GMIDRold(i,j,2)
!     sib%GMIDR(i,j,3)=sib%GMIDRold(i,j,3)
!     sib%GMIDR(i,j,4)=sib%GMIDRold(i,j,4)
!     sib%GMIDR(i,j,5)=sib%GMIDRold(i,j,5)
!     sib%GMIDR(i,j,6)=sib%GMIDRold(i,j,6)
!     sib%GMIDR(i,j,7)=sib%GMIDRold(i,j,7)
!     sib%GMIDR(i,j,8)=sib%GMIDRold(i,j,8)
!     sib%GMIDR(i,j,9)=sib%GMIDRold(i,j,9)
!    enddo
!   enddo

! ELSE
  
  do j=1,2
     do i=1,14
        sib%GMIDR(i,j,1)=sib%GMIDRold(i,j,3)
        
        sib%GMIDR(i,j,2)=sib%GMIDRold(i,j,1)
        
        sib%GMIDR(i,j,3)=sib%GMIDRold(i,j,3)
        
        sib%GMIDR(i,j,4)=sib%GMIDRold(i,j,2)
        
        sib%GMIDR(i,j,5)=(.5*sib%GMIDRold(i,j,2) +    &
                          .5*sib%GMIDRold(i,j,3))
        
        sib%GMIDR(i,j,6)=(.5702*sib%GMIDRold(i,j,3) + &
                          .2200*sib%GMIDRold(i,j,4) + &
                          .2098*sib%GMIDRold(i,j,2))
        
        sib%GMIDR(i,j,7)=(.6372*sib%GMIDRold(i,j,4) + &
                          .1955*sib%GMIDRold(i,j,3) + &
                          .1564*sib%GMIDRold(i,j,2) + &
                          .0109*sib%GMIDRold(i,j,5))
        
        sib%GMIDR(i,j,8)=(.4403*sib%GMIDRold(i,j,8) + &
                          .4365*sib%GMIDRold(i,j,4) + &
                          .0743*sib%GMIDRold(i,j,6) + &
                          .0489*sib%GMIDRold(i,j,5))
        
        sib%GMIDR(i,j,9)=(.8506*sib%GMIDRold(i,j,8) + &
                          .0950*sib%GMIDRold(i,j,5) + &
                          .0399*sib%GMIDRold(i,j,4) + &
                          .0145*sib%GMIDRold(i,j,6))
        
        sib%GMIDR(i,j,10)=sib%GMIDRold(i,j,4)
        
        sib%GMIDR(i,j,11)=sib%GMIDRold(i,j,4)
        
        sib%GMIDR(i,j,12)=sib%GMIDRold(i,j,7)
        
        sib%GMIDR(i,j,13)=(.7114*sib%GMIDRold(i,j,4) + &
                           .1055*sib%GMIDRold(i,j,2) + &
                           .0723*sib%GMIDRold(i,j,3) + &
                           .0526*sib%GMIDRold(i,j,8) + &
                           .0178*sib%GMIDRold(i,j,6) + &
                           .0077*sib%GMIDRold(i,j,5) + &
                           .0327*sib%GMIDRold(i,j,7))
     end do !i
  end do !j

!  ENDIF

  return
end subroutine mapsib2umd

