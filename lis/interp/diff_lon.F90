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
! !ROUTINE: diff_lon
! \label{diff_lon}
! 
! !INTERFACE: 
function diff_lon(elon, slon)

! !USES: 
! none

! !ARGUMENTS: 
   implicit none

   real             :: diff_lon
   real, intent(in) :: elon, slon

! !DESCRIPTION: 
! This function returns the difference of two longitudes,
! where the difference is a positive value representing the number
! of degrees eastward from slon to elon.
!
! In general, to compute the number of degrees between
! longitude lon\_a and longitude lon\_b, you would subtract
! them --- dlon = lon\_a - lon\_b
!
! However, longitudes can be specified as [-180, 180]
! or as [0, 360] and issues arise when the starting
! and ending longitudes span 0 or -180.
!
! To avoid these issues, this function adjusts the input longitudes
! to [0, 360] before subtracting.  Then this function assumes that a
! negative result indicates that the difference has wrapped the globe,
! and it adjusts the difference.
!
! For example, say your starting longitude is 1 deg, and your ending
! longitude is -1 deg.  The difference should be ending - starting,
! which gives -1 - 1 = -2 degrees; where it is actually 358 degrees.
!
!
! USAGE: if you wish to compute:
!
!    dlon = lon\_a - lon\_b
!
! then use this function thusly:
!
!    dlon = diff\_lon(lon\_a, lon\_b)
!
! I.e.; the order of the arguments to diff\_lon are in the same
! order that they would be in if you simply performed the subtraction.
!
! !REVISION HISTORY:
!  05 Nov 2014: James Geiger; Initial specification
! 
!EOP      
   
   real :: slon360, elon360

   if ( slon == elon ) then
      diff_lon = 0.0
   else
      ! Adjust longitudes to [0, 360].
      if ( slon < 0 ) then
         slon360 = slon + 360.0
      else
         slon360 = slon
      endif
      if ( elon < 0 ) then
         elon360 = elon + 360.0
      else
         elon360 = elon
      endif

      diff_lon = elon360 - slon360

      if ( diff_lon <= 0 ) then
         diff_lon = diff_lon + 360.0
      endif
   endif

end function diff_lon
