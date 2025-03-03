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
! !ROUTINE: LVT_getADCTimeIndex
! \label{LVT_getADCTimeIndex}
!
! !INTERFACE:
subroutine LVT_getADCTimeIndex(tind)
! 
! !USES:
  use LVT_coreMod, only : LVT_rc

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine computes the time index for the seasonal cycle,
!  depending on the current LVT clock. If the seasonal cycle interval is
!  set to 1, then the seasonal cycle will be computed on a monthly basis. 
!  If the seaonal cycle interval is set to 2, then the seasonal cycle will
!  be computed on a 3-month basis (DJF, MAM, JJA and SON). 
!
!  If the stats writing interval in LVT is set to 2592000 (month), then 
!  month switch is detected by change in the month value. So as a result, 
!  the index is defined as the current month subtracted by 1. On the
!  other hand, if the stats writing interval is less than a month, then
!  the index is same as the current month. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  June 1, 2011 : Sujay Kumar, Initial specification
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer    :: interval
  integer    :: tind
!EOP
  
  tind = nint((LVT_rc%hr*3600.0+LVT_rc%mn*60+LVT_rc%ss)/&
       real(LVT_rc%statswriteint))+1

end subroutine LVT_getADCTimeIndex
