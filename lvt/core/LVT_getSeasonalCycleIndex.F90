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
! !ROUTINE: LVT_getSeasonalCycleTimeIndex
! \label{LVT_getSeasonalCycleTimeIndex}
!
! !INTERFACE:
subroutine LVT_getSeasonalCycleTimeIndex(interval,tind)
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
  if(interval.eq.1) then 
     if(LVT_rc%statswriteint.eq.2592000) then 
        tind = LVT_rc%mo-1
        if(tind.eq.0) then 
           tind = 12
        endif
     elseif(LVT_rc%statswriteint.lt.2592000) then 
        tind = LVT_rc%mo
     endif
  elseif(interval.eq.2) then 
     if(LVT_rc%statswriteint.eq.2592000) then 
        if(LVT_rc%mo.eq.1.or.LVT_rc%mo.eq.2.or.LVT_rc%mo.eq.3) then
           !DJF
           tind = 1
        elseif(LVT_rc%mo.eq.4.or.LVT_rc%mo.eq.5.or.LVT_rc%mo.eq.6) then
           !MAM
           tind = 2
        elseif(LVT_rc%mo.eq.7.or.LVT_rc%mo.eq.8.or.LVT_rc%mo.eq.9) then
           !JJA
           tind = 3
        elseif(LVT_rc%mo.eq.10.or.LVT_rc%mo.eq.11.or.LVT_rc%mo.eq.12) then
           !SON
           tind = 4
        endif
     elseif(LVT_rc%statswriteint.lt.2592000) then 
        if(LVT_rc%mo.eq.12.or.LVT_rc%mo.eq.1.or.LVT_rc%mo.eq.2) then
           !DJF
           tind = 1
        elseif(LVT_rc%mo.eq.3.or.LVT_rc%mo.eq.4.or.LVT_rc%mo.eq.5) then
           !MAM
           tind = 2
        elseif(LVT_rc%mo.eq.6.or.LVT_rc%mo.eq.7.or.LVT_rc%mo.eq.8) then
           !JJA
           tind = 3
        elseif(LVT_rc%mo.eq.9.or.LVT_rc%mo.eq.10.or.LVT_rc%mo.eq.11) then
           !SON
           tind = 4
        endif

     endif
  elseif(interval.eq.21) then
     if(LVT_rc%statswriteint.eq.2592000) then
        if(LVT_rc%mo.eq.2.or.LVT_rc%mo.eq.3.or.LVT_rc%mo.eq.4) then
           !JFM
           tind = 1
        elseif(LVT_rc%mo.eq.5.or.LVT_rc%mo.eq.6.or.LVT_rc%mo.eq.7) then
           !AMJ
           tind = 2
        elseif(LVT_rc%mo.eq.8.or.LVT_rc%mo.eq.9.or.LVT_rc%mo.eq.10) then
           !JAS
           tind = 3
        elseif(LVT_rc%mo.eq.11.or.LVT_rc%mo.eq.12.or.LVT_rc%mo.eq.1) then
           !OND
           tind = 4
        endif
     elseif(LVT_rc%statswriteint.lt.2592000) then
        if(LVT_rc%mo.eq.1.or.LVT_rc%mo.eq.2.or.LVT_rc%mo.eq.3) then
           !JFM
           tind = 1
        elseif(LVT_rc%mo.eq.4.or.LVT_rc%mo.eq.5.or.LVT_rc%mo.eq.6) then
           !AMJ
           tind = 2
        elseif(LVT_rc%mo.eq.7.or.LVT_rc%mo.eq.8.or.LVT_rc%mo.eq.9) then
           !JAS
           tind = 3
        elseif(LVT_rc%mo.eq.10.or.LVT_rc%mo.eq.11.or.LVT_rc%mo.eq.12) then
           !OND
           tind = 4
        endif
     endif
  endif
end subroutine LVT_getSeasonalCycleTimeIndex
