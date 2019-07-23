!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_gdasLSWG
!  \label{get_gdasLSWG}
!
! !REVISION HISTORY:
!  20 Oct 2009; Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine get_gdasLSWG(n, findex, metdata1, metdata2)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use LIS_metforcingMod, only : LIS_forc
  use LIS_timeMgrMod,     only : LIS_calendar, LIS_get_nstep, LIS_tick
  use LIS_logMod,         only : LIS_logunit, LIS_endrun, LIS_verify
  use gdasLSWG_forcingMod,  only : gdasLSWG_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  real                :: metdata1(LIS_rc%met_nf(findex), LIS_rc%ngrid(n))
  real                :: metdata2(LIS_rc%met_nf(findex), LIS_rc%ngrid(n))
!  
! !DESCRIPTION:
!  This routine issues the calls to open and read the GDAS (LSWG) data. 
!  The data interval is assumed to be 6 hours. 
! 
!EOP
  
  integer             :: yr1,mo1,da1,hr1,mn1,ss1
  integer             :: yr2,mo2,da2,hr2,mn2,ss2
  real*8              :: time1, time2, timenow
  type(ESMF_Time)     :: btime1, btime2
  integer             :: movetime
  integer             :: doy1, doy2
  real                :: gmt1, gmt2,ts1,ts2
  integer             :: status
  integer             :: nstep
  character*100       :: narrfile

  gdasLSWG_struc(n)%findtime1 = 0 
  gdasLSWG_struc(n)%findtime2 = 0 
  movetime = 0 

  nstep=LIS_get_nstep(LIS_rc,n)
  
  if(nstep.eq.0.or.nstep.eq.1.or.LIS_rc%rstflag(n).eq.1) then
     gdasLSWG_struc(n)%findtime1=1
     gdasLSWG_struc(n)%findtime2=1
     movetime=0        ! movetime is not properly set at time-step = 1
     LIS_rc%rstflag(n) = 0
  endif

!-----------------------------------------------------------------
! Determine required NARR data times 
! (previous assimilation, current & future assimilation hours) 
!-----------------------------------------------------------------
! Current time
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0 
  ts1 = 0 
  call LIS_tick(timenow, doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
!previous hour (nearest 6 hour interval)
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = 6*(int(real(LIS_rc%hr)/6.0))
  mn1 = 0 
  ss1 = 0 
  ts1 = 0 

  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  call ESMF_TimeSet(btime1, yy=yr1, &
       mm = mo1, dd = da1, h=hr1, m = 0 , s = 0, calendar=LIS_calendar, &
       rc=status)
  call LIS_verify(status, 'ESMF_TimeSet : get_gdasLSWG')

!Next hour (nearest 6 hour interval)
  yr2 = LIS_rc%yr
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 6*(int(real(LIS_rc%hr)/6.0))
  mn2 = 0 
  ss2 = 0 
  ts2 = 6*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  call ESMF_TimeSet(btime2, yy=yr2, &
       mm = mo2, dd = da2, h=hr2, m = 0 , s = 0, calendar=LIS_calendar, &
       rc=status)
  call LIS_verify(status, 'ESMF_TimeSet : get_gdasLSWG')

  if(timenow.ge.gdasLSWG_struc(n)%st_real.and.&
       gdasLSWG_struc(n)%startRead) then 
     gdasLSWG_struc(n)%time1 = time1
     gdasLSWG_struc(n)%time2 = time2
     gdasLSWG_struc(n)%btime1 = btime1
     gdasLSWG_struc(n)%btime2 = btime2
     movetime = 0 
  endif
  
  if(timenow.ge.gdasLSWG_struc(n)%time2) then 
     movetime = 1
     gdasLSWG_struc(n)%findtime2 = 1
  endif
  

!-----------------------------------------------------------------
! Reading first bookmark
!-----------------------------------------------------------------

  if(gdasLSWG_struc(n)%findtime1.eq.1) then 
    
     call retgdasLSWG(n, 1, findex, metdata1)
     
     gdasLSWG_struc(n)%time1 = time1
     gdasLSWG_struc(n)%btime1 = btime1
  endif

  if(movetime.eq.1) then 
     gdasLSWG_struc(n)%time1 = gdasLSWG_struc(n)%time2
     gdasLSWG_struc(n)%btime1 = gdasLSWG_struc(n)%btime2
     gdasLSWG_struc(n)%findtime2 = 1
     
! need to transfer data...
     metdata1 = metdata2
  endif
  
  if(gdasLSWG_struc(n)%findtime2.eq.1) then 
     call retgdasLSWG(n, 2, findex, metdata2)
     gdasLSWG_struc(n)%time2 = time2
     gdasLSWG_struc(n)%btime2 = btime2
  endif

end subroutine get_gdasLSWG

