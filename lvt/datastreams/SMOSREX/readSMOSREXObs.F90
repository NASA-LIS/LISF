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
! !ROUTINE: readSMOSREXObs
! \label{readSMOSREXObs}
!
! !INTERFACE: 
subroutine readSMOSREXObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_endrun, LVT_verify
  use SMOSREX_obsMod,      only : smosrexobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,     intent(in)  :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP

  real                :: smc(LVT_rc%lnc, LVT_rc%lnr)
  real                :: nsmc(LVT_rc%lnc, LVT_rc%lnr)
  integer             :: t, st, et, c,r
  real                :: col, row
  integer             :: status
  integer             :: stn_col, stn_row
  integer             :: yr, mo, da, hr, mn, ss
  real                :: lat, lon
  real*8              :: lis_prevtime
  real                :: gmt
  integer             :: doy
  type(ESMF_Time)     :: scantime1, scantime2
  
  smc = 0 
  nsmc = 0 
  
  call ESMF_TimeSet(scantime1, yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'scantime1 set failed')
  
  et = nint((scantime1 - smosrexobs(source)%starttime)/smosrexobs(source)%timestep)+1
  
  yr = LVT_rc%dyr(source)
  mo = LVT_rc%dmo(source)
  da = LVT_rc%dda(source)
  hr = LVT_rc%dhr(source)
  mn = LVT_rc%dmn(source)
  ss = LVT_rc%dss(source)
  
  call LVT_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,(-1)*LVT_rc%ts)

  call ESMF_TimeSet(scantime2, yy=yr, &
       mm = mo, &
       dd = da, &
       h = hr, &
       m = mn, &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status)

  st = nint((scantime2 - smosrexobs(source)%starttime)/smosrexobs(source)%timestep) + 1

  do t=st+1,et !to avoid double counting from last cycle
     lat = 43.36
     lon = 1.3015
     call latlon_to_ij(LVT_domain%lvtproj, lat, lon,&
          col, row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if(t.le.87648) then 
        if(smosrexobs(source)%sm(t).gt.0) then 
           smc(stn_col, stn_row) = smc(stn_col,stn_row) +  smosrexobs(source)%sm(t)
           nsmc(stn_col,stn_row) = nsmc(stn_col,stn_row) + 1          
           if(smosrexobs(source)%sm(t).lt.0) then 
              print*, 'Error in smc ',stn_col, stn_row,  t, smosrexobs(source)%sm(t)
              stop
           endif
        endif
     else
        smc(stn_col, stn_row) = LVT_rc%udef
        nsmc(stn_col, stn_row) = 0 
     endif
  enddo
 
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nsmc(c,r).gt.0) then 
           smc(c,r) = smc(c,r)/nsmc(c,r)
        else
           smc(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source, smc,vlevel=1,units="m3/m3")

end subroutine readSMOSREXObs
