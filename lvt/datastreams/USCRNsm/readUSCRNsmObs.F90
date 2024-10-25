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
! !ROUTINE: readUSCRNsmObs
! \label{readUSCRNsmObs}
!
! !INTERFACE: 
subroutine readUSCRNsmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_endrun, LVT_verify
  use USCRNsm_obsMod,      only : uscrnsmobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)       :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for USCRN station data, 
! to process the soil moisture and soil temperature measurements.  
! At the start of the program, an entire year's worth of 
! data is read and stored based on the time location and the station 
! it corresponds to. 
!
! At future times, the read routine simply indexes into the right
! location. Depending upon the frequency of computing output statistics,
! the routine also computes time average (between different LIS 
! output intervals). The routine also interpolates the soil moisture
! and soil temperature data to the model's vertical resolution. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Jan 2017   Sujay Kumar  Initial Specification
!  20 Jul 2018   Mahdi Navari Bug fix: issues with leap years and 
!  	shifting the time between 36 to 12 hours. Please see the ticket #708
!
!EOP

  integer             :: stn_col, stn_row
  real                :: col, row
  character*10        :: wbanno, crx_vn,sur_temp_daily_type
  integer             :: id, yr, mo, da, hr, mn, ss, doy, pentad
  real                :: c1smv, c1tmp, c2smv, c2tmp
  real                :: c3smv, c3tmp, c4smv, c4tmp, c5smv, c5tmp
  integer             :: ftn, i, ios
  integer             :: c,t, st, et, r,kk
  integer             :: status
  logical             :: readflag, file_exists
  character*100       :: uscrnsmfilename 
  real                :: gmt,lat, lon
  real                :: t_daily_max, t_daily_min, t_daily_mean, t_daily_avg
  real                :: p_daily_calc,solrad_daily
  real                :: sur_temp_daily_max,sur_temp_daily_min,sur_temp_daily_avg
  real                :: rh_daily_max,rh_daily_min,rh_daily_avg
  real*8              :: lis_prevtime
  type(ESMF_Time)     :: uscrnsmtime, uscrnsmtime1, uscrnsmtime2, starttime
  type(ESMF_Time)     :: initTime
  type(ESMF_TimeInterval) :: dayInterval
  real                :: sm_tp(5), soilt_tp(5)  
  real                :: depth(5)
  real                :: sf_wt(5), rz_wt(5)
  real                :: sf_wtsum, rz_wtsum
  real                :: sfsm, rzsm, sfst, rzst

  real                :: stc(LVT_rc%lnc, LVT_rc%lnr)
  real                :: smc(LVT_rc%lnc, LVT_rc%lnr)
  real                :: roott(LVT_rc%lnc, LVT_rc%lnr)
  real                :: rootsm(LVT_rc%lnc, LVT_rc%lnr)
  real                :: dummy(LVT_rc%lnc, LVT_rc%lnr)
  real                :: tmp

  smc     = LVT_rc%udef
  stc     = LVT_rc%udef
  roott   = LVT_rc%udef
  rootsm  = LVT_rc%udef

  yr = LVT_rc%dyr(source)
  mo = LVT_rc%dmo(source)
  da = LVT_rc%dda(source)
  hr = LVT_rc%dhr(source)
  mn = LVT_rc%dmn(source)
  ss = LVT_rc%dss(source) 

  depth(1) = 0.05
  depth(2) = 0.10
  depth(3) = 0.20
  depth(4) = 0.50
  depth(5) = 1.00

  if(((uscrnsmobs(source)%yr.ne.LVT_rc%dyr(source))).or.&
       LVT_rc%resetFlag(source))then 
     
     LVT_rc%resetFlag(source) = .false. 

     call ESMF_TimeSet(uscrnsmobs(source)%startTime,  yy=LVT_rc%dyr(source), &
          mm = 1, &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting USCRNsm start time')

     uscrnsmobs(source)%yr = LVT_rc%dyr(source)
     uscrnsmobs(source)%sm = LVT_rc%udef
     uscrnsmobs(source)%rootsm = LVT_rc%udef
     uscrnsmobs(source)%soilt = LVT_rc%udef
     uscrnsmobs(source)%roott = LVT_rc%udef

     ss = 0 
     do i=1,uscrnsmobs(source)%nstns
        call create_USCRNsm_filename(uscrnsmobs(source)%odir, &
             LVT_rc%dyr(source), &
             uscrnsmobs(source)%stnid(i), uscrnsmfilename)
        
        inquire(file=trim(uscrnsmfilename), exist=file_exists) 
        
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading USCRNsm file, ',&
                trim(uscrnsmfilename)
           
           ftn=LVT_getNextUnitNumber()
           open(ftn,file=trim(uscrnsmfilename),form='formatted')
           readflag = .true.         
           do while(readflag) 
              
              read(ftn,fmt='(a5,1X,I4.4,I2.2,I2.2,1X,a6,1X,F7.2,1X,F7.2,5(1x,F7.2),1x,F8.2,1x,a1,16(1x,F7.2))',&
                   iostat=ios) wbanno, yr,mo,da,&
                   crx_vn, lat, lon, &
                   t_daily_max, t_daily_min, t_daily_mean, t_daily_avg,&
                   p_daily_calc, solrad_daily,sur_temp_daily_type, &
                   sur_temp_daily_max,sur_temp_daily_min,sur_temp_daily_avg,&
                   rh_daily_max,rh_daily_min,rh_daily_avg, &
                   c1smv, c2smv, c3smv, c4smv, c5smv, &
                   c1tmp, c2tmp, c3tmp, c4tmp, c5tmp

              if(ios.eq.0) then 

                 sm_tp(1) = c1smv
                 sm_tp(2) = c2smv
                 sm_tp(3) = c3smv
                 sm_tp(4) = c4smv
                 sm_tp(5) = c5smv
                 
                 soilt_tp(1) = c1tmp
                 soilt_tp(2) = c2tmp
                 soilt_tp(3) = c3tmp
                 soilt_tp(4) = c4tmp
                 soilt_tp(5) = c5tmp

                 call ESMF_TimeSet(uscrnsmtime, yy=yr,&
                      mm = mo, dd=da, calendar=LVT_calendar, &
                      rc=status)
                 call LVT_verify(status, 'Error in TimeSet:readUSCRNsmobs')
                 
                 t = nint((uscrnsmtime-uscrnsmobs(source)%starttime)/&
                      uscrnsmobs(source)%timestep) + 1
                          
                 call compute_vinterp_weights(&
                      3,LVT_rc%lis_sf_d,&
                      LVT_rc%lis_rz_d,&
                      depth(1:3),sf_wt,rz_wt)
                 
                 sfsm = 0 
                 rzsm = 0 
                 sf_wtsum = 0 
                 rz_wtsum = 0 
                 do kk=1,3
                    if(sf_wt(kk).ne.0.and.sm_tp(kk).lt.0.01) then 
                       sfsm = LVT_rc%udef
                       exit
                    else
                       sfsm = sfsm + sf_wt(kk)*sm_tp(kk)
                       sf_wtsum = sf_wtsum + sf_wt(kk)
                    endif
                    if(rz_wt(kk).ne.0.and.sm_tp(kk).lt.0.01) then 
                       rzsm = LVT_rc%udef
                       exit
                    else
                       rzsm = rzsm + rz_wt(kk)*sm_tp(kk)
                       rz_wtsum = rz_wtsum + rz_wt(kk)
                    endif
                 enddo
                 if(sfsm.ne.LVT_rc%udef.and.&
                      abs(sf_wtsum-1.0).lt.0.001) then 
                    uscrnsmobs(source)%sm(i,t) = sfsm
                 else
                    uscrnsmobs(source)%sm(i,t) = LVT_rc%udef
                 endif
                 if(rzsm.ne.LVT_rc%udef.and.&
                      abs(rz_wtsum-1.0).lt.0.001) then 
                    uscrnsmobs(source)%rootsm(i,t) = rzsm
                 else
                    uscrnsmobs(source)%rootsm(i,t) = LVT_rc%udef
                 endif

                 call compute_vinterp_weights(&
                      3,LVT_rc%lis_sf_d,&
                      LVT_rc%lis_rz_d,&
                      depth(1:3),sf_wt,rz_wt)
                 
                 sfst = 0 
                 rzst = 0 
                 sf_wtsum = 0 
                 rz_wtsum = 0 

                 do kk=1,3
                    if(sf_wt(kk).ne.0.and.soilt_tp(kk).lt.0.01) then 
                       sfst = LVT_rc%udef
                       exit
                    else
                       sfst = sfst + sf_wt(kk)*soilt_tp(kk)
                       sf_wtsum = sf_wtsum + sf_wt(kk)
                    endif
                    if(rz_wt(kk).ne.0.and.soilt_tp(kk).lt.0.01) then 
                       rzst = LVT_rc%udef
                       exit
                    else
                       rzst = rzst + rz_wt(kk)*soilt_tp(kk)
                       rz_wtsum = rz_wtsum + rz_wt(kk)
                    endif
                 enddo
                 if(sfst.ne.LVT_rc%udef.and.&
                      abs(sf_wtsum-1).lt.0001) then 
                    uscrnsmobs(source)%soilt(i,t) = sfst                    
                 else
                    uscrnsmobs(source)%soilt(i,t) = LVT_rc%udef
                 endif
                 if(rzst.ne.LVT_rc%udef.and.&
                      abs(rz_wtsum-1).lt.0.001) then 
                    uscrnsmobs(source)%roott(i,t) = rzst
                 else
                    uscrnsmobs(source)%roott(i,t) = LVT_rc%udef
                 endif

              endif
              if(ios.ne.0) then 
                 readflag = .false. 
              endif
           enddo

           call LVT_releaseUnitNumber(ftn)

        endif
     enddo   
  endif
                 
  call ESMF_TimeSet(uscrnsmtime1, yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'uscrnsmtime1 set failed')
!MN
!  t = nint((uscrnsmtime1 - uscrnsmobs(source)%starttime)/&
!       uscrnsmobs(source)%timestep)+1

  tmp = nint((uscrnsmtime1 - uscrnsmobs(source)%starttime)/&
       (uscrnsmobs(source)%timestep/24))+1
  t = ceiling(tmp/24)

  do i=1,uscrnsmobs(source)%nstns
     call latlon_to_ij(LVT_domain%lvtproj, &
          uscrnsmobs(source)%stnlat(i), uscrnsmobs(source)%stnlon(i),&
          col, row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
          stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
        if(uscrnsmobs(source)%soilt(i,t).ne.LVT_rc%udef) then 
           stc(stn_col, stn_row) = &
                uscrnsmobs(source)%soilt(i,t)
           if(uscrnsmobs(source)%soilt(i,t).lt.0) then 
              print*, '[ERR] Error in stc ',stn_col, stn_row, i, t, &
                   uscrnsmobs(source)%soilt(i,t)
              stop
           endif
        endif
        
        if(uscrnsmobs(source)%roott(i,t).ne.LVT_rc%udef) then 
           roott(stn_col, stn_row) =&
                uscrnsmobs(source)%roott(i,t)
           if(uscrnsmobs(source)%roott(i,t).lt.0) then 
              print*, '[ERR] Error in roott ',stn_col, stn_row, i, t,&
                   uscrnsmobs(source)%roott(i,t)
              stop
           endif
        endif
        if(uscrnsmobs(source)%sm(i,t).ne.LVT_rc%udef) then 
           if(uscrnsmobs(source)%sm(i,t).lt.0) then 
              print*, '[ERR] Error in sm ',stn_col, stn_row, i, t, &
                   uscrnsmobs(source)%sm(i,t), LVT_rc%udef
              stop
           endif
           smc(stn_col, stn_row) = &
                uscrnsmobs(source)%sm(i,t)
        endif
        if(uscrnsmobs(source)%rootsm(i,t).ne.LVT_rc%udef) then 
           rootsm(stn_col, stn_row) = &
                uscrnsmobs(source)%rootsm(i,t)
           if(uscrnsmobs(source)%rootsm(i,t).lt.0) then 
              print*, '[ERR] Error in rootsm ',stn_col, stn_row, i, t, &
                   uscrnsmobs(source)%rootsm(i,t)
              stop
           endif
        endif
     endif
  enddo
       
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(smc(c,r).gt.0.5) then !ignore
           smc(c,r) = LVT_rc%udef
        endif
        if(rootsm(c,r).gt.0.5) then !ignore
           rootsm(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  dummy = LVT_rc%udef

  call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp, source, stc,&
       vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, smc,&
       vlevel=1,units="m3/m3")
  do c=2,LVT_rc%nstlayers
     call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp, source, &
          dummy,vlevel=c,units="K")
  enddo

  do c=2,LVT_rc%nsmlayers
     call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
          dummy,vlevel=c,units="m3/m3")
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_roottemp, source, &
       roott,vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_rootmoist, source, &
       rootsm,vlevel=1,units="m3/m3")
 
end subroutine readUSCRNsmObs

!BOP
! 
! !ROUTINE: create_USCRNsm_filename
! \label{create_USCRNsm_filename}
!
! !INTERFACE: 
subroutine create_USCRNsm_filename(odir, yr, stnid, uscrnsmname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  character(len=*), intent(in)  :: stnid
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: uscrnsmname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the USCRN station
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      USCRNsm base directory
!   \item[yr]        year of data
!   \item[stnid]     Station ID 
!   \item[uscrnsmname]  Name of the USCRNsm file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4       :: fyr
  
  write(unit=fyr, fmt='(i4.4)') yr
  
  uscrnsmname = trim(odir)//'/'//trim(fyr)//'/CRND0103-'//trim(fyr)//&
       '-'//trim(stnid)//'.txt'

end subroutine create_USCRNsm_filename
