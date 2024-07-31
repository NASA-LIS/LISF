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
! !ROUTINE: readSCANObs
! \label{readSCANObs}
!
! !INTERFACE: 
subroutine readSCANObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_endrun, LVT_verify
  use SCAN_obsMod,      only : scanobs
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
! This subroutine provides the data reader for SCAN station data, 
! to process the soil moisture and soil temperature measurements.  
! At the start of the program, an entire year's worth of 
! data is read and stored
! based on the time location and the station it corresponds to. 
! At future times, the read routine simply indexes into the right
! location. Depending upon the frequency of computing output statistics,
! the routine also computes time average (between different LIS 
! output intervals). The routine also interpolates the soil moisture
! and soil temperature data to the model's vertical resolution. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP

  integer             :: stn_col, stn_row
  real                :: col, row
  integer             :: id, yr, mo, da, hr, mn, ss, doy, pentad
  real                :: tair
  integer             :: sc_yr,sc_mo,sc_da,sc_hr,sc_mn,sc_ss
  real                :: swd, t2m, c1smv, c1tmp, c2smv, c2tmp
  real                :: c3smv, c3tmp, c4smv, c4tmp, c5smv, c5tmp
  integer             :: ftn, i, ios
  integer             :: c,t, st, et, r,kk
  integer             :: status
  logical             :: readflag, file_exists
  character*100       :: scanfilename 
  real                :: gmt
  real*8              :: lis_prevtime
  type(ESMF_Time)     :: scantime, scantime1, scantime2, starttime
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

  depth(1) = 2*0.0254
  depth(2) = 4*0.0254
  depth(3) = 8*0.0254
  depth(4) = 20*0.0254
  depth(5) = 40*0.0254

  if((scanobs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(scanobs(source)%startTime,  yy=LVT_rc%dyr(source), &
          mm = 1, &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting scan start time')

     scanobs(source)%yr = LVT_rc%dyr(source)
     scanobs(source)%sm = LVT_rc%udef
     scanobs(source)%rootsm = LVT_rc%udef
     scanobs(source)%soilt = LVT_rc%udef
     scanobs(source)%roott = LVT_rc%udef

     ss = 0 
     do i=1,scanobs(source)%nstns
        call create_SCAN_filename(scanobs(source)%odir, &
             LVT_rc%dyr(source), &
             scanobs(source)%stnid(i), scanfilename)
        
        inquire(file=trim(scanfilename), exist=file_exists) 
        
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading SCAN file ',&
                trim(scanfilename)
           
           ftn=LVT_getNextUnitNumber()
           open(ftn,file=trim(scanfilename),form='formatted')
           readflag = .true.         
           read(ftn,*) 
           read(ftn,*) 
           read(ftn,*) 
           read(ftn,*) 

           do while(readflag) 
              if(scanobs(source)%stnid(i).lt.1000) then
                 read(ftn,fmt='(I3.3,1X,I4.4,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,1X,10(F8.3,1X))',&
                      iostat=ios) id, yr,mo,da,hr,mn,&
                      c1smv, c2smv, c3smv, c4smv, c5smv, &
                      c1tmp, c2tmp, c3tmp, c4tmp, c5tmp 
              else                
                 read(ftn,fmt='(I4.4,1X,I4.4,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,1X,10(F8.3,1X))',&
                      iostat=ios) id, yr,mo,da,hr,mn,&
                      c1smv, c2smv, c3smv, c4smv, c5smv, &
                      c1tmp, c2tmp, c3tmp, c4tmp, c5tmp
              endif
              if(ios.eq.0.and.id.ne.0) then 

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

                 call ESMF_TimeSet(scantime, yy=yr,&
                      mm = mo, dd=da, h=hr, m =mn, s=ss,&
                      calendar=LVT_calendar, &
                      rc=status)
                 call LVT_verify(status, 'Error in TimeSet:readSCANobs')
                 
                 t = nint((scantime-scanobs(source)%starttime)/&
                      scanobs(source)%timestep) + 1

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
                       sfsm = sfsm + sf_wt(kk)*sm_tp(kk)/100.0
                       sf_wtsum = sf_wtsum + sf_wt(kk)
                    endif
                    if(rz_wt(kk).ne.0.and.sm_tp(kk).lt.0.01) then 
                       rzsm = LVT_rc%udef
                       exit
                    else
                       rzsm = rzsm + rz_wt(kk)*sm_tp(kk)/100.0
                       rz_wtsum = rz_wtsum + rz_wt(kk)
                    endif
                 enddo
                 if(t.gt.0.and.t.le.8784) then 
                    if(sfsm.ne.LVT_rc%udef.and.&
                         abs(sf_wtsum-1.0).lt.0.001) then 
                       scanobs(source)%sm(i,t) = sfsm
                    else
                       scanobs(source)%sm(i,t) = LVT_rc%udef
                    endif
                    if(rzsm.ne.LVT_rc%udef.and.&
                         abs(rz_wtsum-1.0).lt.0.001) then 
                       scanobs(source)%rootsm(i,t) = rzsm
                    else
                       scanobs(source)%rootsm(i,t) = LVT_rc%udef
                    endif
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
                 if(t.gt.0.and.t.le.8784) then 
                    if(sfst.ne.LVT_rc%udef.and.&
                         abs(sf_wtsum-1).lt.0001) then 
                       scanobs(source)%soilt(i,t) = sfst                    
                    else
                       scanobs(source)%soilt(i,t) = LVT_rc%udef
                    endif
                    if(rzst.ne.LVT_rc%udef.and.&
                         abs(rz_wtsum-1).lt.0.001) then 
                       scanobs(source)%roott(i,t) = rzst
                    else
                       scanobs(source)%roott(i,t) = LVT_rc%udef
                    endif
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
                 
  call ESMF_TimeSet(scantime1, yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'scantime1 set failed')

  t = nint((scantime1 - scanobs(source)%starttime)/&
       scanobs(source)%timestep)+1
  do i=1,scanobs(source)%nstns
     call latlon_to_ij(LVT_domain%lvtproj, &
          scanobs(source)%stnlat(i), scanobs(source)%stnlon(i),&
          col, row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
          stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
        if(scanobs(source)%soilt(i,t).ne.LVT_rc%udef) then 
           stc(stn_col, stn_row) = &
                scanobs(source)%soilt(i,t)
           if(scanobs(source)%soilt(i,t).lt.0) then 
              print*, 'Error in stc ',stn_col, stn_row, i, t, &
                   scanobs(source)%soilt(i,t)
!              stop
           endif
        endif
        
        if(scanobs(source)%roott(i,t).ne.LVT_rc%udef) then 
           roott(stn_col, stn_row) =&
                scanobs(source)%roott(i,t)
           if(scanobs(source)%roott(i,t).lt.0) then 
              print*, 'Error in roott ',stn_col, stn_row, i, t,&
                   scanobs(source)%roott(i,t)
!              stop
           endif
        endif
        if(scanobs(source)%sm(i,t).ne.LVT_rc%udef) then 
           if(scanobs(source)%sm(i,t).lt.0) then 
              print*, 'Error in sm ',stn_col, stn_row, i, t, &
                   scanobs(source)%sm(i,t), LVT_rc%udef
              stop
           endif
           smc(stn_col, stn_row) = &
                scanobs(source)%sm(i,t)
        endif
        if(scanobs(source)%rootsm(i,t).ne.LVT_rc%udef) then 
           rootsm(stn_col, stn_row) = &
                scanobs(source)%rootsm(i,t)
           if(scanobs(source)%rootsm(i,t).lt.0) then 
              print*, 'Error in rootsm ',stn_col, stn_row, i, t, &
                   scanobs(source)%rootsm(i,t)
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
 
end subroutine readSCANObs

!BOP
! 
! !ROUTINE: create_SCAN_filename
! \label{create_SCAN_filename}
!
! !INTERFACE: 
subroutine create_SCAN_filename(odir, yr, stnid, scanname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: stnid
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: scanname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the SCAN station
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      SCAN base directory
!   \item[yr]        year of data
!   \item[stnid]     Station ID 
!   \item[scanname]  Name of the SCAN file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4       :: fyr
  character*4       :: fstnid
  
  write(unit=fyr, fmt='(i4.4)') yr
  if(stnid.ge.1000) then 
     write(unit=fstnid, fmt='(i4)') stnid
  elseif(stnid.ge.100) then 
     write(unit=fstnid, fmt='(i3)') stnid
  else
     write(unit=fstnid, fmt='(i2)') stnid
  endif
  
  scanname = trim(odir)//'/'//trim(fyr)//'/'//trim(fstnid)//&
       '-'//trim(fyr)//'CY.csv'
  
  
end subroutine create_SCAN_filename
