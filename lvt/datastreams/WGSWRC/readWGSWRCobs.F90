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
! !ROUTINE: readWGSWRCObs
! \label{readWGSWRCObs}
!
! !INTERFACE: 
subroutine readWGSWRCObs(source)
! 
! !USES: 
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod,  only : LVT_logSingleVar
  use LVT_obsDataMod,   only : LVT_obsData
  use LVT_timeMgrMod,   only : LVT_tick, LVT_date2time, LVT_calendar, &
       LVT_doy2moda
  use LVT_logMod,       only : LVT_getNextUnitNumber, LVT_releaseUnitNumber, &
       LVT_logunit, LVT_endrun, LVT_verify
  use WGSWRC_obsMod,      only : WGSWRCobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)  :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 May 2011: Sujay Kumar, Initial Specification
! 
!EOP

  character*200         :: filename
  integer               :: ios
  integer               :: yr,doy,mo,da,hr,mn,ss
  real*8                :: lis_prevtime
  real                  :: gmt
  integer               :: st,et
  logical               :: file_exists
  character*20          :: name1, name2
  integer               :: stnid
  integer               :: sm(7)
  real                  :: st1, st2, st3, st4, st5,st6,st7
  character*1           :: st_f, sm_f
  real                  :: col,row
  integer               :: stn_col, stn_row
  integer               :: status
  type(ESMF_Time)       :: obstime,obstime1,obstime2
  real                  :: offset
  integer               :: ftn
  integer               :: k,c,r,kk 
  integer               :: tind
  logical               :: data_flg
  real                  :: smdepth(6),sf_wt(6),rz_wt(6)
  real                  :: sfsm, rzsm
  real                  :: smc(LVT_rc%lnc,LVT_rc%lnr)
  integer               :: nsmc(LVT_rc%lnc,LVT_rc%lnr)
  real                  :: rzsmc(LVT_rc%lnc,LVT_rc%lnr)
  integer               :: nrzsmc(LVT_rc%lnc,LVT_rc%lnr)
  real                  :: stc(LVT_rc%lnc,LVT_rc%lnr)
  integer               :: nstc(LVT_rc%lnc,LVT_rc%lnr)

  smc = -9999.0
  stc = -9999.0
!soil moisture depths in meters
  smdepth(1) = 0.05
  smdepth(2) = 0.10
  smdepth(3) = 0.15
  smdepth(4) = 0.20
  smdepth(5) = 0.30
  smdepth(6) = 0.50
  
  if((LVT_rc%dyr(source).ne.WGSWRCobs(source)%yr).or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 
     WGSWRCobs(source)%yr = LVT_rc%dyr(source)
     call ESMF_TimeSet(WGSWRCobs(source)%startTime,&
          yy=LVT_rc%dyr(source),&
          mm=LVT_rc%dmo(source),&
          dd=1,&
          h=0,&
          m=0,&
          calendar=LVT_calendar,&
          rc=status)

     call LVT_verify(status,'error in setting WGSWRC start time')

     do k = 1, WGSWRCobs(source)%n_stns
        call create_WGSWRCsmobs_filename(WGSWRCobs(source)%odir,&
             WGSWRCobs(source)%stn_name(k), &
             LVT_rc%dyr(source), &
             filename)
        inquire(file=trim(filename), exist=file_exists)
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading ',trim(filename)
           ftn = LVT_getNextUnitNumber()
           open(ftn,file=trim(filename),form='formatted')
           ios = 0 
           do while(ios.eq.0) 
              read(ftn,200,iostat=ios) stnid, yr, doy, hr, mn, &
                   (sm(kk),kk=1,7)
              call LVT_doy2moda(yr,doy,mo,da)

              if(ios.ne.0) exit
              call ESMF_TimeSet(obstime, yy=yr,mm=mo,dd=da,h=hr,m=mn,&
                   calendar=LVT_calendar,rc=status)
              call LVT_verify(status, 'ESMF_TimeSet in readWGSWRCObs')
             
              tind = nint((obstime-WGSWRCobs(source)%startTime)/WGSWRCobs(source)%timestep)+1

              if(tind.gt.0) then 
                 data_flg = .true. 
                 do kk=1,6
                    if(float(sm(kk))/100.0.gt.0.5) then 
                       data_flg = .false. 
                    endif
                 enddo
                 
                 if(data_flg) then 
                    
                    call compute_vinterp_weights(LVT_rc%nsmlayers,&
                         6,LVT_rc%lis_sf_d,&
                         LVT_rc%lis_rz_d,&
                         LVT_rc%smthick(:),&
                         smdepth,sf_wt,rz_wt)
                    sfsm = 0
                    rzsm = 0
                    do kk=1,6
                       sfsm = sfsm +sf_wt(kk)*float(sm(kk))/100.0
                       rzsm = rzsm +rz_wt(kk)*float(sm(kk))/100.0
                    enddo
                    
                    WGSWRCobs(source)%sm(:,:,tind) = sfsm
                    WGSWRCobs(source)%rootsm(:,:,tind) = rzsm
                    
                 else
                    WGSWRCobs(source)%sm(:,:,tind) = LVT_rc%udef
                    WGSWRCobs(source)%rootsm(:,:,tind) = LVT_rc%udef
                 endif
              endif
           enddo
           call LVT_releaseUnitNumber(ftn)
        endif
     enddo
!soil temperature
     do k = 1, WGSWRCobs(source)%n_stns
        call create_WGSWRCstobs_filename(WGSWRCobs(source)%odir,&
             WGSWRCobs(source)%stn_name(k), &
             LVT_rc%dyr(source), &
             filename)
        inquire(file=trim(filename), exist=file_exists)
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading ',trim(filename)
           ftn = LVT_getNextUnitNumber()
           open(ftn,file=trim(filename),form='formatted')
           ios = 0 
           do while(ios.eq.0) 
              read(ftn,201,iostat=ios) stnid, yr, doy, hr, mn, &
                   st1, st2, st3, st4, st5,st6,st7
              call LVT_doy2moda(yr,doy,mo,da)
              if(ios.ne.0) exit
              call ESMF_TimeSet(obstime, yy=yr,mm=mo,dd=da,h=hr,m=mn,&
                   calendar=LVT_calendar,rc=status)
              call LVT_verify(status, 'ESMF_TimeSet in readWGSWRCObs')
              
              tind = nint((obstime-WGSWRCobs(source)%startTime)/WGSWRCobs(source)%timestep)+1
              if(tind.gt.0) then
                 if((st1+273.15).le.350) then 
                    WGSWRCobs(source)%st(:,:,tind) = st1 + 273.15
                 endif
              endif
           enddo
           call LVT_releaseUnitNumber(ftn)
        endif
     enddo
  endif
  
  call ESMF_TimeSet(obstime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'obstime1 set failed')
  
  yr = LVT_rc%dyr(source)
  mo = LVT_rc%dmo(source)
  da = LVT_rc%dda(source)
  hr = LVT_rc%dhr(source)
  mn = LVT_rc%dmn(source)
  ss = LVT_rc%dss(source)
  
  call LVT_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,(-1)*LVT_rc%ts)
  
  call ESMF_TimeSet(obstime2, yy=yr, &
       mm = mo, &
       dd = da, &
       h = hr, &
       m = mn, &
       calendar = LVT_calendar, &
       rc=status)

  call LVT_verify(status, 'error in timeset: readARMobs')

  et = nint((obstime1-WGSWRCobs(source)%starttime)/WGSWRCobs(source)%timestep)+1
  st = nint((obstime2-WGSWRCobs(source)%starttime)/WGSWRCobs(source)%timestep)+1

  smc = 0.0
  nsmc = 0
  do tind=st,et
     if(WGSWRCobs(source)%sm(1,1,tind).ne.-9999.0) then 
        smc(:,:) = smc(:,:)+ WGSWRCobs(source)%sm(:,:,tind)
        nsmc(:,:) = nsmc(:,:) + 1
     endif
  enddo

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nsmc(c,r).ne.0) then 
           smc(c,r) = smc(c,r)/nsmc(c,r)
        else
           smc(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  rzsmc = 0.0
  nrzsmc = 0
  do tind=st,et
     if(WGSWRCobs(source)%sm(1,1,tind).ne.-9999.0) then 
        rzsmc(:,:) = rzsmc(:,:)+ WGSWRCobs(source)%rootsm(:,:,tind)
        nrzsmc(:,:) = nrzsmc(:,:) + 1
     endif
  enddo

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nrzsmc(c,r).ne.0) then 
           rzsmc(c,r) = rzsmc(c,r)/nrzsmc(c,r)
        else
           rzsmc(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  stc = 0.0
  nstc = 0
  do tind=st,et
     if(WGSWRCobs(source)%st(1,1,tind).ne.-9999.0) then 
        stc(:,:) = stc(:,:)+ WGSWRCobs(source)%st(:,:,tind)
        nstc(:,:) = nstc(:,:) + 1
     endif
  enddo

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nstc(c,r).ne.0) then 
           stc(c,r) = stc(c,r)/nstc(c,r)
        else
           stc(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_obsData(i)%soilmoist_obs,source, smc,vlevel=1)
  call LVT_logSingleDataStreamVar(LVT_obsData(i)%soiltemp_obs,source, stc,vlevel=1)
  call LVT_logSingleDataStreamVar(LVT_obsData(i)%rootmoist_obs,source, rzsmc,vlevel=1)
 
200 format (2X,I4.4,1X,I4.4,I5.5,I5.5,I5.5,7I5.5)
201 format (2X,I4.4,1X,I4.4,I5.5,I5.5,I5.5,7F7.1)
end subroutine readWGSWRCObs

!BOP
! 
! !ROUTINE: create_WGSWRCsmobs_filename
! \label(create_WGSWRCsmobs_filename)
!
! !INTERFACE:
subroutine create_WGSWRCsmobs_filename(odir, stn, yr, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in) :: odir
  character(len=*), intent(in) :: stn
  integer  ,        intent(in) :: yr
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  character(len=*)             :: filename

  character*4 :: fyr

  write(unit=fyr, fmt='(i4.4)') yr
  
  filename = trim(odir)//'/'//trim(stn)//'sm_'//trim(fyr)//'.dat'
  
end subroutine create_WGSWRCsmobs_filename

!BOP
! 
! !ROUTINE: create_WGSWRCstobs_filename
! \label(create_WGSWRCstobs_filename)
!
! !INTERFACE:
subroutine create_WGSWRCstobs_filename(odir, stn, yr, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in) :: odir
  character(len=*), intent(in) :: stn
  integer  ,        intent(in) :: yr
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  character(len=*)             :: filename

  character*4 :: fyr

  write(unit=fyr, fmt='(i4.4)') yr
  
  filename = trim(odir)//'/'//trim(stn)//'st_'//trim(fyr)//'.dat'
  
end subroutine create_WGSWRCstobs_filename
