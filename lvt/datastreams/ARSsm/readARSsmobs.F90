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
! !ROUTINE: readARSsmObs
! \label{readARSsmObs}
!
! !INTERFACE: 
subroutine readARSsmObs(source)
! 
! !USES: 
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_tick, LVT_date2time, LVT_calendar, &
       LVT_doy2moda
  use LVT_logMod,       only : LVT_getNextUnitNumber, LVT_releaseUnitNumber, &
       LVT_logunit, LVT_endrun, LVT_verify
  use ARSsm_obsMod,      only : ARSsmobs

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
!  29 Mar 2012   Sujay Kumar  Initial Specification
! 
!EOP

  character*200         :: filename
  integer               :: ios
  integer               :: yr,doy,mo,da,hr,mn,ss
  real                  :: gmt
  integer               :: st,et
  logical               :: file_exists
  character*20          :: name1, name2
  integer               :: stnid
  character*1           :: st_f, sm_f
  integer               :: stn_col, stn_row
  integer               :: status
  type(ESMF_Time)       :: obstime,obstime1
  real                  :: offset
  integer               :: ftn
  integer               :: k,c,r,kk 
  integer               :: tind
  logical               :: data_flg
  real                  :: sfsm, sfsm_std
  real                  :: smc(LVT_rc%lnc,LVT_rc%lnr)
  real                  :: smstd(LVT_rc%lnc,LVT_rc%lnr)
  real                  :: timenow

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)

  smc = LVT_rc%udef
  smstd = LVT_rc%udef

  if(LVT_rc%dyr(source).ne.ARSsmobs(source)%yr.or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     ARSsmobs(source)%sm = -9999.0
     ARSsmobs(source)%sm_std = -9999.0

     ARSsmobs(source)%yr = LVT_rc%dyr(source)
     call ESMF_TimeSet(ARSsmobs(source)%startTime,&
          yy=LVT_rc%dyr(source),&
          mm=LVT_rc%dmo(source),&
          dd=1,&
          h=0,&
          m=0,&
          calendar=LVT_calendar,&
          rc=status)
     call LVT_verify(status,'error in setting ARSsm start time')

     do k = 1, ARSsmobs(source)%n_stns

        call create_ARSsmobs_filename(ARSsmobs(source)%odir,&
             ARSsmobs(source)%stn_name(k), &
             LVT_rc%dyr(source), &
             filename)
        
        inquire(file=trim(filename), exist=file_exists)
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading ',trim(filename)
           ftn = LVT_getNextUnitNumber()
           open(ftn,file=trim(filename),form='formatted')
           ios = 0 
           do while(ios.eq.0) 
              read(ftn,*,iostat=ios) yr, mo,da, hr, mn, sfsm, sfsm_std
!              call LVT_doy2moda(yr,doy,mo,da)

              if(ios.ne.0) exit
              call ESMF_TimeSet(obstime, yy=yr,mm=mo,dd=da,h=hr,m=mn,&
                   calendar=LVT_calendar,rc=status)
              call LVT_verify(status, 'ESMF_TimeSet in readARSsmObs')
             
              tind = nint((obstime-ARSsmobs(source)%startTime)/&
                   ARSsmobs(source)%timestep)+1

              if(tind.gt.0.and.sfsm.gt.0.and.sfsm.lt.0.5) then 
                 ARSsmobs(source)%sm(k,tind) = sfsm
                 if(sfsm_std.gt.0) then 
                    ARSsmobs(source)%sm_std(k,tind) = sfsm_std
                 else
                    ARSsmobs(source)%sm_std(k,tind) = LVT_rc%udef
                 endif
              endif
           enddo
           call LVT_releaseUnitNumber(ftn)
           write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)
        endif
     enddo
  end if

  if(mod(timenow, 1800.0).eq.0.0) then 
     call ESMF_TimeSet(obstime1,yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
          h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'obstime1 set failed')
     yr = LVT_rc%dyr(source)
     mo = LVT_rc%dmo(source)
     da = LVT_rc%dda(source)
     hr = LVT_rc%dhr(source)
     mn = LVT_rc%dmn(source)
     ss = LVT_rc%dss(source)
  
     tind = nint((obstime1-ARSsmobs(source)%starttime)/&
          ARSsmobs(source)%timestep)+1

     do r=1, LVT_rc%lnr
        do c=1,LVT_rc%lnc
           do k=1,ARSsmobs(source)%n_stns             
              if(c.eq.ARSsmobs(source)%stn_col(k).and.&
                   r.eq.ARSsmobs(source)%stn_row(k).and.&
                   ARSsmobs(source)%sm(k,tind).ne.-9999.0) then 
                 smc(c,r) = ARSsmobs(source)%sm(k,tind)
                 if(c.eq.ARSsmobs(source)%stn_col(k).and.&
                      r.eq.ARSsmobs(source)%stn_row(k).and.&
                      ARSsmobs(source)%sm_std(k,tind).ne.-9999.0) then 
                    smstd(c,r) = ARSsmobs(source)%sm_std(k,tind)
                 endif
              endif
           enddo
        enddo
     enddo
  else
     smc = LVT_rc%udef
     smstd = LVT_rc%udef
  endif

!  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist,source,smc,vlevel=1, &
!       units="m3/m3",stdev=smstd)

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist,source,smc,vlevel=1, &
       units="m3/m3")
 
end subroutine readARSsmObs

!BOP
! 
! !ROUTINE: create_ARSsmobs_filename
! \label(create_ARSsmobs_filename)
!
! !INTERFACE:
subroutine create_ARSsmobs_filename(odir, stn, yr, filename)
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
  
  filename = trim(odir)//'/'//trim(stn)//'_'//trim(fyr)//'.txt'
  
end subroutine create_ARSsmobs_filename

