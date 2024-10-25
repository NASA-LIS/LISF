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
! !ROUTINE: readGRDCObs
! \label{readGRDCObs}
!
! !INTERFACE: 
  subroutine readGRDCObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,        only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,     only : LVT_calendar, LVT_tick
  use LVT_logMod,         only : LVT_verify, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_logunit
  use GRDC_obsMod, only : GRDCobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!  The Global Runoff Data Centre (GRDC), a repository for 
!  the world's river discharge data. The Global Runoff Database 
!  is a unique collection of river discharge data collected at 
!  daily or monthly intervals from more than 9,300 stations in 
!  160 countries. This adds up to around 400,000 station-years 
!  with an average record length of 43 years. The GRDC provides 
!  discharge data and data products for non-commercial applications.
!
!  The GRDC operates under the auspices of the World Meteorological 
!  Organisation (WMO), and the German Federal Institute of Hydrology 
!  (Bundesanstalt für Gewässerkunde or BfG) hosts the GRDC in Koblenz. 
!
!  http://www.bafg.de/GRDC/EN/Home/homepage_node.html
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  11 May 2011: Sujay Kumar, Initial Specification
! 
!EOP
!--------------------------------------------------------------------------



  character*100       :: filename
  character*20        :: datestring, qstring
  integer             :: stn_col, stn_row, c,r
  real                :: col,row
  integer             :: tind
  real                :: offset
  real                :: q_orig, q_calc, flag
  real                :: discharge
  integer             :: i,t,kk
  integer             :: ftn
  integer             :: iloc
  integer             :: ios,ios1,ios2,status
  integer             :: yr, mo, da, cmo
  character*100       :: line
  logical             :: file_exists
  type(ESMF_Time)     :: obsTime,obsTime1
  real                :: q(LVT_rc%lnc,LVT_rc%lnr)

  q = LVT_rc%udef

  ! Every new year, read the data for each station and store it in memory:

  if((LVT_rc%dyr(source).ne.GRDCobs(source)%yr).or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 
     GRDCobs(source)%yr = LVT_rc%dyr(source)
     GRDCobs(source)%q = LVT_rc%udef

     call ESMF_TimeSet(GRDCobs(source)%startTime, yy=LVT_rc%dyr(source), &
          mm= LVT_rc%dmo(source), &
          dd = 1, &
          h = 0, &
          m = 0, & 
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, "error in setting GRDC start time")

     do i = 1, GRDCobs(source)%n_stns
        call create_GRDCobs_filename(GRDCobs(source)%odir, &
             GRDCobs(source)%use_daily, &
             GRDCobs(source)%stn_id(i), &
             LVT_rc%dyr(source),filename)

        inquire(file=trim(filename), exist=file_exists) 
        if(file_exists) then            
           write(LVT_logunit,*) '[INFO] reading GRDC streamflow file: ',&
                trim(filename)
           ftn = LVT_getNextUnitNumber()
           open(ftn,file=trim(filename),form='formatted')
           ios = 0 
           do while(ios.eq.0) 
              read(ftn,'(a)',iostat=ios) line
              if(ios.ne.0) exit

              iloc = index(line,"-")
              read(line(1:iloc-1),*) yr
              line = line(iloc+1:len(line))
              
              iloc = index(line,"-")
              read(line(1:iloc-1),*) mo
              line = line(iloc+1:len(line))

              iloc = index(line,";")
              read(line(1:iloc-1),*) da
              line = line(iloc+1:len(line))

              iloc = index(line,";")
              line = line(iloc+1:len(line))

              if( GRDCobs(source)%version == 1 ) then ! Prior to 2017
              ! Previous line format (example):
              !  1996-01-01;--:--;     57.837;     57.837; -999

                iloc = index(line,";")
                read(line(1:iloc-1),*) q_orig
                line = line(iloc+1:len(line))

                iloc = index(line,";")
                read(line(1:iloc-1),*) q_calc
                line = line(iloc+1:len(line))

                read(line,*) flag

              elseif( GRDCobs(source)%version == 2 ) then ! For 2017 and onward
              ! New line format:
              !  1981-01-01;--:--;  57490.000
              ! (-999 undefined flag now incorporated in above column)

                read(line,*) q_calc   ! Technically, this is actually q_orig ...

              endif

              if(da.eq.0) da = 1
              call ESMF_TimeSet(obsTime, yy=yr, &
                   mm=mo, dd=da, h=0, calendar=LVT_calendar, &
                   rc=status)
              call LVT_verify(status, 'ESMF_TimeSet in readGRDCobs')
              
              ! Read daily obs
              if(GRDCobs(source)%use_daily.eq.1) then
                 t = nint((obstime-GRDCobs(source)%starttime)/&
                      GRDCobs(source)%timestep)+1
                 
                 if(t.ge.1.and.t.le.366.and. &
                      flag.ne.99) then 
!                    if(q_calc.gt.0) then   ! Original
                    if(q_calc.ge.0) then 
                       GRDCobs(source)%q(i,t) = q_calc
                    endif
                 endif
              ! Read monthly GRDC obs
              else
                 t = mo
                 if(t.ge.1.and.t.le.12.and. &
                      flag.ne.99) then 
!                    if(q_calc.gt.0) then   ! Original
                    if(q_calc.ge.0) then 
                       GRDCobs(source)%q(i,t) = q_calc
                    endif
                 endif
              endif
           end do
           call LVT_releaseUnitNumber(ftn)
        else
           write(LVT_logunit,*) '[WARN] unable to read GRDC streamflow file: ',&
                trim(filename)
        end if
     enddo
  endif

  ! Daily matching of the station with LIS
  if(GRDCobs(source)%use_daily.eq.1) then
     call ESMF_TimeSet(obstime1,yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
          h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'obstime1 set failed')
     
     offset = (obstime1-GRDCobs(source)%starttime)/GRDCobs(source)%timestep
     
     if((nint(offset)-offset).eq.0) then !only when LIS time matches the observation
        
        tind = nint((obstime1-GRDCobs(source)%starttime)/&
             GRDCobs(source)%timestep)+1
        
        ! Locate closest station lat/lon to model gridcell
        do i=1,GRDCobs(source)%n_stns
           call latlon_to_ij(LVT_domain%lvtproj, &
                GRDCobs(source)%stnlat(i),GRDCobs(source)%stnlon(i),&
                col,row)
           stn_col = nint(col)
           stn_row = nint(row)
           
           if(stn_col.gt.0.and.stn_row.gt.0.and.tind.ge.0.and.&
                stn_col.le.LVT_rc%lnc.and.stn_row.le.LVT_rc%lnr) then 
              if(GRDCobs(source)%q(i,tind).ne.LVT_rc%udef) then 
                 q(stn_col, stn_row) = &
                      GRDCobs(source)%q(i,tind)
              endif
           endif
        enddo
     endif
     
  ! Monthly
  else
     if(LVT_rc%d_nmo(source).ne.GRDCobs(source)%mo) then 
        cmo = GRDCobs(source)%mo
        tind = cmo

        ! Locate closest station lat/lon to model gridcell
        do i=1,GRDCobs(source)%n_stns
           call latlon_to_ij(LVT_domain%lvtproj, &
                GRDCobs(source)%stnlat(i),GRDCobs(source)%stnlon(i),&
                col,row)
           stn_col = nint(col)
           stn_row = nint(row)
           
           if(stn_col.gt.0.and.stn_row.gt.0.and.tind.ge.0.and.&
                stn_col.le.LVT_rc%lnc.and.stn_row.le.LVT_rc%lnr) then 
              if(GRDCobs(source)%q(i,tind).ne.LVT_rc%udef) then 
                 q(stn_col, stn_row) = &
                      GRDCobs(source)%q(i,tind)
              endif
           endif
        enddo
        GRDCobs(Source)%mo = LVT_rc%d_nmo(source)
     endif
  endif
  call LVT_logSingleDataStreamVar(LVT_MOC_streamflow,source, q,vlevel=1,units="m3/s")
     
end subroutine readGRDCObs

!BOP
! 
! !ROUTINE: create_GRDCobs_filename
! \label(create_GRDCobs_filename)
!
! !INTERFACE:
subroutine create_GRDCobs_filename(odir, use_daily, stn_id,yr,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*),     intent(in) :: odir
  integer,              intent(in) :: use_daily
  character(len=*),     intent(in) :: stn_id
  integer,              intent(in) :: yr
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
  
  character(len=*)  :: filename
  character*4       :: fyr

  write(unit=fyr,fmt='(i4.4)') yr
  
  if(use_daily.eq.1) then 
     filename = trim(odir)//'/daily/'//trim(stn_id)//&
          '_'//trim(fyr)//'.dat'
  else
     filename = trim(odir)//'/monthly/'//trim(stn_id)//&
          '_'//trim(fyr)//'.dat'
  endif

end subroutine create_GRDCobs_filename
