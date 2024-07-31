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
! !ROUTINE: readSURFRADObs
! \label{readSURFRADObs}
!
! !INTERFACE:
subroutine readSURFRADObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod, only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use LVT_logMod,       only : LVT_verify
  use SURFRAD_obsMod, only : SURFRADobs
  use map_utils
  
  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! This subroutine provides the data reader for the SURFRAD station data.
! This routine is called by LVT at each timestep and returns the
! observations corresponding to the current time. The SURFRAD data 
! which is provided at either a 1 minute interval or at 3 minute
! interval is temporally aggregated to the LIS output frequency 
! before mapping to the observational data structures of LVT. 
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! !ARGUMENTS:


    integer :: station_index
    integer, parameter :: num_stations = 7
    
!
!EOP

! read the data and then use the log-method (see below) to map the relevant variable(s)
! to the LVT data structures.
!
    real  		:: time
    integer 		:: yr, mo, da, hr, mn, ss, doy
    real     		:: gmt
    real*8   		:: lis_prevtime
    type(ESMF_Time)     :: initTime
    type(ESMF_Time)     :: reftime
    type(ESMF_TimeInterval) :: dayInterval
    real     		:: col, row
    integer 		:: c,r,stn_col, stn_row, t, st, et
    real          	:: swd(LVT_rc%lnc, LVT_rc%lnr)
    integer          	:: nswd(LVT_rc%lnc, LVT_rc%lnr)
    real          	:: lwd(LVT_rc%lnc, LVT_rc%lnr)
    integer          	:: nlwd(LVT_rc%lnc, LVT_rc%lnr)
    real          	:: wspd(LVT_rc%lnc, LVT_rc%lnr)
    integer          	:: nwspd(LVT_rc%lnc, LVT_rc%lnr)
    real          	:: psurf(LVT_rc%lnc, LVT_rc%lnr)
    integer          	:: npsurf(LVT_rc%lnc, LVT_rc%lnr)
    integer          	:: rc
    integer 		:: i,j
    
    type(ESMF_Time) :: surftime1, surftime2

    swd = LVT_rc%udef
    lwd =  LVT_rc%udef
    wspd =  LVT_rc%udef
    psurf =  LVT_rc%udef
    
    if(LVT_rc%dyr(source).ge.2009.and.LVT_rc%dmo(source).ge.1.and.&
         LVT_rc%dda(source).ge.1) then 
       call ESMF_TimeIntervalSet(SURFRADobs(source)%timestep, s=1*60, rc=rc)
       call LVT_verify(rc, 'error in setting timestep (SURFRADobs)')
    endif

    time = LVT_rc%dhr(source)*3600+LVT_rc%dmn(source)*60+LVT_rc%dss(source)
    if((mod(time,86400.0).eq.0.0) .or. (SURFRADobs(source)%day .ne. &
         LVT_rc%dda(source)).or.&
         LVT_rc%resetFlag(source)) then

       LVT_rc%resetFlag(source) = .false. 
       call ESMF_TimeSet(SURFRADobs(source)%starttime, yy=LVT_rc%dyr(source), &
            mm = LVT_rc%dmo(source), &
            dd = LVT_rc%dda(source), &
            h = 0, &
            m = 0, &
            calendar = LVT_calendar, &
            rc=rc)
       SURFRADobs(source)%day = LVT_rc%dda(source)
       reftime = SURFRADobs(source)%starttime

       call read_station_data(source,1,'bon',LVT_rc%dyr(source),LVT_rc%dmo(source), &
            LVT_rc%dda(source), reftime)
       call read_station_data(source,2,'tbl',LVT_rc%dyr(source),LVT_rc%dmo(source), &
            LVT_rc%dda(source), reftime)
       call read_station_data(source,3,'dra',LVT_rc%dyr(source),LVT_rc%dmo(source), &
            LVT_rc%dda(source), reftime)
       call read_station_data(source,4,'fpk',LVT_rc%dyr(source),LVT_rc%dmo(source), &
            LVT_rc%dda(source), reftime)
       call read_station_data(source,5,'gwn',LVT_rc%dyr(source),LVT_rc%dmo(source), &
            LVT_rc%dda(source), reftime)
       call read_station_data(source,6,'psu',LVT_rc%dyr(source),LVT_rc%dmo(source), &
            LVT_rc%dda(source), reftime)
       call read_station_data(source,7,'sxf',LVT_rc%dyr(source),LVT_rc%dmo(source), &
            LVT_rc%dda(source), reftime)

    endif
	
    call ESMF_TimeSet(surftime1, yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=rc)
    call LVT_verify(rc)
    t = nint((surftime1 - &
         SURFRADobs(source)%starttime)/SURFRADobs(source)%timestep) + 1
    do c=1,7
       call latlon_to_ij(LVT_domain%lvtproj,SURFRADobs(source)%stnlat(c),&
            SURFRADobs(source)%stnlon(c),col,row)
       stn_col = nint(col)
       stn_row = nint(row)
       
       if(SURFRADobs(source)%dw_psp(c,t).ne.LVT_rc%udef) then 
          swd(stn_col,stn_row) = &
               SURFRADobs(source)%dw_psp(c,t)
       endif
       
       if(SURFRADobs(source)%dw_pir(c,t).ne.LVT_rc%udef) then 
          lwd(stn_col,stn_row) = &
               SURFRADobs(source)%dw_pir(c,t)
       endif
       
       if(SURFRADobs(source)%windspd(c,t).ne.LVT_rc%udef) then 
          wspd(stn_col,stn_row) = &
               SURFRADobs(source)%windspd(c,t)
       endif
       
       if(SURFRADobs(source)%pres(c,t).ne.LVT_rc%udef) then 
          psurf(stn_col,stn_row) = &
               SURFRADobs(source)%pres(c,t)
       endif
       
    enddo

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(swd(c,r).lt.-100) then 
             swd(c,r) = LVT_rc%udef
          endif
       enddo
    enddo
    call LVT_logSingleDataStreamVar(LVT_MOC_swdownforc, source, swd,&
         vlevel=1,units="W/m2")
    
    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(lwd(c,r).lt.-100) then 
             lwd(c,r) = LVT_rc%udef
          endif
       enddo
    enddo
    call LVT_logSingleDataStreamVar(LVT_MOC_lwdownforc, source,lwd,&
         vlevel=1,units="W/m2")

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(wspd(c,r).lt.-100) then 
             wspd(c,r) = LVT_rc%udef
          endif
       enddo
    enddo
    call LVT_logSingleDataStreamVar(LVT_MOC_windforc, source, wspd,&
         units="m/s",vlevel=1)

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(wspd(c,r).gt.0) then 
             wspd(c,r) = wspd(c,r)*86.4 
          endif
       enddo
    enddo
    call LVT_logSingleDataStreamVar(LVT_MOC_windforc, source, wspd,&
         units="km/day",vlevel=1)

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(psurf(c,r).lt.-100) then 
             psurf(c,r) = LVT_rc%udef
          else
             psurf(c,r) = psurf(c,r)*100
          endif
       enddo
    enddo
    call LVT_logSingleDataStreamVar(LVT_MOC_psurfforc, source, psurf,&
         units="Pa",vlevel=1)

  end subroutine readSURFRADObs

!BOP
! 
! !ROUTINE: read_station_data
!  \label{read_station_data}
!
! !INTERFACE:
subroutine read_station_data(source, station_index, station_name, &
     year, month, day, &
     reftime)
! 
! !USES:
    use ESMF
    use LVT_timeMgrMod,   	only : LVT_calendar
    use SURFRAD_obsMod, 	only : SURFRADobs
    use LVT_logMod,     	only : LVT_logunit, LVT_getNextUnitNumber, LVT_releaseUnitNumber

    
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine reads and stores the surfrad data for the current day, for a 
!  single station.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUEMENTS:
    integer 			source 
    integer 			station_index
    integer 			year,month, day
    type(ESMF_Time)             :: reftime

!EOP
    integer, parameter 		:: nlines = 1440
    integer 			ftn, icount
    character*3			station_name
    character*24		file_name
    character*100 		directory
    real 			swnet_obs, lwnet_obs
    logical 			file_exists
    
    integer i
    integer dyear,dday,dmonth,jday,dhour,dminute
    character *80 station

    real    latitude,longitude,dt,zen,direct2, elevation
    real    uw_psp, dw_psp
    real    diffuse,dw_pir,dw_casetemp
    real    dw_dometemp,uw_pir,uw_casetemp
    real    uw_dometemp,uvb,par
    real    netsolar,netir,totalnet,temp
    real    rh,windspd,winddir,pressure

    integer    qc_direct,qc_netsolar,qc_netir
    integer    qc_dwpsp,qc_uwpsp,qc_diffuse
    integer    qc_dwpir,qc_dwcasetemp
    integer    qc_dwdometemp,qc_uwpir
    integer    qc_uwcasetemp,qc_uwdometemp
    integer    qc_uvb,qc_par
    integer    qc_totalnet,qc_temp
    integer    qc_rh,qc_windspd,qc_winddir
    integer    qc_pressure
    integer    ios
    logical    readflag
    type(ESMF_Time) :: surftime1
    integer         :: rc    
    
    call get_file_name(year,month,day,station_name,file_name)
    call get_directory(year, station_index, directory)
 
    inquire(file=trim(SURFRADobs(source)%odir)//trim(directory)//file_name, &
         exist=file_exists) 
    if(file_exists) then
       
       ftn = LVT_getNextUnitNumber()
       
       write(LVT_logunit,*) '[INFO] Reading file ',trim(SURFRADobs(source)%odir)//&
            trim(directory)//file_name
       open(ftn,file=(trim(SURFRADobs(source)%odir)//trim(directory)//file_name), &
            status='old')
       
       read (ftn,10) station
       
10     format(1x,a)
       read(ftn,*) latitude, longitude, elevation
       
       icount = 0
       
       readflag = .true. 
       do while(readflag) 
          read(ftn,30,iostat=ios) year,jday,dmonth,dday,dhour,dminute,dt,&
               zen,dw_psp,qc_dwpsp,uw_psp,qc_uwpsp,direct2,&
               qc_direct,diffuse,qc_diffuse,dw_pir,&
               qc_dwpir,dw_casetemp,qc_dwcasetemp,dw_dometemp,&
               qc_dwdometemp,uw_pir,qc_uwpir,uw_casetemp,&
               qc_uwcasetemp,uw_dometemp,qc_uwdometemp,uvb,&
               qc_uvb,par,qc_par,netsolar,qc_netsolar,netir,&
               qc_netir,totalnet,qc_totalnet,temp,qc_temp,rh,&
               qc_rh,windspd,qc_windspd,winddir,qc_winddir,&
               pressure,qc_pressure
          icount = icount + 1
          if(ios.ne.0) then 
             readflag = .false. 
          else
             call ESMF_TimeSet(surftime1, yy=year, &
                  mm = dmonth, &
                  dd = dday, &
                  h = dhour, &
                  m = dminute, &
                  calendar = LVT_calendar, &
                  rc=rc)
             !print *, dw_psp
             if(qc_dwpsp > 0) then!check quality control value
                dw_psp = -9999.0
             endif
             
             if(qc_dwpir >0) then 
                dw_pir = -9999.0
             endif
             
             if(qc_windspd > 0) then 
                windspd = -9999.0
             endif
             
             if(qc_pressure > 0) then
                pressure = -9999.0
             endif

             i = nint((surftime1 - reftime)/SURFRADobs(source)%timestep) + 1
             SURFRADobs(source)%dw_psp(station_index, i) = dw_psp
             SURFRADobs(source)%dw_pir(station_index, i) = dw_pir
             SURFRADobs(source)%windspd(station_index, i) = windspd
             SURFRADobs(source)%pres(station_index, i) = pressure
          endif
          
30        format(1x,i4,1x,i3,4(1x,i2),1x,f6.3,1x,f6.2,20(1x,f7.1,1x,i1))

       enddo
 
!    40    print *,'end of file reached, ',icount,' records read'

       call LVT_releaseUnitNumber(ftn)
	
    else
       SURFRADobs(source)%dw_psp(station_index, :) = -9999.0
       SURFRADobs(source)%dw_pir(station_index, :) = -9999.0
       SURFRADobs(source)%windspd(station_index, :) =  -9999.0
       SURFRADobs(source)%pres(station_index, :) =  -9999.0

       write(LVT_logunit,*) '[ERR] ',trim(SURFRADobs(source)%odir)//trim(directory)//file_name, &
            &'does not exist'
    endif
    
end subroutine read_station_data

!BOP
! 
! !ROUTINE: get_directory
!  \label{get_directory}
!
! !INTERFACE:
subroutine get_directory(year, station_index, directory)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION
!	This subroutine returns the directory which contains the current 
!	datafile given a station and a year
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
! !ARGUMENTS    
    integer year
    integer station_index
    character*100 directory
!EOP
!    print *,'get_directory ', station_index
    character*4 yearC
    character *18 stationC!char len...

    
    WRITE(yearC, '(i4)' ) year

    select case(station_index)
    
        case (1)
            stationC = 'Bondville_IL'
        case (2)
            stationC = 'Boulder_CO'
        case (3)
            stationC = 'Desert_Rock_NV'
        case (4)
            stationC = 'Fort_Peck_MT'
        case (5)
            stationC = 'Goodwin_Creek_MS'
        case (6)
            stationC = 'Penn_State_PA'
        case (7)
            stationC = 'Sioux_Falls_SD'
        case default
            print *, 'bad station index'
	    
    end select

    directory = '/' // trim(stationC) // '/' // trim(yearC) // '/'
    
end subroutine get_directory


!BOP
! 
! !ROUTINE: get_file_name
!  \label{get_file_name}
!
! !INTERFACE:
subroutine get_file_name(year,month,day,station_name, file_name)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION
!	This subroutine returns the file name of the data associated with the 
!given the station, year, month and day
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
!
! !INTERFACE:

! !ARGUMENTS
    integer year,month,day
    character*24 file_name
!EOP    
    
    integer days_passed, yr
    character*3 dayC, station_name
    character*2 yearC

    call get_days_passed(year,month,day,days_passed)
    
    yr = year
    yr = modulo(year,100)!to get the last 2 digits, 2000 is an exception

    if(yr == 0) then
        yearC = '00'
    else if(yr == modulo(year,10)) then
        WRITE(yearC, '(i1)' ) yr
        yearC = '0' // trim(yearC)
    else
         WRITE(yearC, '(i2)' ) yr
    endif

    write(dayC,'(i3.3)') days_passed

    file_name = trim(station_name) // trim(yearC) // trim(dayC) // '.dat'

end subroutine get_file_name


!BOP
! 
! !ROUTINE: get_days_passed
! \label{get_days_passed}
!
! !INTERFACE:
subroutine get_days_passed(year,month,day,days_passed)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION
!  This subroutine calculated the number of days which have passed in 
!  the year give the month, year, and day. Used by get_file_name
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS
    integer year, month, day, days_passed
!EOP
    integer feb
    if(((modulo(year,4) == 0) .AND. (modulo(year,100) /= 0)) .OR. &
    	& (modulo(year,400) == 0)) then
	feb = 29
    else
	feb = 28
    endif    
    
    select case(month)
    	
	case(1) 
    		days_passed = day
	case(2)
		days_passed = day + 31
	case(3) 
		days_passed = day + 31 + feb
	case(4)
		days_passed = day + 62 + feb
	case(5) 
		days_passed = day + 92 + feb
	case(6)
		days_passed = day + 123 + feb
	case(7) 
		days_passed = day + 153 + feb
	case(8) 
		days_passed = day + 184 + feb
	case(9) 
		days_passed = day + 215 + feb
	case(10) 
		days_passed = day + 245 + feb
	case(11) 
		days_passed = day + 276 + feb
	case(12) 
		days_passed = day + 306 + feb
	case default
		print *, 'non gregorian calender'

     end select
   !find skip years, calc month and date from that?
end subroutine get_days_passed

