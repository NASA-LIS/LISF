!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readSMOSL2smObs
! \label{readSMOSL2smObs}
!
! !INTERFACE: 
subroutine readSMOSL2smObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMOSL2sm_obsMod, only : SMOSL2smobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer,   intent(in):: source 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!  11 Dec 2014: Sujay Kumar, Initial Specification
! 
!EOP

  real                    :: timenow
  logical                 :: alarmCheck
  logical                 :: file_exists
  integer                 :: c,r,i,j
  integer                 :: status
  integer                 :: yr, mo, da, hr, mn, ss
  real                    :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: smobs_av(LVT_rc%lnc, LVT_rc%lnr)
  integer                 :: count_av(LVT_rc%lnc, LVT_rc%lnr)
  type(ESMF_Time)         :: cTime, stTime, enTime

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  
  smobs= LVT_rc%udef
  
  if(SMOSL2smobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then
     
     LVT_rc%resetFlag(source) = .false. 
     SMOSL2smobs(source)%startmode = .false.

     if(LVT_rc%smoothObs.eq.1) then 
        smobs_av = 0.0
        count_av = 0 
        call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h  = LVT_rc%dhr(source), &
             m  = LVT_rc%dmn(source), & 
             s  = LVT_rc%dss(source), &
             calendar = LVT_calendar, & 
             rc = status)
        call LVT_verify(status)
        
        stTime = cTime - LVT_obsSmTwL
        enTime = cTime + LVT_obsSmTwL

        cTime = stTime
        do while (cTime.le.enTime) 
           call ESMF_TimeGet(cTime, yy = yr, mm=mo, dd=da,&
                h = hr, m = mn, s =ss, calendar=LVT_calendar, &
                rc=status)
           call LVT_verify(status)

           smobs = LVT_rc%udef
           call processSMOSfiles(source, smobs, yr, mo, da, hr, mn, ss)
           
           cTime = cTime + LVT_obsSmTwI
           
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(smobs(c,r).ne.LVT_rc%udef) then 
                    smobs_av(c,r) = smobs_av(c,r) + smobs(c,r)
                    count_av(c,r) = count_av(c,r) + 1
                 endif
              enddo
           enddo
        enddo

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(count_av(c,r).gt.0) then 
                 smobs(c,r) = smobs_av(c,r)/count_av(c,r)
              else
                 smobs(c,r) = LVT_rc%udef
              endif
           enddo
        enddo
     else
        call processSMOSfiles(source, smobs,LVT_rc%dyr(source), &
             LVT_rc%dmo(source), LVT_rc%dda(source), &
             LVT_rc%dhr(source), LVT_rc%dmn(source), LVT_rc%dss(source))
        
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       smobs,vlevel=1,units="m3/m3")
 
end subroutine readSMOSL2smObs

subroutine processSMOSfiles(source, smobs, iyr, imo, ida, ihr, imn, iss)

  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_logMod
  use SMOSL2sm_obsMod, only : SMOSL2smobs
  
  implicit none

  integer                 :: source
  real                    :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  integer                 :: iyr, imo, ida, ihr, imn, iss

  character*100           :: fname
  character*100           :: smos_filename
  character*8             :: yyyymmdd
  character*200           :: list_files
  integer                 :: sind
  integer                 :: yr,mo,da,hr,mn,ss
  integer                 :: ftn
  type(ESMF_Time)         :: currTime, prevTime
  type(ESMF_Time)         :: obsStartTime,obsEndTime,obsTime
  type(ESMF_TimeInterval) :: dayInt,dt
  integer                 :: status
  integer                 :: ierr

! dump the list of files for the current date to a file (note that
! we assume a flat organization of the files under the SMOS observation
! directory.

  call ESMF_TimeintervalSet(dayInt, d=1,rc=status)
  call LVT_verify(status, 'ESMF_TimeIntervalSet failed in readSMOSL2smObs')

  write(yyyymmdd,'(i4.4,2i2.2)') iyr, imo, ida

  list_files = 'ls '//trim(SMOSL2smobs(source)%odir)//'/*'//trim(yyyymmdd)&
       //'*.DBL > SMOS_filelist.dat'
  call system(trim(list_files))
  
  ftn = LVT_getNextUnitNumber()
  open(ftn,file="./SMOS_filelist.dat",status='old',iostat=ierr)
  do while(ierr.eq.0)
     read(ftn,'(a)',iostat=ierr) smos_filename
     if(ierr.ne.0) then
        exit
     else
        ! check first if the observation time is within the LVT time window (daily)
        call ESMF_TimeSet(currTime, yy=iyr,&
             mm = imo, dd=ida, h = ihr, &
             m = imn, s = iss, calendar = LVT_calendar,&
             rc=status)
        call LVT_verify(status, 'ESMF_TimeSet failed in readSMOSL2smObs (1)')
        prevTime = currTime
        currTime = currTime + dayInt
        sind = index(smos_filename, "SMUDP2_")
        sind = sind+7
        read(smos_filename(sind:sind+3),'(i4)') yr
        read(smos_filename(sind+4:sind+5),'(i2)') mo
        read(smos_filename(sind+6:sind+7),'(i2)') da
        read(smos_filename(sind+9:sind+10),'(i2)') hr
        read(smos_filename(sind+11:sind+12),'(i2)') mn
        read(smos_filename(sind+13:sind+14),'(i2)') ss
        
        call ESMF_TimeSet(obsStartTime, yy=yr,&
             mm = mo, dd=da, h = hr, &
             m = mn, s = ss, calendar = LVT_calendar,&
             rc=status)
        call LVT_verify(status, 'ESMF_TimeSet failed in readSMOSL2smObs (2)')
        read(smos_filename(sind+16:sind+19),'(i4)') yr
        read(smos_filename(sind+20:sind+21),'(i2)') mo
        read(smos_filename(sind+22:sind+23),'(i2)') da
        read(smos_filename(sind+25:sind+26),'(i2)') hr
        read(smos_filename(sind+27:sind+28),'(i2)') mn
        read(smos_filename(sind+29:sind+30),'(i2)') ss
        
        call ESMF_TimeSet(obsEndTime, yy=yr,&
             mm = mo, dd=da, h = hr, &
             m = mn, s = ss, calendar = LVT_calendar,&
             rc=status)
        call LVT_verify(status, 'ESMF_TimeSet failed in readSMOSL2smObs (3)')
        
        dt = (obsEndTime - obsStartTime)/2
        obsTime = (obsStartTime + dt)
        if((obsTime.gt.prevTime).and.(obsTime.le.currTime)) then
           call readSingle_SMOSL2_data(smos_filename, smobs)
        endif
     endif
  enddo

  call LVT_releaseUnitNumber(ftn)

end subroutine processSMOSfiles

!BOP
! 
! !ROUTINE: readSingle_SMOSL2_data
! \label(readSingle_SMOSL2_data)
!
! !INTERFACE:
subroutine readSingle_SMOSL2_data(fname, smobs)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: smobs(LVT_rc%lnc,LVT_rc%lnr)


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads a given SMOS L2 soil moisture data, applies the
!  quality control flags and then grids the data into the LIS/LVT grid. 
!   
!  The data format is described in the SMOS Level 2 auxilary data products
!  specifications document (pages 78 onwards)
!  https://earth.esa.int/documents/10174/127856/SO-TN-IDR-GS-0006_L2+Spec_v6.1_2012-02-09.pdf/ea56dca3-fe04-4196-a200-ab34da105a9c?version=1.0
!  
!  The code here applies available dataflags in the data for inline
!  quality control. Data is excluded when the following criteria is met
!   * RFI probabilty is > 50
!   * vegetation optical thickness is > 0.8
!   * data quality is reported to be poor
!   * when reported uncertainty is high
!
!  For grid points with multiple data values during a day, data at the
!  latest observation time is chosen. The observation times in the 
!  SMOS data is reported as the number of days since Jan 1, 2000. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOSL2 file
!  \item[smobs]    soil moisture data processed to the LVT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer                          :: ftn
  integer                          :: totpts
  integer                          :: i,c,r
  real                             :: col,row
  integer                          :: grid_point_id
  real                             :: lat,lon,alt,sm,sm_dqx
  real                             :: opt_thick, opt_thick_dqx
  real                             :: tsurf,tsurf_dqx,dummy(26)
  integer*2                        :: conflags
  character*1                      :: gqx,chi2,chi2p,rest_conf(36)
  character*1                      :: science(6)
  character*1                      :: proc(4),dgg_current(13),rfi_prob
  integer                          :: acday,acsec,acms
  integer                          :: yr,mo,da,da1,hr,mn,ss
  logical                          :: file_exists
  integer                          :: stc, enc, str, enr, c1,r1
  integer                          :: status
  type(ESMF_Time)                  :: acTime
  type(ESMF_Time)                  :: smTime(LVT_rc%lnc,LVT_rc%lnr)

  inquire(file=fname,exist=file_exists)
  
  if(file_exists) then 
     call ESMF_TimeSet(acTime, yy=2000,&
          mm = 1, dd=1, h =0, &
          m = 0, s = 0, calendar = LVT_calendar,&
          rc=status)
     call LVT_verify(status, 'ESMF_TimeSet failed in readSingle_SMOSL2_data (3)')
     smTime = acTime

     ftn = LVT_getNextUnitNumber()
     write(LVT_logunit,*) '[INFO] Reading '//trim(fname)
     open(ftn,file=fname,form='unformatted',access='stream',&
          convert='little_endian')
     read(ftn) totpts

     do i=1,totpts
        read(ftn)grid_point_id,lat,lon,alt,acday,acsec,acms
        read(ftn)SM, SM_DQX, opt_thick,opt_thick_dqx,tsurf,tsurf_dqx,dummy 
        read(ftn)conflags,gqx,chi2,chi2p,rest_conf  !41 bytes total
        read(ftn)science  ! 6 bytes
        read(ftn)proc     ! 4 bytes
        read(ftn)dgg_current,rfi_prob   !  14 bytes total

        call latlon_to_ij(LVT_domain%lvtproj,lat,lon,&
             col,row)
        c = nint(col)
        r = nint(row)
        stc = max(1,c-2)
        enc = min(LVT_rc%lnc,c+2)
        str = max(1,r-2)
        enr = min(LVT_rc%lnr,r+2)
        do c1=stc,enc
           do r1=str,enr
              if((sm.gt.0.001.and.sm.lt.1.00).and.&
                   (sm_dqx.le.0.1).and.&            !high uncertainty
                   (ichar(gqx).lt.10).and.&         !poor quality
                   (ichar(rfi_prob).lt.50).and.&    !RFI
                   (opt_thick.lt.0.8).and.&
!             (tsurf.gt.273.15).and.&
                   (c1.ge.1.and.c1.le.LVT_rc%lnc).and.&
                   (r1.ge.1.and.r1.le.LVT_rc%lnr)) then 

                 if(smobs(c1,r1).gt.0) then 
!              print*, 'data already there',c,r
                 else
                    if(smobs(c1,r1).eq.-9999.0) then !if data is not present already
                       call SMOS_julhr_date( hr, da, mo, yr, acday*24 ) 
                       call LVT_seconds2time(acsec,da1, hr, mn, ss)
                       
                       call ESMF_TimeSet(acTime, yy=yr,&
                            mm = mo, dd=da, h = hr, &
                            m = mn, s = ss, calendar = LVT_calendar,&
                            rc=status)
                       call LVT_verify(status, 'ESMF_TimeSet failed in readSingle_SMOSL2_data (4)')
                       if(acTime > smTime(c1,r1)) then 
                          smTime(c1,r1) = acTime
                          smobs(c1,r1) = sm
                       endif
                    endif
                 endif
              endif
           enddo
        enddo
     enddo
     call LVT_releaseUnitNumber(ftn)
  endif

end subroutine readSingle_SMOSL2_data

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
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
!BOP
!
! 
! !INTERFACE:    
subroutine SMOS_julhr_date( hour, day, month, year, julhr ) 

  implicit none 
! !ARGUMENTS:   
  integer,  intent(out)       :: day  
  integer,  intent(out)       :: hour 
  integer,  intent(out)       :: month
  integer,  intent(out)       :: year 
  integer,  intent(in)        :: julhr
!
! !DESCRIPTION:
!     uses the julian hour to determine the 2-digit zulu time, day,  
!     month, and 4-digit year.    
!
!     
!    \subsubsection{Method}
!     - determine the current zulu hour  \newline
!     - determine the total number of elapsed days  \newline
!     - count forward to the current day/month/year  \newline
!
!    The arguments are: 
!    \begin{description}
!    \item[hour]
!     the zulu time of the julian hour   
!    \item[day]
!     day of the month (1..31)   
!    \item[month]
!     month of the year (1..12)  
!    \item[year]
!     four digit year
!    \item[julhr]
!     the julian hour being processed
!    \end{description}
!
!EOP  
  integer                     :: dypmon(12)   
  integer                     :: elapdy   
  integer                     :: i
  logical                     :: done 
  
  data dypmon   /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  
!     ------------------------------------------------------------------
!     initialize done flag to false.
!     ------------------------------------------------------------------

  done = .false.
      
!    ------------------------------------------------------------------
!     extract the zulu hour from the julian hour.   
!     ------------------------------------------------------------------
          
  hour = mod(julhr, 24) 

!    ------------------------------------------------------------------
!     determine the number of days that have elapsed since dec 31, 1967 
!     (julian hour 0).  
!     ------------------------------------------------------------------
    
  elapdy = julhr / 24   
          
!     ------------------------------------------------------------------
!     initialize the starting day, month, and year values.  
!     ------------------------------------------------------------------
    
  if (elapdy .gt. 0) then        
     year   = 2000
     day    = 2
     month  = 1
     elapdy = elapdy - 1      
  else       
     print*, 'error in SMOS_julhr_date'
     stop
  endif
      
!    ------------------------------------------------------------------
!     loop through the elapsed days to determine the year.  
!     ------------------------------------------------------------------
    
  do while(.not.done) 

     dypmon(2) = 28  

!    ------------------------------------------------------------------
!     determine if year value is a leap year.  leap years occur every 
!     4 years, with the exception of century years not evenly 
!     divisible by 400.   
!     ------------------------------------------------------------------
    
     if ((mod(year, 4)   .eq. 0) .and.   &
          ((mod(year, 100) .ne. 0) .or. &
          (mod(year, 400) .eq. 0))) then  
        
        dypmon(2) = 29
        
     endif

!     ------------------------------------------------------------------
!     if the elapsed number of days is more than a year's worth,  
!     subtract the appropriate number of days, and increment the year 
!     value.  
!     ------------------------------------------------------------------
    
     if (dypmon(2) .eq. 28) then      
        if (elapdy .ge. 365) then         
           year = year + 1 
           elapdy = elapdy - 365           
        else          
           done = .true.   
        endif
     else     
        if (elapdy .ge. 366) then         
           year = year + 1 
           elapdy = elapdy - 366           
        else          
           done = .true.           
        endif
     endif

!     ------------------------------------------------------------------
!     if the elapsed number of days is less than a year's worth, then   
!     exit loop.
!     ------------------------------------------------------------------    
  enddo

!     ------------------------------------------------------------------
!     count the days and months elapsed in the current year.
!     ------------------------------------------------------------------
    
  do i = 1, elapdy      
     day = day + 1
     if (day .gt. dypmon(month)) then        
        day = 1   
        month = month + 1                
        if(month.gt.12) then 
           month = 1
           year = year + 1
        endif
     endif     
  enddo

  return
  
end subroutine SMOS_julhr_date
