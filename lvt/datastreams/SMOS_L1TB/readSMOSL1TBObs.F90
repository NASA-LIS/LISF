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
! !ROUTINE: readSMOSL1TBObs
! \label{readSMOSL1TBObs}
!
! !INTERFACE: 
subroutine readSMOSL1TBObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMOSL1TB_obsMod, only : SMOSL1TBobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in) :: source
! 
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
  real                    :: smobs(LVT_rc%lnc,LVT_rc%lnr)

  call ESMF_TimeintervalSet(dayInt, d=1,rc=status)
  call LVT_verify(status, 'ESMF_TimeIntervalSet failed in readSMOSL1TBObs')

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  
  smobs= LVT_rc%udef
  
  if(SMOSL1TBobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then
     
     LVT_rc%resetFlag(source) = .false. 
     SMOSL1TBobs(source)%startmode = .false.

! dump the list of files for the current date to a file (note that
! we assume a flat organization of the files under the SMOS observation
! directory.

     write(yyyymmdd,'(i4.4,2i2.2)') LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source)
     list_files = 'ls '//trim(SMOSL1TBobs(source)%odir)//'/*'//trim(yyyymmdd)&
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
           call ESMF_TimeSet(currTime, yy=LVT_rc%dyr(source),&
                mm = LVT_rc%dmo(source), dd=LVT_rc%dda(source), h = LVT_rc%dhr(source), &
                m = LVT_rc%dmn(source), s = LVT_rc%dss(source), calendar = LVT_calendar,&
                rc=status)
           call LVT_verify(status, 'ESMF_TimeSet failed in readSMOSL1TBObs (1)')
           prevTime = currTime
           currTime = currTime + dayInt

           sind = index(smos_filename, "BWLF1C_")
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
           call LVT_verify(status, 'ESMF_TimeSet failed in readSMOSL1TBObs (2)')
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
           call LVT_verify(status, 'ESMF_TimeSet failed in readSMOSL1TBObs (3)')

           dt = (obsEndTime - obsStartTime)/2
           obsTime = (obsStartTime + dt)
           if((obsTime.gt.prevTime).and.(obsTime.le.currTime)) then
              call ESMF_TimeGet(obsTime, yy=yr,&
                   mm = mo, dd=da, h = hr, &
                   m = mn, s = ss, calendar = LVT_calendar,&
                   rc=status)
              call read_SMOSL1_data(smos_filename, yr,mo,da,hr,mn,ss,smobs)
           endif
        endif
     enddo

     call LVT_releaseUnitNumber(ftn)

!     do r=1,LVT_rc%lnr
!        do c=1,LVT_rc%lnc
!           if(smobs(c,r).ne.-9999.0) then
!              SMOSL1TBobs(source)%smobs(c,r) = smobs(c,r)
!           endif
!        enddo
!     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist,source,&
       smobs,vlevel=1,units="m3/m3")
 
end subroutine readSMOSL1TBObs


!BOP
! 
! !ROUTINE: read_SMOSL1_data
! \label(read_SMOSL1_data)
!
! !INTERFACE:
subroutine read_SMOSL1_data(fname,yr,mo,da,hr,mn,ss, smobs)
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
  integer                       :: yr,mo,da,hr,mn,ss

  character (len=*)             :: fname
  real                          :: smobs(LVT_rc%lnc,LVT_rc%lnr)


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads a given SMOS L1 TB data, applies the
!  quality control flags and then grids the data into the LIS/LVT grid. 
!   
!  BWLF1C - Level 1C browse full polarization land science measurements
!  product
!  The data format is described in the SMOS Level 1 auxilary data products
!  specifications document (pages 313 onwards)
!  http://www.cesbio.ups-tlse.fr/SMOS_blog/wp-content/uploads/DOCS/SO-TN-IDR-GS-0005_L1_Spec_v5.7_2009-04-02.pdf
!
!Polarisation flags:
! [XXXX:XXXX:XXXX:XX00]representsHHpolarisation
! [XXXX:XXXX:XXXX:XX01]representsVVpolarisation
! SUN FOV flags
! [X X X X:X X X X:X X X X:X 0 X X] means that no Direct Sun correction has been performed during image reconstruction of this pixel
! [XXXX:XXXX:XXXX:X1XX]meansthatDirectSuncorrectionhasbeenperformedduringimagereconstructionof this pixel
!  
!SUN GLINT FOV flag:
! [X X X X:X X X X:X X X X:0 X X X] means that no Reflected Sun correction has been performed during image
!reconstruction of this pixel
! [X X X X:X X X X:X X X X:1 X X X] means that Reflected Sun correction has been performed during image reconstruction of this pixel. Sun correction is based on the Sea bistatic coefficients defined in the AUX_BSCAT ADF and computed for a fixed wind speed of 7 m/s and wind direction of 0 deg North
! MOON FOV flag:
! [X X X X:X X X X:X X X 0:X X X X] means that no Direct Moon correction has been performed during image
!reconstruction of this pixel
! [X X X X:X X X X:X X X 1:X X X X] means that Direct Moon correction has been performed during image reconstruction of this pixel
!
! SINGLE_SNAPSHOT flag:
![X X X X:X X X X:X X 0 X:X X X X] means that this scene has been combined with an adjacent scene in opposite polarisation during image reconstruction to account for crosspolarisation leakage
! [X X X X:X X X X:X X 1 X:X X X X] means that this scene has not been combined with an adjacent scene in opposite polarisation during image reconstruction to account for crosspolarisation leakage (it has been processed with only co- polar antenna patterns information)
! RFI Mitigation flag:
! [XXXX:XXXX:X0XX:XXXX]meansthatnoRFIMitigationhasbeenperformedduringimagereconstructionofthis
!pixel measurement
! [X X X X:X X X X:X 1 X X:X X X X] means that RFI Mitigation has been performed during image reconstruction of this pixel measurement
! SUN POINT flag:
! [XXXX:XXXX:0XXX:XXXX]meansthatthispixelisnotlocatedinazone(seebelow)whereaSunaliaswas
!reconstructed
! [XXXX:XXXX:1XXX:XXXX]meansthatthispixelislocatedinazone(circlearoundSunaliaspositionwithradius configurable through Sun_Point_Flag_Size field in AUX CNFL1P) where a Sun alias was reconstructed (if Sun removal is active, measurement may be degraded)
!
!
!SUN POINT flag:
! [XXXX:XXXX:0XXX:XXXX]meansthatthispixelisnotlocatedinazone(seebelow)whereaSunaliaswas
!reconstructed
! [XXXX:XXXX:1XXX:XXXX]meansthatthispixelislocatedinazone(circlearoundSunaliaspositionwithradius configurable through Sun_Point_Flag_Size field in AUX CNFL1P) where a Sun alias was reconstructed (if Sun removal is active, measurement may be degraded)
! SUN GLINT_AREA flag:
! [XXXX:XXX0:XXXX:XXXX]meansthatthispixelislocatedinazonewherenoSunreflectionhasbeendetected
! [X X X X:X X X 1:X X X X:X X X X] means that this pixel is located in a zone where Sun reflection has been detected using the bi-static scattering coefficient threshold defined in the configuration file.
! MOON POINT flag:
! [XXXX:XX0X:XXXX:XXXX]meansthatthispixelislocatedinazonewherenoMoonaliaswasreconstructed
! [X X X X:X X 1 X:X X X X:X X X X] means that this pixel is located in a zone where a Moon alias was reconstructed (after Moon removal, measurement may be degraded)
!
! AF FOV flag:
![XXXX:X0XX:XXXX:XXXX]meansthatthepixelisnotinsidetheexclusivezoneofAliasfree(delimitedbythesix aliased unit circles)
! [X X X X:X 1 X X:X X X X:X X X X] means that the pixel is inside the exclusive zone of Alias free (delimited by the six aliased unit circles)
! BORDER FOV flag:
! [X X X 0:X X X X:X X X X:X X X X] means that the pixel is far from the border delimiting the Extended Alias free zone
! SUN TAILS flag:
!and from the unit circle replicas borders (also known as “suspenders and belts”)
! [XXX1:XXXX:XXXX:XXXX]meansthatthepixelisclosetotheborderdelimitingtheExtendedAliasfreezoneor to the unit circle replicas borders (also known as “suspenders and belts”). Distance threshold is configurable through FOV_Border_Flag_Size field in AUX CNFL1P
! [XX0X:XXXX:XXXX:XXXX]meansthatthispixelislocatedinazonewithnopotentialproblemswithSunaliases
! [X X 1 X:X X X X:X X X X:X X X X] means that this pixel is located in the hexagonal alias directions centred on a Sun alias (if Sun is not removed, measurement may be degraded in these directions)
! RFI flags:
! [X 0 X X:X X X X:X X X X:X X X X] means that the measurement is not affected by strong RFI as detected in the L1b
!processing
! [X 1 X X:X X X X:X X X X:X X X X] means that the measurement is affected by strong RFI as detected in the L1b processing (max RFI self estimated BT above a defined threshold)
! [0 X X X:X X X X:X X X X:X X X X] means that the measurement is not affected by any point source RFI as identified in the AUX RFI list and it does not exceed the threshold defined in BT_Dual/Full_RFI_Pixel_Flag_Threshold fields in AUX CNFL1P
! [1XXX:XXXX:XXXX:XXXX]meansthatthemeasurementisaffectedbypointsourceRFIasidentifiedintheAUX RFI list (flag is set in a circle around the RFI position, with a radius dependant on the RFI expected BT defined in the AUX RFILST), or it has exceeded the threshold defined in BT_Dual/Full_RFI_Pixel_Flag_Threshold fields in AUX CNFL1P, or is negative for Dual Polarisation values
!
![XXXX:0XXX:XXXX:XXXX]meansthatthepixelisnotaffectedbyanytailsfrompointsourceRFIasidentifiedin the AUX RFILST ADF
! [XXXX:1XXX:XXXX:XXXX]meansthatthemeasurementisaffectedbythetailsofapointsourceRFIas identified in the AUX RFI list (tail width is dependant on the RFI expected BT defined in the AUX RFILST)
!  For grid points with multiple data values during a day, data at the
!  latest observation time is chosen. The observation times in the 
!  SMOS data is reported as the number of days since Jan 1, 2000. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOSL1 file
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
  integer                          :: i,k,c,r
  real                             :: col,row
  integer                          :: grid_point_id
  real                             :: lat,lon,alt,sm,sm_dqx
  integer*1                        :: mask,bt_counter
  integer*2                        :: rad_acc_pix
  integer*2                        :: flags
  real                             :: bt_value(4)
  integer*2                        :: az_angle, ft_axis1, ft_axis2
  logical                          :: file_exists
  integer                          :: stc, enc, str, enr, c1,r1
  integer                          :: status
  type(ESMF_Time)                  :: acTime
  logical                          :: rfi_pass
  type(ESMF_Time)                  :: smTime(LVT_rc%lnc,LVT_rc%lnr)

  inquire(file=fname,exist=file_exists)
  
  if(file_exists) then 
     call ESMF_TimeSet(acTime, yy=2000,&
          mm = 1, dd=1, h =0, &
          m = 0, s = 0, calendar = LVT_calendar,&
          rc=status)
     call LVT_verify(status, 'ESMF_TimeSet failed in read_SMOSL1_data (3)')
     smTime = acTime

     ftn = LVT_getNextUnitNumber()
     write(LVT_logunit,*) '[INFO] Reading '//trim(fname)
     open(ftn,file=fname,form='unformatted',access='stream',&
          convert='little_endian')
     read(ftn) totpts

     do i=1,totpts
        read(ftn) grid_point_id,lat, lon, alt, mask, bt_counter
        do k=1,bt_counter
           read(ftn) flags, bt_value(k), rad_acc_pix, az_angle, ft_axis1, &
                ft_axis2
        enddo
!        write(*,fmt='(16(I1.1,1X),F7.2)') ibits(flags,0,1),ibits(flags,1,1), &
!             ibits(flags,2,1),ibits(flags,3,1),&
!             ibits(flags,4,1),ibits(flags,5,1),ibits(flags,6,1),&
!             ibits(flags,7,1),ibits(flags,8,1),&
!             ibits(flags,9,1),ibits(flags,10,1), ibits(flags,11,1),&
!             ibits(flags,12,1),ibits(flags,13,1),ibits(flags,14,1),&
!             ibits(flags,15,1),bt_value(1)
        call latlon_to_ij(LVT_domain%lvtproj,lat,lon,&
             col,row)
        c = nint(col)
        r = nint(row)
        stc = max(1,c-2)
        enc = min(LVT_rc%lnc,c+2)
        str = max(1,r-2)
        enr = min(LVT_rc%lnr,r+2)

        rfi_pass = .true.
!        if(ibits(flags,1,1).eq.0) then 
!           rfi_pass = .true. 
!        else
!           rfi_pass = .false. 
!        endif
!        if(.not.rfi_pass) then 
!           if(ibits(flags,9,1).eq.1) then  !RFI mitigation done
!              rfi_pass = .true. 
!           endif
!        endif
        do c1=stc,enc
           do r1=str,enr
              if(((ibits(flags,15,1).eq.0.and.ibits(flags,14,1).eq.0)).and.& !only HH 
                   rfi_pass.and.&
                   (c1.ge.1.and.c1.le.LVT_rc%lnc).and.&
                   (r1.ge.1.and.r1.le.LVT_rc%lnr)) then 

                 if(smobs(c1,r1).gt.0) then 
!              print*, 'data already there',c,r
                 else
                    if(smobs(c1,r1).eq.-9999.0) then !if data is not present already
                       call ESMF_TimeSet(acTime, yy=yr,&
                            mm = mo, dd=da, h = hr, &
                            m = mn, s = ss, calendar = LVT_calendar,&
                            rc=status)
                       call LVT_verify(status, 'ESMF_TimeSet failed in read_SMOSL1_data (4)')
                       if(acTime > smTime(c1,r1)) then 
                          smTime(c1,r1) = acTime
                          smobs(c1,r1) = bt_value(1)
                       endif
                    endif
                 endif
              endif
           enddo
        enddo

     enddo
     call LVT_releaseUnitNumber(ftn)

!     open(100,file='smobs.bin',form='unformatted')
!     write(100) smobs
!     close(100)
!     stop
  endif

end subroutine read_SMOSL1_data

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
subroutine SMOS_julhr_date1( hour, day, month, year, julhr ) 

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
  
end subroutine SMOS_julhr_date1
