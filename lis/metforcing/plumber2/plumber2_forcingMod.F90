!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module plumber2_forcingMod
!BOP
! !MODULE: plumber2_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various PLUMBER2 stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt plumber2\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 15 Sep 2021: Mark Beauharnois, Derived from Bondville and GSWP2 readers
! 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_plumber2 !defines the native resolution of 
                                       !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: plumber2_struc
!EOP

  type, public :: plumber2_type_dec

     real                     :: ts
     character(len=LIS_CONST_PATH_LEN) :: plumber2file
     real                     :: undef
     real*8                   :: starttime,plumber2time1,plumber2time2
     integer                  :: sYear, sMon, sDay, sHour, sMin, sSec
     integer                  :: eYear, eMon, eDay, eHour, eMin, eSec
     integer                  :: findtime1,findtime2
     integer                  :: utcoffset_sec
     integer                  :: nstns
     character*6,allocatable  :: stnid(:)
     character*3              :: veg_type
     real, allocatable        :: stnlat(:),stnlon(:)
     real, allocatable        :: stnwt(:,:)
     real, allocatable        :: metdata1(:,:) 
     real, allocatable        :: metdata2(:,:) 

     real*8                   :: ringtime

     logical                  :: reset_flag

     integer                  :: read_index  !! PLUMBER2 'time' index
  end type plumber2_type_dec
  
  type(plumber2_type_dec), allocatable :: plumber2_struc(:)

contains

#include "LIS_misc.h"

!BOP
!
! !ROUTINE: init_plumber2
! \label{init_plumber2}
! 
! !INTERFACE:
  subroutine init_plumber2(findex)
! !USES:
    use LIS_coreMod,only     : LIS_rc
    use LIS_logMod, only     : LIS_logunit,LIS_endrun,LIS_verify
    use LIS_timeMgrMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    use netcdf
#endif
    
    implicit none
    integer, intent(in) :: findex
    logical :: file_exists
    character*256 :: startString
    character*4 :: yrStr
    character*2 :: moStr,daStr,hhStr,mnStr,ssStr
    integer :: pYear, pMon, pDay, pHour, pMin, pSec
    integer :: utcYr,utcMo,utcDa,utcHr,utcMn
    integer :: yr1,mo1,da1,hr1,mn1,ss1
    real    :: pLongitude
    real*8  :: utcoffset_seconds
    
! !DESCRIPTION:
!  This routines reads the runtime configurations for using the
!  PLUMBER2 station data. Using the metadata provided for the
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used.
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_plumber2](\ref{readcrd_plumber2}) \newline
!     reads the runtime options specified for PLUMBER2 station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP

    real    :: gmt,gmt1
    integer :: i,n,doy,doy1,yrdays
    integer :: ftn,timeId,longId
    real*8  :: localtime,gmttime,timeDelta,timenow

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the PLUMBER2 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    ! MCB 20211101
    ! Force ONLY 1 nest for PLUMBER2 simulations
    if (LIS_rc%nnest .ne. 1) then
      write(LIS_logunit,*) '[ERR] Only 1 nest allowed for PLUMBER2'
      write(LIS_logunit,*) '[ERR]   simulations at this time'
      write(LIS_logunit,*) '[ERR]   LIS forecast run-time ending.'
      call LIS_endrun()
    endif

    allocate(plumber2_struc(LIS_rc%nnest))

    ! allocate station ID array before reading config entries
    do n=1, LIS_rc%nnest
       plumber2_struc(n)%nstns = 1
       allocate(plumber2_struc(n)%stnid(plumber2_struc(n)%nstns))
    enddo

    call readcrd_plumber2()

    !! MCB 20211101 
    !! Check to make sure user specified correct starting date and time
    !! For this PLUMBER2 case.  We do this by reading the 'units' attribute
    !! of the PLUMBER2 input file "time" variable, and the site longitude value.
    !! Then, using the longitude, figure out the UTC offset in seconds.  The date/time
    !! in UTC is calculated, and then compared to what the user specified in
    !! the LIS config file.

    inquire(file=plumber2_struc(1)%plumber2file,exist=file_exists)
    if(file_exists) then
       write(LIS_logunit,*)'[INFO] Checking PLUMBER2 start...'
       call LIS_verify(nf90_open(path=trim(plumber2_struc(1)%plumber2file),&
          mode=NF90_SHARE, &
          ncid=ftn), &
          'nf90_open failed for PLUMBER2 forcing file in init_plumber2')

       ! Query the 'time' variable
       call LIS_verify(nf90_inq_varid(ftn,'time',timeId), &
          'nf90_inq_varid failed for TIME in init_plumber2')
       call LIS_verify(nf90_get_att(ftn, timeId,'units',startString), &
          'nf90_get_att failed for TIME:units in init_plumber2')

       write(LIS_logunit,*) '[INFO] PLUMBER2 startstring: ',trim(startString)
       
       !! WARNING, WARNING !!
       !! Substring positions are HARD-CODED for extraction of date/time
       !! information from the 'units' attribute of the NetCDF 'time'
       !! variable.  This WILL barf if the 'units' attribute
       !! of different PLUMBER2 input files varies AT ALL positionally.
       !! Example: seconds since 2000-01-01 00:00:00
       !!          ^
       !!          |
       !!          ASSUMES Position #1
       !!
       yrStr = startString(15:18)
       moStr = startString(20:21)
       daStr = startString(23:24)
       hhStr = startString(26:27)
       mnStr = startString(29:30)
       ssStr = startString(32:33)

       read(yrStr, '(i4)') pYear
       read(moStr, '(i2)') pMon
       read(daStr, '(i2)') pDay
       read(hhStr, '(i2)') pHour
       read(mnStr, '(i2)') pMin
       read(ssStr, '(i2)') pSec
       write(LIS_logunit,*) '[INFO] pYear: ',pYear
       write(LIS_logunit,*) '[INFO] pMon : ',pMon
       write(LIS_logunit,*) '[INFO] pDay : ',pDay
       write(LIS_logunit,*) '[INFO] pHour: ',pHour
       write(LIS_logunit,*) '[INFO] pMin : ',pMin
       write(LIS_logunit,*) '[INFO] pSec : ',pSec

       !! Get the PLUMBER2 longitude

       call LIS_verify(nf90_inq_varid(ftn,'longitude',longId), &
          'nf90_inq_varid failed for longitude in init_plumber2')
       call LIS_verify(nf90_get_var(ftn,longId,pLongitude), &
          'nf90_get_var failed for longitude in init_plumber2')

       write(LIS_logunit,*) '[INFO] pLongitude: ',pLongitude

       ! Craig Ferguson's offset calculation...
       if (pLongitude .le. 0.0) then
            utcoffset_seconds = &
              ceiling((pLongitude - 7.5)/15.0) * 3600.0
       else
            utcoffset_seconds = &
              floor((pLongitude + 7.5)/15.0) * 3600.0
       endif

       write(LIS_logunit,*) '[INFO] utcoffset_seconds: ',utcoffset_seconds

       call LIS_date2time(localtime,doy,gmt,pYear,pMon,pDay,pHour,pMin,pSec)

       write(LIS_logunit,*) '[INFO] localtime: ',localtime
       write(LIS_logunit,*) '[INFO] doy      : ',doy
       write(LIS_logunit,*) '[INFO] gmt      : ',gmt

       ! 20211123 CRF, MCB
       ! Use what was entered as a start year in the LIS config file to
       ! determine the fractional time delta
       if((mod(LIS_rc%yr,4).eq.0.and.mod(LIS_rc%yr,100).ne.0) &     !correct for leap year
         .or.(mod(LIS_rc%yr,400).eq.0))then             !correct for y2k
         yrdays=366
       else
         yrdays=365
       endif

       timeDelta = (abs(utcoffset_seconds)/86400.d0) / float(yrdays)

       if (utcoffset_seconds .lt. 0.0) then
          gmttime = localtime + timeDelta
       else
          gmttime = localtime - timeDelta
       endif

       write(LIS_logunit,*) '[INFO] gmttime  : ',gmttime

       call LIS_time2date(gmttime,doy,gmt,utcYr,utcMo,utcDa,utcHr,utcMn)

       write(LIS_logunit,*) '[INFO] UTC Year : ',utcYr
       write(LIS_logunit,*) '[INFO] UTC Mon  : ',utcMo
       write(LIS_logunit,*) '[INFO] UTC Day  : ',utcDa
       write(LIS_logunit,*) '[INFO] UTC Hour : ',utcHr
       write(LIS_logunit,*) '[INFO] UTC Min  : ',utcMn

       ! What is LIS start time?  Find it and compare to what is
       ! expected for the PLUMBER2 dataset
       yr1 = LIS_rc%yr
       mo1 = LIS_rc%mo
       da1 = LIS_rc%da
       hr1 = LIS_rc%hr
       mn1 = LIS_rc%mn
       ss1 = 0

       call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,0.0)

       write(LIS_logunit,*) '[INFO] LIS Start time: ',timenow
 
       if(gmttime .ne. timenow) then
         write(LIS_logunit,*) '[ERR] Start date/time in PLUMBER2 input file'
         write(LIS_logunit,*) '[ERR]   ("time" variable, "units" attribute)'
         write(LIS_logunit,*) '[ERR]   does NOT match LIS start date/time  '
         write(LIS_logunit,*) '[ERR]   LIS simulation aborting...'
         call LIS_endrun()
       endif

    else
       write(LIS_logunit,*) &
       trim(plumber2_struc(1)%plumber2file)//' does not exist'
       call LIS_endrun()
    endif

    call LIS_verify(nf90_close(ftn),&
       'nf90_close failed for PLUMBER2 forcing file in init_plumber2')


       do n=1, LIS_rc%nnest
          !! MCB next line commented for PLUMBER2 because 'ts' is now
          !! set in the 'readcrd_plumber2' routine for automation
          !! plumber2_struc(n)%ts = 1800
          !!call LIS_update_timestep(LIS_rc, n, plumber2_struc(n)%ts)

          plumber2_struc(n)%reset_flag = .false.
       enddo

       ! MCB NOTE: For PLUMBER2 the following met forcing variables are
       !           required:
       !  1. prec rain (total precip)
       !  2. psurf
       !  3. qair
       !  4. tair
       !  5. swdown
       !  6. lwdown
       !  7. wind u
       !  8. wind v
       !  9. lai    (currently not used)
       !  
       LIS_rc%met_nf(findex) = 8 !number of met variables from PLUMBER2

       do n = 1,LIS_rc%nnest

          call LIS_registerAlarm("PLUMBER2 forcing alarm",&
               plumber2_struc(n)%ts, plumber2_struc(n)%ts)

          allocate(plumber2_struc(n)%metdata1(LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(plumber2_struc(n)%metdata2(LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

          plumber2_struc(n)%metdata1 = 0
          plumber2_struc(n)%metdata2 = 0

          plumber2_struc(n)%undef = -9999.0

          plumber2_struc(n)%read_index = 0

          allocate(plumber2_struc(n)%stnlat(plumber2_struc(n)%nstns))
          allocate(plumber2_struc(n)%stnlon(plumber2_struc(n)%nstns))

          allocate(plumber2_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               plumber2_struc(n)%nstns))

          !      MCB Note: No spatial interpolation for PLUMBER2 test case(s)
          !      call compute_stnwts(plumber2_struc(n)%nstns,LIS_rc%gridDesc,&
          !           plumber2_struc(n)%stnlat,plumber2_struc(n)%stnlon,&
          !           LIS_rc%lnc(n)*LIS_rc%lnr(n),plumber2_struc(n)%stnwt)

       enddo

  end subroutine init_plumber2

end module plumber2_forcingMod
