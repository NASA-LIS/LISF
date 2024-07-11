!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: USAF_fldbld_radflux_galwem
! \label{USAF_fldbld_radflux_galwem}
!
! !REVISION HISTORY:
! 11 Aug 2016  Initial specification based on AGRMET_fldbld_precip_gfs
!              ...........................................James Geiger/NASA
! 29 May 2024  Added reading in GALWEM radiation fluxes
!              ...........................................K. Arsenault/SAIC
!
! !INTERFACE:
subroutine USAF_fldbld_radflux_galwem(n,julhr,fg_swdata,fg_lwdata,rc)

! !USES:
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_coreMod,       only : LIS_rc, LIS_masterproc
  use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_alert, &
                                LIS_verify, LIS_endrun
  use LIS_timeMgrMod,    only : LIS_julhr_date

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer             :: julhr
  real, intent(out)   :: fg_swdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real, intent(out)   :: fg_lwdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
  integer, intent(out):: rc

!
! !DESCRIPTION:
!  This routine reads in the GALWEM downward short- 
!   and long-wave radiation flux data and interpolates 
!   to the USAF/AGRMET grid.
!
! The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fc\_hr]
!    forecast hour or the difference between reference and valid time
!  \item[fg\_swdata]
!    array of galwem shortwave radiation data
!  \item[fg\_lwdata]
!    array of galwem longwave radiation data
!
!  \item[ftn]
!    file unit number
!  \item[igrib]
!    GRIB message handle
!  \item[yr1, mo1, da1, hr1]
!    date/time
!  \item[avnfile]
!    name of the 3 hour forecast file
!  \item[avnfile2]
!   name of the 6 hour forecast file
!  \item[message]
!    error message
!  \item[iginfo]
!   array of grid information from a grib record
!  \item[gridres]
!    the resolution, in degrees, of the first
!    guess grid
!  \item[fg\_swdown]
!    SWdown radiation flux data to be interpolated to lis grid
!  \item[fg\_lwdown]
!    LWdown radiation flux data to be interpolated to lis grid
!  \item[fg\_swdown1]
!    3 hour forecast swdown radiation flux data
!  \item[fg\_lwdown1]
!    3 hour forecast lwdown radiation flux data
!  \item[alert\_number]
!    number of alerts that occur in the program
!  \item[ifguess]
!    east/west dimension of first guess grid
!  \item[jfguess]
!    north/south dimension of first guess grid
!  \item[center]
!    meteorological center that produced the first
!    guess data (7-NCEP, 57-AFGWC, 58-FNMOC)
!  \item[ierr]
!    error code
!  \item[yr\_2d]
!    2 digit year for comparison with GRIB header
!  \item[found]
!    logical flag set true when an input file with the correct valid
!    time is found
!  \item[file\_julhr]
!    julian hour used to determine names of forecast files from previous cycles
!  \item[getsixhr]
!    indicates whether to get data from the 6 hour forecast file
!  \item[dataDate]
!    date of values in the GRIB message
!  \item[dataTime]
!    time of values in the GRIB message
!  \item[gtype]
!    type of grid determined by querying GRIB message
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!    converts the julian hour to a date format
!  \item[AGRMET_getGALWEMfilename](\ref{AGRMET_getGALWEMfilename}) \newline
!    generates the first guess GALWEM filename
!  \item[AGRMET\_fldbld\_read\_radflux\_galwem]
!   (\ref{USAF_fldbld_read_radflux_galwem}) \newline
!    read GALWEM radiation flux data in grib format
!  \item[interp\_galwem\_first\_guess](\ref{interp_galwem_first_guess}) \newline
!    interpolate first guess data to the USAF/AGRMET grid
!  \item[LIS\_abort](\ref{LIS_abort}) \newline
!    abort in case of error
!  \end{description}
!EOP
  character(255) :: message(20)
  integer                 :: ftn, igrib
  character*250           :: avnfile
  integer                 :: yr1, mo1, da1, hr1
  integer                 :: fc_hr
  integer                 :: iginfo      ( 2 )
  real                    :: gridres
  integer, save           :: alert_number = 1
  real, allocatable       :: fg_swdown1   ( : , : )
  real, allocatable       :: fg_lwdown1   ( : , : )
  integer                 :: ifguess, jfguess
  integer                 :: center
  integer                 :: ierr
  logical*1               :: found
  logical                 :: first_time
  integer                 :: yr_2d
  integer                 :: file_julhr
  integer                 :: dataDate, dataTime
  character*100           :: gtype
  logical                 :: found_inq

  ! External functions
  external :: AGRMET_getGALWEMfilename
  external :: interp_galwem_first_guess
  external :: USAF_fldbld_read_radflux_galwem

  ! ---------------------------------------------------

  ! Initialize return code to "no error".  We will change it below if 
  ! necessary.
  rc = 0

  ! Will search previous GALWEM cycles every six hours, up to 24 hours,
  !  until we find an acceptable file.
  fc_hr = 0           ! Incremented below
  file_julhr = julhr  ! Decremented below
  call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)

  if (hr1 == 1 .or. hr1 == 7 .or. hr1 == 13 .or. hr1 == 19) then
     fc_hr = 1
  else if (hr1 == 2 .or. hr1 == 8 .or. hr1 == 14 .or. hr1 == 20) then
     fc_hr = 2
  else if (hr1 == 3 .or. hr1 == 9 .or. hr1 == 15 .or. hr1 == 21) then
     fc_hr = 3
  else if (hr1 == 4 .or. hr1 == 10 .or. hr1 == 16 .or. hr1 == 22) then
     fc_hr = 4
  else if (hr1 == 5 .or. hr1 == 11 .or. hr1 == 17 .or. hr1 == 23) then
     fc_hr = 5
  end if
  file_julhr = file_julhr - fc_hr

  found = .FALSE.
  first_time = .true.
  do while ( .not. found )

     ! Make sure we start with the previous GALWEM cycle.
     if ( (.not. first_time) .or. &
          (first_time .and. fc_hr < 6)) then
        fc_hr = fc_hr + 6
        if (fc_hr > 30) exit ! Give up

        file_julhr = file_julhr - 6
        call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)
     end if
     first_time = .false.

     yr_2d = mod(yr1,100)
     if(yr_2d.eq.0) yr_2d = 100

     call AGRMET_getGALWEMfilename(avnfile, agrmet_struc(n)%agrmetdir,&
              agrmet_struc(n)%galwemraddir, agrmet_struc(n)%use_timestamp,&
              agrmet_struc(n)%galwem_res, yr1,mo1,da1,hr1,fc_hr)

!     ------------------------------------------------------------------
!     open first guess grib data using library utility.  just read
!     the first file only, as all data will be of the same type
!     (avn or nogaps) because the search script ensures that it is.
!     ------------------------------------------------------------------
     inquire(file=trim(avnfile),exist=found_inq)
     if (.not. found_inq) then
        write(LIS_logunit,*) '[WARN] Cannot find file '//trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = '  Cannot find file ' // trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        cycle
     end if

#if (defined USE_GRIBAPI)

     ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
     ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
     ! writing error messages to stdout/stderr, which may lead to runtime
     ! problems.

     call grib_open_file(ftn,trim(avnfile),'r',ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) &
             '[WARN] Failed to open GALWEM radiation file ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = '  Cannot open file ' // trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        cycle
     end if

     call grib_new_from_file(ftn,igrib,ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) &
             '[WARN] Failed file read check '//trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = '  Failed to read field from ' // trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'centre',center,ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'error in grib_get: centre in ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = '  Failed to read centre from ' // trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     end if

     call grib_get(igrib,'gridType',gtype,ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'error in grid_get: gridType in ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = '  Failed to read gridType from ' // trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'Ni',iginfo(1),ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'error in grid_get:Ni in ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = '  Failed to read Ni from ' // trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'Nj',iginfo(2),ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'error in grid_get:Nj in ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = '  Failed to read Nj from ' // trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) &
             'error in grid_get:jDirectionIncrementInDegrees in ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = &
                '  Failed to read jDirectionIncrementInDegrees from ' // &
                trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'dataDate',dataDate,ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'error in grid_get:dataDate in ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = &
                '  Failed to read dataDate from ' // &
                trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'dataTime',dataTime,ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'error in grid_get:dataTime in ' // &
             trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = &
                '  Failed to read dataTime from ' // &
                trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     if ( yr1*10000+mo1*100+da1 /= dataDate .or. &
          hr1*100 /= dataTime ) then
        write(LIS_logunit,*) &
             '[WARN] Bad time found in ', trim(avnfile)
        write(LIS_logunit,*) &
             '[WARN] Found: ', dataDate, dataTime
        write(LIS_logunit,*) &
             '[WARN] Expected ', (yr1*10000+mo1*100+da1), &
             (hr1*100)
        flush(LIS_logunit)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = &
                '  Bad date and time from ' // &
                trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     if ( gtype /= "regular_ll" ) then
        write(LIS_logunit,*) &
             '[WARN] Did not find lat/long data in ', trim(avnfile)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_radflux_galwem.'
           message(3) = &
                '  Lat/lon grid not found in ' // &
                trim(avnfile)
           call LIS_alert( 'USAF_fldbld_radflux_galwem', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_close_file(ftn)
        cycle
     endif

     call grib_release(igrib,ierr)
     call grib_close_file(ftn)

#else
     write(LIS_logunit,*) '[ERR]: USAF_fldbld_radflux_galwem requires GRIB-API'
     write(LIS_logunit,*) '[ERR]: please recompile LIS'
     call LIS_endrun
#endif

     write(LIS_logunit,*) &
          '[INFO] Using NWP Radiation fields from ',trim(avnfile)
     rc = 0

     write(LIS_logunit,*)'[INFO] FIRST GUESS DATA IS ON A ', gridres,&
          ' DEGREE LAT/LON GRID'
     ifguess = iginfo(1)
     jfguess = iginfo(2)

     if (center .eq. 57) then
        write(LIS_logunit,*) &
             '[INFO] RADIATION DATA IS FROM UK UM (GALWEM) MODEL'
     else
        write(LIS_logunit,*)'[INFO] UNKNOWN SOURCE FOR RADIATION'
     end if

! ------------------------------------------------------------------
! allocate first guess grid-specific variables.
! ------------------------------------------------------------------
     allocate ( fg_swdown1 (ifguess, jfguess) )
     allocate ( fg_lwdown1 (ifguess, jfguess) )

! ------------------------------------------------------------------
! read in first guess data for this julian hour.
! ------------------------------------------------------------------
     call USAF_fldbld_read_radflux_galwem(avnfile, ifguess, jfguess,&
          fg_swdown1, fg_lwdown1, rc)

     if (rc /= 0) then
        deallocate ( fg_swdown1 )
        deallocate ( fg_lwdown1 )
        cycle
     end if

     ! Interpolate
     call interp_galwem_first_guess(n, ifguess, jfguess, .true., &
          fg_swdown1, fg_swdata)
     call interp_galwem_first_guess(n, ifguess, jfguess, .true., &
          fg_lwdown1, fg_lwdata)
     deallocate ( fg_swdown1 )
     deallocate ( fg_lwdown1 )

     ! At this point, we are done.
     found = .true.
     if (found) exit

  enddo ! Loop through cycles and forecast hours

  ! Give up if no acceptable GALWEM file was found
  if (.not. found) then
     write(LIS_logunit,*)'[WARN] No matching GALWEM Radiation file found!'
     rc = 1
     return
  end if

end subroutine USAF_fldbld_radflux_galwem


!BOP
!
! !ROUTINE: USAF_fldbld_read_radflux_galwem
!  \label{USAF_fldbld_read_radflux_galwem}
!
! !REVISION HISTORY:
! 11 Aug 2016  Initial specification based on AGRMET_fldbld_read_precip_gfs
!              ...........................................James Geiger/NASA
! 29 May 2024  Added reading in GALWEM radiation fluxes
!              ...........................................K. Arsenault/SAIC
!
! !INTERFACE:
subroutine USAF_fldbld_read_radflux_galwem(fg_filename, &
     ifguess, jfguess, fg_swdown, fg_lwdown, rc )

  ! !USES:
  use LIS_coreMod, only : LIS_masterproc
  use LIS_logMod,  only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none

! !ARGUMENTS:
  character(len=*),  intent(in) :: fg_filename
  integer,        intent(in)    :: ifguess
  integer,        intent(in)    :: jfguess
  real,           intent(out)   :: fg_swdown ( ifguess,jfguess )
  real,           intent(out)   :: fg_lwdown ( ifguess,jfguess )
  integer,        intent(out) :: rc
!
! !DESCRIPTION:
!
!     Read UK Unified Model (GALWEM) radiation data in GRIB-2 format.
!
!     \textbf{Method} \newline
!
!     - open file and allocate variables. \newline
!     - read in each grib record. \newline
!     - retrieve section 1 information. \newline
!     - if data is what we need, store it in the proper arrays. \newline
!       also, keep track of what has been read in using \newline
!       counter variables. \newline
!     - check forecast hour of data, if it is not the 3 or \newline
!       or 6 hour forecast, then send an alert message \newline
!       to warn of possible degradation. \newline
!     - check counter variables.  if data is missing abort. \newline
!
!
!     \begin{description}
!      \item[fg\_filename]
!        name, including path, of the first guess file
!        being read in
!      \item[ifguess]
!        east-west dimension of first guess grid
!      \item[jfguess]
!        north-south dimension of first guess grid
!      \item[fg\_swdown]
!        SWdown radiation flux read in from galwem file
!      \item[fg\_lwdown]
!        LWdown radiation flux read in from galwem file
!       \item[alert\_number]
!        counts number of alert messages sent
!       \item[cstat]
!        I/O status, character
!       \item[message]
!        Error message
!       \item[count\_swdown]
!        counts number of SWdown radiation flux
!        levels read in from first guess file
!       \item[count\_lwdown]
!        counts number of LWdown radiation flux
!        levels read in from first guess file
!       \item[i,j,k]
!        looping and indexing variables
!       \item[ierr,istat1]
!        error status
!       \item[dum1d]
!        dummy array
!       \end{description}
!
!EOP
  character(255) :: message(20)
  integer, save :: alert_number = 1
  integer                       :: count_swdown
  integer                       :: count_lwdown
  integer                       :: ierr
  integer :: k
  integer                       :: ftn, igrib, nvars
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, forecasttime_val
  real,           allocatable   :: dum1d    ( : )
  logical                       :: found_inq
  integer :: productDefinitionTemplateNumber

! ------------------------------------------------------------------

  rc = 0

! ------------------------------------------------------------------
!   read in grib file.
! ------------------------------------------------------------------

  ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
  ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
  ! writing error messages to stdout/stderr, which may lead to runtime
  ! problems.

  inquire(file=trim(fg_filename),exist=found_inq)
  if (.not. found_inq) then
     write(LIS_logunit,*)'[WARN] Cannot find ' // trim(fg_filename)
     if (LIS_masterproc) then
        message(:) = ''
        message(1) = '[WARN] Program: LIS'
        message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
        message(3) = '  Cannot find file ' // trim(fg_filename)
        call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
             alert_number, message)
        alert_number = alert_number + 1
     end if
     rc = 1
     return
  end if

#if (defined USE_GRIBAPI)
  call grib_open_file(ftn, trim(fg_filename), 'r', ierr)
  if (ierr /= 0) then
     write(LIS_logunit,*)'[WARN] Cannot open ' // trim(fg_filename)
     if (LIS_masterproc) then
        message(:) = ''
        message(1) = '[WARN] Program: LIS'
        message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
        message(3) = '  Cannot open file ' // trim(fg_filename)
        call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
             alert_number, message)
        alert_number = alert_number + 1
     end if
     call grib_close_file(ftn)
     rc = 1
     return
  end if

  write(LIS_logunit,*)' '
  write(LIS_logunit,*) &
       '[INFO] Reading GALWEM radiation fluxes '
  write(LIS_logunit,*) trim(fg_filename)

  call grib_count_in_file(ftn, nvars, ierr)
  if (ierr /= 0) then
     write(LIS_logunit,*) &
          '[WARN] Problem counting GRIB messages in ' //trim(fg_filename)
     if (LIS_masterproc) then
        message(:) = ''
        message(1) = '[WARN] Program: LIS'
        message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
        message(3) = '  Problem counting GRIB messages in ' // &
             trim(fg_filename)
        call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
             alert_number, message)
        alert_number = alert_number + 1
     end if
     call grib_close_file(ftn)
     rc = 1
     return
  end if

  allocate ( dum1d   (ifguess*jfguess) )
  count_swdown = 0
  count_lwdown = 0

  do k = 1, nvars

     call grib_new_from_file(ftn, igrib, ierr)
     if (ierr /= 0) then
        write(LIS_logunit,*) &
             '[WARN] Cannot get GRIB message in ' // trim(fg_filename)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
           message(3) = '  Cannot get GRIB message in ' // &
                trim(fg_filename)
           call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        rc = 1
        exit
     end if

     call grib_get(igrib, 'discipline', param_disc_val, ierr)
     if (ierr /= 0) then
        write(LIS_logunit,*) &
             '[WARN] Cannot get discipline in ' // trim(fg_filename)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
           message(3) = '  Cannot read discipline in ' // &
                trim(fg_filename)
           call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_release(igrib, ierr)
        rc = 1
        exit
     end if

     call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
     if (ierr /= 0) then
        write(LIS_logunit,*) &
             '[WARN] Cannot get parameterCategory in ' // &
             trim(fg_filename)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
           message(3) = '  Cannot read parameterCategory in ' // &
                trim(fg_filename)
           call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_release(igrib, ierr)
        rc = 1
        exit
     end if

     call grib_get(igrib, 'parameterNumber', param_num_val, ierr)
     if (ierr /= 0) then
        write(LIS_logunit,*) &
             '[WARN] Cannot get parameterNumber in ' // &
             trim(fg_filename)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
           message(3) = '  Cannot read parameterNumber in ' // &
                trim(fg_filename)
           call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_release(igrib, ierr)
        rc = 1
        exit
     end if

     call grib_get(igrib, 'productDefinitionTemplateNumber', &
          productDefinitionTemplateNumber, ierr)
     if (ierr /= 0) then
        write(LIS_logunit,*) &
             '[WARN] Cannot get productDefinitionTemplateNumber in ' // &
             trim(fg_filename)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
           message(3) = &
                '  Cannot read productDefinitionTemplateNumber in ' // &
                trim(fg_filename)
           call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_release(igrib, ierr)
        rc = 1
        exit
     end if

     call grib_get(igrib, 'forecastTime', forecasttime_val, ierr)
     if (ierr /= 0) then
        write(LIS_logunit,*) &
             '[WARN] Cannot get forecastTime in ' // &
             trim(fg_filename)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
           message(3) = &
                '  Cannot read forecastTime in ' // &
                trim(fg_filename)
           call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        call grib_release(igrib, ierr)
        rc = 1
        exit
     end if

     !  Instantaneous SW Down Radiation flux:
     if ( param_disc_val == 0 .and. param_cat_val == 4 .and. &
          param_num_val  == 7 .and. &
          productDefinitionTemplateNumber == 0) then

        call grib_get(igrib, 'values', dum1d, ierr)
        if (ierr /= 0) then
           write(LIS_logunit,*) &
                '[WARN] Cannot get values in ' // &
                trim(fg_filename)
           if (LIS_masterproc) then
              message(:) = ''
              message(1) = '[WARN] Program: LIS'
              message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
              message(3) = &
                   '  Cannot read values in ' // &
                   trim(fg_filename)
              call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                   alert_number, message)
              alert_number = alert_number + 1
           end if
           call grib_release(igrib,ierr)
           rc = 1
           exit
        end if
        fg_swdown = reshape(dum1d, (/ifguess,jfguess/))
        count_swdown = count_swdown + 1
     end if

     !  Instantaneous LW Down Radiation flux:
     if ( param_disc_val == 0 .and. param_cat_val == 5 .and. &
          param_num_val  == 3 .and. &
          productDefinitionTemplateNumber == 0) then

        call grib_get(igrib, 'values', dum1d, ierr)
        if (ierr /= 0) then
           write(LIS_logunit,*) &
                '[WARN] Cannot get values in ' // &
                trim(fg_filename)
           if (LIS_masterproc) then
              message(:) = ''
              message(1) = '[WARN] Program: LIS'
              message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
              message(3) = &
                   '  Cannot read values in ' // &
                   trim(fg_filename)
              call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                   alert_number, message)
              alert_number = alert_number + 1
           end if
           call grib_release(igrib, ierr)
           rc = 1
           exit
        end if
        fg_lwdown = reshape(dum1d, (/ifguess,jfguess/))
        count_lwdown = count_lwdown + 1
     end if

     ! Done with this GRIB message
     call grib_release(igrib, ierr)
     if (ierr /= 0) then
        write(LIS_logunit,*) &
             '[WARN] Cannot release GRIB message from ' // &
             trim(fg_filename)
        if (LIS_masterproc) then
           message(:) = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
           message(3) = &
                '  Cannot release GRIB message from ' // &
                trim(fg_filename)
           call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        rc = 1
        exit
     end if

     ! Break out of loop if we have what we need.
     if (count_swdown > 0 .and. &
         count_lwdown > 0) exit
  enddo

  call grib_close_file(ftn)

  deallocate ( dum1d )

!------------------------------------------------------------------
! See if we have everything. if not, abort.
!------------------------------------------------------------------

  if (count_swdown == 0) then
     write(LIS_logunit,*)'[WARN] Missing downward SW radiation in ', &
          trim(fg_filename)
     if (LIS_masterproc) then
        message(:) = ''
        message(1) = '[WARN] Program: LIS'
        message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
        message(3) = &
             '  Missing downward SW radiation in  ' // &
             trim(fg_filename)
        call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
             alert_number, message)
        alert_number = alert_number + 1
     end if
     rc = 1
  end if
  if (count_lwdown == 0) then
     write(LIS_logunit,*)'[WARN] Missing downward LW radiation in ', &
          trim(fg_filename)
     if (LIS_masterproc) then
        message(:) = ''
        message(1) = '[WARN] Program: LIS'
        message(2) = '  Routine USAF_fldbld_read_radflux_galwem.'
        message(3) = &
             '  Missing downward LW radiation in  ' // &
             trim(fg_filename)
        call LIS_alert( 'USAF_fldbld_read_radflux_galwem.', &
             alert_number, message)
        alert_number = alert_number + 1
     end if
     rc = 1
  end if

#endif

end subroutine USAF_fldbld_read_radflux_galwem
