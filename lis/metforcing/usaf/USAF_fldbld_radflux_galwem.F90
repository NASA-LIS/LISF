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
!subroutine USAF_fldbld_radflux_galwem(n,julhr,fc_hr,fg_swdata,fg_lwdata,rc)
subroutine USAF_fldbld_radflux_galwem(n,julhr,fg_swdata,fg_lwdata,rc)

! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_masterproc
  use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_alert, &
                                LIS_verify, LIS_endrun
  use LIS_timeMgrMod,    only : LIS_julhr_date
  use AGRMET_forcingMod, only : agrmet_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer             :: julhr
!  integer             :: fc_hr
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
!  \item[found2]
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
  integer                 :: ftn, igrib
  character*250           :: avnfile, avnfile2
  integer                 :: yr1, mo1, da1, hr1
  integer                 :: fc_hr
  character*255           :: message     ( 20 )
  integer                 :: iginfo      ( 2 )
  real                    :: gridres
  integer                 :: alert_number
  real, allocatable       :: fg_swdown    ( : , : )
  real, allocatable       :: fg_swdown1   ( : , : )
  real, allocatable       :: fg_lwdown    ( : , : )
  real, allocatable       :: fg_lwdown1   ( : , : )
  integer                 :: ifguess, jfguess
  integer                 :: center
  integer                 :: ierr
  logical*1               :: found, found2
  logical                 :: first_time
  integer                 :: yr_2d
  integer                 :: file_julhr
!  integer                 :: getsixhr
  integer                 :: dataDate, dataTime
  character*100           :: gtype
  logical                 :: found_inq
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
        !if (fc_hr > 24) exit ! Give up
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
        cycle
     end if

#if (defined USE_GRIBAPI)

     ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
     ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
     ! writing error messages to stdout/stderr, which may lead to runtime
     ! problems.

     if (found_inq) then
        call grib_open_file(ftn,trim(avnfile),'r',ierr)
     else
        ierr = 1
     end if
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) '[WARN] (1) Failed to open first guess - ',trim(avnfile)
     else
   !    ------------------------------------------------------------------
   !    read in the first grib record, unpack the header and extract
   !    section 1 and section 2 information.
   !    ------------------------------------------------------------------
        call grib_new_from_file(ftn,igrib,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) '[WARN] (1) Failed file read check - '//trim(avnfile)
        endif

        call grib_get(igrib,'centre',center,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grib_get: centre in ' // &
                                'USAF_fldbld_radflux_galwem'
        endif

        call grib_get(igrib,'gridType',gtype,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: gridtype in ' // &
                                'USAF_fldbld_radflux_galwem'
        endif

        call grib_get(igrib,'Ni',iginfo(1),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:Ni in ' // &
                                'USAF_fldbld_radflux_galwem'
        endif

        call grib_get(igrib,'Nj',iginfo(2),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:Nj in ' // &
                                'USAF_fldbld_radflux_galwem'
        endif

        call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) &
              'error in grid_get:jDirectionIncrementInDegrees in ' // &
              'USAF_fldbld_radflux_galwem'
        endif

        call grib_get(igrib,'dataDate',dataDate,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:dataDate in ' // &
                                'USAF_fldbld_radflux_galwem'
        endif

        call grib_get(igrib,'dataTime',dataTime,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:dataTime in ' // &
                                'USAF_fldbld_radflux_galwem'
        endif

        if ( yr1*10000+mo1*100+da1 == dataDate .and. hr1*100 == dataTime ) then
           found = .TRUE.
           if ( gtype /= "regular_ll" ) then
              message(1) = 'program: LIS'
              message(2) = '  Subroutine: USAF_fldbld_radflux_galwem'
              message(3) = '  First guess source is not a lat/lon grid'
              message(4) = '  USAF_fldbld_radflux_galwem expects lat/lon data'
              call lis_abort(message)
           endif
        endif

        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
     endif
#else
     write(LIS_logunit,*) '[ERR]: USAF_fldbld_radflux_galwem requires GRIB-API'
     write(LIS_logunit,*) '[ERR]: please recompile LIS'
     call LIS_endrun
#endif

     ! At this point, we have everything we need.
     found = .true.
     if (found) exit

  enddo ! Loop through cycles and forecast hours

  ! Give up if no acceptable GALWEM file was found
  if (.not. found) then
     write(LIS_logunit,*)'[WARN] No matching GALWEM Radiation file found!'
     rc = 1
     return
  end if

  write(LIS_logunit,*) &
       '[INFO] Using NWP Radiation fields from ',trim(avnfile)
  rc = 0


  if( found ) then
     write(LIS_logunit,*)'- FIRST GUESS DATA IS ON A ', gridres,&
          ' DEGREE LAT/LON GRID'
     ifguess = iginfo(1)
     jfguess = iginfo(2)

     if (center .eq. 57) then
        write(LIS_logunit,*)'- FIRST GUESS DATA IS FROM UK UM (GALWEM) MODEL'
     else
        write(LIS_logunit,*)'- UNKNOWN SOURCE FOR FIRST GUESS DATA'
     end if

!    ------------------------------------------------------------------
!      allocate first guess grid-specific variables.
!    ------------------------------------------------------------------
     allocate ( fg_swdown1 (ifguess, jfguess) )
     allocate ( fg_lwdown1 (ifguess, jfguess) )

!    ------------------------------------------------------------------
!      read in first guess data for this julian hour.
!    ------------------------------------------------------------------
     alert_number = 0

     call USAF_fldbld_read_radflux_galwem(avnfile, ifguess, jfguess,&
                                           fg_swdown1, fg_lwdown1, alert_number)

     !write(LIS_logunit,*)'EMK: maxval(fg_swdown1) = ', maxval(fg_swdown1)
     !write(LIS_logunit,*)'EMK: maxval(fg_lwdown1) = ', maxval(fg_lwdown1)

     !allocate ( fg_swdown (ifguess, jfguess) )
     !allocate ( fg_lwdown (ifguess, jfguess) )

! -----------------------------------------------------------------------
!    sometimes the subtraction of 3 hour radiation from 6 hour rad fluxes causes
!    slightly negative radiation values, set back to 0
! -----------------------------------------------------------------------
     !where (fg_swdown .lt. 0)
     !   fg_swdown=0
     !endwhere
     !where (fg_lwdown .lt. 0)
     !   fg_lwdown=0
     !endwhere
     !fg_swdown=0
     !fg_lwdown=0

     ! Spatially interpolate native grid to target domain grid:
     !call interp_galwem_first_guess(n, ifguess, jfguess, .true., &
     !                               fg_swdown, fg_swdata)
     !call interp_galwem_first_guess(n, ifguess, jfguess, .true., &
     !                               fg_lwdown, fg_lwdata)
     !deallocate ( fg_swdown, fg_swdown1 )
     !deallocate ( fg_lwdown, fg_lwdown1 )

     call interp_galwem_first_guess(n, ifguess, jfguess, .true., &
                                    fg_swdown1, fg_swdata)
     call interp_galwem_first_guess(n, ifguess, jfguess, .true., &
                                    fg_lwdown1, fg_lwdata)
     deallocate ( fg_swdown1 )
     deallocate ( fg_lwdown1 )

  else
        write(LIS_logunit,*)'[WARN] ** GALWEM RADIATION FLUX data not available **'

        message(1) = 'Program:  LIS'
        message(2) = '  Routine: USAF_fldbld_radflux_galwem.'
        message(3) = '  galwem radiation flux data not available, ' // &
                     'possible degradation.'
        alert_number = alert_number + 1
        if(LIS_masterproc) then
           call LIS_alert( 'fldbld              ', alert_number, message )
        endif
  endif

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
subroutine USAF_fldbld_read_radflux_galwem(fg_filename, ifguess, jfguess,&
                                            fg_swdown, fg_lwdown, alert_number )
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
  integer,        intent(inout) :: alert_number
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
  character*9                   :: cstat
  character*255                 :: message  ( 20 )
  integer                       :: count_swdown
  integer                       :: count_lwdown
  integer                       :: i
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: j
  integer                       :: k, c, r
  integer                       :: ftn, igrib, nvars
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, forecasttime_val
  real,           allocatable   :: dum1d    ( : )
  logical                       :: found_inq

! ------------------------------------------------------------------

! ------------------------------------------------------------------
!   read in grib file.
! ------------------------------------------------------------------

  ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
  ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
  ! writing error messages to stdout/stderr, which may lead to runtime
  ! problems.

  inquire(file=trim(fg_filename),exist=found_inq)
  if (.not. found_inq) then
     ierr = 1
  else
     ierr = 0
  end if
  call LIS_verify(ierr,'[WARN] FILE NOT FOUND - '//trim(fg_filename))

#if (defined USE_GRIBAPI)
  call grib_open_file(ftn,trim(fg_filename),'r',ierr)
  call LIS_verify(ierr,'[WARN] (3) Failed to open in read routine - '//trim(fg_filename))

  if ( ierr == 0 ) then
     allocate ( dum1d   (ifguess*jfguess) )
     count_swdown = 0
     count_lwdown = 0

     write(LIS_logunit,*)' '
     write(LIS_logunit,*)'[INFO] Reading first guess GALWEM radiation fluxes '
     write(LIS_logunit,*) trim(fg_filename)

     call grib_count_in_file(ftn,nvars,ierr)
     call LIS_verify(ierr, 'error in grib_count_in_file in ' // &
                           'USAF_fldbld_read_radflux_galwem')

     do k = 1, nvars

        call grib_new_from_file(ftn,igrib,ierr)
        call LIS_verify(ierr, 'failed to read - '//trim(fg_filename))

        call grib_get(igrib,'discipline',param_disc_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                              'USAF_fldbld_read_radflux_galwem')

        call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterCategory in ' // &
                              'USAF_fldbld_read_radflux_galwem')

        call grib_get(igrib,'parameterNumber',param_num_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                              'USAF_fldbld_read_radflux_galwem')
        
        call grib_get(igrib,'forecastTime',forecasttime_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: forecastTime in ' // &
                              'USAF_fldbld_read_radflux_galwem')

        !  SW Down Radiation flux:
        if ( param_disc_val == 0 .and. param_cat_val    == 4 .and. &
             param_num_val  == 7 ) then
!             param_num_val  == 7 .and. forecasttime_val == 0 ) then

           call grib_get(igrib,'values',dum1d,ierr)
           call LIS_verify(ierr, 'error in grib_get: SWdown values in ' // &
                                 'USAF_fldbld_read_radflux_galwem')
           if ( ierr == 0 ) then 
!     ------------------------------------------------------------------
!       store SW radiation flux to the fg_swdown array.
!     ------------------------------------------------------------------
              fg_swdown = reshape(dum1d, (/ifguess,jfguess/))
              count_swdown = count_swdown + 1
              ! found swdown; clean up and break out of do-loop.
!              call grib_release(igrib,ierr)
!              exit
           else
              write(cstat,'(i9)',iostat=istat1) ierr
              message(1) = 'Program: LIS'
              message(2) = '  Subroutine:  AGRMET_fldbld_read.'
              message(3) = '  Error reading first guess file:'
              message(4) = '  ' // trim(fg_filename)
              if( istat1 .eq. 0 )then
                 message(5) = '  Status = ' // trim(cstat)
              endif
              if ( allocated(dum1d) )   deallocate(dum1d)
              call LIS_abort( message)
           endif
        endif ! SWdown

        !  LW Down Radiation flux:
        if ( param_disc_val == 0 .and. param_cat_val    == 5 .and. &
             param_num_val  == 3 ) then
!             param_num_val  == 3 .and. forecasttime_val == 0 ) then

           call grib_get(igrib,'values',dum1d,ierr)
           call LIS_verify(ierr, 'error in grib_get: LWdown values in ' // &
                                 'USAF_fldbld_read_radflux_galwem')
           if ( ierr == 0 ) then
!     ------------------------------------------------------------------
!       store LW radiation flux to the fg_lwdown array.
!     ------------------------------------------------------------------
              fg_lwdown = reshape(dum1d, (/ifguess,jfguess/))
              count_lwdown = count_lwdown + 1
           else
              write(cstat,'(i9)',iostat=istat1) ierr
              message(1) = 'Program: LIS'
              message(2) = '  Subroutine:  USAF_fldbld_read_radflux_galwem.'
              message(3) = '  Error reading first guess file:'
              message(4) = '  ' // trim(fg_filename)
              if( istat1 .eq. 0 )then
                 message(5) = '  Status = ' // trim(cstat)
              endif
              if ( allocated(dum1d) )   deallocate(dum1d)
              call LIS_abort( message)
           endif
        endif ! LWdown

        call grib_release(igrib,ierr)
     enddo

     call grib_close_file(ftn)

     deallocate ( dum1d )

!     ------------------------------------------------------------------
!     see if we have everything. if not, abort.
!     ------------------------------------------------------------------

     if ( count_swdown == 0 .or. count_lwdown == 0 ) then
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  USAF_fldbld_read_radflux_galwem.'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(fg_filename)
        message(5) = '  File does not contain all data that fldbld needs.'
        if ( allocated(dum1d) ) deallocate(dum1d)
        call LIS_abort( message)
     endif

  else
     write(cstat,'(i9)',iostat=istat1) ierr
     message(1) = 'Program: LIS'
     message(2) = '  Subroutine:  USAF_fldbld_read_radflux_galwem.'
     message(3) = '  Error opening first guess file:'
     message(4) = '  ' // trim(fg_filename)
     if ( istat1 == 0 ) then
        message(5) = '  Status = ' // trim(cstat)
     endif
     call LIS_abort(message)
  endif
#endif

end subroutine USAF_fldbld_read_radflux_galwem
