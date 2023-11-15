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
! !ROUTINE: AGRMET_fldbld_precip_gfs
! \label{AGRMET_fldbld_precip_gfs}
!
! !REVISION HISTORY:
! 05 May 2013  Initial version, forked from AGRMET_fldbld
!              ...................................Ryan Ruhge/16WS/WXE/SEMS
! 30 Jan 2014  Changed the grib reader to grib api library
!              ...................................James Geiger/NASA
! 11 Aug 2016  Renamed to AGRMET_fldbld_precip_gfs
!              ...................................James Geiger/NASA
!
! !INTERFACE:
subroutine AGRMET_fldbld_precip_gfs(n,findex,julhr,fc_hr,gfsdata)
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
  integer, intent(in) :: findex
  integer             :: fc_hr
  real, intent(out)   :: gfsdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
!
! !DESCRIPTION:
!  This routine interpolates the precipitation
!  data to the AGRMET grid.
!
! The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \item[avnfile]
!    name of the 3 hour forecast file
!  \item[avnfile2]
!   name of the 6 hour forecast file
!  \item[nunit]
!    unix file descriptor
!  \item[nunit2]
!    unix file descriptor
!  \item[ksec2]
!   array of section 2 information from a grib record
!   see comments for more information.
!  \item[message]
!    error message
!  \item[iginfo]
!   array of grid information from a grib record
!  \item[ginfo]
!   array of grid information from a grib record
!  \item[gridres]
!    the resolution, in degrees, of the first
!    guess grid
!  \item[fg\_prec]
!   precipitation data to be interpolated to lis grid
!  \item[fg\_prec1]
!   3 hour forecast precipitation data
!  \item[fg\_prec2]
!   6 hour forecast precipitation data
!  \item[alert\_number]
!    number of alerts that occur in the program
!  \item[ifguess]
!    east/west dimension of first guess grid
!  \item[jfguess]
!    north/south dimension of first guess grid
!  \item[gfsdata]
!   array of gfs rain rate data
!  \item[istat]
!    io error stat variable
!  \item[nmlfil]
!    name of first guess namelist file
!  \item[ksec1]
!   array of grib section 1 information
!   see comments for more specific information
!  \item[avnflag]
!    flag used in grid calculations to account for the
!    fact that point (1,1) on the avn/nogaps
!    grid is at the north/south pole.
!  \item[center]
!    meteorological center that produced the first
!    guess data (7-NCEP, 58-FNMOC
!  \item[ierr]
!    error code
!  \item[fc\_hr]
!    forecast hour or the difference between reference and valid time
!  \item[yr\_2d]
!    2 digit year for comparison with GRIB header
!  \item[found]
!    logical flag set true when an input file with the correct valid time is found
!  \item[found2]
!    logical flag set true when an input file with the correct valid time is found
!  \item[file\_julhr]
!    julian hour used to determine names of forecast files from previous cycles
!  \item[getsixhr]
!   indicates whether to get data from the 6 hour forecast file
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!    converts the julian hour to a date format
!  \item[getAVNfilename](\ref{getAVNfilename}) \newline
!    generates the first guess AVN filename
!  \item[AGRMET\_fldbld\_read\_precip\_gfs]
!  (\ref{AGRMET_fldbld_read_precip_gfs}) \newline
!   read AVN or NOGAPS data in grib format
!  \item[AGRMET\_fg2lis\_precip](\ref{AGRMET_fg2lis_precip}) \newline
!   interpolate first guess data to the AGRMET grid
!  \item[LIS\_abort](\ref{LIS_abort}) \newline
!   abort in case of error
!  \end{description}
!EOP
  integer                 :: ftn, igrib
  !character*120           :: avnfile, avnfile2
  character*255           :: avnfile, avnfile2
  integer                 :: yr1, mo1, da1, hr1
  integer                 :: julhr
  integer                 :: nunit, nunit2
  integer                 :: ksec2       ( 10 )
  character*255           :: message     ( 20 )
  integer                 :: iginfo      ( 40 )
  real                    :: ginfo       ( 40 )
  real                    :: gridres
  integer                 :: alert_number
  real, allocatable       :: fg_prec    ( : , : )
  real, allocatable       :: fg_prec1   ( : , : )
  real, allocatable       :: fg_prec2   ( : , : )
  integer                 :: ifguess, jfguess
  integer                 :: ksec1       ( 100 )
  integer                 :: avnflag
  integer                 :: center
  integer                 :: ierr
  logical*1               :: found, found2
  integer                 :: yr_2d
  integer                 :: file_julhr
  integer                 :: getsixhr
  integer                 :: dataDate, dataTime
  character*100           :: gtype
  logical :: found_inq

  alert_number = 0 ! EMK BUG FIX

  found = .false.
  found2 = .false.
  call LIS_julhr_date(julhr,yr1,mo1,da1,hr1)
  file_julhr = julhr

!     ------------------------------------------------------------------
!     Need to process the current and previous 6 hour instances
!     Search for an analysis or forecast file for upto 24 hours with
!     the needed valid time
!     ------------------------------------------------------------------

  do while( ((.not.found) .or. (.not.found2)) .and. (fc_hr <= 12))
     found = .false.

!    --------------------------------------------------------------------
!    determine if we need the 6 hour forecast file
!    --------------------------------------------------------------------
     if (mod(fc_hr,6) .eq. 0) then
        getsixhr=1
        found2=.false.
     else
        getsixhr=0
        found2=.true.
     endif

     yr_2d = mod(yr1,100)
     if(yr_2d.eq.0) yr_2d = 100
     ! EMK...Added support for new GFS filename convention
     call getAVNfilename(avnfile, agrmet_struc(n)%agrmetdir,&
          agrmet_struc(n)%gfsdir, agrmet_struc(n)%use_timestamp,&
          agrmet_struc(n)%gfs_timestamp, &
          agrmet_struc(n)%gfs_filename_version, &
          yr1, mo1, da1, hr1, fc_hr)
     if (getsixhr.eq.1) then
        call getAVNfilename(avnfile2, agrmet_struc(n)%agrmetdir,&
             agrmet_struc(n)%gfsdir, agrmet_struc(n)%use_timestamp,&
             agrmet_struc(n)%gfs_timestamp, &
             agrmet_struc(n)%gfs_filename_version, &
             yr1, mo1, da1, hr1, fc_hr-3)
     endif

!     ------------------------------------------------------------------
!     open first guess grib data using library utility.  just read
!     the first file only, as all data will be of the same type
!     (avn or nogaps) because the search script ensures that it is.
!     ------------------------------------------------------------------

#if (defined USE_GRIBAPI)

     ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
     ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
     ! writing error messages to stdout/stderr, which may lead to runtime
     ! problems.
     !call grib_open_file(ftn,trim(avnfile),'r',ierr)
     inquire(file=trim(avnfile),exist=found_inq)
     if (found_inq) then
        call grib_open_file(ftn,trim(avnfile),'r',ierr)
     else
        ierr = 1
     end if
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'Failed to open - ', trim(avnfile)
     else
   !    ------------------------------------------------------------------
   !    read in the first grib record, unpack the header and extract
   !    section 1 and section 2 information.
   !    ------------------------------------------------------------------
        call grib_new_from_file(ftn,igrib,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'failed to read - '//trim(avnfile)
        endif

        call grib_get(igrib,'centre',center,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grib_get: centre in ' // &
                                'AGRMET_fldbld_precip_gfs'
        endif

        call grib_get(igrib,'gridType',gtype,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: gridtype in ' // &
                                'AGRMET_fldbld_precip_gfs'
        endif

        call grib_get(igrib,'Ni',iginfo(1),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:Ni in ' // &
                                'AGRMET_fldbld_precip_gfs'
        endif

        call grib_get(igrib,'Nj',iginfo(2),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:Nj in ' // &
                                'AGRMET_fldbld_precip_gfs'
        endif

        call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) &
              'error in grid_get:jDirectionIncrementInDegrees in ' // &
              'AGRMET_fldbld_precip_gfs'
        endif

        call grib_get(igrib,'dataDate',dataDate,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:dataDate in ' // &
                                'AGRMET_fldbld_precip_gfs'
        endif

        call grib_get(igrib,'dataTime',dataTime,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:dataTime in ' // &
                                'AGRMET_fldbld_precip_gfs'
        endif

        if ( yr1*10000+mo1*100+da1 == dataDate .and. hr1*100 == dataTime ) then
           found = .TRUE.

   !    ------------------------------------------------------------------
   !    Ensure data is on a lat/lon grid by ensuring section 2 octet 6
   !    is a "0".  This is stored in array ksec2(4).  FLDBLD can not
   !    handle any other grid type, so abort otherwise.
   !    See GRIB utility and GRIB manual for more information.
   !    ------------------------------------------------------------------

           if ( gtype /= "regular_ll" ) then
              message(1) = 'program: LIS'
              message(2) = '  Subroutine: agrmet_sfcalc'
              message(3) = '  First guess source is not a lat/lon grid'
              message(4) = '  agrmet_sfcalc expects lat/lon data'
              call lis_abort(message)
           endif
        endif

        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
     endif
#else
  write(LIS_logunit,*) 'ERR: AGRMET_fldbld_precip_gfs requires GRIB-API'
  write(LIS_logunit,*) 'ERR: please recompile LIS'
  call LIS_endrun
#endif

!    -------------------------------------------------------------------
!    if we need to get the six hour forecast file, repeat steps above
!    for that file
!    -------------------------------------------------------------------
     if (getsixhr.eq.1) then
#if (defined USE_GRIBAPI)

        ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
        ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
        ! writing error messages to stdout/stderr, which may lead to runtime
        ! problems.
        !call grib_open_file(ftn,trim(avnfile),'r',ierr)
        inquire(file=trim(avnfile2),exist=found_inq)
        if (found_inq) then
           call grib_open_file(ftn,trim(avnfile2),'r',ierr)
        else
           ierr = 1
        end if
        if(ierr.ne.0) then
           write(LIS_logunit,*) 'Failed to open - ', trim(avnfile2)
        else
   !     ------------------------------------------------------------------
   !     read in the first grib record, unpack the header and extract
   !     section 1 and section 2 information.
   !     ------------------------------------------------------------------
           call grib_new_from_file(ftn,igrib,ierr)
           call LIS_verify(ierr, 'failed to read - '//trim(avnfile2))

           call grib_get(igrib,'centre',center,ierr)
           call LIS_verify(ierr, 'error in grib_get: centre in ' // &
                                 'AGRMET_fldbld_precip_gfs')

           call grib_get(igrib,'gridType',gtype,ierr)
           call LIS_verify(ierr, 'error in grid_get: gridtype in ' // &
                                 'AGRMET_fldbld_precip_gfs')

           call grib_get(igrib,'Ni',iginfo(1),ierr)
           call LIS_verify(ierr, 'error in grid_get:Ni in ' // &
                                 'AGRMET_fldbld_precip_gfs')

           call grib_get(igrib,'Nj',iginfo(2),ierr)
           call LIS_verify(ierr, 'error in grid_get:Nj in ' // &
                                 'AGRMET_fldbld_precip_gfs')

           call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
           call LIS_verify(ierr, &
                   'error in grid_get:jDirectionIncrementInDegrees in ' // &
                   'AGRMET_fldbld_precip_gfs')

           call grib_get(igrib,'dataDate',dataDate,ierr)
           call LIS_verify(ierr, 'error in grid_get:dataDate in ' // &
                                 'AGRMET_fldbld_precip_gfs')

           call grib_get(igrib,'dataTime',dataTime,ierr)
           call LIS_verify(ierr, 'error in grid_get:dataTime in ' // &
                                 'AGRMET_fldbld_precip_gfs')

           if ( yr1*10000+mo1*100+da1 == dataDate .and. &
                hr1*100 == dataTime ) then
              found2 = .TRUE.

   !    ------------------------------------------------------------------
   !    Ensure data is on a lat/lon grid by ensuring section 2 octet 6
   !    is a "0".  This is stored in array ksec2(4).  FLDBLD can not
   !    handle any other grid type, so abort otherwise.
   !    See GRIB utility and GRIB manual for more information.
   !    ------------------------------------------------------------------
              if ( gtype /= "regular_ll" ) then
                 message(1) = 'program: LIS'
                 message(2) = '  Subroutine: agrmet_sfcalc'
                 message(3) = '  First guess source is not a lat/lon grid'
                 message(4) = '  agrmet_sfcalc expects lat/lon data'
                 call lis_abort(message)
              endif
           endif

           call grib_release(igrib,ierr)
           call grib_close_file(ftn)
        endif
#endif
     endif


!    ------------------------------------------------------------------
!    If the correct valid time is not found:
!       Increment forecast hour by 6.
!       Decrement file_julhr by 6 and get the new filename elements.
!    ------------------------------------------------------------------

     if ((.not. found).or.(.not. found2)) then
        fc_hr = fc_hr + 6
        file_julhr = file_julhr - 6
        call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)
     endif
  enddo


!     ------------------------------------------------------------------
!     determine resolution of grid.  this is in section 2 octet 24-25
!     and is stored in millidegrees (iginfo(8)).  also get the
!     i and j dimensions of the grid, section 2, octets 7-8 and 9-10
!     (iginfo(1 and 2).  see grib utility and grib manual for
!     more details.
!     ------------------------------------------------------------------
   if ((found) .and. (found2)) then
     write(LIS_logunit,*)'- FIRST GUESS DATA IS ON A ', gridres,&
          ' DEGREE LAT/LON GRID'
     ifguess = iginfo(1)
     jfguess = iginfo(2)
!     ------------------------------------------------------------------
!     determine starting point (j coordinate) of first guess grid.
!     for avn, j=1 at the north pole, for nogaps, j=1 at the south
!     pole.  this is stored in section 2, octets 11-13 in millidegrees
!     (stored in iginfo(3)). this information is needed by subroutine
!     setup which calculates the location of the agrmet grid with respect
!     to the first guess grid.  see grib utility and grib manual for
!     more details.
!     ------------------------------------------------------------------
! Not used.
!     if ( iginfo(3) .eq. -90000 ) then
!        avnflag =  1
!     else if ( iginfo(3) .eq. 90000 ) then
!        avnflag = -1
!     end if
!     ------------------------------------------------------------------
!     get modeling center (58 is FNMOC, 7 is NCEP), section 1 octet 5.
!     our degribber stores this in ksec1(3).  see grib utility and
!     grib manual for more details.
!     ------------------------------------------------------------------

     if (center .eq. 58) then
        write(LIS_logunit,*)'- FIRST GUESS DATA IS FROM NOGAPS MODEL'
     elseif (center .eq. 7) then
        write(LIS_logunit,*)'- FIRST GUESS DATA IS FROM GFS MODEL'
     end if


!     ------------------------------------------------------------------
!       allocate first guess grid-specific variables.
!     ------------------------------------------------------------------

     allocate ( fg_prec1 (ifguess, jfguess) )
     if (getsixhr.eq.1) allocate ( fg_prec2 (ifguess, jfguess) )


!     ------------------------------------------------------------------
!         read in first guess data for this julian hour.
!     ------------------------------------------------------------------
     alert_number = 0

     call AGRMET_fldbld_read_precip_gfs( avnfile, ifguess, jfguess,&
          fg_prec1, alert_number )
     if (getsixhr.eq.1) &
          call AGRMET_fldbld_read_precip_gfs( avnfile2, ifguess, jfguess,&
          fg_prec2, alert_number )

     allocate ( fg_prec (ifguess, jfguess) )

! -----------------------------------------------------------------------
!    if we need 6 hour precip data, subtract from 3 hour precip data
!    since the data are originally for the 0-6 hour forecast.
! -----------------------------------------------------------------------
     if (getsixhr.eq.1) then
        fg_prec = fg_prec1 - fg_prec2
     else
        fg_prec = fg_prec1
     endif

! -----------------------------------------------------------------------
!    sometimes the subtraction of 3 hour precip from 6 hour precip causes
!    slightly negative precipitation values, set back to 0
! -----------------------------------------------------------------------
     where (fg_prec .lt. 0)
        fg_prec=0
     endwhere

     call AGRMET_fg2lis_precip(n,findex, ifguess, jfguess, &
          fg_prec, gfsdata)

     deallocate ( fg_prec )
     deallocate ( fg_prec1 )
     if (getsixhr.eq.1) deallocate ( fg_prec2 )
  else
     write(LIS_logunit,*)'- ** GFS Precipitation data not available **'
     flush(LIS_logunit)

        message(1) = 'Program:  LIS'
        message(2) = '  Routine:  fldbld_precip.'
        message(3) = '  GFS Precipitation data not available! '
        alert_number = alert_number + 1
        if(LIS_masterproc) then
           call LIS_alert( 'fldbld              ', alert_number, message )
        endif
        call LIS_abort( message) ! EMK...Cannot continue w/o background field!

  endif
end subroutine AGRMET_fldbld_precip_gfs


!BOP
!
! !ROUTINE: AGRMET_fldbld_read_precip_gfs
!  \label{AGRMET_fldbld_read_precip_gfs}
!
! !REVISION HISTORY:
! 05 May 2013  Initial version, forked from AGRMET_fldbld_read
!              ...................................Ryan Ruhge/16WS/WXE/SEMS
! 30 Jan 2014  Changed the grib reader to grib api library
!              ...................................James Geiger/NASA
!
! !INTERFACE:
subroutine AGRMET_fldbld_read_precip_gfs( fg_filename, ifguess, jfguess,&
     fg_prec,&
     alert_number )
! !USES:
  use LIS_coreMod, only : LIS_masterproc
  use LIS_logMod, only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  character(len=*),  intent(in) :: fg_filename
  integer,        intent(in)    :: ifguess
  integer,        intent(in)    :: jfguess
  real,           intent(out)   :: fg_prec     ( ifguess,jfguess )
  integer,        intent(inout) :: alert_number
!
! !DESCRIPTION:
!
!     to read gfs or navgem data in grib format.
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
!      \item[fg\_prec]
!        total precipitation read in from gfs file
!       \item[alert\_number]
!        counts number of alert messages sent
!       \item[cstat]
!        I/O status, character
!       \item[message]
!        Error message
!       \item[count\_prec]
!        counts number of precipitation
!        levels read in from first guess file:q
!       \item[file\_age]
!        stores forecast hour (0 for analysis) of first
!        guess data
!       \item[i,j,k]
!        looping and indexing variables
!       \item[ierr,istat1]
!        error status
!       \item[ksec1]
!        array of grib section 1 information
!       \item[nunit]
!        unix file descriptor
!       \item[dum1d]
!        dummy array
!       \end{description}
!
!EOP
  character*9                   :: cstat
  character*255                 :: message     ( 20 )
  integer                       :: count_prec
  integer                       :: file_age
  integer                       :: i
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: j
  integer                       :: k, c, r
  integer                       :: ksec1       ( 100 )
  integer                       :: nunit
  integer                       :: ftn, igrib, nvars
  integer                       :: pds7_val, pds8_val, pds9_val
  ! EMK...For GRIB 2
  integer :: editionNumber_val
  integer :: param_disc_val, param_cat_val, &
       param_num_val
  real,           allocatable   :: dum1d       ( : )
  logical :: found_inq

!     ------------------------------------------------------------------
!     executable code begins here ...
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     read in grib file.
!     ------------------------------------------------------------------


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
  call LIS_verify(ierr,'Failed to open - '//trim(fg_filename))

#if (defined USE_GRIBAPI)
  call grib_open_file(ftn,trim(fg_filename),'r',ierr)
  call LIS_verify(ierr,'Failed to open - '//trim(fg_filename))

  if ( ierr .eq. 0 ) then
     allocate ( dum1d   (ifguess*jfguess) )
     count_prec = 0

     write(LIS_logunit,*)'[INFO] Reading first guess precip ', trim(fg_filename)

     call grib_count_in_file(ftn,nvars,ierr)
     call LIS_verify(ierr, 'error in grib_count_in_file in AGRMET_fldbld_read')

     do k = 1, nvars

!     ------------------------------------------------------------------
!       read a single grib record from unix file descriptor nunit,
!       unpack header data, and place data into dum1d array.
!
!       nunit - is the unix file descriptor passed back from call
!               to copen.
!       dum1d - contains the degribbed data.
!       ierr  - is an error stat.  0 - no error.   1 - end of file
!                                  2 - read error      indicator.
!
!       note: parameters are identified by grib section 1 information
!
!       octet 9  - ksec1(7) - variable type
!
!       octet 10 - ksec1(8) - level type, ex. isobaric is 100,
!                                             AGL      is 105.
!
!       octet 11 - ksec1(9) - height, pressure, etc. of level.
!
!       see grib utility code and grib manual for more details.
!     -------------------------------------------------------------------

        call grib_new_from_file(ftn,igrib,ierr)
        call LIS_verify(ierr, 'failed to read - '//trim(fg_filename))

        ! EMK...Allow for GRIB2
        call grib_get(igrib,'editionNumber',editionNumber_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: editionNumber in ' //&
                              'AGRMET_fldbld_read_precip_gfs')

        if (editionNumber_val .eq. 1) then
           call grib_get(igrib,'indicatorOfParameter',pds7_val,ierr)
           call LIS_verify(ierr, &
                'error in grib_get: indicatorOfParameter in ' //&
                'AGRMET_fldbld_read_precip_gfs')

           call grib_get(igrib,'level',pds9_val,ierr)
           call LIS_verify(ierr, 'error in grib_get: level in ' // &
                'AGRMET_fldbld_read_precip_gfs')

           call grib_get(igrib,'indicatorOfTypeOfLevel',pds8_val,ierr)
           call LIS_verify(ierr, 'error in grib_get: level in ' // &
                'AGRMET_fldbld_read_precip_gfs')

           call grib_get(igrib, 'values',dum1d,ierr)
           call LIS_verify(ierr, 'error in grib_get: values in ' // &
                'AGRMET_fldbld_read_precip_gfs')
           
        else ! GRIB 2
            call grib_get(igrib,'discipline',param_disc_val,ierr)
            call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                              'AGRMET_fldbld_read_precip_gfs')
            
            call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
            call LIS_verify(ierr, &
                 'error in grib_get: parameterCategory in ' // &
                 'AGRMET_fldbld_read_precip_gfs')

            call grib_get(igrib,'parameterNumber',param_num_val,ierr)
            call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                 'AGRMET_fldbld_read_precip_gfs')
            
            if ( param_disc_val == 0 .and. param_cat_val    == 1 .and. &
                 param_num_val  == 8 ) then
               call grib_get(igrib, 'values',dum1d,ierr)
               call LIS_verify(ierr, 'error in grib_get: values in ' // &
                    'AGRMET_fldbld_read_precip_gfs')
            end if

        end if

        if(ierr.eq.0) then
!     ------------------------------------------------------------------
!       store total precipitation to the fg_prec array.
!     ------------------------------------------------------------------

           if (editionNumber_val .eq. 1) then ! GRIB1
              if ( (pds7_val .eq. 61) .and. &
                   (pds8_val .eq. 1)  .and. &
                   (pds9_val .eq. 0)  ) then
                 do r=1,jfguess
                    do c=1,ifguess
                       fg_prec(c,r) = dum1d(c+(r-1)*ifguess)
                    enddo
                 enddo
                 count_prec = count_prec + 1
              end if
           else ! GRIB2
              if ( param_disc_val == 0 .and. param_cat_val    == 1 .and. &
                   param_num_val  == 8 ) then
                 do r=1,jfguess
                    do c=1,ifguess
                       fg_prec(c,r) = dum1d(c+(r-1)*ifguess)
                    enddo ! c
                 enddo ! r
                 count_prec = count_prec + 1
              end if
           end if ! editionNumber

        elseif(ierr.eq.1) then
           exit
        else
           write(cstat,'(i9)',iostat=istat1) ierr
           message(1) = 'Program: LIS'
           message(2) = '  Subroutine:  AGRMET_fldbld_read_precip_gfs'
           message(3) = '  Error reading first guess file:'
           message(4) = '  ' // trim(fg_filename)
           if( istat1 .eq. 0 )then
              message(5) = '  Status = ' // trim(cstat)
           endif
           if ( allocated(dum1d) )   deallocate(dum1d)
           call LIS_abort( message)
        endif
        call grib_release(igrib,ierr)
     enddo

     call grib_close_file(ftn)

     deallocate ( dum1d )

!     ------------------------------------------------------------------
!     check if first guess data is forecast data (more than six hours)
!     if forecast data more than six hours is used, send alert
!     message.  file age is determined from section 1, octets 19 and/or
!     20 depending on the time range indicator, octet 21.  if the
!     time range indicator (stored in ksec1(19)) is zero, then
!     the forecast hour is in octet 19 (stored in ksec1(17).  if
!     the time range indicator is 10, then the forecast hour is
!     gribbed across octets 19 and 20 and the degrib utility stores
!     it in ksec1(18).  see degrib utility and grib manual for
!     more details.
!
!     note, for nogaps, the 0 hour forecast is gribbed as a 1 hour
!     forecast (i don't know why, but this violates the grib standard).
!     so this is handled with an if statement below.
!     ------------------------------------------------------------------
!     jvg: disabling the check below since it only seems to be applicable
!     to nogaps
#if 0
     if (ksec1(19) .eq. 10) then
        file_age = ksec1(18)
     elseif (ksec1(19) .eq. 0) then
        file_age = ksec1(17)
     elseif (ksec1(19) .eq. 1) then
        file_age = 0
     end if

     if( file_age .gt. 6 )then

        write(LIS_logunit,*)'- ** FORECAST DATA BEING USED, POSSIBLE ' // &
                            'DEGRADATION **'

        message(1) = 'Program:  LIS'
        message(2) = '  Routine:  fg_read.'
        message(3) = '  Forecast data being used, possible degradation.'
        message(4) = '  First guess file used is:'
        message(5) = '  ' // trim(fg_filename)
        alert_number = alert_number + 1
        if(LIS_masterproc) then
           call LIS_alert( 'fldbld              ', alert_number, message )
        endif
     end if
#endif

!     ------------------------------------------------------------------
!     see if we have everything. if not, abort.
!     ------------------------------------------------------------------

     if (count_prec.eq. 0) then
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  AGRMET_fldbld_read_precip_gfs'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(fg_filename)
        message(5) = '  File does not contain all data that fldbld needs.'
        if ( allocated(dum1d) )   deallocate(dum1d)
        call LIS_abort( message)
     endif
  else
     write(cstat,'(i9)',iostat=istat1) ierr
     message(1) = 'Program: LIS'
     message(2) = '  Subroutine:  AGRMET_fldbld_read_precip_gfs'
     message(3) = '  Error opening first guess file:'
     message(4) = '  ' // trim(fg_filename)
     if( istat1 .eq. 0 )then
        message(5) = '  Status = ' // trim(cstat)
     endif
     call LIS_abort( message)
  endif
#endif
end subroutine AGRMET_fldbld_read_precip_gfs
