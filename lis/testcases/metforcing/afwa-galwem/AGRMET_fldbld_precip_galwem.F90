!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: AGRMET_fldbld_precip_galwem
! \label{AGRMET_fldbld_precip_galwem}
!
! !REVISION HISTORY:
! 11 Aug 2016  Initial specification based on AGRMET_fldbld_precip_gfs
!              ...........................................James Geiger/NASA
!
! !INTERFACE:
subroutine AGRMET_fldbld_precip_galwem(n,julhr,fc_hr,fg_data)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_masterproc
  use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_alert, &
                                LIS_verify, LIS_endrun
  use LIS_timeMgrMod,    only : LIS_julhr_date
  use AGRMET_forcingMod, only : agrmet_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif
!<debug -- jim testing>
  use LIS_mpiMod
  use LIS_historyMod
!</debug -- jim testing>

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer             :: fc_hr
  real, intent(out)   :: fg_data(LIS_rc%lnc(n), LIS_rc%lnr(n))
!
! !DESCRIPTION:
!  This routine interpolates the precipitation
!  data to the AGRMET grid.
!
! The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fc\_hr]
!    forecast hour or the difference between reference and valid time
!  \item[fg\_data]
!   array of galwem rain rate data
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
!   indicates whether to get data from the 6 hour forecast file
!  \item[dataDate]
!   date of values in the GRIB message
!  \item[dataTime]
!   time of values in the GRIB message
!  \item[gtype]
!   type of grid determined by querying GRIB message
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!    converts the julian hour to a date format
!  \item[AGRMET_getGALWEMfilename](\ref{AGRMET_getGALWEMfilename}) \newline
!    generates the first guess GALWEM filename
!  \item[AGRMET\_fldbld\_read\_precip\_galwem]
!   (\ref{AGRMET_fldbld_read_precip_galwem}) \newline
!   read GALWEM precipitation data in grib format
!  \item[interp\_galwem\_first\_guess](\ref{interp_galwem_first_guess}) \newline
!   interpolate first guess data to the AGRMET grid
!  \item[LIS\_abort](\ref{LIS_abort}) \newline
!   abort in case of error
!  \end{description}
!EOP
  integer                 :: ftn, igrib
  character*120           :: avnfile, avnfile2
  integer                 :: yr1, mo1, da1, hr1
  integer                 :: julhr
  character*100           :: message     ( 20 )
  integer                 :: iginfo      ( 2 )
  real                    :: gridres
  integer                 :: alert_number
  real, allocatable       :: fg_prec    ( : , : )
  real, allocatable       :: fg_prec1   ( : , : )
  real, allocatable       :: fg_prec2   ( : , : )
  integer                 :: ifguess, jfguess
  integer                 :: center
  integer                 :: ierr
  logical*1               :: found, found2
  integer                 :: yr_2d
  integer                 :: file_julhr
  integer                 :: getsixhr
  integer                 :: dataDate, dataTime
  character*100           :: gtype



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
     call AGRMET_getGALWEMfilename(avnfile, agrmet_struc(n)%agrmetdir,&
          agrmet_struc(n)%galwemdir,agrmet_struc(n)%use_timestamp,&
          yr1,mo1,da1,hr1,fc_hr)
     if (getsixhr.eq.1) then
        call AGRMET_getGALWEMfilename(avnfile2, agrmet_struc(n)%agrmetdir,&
             agrmet_struc(n)%galwemdir,agrmet_struc(n)%use_timestamp,&
             yr1,mo1,da1,hr1,fc_hr-3)
     endif
!     ------------------------------------------------------------------
!     open first guess grib data using library utility.  just read
!     the first file only, as all data will be of the same type
!     (avn or nogaps) because the search script ensures that it is.
!     ------------------------------------------------------------------

#if (defined USE_GRIBAPI)
     call grib_open_file(ftn,trim(avnfile),'r',ierr)
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
                                'AGRMET_fldbld_precip_galwem'
        endif

        call grib_get(igrib,'gridType',gtype,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: gridtype in ' // &
                                'AGRMET_fldbld_precip_galwem'
        endif

        call grib_get(igrib,'Ni',iginfo(1),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:Ni in ' // &
                                'AGRMET_fldbld_precip_galwem'
        endif

        call grib_get(igrib,'Nj',iginfo(2),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:Nj in ' // &
                                'AGRMET_fldbld_precip_galwem'
        endif

        call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) &
              'error in grid_get:jDirectionIncrementInDegrees in ' // &
              'AGRMET_fldbld_precip_galwem'
        endif

        call grib_get(igrib,'dataDate',dataDate,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:dataDate in ' // &
                                'AGRMET_fldbld_precip_galwem'
        endif

        call grib_get(igrib,'dataTime',dataTime,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get:dataTime in ' // &
                                'AGRMET_fldbld_precip_galwem'
        endif

        if ( yr1*10000+mo1*100+da1 == dataDate .and. hr1*100 == dataTime ) then
           found = .TRUE.
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
  write(LIS_logunit,*) 'ERR: AGRMET_fldbld_precip_galwem requires GRIB-API'
  write(LIS_logunit,*) 'ERR: please recompile LIS'
  call LIS_endrun
#endif

!    -------------------------------------------------------------------
!    if we need to get the six hour forecast file, repeat steps above
!    for that file
!    -------------------------------------------------------------------
     if (getsixhr.eq.1) then
#if (defined USE_GRIBAPI)
        call grib_open_file(ftn,trim(avnfile2),'r',ierr)
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
                                 'AGRMET_fldbld_precip_galwem')

           call grib_get(igrib,'gridType',gtype,ierr)
           call LIS_verify(ierr, 'error in grid_get: gridtype in ' // &
                                 'AGRMET_fldbld_precip_galwem')

           call grib_get(igrib,'Ni',iginfo(1),ierr)
           call LIS_verify(ierr, 'error in grid_get:Ni in ' // &
                                 'AGRMET_fldbld_precip_galwem')

           call grib_get(igrib,'Nj',iginfo(2),ierr)
           call LIS_verify(ierr, 'error in grid_get:Nj in ' // &
                                 'AGRMET_fldbld_precip_galwem')

           call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
           call LIS_verify(ierr, &
                   'error in grid_get:jDirectionIncrementInDegrees in ' // &
                   'AGRMET_fldbld_precip_galwem')

           call grib_get(igrib,'dataDate',dataDate,ierr)
           call LIS_verify(ierr, 'error in grid_get:dataDate in ' // &
                                 'AGRMET_fldbld_precip_galwem')

           call grib_get(igrib,'dataTime',dataTime,ierr)
           call LIS_verify(ierr, 'error in grid_get:dataTime in ' // &
                                 'AGRMET_fldbld_precip_galwem')

           if ( yr1*10000+mo1*100+da1 == dataDate .and. &
                hr1*100 == dataTime ) then
              found2 = .TRUE.
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


   if ((found) .and. (found2)) then
     write(LIS_logunit,*)'- FIRST GUESS DATA IS ON A ', gridres,&
          ' DEGREE LAT/LON GRID'
     ifguess = iginfo(1)
     jfguess = iginfo(2)

     if (center .eq. 57) then
        write(LIS_logunit,*)'- FIRST GUESS DATA IS FROM UK UM (GALWEM) MODEL'
     else
        write(LIS_logunit,*)'- UNKNOWN SOURCE FOR FIRST GUESS DATA'
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

     call AGRMET_fldbld_read_precip_galwem(avnfile, ifguess, jfguess,&
                                           fg_prec1, alert_number)
     if (getsixhr.eq.1) &
          call AGRMET_fldbld_read_precip_galwem(avnfile2, ifguess, jfguess,&
                                                fg_prec2, alert_number)

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

     call interp_galwem_first_guess(n, ifguess, jfguess, .true., &
                                    fg_prec, fg_data)

     deallocate ( fg_prec )
     deallocate ( fg_prec1 )
     if (getsixhr.eq.1) deallocate ( fg_prec2 )
  else
        write(LIS_logunit,*)'- ** GALWEM Precipitation data not available **'

        message(1) = 'Program:  LIS'
        message(2) = '  Routine:  fldbld_precip.'
        message(3) = '  galwem Precipitation data not available, ' // &
                     'possible degradation.'
        alert_number = alert_number + 1
        if(LIS_masterproc) then
           call LIS_alert( 'fldbld              ', alert_number, message )
        endif
  endif
end subroutine AGRMET_fldbld_precip_galwem


!BOP
!
! !ROUTINE: AGRMET_fldbld_read_precip_galwem
!  \label{AGRMET_fldbld_read_precip_galwem}
!
! !REVISION HISTORY:
! 11 Aug 2016  Initial specification based on AGRMET_fldbld_read_precip_gfs
!              ...........................................James Geiger/NASA
!
! !INTERFACE:
subroutine AGRMET_fldbld_read_precip_galwem(fg_filename, ifguess, jfguess,&
                                            fg_prec, alert_number )
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
!     to read UK Unified Model (GALWEM) data in GRIB-2 format.
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
!        total precipitation read in from galwem file
!       \item[alert\_number]
!        counts number of alert messages sent
!       \item[cstat]
!        I/O status, character
!       \item[message]
!        Error message
!       \item[count\_prec]
!        counts number of precipitation
!        levels read in from first guess file:q
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
  character*100                 :: message     ( 20 )
  integer                       :: count_prec
  integer                       :: i
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: j
  integer                       :: k, c, r
  integer                       :: ftn, igrib, nvars
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, forecasttime_val
  real,           allocatable   :: dum1d       ( : )
!<debug -- jim testing>
character(len=40) :: jim_name
!</debug -- jim testing>


!     ------------------------------------------------------------------
!     executable code begins here ...
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     read in grib file.
!     ------------------------------------------------------------------

#if (defined USE_GRIBAPI)
  call grib_open_file(ftn,trim(fg_filename),'r',ierr)
  call LIS_verify(ierr,'Failed to open - '//trim(fg_filename))

  if ( ierr == 0 ) then
     allocate ( dum1d   (ifguess*jfguess) )
     count_prec = 0

     write(LIS_logunit,*)' '
     write(LIS_logunit,*)'[MSG] Reading first guess precip ', trim(fg_filename)

     call grib_count_in_file(ftn,nvars,ierr)
     call LIS_verify(ierr, 'error in grib_count_in_file in ' // &
                           'AGRMET_fldbld_read_precip_galwem')

!<debug -- jim testing>
if ( LIS_masterproc ) then
!20160411_CY.00_FH.000_DF.GR2_sp.bin
!--------1---------2---------3---------4---------5---------6---------7---------8----5
!./input/AFWA/20160411/GALWEM/PS.557WW_SC.U_DI.F_GP.GALWEM-GD_GR.C17KM_AR.GLOBAL_DD.20160411_CY.00_FH.003_DF.GR2
jim_name=trim(fg_filename(84:))
jim_name=trim(jim_name)//'_prec.bin'
open(unit=666,file=trim(jim_name),access='direct',recl=ifguess*jfguess*4)
endif
!</debug -- jim testing>
     do k = 1, nvars

        call grib_new_from_file(ftn,igrib,ierr)
        call LIS_verify(ierr, 'failed to read - '//trim(fg_filename))

        call grib_get(igrib,'discipline',param_disc_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                              'AGRMET_fldbld_read_precip_galwem')

        call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterCategory in ' // &
                              'AGRMET_fldbld_read_precip_galwem')

        call grib_get(igrib,'parameterNumber',param_num_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                              'AGRMET_fldbld_read_precip_galwem')
        
        call grib_get(igrib,'forecastTime',forecasttime_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: forecastTime in ' // &
                              'AGRMET_fldbld_read_precip_galwem')

        if ( param_disc_val == 0 .and. param_cat_val    == 1 .and. &
             param_num_val  == 8 .and. forecasttime_val == 0 ) then

           call grib_get(igrib,'values',dum1d,ierr)
           call LIS_verify(ierr, 'error in grib_get: values in ' // &
                                 'AGRMET_fldbld_read_precip_galwem')

           if ( ierr == 0 ) then 
!     ------------------------------------------------------------------
!       store total precipitation to the fg_prec array.
!     ------------------------------------------------------------------
              fg_prec = reshape(dum1d, (/ifguess,jfguess/))
              count_prec = count_prec + 1
!<debug -- jim testing>
if ( LIS_masterproc ) then
write(666,rec=1) fg_prec
endif
!</debug -- jim testing>

              ! found precip; clean up and break out of do-loop.
              call grib_release(igrib,ierr)
              exit
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
        endif
        call grib_release(igrib,ierr)
     enddo
!<debug -- jim testing>
if ( LIS_masterproc ) then
flush(666)
close(666)
endif
!</debug -- jim testing>

     call grib_close_file(ftn)

     deallocate ( dum1d )

!     ------------------------------------------------------------------
!     see if we have everything. if not, abort.
!     ------------------------------------------------------------------

     if ( count_prec == 0 ) then
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  AGRMET_fldbld_read_precip_galwem.'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(fg_filename)
        message(5) = '  File does not contain all data that fldbld needs.'
        if ( allocated(dum1d) ) deallocate(dum1d)
        call LIS_abort( message)
     endif
  else
     write(cstat,'(i9)',iostat=istat1) ierr
     message(1) = 'Program: LIS'
     message(2) = '  Subroutine:  AGRMET_fldbld_read_precip_galwem.'
     message(3) = '  Error opening first guess file:'
     message(4) = '  ' // trim(fg_filename)
     if ( istat1 == 0 ) then
        message(5) = '  Status = ' // trim(cstat)
     endif
     call LIS_abort(message)
  endif
#endif
end subroutine AGRMET_fldbld_read_precip_galwem
