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
! !ROUTINE: AGRMET_fldbld_galwem
! \label{AGRMET_fldbld_galwem}
!
! !REVISION HISTORY:
! 14 Jun 2016  Initial specification based on AGRMET_fldbld_gfs
!              ...........................................James Geiger/NASA
!
! !INTERFACE:    
subroutine AGRMET_fldbld_galwem(n,order,julhr)
! !USES: 
  use LIS_coreMod,       only : LIS_rc
  use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_verify
  use LIS_timeMgrMod,    only : LIS_julhr_date
  use AGRMET_forcingMod, only : agrmet_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: order
!
! !DESCRIPTION: 
!  This routine interpolates the UK Unified Model (GALWEM) first guess height,
!  temperature, moisture, and wind data to the AGRMET grid.
!
! The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[order]
!    specifies which forcing bookend is being processed
!    1 = right/current bookend
!    2 = left/previous bookend
!  \item[avnfile]
!    name of the first guess file
!  \item[message]
!    error message
!  \item[iginfo]
!   array of grid information from a grib record 
!  \item[gridres]
!    the resolution, in degrees, of the first 
!    guess grid
!  \item[alert\_number]
!    number of alerts that occur in the program
!  \item[ifguess]
!    east/west dimension of first guess grid
!  \item[jfguess]
!    north/south dimension of first guess grid
!  \item[kprs]
!    number of isobaric levels for first guess data
!  \item[prslvls]
!    isobaric levels for first guess data
!  \item[center]
!    meteorological center that produced the first 
!    guess data (7-NCEP, 57-AFGWC, 58-FNMOC)
!  \item[ierr]
!    error code
!  \item[fc\_hr]
!    forecast hour or the difference between reference and valid time
!  \item[yr\_2d]
!    2 digit year for comparison with GRIB header
!  \item[found]
!    logical flag set true when an input file with the correct valid
!    time is found
!  \item[file\_julhr]
!    julian hour used to determine names of forecast files from previous cycles
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!    converts the julian hour to a date format
!  \item[AGRMET_getGALWEMfilename](\ref{AGRMET_getGALWEMfilename}) \newline
!    generates the first guess GALWEM filename
!  \item[AGRMET\_fldbld\_read\_galwem](\ref{AGRMET_fldbld_read_galwem}) \newline
!   read GALWEM data in grib format
!  \item[interp\_galwem\_first\_guess](\ref{interp_galwem_first_guess}) \newline
!   interpolate first guess data to the LIS grid
!  \item[lis\_abort](\ref{LIS_abort}) \newline
!   abort in case of error
!  \end{description}
!EOP
  integer                 :: ftn, igrib
  character*120           :: avnfile
  integer                 :: yr1, mo1, da1, hr1
  integer                 :: julhr
  character*100           :: message     ( 20 )
  integer                 :: iginfo      ( 40 )
  real                    :: gridres
  integer                 :: alert_number
  integer                 :: ifguess, jfguess
  integer                 :: kprs 
  integer                 :: prslvls      (30) 
  integer                 :: center
  integer                 :: ierr
  logical*1               :: found
  integer                 :: yr_2d
  integer                 :: file_julhr
  character*100           :: gtype
  integer                 :: fc_hr
  integer                 :: dataDate, dataTime

  data prslvls / 1000,975,950,925,900,850,800,750,700,650,600,550,500,&
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
  kprs = 13
  call LIS_julhr_date(julhr,yr1,mo1,da1,hr1)
  file_julhr = julhr
  
!     ------------------------------------------------------------------
!     Need to process the current and previous 6 hour instances
!     Search for an analysis or forecast file for upto 24 hours with 
!     the needed valid time
!     ------------------------------------------------------------------       

  fc_hr = 0 
  found = .FALSE. 
  do while( (.not.found) .and. (fc_hr <= 24))
     yr_2d = mod(yr1,100)
     if(yr_2d.eq.0) yr_2d = 100 
     call AGRMET_getGALWEMfilename(avnfile, agrmet_struc(n)%agrmetdir,&
          agrmet_struc(n)%galwemdir, agrmet_struc(n)%use_timestamp,&
          yr1,mo1,da1,hr1,fc_hr)

!     ------------------------------------------------------------------
!     open first guess grib data using library utility.  just read
!     the first file only, as all data will be of the same type
!     (avn or nogaps) because the search script ensures that it is.
!     ------------------------------------------------------------------
#if (defined USE_GRIBAPI) 
     
     call grib_open_file(ftn,trim(avnfile),'r',ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) 'Failed to open - '//trim(avnfile)
     else

!     ------------------------------------------------------------------
!     read in the first grib record, unpack the header and extract
!     section 1 and section 2 information.
!     ------------------------------------------------------------------   
        call grib_new_from_file(ftn,igrib,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'failed to read - '//trim(avnfile)
        endif

        call grib_get(igrib,'centre',center,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grib_get: ' // &
                                'centre in AGRMET_fldbld_galwem'
        endif
        
        call grib_get(igrib,'gridType',gtype,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: ' // &
                                'gridtype in AGRMET_fldbld_galwem'
        endif

        call grib_get(igrib,'Ni',iginfo(1),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: Ni in AGRMET_fldbld_galwem'
        endif

        call grib_get(igrib,'Nj',iginfo(2),ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: Nj in AGRMET_fldbld_galwem'
        endif

        call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: ' // &
                                'jDirectionIncrementInDegrees ' // &
                                'in AGRMET_fldbld_galwem'
        endif

        call grib_get(igrib,'dataDate',dataDate,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: ' // &
                                'dataDate in AGRMET_fldbld_galwem'
        endif

        call grib_get(igrib,'dataTime',dataTime,ierr)
        if ( ierr /= 0 ) then
           write(LIS_logunit,*) 'error in grid_get: ' // &
                                'dataTime in AGRMET_fldbld_galwem'
        endif

        if ( yr1*10000+mo1*100+da1 == dataDate .and. hr1*100 == dataTime ) then
           found = .true.
        endif
    
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
     endif
#endif     

!    ------------------------------------------------------------------   
!    If the correct valid time is not found:
!       Increment forecast hour by 6.
!       Decrement file_julhr by 6 and get the new filename elements.
!    ------------------------------------------------------------------   

     if (.not. found) then
        fc_hr = fc_hr + 6
        file_julhr = file_julhr - 6
        call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)
     endif
  enddo
  
  if (.not. found) then
     message(1) = 'program: LIS'
     message(2) = '  Subroutine: fldbld_galwem'
     message(3) = '  No matching avn file found:'
     call LIS_abort(message)
  else
     
     if(trim(gtype).ne."regular_ll") then  
        message(1) = 'program: LIS'
        message(2) = '  Subroutine: agrmet_sfcalc'
        message(3) = '  First guess source is not a lat/lon grid'
        message(4) = '  agrmet_sfcalc expects lat/lon data'
        call lis_abort(message)
     endif
!     ------------------------------------------------------------------
!     determine resolution of grid.  this is in section 2 octet 24-25
!     and is stored in millidegrees (iginfo(8)).  also get the 
!     i and j dimensions of the grid, section 2, octets 7-8 and 9-10
!     (iginfo(1 and 2).  see grib utility and grib manual for
!     more details.
!     ------------------------------------------------------------------
!     gridres = float(iginfo(8))/1000.0
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
!         read in first guess data for this julian hour.
!     ------------------------------------------------------------------
     alert_number = 0 
     
     if(order.eq.1) then 
        call AGRMET_fldbld_read_galwem(n, avnfile, ifguess, jfguess,  &
                kprs, prslvls, alert_number,                          &
                agrmet_struc(n)%wndwgt, agrmet_struc(n)%minwnd,       &
                agrmet_struc(n)%agr_tmp_c, agrmet_struc(n)%agr_hgt_c, &
                agrmet_struc(n)%agr_rh_c, agrmet_struc(n)%agr_wspd_c, &
                agrmet_struc(n)%agr_pres_c)
     else
        call AGRMET_fldbld_read_galwem(n, avnfile, ifguess, jfguess,  &
                kprs, prslvls, alert_number,                          &
                agrmet_struc(n)%wndwgt, agrmet_struc(n)%minwnd,       &
                agrmet_struc(n)%agr_tmp_p, agrmet_struc(n)%agr_hgt_p, &
                agrmet_struc(n)%agr_rh_p, agrmet_struc(n)%agr_wspd_p, &
                agrmet_struc(n)%agr_pres_p)
     endif

  endif
end subroutine AGRMET_fldbld_galwem


!BOP
! 
! !ROUTINE: getGALWEMfilename
! \label{getGALWEMfilename}
!
! !INTERFACE: 
subroutine AGRMET_getGALWEMfilename(filename,rootdir,dir,use_timestamp, &
                             yr,mo,da,hr,fc_hr)

  implicit none
! !ARGUMENTS: 
  character(*)        :: filename
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: yr,mo,da,hr
  integer, intent(in) :: fc_hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the timestamped GALWEM file
! 
!  The arguments are: 
!  \begin{description}
!   \item[filename]
!    created filename
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[use\_timestamp]
!    flag to indicate whether the directories 
!    should be timestamped or not
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!   \item[fc\_hr]
!    hours from reference time to valid time
!  \end{description}
!EOP
  character(8) :: ftime1
  character(2) :: fhr
  character(3) :: fchr
  
  character(len=54) :: fname1
  character(len=20) :: fname2

  write (UNIT=fhr, FMT='(i2.2)') hr
  write (UNIT=fchr, FMT='(i3.3)') fc_hr

  fname1 = 'PS.557WW_SC.U_DI.F_GP.GALWEM-GD_GR.C17KM_AR.GLOBAL_DD.'

  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(i4, i2.2, i2.2)') yr, mo, da
     filename = trim(rootdir) // '/' // ftime1 // '/' // trim(dir) // '/' // &
                fname1 // ftime1 // '_CY.' // fhr // '_FH.' // fchr // '_DF.GR2'
  else
     filename = trim(rootdir) // '/' // trim(dir) // '/' // &
                fname1 // ftime1 // '_CY.' // fhr // '_FH.' // fchr // '_DF.GR2'
  endif
end subroutine AGRMET_getGALWEMfilename


!BOP
!
! !ROUTINE: AGRMET_fldbld_read_galwem
!  \label{AGRMET_fldbld_read_galwem}
!
! !REVISION HISTORY:
! 14 Jun 2016  Initial specification based on AGRMET_fldbld_read_gfs
!              ...........................................James Geiger/NASA
!    
! 
! !INTERFACE: 
subroutine AGRMET_fldbld_read_galwem(n, fg_filename, ifguess, jfguess,   &
                                     kprs, prslvls, alert_number,        &
                                     wndwgt, minwnd,                     &
                                     agr_tmp, agr_hgt, agr_rh, agr_wspd, &
                                     agr_pres)
! !USES:
!<debug -- jim testing>
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_mpiMod
  use LIS_historyMod
!</debug -- jim testing>
  !use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer,        intent(in)    :: n
  character(len=*),  intent(in) :: fg_filename
  integer,        intent(in)    :: ifguess
  integer,        intent(in)    :: jfguess
  integer,        intent(in)    :: kprs
  integer,        intent(in)    :: prslvls(30)
  integer,        intent(inout) :: alert_number
  real,           intent(in)    :: minwnd
  real,           intent(in)    :: wndwgt
  real,           intent(out)   :: agr_hgt (LIS_rc%lnc(n),LIS_rc%lnr(n),kprs)
  real,           intent(out)   :: agr_rh  (LIS_rc%lnc(n),LIS_rc%lnr(n),kprs)
  real,           intent(out)   :: agr_tmp (LIS_rc%lnc(n),LIS_rc%lnr(n),kprs)
  real,           intent(out)   :: agr_wspd(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,           intent(out)   :: agr_pres(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:  
!
!     to read UK Unified Model (GALWEM) data in GRIB-2 format.
!     
!     \textbf{Method} \newline
!     
!     - open file and allocate variables. \newline
!     - read in each grib record. \newline
!     - if data is what we need, store it in the proper arrays. \newline
!       also, keep track of what has been read in using \newline
!       counter variables. \newline
!     - check forecast hour of data, if it is not the analysis \newline
!       or zero hour forecast, then send an alert message \newline
!       to warn of possible degradation. \newline
!     - check counter variables.  if data is missing abort. \newline
!     - calculate surface wind speed from u and v components. \newline
!       multiply by wind weighting factor. restrict lowest \newline
!       wind speed to value of minwnd.
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
!      \item[kprs]
!        number of isobaric levels where fldbld
!        needs first guess data
!      \item[prslvls]
!        the isobaric levels where fldbld needs
!        first guess data
!      \item[alert\_number]
!        counts number of alert messages sent
!      \item[wndwgt]
!        adjustment factor for first guess
!        winds.
!      \item[minwnd]
!        minimum allowable wind speed on the
!        agrmet grid
!      \item[agr\_hgt]
!        isobaric heights
!      \item[agr\_rh]
!        isobaric relative humidity
!      \item[agr\_tmp]
!        isobaric temperatures
!      \item[agr\_wspd]
!        surface wind speeds
!      \item[agr\_pres]
!        surface pressure
!      \item[cstat]
!        I/O status, character
!      \item[message]
!        Error message
!      \item[count\_hgt]
!        counts number of isobaric height levels
!        read in from first guess file
!      \item[count\_rh]
!        counts number of isobaric relative humidity levels
!        read in from first guess file
!      \item[count\_tmp]
!        counts number of isobaric temperature levels
!        read in from first guess file
!      \item[count\_uwnd]
!        counts number of isobaric u wind levels
!        read in from first guess file
!      \item[count\_vwnd]
!        counts number of isobaric v wind levels
!        read in from first guess file
!      \item[i,j,kk]
!        looping and indexing variables
!      \item[ierr,istat1]
!        error status
!      \item[ftn]
!        unix file descriptor
!      \item[dum1d]
!        dummy array
!      \item[fg\_hgt]
!        array of first guess isobaric heights
!      \item[fg\_rh]
!        array of first guess isobaric relative humidity
!      \item[fg\_tmp]
!        array of first guess isobaric temperatures
!      \item[fg\_wspd]
!        array of first guess surface wind speeds
!      \item[fg\_uwnd]
!        array of first guess surface u winds
!      \item[fg\_vwnd]
!        array of first guess surface v winds
!      \item[nvars]
!        total number of messages in the GRIB2 file
!        array of first guess surface v winds
!     \end{description}
!EOP
  character*9                   :: cstat
  character*100                 :: message     ( 20 )
  character(len=4)              :: grib_msg
  character(len=4)              :: AGRMET_check_galwem_message
  integer                       :: count_hgt
  integer                       :: count_rh
  integer                       :: count_tmp
  integer                       :: count_uwnd
  integer                       :: count_vwnd
  integer                       :: i
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: j,c,r
  integer                       :: igrib
  integer                       :: ftn
  integer                       :: kk, nvars
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, surface_val, level_val
  
  real, allocatable :: dum1d   ( : )
  real, allocatable :: fg_hgt  ( : , : , : )
  real, allocatable :: fg_rh   ( : , : , : )
  real, allocatable :: fg_tmp  ( : , : , : )
  real, allocatable :: fg_wspd ( : , : )
  real, allocatable :: fg_pres ( : , : )
  real, allocatable :: fg_uwnd ( : , : )
  real, allocatable :: fg_vwnd ( : , : )
!<debug -- jim testing>
  character(len=100) :: jim_name
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

  if ( ierr .eq. 0 ) then 
     allocate ( fg_hgt  (ifguess, jfguess, kprs) )
     allocate ( fg_rh   (ifguess, jfguess, kprs) )
     allocate ( fg_tmp  (ifguess, jfguess, kprs) )
     allocate ( fg_wspd (ifguess, jfguess) )
     allocate ( fg_pres (ifguess, jfguess) )

     allocate ( dum1d   (ifguess*jfguess) )
     allocate ( fg_uwnd (ifguess, jfguess) )
     allocate ( fg_vwnd (ifguess, jfguess) )
     
     count_hgt  = 0
     count_rh   = 0
     count_tmp  = 0
     count_uwnd = 0
     count_vwnd = 0
     
     write(LIS_logunit,*)' '
     write(LIS_logunit,*)'- READING ', trim(fg_filename)

     call grib_count_in_file(ftn,nvars,ierr)
     call LIS_verify(ierr, 'error in grib_count_in_file in ' // &
                           'AGRMET_fldbld_read_galwem')

!<debug -- jim testing>
     if ( LIS_masterproc ) then
        kk=index(fg_filename, 'GLOBAL_DD')
        jim_name=fg_filename(kk+10:)
        open(unit=666,file=trim(jim_name)//'_sp.bin',access='direct',recl=ifguess*jfguess*4)
        open(unit=667,file=trim(jim_name)//'_t.bin',access='direct',recl=ifguess*jfguess*4)
        open(unit=668,file=trim(jim_name)//'_gh.bin',access='direct',recl=ifguess*jfguess*4)
        open(unit=669,file=trim(jim_name)//'_r.bin',access='direct',recl=ifguess*jfguess*4)
        open(unit=670,file=trim(jim_name)//'_10u.bin',access='direct',recl=ifguess*jfguess*4)
        open(unit=671,file=trim(jim_name)//'_10v.bin',access='direct',recl=ifguess*jfguess*4)
     endif
!</debug -- jim testing>

     do kk=1,nvars
 
        call grib_new_from_file(ftn,igrib,ierr)
        call LIS_verify(ierr, 'failed to read - '//trim(fg_filename))
        
        call grib_get(igrib,'discipline',param_disc_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                              'AGRMET_fldbld_read_galwem')

        call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterCategory in ' // &
                              'AGRMET_fldbld_read_galwem')

        call grib_get(igrib,'parameterNumber',param_num_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: parameterNumber in ' // &
                              'AGRMET_fldbld_read_galwem')
        
        call grib_get(igrib,'typeOfFirstFixedSurface',surface_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: level in ' // &
                              'AGRMET_fldbld_read_galwem')

        call grib_get(igrib,'level',level_val,ierr)
        call LIS_verify(ierr, 'error in grib_get: level in ' // &
                              'AGRMET_fldbld_read_galwem')

        grib_msg = AGRMET_check_galwem_message(param_disc_val, param_cat_val, &
                                        param_num_val, surface_val)

        if ( grib_msg /= 'none' ) then

           call grib_get(igrib,'values',dum1d,ierr)
        
           if ( ierr == 0 ) then 

              select case (grib_msg)
              case('sp') ! surface pressure
                 fg_pres = reshape(dum1d, (/ifguess,jfguess/))
!<debug -- jim testing>
                 if ( LIS_masterproc ) then
                    write(666,rec=1) fg_pres
                 endif
!</debug -- jim testing>
              case('t') ! temperature
                 do i = 1, kprs
                    if ( level_val == prslvls(i) ) then
                       fg_tmp(:,:,i) = reshape(dum1d, (/ifguess,jfguess/))
                       count_tmp = count_tmp + 1
!<debug -- jim testing>
                       if ( LIS_masterproc ) then
                          write(667,rec=i) fg_tmp(:,:,i)
                       endif
!</debug -- jim testing>
                       exit
                    endif
                 enddo
              case('gh') ! geopotential height
                 do i = 1, kprs
                    if ( level_val == prslvls(i) ) then
                       fg_hgt(:,:,i) = reshape(dum1d, (/ifguess,jfguess/))
                       count_hgt = count_hgt + 1
!<debug -- jim testing>
                       if ( LIS_masterproc ) then
                          write(668,rec=i) fg_hgt(:,:,i)
                       endif
!</debug -- jim testing>
                       exit
                    end if
                 enddo
              case('r') ! relative humidity
                 ! according to WMO standard,
                 ! relative humidity is gribbed in percent.
                 ! agrmet expects decimal, so divide by 100.
                 do i = 1, kprs
                    if ( level_val == prslvls(i) ) then
                       dum1d = dum1d / 100.0
                       fg_rh(:,:,i) = reshape(dum1d, (/ifguess,jfguess/))
                       count_rh = count_rh + 1
!<debug -- jim testing>
                       if ( LIS_masterproc ) then
                          write(669,rec=i) fg_rh(:,:,i)
                       endif
!</debug -- jim testing>
                       exit
                    end if
                 enddo
              case('10u') ! 10m u-wind
                 fg_uwnd = reshape(dum1d, (/ifguess,jfguess/))
                 count_uwnd = 1
!<debug -- jim testing>
                 if ( LIS_masterproc ) then
                    write(670,rec=1) fg_uwnd
                 endif
!</debug -- jim testing>
              case('10v') ! 10m v-wind
                 fg_vwnd = reshape(dum1d, (/ifguess,jfguess/))
                 count_vwnd = 1
!<debug -- jim testing>
                 if ( LIS_masterproc ) then
                    write(671,rec=1) fg_vwnd
                 endif
!</debug -- jim testing>
              case default
                 write(cstat,'(i9)',iostat=istat1) ierr
                 message(1) = 'Program: LIS'
                 message(2) = '  Subroutine:  AGRMET_fldbld_read_galwem.'
                 message(3) = '  Error reading first guess file:'
                 message(4) = '  ' // trim(fg_filename)
                 if( istat1 .eq. 0 )then
                    message(5) = '  Status = ' // trim(cstat)
                 endif
                 call LIS_abort( message)
              end select

           endif
        endif
        call grib_release(igrib,ierr)
     enddo
!<debug -- jim testing>
     if ( LIS_masterproc ) then
        close(666)
        close(667)
        close(668)
        close(669)
        close(670)
        close(671)
     endif
!</debug -- jim testing>
     
     call grib_close_file(ftn)

!     ------------------------------------------------------------------
!     see if we have everything. if not, abort.
!     ------------------------------------------------------------------
     
     if ( .not. ( (count_tmp .lt. kprs) .or. (count_uwnd .eq. 0) .or. &
                  (count_hgt .lt. kprs) .or. (count_vwnd .eq. 0) .or. &
                  (count_rh  .lt. kprs) ) ) then 
        
!     ------------------------------------------------------------------
!     convert from u and v component winds to wind speed. 
!     contrain the speeds between 50 m/s and the user
!     specified minimum.
!     ------------------------------------------------------------------
  
        fg_wspd = sqrt((fg_uwnd * fg_uwnd) + (fg_vwnd * fg_vwnd)) * wndwgt
        
        fg_wspd = min( max( fg_wspd, minwnd ), 50.0 ) 

        call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
                                       fg_pres, agr_pres)
        do i = 1, kprs
           call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
                                          fg_tmp(:,:,i), agr_tmp(:,:,i))
           call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
                                          fg_hgt(:,:,i), agr_hgt(:,:,i))
           call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
                                          fg_rh(:,:,i), agr_rh(:,:,i))
        enddo
        call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
                                       fg_wspd, agr_wspd)
!<debug -- jim testing>
        if ( LIS_masterproc ) then
           kk=index(fg_filename, 'GLOBAL_DD')
           jim_name=fg_filename(kk+10:)
           open(unit=666,file=trim(jim_name)//'_lis.bin',access='direct',recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
        endif
        kk = 1
        call LIS_writevar_bin(666,n,agr_pres,kk)
        do i = 1, kprs
           kk = kk + 1
           call LIS_writevar_bin(666,n,agr_tmp(:,:,i),kk)
        enddo
        do i = 1, kprs
           kk = kk + 1
           call LIS_writevar_bin(666,n,agr_hgt(:,:,i),kk)
        enddo
        do i = 1, kprs
           kk = kk + 1
           call LIS_writevar_bin(666,n,agr_rh(:,:,i),kk)
        enddo
        kk = kk + 1
        call LIS_writevar_bin(666,n,agr_wspd,kk)
        if ( LIS_masterproc ) then
           close(666)
        endif
!</debug -- jim testing>

     else
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  AGRMET_fldbld_read_galwem.'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(fg_filename)
        message(5) = '  File does not contain all data that fldbld needs.'
        call LIS_abort( message)
     endif

     deallocate ( dum1d )
     deallocate ( fg_uwnd )
     deallocate ( fg_vwnd )

     deallocate ( fg_hgt  )
     deallocate ( fg_rh   )
     deallocate ( fg_tmp  )
     deallocate ( fg_wspd )
     deallocate ( fg_pres )
  else
     write(cstat,'(i9)',iostat=istat1) ierr
     message(1) = 'Program: LIS'
     message(2) = '  Subroutine:  AGRMET_fldbld_read_galwem.'
     message(3) = '  Error opening first guess file:'
     message(4) = '  ' // trim(fg_filename)
     if( istat1 .eq. 0 )then
        message(5) = '  Status = ' // trim(cstat)
     endif
     call LIS_abort( message)
  endif
#endif

end subroutine AGRMET_fldbld_read_galwem

!BOP
!
! !ROUTINE: AGRMET_check_galwem_message
! \label{AGRMET_check_galwem_message}
!
! !REVISION HISTORY:
! 14 Jun 2016 James Geiger; Initial specification
!
! !INTERFACE:    
function AGRMET_check_galwem_message(param_disc_val, param_cat_val, &
                              param_num_val, surface_val)
! !USES: 
! none

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: param_disc_val, param_cat_val, &
                          param_num_val, surface_val
   character(len=4)    :: AGRMET_check_galwem_message
!
! !DESCRIPTION: 
!  This function compares given grib id values against desired values
!  and returns a string indicting with variable matches those given
!  grib id values.  Returns 'none' if there is no match.
!
! The arguments are:
! \begin{description}
!    \item[param\_disc\_val]
!       value of discipline key from GRIB2 message
!    \item[param\_cat\_val]
!      value of parameterCategory key from GRIB2 message
!    \item[param\_num\_val]
!       value of parameterNumber key from GRIB2 message
!    \item[surface\_val]
!       value of typeOfFirstFixedSurface key from GRIB2 message
! \end{description}
!EOP
   if     ( param_disc_val == 0 .and. &
            param_cat_val  == 3 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 1 ) then
      AGRMET_check_galwem_message = 'sp'
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 0 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 100 ) then
      AGRMET_check_galwem_message = 't'
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 3 .and. &
            param_num_val  == 5 .and. &
            surface_val    == 100 ) then
      AGRMET_check_galwem_message = 'gh'
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 1 .and. &
            param_num_val  == 1 .and. &
            surface_val    == 100 ) then
      AGRMET_check_galwem_message = 'r'
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 2 .and. &
            surface_val    == 103 ) then
      AGRMET_check_galwem_message = '10u'
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 3 .and. &
            surface_val    == 103 ) then
      AGRMET_check_galwem_message = '10v'
   else 
      AGRMET_check_galwem_message = 'none'
   endif
 end function AGRMET_check_galwem_message


!BOP
!
! !ROUTINE: interp_galwem_first_guess
! \label{interp_galwem_first_guess}
!
! !REVISION HISTORY: 
!  01 Jul 2016 James Geiger; Initial Specification
!
! !INTERFACE: 
subroutine interp_galwem_first_guess(n, ifguess, jfguess, prec_flag, &
                                     fg_field, agr_field)
  
! !USES: 
   use LIS_coreMod,       only : LIS_rc, LIS_domain
   use LIS_logMod,        only : LIS_logunit, LIS_endrun
   use AGRMET_forcingMod, only : agrmet_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in)  :: n
   integer, intent(in)  :: ifguess
   integer, intent(in)  :: jfguess
   logical, intent(in)  :: prec_flag
   real,    intent(in)  :: fg_field  ( ifguess,jfguess )
   real,    intent(out) :: agr_field ( LIS_rc%lnc(n),LIS_rc%lnr(n) )
!     
! !DESCRIPTION:    
! 
! This routine interpolates the GALWEM first guess data to the LIS/AGRMET grid.
! 
! The arguments are: 
! \begin{description}
!  \item[n]
!    nest index
!  \item[ifguess]
!    i-dimension of the first guess grid
!  \item[jfguess]
!    j-dimension of the first guess grid
!  \item[fg\_field]
!   first guess field input array
!  \item[agr\_field]
!    agrmet grid output array
!  \end{description}
!EOP

   integer   :: mi, mo
   integer   :: k 
   integer   :: i,j
   integer   :: iret
   integer   :: midway
   character(len=50) :: method
   real, allocatable, dimension(:,:)    :: var
   logical*1, allocatable, dimension(:) :: lb
   logical*1, allocatable, dimension(:) :: lo

   allocate(var(ifguess,jfguess))
   allocate(lb(ifguess*jfguess))
   allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

   mi = ifguess * jfguess
   mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
   lb = .true. 
   
   ! translate from (0,360) to (-180,180).
   midway = ifguess/2
   do j = 1, jfguess
      do i = 1, ifguess
         if ( (i+midway) < ifguess ) then
            var(i,j) = fg_field((i+midway),j)
         else
            var(i,j) = fg_field((i-midway+1),j)
         endif
      enddo
   enddo


  method = agrmet_struc(n)%fg_galwem_interp

  ! When budget-bilinear is selected, it is applied only to precip,
  ! all other fields use bilinear.
  if ( method == 'budget-bilinear' ) then
     if ( .not. prec_flag ) then
        method = 'bilinear'
     endif
  endif

  if ( method == 'bilinear' ) then
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                      &
             var,lo,agr_field,mi,mo,                                    &
             LIS_domain(n)%lat,LIS_domain(n)%lon,                       &
             agrmet_struc(n)%w11_1_galwem,agrmet_struc(n)%w12_1_galwem, &
             agrmet_struc(n)%w21_1_galwem,agrmet_struc(n)%w22_1_galwem, &
             agrmet_struc(n)%n11_1_galwem,agrmet_struc(n)%n12_1_galwem, &
             agrmet_struc(n)%n21_1_galwem,agrmet_struc(n)%n22_1_galwem, &
             LIS_rc%udef, iret)
  elseif ( method == 'budget-bilinear' ) then
     call conserv_interp(LIS_rc%gridDesc(n,:),lb,                       &
             var,lo,agr_field,mi,mo,                                    &
             LIS_domain(n)%lat,LIS_domain(n)%lon,                       &
             agrmet_struc(n)%w11_2_galwem,agrmet_struc(n)%w12_2_galwem, &
             agrmet_struc(n)%w21_2_galwem,agrmet_struc(n)%w22_2_galwem, &
             agrmet_struc(n)%n11_2_galwem,agrmet_struc(n)%n12_2_galwem, &
             agrmet_struc(n)%n21_2_galwem,agrmet_struc(n)%n22_2_galwem, &
             LIS_rc%udef, iret)
  elseif ( method == 'neighbor' ) then
      call neighbor_interp(LIS_rc%gridDesc(n,:),lb,                     &
                           var,lo,agr_field,mi,mo,                      &
                           LIS_domain(n)%lat, LIS_domain(n)%lon,        &
                           agrmet_struc(n)%n11_1_galwem,LIS_rc%udef,iret)
  elseif ( method == 'average' ) then
     call upscaleByAveraging(mi,mo,LIS_rc%udef, &
                             agrmet_struc(n)%n11_1_galwem,lb,var,lo,agr_field)
  else
     write(LIS_logunit,*) 'ERR: Unexpected interpolation method'
     write(LIS_logunit,*) '     in interp_galwem_first_guess'
     write(LIS_logunit,*) '     ', trim(method)
     call LIS_endrun
  endif

   deallocate(var)
   deallocate(lb)
   deallocate(lo)

end subroutine interp_galwem_first_guess

!BOP
!
! !ROUTINE: galwem_reset_interp_input
!  \label{galwem_reset_interp_input}
!
! !REVISION HISTORY:
!  12 Jun 2017: James Geiger: Initial specification
!
! !INTERFACE:
subroutine galwem_reset_interp_input(n, findex, gridDesci)
! !USES:

   use LIS_coreMod,        only : LIS_rc, LIS_howtoTransform
   use LIS_logMod,         only : LIS_logunit, LIS_endrun
   use agrmet_forcingMod,  only : agrmet_struc

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: findex
   real, intent(in)    :: gridDesci(50)

! !DESCRIPTION:
! Resets the neighbours and weights arrays used for spatially
! interpolating the first guess forcing data to the LIS running domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing index
!  \item[gridDesci]
!    array of magic numbers describing the first guess forcing domain
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!   \item[LIS\_howtoTransform](\ref{LIS_howtoTransform}) \newline
!    determines whether LIS' running domain is at a finer resolution
!    than the given resolution
!  \end{description}
!EOP

   integer :: rc
   integer :: mi
   character(len=16) :: howtoTransform

   deallocate(agrmet_struc(n)%n11_1_galwem, stat=rc)
   deallocate(agrmet_struc(n)%n12_1_galwem, stat=rc)
   deallocate(agrmet_struc(n)%n21_1_galwem, stat=rc)
   deallocate(agrmet_struc(n)%n22_1_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w11_1_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w12_1_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w21_1_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w22_1_galwem, stat=rc)

   deallocate(agrmet_struc(n)%n11_2_galwem, stat=rc)
   deallocate(agrmet_struc(n)%n12_2_galwem, stat=rc)
   deallocate(agrmet_struc(n)%n21_2_galwem, stat=rc)
   deallocate(agrmet_struc(n)%n22_2_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w11_2_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w12_2_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w21_2_galwem, stat=rc)
   deallocate(agrmet_struc(n)%w22_2_galwem, stat=rc)

   howtoTransform = LIS_howtoTransform(n,max(gridDesci(9),gridDesci(10)))
!<debug -- jim testing>
   howtoTransform = 'interpolate'
!</debug -- jim testing>

   if ( howtoTransform == 'interpolate') then

      agrmet_struc(n)%fg_galwem_interp = LIS_rc%met_interp(findex)

      write(LIS_logunit,*) 'MSG: The GALWEM forcing resolution is coarser ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Interpolating with the ' // &
                           trim(agrmet_struc(n)%fg_galwem_interp) // ' method.'

      if ( agrmet_struc(n)%fg_galwem_interp == 'bilinear' ) then

         allocate(agrmet_struc(n)%n11_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n12_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n21_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n22_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w11_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w12_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w21_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w22_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,        &
                 agrmet_struc(n)%n11_1_galwem,agrmet_struc(n)%n12_1_galwem, &
                 agrmet_struc(n)%n21_1_galwem,agrmet_struc(n)%n22_1_galwem, &
                 agrmet_struc(n)%w11_1_galwem,agrmet_struc(n)%w12_1_galwem, &
                 agrmet_struc(n)%w21_1_galwem,agrmet_struc(n)%w22_1_galwem)

      elseif ( agrmet_struc(n)%fg_galwem_interp == 'budget-bilinear' ) then

         allocate(agrmet_struc(n)%n11_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n12_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n21_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n22_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w11_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w12_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w21_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w22_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,        &
                 agrmet_struc(n)%n11_1_galwem,agrmet_struc(n)%n12_1_galwem, &
                 agrmet_struc(n)%n21_1_galwem,agrmet_struc(n)%n22_1_galwem, &
                 agrmet_struc(n)%w11_1_galwem,agrmet_struc(n)%w12_1_galwem, &
                 agrmet_struc(n)%w21_1_galwem,agrmet_struc(n)%w22_1_galwem)

         allocate(agrmet_struc(n)%n11_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%n12_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%n21_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%n22_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w11_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w12_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w21_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w22_2_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

         call conserv_interp_input(n,gridDesci,         &
                 agrmet_struc(n)%n11_2_galwem,agrmet_struc(n)%n12_2_galwem, &
                 agrmet_struc(n)%n21_2_galwem,agrmet_struc(n)%n22_2_galwem, &
                 agrmet_struc(n)%w11_2_galwem,agrmet_struc(n)%w12_2_galwem, &
                 agrmet_struc(n)%w21_2_galwem,agrmet_struc(n)%w22_2_galwem)
      endif
   elseif ( howtoTransform == 'neighbor') then
      agrmet_struc(n)%fg_galwem_interp = 'neighbor'

      write(LIS_logunit,*) 'MSG: The GALWEM forcing resolution is comparable ' // &
                           'to the running domain.'
      write(LIS_logunit,*) '     Interpolating with the ' // &
                           trim(agrmet_struc(n)%fg_galwem_interp) // ' method.'

      allocate(agrmet_struc(n)%n11_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

      call neighbor_interp_input(n,gridDesci,agrmet_struc(n)%n11_1_galwem)
   elseif ( howtoTransform == 'upscale' ) then
      agrmet_struc(n)%fg_galwem_interp = LIS_rc%met_upscale(findex)

      write(LIS_logunit,*) 'MSG: The GALWEM forcing resolution is finer ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Upscaling with the ' // &
                           trim(agrmet_struc(n)%fg_galwem_interp) // ' method.'


      select case( agrmet_struc(n)%fg_galwem_interp )
      case( 'average' )
         mi = gridDesci(2) * gridDesci(3)
         allocate(agrmet_struc(n)%n11_1_galwem(mi))

         call upscaleByAveraging_input(gridDesci,                   &
                                       LIS_rc%gridDesc(n,:),        &
                                       mi,                          &
                                       LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                                       agrmet_struc(n)%n11_1_galwem)
      case default
         write(LIS_logunit,*) 'The specified spatial interpolation option '
         write(LIS_logunit,*) 'is not supported for GALWEM.'
         write(LIS_logunit,*) 'LIS is stopping.'
         call LIS_endrun()
      end select
   else
      write(LIS_logunit,*) 'Unexpected spatial transformation method '
      write(LIS_logunit,*) howtoTransform
      write(LIS_logunit,*) 'LIS is stopping.'
      call LIS_endrun()
   endif
end subroutine galwem_reset_interp_input

