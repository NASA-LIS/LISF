!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: AGRMET_getsfc
! \label{AGRMET_getsfc}
!
! !REVISION HISTORY: 
! 
!     09 Jul 94  Initial version...............Capt Bertone,SYSM/AGROMET  
!     06 Dec 95  Changed the second call of offhr2 to call original   
!                routine offhr1........................SSgt Erkkila/SYSM  
!     12 Sep 97  Eliminated calls to onhour, offhr1 and offhr2.  
!                Converted from 6-hour to 1-hour retrieval of all
!                observations. Brought code up to AFGWC & SYSM software   
!                standards........Mr Moore(AGROMET), SSgt McCormick/SYSM
!     28 Feb 99  Modified to process data by hemisphere instead of by 
!                box, and to pull data from CDMS, to port the program
!                from System 5 to GTWAPS (workstation).....Mr Moore/DNXM
!     21 Ful 99  Ported to IBM SP-2.  Fixed some bugs.  Revamped
!                observation quality control logic. Removed hardwired
!                AGRMET grid dimensions....................Mr Gayno/DNXM
!     30 Dec 99  Fixed error in call to agr_alert - the first argument
!                was specified incorrectly.................Mr Gayno/DNXM
!     19 Jul 04  Added tracking of CDMS access time.....Mr Lewiston/DNXM
!     11 Mar 10  Changed program names in messages to LIS.  Left argument
!                to LIS_alert with old program name to keep alert numbers
!                unique........................Chris Franks/16WS/WXE/SEMS
!     30 Jul 10  Limit the number of obs read to the size of the array
!                ..............................Chris Franks/16WS/WXE/SEMS
!      9 Sep 10  Modified to read data from JMOBS and to calculate RH from
!                temp & dew pt if RH not provided.Chris Franks/16WS/WXE/SEMS
!     25 Jun 20  Modified to check valid times of surface obs.
!                ..............................Eric Kemp/NASA/SSAI
!     16 Dec 21  Replaced julhr with YYYYMMDDHH in log...Eric Kemp/NASA/SSAI
!     07 Nov 23  Added new sfcobs format.......Eric Kemp/NASA/SSAI
! 
! !INTERFACE: 
subroutine AGRMET_getsfc( n, julhr, t2mObs, rh2mObs, spd10mObs, &
     ri, rj, obstmp, obsrlh, obsspd, obscnt, &
     isize, minwnd, alert_number, imax, jmax, rootdir,cdmsdir,&
     use_timestamp, use_wigos_sfcobs)
! !USES: 
  use ESMF ! EMK Patch for DTG check
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_coreMod,    only  : LIS_domain, LIS_masterproc
  use LIS_timeMgrMod, only  : LIS_julhr_date
  use LIS_logMod,     only  : LIS_logunit, LIS_alert, &
       LIS_getNextUnitNumber,  LIS_releaseUnitNumber
  use map_utils,      only  : latlon_to_ij
  use USAF_bratsethMod, only: USAF_ObsData, USAF_assignObsData


  implicit none
! !ARGUMENTS:   
  integer,    intent(in)         :: n 
  integer,    intent(in)         :: julhr
  integer,    intent(in)         :: imax
  integer,    intent(in)         :: jmax
  integer,    intent(in)         :: isize  
  type(USAF_ObsData), intent(inout) :: t2mObs
  type(USAF_ObsData), intent(inout) :: rh2mObs
  type(USAF_ObsData), intent(inout) :: spd10mObs
  real,       intent(out)        :: ri       ( isize ) 
  real,       intent(out)        :: rj       ( isize ) 
  real,       intent(out)        :: obsrlh   ( isize )   
  real,       intent(out)        :: obsspd   ( isize ) 
  real,       intent(out)        :: obstmp   ( isize )     
  integer,    intent(out)        :: obscnt
  integer,    intent(inout)      :: alert_number
  real,       intent(in)         :: minwnd
  character(len=*)               :: rootdir
  character(len=*)               :: cdmsdir
  integer,    intent(in)         :: use_timestamp
  logical, intent(in) :: use_wigos_sfcobs
!    
! !DESCRIPTION: 
!    
!    Retrieve observed temperature, relative humidity, and
!    wind speed values from the JMOBS database.
!    
!    \textbf{Method} \newline
!     - Set observation counter to zero. \newline
!     - Retrieve all observations for this hemisphere and time from 
!       database. \newline
!     - If there is a read error, send an alert message.  The
!       observation count will remain at zero and no Barnes Analysis
!       will be performed. \newline
!     - Check lat/lon of each observation for missing data and bad
!       values.  Observations that fail check are not stored
!       in the ob arrays passed to routine barnes. \newline
!     - Convert from lat/lon to i and j and hemisphere on the 
!       AGRMET model grid. \newline
!       - Check for observations located outside of grid. \newline
!       - Check for observations converted to other hemisphere. \newline
!       - Observations that fail either check are not stored to
!         the ob arrays. \newline
!     - Check for observations with missing/unreported temperatures,
!       wind speeds and relative humidities. Don't store to 
!       ob arrays if it fails check. \newline
!     - Descale temperature, wind speed and relative humidity for
!       each observation.  Range check data.  If temperature or
!       relative humidity is outside the range, store to ob arrays
!       arrays as minus one so the Barnes Analysis will ignore it.
!       If wind speed is out of range, constrain it to between
!       75 ms-1 and the value of minwnd.  If any of the fields
!       are missing/unreported, store as minus 1. \newline
!     - If all three fields are missing, unreported and/or out of 
!       range, don't store this observation to the ob array. \newline
!     - If the observation passes all the above checks, store
!       the information to the final ob arrays for use in 
!       routine barnes and increment the observation counter. \newline
!    
! The arguments and variables are:
!  \begin{description}
!  \item[julhr]
!    current julian hour
!  \item[ri]
!    Array of observation locations on the AGRMET grid
!   (i dimension)
!  \item[rj]
!    Array of observation locations on the AGRMET grid
!   (j dimension)
!  \item[obstmp]
!   Array of temperature observations that
!   have passed all quality control checks
!  \item[obsrlh]
!   Array of relative humidity observations that
!   have passed all quality control checks
!  \item[obsspd]
!   Array of wind speed observations that
!   have passed all quality control checks
!  \item[obscnt]
!   Number of observations that have passed all
!   quality control checks
!  \item[isize]
!   Max number of observations allowed for 
!   a hemisphere 
!  \item[minwnd]
!    minimum allowable windspeed on the agrmet grid   
!  \item[alert\_number]
!    incremental number of alert message
!  \item[imax]
!    east/west dimension of agrmet grid
!  \item[jmax]
!    north/south dimension of agrmet grid
!  \item[cdmsdir]
!    full path of the JMOBS directory
!  \item[cjulhr]
!   character equivalent of variable julhr
!  \item[hem]
!   Hemisphere returned from lltops routine
!  \item[hemi]
!   Hemisphere (1 = north, 2 = south)
!  \item[ierr1,ierr2,ierr3,istat1]
!   I/O error status
!  \item[ilat]
!   scaled observation latitude (integer)
!  \item[ilon]
!   scaled observation longitude (integer)
!  \item[irecord]
!   loop index for raw observation arrays
!  \item[irelh]
!   scaled relative humidity observations (integer)
!  \item[ispd]
!   Scaled wind speed observations (integer)
!  \item[itmp]
!   Scaled temperature observations (integer)
!  \item[norsou]
!   Character string for northern/southern hemisphere
!  \item[nsize]
!   Number of observations returned from database
!  \item[rlat]
!   descaled observation latitude
!  \item[rlon]
!   descaled observation longitude
!  \item[rrelh]
!   descaled observation relative humidity
!  \item[rrspd]
!   descaled observation wind speed
!  \item[rtmp]
!   descaled observation temperature
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[julhr\_date](\ref{LIS_julhr_date}) \newline
!   convert julian hr to a date format
!  \item[getsfcobsfilename](\ref{getsfcobsfilename}) \newline
!   generates the surface OBS filename 
!  \item[lltops](\ref{lltops}) \newline
!   convert latlon to points on the AGRMET grid
!  \item[LIS\_alert](\ref{LIS_alert}) \newline
!   send an alert message with problems in OBS processing
!  \end{description}
!EOP

  character*100                  :: sfcobsfile
  integer                        :: hemi
  integer, allocatable           :: idpt     ( : )
  integer                        :: ierr1
  integer                        :: ierr2
  integer                        :: ierr3
  integer, allocatable           :: ilat     ( : )
  integer, allocatable           :: ilon     ( : )

  integer                        :: irecord
  integer, allocatable           :: irelh    ( : )
  integer                        :: istat1
  integer, allocatable           :: ispd     ( : )
  integer, allocatable           :: itmp     ( : )
  integer                        :: nsize
  real                           :: rdpt
  real                           :: rigrid
  real                           :: rjgrid
  character*6                    :: cjulhr
  character*10                   :: date10 ! EMK replace cjulhr in log
  character*14, allocatable      :: dtg      ( : )
  character*100                  :: message  ( 20 )
  !character*8, allocatable       :: netyp    ( : )
  character*9, allocatable       :: netyp    ( : )
  character*8                    :: norsou   ( 2 )
!  character*8, allocatable       :: platform ( : )
  !character*9, allocatable       :: platform ( : ) ! EMK BUG FIX
  character*32, allocatable       :: platform ( : ) ! EMK WIGOS
  !character*8, allocatable       :: rptyp    ( : )
  character*9, allocatable       :: rptyp    ( : ) ! WIGOS
  logical, allocatable           :: skip (:) ! EMK
  real                           :: rlat   
  real                           :: rlon   
  real                           :: rrelh
  real                           :: rspd
  real                           :: rtmp
  integer                        :: yr, mo, da, hr
  integer                        :: i 

  ! EMK Patch for DTG check
  character*14 :: yyyymmddhhmmss
  integer :: dtg_year, dtg_month, dtg_day, dtg_hour, dtg_minute, dtg_second
  integer :: diff_year, diff_month, diff_day, diff_hour, diff_minute, &
       diff_second
  type(ESMF_Time) :: dtg_time, lis_time
  type(ESMF_TimeInterval) :: time_interval
  integer :: rc

  real,        external          :: AGRMET_calcrh_dpt
  !character(len=10) :: net10, platform10
  character(len=8) :: netyp8
  character(len=9) :: platform9
  character(len=32) :: net32, platform32
  character(len=8) :: rptyp8
  character*32, parameter :: blank32 = "                                "
  integer :: iunit
  logical :: found_file

  data norsou  / 'NORTHERN', 'SOUTHERN' /

  message = ''

!     ------------------------------------------------------------------
!     Executable code begins here. Intialize observation counter to 0.
!     ------------------------------------------------------------------

  obscnt =  0
  obstmp = -1
  obsrlh = -1
  obsspd = -1
 
  allocate ( idpt     ( isize ) )
  allocate ( ilat     ( isize ) )
  allocate ( ilon     ( isize ) )
  allocate ( irelh    ( isize ) )
  allocate ( ispd     ( isize ) )
  allocate ( itmp     ( isize ) )
  allocate ( dtg      ( isize ) )
  allocate ( netyp    ( isize ) )
  allocate ( platform ( isize ) )
  allocate ( rptyp    ( isize ) )
  allocate ( skip     ( isize ) ) ! EMK

!     ------------------------------------------------------------------
!     Retrieve observations for this julhr and hemi from database.
!     ------------------------------------------------------------------

  call LIS_julhr_date( julhr, yr, mo, da, hr)

  ! EMK Patch for DTG check
  call ESMF_TimeSet(lis_time, yy=yr, mm=mo, dd=da, h=hr, m=0, s=0, rc=rc)

  do hemi=1,2

     ! EMK...For sanity, initialize arrays with missing values
     idpt(:) =  -99999999
     ilat(:) =  -99999999
     ilon(:) =  -99999999
     irelh(:) = -99999999
     ispd(:) =  -99999999
     itmp(:) =  -99999999
     netyp(:) = '-99999999'
     platform(:) = '-99999999'

     skip(:) = .false. ! EMK

     call getsfcobsfilename(sfcobsfile, rootdir, cdmsdir, &
          use_timestamp,hemi, yr, mo, da, hr, use_wigos_sfcobs)

     inquire(file=trim(sfcobsfile), exist=found_file)
     if (.not. found_file) then
        write(LIS_logunit,*) '[WARN] Cannot find ', trim(sfcobsfile)
        message(1) = '[WARN] Program:  LIS'
        message(2) = '  Routine: AGRMET_getsfc'
        message(3) = '  Cannot find file ' // trim(sfcobsfile)
        if (LIS_masterproc) then
           call LIS_alert('LIS.AGRMET_getsfc', &
                alert_number, message)
           alert_number = alert_number + 1
        end if
        message = ''
        if (use_wigos_sfcobs) exit ! These files are global
        cycle
     end if

     write(LIS_logunit,*)'[INFO] Opening: ',trim(sfcobsfile)
     iunit = LIS_getNextUnitNumber()
     open(iunit,file=trim(sfcobsfile),status='old', iostat=ierr1)

     if (ierr1 .ne. 0) then
        write(date10,'(i4, i2.2, i2.2, i2.2)', iostat=istat1) &
             yr, mo, da, hr
        write(LIS_logunit,*)' '
        write(LIS_logunit,*) &
             '[WARN] ROUTINE AGRMET_GETSFC: ERROR OPENING ', &
             trim(sfcobsfile)
        write(LIS_logunit,*)'[WARN] ISTAT IS ', ierr1
        message(1) = 'program:  LIS'
        message(2) = '  routine:  AGRMET_getsfc'
        message(3) = '  error opening file ' // trim(sfcobsfile)
        alert_number = alert_number + 1
        if(LIS_masterproc) then
           call lis_alert( 'LIS.AGRMET_getsfc', alert_number, message )
        endif
        message = ''
        if (use_wigos_sfcobs) exit ! New WIGOS version is global
        cycle
     endif

     if(ierr1.eq.0) then
        read(iunit,*, iostat=ierr2) nsize

        if (ierr2 .ne. 0) then
           write(LIS_logunit,*) &
                '[WARN] Problem reading total report count from ', &
                trim(sfcobsfile)
           message(1) = 'program:  LIS'
           message(2) = '  routine:  AGRMET_getsfc'
           message(3) = '  Cannot read total report count from ' // &
                trim(sfcobsfile)
           alert_number = alert_number + 1
           if(LIS_masterproc) then
              call lis_alert( 'LIS.AGRMET_getsfc', alert_number, &
                   message )
           endif
           message = ''
           close(iunit)
           call LIS_releaseUnitNumber(iunit)
           if (use_wigos_sfcobs) exit
           cycle
        end if

!     ------------------------------------------------------------------
!     If the number of obs in the file is greater than the array size
!     write an alert to the log and set back the number to read to
!     prevent a segfault.
!     ------------------------------------------------------------------

         if ( nsize .GT. isize ) then

           write(LIS_logunit,*)' '
           write(LIS_logunit,*)"******************************************************"
           write(LIS_logunit,*)"* NUMBER OF SURFACE OBSERVATIONS EXCEEDS ARRAY SIZE."
           write(LIS_logunit,*)"* NUMBER OF SURFACE OBS IS ", nsize
           write(LIS_logunit,*)"* ARRAY SIZE IS ", isize
           write(LIS_logunit,*)"* OBSERVATIONS BEYOND ARRAY SIZE WILL BE IGNORED."
           write(LIS_logunit,*)"******************************************************"
           write(LIS_logunit,*)' '

           nsize = isize

        end if

        do i=1,nsize
           if (.not. use_wigos_sfcobs) then
              !read(iunit, *, iostat=ierr3) platform(i), dtg(i), &
              !     itmp(i), idpt(i), irelh(i), ilat(i), ilon(i), &
              !     ispd(i), rptyp(i), netyp(i)
              read(iunit, *, iostat=ierr3) platform9, dtg(i), &
                   itmp(i), idpt(i), irelh(i), ilat(i), ilon(i), &
                   ispd(i), rptyp8, netyp8
              if (ierr3 == 0) then
                 platform(i) = platform9
                 rptyp(i) = rptyp8
                 netyp(i) = netyp8
              end if
           else
              read(iunit, 1000, iostat=ierr3) platform(i), dtg(i), &
                   itmp(i), idpt(i), irelh(i), ilat(i), ilon(i), &
                   ispd(i), rptyp(i), netyp(i)
1000          format(a32, 1x, a14, 1x, i9, 1x, i9, 1x, i9, 1x, &
                   i9, 1x, i9, 1x, i9, 1x, a9, 1x, a9)
           end if

           !if (ierr3 .ne. 0) skip(i) = .true. ! EMK
           ! EMK Patch Skip report if problem occurred reading it
           if (ierr3 .ne. 0) then
              write(LIS_logunit,*) &
                   '[WARN] Problem reading report ', i, &
                   ' from sfcobs file, skipping line'
              skip(i) = .true.
              cycle
           end if

           ! EMK Patch DTG check.  Skip report if more than +/- 30 min from
           ! expected valid time.
           yyyymmddhhmmss = dtg(i)
           read(yyyymmddhhmmss(1:4), '(i4.4)') dtg_year
           read(yyyymmddhhmmss(5:6), '(i2.2)') dtg_month
           read(yyyymmddhhmmss(7:8), '(i2.2)') dtg_day
           read(yyyymmddhhmmss(9:10), '(i2.2)') dtg_hour
           read(yyyymmddhhmmss(11:12), '(i2.2)') dtg_minute
           read(yyyymmddhhmmss(13:14), '(i2.2)') dtg_second
           call ESMF_TimeSet(dtg_time, &
                yy=dtg_year, mm=dtg_month, dd=dtg_day, &
                h=dtg_hour, m=dtg_minute, s=dtg_second, rc=rc)
           time_interval = lis_time - dtg_time
           time_interval = ESMF_TimeIntervalAbsValue(time_interval)
           call ESMF_TimeIntervalGet(time_interval, &
                yy=diff_year, mm=diff_month, d=diff_day, &
                h=diff_hour, m=diff_minute, s=diff_second, rc=rc)
           if (diff_year > 0) then
              skip(i) = .true.
           else if (diff_month > 0) then
              skip(i) = .true.
           else if (diff_day > 0) then
              skip(i) = .true.
           else if (diff_hour > 0) then
              skip(i) = .true.
           else if (diff_minute > 30) then
              skip(i) = .true.
           else if (diff_minute == 30 .and. diff_second > 0) then
              skip(i) = .true.
           end if

        enddo
        close(iunit)
        call LIS_releaseUnitNumber(iunit)

        !if(ierr2.eq.0.and.ierr3.eq.0) then
        ! EMK Patch...Allow observation storage even if problem occurred
        ! reading a particular report
        if (ierr2 .eq. 0) then

!     ------------------------------------------------------------------
!         Retrieve and descale observation latitude and longitude.
!         check to make sure lat/lon is not missing (-99999999) or
!         or unreported (-99999998).  If it is, skip to the next
!         observation.  also ensure there are no goofy values.
!     ------------------------------------------------------------------

           do irecord=1, nsize

              if (skip(irecord)) cycle ! EMK

              if(ilat(irecord).gt.-99999998) then
                 rlat = float(ilat(irecord))/100.0
                 if((rlat.gt.90.0) .or. rlat.lt.-90.0) cycle
              else
                 cycle
              endif
              if (ilon(irecord) .gt. -99999998) then
                 rlon = float(ilon(irecord)) / 100.0
                 if ((rlon .gt. 180.0) .or. (rlon .lt. -180.0)) cycle
              else
                 cycle
              end if
!     ------------------------------------------------------------------
!         Only process the observations for the current hemisphere.
!     ------------------------------------------------------------------

              ! EMK...Use all obs if we are using a WIGOS global file
              if( ((hemi .eq. 1) .and. (rlat .ge. 0.0)) .or. &
                   ((hemi .eq. 2) .and. (rlat .lt. 0.0)) .or. &
                   use_wigos_sfcobs) then
!     ------------------------------------------------------------------
!         Convert point's lat/lon to i/j coordinates.
!         check for values of rigrid and rjgrid outside of
!         the AGRMET model grid, or in the other hemisphere.
!     ------------------------------------------------------------------

                 call latlon_to_ij(LIS_domain(n)%lisproj, rlat, rlon, &
                      rigrid, rjgrid)

                 !if(rigrid.ge.1.and.rigrid.le.imax.and. &
!                      rjgrid.ge.1.and.rjgrid.le.jmax) then 
! EMK TEST
                 if (.true.) then
!     ------------------------------------------------------------------
!         Make sure the observation is actually reporting temp, wind, 
!         and rh or dew point before we store it.  
!         Initial tests in july 99 showed that, on average, 5% of the 
!         observations had to be discared for this reason.
!     ------------------------------------------------------------------

                    if ( ( (itmp(irecord)  .gt. -99999998) .and. &
                         (ispd(irecord)  .gt. -99999998) ) .and. &
                         ( (irelh(irecord) .gt. -99999998) .or. &
                         (idpt(irecord) .gt. -99999998) ) ) then

                       ! EMK...Make sure dew point .le. temperature
                       if ( (itmp(irecord)  .gt. -99999998) .and. &
                            (idpt(irecord) .gt. -99999998) ) then
                          if (idpt(irecord) .gt. itmp(irecord)) cycle
                       end if

!     ------------------------------------------------------------------
!         Gross error check and set temperature.
!         Ensure temperature values are between 200.0 and 350.0 K.
!         a "-99999999" indicates missing data and a "-99999998" 
!         indicates a non report.  Bad data is given a value of -1, 
!         so the Barnes scheme will ignore it.
!     ------------------------------------------------------------------

                       if (itmp(irecord) .gt. -99999998) then
                          rtmp = float (itmp(irecord)) / 100.0
                          if ( (rtmp .lt. 200.0) .or. (rtmp .gt. 350.0) ) then
                             rtmp = -1.0
                          end if
                       else
                          rtmp = -1.0
                       end if
        
!     ------------------------------------------------------------------
!         Gross error check and set relative humidity values.
!         ensure RH values are between 0.0 and 1.0.
!         JMOBS RH is in percent.  Missing or bad values are given a
!         value of -1, so the Barnes scheme will ignore it.  
!     ------------------------------------------------------------------

                       if (irelh(irecord) .gt. -99999998) then
                          rrelh = float(irelh(irecord)) / 10000.0
                          if ( (rrelh .gt. 1.0) .or. (rrelh .lt. 0.0) ) then   
                             rrelh = -1.0
                          end if
                       else if ( (idpt(irecord) .gt. -99999998) .and. &
                            (rtmp .ne. -1.0) ) then
                          rdpt = float(idpt(irecord)) / 100.0
                          rrelh = AGRMET_calcrh_dpt(rtmp, rdpt)
                          if ( (rrelh .gt. 1.0) .or. (rrelh .lt. 0.0) ) then   
                             rrelh = -1.0
                          end if
                       else
                          rrelh = -1.0
                       end if
           
!     ------------------------------------------------------------------
!         Constrain surface wind speed to a value between
!         75 ms-1 and minwind (from control file).
!         JMOBS winds are scaled by 10.  Bad values are given a value
!         of -1, so the Barnes scheme will ignore it.
!     ------------------------------------------------------------------

                       if (ispd(irecord) .gt. -99999998) then
                          rspd = float(ispd(irecord)) / 10.0
                          rspd = max( min( rspd, 75.0 ), minwnd )
                       else
                          rspd = -1.0
                       end if
                       
!     ------------------------------------------------------------------
!         If all the data is out of range and/or missing, don't store
!         this observation.
!     ------------------------------------------------------------------
           
                       if ( (nint(rtmp)  .gt. 0)  .and. &
                            (nint(rrelh) .gt. 0)  .and.&
                            (nint(rspd)  .gt. 0) ) then 
                          
!     ------------------------------------------------------------------
!         If we make it here, the observation probably has some good
!         data, so increment the observation counter, store the
!         data into the observation arrays.
!     ------------------------------------------------------------------

                          ! EMK...Add to data structures.  Handle reformated
                          ! CDMS data that is missing platform and network
                          if (rtmp .gt. 0) then
                             net32 = blank32
                             net32 = netyp(irecord)
                             platform32 = blank32
                             platform32 = platform(irecord)
                             if (trim(net32) .eq. 'NULL') then
                                net32 = 'CDMS'
                             end if
                             if (trim(platform32) .eq. '-99999999') then
                                platform32 = '00000000'
                             end if
                             call USAF_assignObsData(t2mObs,net32, &
                                  platform32,rtmp,rlat,rlon, &
                                  agrmet_struc(n)%bratseth_t2m_stn_sigma_o_sqr,&
                                  0.)

                          end if
                          if (rrelh .gt. 0) then
                             net32 = blank32
                             net32 = netyp(irecord)
                             platform32 = blank32
                             platform32 = platform(irecord)
                             if (trim(net32) .eq. 'NULL') then
                                net32 = 'CDMS'
                             end if
                             if (trim(platform32) .eq. '-99999999') then
                                platform32 = '00000000'
                             end if
                             call USAF_assignObsData(rh2mObs,net32, &
                                  platform32,rrelh,rlat,rlon, &
                                  agrmet_struc(n)%bratseth_t2m_stn_sigma_o_sqr, &
                                  0.)
                          end if
                          if (rspd .gt. 0) then
                             net32 = blank32
                             net32 = netyp(irecord)
                             platform32 = blank32
                             platform32 = platform(irecord)
                             if (trim(net32) .eq. 'NULL') then
                                net32 = 'CDMS'
                             end if
                             if (trim(platform32) .eq. '-99999999') then
                                platform32 = '00000000'
                             end if
                             call USAF_assignObsData(spd10mObs,net32, &
                                  platform32,rspd,rlat,rlon, &
                                  agrmet_struc(n)%bratseth_spd10m_stn_sigma_o_sqr, &
                                  0.)
                          end if

                          obscnt         = obscnt + 1

                          ri(obscnt)     = rigrid
                          rj(obscnt)     = rjgrid
                          obstmp(obscnt) = rtmp
                          obsrlh(obscnt) = rrelh
                          obsspd(obscnt) = rspd
                       endif
!     ------------------------------------------------------------------
!         If we have reached the number of obs that our hard-wired
!         arrays can hold, then exit the loop.
!     ------------------------------------------------------------------

                       if ( obscnt .eq. isize ) then
                          write(LIS_logunit,*) &
                               '[WARN] ROUTINE GETSFC: REACHED MAXIMUM NUMBER OF SFC OBS TO SAVE IN MEMORY'
                          message(1) = 'program:  LIS'
                          message(2) = '  routine:  AGRMET_getsfc'
                          message(3) = '  reached maximum number of sfc obs to save memory in memory'
                          alert_number = alert_number + 1
                          if(LIS_masterproc) then
                             call lis_alert( 'LIS.AGRMET_getsfc', &
                                  alert_number, &
                                  message )
                          endif
                          message = ''
                          exit ! Jump out of hemi do loop
                       end if
                    endif
                 endif
              endif
           enddo
        endif
     endif

     if (use_wigos_sfcobs) exit ! New sfcobs file is global
  enddo

  deallocate ( idpt )
  deallocate ( ilat )
  deallocate ( ilon )
  deallocate ( irelh )
  deallocate ( ispd )
  deallocate ( itmp )
  deallocate ( dtg )
  deallocate ( netyp )
  deallocate ( platform )
  deallocate ( rptyp )

  deallocate(skip) ! EMK

  write(LIS_logunit,*)'[INFO] Stored ', obscnt, ' sfc obs'
end subroutine AGRMET_getsfc

!BOP
! 
! !ROUTINE: getsfcobsfilename
! \label{getsfcobsfilename}
!
! !INTERFACE: 
subroutine getsfcobsfilename(filename, rootdir, dir, &
     use_timestamp, hemi, yr, mo, da, hr, use_wigos_sfcobs)

  implicit none

! !ARGUMENTS: 
  character(len=*), intent(inout) :: filename
  character(len=*), intent(in) :: rootdir
  character(len=*), intent(in) :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi
  integer, intent(in) :: yr, mo, da, hr
  logical, intent(in) :: use_wigos_sfcobs
!
! !DESCRIPTION: 
!  This routines generates the name of the surface obs file
! 
!  The arguments are: 
!  \begin{description}
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[use\_timestamp]
!    flag to indicate whether the directories 
!    should be timestamped or not
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[filename]
!    created filename
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP
  character*2 :: fhr
  character*4 :: fyr
  character*2 :: fmo
  character*2 :: fda
  character*2 :: fhemi 

  write(unit=fhr,fmt='(i2.2)') hr
  write(unit=fyr,fmt='(i4.4)') yr
  write(unit=fmo,fmt='(i2.2)') mo
  write(unit=fda,fmt='(i2.2)') da
  if(hemi.eq.1) then 
     write(unit=fhemi,fmt='(a2)') 'nh'
  else
     write(unit=fhemi,fmt='(a2)') 'sh'
  endif

  if(use_timestamp.eq.1) then
     if (use_wigos_sfcobs) then
        filename = trim(rootdir) // '/' // &
             trim(fyr) // trim(fmo) // trim(fda) // &
             '/' // trim(dir) // '/sfcobs_01hr_' // &
             trim(fyr) // trim(fmo) // trim(fda) // trim(fhr) // ".txt"
     else
        filename = trim(rootdir) // '/' // &
             trim(fyr) // trim(fmo) // trim(fda) // &
             '/' // trim(dir) // '/sfcobs_' // &
             trim(fhemi) // '.01hr.' // &
             trim(fyr) // trim(fmo) // trim(fda) // trim(fhr)
     end if
  else
     if (use_wigos_sfcobs) then
        filename = trim(rootdir) // &
             '/' // trim(dir) // '/sfcobs_01hr_' // &
             trim(fyr) // trim(fmo) // trim(fda) // trim(fhr) // ".txt"
     else
        filename = trim(rootdir) // &
             '/' // trim(dir) // '/sfcobs_' // &
             trim(fhemi) // '.01hr.' // &
             trim(fyr) // trim(fmo) // trim(fda) // trim(fhr)
     end if
  endif
end subroutine getsfcobsfilename

