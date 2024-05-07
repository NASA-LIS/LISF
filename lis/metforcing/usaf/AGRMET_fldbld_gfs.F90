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
! !ROUTINE: AGRMET_fldbld_gfs
! \label{AGRMET_fldbld_gfs}
!
! !REVISION HISTORY:
! 12 Aug 1988  Initial version.............................MSgt Neill/SDDC
! 04 Aug 1999  Ported to IBM SP-2. Major changes to read first guess data
!              directly from GRIB files and process AVN as well as NOGAPS.
!              Simplified driver by placing read of first guess data,
!              interpolation of data to AGRMET grid, and writing of the
!              interpolated data into separate subroutines.  Made grid-
!              specific variables dynamically allocatable....Mr Gayno/DNXM
! 29 Jul 2005  Incorporated into LIS......................Sujay Kumar/NASA
! 26 Dec 2007  Simplified filename creation.............Marv Freimund/DNXM
! 12 Aug 2008  Added search for GFS data up to 24 hours before
!              cycle time......................Chris Franks/2WXG/WEA(SEMS)
! 04 Mar 2009  Added ability to add timestamp only to GFS
!              directories.......................Ryan Ruhge/16WS/WXE(SEMS)
! 11 Mar 2010  Changed program names in messages to LIS.
!              ................................Chris Franks/16WS/WXE(SEMS)
! 25 Jan 2012  Changed the grib reader to grib api library Sujay Kumar/NASA
! 14 Jun 2016  Renamed from AGRMET_fldbld James Geiger/NASA
! 09 Jun 2017  Refactor interpolation James Geiger/NASA
! 14 Jun 2017  Added GFS terrain height and 2-m T, RH; use of total LIS
!              domain..........................Eric Kemp/GSFC
! 16 Oct 2017  Changed logic to roll-back to earlier GFS cycle if file is
!              missing or incomplete...........Eric Kemp/GSFC
! 05 Mar 2020  Added support for new GFS filename convention.
!              ................................Eric Kemp/GSFC
! !INTERFACE:    
subroutine AGRMET_fldbld_gfs(n,order,julhr,rc)
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
  integer, intent(in) :: julhr
  integer, intent(out) :: rc
!
! !DESCRIPTION: 
!  This routine interpolates the GFS first guess height, temperature, moisture, 
!  and wind data to the AGRMET grid. 
!
! The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[avnfile]
!    name of the first guess file
!  \item[nunit]
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
!  \item[file\_julhr]
!    julian hour used to determine names of forecast files from previous cycles
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!    converts the julian hour to a date format
!  \item[getAVNfilename](\ref{getAVNfilename}) \newline
!    generates the first guess AVN filename
!  \item[AGRMET\_fldbld\_read\_gfs](\ref{AGRMET_fldbld_read_gfs}) \newline
!   read AVN or NOGAPS data in grib format
!  \item[AGRMET\_fg2lis](\ref{AGRMET_fg2lis}) \newline
!   interpolate first guess data to the LIS grid
!  \item[lis\_abort](\ref{LIS_abort}) \newline
!   abort in case of error
!  \end{description}
!EOP
  integer                 :: ftn, igrib
  !character*120           :: avnfile
  character*255           :: avnfile ! EMK
  integer                 :: yr1, mo1, da1, hr1

  integer                 :: nunit
  integer                 :: ksec2       ( 10 )
  character*255           :: message     ( 20 )
  integer                 :: iginfo      ( 40 )
  real                    :: ginfo       ( 40 )
  real                    :: gridres
  integer                 :: alert_number
  real, allocatable       :: fg_hgt     ( : , : , : )
  real, allocatable       :: fg_rh      ( : , : , : )
  real, allocatable       :: fg_tmp     ( : , : , : )
  real, allocatable       :: fg_hgt_sfc ( : , : )
  real, allocatable       :: fg_rh_sfc  ( : , : )
  real, allocatable       :: fg_tmp_sfc ( : , : )
  real, allocatable       :: fg_wspd    ( : , : )
  real, allocatable       :: fg_pres    ( : , : )
  integer                 :: ifguess, jfguess
  integer                 :: kprs 
  integer                 :: prslvls      (30) 
  integer                 :: ksec1       ( 100 )
  integer                 :: avnflag
  integer                 :: center
  integer                 :: ierr
  logical*1               :: found
  integer                 :: yr_2d
  integer                 :: file_julhr
  character*100           :: gtype
  integer                 :: fc_hr
  integer                 :: dataDate, dataTime
  integer                 :: i
  integer :: rc2
  logical :: first_time, found_inq

  ! FIXME...This should be moved into a module for all of LIS to 
  ! reference.
  data prslvls / 1000,975,950,925,900,850,800,750,700,650,600,550,500,&
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
  kprs = 13

  rc = 0 

  ! Will search GFS cycles every six hours, up to 24 hours,
  ! until we find an acceptable file.
  call LIS_julhr_date(julhr,yr1,mo1,da1,hr1)
  fc_hr = 0           
  file_julhr = julhr  
  found = .FALSE. 
  first_time = .true.
  do while (.not. found)

     ! Allow use of current GFS cycle for legacy runs.
     if (.not. first_time) then
        fc_hr = fc_hr + 6
        if (fc_hr .gt. 24) exit ! Give up
        file_julhr = file_julhr - 6
        call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)
     end if
     first_time = .false.
     
     yr_2d = mod(yr1,100)
     if(yr_2d.eq.0) yr_2d = 100 
     call getAVNfilename(avnfile, agrmet_struc(n)%agrmetdir,&
          agrmet_struc(n)%gfsdir, agrmet_struc(n)%use_timestamp,&
          agrmet_struc(n)%gfs_timestamp, &
          agrmet_struc(n)%gfs_filename_version, &
          yr1, mo1, da1, hr1, fc_hr)

     ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
     ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
     ! writing messages to stdout/stderr, which may lead to runtime problems.
     inquire(file=trim(avnfile),exist=found_inq)
     if (.not. found_inq) then
        write(LIS_logunit,*) '[WARN] Cannot find file '//trim(avnfile)
        cycle
     end if

     ! Open first guess grib data using library utility.  just read
     ! the first file only, as all data will be of the same type
     ! (avn or nogaps) because the search script ensures that it is.

#if (defined USE_GRIBAPI) 
     
     call grib_open_file(ftn,trim(avnfile),'r',ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] Failed to open - '//trim(avnfile)
        cycle
     end if
     
     ! Read in the first grib record, unpack the header and extract
     ! section 1 and section 2 information.
     call grib_new_from_file(ftn,igrib,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] failed to read - '//trim(avnfile)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'centre',center,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
                '[WARN] in grib_get: centre in AGRMET_fldbld_gfs'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif
        
     call grib_get(igrib,'gridType',gtype,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
                '[WARN] in grid_get: gridtype in AGRMET_fldbld_gfs'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     if(trim(gtype).ne."regular_ll") then  
        write(LIS_logunit,*)'[WARN] GFS file not on lat/lon grid!'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'Ni',iginfo(1),ierr)
     if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) '[WARN] in grid_get:Ni in AGRMET_fldbld_gfs'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'Nj',iginfo(2),ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get:Nj in AGRMET_fldbld_gfs'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grid_get:jDirectionIncrementInDegrees in ' // &
             'AGRMET_fldbld_gfs'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'dataDate',dataDate,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
                '[WARN] in grid_get:dataDate in AGRMET_fldbld_gfs'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'dataTime',dataTime,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grid_get:dataTime in AGRMET_fldbld_gfs'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     if ( yr1*10000+mo1*100+da1 .ne. dataDate) then
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     end if


     if (hr1*100 .ne. dataTime ) then
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     end if

     ! Here we tentatively have a file we can use.   Close it for now, and
     ! prepare to pull the appropriate variables.
     call grib_release(igrib,ierr)
     call grib_close_file(ftn)

     write(LIS_logunit,*)'[INFO] FIRST GUESS DATA IS ON A ', gridres,&
          ' DEGREE LAT/LON GRID'
     ifguess = iginfo(1)
     jfguess = iginfo(2)

     ! FIXME...Reject file if wrong center?
     if (center .eq. 58) then
        write(LIS_logunit,*)'- FIRST GUESS DATA IS FROM NOGAPS MODEL'
     elseif (center .eq. 7) then
        write(LIS_logunit,*)'- FIRST GUESS DATA IS FROM GFS MODEL'
     else
        write(LIS_logunit,*)'[WARN] UNKNOWN SOURCE FOR FIRST GUESS DATA'
     end if

     ! Allocate memory to store GRIB data
     allocate ( fg_hgt  (ifguess, jfguess, kprs) )
     allocate ( fg_rh   (ifguess, jfguess, kprs) )
     allocate ( fg_tmp  (ifguess, jfguess, kprs) )
     allocate ( fg_hgt_sfc  (ifguess, jfguess) )
     allocate ( fg_rh_sfc   (ifguess, jfguess) )
     allocate ( fg_tmp_sfc  (ifguess, jfguess) )
     allocate ( fg_wspd (ifguess, jfguess) )
     allocate ( fg_pres (ifguess, jfguess) )

     alert_number = 0 
     
     call AGRMET_fldbld_read_gfs( avnfile, ifguess, jfguess,&
          agrmet_struc(n)%wndwgt, agrmet_struc(n)%minwnd, &
          fg_hgt, fg_rh, fg_tmp, fg_hgt_sfc, fg_rh_sfc, fg_tmp_sfc, &
          fg_wspd, fg_pres, kprs, prslvls, alert_number, rc2 )

     ! Check for error
     if (rc2 .ne. 0) then
        deallocate (fg_hgt)
        deallocate (fg_rh)
        deallocate (fg_tmp)
        deallocate (fg_hgt_sfc)
        deallocate (fg_rh_sfc)
        deallocate (fg_tmp_sfc)
        deallocate (fg_wspd)
        deallocate (fg_pres)
        cycle
     end if

     ! At this point, we have the data from the GRIB file.  Interpolate
     ! and wrap up
     if(order.eq.1) then 
        agrmet_struc(n)%agr_bgrd_src_c = "GFS"
        do i = 1, kprs
           call AGRMET_fg2lis(n, ifguess, jfguess, &
                          fg_hgt(:,:,i), agrmet_struc(n)%agr_hgt_c(:,:,i))
           call AGRMET_fg2lis(n, ifguess, jfguess, &
                          fg_rh(:,:,i), agrmet_struc(n)%agr_rh_c(:,:,i))
           call AGRMET_fg2lis(n, ifguess, jfguess, &
                          fg_tmp(:,:,i), agrmet_struc(n)%agr_tmp_c(:,:,i))
        enddo ! i
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_wspd, agrmet_struc(n)%agr_wspd_c)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_pres, agrmet_struc(n)%agr_pres_c)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_hgt_sfc, agrmet_struc(n)%agr_hgt_sfc_c)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_rh_sfc, agrmet_struc(n)%agr_rh_sfc_c)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_tmp_sfc, agrmet_struc(n)%agr_tmp_sfc_c)
     else
        agrmet_struc(n)%agr_bgrd_src_p = "GFS"
        do i = 1, kprs
           call AGRMET_fg2lis(n, ifguess, jfguess, &
                          fg_hgt(:,:,i), agrmet_struc(n)%agr_hgt_p(:,:,i))
           call AGRMET_fg2lis(n, ifguess, jfguess, &
                          fg_rh(:,:,i), agrmet_struc(n)%agr_rh_p(:,:,i))
           call AGRMET_fg2lis(n, ifguess, jfguess, &
                          fg_tmp(:,:,i), agrmet_struc(n)%agr_tmp_p(:,:,i))
        enddo ! i
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_wspd, agrmet_struc(n)%agr_wspd_p)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_pres, agrmet_struc(n)%agr_pres_p)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_hgt_sfc, agrmet_struc(n)%agr_hgt_sfc_p)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_rh_sfc, agrmet_struc(n)%agr_rh_sfc_p)
        call AGRMET_fg2lis(n, ifguess, jfguess, &
                           fg_tmp_sfc, agrmet_struc(n)%agr_tmp_sfc_p)
     endif
     
     ! Clean up
     deallocate ( fg_hgt  )
     deallocate ( fg_rh   )
     deallocate ( fg_tmp  )
     deallocate ( fg_hgt_sfc  )
     deallocate ( fg_rh_sfc   )
     deallocate ( fg_tmp_sfc  )
     deallocate ( fg_wspd )         
     deallocate ( fg_pres )

     found = .true.
     exit

  end do
#endif     

  ! Give up if no acceptable GFS files were found
  if (.not. found) then
     write(LIS_logunit,*)'[WARN] No matching GFS file found!'
     rc = 1
     return
  end if

  write(LIS_logunit,*) &
       '[INFO] Using NWP fields from ',trim(avnfile)
  rc = 0

end subroutine AGRMET_fldbld_gfs

!BOP
! 
! !ROUTINE: getAVNfilename
! \label{getAVNfilename}
!
! !INTERFACE: 
subroutine getAVNfilename(filename, rootdir, dir, &
     use_timestamp, gfs_timestamp, gfs_filename_version, &
     yr, mo, da, hr, fc_hr)

   use LIS_logMod, only: LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS: 
  character(*)        :: filename
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: gfs_timestamp
  integer, intent(in) :: gfs_filename_version
  integer, intent(in) :: yr, mo, da, hr
  integer, intent(in) :: fc_hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the timestamped AVN file 
! 
!  The arguments are: 
!  \begin{description}
!   \item[dir]
!    full path to the directory containing the data
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
!   \item[fc\_hr]
!    hours from reference time to valid time
!  \end{description}
!EOP
  character(10) :: ftime1
  character( 2) :: fhr
  character( 4) :: fchr
  character(255) :: basename
  character(8) :: cyyyymmdd
  character(3) :: chhh

  write (UNIT=fhr, FMT='(i2.2)') hr
  write (UNIT=fchr, FMT='(i4.4)') fc_hr

  select case (gfs_filename_version)
  case (1)
     ! Longtime convention used by AGRMET and LIS
     basename = 'MT.avn_CY.' // &
          fhr // '_fh.'//fchr//'_tl.press_gr.0p5deg'
  case (2)
     ! New convention to be adopted in 2020
     write(cyyyymmdd, '(I4.4,I2.2,I2.2)') yr, mo, da
     write(chhh, '(I3.3)') fc_hr
     basename = 'PS.NCEP_SC.U_DI.A_DC.GRID_GP.GFS_SP.SIMPLE_GR' // &
          '.C0P5DEG_AR.GLOBAL_PA.GFS_DD.' // cyyyymmdd // '_CY' // &
          '.' // fhr // '_FH.' // chhh // '_DF.GR2'
  case default
     write(LIS_logunit,*) &
          '[ERR] Invalid GFS filename version!'
     write(LIS_logunit,*) &
          '[ERR] Available options are 1 and 2'
     write(LIS_logunit,*) &
          '[ERR] Found ', gfs_filename_version
     call LIS_endrun()
  end select

  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     filename = trim(rootdir) // ftime1 // trim(dir) // &
          "/" // trim(basename)
  elseif (gfs_timestamp .eq. 1) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/' 
     filename = trim(rootdir) // trim(dir) // '/' // ftime1 // &
          "/" // trim(basename)
  else
     filename = trim(rootdir) //  '/'  // trim(dir) // &
          "/" // trim(basename)
  endif

end subroutine getAVNfilename


!BOP
! 
! !ROUTINE: getDataLevelsFilename
! \label{getDataLevelsFilename}
!
! !INTERFACE: 
subroutine getDataLevelsFilename(filename, dir, hr)

  implicit none
! !ARGUMENTS: 
  character(*)        :: filename
  character(*)        :: dir
  integer, intent(in) :: hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the timestamped AVN file 
! 
!  The arguments are: 
!  \begin{description}
!   \item[filename]
!    created filename
!   \item[dir]
!    full path to the directory containing the data
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP
  character(2) :: fhr

  write (UNIT=fhr, FMT='(i2.2)') hr

  filename = trim(dir) // '/data_levels.' // trim(fhr) // 'z'

end subroutine getDataLevelsFilename


!BOP
!
! !ROUTINE: AGRMET_fldbld_read_gfs
!  \label{AGRMET_fldbld_read_gfs}
!
! !REVISION HISTORY:
!    
!     05 aug 99  initial version...........................mr gayno/dnxm
!     31 mar 00  changed variable "fg_filename" to character*120 to
!                be compatable with routine copen.c........mr gayno/dnxm
!     07 jan 02  modified for new avn file naming convention............
!                ..........................................mr gayno/dnxm
!     01 aug 02  removed hard-wired references to isobaric levels of
!                first guess source........................mr gayno/dnxm
!     20 aug 05 Sujay Kumar, Initial Specification in LIS
!     11 Mar 10 Changed program names in messages to LIS.  Left argument
!               to LIS_alert with old program name to keep alert numbers
!               unique........................Chris Franks/16WS/WXE/SEMS
!     17 Mar 10 Account for GFS changing the GRIB time range indicator 
!               for 0 hour files to 1.........Chris Franks/16WS/WXE/SEMS
!     25 Jan 12 Switched to grib api for reading AVN/GFS files .........
!               ........................................Sujay Kumar/NASA
!   14 Jun 2016 Renamed from AGRMET_fldbld_read James Geiger/NASA
!   14 Jun 2017 Added GFS ground height, 2-m T, RH; also added dynamic
!               GRIB2 support...........................Eric Kemp/GSFC
!   28 Jul 2017 Added standard out messages clarifying if a GFS GRIB file
!               could not be successfully opened........Eric Kemp/GSFC
!   16 Oct 2017 Refactored to return if a read error occurs of if not all
!               fields are available. Return code now included in argument
!               list....................................Eric Kemp/GSFC
! 
! !INTERFACE: 
subroutine AGRMET_fldbld_read_gfs( fg_filename, ifguess, jfguess,&
     wndwgt, minwnd, fg_hgt, fg_rh,&
     fg_tmp, fg_hgt_sfc, fg_rh_sfc, fg_tmp_sfc, fg_wspd, fg_pres, &
     kprs, prslvls, alert_number, rc )
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
  real,           intent(in)    :: minwnd
  real,           intent(in)    :: wndwgt
  integer,        intent(in)    :: kprs
  real,           intent(out)   :: fg_hgt      ( ifguess,jfguess,kprs )
  real,           intent(out)   :: fg_rh       ( ifguess,jfguess,kprs )
  real,           intent(out)   :: fg_tmp      ( ifguess,jfguess,kprs )
  real,           intent(out)   :: fg_hgt_sfc  ( ifguess,jfguess)
  real,           intent(out)   :: fg_rh_sfc   ( ifguess,jfguess)
  real,           intent(out)   :: fg_tmp_sfc  ( ifguess,jfguess)

  real,           intent(out)   :: fg_wspd     ( ifguess,jfguess )
  real,           intent(out)   :: fg_pres     ( ifguess,jfguess )
  integer,        intent(in)    :: prslvls     ( 30 )
  integer,        intent(inout) :: alert_number
  integer, intent(out) :: rc
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
!     - check forecast hour of data, if it is not the analysis \newline
!       or zero hour forecast, then send an alert message \newline
!       to warn of possible degradation. \newline
!     - check counter variables.  if data is missing abort. \newline
!     - if nogaps, convert dewpoint depression to relative \newline
!       humidity. \newline
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
!      \item[wndwgt]
!        adjustment factor for first guess
!        winds.
!      \item[minwnd]
!        minimum allowable wind speed on the
!        agrmet grid 
!      \item[fg\_hgt]
!        array of first guess isobaric heights
!      \item[fg\_rh]
!        array of first guess isobaric relative humidity
!      \item[fg\_tmp]
!        array of first guess isobaric temperatures
!      \item[fg\_wspd]
!        array of first guess surface wind speeds
!      \item[kprs]
!        number of isobaric levels where fldbld
!        needs first guess data
!       \item[prslvls]
!        the isobaric levels where fldbld needs
!        first guess data
!       \item[alert\_number]
!        counts number of alert messages sent
!       \item[cstat]
!        I/O status, character
!       \item[message]
!        Error message
!       \item[count\_dpd]
!        counts number of isobaric dewpoint 
!        depression levels read in from first guess file
!       \item[count\_hgt]
!        counts number of isobaric height levels
!        read in from first guess file
!       \item[count\_rh]
!        counts number of isobaric relative humidity levels
!        read in from first guess file
!       \item[count\_tmp]
!        counts number of isobaric temperature levels
!        read in from first guess file
!       \item[count\_uwnd]
!        counts number of isobaric u wind levels
!        read in from first guess file
!       \item[count\_vwnd]
!        counts number of isobaric v wind levels
!        read in from first guess file
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
!       \item[AGRMET\_calcrh]
!        function to calculate relative humidity
!       \item[dum1d]
!        dummy array
!       \item[fg\_dpd]
!        array of first guess dewpoint depressions
!       \item[fg\_uwnd]
!        array of first guess surface u winds
!       \item[fg\_vwnd]
!        array of first guess surface v winds
!       \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_calcrh](\ref{AGRMET_calcrh}) \newline
!   calculate relative humdities
!  \end{description}
!EOP
  character*9                   :: cstat
  character*255                 :: message     ( 20 )
  integer                       :: count_dpd
  integer                       :: count_hgt
  integer                       :: count_rh
  integer                       :: count_tmp
!  integer                       :: count_uwnd
!  integer                       :: count_vwnd
  integer :: count_hgt_sfc
  integer :: count_rh_sfc
  integer :: count_tmp_sfc
  integer :: count_uwnd_sfc
  integer :: count_vwnd_sfc
  integer :: count_pres_sfc
  integer                       :: file_age
  integer                       :: i
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: j,c,r
  integer                       :: igrib
  integer                       :: ftn
  integer                       :: kk, nvars
  integer                       :: pds7_val, pds8_val, pds9_val
  integer                       :: k
  integer                       :: ksec1       ( 100 )
  integer                       :: nunit
  
  real,           external      :: AGRMET_calcrh
  real,           allocatable   :: dum1d       ( : )
  real,           allocatable   :: fg_dpd      ( : , : , : )
  real,           allocatable   :: fg_uwnd     ( : , : )
  real,           allocatable   :: fg_vwnd     ( : , : )
  ! EMK...For GRIB 2
  integer :: editionNumber_val
  integer :: param_disc_val, param_cat_val, &
       param_num_val
  integer :: typeOfFirstFixedSurface_val, level_val
  character(len=4) :: check_gfs_grib1_message
  character(len=4) :: check_gfs_grib2_message
  character(len=4) :: grib_msg
  logical :: found
  integer :: plevel
  integer, external :: set_plevel
  logical :: found_inq

  ! Executable code begins here ...
  rc = 0

  ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
  ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
  ! writing messages to stdout/stderr, which may lead to runtime problems.
  inquire(file=trim(fg_filename),exist=found_inq)
  if (.not. found_inq) then
     write(LIS_logunit,*) '[WARN] Cannot find file '//trim(fg_filename)
     rc = 1
     return
  end if

#if (defined USE_GRIBAPI)

  ! If a problem occurs here, we can just return immediately since no
  ! memory has been allocated yet.
  call grib_open_file(ftn,trim(fg_filename),'r',ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] Failed to open - '//trim(fg_filename)
     rc = 1
     return
  end if

  allocate ( dum1d   (ifguess*jfguess) )
  allocate ( fg_dpd  (ifguess, jfguess, kprs) )
  allocate ( fg_uwnd (ifguess, jfguess) )
  allocate ( fg_vwnd (ifguess, jfguess) )
       
  ! From this point, we must deallocate memory before returning. 
  ! Unfortunately this means using a GOTO statement if a problem is
  ! encountered, but such is life.
  count_dpd  = 0
  count_hgt  = 0
  count_rh   = 0
  count_tmp  = 0
!  count_uwnd = 0
!  count_vwnd = 0
  count_hgt_sfc = 0
  count_rh_sfc = 0
  count_tmp_sfc = 0
  count_uwnd_sfc = 0
  count_vwnd_sfc = 0
  count_pres_sfc = 0

  write(LIS_logunit,*)' '
  write(LIS_logunit,*)'[INFO] READING ', trim(fg_filename)

  call grib_count_in_file(ftn,nvars,ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] in grib_count_in_file in ' // &
          'AGRMET_fldbld_read_gfs'
     goto 100
  end if

  ! Tentatively loop through every field in GRIB file looking for the variables
  ! we want. The code below will exit the loop early if a problem is found *or*
  ! once all the required variables are found and read in.

  do kk=1,nvars
 
     call grib_new_from_file(ftn,igrib,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] failed to read - '//trim(fg_filename)
        goto 100
     end if

     ! See which version of GRIB this is
     call grib_get(igrib,'editionNumber',editionNumber_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grib_get: editionNumber in ' // &
             'AGRMET_fldbld_read_gfs'
        call grib_release(igrib,ierr)
        goto 100
     end if

     ! See if we have one of the right fields. The GRIB parameters to be 
     ! checked depend on if the file is GRIB1 or GRIB2
     pds7_val = 0
     pds8_val = 0
     pds9_val = 0
     param_disc_val = 0
     param_cat_val = 0
     param_num_val = 0
     typeOfFirstFixedSurface_val = 0
     level_val = 0          
     if (editionNumber_val .eq. 1) then ! GRIB1
        call grib_get(igrib,'indicatorOfParameter',pds7_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: indicatorOfParameter in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        
        call grib_get(igrib,'level',pds9_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: level in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        
        call grib_get(igrib,'indicatorOfTypeOfLevel',pds8_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: indicatorOfTypeOfLevel in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if

        grib_msg = check_gfs_grib1_message(editionNumber_val,pds7_val, &
             pds8_val,pds9_val)
        
     else if (editionNumber_val .eq. 2) then ! GRIB2
        call grib_get(igrib,'discipline',param_disc_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: discipline in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        
        call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: parameterCategory in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        
        call grib_get(igrib,'parameterNumber',param_num_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: parameterNumber in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        
        call grib_get(igrib,'typeOfFirstFixedSurface', &
             typeOfFirstFixedSurface_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: typeOfFirstFixedSurface in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        
        call grib_get(igrib,'level',level_val,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: level in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        
        grib_msg = check_gfs_grib2_message(editionNumber_val,param_disc_val, &
             param_cat_val,param_num_val,typeOfFirstFixedSurface_val,level_val)

     else
	write(LIS_logunit,*)'[WARN] Unknown GRIB edition number ',editionNumber_val
	call grib_release(igrib,ierr)
	goto 100
     end if ! editionNumber_val

     ! Handle surface pressure
     if (trim(grib_msg) .eq. 'sp') then
        call grib_get(igrib, 'values',dum1d,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: values in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        do r=1,jfguess
           do c=1,ifguess
              fg_pres (c,r) = dum1d(c+(r-1)*ifguess)
           enddo ! c
        enddo ! r
        count_pres_sfc = count_pres_sfc + 1
     endif
     
     ! Handle 2-m temperature
     if (trim(grib_msg) .eq. '2t') then
        call grib_get(igrib, 'values',dum1d,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: values in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        do r=1,jfguess
           do c=1,ifguess
              fg_tmp_sfc (c,r) = dum1d(c+(r-1)*ifguess)
           enddo ! c
        enddo ! r
        count_tmp_sfc = count_tmp_sfc + 1
     end if
     
     ! Handle 2-m RH
     if (trim(grib_msg) .eq. '2rh') then
        call grib_get(igrib, 'values',dum1d,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: values in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        do r=1,jfguess
           do c=1,ifguess
              fg_rh_sfc(c,r) = dum1d(c+(r-1)*ifguess)*0.01
           enddo ! c
        enddo ! r
        count_rh_sfc = count_rh_sfc + 1
     end if
     
     ! Handle surface height
     if (trim(grib_msg) .eq. 'sfch') then
        call grib_get(igrib, 'values',dum1d,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: values in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        do r=1,jfguess
           do c=1,ifguess
              fg_hgt_sfc (c,r) = dum1d(c+(r-1)*ifguess)
           enddo ! c
        enddo ! r
        count_hgt_sfc = count_hgt_sfc + 1
     end if
     
     ! Handle isobaric temperature
     if (trim(grib_msg) .eq. 't') then
        plevel = set_plevel(editionNumber_val,pds9_val,level_val)
        do i = 1, kprs
           if (plevel .eq. prslvls(i)) then
              call grib_get(igrib, 'values',dum1d,ierr)
              if ( ierr .ne. 0 ) then
                 write(LIS_logunit,*) &
                      '[WARN] in grib_get: values in ' // &
                      'AGRMET_fldbld_read_gfs'
                 call grib_release(igrib,ierr)
                 goto 100
              end if
              do r=1,jfguess
                 do c=1,ifguess
                    fg_tmp (c,r,i) = dum1d(c+(r-1)*ifguess)
                 enddo ! c
              enddo ! r
              count_tmp      = count_tmp + 1
              exit
           end if
        enddo ! i
     end if
     
     ! Handle geopotential heights
     if (trim(grib_msg) .eq. 'gh') then
        plevel = set_plevel(editionNumber_val,pds9_val,level_val)
        do i = 1, kprs
           if (plevel .eq. prslvls(i)) then
              call grib_get(igrib, 'values',dum1d,ierr)
              if ( ierr .ne. 0 ) then
                 write(LIS_logunit,*) &
                      '[WARN] in grib_get: values in ' // &
                      'AGRMET_fldbld_read_gfs'
                 call grib_release(igrib,ierr)
                 goto 100
              end if
              do r=1,jfguess
                 do c=1,ifguess
                    fg_hgt (c,r,i) = dum1d(c+(r-1)*ifguess)
                 enddo
              enddo
              count_hgt      = count_hgt + 1
              exit
           end if
        enddo
     end if
     
     ! Handle isobaric RH.
     if (trim(grib_msg) .eq. 'r') then
        plevel = set_plevel(editionNumber_val,pds9_val,level_val)
        do i = 1, kprs
           if (plevel .eq. prslvls(i)) then
              call grib_get(igrib, 'values',dum1d,ierr)
              if ( ierr .ne. 0 ) then
                 write(LIS_logunit,*) &
                      '[WARN] in grib_get: values in ' // &
                      'AGRMET_fldbld_read_gfs'
                 call grib_release(igrib,ierr)
                 goto 100
              end if
              do r=1,jfguess
                 do c=1,ifguess
                    fg_rh (c,r,i) = dum1d(c+(r-1)*ifguess)/100.0
                 enddo
              enddo
              count_rh      = count_rh + 1
              exit
           end if
        enddo
     end if
     
     ! Handle dewpoint depression
     if (trim(grib_msg) .eq. 'dd') then
        plevel = set_plevel(editionNumber_val,pds9_val,level_val)
        do i = 1, kprs
           if (plevel .eq. prslvls(i)) then
              call grib_get(igrib, 'values',dum1d,ierr)
              if ( ierr .ne. 0 ) then
                 write(LIS_logunit,*) &
                      '[WARN] in grib_get: values in ' // &
                      'AGRMET_fldbld_read_gfs'
                 call grib_release(igrib,ierr)
                 goto 100
              end if
              do r=1,jfguess
                 do c=1,ifguess
                    fg_dpd (c,r,i) = dum1d(c+(r-1)*ifguess)
                 enddo
              enddo
              count_dpd      = count_dpd + 1
              exit
           end if
        enddo
     end if
  
     ! Handle 10-m u-wind component
     if (trim(grib_msg) .eq. '10u') then
        call grib_get(igrib, 'values',dum1d,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: values in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        do r=1,jfguess
           do c=1,ifguess
              fg_uwnd (c,r) = dum1d(c+(r-1)*ifguess)
           enddo
        enddo
        count_uwnd_sfc = count_uwnd_sfc + 1
     end if

     ! Handle 10-m v-wind component     
     if (trim(grib_msg) .eq. '10v') then
        call grib_get(igrib, 'values',dum1d,ierr)
        if ( ierr .ne. 0 ) then
           write(LIS_logunit,*) &
                '[WARN] in grib_get: values in ' // &
                'AGRMET_fldbld_read_gfs'
           call grib_release(igrib,ierr)
           goto 100
        end if
        do r=1,jfguess
           do c=1,ifguess
              fg_vwnd (c,r) = dum1d(c+(r-1)*ifguess)
           enddo
        enddo
        count_vwnd_sfc = count_vwnd_sfc + 1
     end if
     
     call grib_release(igrib,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_release in ' // &
             'AGRMET_fldbld_read_gfs'
        goto 100
     end if

     ! See if we have everything
!     if ( count_tmp .eq. kprs .and. count_uwnd .eq. 1 .and. &
!          count_hgt .eq. kprs .and. count_vwnd .eq. 1 .and. &
!          (count_rh .eq. kprs .or. count_dpd .eq. kprs) ) then
     if ( count_tmp .eq. kprs .and. count_tmp_sfc .eq. 1 .and. &
          count_hgt .eq. kprs .and. count_hgt_sfc .eq. 1 .and. &
          (count_rh  .eq. kprs .or. count_dpd .eq. kprs) .and. &
          count_rh_sfc .eq.  1 .and. &
          count_uwnd_sfc .eq. 1 .and. count_vwnd_sfc .eq. 1 .and. &
          count_pres_sfc .eq. 1) then
        exit
     end if
  enddo ! k
  
  ! Make sure we have all required fields.  This may require converting
  ! dew point depression to RH
  found = .false.
!  if ( count_tmp .eq. kprs .and. count_uwnd .eq. 1 .and. &
!       count_hgt .eq. kprs .and. count_vwnd .eq. 1 ) then
  if ( count_tmp .eq. kprs .and. count_tmp_sfc .eq. 1 .and. &
       count_hgt .eq. kprs .and. count_hgt_sfc .eq. 1 .and. &
       count_uwnd_sfc .eq. 1 .and. count_vwnd_sfc .eq. 1 .and. &
       count_pres_sfc .eq. 1) then

     if (count_rh .eq. kprs) then
        found = .true.
     else
        if  (count_dpd .eq. kprs)  then           
           do k = 1, kprs
              do j = 1, jfguess
                 do i = 1, ifguess
                    fg_rh(i,j,k) = &
                         AGRMET_calcrh ( fg_tmp(i,j,k), fg_dpd(i,j,k) )
                 enddo ! i
              enddo ! j
           enddo ! k
           found = .true.
        end if
     end if
  end if
  
  ! Handle missing data
  if (.not. found) then
     write(LIS_logunit,*) &
          '[WARN] Missing data from GFS GRIB file!'
     write(LIS_logunit,*) &
          'count_tmp = ',count_tmp,', should be ',kprs
     write(LIS_logunit,*) &
          'count_tmp_sfc = ',count_tmp_sfc,', should be ',1
     write(LIS_logunit,*) &
          'count_hgt = ',count_hgt,', should be ',kprs
     write(LIS_logunit,*) &
          'count_hgt_sfc = ',count_hgt_sfc,', should be ',1
     write(LIS_logunit,*) &
          'count_rh = ',count_rh,', should be ',kprs
     write(LIS_logunit,*) &
          'OR: count_dpd = ',count_dpd,', should be ',kprs
     write(LIS_logunit,*) &
          'count_rh_sfc = ',count_rh_sfc,', should be ',1
     write(LIS_logunit,*) &
          'count_uwnd_sfc = ',count_uwnd_sfc,', should be ',1
     write(LIS_logunit,*) &
          'count_vwnd_sfc = ',count_vwnd_sfc,', should be ',1
     write(LIS_logunit,*) &
          'count_pres_sfc = ',count_pres_sfc,', should be ',1
     goto 100
  end if

  ! Convert from u and v component winds to wind speed. 
  ! Constrain the speeds between 50 m/s and the user
  ! specified minimum.
  fg_wspd = &
       sqrt ((fg_uwnd * fg_uwnd) + (fg_vwnd * fg_vwnd)) * wndwgt
  
  fg_wspd = min( max( fg_wspd, minwnd ), 50.0 ) 
  
  ! Clean up
  deallocate ( dum1d )     
  deallocate ( fg_dpd  )
  deallocate ( fg_uwnd )
  deallocate ( fg_vwnd )
  rc = 0
  return

  ! Jump down here to clean up memory before returning after finding a
  ! problem.
  100 continue
  call grib_close_file(ftn)
  deallocate ( dum1d )
  deallocate ( fg_dpd  )
  deallocate ( fg_uwnd )
  deallocate ( fg_vwnd )
  rc = 1
  return

#else
  
  ! No GRIB_API
  rc = 1
  return

#endif

end subroutine AGRMET_fldbld_read_gfs


!BOP
!
! !ROUTINE: AGRMET_fg2lis
! \label{AGRMET_fg2lis}
!
! !REVISION HISTORY: 
! 28 Jul 08    Sujay Kumar; Initial Specification
! 09 Jun 2017  James Geiger; Refactor interpolation
!
! !INTERFACE: 
subroutine AGRMET_fg2lis(n, ifguess, jfguess, fg_field, agr_field)
  
! !USES: 
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer,     intent(in)  :: n
  integer,     intent(in)  :: ifguess
  integer,     intent(in)  :: jfguess
  real,        intent(in)  :: fg_field  ( ifguess, jfguess )
  real,        intent(out) :: agr_field ( LIS_rc%lnc(n),LIS_rc%lnr(n))
!     
! !DESCRIPTION:    
! 
!     interpolate first guess data to the agrmet grid.
!     
!     \textbf{Method} \newline
!       interpolate geopotential heights, temperature, relative
!       humidity, and wind speed to the agrmet grid using a
!       four point interpolation.
!       
!     The arguments are: 
!     \begin{description}
!      \item[n]
!        nest index
!      \item[ifguess]
!        i-dimension of the first guess grid
!      \item[jfguess]
!        j-dimension of the first guess grid
!      \item[fg\_field]
!       array of first guess field
!      \item[agr\_field]
!        field on agrmet grid
!      \item[i,j,k]
!        loop indices
!      \end{description}
!EOP
  integer                              :: mi, mo
  integer                              :: i,j
  integer                              :: iret
  integer                              :: midway
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
  ! translate from N-S to S-N
  midway = ifguess/2
  do j=1,jfguess
     do i=1,ifguess
        if((i+midway) <= ifguess) then
           var(i,j) = fg_field((i+midway),(jfguess-j)+1)
        else
           var(i,j) = fg_field((i-midway),(jfguess-j)+1)
        endif
     enddo
  enddo

  ! This routine does not process precip.  See AGRMET_fg2lis_precip.
  ! When budget-bilinear is selected, it is applied only to precip,
  ! all other fields use bilinear.  Handle this case with the following
  ! compound conditional expression.
  if ( agrmet_struc(n)%fg_gfs_interp == 'bilinear' .or. &
       agrmet_struc(n)%fg_gfs_interp == 'budget-bilinear' ) then
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                &
             var,lo,agr_field,mi,mo,                              &
             LIS_domain(n)%lat,LIS_domain(n)%lon,                 &
             agrmet_struc(n)%w11_1_gfs,agrmet_struc(n)%w12_1_gfs, &
             agrmet_struc(n)%w21_1_gfs,agrmet_struc(n)%w22_1_gfs, &
             agrmet_struc(n)%n11_1_gfs,agrmet_struc(n)%n12_1_gfs, &
             agrmet_struc(n)%n21_1_gfs,agrmet_struc(n)%n22_1_gfs, &
             LIS_rc%udef, iret)
  elseif ( agrmet_struc(n)%fg_gfs_interp == 'neighbor' ) then
      call neighbor_interp(LIS_rc%gridDesc(n,:),lb,                    &
                           var,lo,agr_field,mi,mo,                     &
                           LIS_domain(n)%lat, LIS_domain(n)%lon,       &
                           agrmet_struc(n)%n11_1_gfs,LIS_rc%udef,iret)
  elseif ( agrmet_struc(n)%fg_gfs_interp == 'average' ) then
     call upscaleByAveraging(mi,mo,LIS_rc%udef, &
                             agrmet_struc(n)%n11_1_gfs,lb,var,lo,agr_field)
  else
     write(LIS_logunit,*) 'ERR: Unexpected interpolation method'
     write(LIS_logunit,*) '     in AGRMET_fg2lis'
     write(LIS_logunit,*) '     ', trim(agrmet_struc(n)%fg_gfs_interp)
     call LIS_endrun
  endif

  deallocate(var)
  deallocate(lb)
  deallocate(lo)

end subroutine AGRMET_fg2lis


!BOP
!
! !ROUTINE: AGRMET_fg2lis_precip
! \label{AGRMET_fg2lis_precip}
!
! !REVISION HISTORY: 
!     10 May 13 Initial specification, forked form AGRMET_fg2lis
!               ..................................Ryan Ruhge/16WS/WXE/SEMS
!
! !INTERFACE: 
subroutine AGRMET_fg2lis_precip(n, findex, ifguess, jfguess, fg_prec,agr_prec)
  
! !USES: 
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer,     intent(in)    :: n
  integer,     intent(in)    :: findex
  integer,     intent(in)    :: ifguess
  integer,     intent(in)    :: jfguess

  real,        intent(in)    :: fg_prec  ( ifguess, jfguess )  

  real,        intent(out)   :: agr_prec ( LIS_rc%lnc(n),LIS_rc%lnr(n))
!     
! !DESCRIPTION:    
! 
!     interpolate first guess data to the agrmet grid.
!     
!     \textbf{Method} \newline
!       interpolate precipitation to the agrmet grid using a
!       four point interpolation.
!       
!     The arguments are: 
!     \begin{description}
!      \item[fg\_prec]
!       array of first guess precipitation
!      \item[ifguess]
!        i-dimension of the first guess grid
!      \item[jfguess]
!        j-dimension of the first guess grid
!      \item[agr\_prec]
!        precipitation on the agrmet grid
!      \item[i,j,k]
!        loop indices
!      \end{description}
!EOP
  integer           :: ibi, ibo
  integer           :: mi, mo
  integer           :: k 
  integer           :: i,j
  integer           :: iret
  real, allocatable, dimension(:,:)    :: var
  logical*1, allocatable, dimension(:) :: lb
  logical*1, allocatable, dimension(:) :: lo

  allocate(var(ifguess,jfguess))
  allocate(lb(ifguess*jfguess))
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  mi = ifguess * jfguess
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  ibi = 1
  lb = .true. 
  
  do i=1,ifguess
     do j=1,jfguess
        if((i+ifguess/2) <= ifguess) then
           var(i,j) = fg_prec((i+ifguess/2),(jfguess-j)+1)
        else
           var(i,j) = fg_prec((i-ifguess/2),(jfguess-j)+1)
        endif
     enddo
  enddo

  if ( agrmet_struc(n)%fg_gfs_interp == 'bilinear' ) then
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                &
             var,lo,agr_prec,mi,mo,                               &
             LIS_domain(n)%lat,LIS_domain(n)%lon,                 &
             agrmet_struc(n)%w11_1_gfs,agrmet_struc(n)%w12_1_gfs, &
             agrmet_struc(n)%w21_1_gfs,agrmet_struc(n)%w22_1_gfs, &
             agrmet_struc(n)%n11_1_gfs,agrmet_struc(n)%n12_1_gfs, &
             agrmet_struc(n)%n21_1_gfs,agrmet_struc(n)%n22_1_gfs, &
             LIS_rc%udef, iret)
  elseif ( agrmet_struc(n)%fg_gfs_interp == 'budget-bilinear' ) then
     call conserv_interp(LIS_rc%gridDesc(n,:),lb,                 &
             var,lo,agr_prec,mi,mo,                               &
             LIS_domain(n)%lat,LIS_domain(n)%lon,                 &
             agrmet_struc(n)%w11_2_gfs,agrmet_struc(n)%w12_2_gfs, &
             agrmet_struc(n)%w21_2_gfs,agrmet_struc(n)%w22_2_gfs, &
             agrmet_struc(n)%n11_2_gfs,agrmet_struc(n)%n12_2_gfs, &
             agrmet_struc(n)%n21_2_gfs,agrmet_struc(n)%n22_2_gfs, &
             LIS_rc%udef, iret)
  elseif ( agrmet_struc(n)%fg_gfs_interp == 'neighbor' ) then
      call neighbor_interp(LIS_rc%gridDesc(n,:),lb,                    &
                           var,lo,agr_prec,mi,mo,                      &
                           LIS_domain(n)%lat, LIS_domain(n)%lon,       &
                           agrmet_struc(n)%n11_1_gfs,LIS_rc%udef,iret)
  elseif ( agrmet_struc(n)%fg_gfs_interp == 'average' ) then
     call upscaleByAveraging(mi,mo,LIS_rc%udef, &
                             agrmet_struc(n)%n11_1_gfs,lb,var,lo,agr_prec)
  else
     write(LIS_logunit,*) 'ERR: Unexpected interpolation method'
     write(LIS_logunit,*) '     in AGRMET_fg2lis_precip'
     write(LIS_logunit,*) '     ', trim(agrmet_struc(n)%fg_gfs_interp)
     call LIS_endrun
  endif

  deallocate(var)
  deallocate(lb)
  deallocate(lo)

end subroutine AGRMET_fg2lis_precip


!BOP
!
! !ROUTINE: gfs_reset_interp_input
!  \label{gfs_reset_interp_input}
!
! !REVISION HISTORY:
!  12 Jun 2017: James Geiger: Initial specification
!
! !INTERFACE:
subroutine gfs_reset_interp_input(n, findex, gridDesci)
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

   deallocate(agrmet_struc(n)%n11_1_gfs, stat=rc)
   deallocate(agrmet_struc(n)%n12_1_gfs, stat=rc)
   deallocate(agrmet_struc(n)%n21_1_gfs, stat=rc)
   deallocate(agrmet_struc(n)%n22_1_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w11_1_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w12_1_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w21_1_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w22_1_gfs, stat=rc)

   deallocate(agrmet_struc(n)%n11_2_gfs, stat=rc)
   deallocate(agrmet_struc(n)%n12_2_gfs, stat=rc)
   deallocate(agrmet_struc(n)%n21_2_gfs, stat=rc)
   deallocate(agrmet_struc(n)%n22_2_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w11_2_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w12_2_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w21_2_gfs, stat=rc)
   deallocate(agrmet_struc(n)%w22_2_gfs, stat=rc)

   howtoTransform = LIS_howtoTransform(n,max(gridDesci(9),gridDesci(10)))

   if ( howtoTransform == 'interpolate') then

      agrmet_struc(n)%fg_gfs_interp = LIS_rc%met_interp(findex)

      write(LIS_logunit,*) '[INFO] The GFS forcing resolution is coarser ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Interpolating with the ' // &
                           trim(agrmet_struc(n)%fg_gfs_interp) // ' method.'

      if ( agrmet_struc(n)%fg_gfs_interp == 'bilinear' ) then

         allocate(agrmet_struc(n)%n11_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n12_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n21_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n22_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w11_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w12_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w21_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w22_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,        &
                 agrmet_struc(n)%n11_1_gfs,agrmet_struc(n)%n12_1_gfs, &
                 agrmet_struc(n)%n21_1_gfs,agrmet_struc(n)%n22_1_gfs, &
                 agrmet_struc(n)%w11_1_gfs,agrmet_struc(n)%w12_1_gfs, &
                 agrmet_struc(n)%w21_1_gfs,agrmet_struc(n)%w22_1_gfs)

      elseif ( agrmet_struc(n)%fg_gfs_interp == 'budget-bilinear' ) then

         allocate(agrmet_struc(n)%n11_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n12_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n21_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%n22_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w11_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w12_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w21_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(agrmet_struc(n)%w22_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,        &
                 agrmet_struc(n)%n11_1_gfs,agrmet_struc(n)%n12_1_gfs, &
                 agrmet_struc(n)%n21_1_gfs,agrmet_struc(n)%n22_1_gfs, &
                 agrmet_struc(n)%w11_1_gfs,agrmet_struc(n)%w12_1_gfs, &
                 agrmet_struc(n)%w21_1_gfs,agrmet_struc(n)%w22_1_gfs)

         allocate(agrmet_struc(n)%n11_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%n12_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%n21_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%n22_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w11_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w12_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w21_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(agrmet_struc(n)%w22_2_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

         call conserv_interp_input(n,gridDesci,         &
                 agrmet_struc(n)%n11_2_gfs,agrmet_struc(n)%n12_2_gfs, &
                 agrmet_struc(n)%n21_2_gfs,agrmet_struc(n)%n22_2_gfs, &
                 agrmet_struc(n)%w11_2_gfs,agrmet_struc(n)%w12_2_gfs, &
                 agrmet_struc(n)%w21_2_gfs,agrmet_struc(n)%w22_2_gfs)
      endif
   elseif ( howtoTransform == 'neighbor') then
      agrmet_struc(n)%fg_gfs_interp = 'neighbor'

      write(LIS_logunit,*) '[INFO] The GFS forcing resolution is comparable ' // &
                           'to the running domain.'
      write(LIS_logunit,*) '     Interpolating with the ' // &
                           trim(agrmet_struc(n)%fg_gfs_interp) // ' method.'

      allocate(agrmet_struc(n)%n11_1_gfs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

      call neighbor_interp_input(n,gridDesci,agrmet_struc(n)%n11_1_gfs)
   elseif ( howtoTransform == 'upscale' ) then
      agrmet_struc(n)%fg_gfs_interp = LIS_rc%met_upscale(findex)

      write(LIS_logunit,*) '[INFO] The GFS forcing resolution is finer ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Upscaling with the ' // &
                           trim(agrmet_struc(n)%fg_gfs_interp) // ' method.'

      select case( agrmet_struc(n)%fg_gfs_interp )
      case( 'average' )
         mi = gridDesci(2) * gridDesci(3)
         allocate(agrmet_struc(n)%n11_1_gfs(mi))

         call upscaleByAveraging_input(gridDesci,                   &
                                       LIS_rc%gridDesc(n,:),        &
                                       mi,                          &
                                       LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                                       agrmet_struc(n)%n11_1_gfs)
      case default
         write(LIS_logunit,*) 'The specified spatial interpolation option '
         write(LIS_logunit,*) 'is not supported for GFS.'
         write(LIS_logunit,*) 'LIS is stopping.'
         call LIS_endrun()
      end select
   else
      write(LIS_logunit,*) 'Unexpected spatial transformation method '
      write(LIS_logunit,*) howtoTransform
      write(LIS_logunit,*) 'LIS is stopping.'
      call LIS_endrun()
   endif
end subroutine gfs_reset_interp_input

! EMK...Determine if GRIB1 field is required
function check_gfs_grib1_message(editionNumber,pds7_val,pds8_val,pds9_val)
   implicit none
   integer, intent(in) :: editionNumber
   integer, intent(in) :: pds7_val
   integer, intent(in) :: pds8_val
   integer, intent(in) :: pds9_val
   character(len=4)    :: check_gfs_grib1_message

   check_gfs_grib1_message = 'none'
   if (editionNumber .ne. 1) return
   
   if (pds7_val.eq.1 .and. pds8_val.eq.1 .and. pds9_val.eq.0) then
      check_gfs_grib1_message = 'sp' ! surface pressure
   else if (pds7_val.eq.11 .and. pds8_val.eq.105 .and. pds9_val.eq.2) then
      check_gfs_grib1_message = '2t' ! 2-m temperature
   else if (pds7_val.eq.52 .and. pds8_val.eq.105 .and. pds9_val.eq.2) then
      check_gfs_grib1_message = '2rh' ! 2-m RH
   else if (pds7_val.eq.7 .and. pds8_val.eq.1) then
      check_gfs_grib1_message = 'sfch' ! surface height
   else if (pds7_val .eq. 11 .and. pds8_val .eq. 100) then
      check_gfs_grib1_message = 't' ! isobaric temperature
   else if (pds7_val .eq. 7 .and. pds8_val .eq. 100) then
      check_gfs_grib1_message = 'gh' ! isobaric geopotential heights
   else if (pds7_val .eq. 52 .and. pds8_val .eq. 100) then
      check_gfs_grib1_message = 'r' ! isobaric RH
   else if (pds7_val .eq. 18 .and. pds8_val .eq. 100) then
      check_gfs_grib1_message = 'dd' ! isobaric dewpoint depression
   else if (pds7_val.eq.33 .and. pds8_val.eq.105 .and. pds9_val.eq.10) then
      check_gfs_grib1_message = '10u' ! 10-m U wind
   else if (pds7_val.eq.34 .and. pds8_val.eq.105 .and. pds9_val.eq.10) then
      check_gfs_grib1_message = '10v' ! 10-m U wind
   end if   
end function check_gfs_grib1_message

! EMK...Determine if GRIB2 field is required
function check_gfs_grib2_message(editionNumber,param_disc_val,param_cat_val,&
     param_num_val,typeOfFirstFixedSurface_val,level_val)
   implicit none
   integer, intent(in) :: editionNumber
   integer, intent(in) :: param_disc_val
   integer, intent(in) :: param_cat_val
   integer, intent(in) :: param_num_val
   integer, intent(in) :: typeOfFirstFixedSurface_val
   integer, intent(in) :: level_val
   character(len=4)    :: check_gfs_grib2_message

   check_gfs_grib2_message = 'none'
   if (editionNumber .ne. 2) return

   if (param_disc_val == 0 .and. param_cat_val == 3 .and. &
        param_num_val == 0 .and. typeOfFirstFixedSurface_val == 1) then
      check_gfs_grib2_message = 'sp' ! surface pressure
   else if (param_disc_val == 0 .and. param_cat_val == 0 .and. &
        param_num_val == 0 .and. &
        typeOfFirstFixedSurface_val == 103 .and. &
        level_val == 2) then
      check_gfs_grib2_message = '2t' ! 2-m temperature
   else if (param_disc_val == 0 .and. param_cat_val == 1 .and. &
        param_num_val == 1 .and. &
        typeOfFirstFixedSurface_val == 103 .and. &
        level_val == 2) then
      check_gfs_grib2_message = '2rh' ! 2-m RH
   else if ( param_disc_val == 0 .and. param_cat_val == 3 .and. &
        param_num_val == 5 .and. &
        typeOfFirstFixedSurface_val == 1) then
      check_gfs_grib2_message = 'sfch' ! surface height
   else if (param_disc_val == 0 .and. param_cat_val == 0 .and. &
        param_num_val == 0 .and.&
        typeOfFirstFixedSurface_val == 100) then
      check_gfs_grib2_message = 't' ! isobaric temperature
   else if (param_disc_val == 0 .and. param_cat_val == 3 .and. &
        param_num_val == 5 .and.&
        typeOfFirstFixedSurface_val == 100) then
      check_gfs_grib2_message = 'gh' ! isobaric geopotential height
   else if (param_disc_val == 0 .and. param_cat_val == 1 .and. &
        param_num_val == 1 .and.&
        typeOfFirstFixedSurface_val == 100) then
      check_gfs_grib2_message = 'r' ! isobaric relative humidity
   else if (param_disc_val == 0 .and. param_cat_val == 2 .and. &
        param_num_val == 2 .and. &
        typeOfFirstFixedSurface_val == 103 .and. &
        level_val == 10) then
      check_gfs_grib2_message = '10u' ! 10-meter U wind
   else if (param_disc_val == 0 .and. param_cat_val == 2 .and. &
        param_num_val == 3 .and. &
        typeOfFirstFixedSurface_val == 103 .and. &
        level_val == 10) then
      check_gfs_grib2_message = '10v' ! 10-meter V wind
   end if
end function check_gfs_grib2_message

! EMK...Return pressure level based on GRIB metadata
integer function set_plevel(editionNumber,pds9,level)

   ! Imports
   use LIS_coreMod, only: LIS_masterproc
   use LIS_logmod, only: LIS_logunit,LIS_abort, &
      LIS_alert,LIS_endrun
   use LIS_mpiMod

   ! Defaults
   implicit none

   ! Arguments
   integer, intent(in) :: editionNumber
   integer, intent(in) :: pds9
   integer, intent(in) :: level

   ! Locals
   integer :: plevel
   integer :: ierr
   character(len=255) :: messages(20)

   if (editionNumber == 1) then
      plevel = pds9
   else if (editionNumber == 2) then
      plevel = level
   else
      write(LIS_logunit,*) &
        '[ERR] Unknown GRIB edition ',editionNumber
      write(LIS_logunit,*) &
        'ABORTING...'
      flush(LIS_logunit)
      messages(:) = ''
      messages(1) = '[ERR] Program: LIS'
      messages(2) = '  Routine: set_plevel'
      messages(3) = '  Unknown GRIB edition'
#if (defined SPMD)
     call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif            
      if (LIS_masterproc) then
         call LIS_alert('LIS.set_plevel',1,messages)
	 call LIS_abort(messages)
      else 
      	 call sleep(10) ! Make sure LIS_masterproc finishes LIS_abort
	 call LIS_endrun()
      end if
   end if
   set_plevel = plevel
end function set_plevel
