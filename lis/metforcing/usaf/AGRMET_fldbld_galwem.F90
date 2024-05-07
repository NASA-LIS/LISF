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
! 14 Jun 2017  Add GALWEM ground height and 2-m T, RH.....Eric Kemp/GSFC
! 03 Oct 2017  Adjusted logic to first search for 6 to 12 hour GALWEM forecast,
!              instead of 0-6 hr..........................Eric Kemp/GSFC
! 06 Oct 2017  Changed logic to roll-back to earlier GALWEM cycle if file
!              is missing or incomplete...................Eric Kemp/GSFC
! 21 Feb 2020  Added support for 10-km GALWEM.............Eric Kemp/GSFC
! !INTERFACE:    
subroutine AGRMET_fldbld_galwem(n,order,julhr,rc)
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
!  This routine interpolates the UK Unified Model (GALWEM) first guess height,
!  temperature, moisture, surface pressure, and wind data to the AGRMET grid.
!
!EOP
  integer                 :: ftn, igrib
  character*120           :: gribfile
  integer                 :: yr1, mo1, da1, hr1
  character*255           :: message     ( 20 )
  integer                 :: iginfo      ( 40 )
  real                    :: gridres_dlat, gridres_dlon
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
  integer                 :: rc2
  logical :: first_time, found_inq

  ! FIXME...This should be moved into a module for all of LIS to 
  ! reference.
  data prslvls / 1000,975,950,925,900,850,800,750,700,650,600,550,500,&
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
  kprs = 13

  ! Initialize return code to "no error".  We will change it below if 
  ! necessary.
  rc = 0
  
  ! Will search previous GALWEM cycles every six hours, up to 24 hours,
  ! until we find an acceptable file.
  fc_hr = 0           ! Incremented below
  file_julhr = julhr  ! Decremented below
  call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)
  found = .FALSE. 
  first_time = .true.
  do while ( .not. found )

     ! Make sure we start with the previous GALWEM cycle.
     if ( (.not. first_time) .or. &
          (first_time .and. fc_hr < 6)) then
        fc_hr = fc_hr + 6
        if (fc_hr > 24) exit ! Give up
        file_julhr = file_julhr - 6
        call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)
     end if
     first_time = .false.

     yr_2d = mod(yr1,100)
     if(yr_2d.eq.0) yr_2d = 100 
     call AGRMET_getGALWEMfilename(gribfile, agrmet_struc(n)%agrmetdir,&
          agrmet_struc(n)%galwemdir, agrmet_struc(n)%use_timestamp,&
          agrmet_struc(n)%galwem_res, yr1,mo1,da1,hr1,fc_hr)

     ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
     ! using a simple inquire statement.  This avoids ECCODES/GRIB_API 
     ! writing error messages to stdout/stderr, which may lead to runtime
     ! problems.
     inquire(file=trim(gribfile),exist=found_inq)
     if (.not. found_inq) then
        write(LIS_logunit,*) '[WARN] Cannot find file '//trim(gribfile)
        cycle
     end if

     ! Open first guess grib data using library utility.  Just read
     ! the first file only, as all data will be of the same type.
     ! (GALWEM) because the search script ensures that it is.
     ! 
#if (defined USE_GRIBAPI) 
     
     call grib_open_file(ftn,trim(gribfile),'r',ierr)     
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] Failed to open - '//trim(gribfile)
        cycle
     end if

     ! Read in the first grib record, unpack the header and extract
     ! section 1 and section 2 information.
     call grib_new_from_file(ftn,igrib,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] failed to read - '//trim(gribfile)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'centre',center,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grib_get: ' // &
             'centre in AGRMET_fldbld_galwem'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif
        
     call grib_get(igrib,'gridType',gtype,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get: ' // &
             'gridtype in AGRMET_fldbld_galwem'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     if(trim(gtype).ne."regular_ll") then  
        write(LIS_logunit,*)'[WARN] GALWEM file not on lat/lon grid!'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'Ni',iginfo(1),ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get: Ni in AGRMET_fldbld_galwem'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'Nj',iginfo(2),ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get: Nj in AGRMET_fldbld_galwem'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'jDirectionIncrementInDegrees',gridres_dlat,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get: ' // &
             'jDirectionIncrementInDegrees ' // &
             'in AGRMET_fldbld_galwem'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     ! EMK...Added dlon
     call grib_get(igrib, 'iDirectionIncrementInDegrees', gridres_dlon, ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get: ' // &
             'iDirectionIncrementInDegrees ' // &
             'in AGRMET_fldbld_galwem'
        call grib_release(igrib, ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'dataDate',dataDate,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get: ' // &
             'dataDate in AGRMET_fldbld_galwem'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     call grib_get(igrib,'dataTime',dataTime,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grid_get: ' // &
             'dataTime in AGRMET_fldbld_galwem'
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     endif

     if ( yr1*10000+mo1*100+da1 .ne. dataDate) then
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     end if
     if (  hr1*100 .ne. dataTime ) then
        call grib_release(igrib,ierr)
        call grib_close_file(ftn)
        cycle
     end if

     ! Here we tentatively have a file we can use. Close it for now, and
     ! prepare to pull the appropriate variables.
     call grib_release(igrib,ierr)
     call grib_close_file(ftn)

     ! EMK...Added delta lon
     write(LIS_logunit,*)'[INFO] GALWEM FIRST GUESS DELTA LAT IS ', &
          gridres_dlat,' DEGREES'
     write(LIS_logunit,*)'[INFO] GALWEM FIRST GUESS DELTA LON IS ', &
          gridres_dlon,' DEGREES'

     ifguess = iginfo(1)
     jfguess = iginfo(2)
     
     ! FIXME...Reject file if wrong center?
     if (center .eq. 57) then
        write(LIS_logunit,*)'[INFO] FIRST GUESS DATA IS FROM UK UM (GALWEM) MODEL'
     else
        write(LIS_logunit,*)'[WARN] UNKNOWN SOURCE FOR FIRST GUESS DATA'
     end if

     ! Read in first guess data for this julian hour.
     ! EMK 14 June 2017...Include GALWEM terrain height and 2-m T, RH.     
     ! EMK 06 Oct 2017...Added rc2 (return code) argument
     if (order.eq.1) then 
        call AGRMET_fldbld_read_galwem(n, gribfile, ifguess, jfguess,  &
             kprs, prslvls,                                        &
             agrmet_struc(n)%wndwgt, agrmet_struc(n)%minwnd,       &
             agrmet_struc(n)%agr_tmp_c, agrmet_struc(n)%agr_hgt_c, &
             agrmet_struc(n)%agr_rh_c,                             &
             agrmet_struc(n)%agr_tmp_sfc_c,                        &
             agrmet_struc(n)%agr_hgt_sfc_c,                        &
             agrmet_struc(n)%agr_rh_sfc_c,                         &
             agrmet_struc(n)%agr_wspd_c,                           &
             agrmet_struc(n)%agr_pres_c, &
             rc2)
        if (rc2 .eq. 0) then
           agrmet_struc(n)%agr_bgrd_src_c = "GALWEM"
        end if
     else
        call AGRMET_fldbld_read_galwem(n, gribfile, ifguess, jfguess,  &
             kprs, prslvls,                                        &
             agrmet_struc(n)%wndwgt, agrmet_struc(n)%minwnd,       &
             agrmet_struc(n)%agr_tmp_p, agrmet_struc(n)%agr_hgt_p, &
             agrmet_struc(n)%agr_rh_p,                             &
             agrmet_struc(n)%agr_tmp_sfc_p,                        &
             agrmet_struc(n)%agr_hgt_sfc_p,                        &
             agrmet_struc(n)%agr_rh_sfc_p,                         &
             agrmet_struc(n)%agr_wspd_p,                           &
             agrmet_struc(n)%agr_pres_p, &
             rc2)
        if (rc2 .eq. 0) then
           agrmet_struc(n)%agr_bgrd_src_p = "GALWEM"
        end if
     endif

     ! Check return code from subroutine.  If a field was missing, roll
     ! back.
     if (rc2 .ne. 0) cycle

     ! At this point, we have everything we need.
     found = .true. 
     if (found) exit
#endif     

  enddo ! Loop through cycles and forecast hours
  
  ! Give up if no acceptable GALWEM file was found
  if (.not. found) then
     write(LIS_logunit,*)'[WARN] No matching GALWEM file found!'
     rc = 1
     return
  end if

  write(LIS_logunit,*) &
       '[INFO] Using NWP fields from ',trim(gribfile)
  rc = 0

end subroutine AGRMET_fldbld_galwem


!BOP
! 
! !ROUTINE: AGRMET_getGALWEMfilename
! \label{AGRMET_getGALWEMfilename}
!
! !INTERFACE: 
!EMK...Added support for 10-km GALWEM
subroutine AGRMET_getGALWEMfilename(filename,rootdir,dir,use_timestamp, &
     nominal_res_km,yr,mo,da,hr,fc_hr)

  use LIS_logMod, only: LIS_logunit, LIS_endrun
  implicit none
! !ARGUMENTS: 
  character(*)        :: filename
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: nominal_res_km
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

  ! EMK...Support GALWEM 17 or 10-km data
  if (nominal_res_km == 17) then
     fname1 = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GD_GR.C17KM_AR.GLOBAL_DD.'
  else if (nominal_res_km == 10) then
     fname1 = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GD_GR.C10KM_AR.GLOBAL_DD.'
  else
     write(LIS_logunit)'[ERR] Invalid nominal resolution for GALWEM!'
     write(LIS_logunit)'[ERR] Found ', nominal_res_km
     write(LIS_logunit)'[ERR] Only supports 17 and 10'
     call LIS_endrun()
  end if

  ! EMK...Make sure ftime1 is always initialized
  write (UNIT=ftime1, FMT='(i4, i2.2, i2.2)') yr, mo, da

  if (use_timestamp .eq. 1) then 
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
! 14 Jun 2017  Added GALWEM terrain height, 2-m T, RH on full LIS
!              grid.......................................Eric Kemp/GSFC    
! 29 Sep 2017  Added GRIB2 Product Definition Template Number and type
!              of second fixed layer to ensure we are working with 
!              instantaneous variables on horizontal levels or layers. Also
!              added better checks for requested vs found fields.
!              ...........................................Eric Kemp/GSFC
! 11 Oct 2017  Refactored to return if a read error occurs or if not all
!              fields are available.  Return code now included in argument
!              list.......................................Eric Kemp/GSFC
! !INTERFACE: 
subroutine AGRMET_fldbld_read_galwem(n, fg_filename, ifguess, jfguess,     &
                                     kprs, prslvls,                        &
                                     wndwgt, minwnd,                       &
                                     agr_tmp, agr_hgt, agr_rh,             &
                                     agr_tmp_sfc, agr_hgt_sfc, agr_rh_sfc, &
                                     agr_wspd_sfc,                         &
                                     agr_pres_sfc,rc)
! !USES:
  use LIS_coreMod, only : LIS_rc
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
  real,           intent(in)    :: minwnd
  real,           intent(in)    :: wndwgt
  real,           intent(out)   :: agr_hgt (LIS_rc%lnc(n),LIS_rc%lnr(n),kprs)
  real,           intent(out)   :: agr_rh  (LIS_rc%lnc(n),LIS_rc%lnr(n),kprs)
  real,           intent(out)   :: agr_tmp (LIS_rc%lnc(n),LIS_rc%lnr(n),kprs)
  real,           intent(out)   :: agr_hgt_sfc (LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,           intent(out)   :: agr_rh_sfc  (LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,           intent(out)   :: agr_tmp_sfc (LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,           intent(out)   :: agr_wspd_sfc(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,           intent(out)   :: agr_pres_sfc(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer, intent(out) :: rc
!
! !DESCRIPTION:  
!
!     To read UK Unified Model (GALWEM) data in GRIB-2 format.
!     
!EOP
  character*9                   :: cstat
  character*255                 :: message     ( 20 )
  character(len=4)              :: grib_msg
  character(len=4)              :: AGRMET_check_galwem_message
  integer                       :: count_hgt
  integer                       :: count_rh
  integer                       :: count_tmp
  integer                       :: count_hgt_sfc
  integer                       :: count_rh_sfc
  integer                       :: count_tmp_sfc
  integer                       :: count_uwnd_sfc
  integer                       :: count_vwnd_sfc
  integer                       :: count_pres_sfc
  integer                       :: i
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: j,c,r
  integer                       :: igrib
  integer                       :: ftn
  integer                       :: kk, nvars
  integer                       :: prod_def_tmpl_num
  integer                       :: surface_val_2
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, surface_val, level_val
  
  real, allocatable :: dum1d   ( : )
  real, allocatable :: fg_hgt  ( : , : , : )
  real, allocatable :: fg_rh   ( : , : , : )
  real, allocatable :: fg_tmp  ( : , : , : )

  real, allocatable :: fg_hgt_sfc  ( : , : )
  real, allocatable :: fg_rh_sfc   ( : , : )
  real, allocatable :: fg_tmp_sfc  ( : , : )

  real, allocatable :: fg_wspd_sfc ( : , : )
  real, allocatable :: fg_pres_sfc ( : , : )
  real, allocatable :: fg_uwnd_sfc ( : , : )
  real, allocatable :: fg_vwnd_sfc ( : , : )
  logical :: found_inq

  ! Executable code begins here ... 

  rc = 0 ! Initialize as "no error"

  ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
  ! using a simple inquire statement.  This avoids ECCODES/GRIB_API 
  ! writing error messages to stdout/stderr, which may lead to runtime
  ! problems.
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

  allocate ( fg_hgt  (ifguess, jfguess, kprs) )
  allocate ( fg_rh   (ifguess, jfguess, kprs) )
  allocate ( fg_tmp  (ifguess, jfguess, kprs) )
  allocate ( fg_hgt_sfc  (ifguess, jfguess) )
  allocate ( fg_rh_sfc   (ifguess, jfguess) )
  allocate ( fg_tmp_sfc  (ifguess, jfguess) )
  
  allocate ( fg_wspd_sfc (ifguess, jfguess) )
  allocate ( fg_pres_sfc (ifguess, jfguess) )
  
  allocate ( dum1d   (ifguess*jfguess) )
  allocate ( fg_uwnd_sfc (ifguess, jfguess) )
  allocate ( fg_vwnd_sfc (ifguess, jfguess) )

  ! From this point, we must deallocate memory before returning. 
  ! Unfortunately this means using a GOTO statement if a problem is
  ! encountered, but such is life.
  count_hgt  = 0
  count_rh   = 0
  count_tmp  = 0
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
          'AGRMET_fldbld_read_galwem'
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
     
     call grib_get(igrib,'discipline',param_disc_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grib_get: parameterNumber in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if
     
     ! EMK...Now read GRIB2 Product Definition Template Number to 
     ! ensure we have an instantaneous variable at a horizontal level.
     call grib_get(igrib,'productDefinitionTemplateNumber', &
          prod_def_tmpl_num,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: productDefinitionTemplateNumber in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if
     
     call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterCategory in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if
     
     call grib_get(igrib,'parameterNumber',param_num_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterNumber in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if
     
     call grib_get(igrib,'typeOfFirstFixedSurface',surface_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if
     
     call grib_get(igrib,'level',level_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if
     
     ! EMK...Now read GRIB2 type of second level to ensure we are not
     ! working with a layer-averaged field.
     call grib_get(igrib,'typeOfSecondFixedSurface',surface_val_2,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if

     ! We have enough information to determine what GRIB parameter this
     ! is.
     grib_msg = AGRMET_check_galwem_message(param_disc_val, &
          prod_def_tmpl_num, &
          param_cat_val, &
          param_num_val, surface_val, level_val, surface_val_2)

     ! Skip this field if GRIB parameter is not required.
     if (grib_msg == 'none') then
        call grib_release(igrib,ierr)
        if (ierr .ne. 0) then
           write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
                'AGRMET_fldbld_read_galwem'
           goto 100
        end if
        cycle ! Not a message we are interested in.
     end if

     call grib_get(igrib,'values',dum1d,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: values in ' // &
             'AGRMET_fldbld_read_galwem'
        call grib_release(igrib,ierr)
        goto 100
     end if
        
     select case (grib_msg)
     case('sp') ! surface pressure
        fg_pres_sfc = reshape(dum1d, (/ifguess,jfguess/))
        count_pres_sfc = count_pres_sfc + 1
     case ('2t') ! 2-m temperature
        fg_tmp_sfc = reshape(dum1d, (/ifguess,jfguess/))
        count_tmp_sfc = count_tmp_sfc + 1
     case ('2rh') ! 2-m relative humidity
        fg_rh_sfc = reshape(dum1d, (/ifguess,jfguess/))
        fg_rh_sfc = fg_rh_sfc*0.01 
        count_rh_sfc = count_rh_sfc + 1
     case ('sfch') ! ground height
        fg_hgt_sfc = reshape(dum1d, (/ifguess,jfguess/))
        count_hgt_sfc = count_hgt_sfc + 1
     case('t') ! temperature
        do i = 1, kprs
           if ( level_val == prslvls(i) ) then
              fg_tmp(:,:,i) = reshape(dum1d, (/ifguess,jfguess/))
              count_tmp = count_tmp + 1
              exit
           endif
        enddo
     case('gh') ! geopotential height
        do i = 1, kprs
           if ( level_val == prslvls(i) ) then
              fg_hgt(:,:,i) = reshape(dum1d, (/ifguess,jfguess/))
              count_hgt = count_hgt + 1
              exit
           end if
        enddo
     case('r') ! relative humidity
        ! according to WMO standard,
        ! relative humidity is gribbed in percent.
        ! agrmet expects decimal, so divide by 100.
        do i = 1, kprs
           if ( level_val == prslvls(i) ) then
              dum1d = dum1d * 0.01
              fg_rh(:,:,i) = reshape(dum1d, (/ifguess,jfguess/))
              count_rh = count_rh + 1
              exit
           end if
        enddo
     case('10u') ! 10m u-wind
        fg_uwnd_sfc = reshape(dum1d, (/ifguess,jfguess/))
        count_uwnd_sfc = count_uwnd_sfc + 1
     case('10v') ! 10m v-wind
        fg_vwnd_sfc = reshape(dum1d, (/ifguess,jfguess/))
        count_vwnd_sfc = count_vwnd_sfc + 1
     case default ! Internal error, we shouldn't be here
        write(LIS_logunit,*)'[WARN] Unknown grib_message ',grib_msg
        write(LIS_logunit,*)'Aborting...'
        flush(LIS_logunit)
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
     
     ! Finished with this field
     call grib_release(igrib,ierr)
     if (ierr .ne. 0) then
        write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
             'AGRMET_fldbld_read_galwem'
        goto 100
     end if
     
     ! Jump out of loop early if we have everything
     if ( (count_tmp .eq. kprs)   .and. (count_tmp_sfc .eq. 1) .and. &
          (count_hgt .eq. kprs)   .and. (count_hgt_sfc .eq. 1) .and. &
          (count_rh  .eq. kprs)   .and. (count_rh_sfc  .eq. 1) .and. &
          (count_uwnd_sfc .eq. 1) .and. (count_vwnd_sfc .eq. 1) .and. &
          (count_pres_sfc .eq. 1) ) then 
        exit
     end if
     
  enddo ! Loop through all GRIB file fields
  
  ! Make sure we have all required fields.
  if ( (count_tmp .ne. kprs) .or. (count_tmp_sfc .ne. 1) .or. &
       (count_hgt .ne. kprs) .or. (count_hgt_sfc .ne. 1) .or. &
       (count_rh  .ne. kprs) .or. (count_rh_sfc  .ne. 1) .or. &
       (count_uwnd_sfc .ne. 1) .or. (count_vwnd_sfc .ne. 1) .or. &
       (count_pres_sfc .ne. 1) ) then
     write(LIS_logunit,*)'[WARN] Missing data from GALWEM GRIB file!'
     write(LIS_logunit,*)'count_tmp = ',count_tmp,', should be ',kprs
     write(LIS_logunit,*)'count_tmp_sfc = ',count_tmp_sfc,', should be 1'
     write(LIS_logunit,*)'count_hgt = ',count_hgt,', should be ',kprs
     write(LIS_logunit,*)'count_hgt_sfc = ',count_hgt_sfc,', should be 1'
     write(LIS_logunit,*)'count_rh = ',count_rh,', should be ',kprs
     write(LIS_logunit,*)'count_rh_sfc = ',count_rh_sfc,', should be 1'
     write(LIS_logunit,*)'count_uwnd_sfc = ',count_uwnd_sfc,', should be 1'
     write(LIS_logunit,*)'count_vwnd_sfc = ',count_vwnd_sfc,', should be 1'
     write(LIS_logunit,*)'count_pres_sfc = ',count_pres_sfc,', should be 1'
     goto 100
  end if

  ! Convert from u and v component winds to wind speed.
  ! Constrain the speeds between 50 m/s and the user-specified minimum.
  fg_wspd_sfc = &
       sqrt((fg_uwnd_sfc * fg_uwnd_sfc) + (fg_vwnd_sfc * fg_vwnd_sfc)) * wndwgt
     
  fg_wspd_sfc = min( max( fg_wspd_sfc, minwnd ), 50.0 ) 
     
  ! Interpolate the fields to the LIS grid
  call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
       fg_pres_sfc, agr_pres_sfc)
  do i = 1, kprs
     call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
          fg_tmp(:,:,i), agr_tmp(:,:,i))
     call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
          fg_hgt(:,:,i), agr_hgt(:,:,i))
     call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
          fg_rh(:,:,i), agr_rh(:,:,i))
  enddo
  call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
       fg_wspd_sfc, agr_wspd_sfc)
  call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
       fg_tmp_sfc, agr_tmp_sfc)
  call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
       fg_hgt_sfc, agr_hgt_sfc)
  call interp_galwem_first_guess(n, ifguess, jfguess, .false., &
       fg_rh_sfc, agr_rh_sfc)
  
  ! At this point, we have everything.  Close the file and return.
  call grib_close_file(ftn)
  rc = 0
  return

  ! Jump down here to clean up memory before returning after finding a
  ! problem.
  100 continue
  call grib_close_file(ftn)
  deallocate ( dum1d )
  deallocate ( fg_uwnd_sfc )
  deallocate ( fg_vwnd_sfc )
  deallocate ( fg_hgt  )
  deallocate ( fg_rh   )
  deallocate ( fg_tmp  )
  deallocate ( fg_hgt_sfc  )
  deallocate ( fg_rh_sfc   )
  deallocate ( fg_tmp_sfc  )
  deallocate ( fg_wspd_sfc )
  deallocate ( fg_pres_sfc )
  rc = 1
#endif

end subroutine AGRMET_fldbld_read_galwem

!BOP
!
! !ROUTINE: AGRMET_check_galwem_message
! \label{AGRMET_check_galwem_message}
!
! !REVISION HISTORY:
! 14 Jun 2016 James Geiger; Initial specification
! 14 Jun 2017 Eric Kemp; Added support for 2-m T, 2-m RH, and ground height;
!             Also improved check for 10-m U and V winds.
! 29 Sep 2017 Eric Kemp; Added GRIB2 product definition template number
!             (prod_def_tmpl_num) and type of second level (surface_val_2)
!             to ensure variable is instantaneous at horizontal level
!             (not layer).
! !INTERFACE:    
function AGRMET_check_galwem_message(param_disc_val, prod_def_tmpl_num, &
     param_cat_val, &
     param_num_val, surface_val, level_val, surface_val_2)
! !USES: 
! none

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: param_disc_val, prod_def_tmpl_num, &
        param_cat_val, &
        param_num_val, surface_val, level_val, surface_val_2
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

   ! EMK...Only use instantaneous variables
   if (prod_def_tmpl_num .ne. 0) then
      AGRMET_check_galwem_message = 'none'
      return
   end if
   ! EMK...Only use single level fields, not layers
   if (surface_val_2 .ne. 255) then
      AGRMET_check_galwem_message = 'none'
      return
   end if
   if     ( param_disc_val == 0 .and. &
            param_cat_val  == 3 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 1 ) then
      AGRMET_check_galwem_message = 'sp' ! Surface pressure
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 0 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 100 ) then
      AGRMET_check_galwem_message = 't' ! Isobaric temperature
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 0 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 103 .and. &
            level_val == 2) then
      AGRMET_check_galwem_message = '2t' ! 2-m temperature
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 1 .and. &
            param_num_val  == 1 .and. &
            surface_val    == 103 .and. &
            level_val == 2) then
      AGRMET_check_galwem_message = '2rh' ! 2-m relative humidity
!EMK...Use values provided by Jerry Wegiel 12 Sep 2017
!   elseif ( param_disc_val == 0 .and. &
!            param_cat_val  == 3 .and. &
!            param_num_val  == 5 .and. &
!            surface_val    == 1) then
   elseif ( param_disc_val == 2 .and. &
            param_cat_val  == 0 .and. &
            param_num_val  == 7 .and. &
            surface_val    == 1) then
      AGRMET_check_galwem_message = 'sfch' ! Surface terrain height
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 3 .and. &
            param_num_val  == 5 .and. &
            surface_val    == 100 ) then
      AGRMET_check_galwem_message = 'gh' ! Isobaric geopotential height
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 1 .and. &
            param_num_val  == 1 .and. &
            surface_val    == 100 ) then
      AGRMET_check_galwem_message = 'r' ! Isobaric relative humidity
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 2 .and. &
            surface_val    == 103 .and. &
            level_val == 10) then
      AGRMET_check_galwem_message = '10u' ! 10-meter U wind
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 3 .and. &
            surface_val    == 103 .and. &
            level_val == 10) then
      AGRMET_check_galwem_message = '10v' ! 10-meter V wind
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

   if ( howtoTransform == 'interpolate') then

      agrmet_struc(n)%fg_galwem_interp = LIS_rc%met_interp(findex)

      write(LIS_logunit,*) '[INFO] The GALWEM forcing resolution is coarser ' // &
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

      write(LIS_logunit,*) '[INFO] The GALWEM forcing resolution is comparable ' // &
                           'to the running domain.'
      write(LIS_logunit,*) '     Interpolating with the ' // &
                           trim(agrmet_struc(n)%fg_galwem_interp) // ' method.'

      allocate(agrmet_struc(n)%n11_1_galwem(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

      call neighbor_interp_input(n,gridDesci,agrmet_struc(n)%n11_1_galwem)
   elseif ( howtoTransform == 'upscale' ) then
      agrmet_struc(n)%fg_galwem_interp = LIS_rc%met_upscale(findex)

      write(LIS_logunit,*) '[INFO] The GALWEM forcing resolution is finer ' // &
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

