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

subroutine USAF_fldbld_radflux_gfs(n, julhr, fg_swdata, &
     fg_lwdata, rc)

  ! Imports
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_coreMod,       only : LIS_rc
  use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_alert, &
                                LIS_verify, LIS_endrun
  use LIS_timeMgrMod,    only : LIS_julhr_date
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: n
  integer, intent(in) :: julhr
  real, intent(out)   :: fg_swdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real, intent(out)   :: fg_lwdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
  integer, intent(out):: rc

  ! Locals
  integer                 :: ftn, ftn2
  integer                 :: igrib, igrib2
  character*250           :: avnfile, avnfile2
  integer                 :: yr1, mo1, da1, hr1
  integer                 :: fc_hr
  character*255           :: message     ( 20 )
  integer                 :: iginfo(2), iginfo2(2)
  real                    :: gridres, gridres2
  integer                 :: alert_number
  real, allocatable       :: fg_swdown1(:,:), fg_swdown2(:,:)
  real, allocatable       :: fg_lwdown1(:,:), fg_lwdown2(:,:)
  integer                 :: ifguess, ifguess2
  integer                 :: jfguess, jfguess2
  integer                 :: center, center2
  integer                 :: ierr, ierr2
  logical*1               :: found, found2
  logical                 :: first_time
  integer                 :: yr_2d
  integer                 :: file_julhr
  integer                 :: dataDate, dataDate2
  integer                 :: dataTime, dataTime2
  character*100           :: gtype, gtype2
  logical                 :: found_inq, found_inq2
  integer :: getsixhr

  ! External subroutines
  external :: AGRMET_fg2lis
  external :: getAVNfilename, getAVNfilename2
  external :: USAF_fldbld_read_radflux_gfs

  ! Initialize return code to "no error".  We will change it below if
  ! necessary.
  rc = 0

  ! Will search previous GFS cycles every six hours, up to 30 hours,
  ! until we find an acceptable file.
  fc_hr = 0           ! Incremented below
  file_julhr = julhr  ! Decremented below
  call LIS_julhr_date(file_julhr, yr1, mo1, da1, hr1)

  ! GFS is only available 3-hrly in NCEI archive, so we will use that here.
  ! Also, files contain time averaged radiation, with the time average
  ! *ending* at the forecast hour.  This creates an eastward bias in the
  ! position of the Sun. We account for that here.
  if (hr1 ==  0 .or. hr1 ==  6 .or. hr1 ==  12 .or. hr1 == 18) then
     fc_hr = 9
     file_julhr = file_julhr - 6
  else if (hr1 ==  1 .or. hr1 ==  7 .or. hr1 == 13 .or. hr1 == 19) then
     fc_hr = 9
     file_julhr = file_julhr - 7
  else if (hr1 == 2 .or. hr1 == 8 .or. hr1 == 14 .or. hr1 == 20) then
     fc_hr = 9
     file_julhr = file_julhr - 8
  else if (hr1 == 3 .or. hr1 == 9 .or. hr1 == 15 .or. hr1 == 21) then
     fc_hr = 12
     file_julhr = file_julhr - 9
  else if (hr1 == 4 .or. hr1 == 10 .or. hr1 == 16 .or. hr1 == 22) then
     fc_hr = 12
     file_julhr = file_julhr - 10
  else if (hr1 == 5 .or. hr1 == 11 .or. hr1 == 17 .or. hr1 == 23) then
     fc_hr = 12
     file_julhr = file_julhr - 11
  end if

  ! Some GFS files have 6-hr time averages of radiation instead of 3-hr.
  ! This requires taking a weighted difference to estimate the most
  ! recent 3-hr time average.  (Equivalent to multiplying average values
  ! by 6 or 3 hours to get "accumulations", differencing the
  ! "accumulations", and then dividing by 3 hours.)
  if (mod(fc_hr,12) .eq. 0) then
     getsixhr = 1
     found2 = .false.
  else
     getsixhr = 0
     found2 = .true.
  end if
  call LIS_julhr_date(file_julhr, yr1, mo1, da1, hr1)

  found = .false.
  first_time = .true.
  do while ( .not. found )

     ! Make sure we start with the previous GFS cycle.
     if ( (.not. first_time) .or. &
          (first_time .and. fc_hr < 6)) then
        fc_hr = fc_hr + 6
        if (fc_hr > 30) exit ! Give up
        file_julhr = file_julhr - 6
        call LIS_julhr_date(file_julhr, yr1, mo1, da1, hr1)
     end if
     first_time = .false.

     yr_2d = mod(yr1, 100)
     if (yr_2d == 0) yr_2d = 100

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
     end if

     ! See if the GRIB file exists before calling ECCODES.
     inquire(file=trim(avnfile), exist=found_inq)
     if (.not. found_inq) then
        write(LIS_logunit,*) '[WARN] Cannot find file ' // trim(avnfile)
        cycle
     end if
     if (getsixhr .eq. 1) then
        inquire(file=trim(avnfile2), exist=found_inq2)
        if (.not. found_inq2) then
           write(LIS_logunit,*) &
                '[WARN] Cannot find file ' // trim(avnfile2)
           cycle
        end if
     end if

#if (defined USE_GRIBAPI)
     call grib_open_file(ftn, trim(avnfile), 'r', ierr)
     if ( ierr /= 0 ) then
        write(LIS_logunit,*) '[WARN] Failed to open first guess - ', &
             trim(avnfile)
        flush(LIS_logunit)
        call grib_close_file(ftn)
        cycle
     end if
     if (getsixhr .eq. 1) then
        call grib_open_file(ftn2, trim(avnfile2), 'r', ierr2)
        if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) '[WARN] Failed to open first guess - ', &
               trim(avnfile2)
          call grib_close_file(ftn2)
          flush(LIS_logunit)
          cycle
       end if
    end if

    ! Extract some information
    call grib_new_from_file(ftn, igrib, ierr)
    ierr2 = 0
    if (getsixhr .eq. 1) call grib_new_from_file(ftn2, igrib2, ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) '[WARN] Failed file read check - ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) '[WARN] Failed file read check - ' // &
               trim(avnfile2)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if (getsixhr .eq. 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       end if
       cycle
    endif

    call grib_get(igrib, 'centre', center, ierr)
    ierr2 = 0
    if (getsixhr .eq. 1) call grib_get(igrib2, 'centre', center2, ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) '[WARN] Cannot read: centre in ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) '[WARN] Cannot read: centre in ' // &
               trim(avnfile2)
          flush(LIS_logunit)
       end if
       if (getsixhr == 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       endif
       cycle
    endif

    call grib_get(igrib, 'gridType', gtype, ierr)
    ierr2 = 0
    if (getsixhr .eq. 1) call grib_get(igrib2, 'gridType', gtype2, ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) '[WARN] Cannot read: gridType in ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) '[WARN] Cannot read: gridType in ' // &
               trim(avnfile2)
          flush(LIS_logunit)
       end if
       if (getsixhr == 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       end if
       cycle
    end if

    call grib_get(igrib, 'Ni', iginfo(1), ierr)
    ierr2 = 0
    if (getsixhr == 1) call grib_get(igrib2, 'Ni', iginfo2(1), ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) '[WARN] Cannot read: Ni in ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if (ierr2 /= 0) then
          write(LIS_logunit,*) '[WARN] Cannot read: Ni in ' // &
               trim(avnfile2)
          flush(LIS_logunit)
       end if
       if (getsixhr == 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       end if
       cycle
    end if

    call grib_get(igrib, 'Nj', iginfo(2), ierr)
    ierr2 = 0
    if (getsixhr == 1) call grib_get(igrib2, 'Nj', iginfo2(2), ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) '[WARN] Cannot read: Nj in ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) '[WARN] Cannot read: Nj in ' // &
               'USAF_fldbld_radflux_gfs'
          flush(LIS_logunit)
       end if
       if (getsixhr == 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       endif
       cycle
    end if


    call grib_get(igrib, 'jDirectionIncrementInDegrees', gridres, &
         ierr)
    ierr2 = 0
    if (getsixhr == 1) call grib_get(igrib2, &
         'jDirectionIncrementInDegrees', gridres2, ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read: jDirectionIncrementInDegrees in ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read: jDirectionIncrementInDegrees in ' // &
               trim(avnfile2)
          flush(LIS_logunit)
       end if
       if (getsixhr == 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       endif
       cycle
    end if

    call grib_get(igrib, 'dataDate', dataDate, ierr)
    ierr2 = 0
    if (getsixhr == 1) call grib_get(igrib2, 'dataDate', dataDate2, ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read: dataDate in ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read: dataDate in ' // &
               trim(avnfile2)
          flush(LIS_logunit)
       end if
       if (getsixhr == 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       endif
       cycle
    end if

    call grib_get(igrib, 'dataTime', dataTime, ierr)
    ierr2 = 0
    if (getsixhr == 1) call grib_get(igrib2, 'dataTime', dataTime2, ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       if ( ierr /= 0 ) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read: dataTime in ' // &
               trim(avnfile)
          flush(LIS_logunit)
       end if
       call grib_release(igrib, ierr)
       call grib_close_file(ftn)
       if ( ierr2 /= 0 ) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read: dataTime in ' // &
               trim(avnfile2)
          flush(LIS_logunit)
       end if
       if (getsixhr == 1) then
          call grib_release(igrib2, ierr2)
          call grib_close_file(ftn2)
       end if
       cycle
    end if

    if ( yr1*10000+mo1*100+da1 == dataDate .and. &
         hr1*100 == dataTime ) then
       found = .TRUE.
       if ( gtype /= "regular_ll" ) then
          message(1) = 'program: LIS'
          message(2) = '  Subroutine: USAF_fldbld_radflux_gfs'
          message(3) = '  First guess source is not a lat/lon grid'
          message(4) = '  USAF_fldbld_radflux_gfs expects lat/lon data'
          call lis_abort(message)
       endif
    endif
    call grib_release(igrib, ierr)
    call grib_close_file(ftn)

    if (getsixhr == 1) then
       if ( yr1*10000+mo1*100+da1 == dataDate2 .and. &
            hr1*100 == dataTime2 ) then
          found2 = .TRUE.
          if ( gtype2 /= "regular_ll" ) then
             message(1) = 'program: LIS'
             message(2) = '  Subroutine: USAF_fldbld_radflux_gfs'
             message(3) = '  First guess source is not a lat/lon grid'
             message(4) = '  USAF_fldbld_radflux_gfs expects lat/lon data'
             call lis_abort(message)
          endif
       endif
       call grib_release(igrib2, ierr2)
       call grib_close_file(ftn2)
    end if

    ! Make sure grids match
    if (getsixhr == 1) then
       if ((iginfo(1) .ne. iginfo2(1)) .or. &
            (iginfo(2) .ne. iginfo2(2)) .or. &
            (gridres .ne. gridres2) .or. &
            (center .ne. center2)) then
          write(LIS_logunit,*) '[ERR]: Grid mismatch between ', &
               trim(avnfile), ' and ', trim(avnfile2)
          message(1) = 'program: LIS'
          message(2) = '  Subroutine: USAF_fldbld_radflux_gfs'
          message(3) = '  First guess source is not a lat/lon grid'
          message(4) = '  USAF_fldbld_radflux_gfs expects lat/lon data'
          call lis_abort(message)
       end if
    end if

#else
     write(LIS_logunit,*) '[ERR]: USAF_fldbld_radflux_gfs requires GRIB-API'
     write(LIS_logunit,*) '[ERR]: please recompile LIS'
     call LIS_endrun
#endif

     ! At this point, we have everything we need.
     found = .true.
     if (found) exit
  end do ! Loop through cycles and forecast hours

  ! Give up if no acceptable GFS file was found
  if (.not. found) then
     write(LIS_logunit,*)'[WARN] No matching GFS Radiation file found!'
     rc = 1
     return
  end if

  write(LIS_logunit,*) &
       '[INFO] Using time-averaged NWP Radiation fields from ',trim(avnfile)
  if (getsixhr == 1) then
     write(LIS_logunit,*) &
          '[INFO] Also using NWP Radiation fields from ',trim(avnfile2)
  end if
  rc = 0

  if (center == 7) then
     write(LIS_logunit,*) '[INFO] Data are from GFS model'
  end if

  write(LIS_logunit,*)'[INFO] GFS RADIATION DATA ARE ON A ', gridres,&
          ' DEGREE LAT/LON GRID'
  ifguess = iginfo(1)
  jfguess = iginfo(2)
  allocate ( fg_swdown1 (ifguess, jfguess) )
  allocate ( fg_lwdown1 (ifguess, jfguess) )

  if (getsixhr == 1) then
     ifguess2 = iginfo2(1)
     jfguess2 = iginfo2(2)
     allocate ( fg_swdown2 (ifguess2, jfguess2) )
     allocate ( fg_lwdown2 (ifguess2, jfguess2) )
  end if

  ! Get radiation for this julian hour
  alert_number = 0
  call USAF_fldbld_read_radflux_gfs(avnfile, ifguess, jfguess, &
       fg_swdown1, fg_lwdown1, alert_number)
  if (getsixhr == 1) then
     call USAF_fldbld_read_radflux_gfs(avnfile2, ifguess2, jfguess2, &
          fg_swdown2, fg_lwdown2, alert_number)
  end if

  ! Estimate time-averaged radiation from last 3 hours.
  ! Some GFS files have 6-hr time averages of radiation instead of 3-hr.
  ! This requires taking a weighted difference to estimate the most
  ! recent 3-hr time average.  (Equivalent to multiplying average values
  ! by 6 or 3 hours to get "accumulations", differencing the
  ! "accumulations", and then dividing by 3 hours.)
  if (getsixhr == 1) then
     write(LIS_logunit,*) &
          '[INFO] Estimating 3-hr time averaged radiation...'
     fg_swdown1 = 2*fg_swdown1 - fg_swdown2
     fg_lwdown1 = 2*fg_lwdown1 - fg_lwdown2
  end if

  ! Interpolate to the LIS grid
  call AGRMET_fg2lis(n, ifguess, jfguess, fg_swdown1, fg_swdata)
  call AGRMET_fg2lis(n, ifguess, jfguess, fg_lwdown1, fg_lwdata)

  ! Clean up
  deallocate(fg_swdown1)
  deallocate(fg_lwdown1)
  if (getsixhr == 1) then
     deallocate(fg_swdown2)
     deallocate(fg_lwdown2)
  end if

end subroutine USAF_fldbld_radflux_gfs

subroutine USAF_fldbld_read_radflux_gfs(fg_filename, ifguess, jfguess, &
     fg_swdown, fg_lwdown, alert_number )

  ! Imports
  use LIS_logMod,  only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  ! Defaults
  implicit none

  ! Arguments
  character(len=*),  intent(in) :: fg_filename
  integer,        intent(in)    :: ifguess
  integer,        intent(in)    :: jfguess
  real,           intent(out)   :: fg_swdown (ifguess,jfguess)
  real,           intent(out)   :: fg_lwdown (ifguess,jfguess)
  integer,        intent(inout) :: alert_number

  ! Locals
  character*255                 :: message  ( 20 )
  integer                       :: count_swdown
  integer                       :: count_lwdown
  integer                       :: ierr
  integer                       :: k
  integer                       :: ftn, igrib, nvars
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val
  real,           allocatable   :: dum1d    ( : )
  logical                       :: found_inq

  ! Make sure file exists before opening with GRIB_API
  inquire(file=trim(fg_filename), exist=found_inq)
  if (.not. found_inq) then
     ierr = 1
  else
     ierr = 0
  end if
  call LIS_verify(ierr, '[ERR] FILE NOT FOUND - ' // trim(fg_filename))

#if (defined USE_GRIBAPI)
  call grib_open_file(ftn, trim(fg_filename), 'r', ierr)
  call LIS_verify(ierr,'[ERR] Failed to open in read routine - ' // &
       trim(fg_filename))

  if ( ierr == 0 ) then
     allocate ( dum1d   (ifguess*jfguess) )
     count_swdown = 0
     count_lwdown = 0

     write(LIS_logunit,*)' '
     write(LIS_logunit,*) &
          '[INFO] Reading time-averaged GFS radiation fluxes'
     write(LIS_logunit,*) trim(fg_filename)

     call grib_count_in_file(ftn, nvars, ierr)
     call LIS_verify(ierr, 'error in grib_count_in_file in ' // &
          'USAF_fldbld_read_radflux_gfs')

     do k = 1, nvars
        call grib_new_from_file(ftn, igrib, ierr)
        call LIS_verify(ierr, &
             '[ERR] failed to read - '// trim(fg_filename))

        call grib_get(igrib, 'discipline', param_disc_val, ierr)
        call LIS_verify(ierr, &
             'error in grib_get: parameterNumber in ' // &
             'USAF_fldbld_read_radflux_gfs')

        call grib_get(igrib, 'parameterCategory', param_cat_val, ierr)
        call LIS_verify(ierr, &
             'error in grib_get: parameterCategory in ' // &
             'USAF_fldbld_read_radflux_gfs')

        call grib_get(igrib, 'parameterNumber', param_num_val, ierr)
        call LIS_verify(ierr, &
             'error in grib_get: parameterNumber in ' // &
             'USAF_fldbld_read_radflux_gfs')

        ! SW Down Radiation flux:
        if ( param_disc_val == 0 .and. param_cat_val == 4 .and. &
             param_num_val  == 192 ) then

           call grib_get(igrib, 'values', dum1d, ierr)
           call LIS_verify(ierr, &
                'error in grib_get: SWdown values in ' // &
                'USAF_fldbld_read_radflux_gfs')
           fg_swdown = reshape(dum1d, (/ifguess,jfguess/))
           count_swdown = count_swdown + 1
        end if

        !  LW Down Radiation flux:
        if ( param_disc_val == 0 .and. param_cat_val    == 5 .and. &
             param_num_val  == 192 ) then
           call grib_get(igrib,'values',dum1d,ierr)
           call LIS_verify(ierr, &
                'error in grib_get: LWdown values in ' // &
                'USAF_fldbld_read_radflux_gfs')
           fg_lwdown = reshape(dum1d, (/ifguess,jfguess/))
           count_lwdown = count_lwdown + 1
        end if

        call grib_release(igrib, ierr)
     end do ! k

     call grib_close_file(ftn)

     deallocate( dum1d)

     ! See if we have everything.  If not, abort.
     if ( count_swdown == 0 .or. count_lwdown == 0 ) then
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  USAF_fldbld_read_radflux_gfs.'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(fg_filename)
        message(5) = '  File does not contain all data that fldbld needs.'
        if ( allocated(dum1d) ) deallocate(dum1d)
        call LIS_abort( message)
     endif
  end if
#endif

end subroutine USAF_fldbld_read_radflux_gfs

