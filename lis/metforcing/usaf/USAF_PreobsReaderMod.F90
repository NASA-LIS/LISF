!
! MODULE: USAF_PreobsReaderMod
!
! DESCRIPTION: Contains code for reading USAF preobs files, performing simple
! preprocessing, and adding to gages database.
!
! AUTHOR: Eric Kemp, SSAI/NASA GSFC
!
!------------------------------------------------------------------------------

module USAF_PreobsReaderMod

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: USAF_read_preobs

  integer, parameter :: MISSING = -99999999

contains

  ! Read preobs files, perform simple preprocessing, and store
  ! in database.
  subroutine USAF_read_preobs(preobsdir, presavdir, &
       use_timestamp, &
       year, month, day, hour, use_expanded_station_ids, &
       alert_number)

    ! Imports
    use ESMF
    use LIS_coreMod, only: LIS_masterproc
    use LIS_logMod, only:  LIS_logunit, LIS_alert, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber
    use LIS_mpiMod, only:  LIS_mpi_comm
    use USAF_GagesMod, only: USAF_Gages_t

    ! Defaults
    implicit none

    ! Arguments
    character(*), intent(in) :: preobsdir
    character(*), intent(in) :: presavdir
    integer, intent(in) :: use_timestamp
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    integer, intent(in) :: use_expanded_station_ids
    integer, intent(inout) :: alert_number

    ! Locals
    character(255) :: filename
    integer, allocatable :: twfprc(:)
    integer, allocatable :: duration(:)
    integer, allocatable :: sixprc(:)
    integer, allocatable :: mscprc(:)
    integer, allocatable :: ilat(:)
    integer, allocatable :: ilon(:)
    integer, allocatable :: bsn(:)
    character(32), allocatable :: network(:)
    character(32), allocatable :: plat_id(:)
    character(2), allocatable :: wmocode_id(:)
    character(2), allocatable :: fipscode_id(:)
    integer, allocatable :: pastwx(:)
    integer, allocatable :: preswx(:)
    integer, allocatable :: wmoblk(:)
    integer :: twfprc_tmp
    integer :: duration_tmp
    integer :: sixprc_tmp
    integer :: mscprc_tmp
    integer :: ilat_tmp
    integer :: ilon_tmp
    character(14) :: YYYYMMDDhhmmss_tmp
    character(32) :: network_tmp
    character(32) :: plat_id_tmp
    character(2) :: wmocode_id_tmp, fipscode_id_tmp
    integer :: pastwx_tmp
    integer :: preswx_tmp
    integer :: wmoblk_tmp
    integer :: bsn_tmp
    integer :: nsize
    integer :: nsize_total
    integer :: ihemi
    logical :: found_file
    integer :: ierr
    integer :: stncnt
    integer :: i
    character(10) :: date10
    character(14), allocatable :: YYYYMMDDhhmmss(:)
    integer, allocatable :: amts24(:)
    integer, allocatable :: amts21(:)
    integer, allocatable :: amts18(:)
    integer, allocatable :: amts15(:)
    integer, allocatable :: amts12(:)
    integer, allocatable :: amts09(:)
    integer, allocatable :: amts06(:)
    integer, allocatable :: amts03(:)
    integer, allocatable :: amts02(:)
    integer, allocatable :: amts01(:)
    integer, allocatable :: amts00(:)
    type(USAF_Gages_t) :: obscur, obsprev
    integer :: rc
    character(255) :: presav_filename
    integer :: ipass, j
    logical :: exchanges
    type(ESMF_Time) :: curtime, prevtime, reporttime
    type(ESMF_TimeInterval) :: deltatime, maxdeltatime
    integer :: prevyear, prevmonth, prevday, prevhour
    logical :: file_exists
    integer :: deltahr
    character(10) :: prevdate10
    integer :: yyyy, mm, dd, h, m, s
    character(255) :: timestring
    integer :: iunit
    character(255) :: message(20)

    message = ''

    write(LIS_logunit,*)'---------------------------------------------'

    ! Set time limit allowed for report
    call esmf_timeintervalset(maxdeltatime, m=5, rc=rc)  ! 5 min after

    ! Find the total number of entries in the two hemispheric preobs
    ! files for this date/time.  This will be the upper-limit for how
    ! much memory we allocate for temporary storage.  If the newer
    ! global preobs file is read, logic in the loop will break out
    ! of the loop before the second pass.
    nsize_total = 0
    do ihemi = 1, 2

       call get_preobs_filename(filename, preobsdir, &
            use_timestamp, &
            ihemi, year, month, day, hour, use_expanded_station_ids)

       inquire(file=trim(filename), exist=found_file)
       if (.not. found_file) then
          write(LIS_logunit,*) '[WARN] Cannot find ', trim(filename)
          message(1) = '[WARN] Program:  LIS'
          message(2) = '  Routine: USAF_read_preobs'
          message(3) = '  Cannot find file ' // trim(filename)
          if (LIS_masterproc) then
             alert_number = alert_number + 1
             call LIS_alert('LIS.USAF_read_preobs', &
                  alert_number, message)
          end if
          if (use_expanded_station_ids == 1) exit ! These files are global
          cycle
       end if

       iunit = LIS_getNextUnitNumber()
       open(iunit, file=trim(filename), status='old', iostat=ierr)
       if (ierr .ne. 0) then
          write(LIS_logunit,*) '[WARN] Problem opening ', trim(filename)
          message(1) = '[WARN] Program:  LIS'
          message(2) = '  Routine: USAF_read_preobs'
          message(3) = '  Cannot open file ' // trim(filename)
          if (LIS_masterproc) then
             alert_number = alert_number + 1
             call LIS_alert('LIS.USAF_read_preobs', &
                  alert_number, message)
          end if
          if (use_expanded_station_ids == 1) exit ! These files are global
          cycle
       end if

       nsize = 0
       read(iunit, *, iostat=ierr) nsize
       if (ierr .ne. 0) then
          write(LIS_logunit,*) '[WARN] Problem reading ', trim(filename)
          message(1) = '[WARN] Program:  LIS'
          message(2) = '  Routine: USAF_read_preobs'
          message(3) = '  Problem reading file ' // trim(filename)
          if (LIS_masterproc) then
             alert_number = alert_number + 1
             call LIS_alert('LIS.USAF_read_preobs', &
                  alert_number, message)
          end if
          close(iunit)
          call LIS_releaseUnitNumber(iunit)
          if (use_expanded_station_ids == 1) exit ! These files are global
          cycle
       end if

       if (nsize == 0) then
          write(LIS_logunit,*)'[WARN] No precip obs found in ', &
               trim(filename)
          message(1) = '[WARN] Program:  LIS'
          message(2) = '  Routine: USAF_read_preobs'
          message(3) = '  No precip obs found in ' // trim(filename)
          if (LIS_masterproc) then
             alert_number = alert_number + 1
             call LIS_alert('LIS.USAF_read_preobs', &
                  alert_number, message)
          end if
       else
          write(LIS_logunit,*) '[INFO] Will process ', trim(filename)
       end if

       nsize_total = nsize_total + nsize
       close(iunit)
       call LIS_releaseUnitNumber(iunit)

       if (use_expanded_station_ids == 1) exit ! These files are global
    end do

    if (nsize_total == 0) then
       write(LIS_logunit,*) '[WARN] No precip obs available!'
       return
    end if

    ! We now have an upper limit of how many gage reports to save.
    allocate(twfprc(nsize_total))
    twfprc = MISSING
    allocate(duration(nsize_total))
    duration = MISSING
    allocate(sixprc(nsize_total))
    sixprc = MISSING
    allocate(mscprc(nsize_total))
    mscprc = MISSING
    allocate(ilat(nsize_total))
    ilat = MISSING
    allocate(ilon(nsize_total))
    ilon = MISSING
    allocate(bsn(nsize_total))
    bsn = MISSING
    allocate(network(nsize_total))
    network = "NULL"
    allocate(plat_id(nsize_total))
    plat_id = "NULL"
    allocate(wmocode_id(nsize_total))
    wmocode_id = "??"
    allocate(fipscode_id(nsize_total))
    fipscode_id = "??"
    allocate(pastwx(nsize_total))
    pastwx = MISSING
    allocate(preswx(nsize_total))
    preswx = MISSING
    allocate(wmoblk(nsize_total))
    wmoblk = MISSING
    allocate(YYYYMMDDhhmmss(nsize_total))
    YYYYMMDDhhmmss = "NULL"

    write(date10, "(i4.4,i2.2,i2.2,i2.2)") year, month, day, hour

    ! Need to set current time in ESMF for calculations
    call esmf_timeset(curtime, yy=year, mm=month, dd=day, h=hour, &
         m=0, s=0, rc=rc)

    ! Now begin reading the data in the files, performing simple checks
    ! along the way.
    stncnt = 0
    do ihemi = 1, 2
       call get_preobs_filename(filename, preobsdir, &
            use_timestamp, &
            ihemi, year, month, day, hour, use_expanded_station_ids)

       inquire(file=trim(filename), exist=found_file)
       if (.not. found_file) then
          if (use_expanded_station_ids == 1) exit
          cycle
       end if

       iunit = LIS_getNextUnitNumber()
       open(iunit, file=trim(filename), status='old', iostat=ierr)
       if (ierr .ne. 0) cycle

       nsize = 0
       read(iunit, *, iostat=ierr) nsize
       if (ierr .ne. 0) then
          close(iunit)
          call LIS_releaseUnitNumber(iunit)
          if (use_expanded_station_ids == 1) exit
          cycle
       end if

       ! Start reading the actual obs
       do i = 1, nsize
          if (use_expanded_station_ids == 1) then
             twfprc_tmp = MISSING
             sixprc_tmp = MISSING
             read(iunit, 6001, iostat=ierr) YYYYMMDDhhmmss_tmp, &
                  network_tmp, plat_id_tmp, ilat_tmp, ilon_tmp, &
                  wmocode_id_tmp, fipscode_id_tmp, wmoblk_tmp, &
                  mscprc_tmp, duration_tmp, pastwx_tmp, preswx_tmp
6001         format(1x, a14, 1x, a10, 1x, a32, 1x, i9, 1x, i9, 1x, &
                  a2, 1x, a2, 1x, i9, 1x, i9, 1x, i9, 1x, i9, 1x, i9)
             if (wmocode_id_tmp == "  ") wmocode_id_tmp = "??"
             if (fipscode_id_tmp == "  ") fipscode_id_tmp = "??"
          else
             read(iunit, 6000, iostat=ierr) twfprc_tmp, duration_tmp, &
                  sixprc_tmp, &
                  mscprc_tmp, ilat_tmp, ilon_tmp, network_tmp, &
                  plat_id_tmp, &
                  pastwx_tmp, preswx_tmp, wmoblk_tmp
6000         format(i10, i10, i10, i10, i10, i10, 1x, a10, 2x, a10 &
                  ,i10, i10, i10)
             wmocode_id_tmp = "??"
             fipscode_id_tmp = "??"
             YYYYMMDDhhmmss_tmp = date10 // "0000"
          end if

          if (ierr .ne. 0) then
             write(LIS_logunit,*) &
                  '[WARN] Problem reading report, skipping line...'
             cycle
          end if

          ! Skip if lat/lon is 0 (this is interpreted as missing).
          if (ilat_tmp == 0 .and. ilon_tmp == 0) cycle

          ! Skip reports that are too much after the analysis time
          ! (but allow earlier reports).  This is a crude way of
          ! allowing for Australian reports that are sometimes one or
          ! two hours behind the sub-synoptic times.
          if (use_expanded_station_ids == 1) then
             read(YYYYMMDDhhmmss_tmp, &
                  '(I4.4, I2.2, I2.2, I2.2, I2.2, I2.2)') yyyy, mm, dd, &
                  h, m, s
             call esmf_timeset(reporttime, yy=yyyy, mm=mm, dd=dd, &
                  h=h, m=m, s=s, rc=rc)
             deltatime = reporttime - curtime

             call esmf_timeintervalget(deltatime, timestring=timestring, &
                  rc=rc)

          !    if (deltatime < mindeltatime) cycle
             if (deltatime > maxdeltatime) cycle
          end if

          ! Don't save 1-hr precip.
          ! FIXME...Add special handling in new preobs files to
          ! estimate 3-hr accums if enough 1-hr reports are available.
          if (duration_tmp .eq. 1) then
             mscprc_tmp = MISSING
             duration_tmp = MISSING
          end if

          ! Skip if no useable data in report.
          if (twfprc_tmp == MISSING .and. &
               sixprc_tmp == MISSING .and. &
               mscprc_tmp == MISSING .and. &
               (pastwx_tmp == MISSING .or. &
               preswx_tmp == MISSING)) then
             cycle
          end if

          ! If we don't have FIPS or WMO codes available, try guessing
          ! the FIPS code for countries with special rules based on
          ! the WIGOS issuer number.
          if (fipscode_id_tmp == "??" .and. &
               wmocode_id_tmp == "??") then
             fipscode_id_tmp = get_fips_from_wigos_issuer(plat_id_tmp)
          end if

          ! Set the numeric bsn field.
          if (network_tmp .eq. "WMO") then
             bsn_tmp = set_bsn_wmo(use_expanded_station_ids, &
                  plat_id_tmp, fipscode_id_tmp, wmocode_id_tmp)
          else if (network_tmp .eq. "CDMS") then
             read (plat_id_tmp, '(i32)') bsn_tmp
          else
             bsn_tmp = 0
          end if

          ! Set the country codes if possible.
          if (network_tmp .eq. "CANA") then
             wmocode_id_tmp = "CA"
             fipscode_id_tmp = "CA"
          else if (network_tmp .eq. "FAA") then
             wmocode_id_tmp = "US"
             fipscode_id_tmp = "US"
          end if

          ! Sometimes 6 and 24-hr data are stored in misc.  Copy
          ! over if necessary, then erase.
          if (duration_tmp == 24) then
             if (twfprc_tmp == MISSING) then
                twfprc_tmp = mscprc_tmp
             end if
             mscprc_tmp = MISSING
             duration_tmp = MISSING
          else if (duration_tmp == 6) then
             if (sixprc_tmp == MISSING) then
                sixprc_tmp = mscprc_tmp
             end if
             mscprc_tmp = MISSING
             duration_tmp = MISSING
          end if

          ! If this is the first report, save it.
          if (stncnt == 0) then
             stncnt = stncnt + 1
             twfprc(stncnt) = twfprc_tmp
             duration(stncnt) = duration_tmp
             sixprc(stncnt) = sixprc_tmp
             mscprc(stncnt) = mscprc_tmp
             ilat(stncnt) = ilat_tmp
             ilon(stncnt) = ilon_tmp
             bsn(stncnt) = bsn_tmp
             network(stncnt) = network_tmp
             plat_id(stncnt) = plat_id_tmp
             wmocode_id(stncnt) = wmocode_id_tmp
             fipscode_id(stncnt) = fipscode_id_tmp
             pastwx(stncnt) = pastwx_tmp
             preswx(stncnt) = preswx_tmp
             wmoblk(stncnt) = wmoblk_tmp
             YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
             cycle
          else
             ! Look at the most recent saved gage report, and see if this
             ! new entry is from the same station.
             if (network(stncnt) .eq. network_tmp .and. &
                  trim(plat_id(stncnt)) .eq. trim(plat_id_tmp)) then

                ! For newer file format with date/time with each report,
                ! pick the more recent date/time for duration other
                ! than 1 hour.
                if (YYYYMMDDhhmmss(stncnt) < YYYYMMDDhhmmss_tmp .and. &
                     (twfprc_tmp .ne. MISSING .or. &
                     sixprc_tmp .ne. MISSING .or. &
                     mscprc_tmp .ne. MISSING .or. &
                     (preswx_tmp .ne. MISSING .and. &
                      pastwx_tmp .ne. MISSING))) then
                   twfprc(stncnt) = twfprc_tmp
                   duration(stncnt) = duration_tmp
                   sixprc(stncnt) = sixprc_tmp
                   mscprc(stncnt) = mscprc_tmp
                   ilat(stncnt) = ilat_tmp
                   ilon(stncnt) = ilon_tmp
                   bsn(stncnt) = bsn_tmp
                   network(stncnt) = network_tmp
                   plat_id(stncnt) = plat_id_tmp
                   wmocode_id(stncnt) = wmocode_id_tmp
                   fipscode_id(stncnt) = fipscode_id_tmp
                   pastwx(stncnt) = pastwx_tmp
                   preswx(stncnt) = preswx_tmp
                   wmoblk(stncnt) = wmoblk_tmp
                   YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                   cycle
                end if

                ! Logic for older file format.  Code above assigns
                ! the file date/time as the report date/time.
                ! See if there is any additional information worth
                ! saving.  The AGRMET logic overwrites the older data
                ! with the newer if the newer data is a higher amount.
                ! An exception is that "zero" six-hour precip is ignored
                ! if 24-hr or 12-hr precip is positive, apparently
                ! because "zero" six-hour precip reports in this case are
                ! typically from an asynoptic ob. We follow that logic
                ! here.
                if (YYYYMMDDhhmmss(stncnt) == YYYYMMDDhhmmss_tmp) then
                   if (twfprc(stncnt) < twfprc_tmp) then
                      twfprc(stncnt) = twfprc_tmp
                      YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                   end if

                   if (mscprc(stncnt) .eq. MISSING .and. &
                        mscprc_tmp .ne. MISSING) then
                      ! We don't have a misc precip stored yet, so do it
                      ! here.
                      mscprc(stncnt) = mscprc_tmp
                      duration(stncnt) = duration_tmp
                      YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                   else if (duration_tmp > duration(stncnt) .and. &
                        duration_tmp > 1 .and. &
                        mscprc_tmp .ne. MISSING) then
                      ! Try to capture 3-hr or other longer duration
                      ! precip.
                      mscprc(stncnt) = mscprc_tmp
                      duration(stncnt) = duration_tmp
                      YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                   else if (mscprc_tmp > mscprc(stncnt) .and. &
                        duration_tmp == duration(stncnt)) then
                      ! Somewhat mimics AGRMET.  Save larger accum for
                      ! current duration.
                      mscprc(stncnt) = mscprc_tmp
                      duration(stncnt) = duration_tmp
                      YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                   else if ((mscprc_tmp > mscprc(stncnt)) .and. &
                        (duration_tmp == MISSING .or. &
                        duration_tmp == 0) .and. &
                        (duration(stncnt) == MISSING .or. &
                        duration(stncnt) == 0)) then
                      ! Somewhat mimics AGRMET.  Save larger accum for
                      ! undefined duration, for use in certain regions
                      ! (e.g., Russia)
                      mscprc(stncnt) = mscprc_tmp
                      duration(stncnt) = duration_tmp
                      YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                   end if
                   if ( (twfprc(stncnt) > 0) .or. &
                        (duration(stncnt) == 12 .and. &
                        mscprc(stncnt) > 0) ) then
                      if (sixprc(stncnt) == MISSING .and. &
                           sixprc_tmp == 0) then
                         sixprc(stncnt) = MISSING
                         YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                      else if (sixprc(stncnt) == 0 .and. &
                           sixprc_tmp == MISSING) then
                         sixprc(stncnt) = MISSING
                         YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                      else
                         if (sixprc(stncnt) < sixprc_tmp) then
                            sixprc(stncnt) = sixprc_tmp
                            YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                         end if
                      end if
                   else
                      if (sixprc(stncnt) < sixprc_tmp) then
                         sixprc(stncnt) = sixprc_tmp
                         YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                      end if
                   end if
                end if
             else
                ! This isn't the same station, so save it.
                !if (trim(network_tmp) == "AMIL" .and. &
                !     trim(plat_id_tmp) == "KQTC") then
                !   write(LIS_logunit,*)'EMK: KQTC, twfprc = ', twfprc_tmp
                !end if

                stncnt = stncnt + 1
                twfprc(stncnt) = twfprc_tmp
                duration(stncnt) = duration_tmp
                sixprc(stncnt) = sixprc_tmp
                mscprc(stncnt) = mscprc_tmp
                ilat(stncnt) = ilat_tmp
                ilon(stncnt) = ilon_tmp
                bsn(stncnt) = bsn_tmp
                network(stncnt) = network_tmp
                plat_id(stncnt) = plat_id_tmp
                wmocode_id(stncnt) = wmocode_id_tmp
                fipscode_id(stncnt) = fipscode_id_tmp
                pastwx(stncnt) = pastwx_tmp
                preswx(stncnt) = preswx_tmp
                wmoblk(stncnt) = wmoblk_tmp
                YYYYMMDDhhmmss(stncnt) = YYYYMMDDhhmmss_tmp
                cycle
             end if
          end if
       end do

       if (use_expanded_station_ids == 1) cycle ! These files are global
    end do ! ihemi

    ! Since we combined both the NH and SH files, the resulting list
    ! is unsorted.
    ! Bubble sort the reports by network.
    if (use_expanded_station_ids == 0) then
       ipass = 1
       exchanges = .true.
       do while ((ipass < stncnt) .and. exchanges)
          exchanges = .false.
          do j = 1, stncnt - ipass
             if (network(j) > network(j+1)) then
                call swap_int(twfprc(j), twfprc(j+1))
                call swap_int(duration(j), duration(j+1))
                call swap_int(sixprc(j), sixprc(j+1))
                call swap_int(mscprc(j), mscprc(j+1))
                call swap_int(ilat(j), ilat(j+1))
                call swap_int(ilon(j), ilon(j+1))
                call swap_int(bsn(j), bsn(j+1))
                call swap_char32(network(j), network(j+1))
                call swap_char32(plat_id(j), plat_id(j+1))
                call swap_char2(wmocode_id(j), wmocode_id(j+1))
                call swap_char2(fipscode_id(j), fipscode_id(j+1))
                call swap_int(pastwx(j), pastwx(j+1))
                call swap_int(preswx(j), preswx(j+1))
                call swap_int(wmoblk(j), wmoblk(j+1))
                call swap_char14(YYYYMMDDhhmmss(j), YYYYMMDDhhmmss(j+1))
                exchanges = .true.
             end if
          end do
          ipass = ipass + 1
       end do

       ! Now bubble sort by station ID within the same network
       ipass = 1
       exchanges = .true.
       do while ((ipass < stncnt) .and. exchanges)
          exchanges = .false.
          do j = 1, stncnt - ipass
             if ((network(j) == network(j+1)) .and. &
                  (plat_id(j) > plat_id(j+1))) then
                call swap_int(twfprc(j), twfprc(j+1))
                call swap_int(duration(j), duration(j+1))
                call swap_int(sixprc(j), sixprc(j+1))
                call swap_int(mscprc(j), mscprc(j+1))
                call swap_int(ilat(j), ilat(j+1))
                call swap_int(ilon(j), ilon(j+1))
                call swap_int(bsn(j), bsn(j+1))
                call swap_char32(network(j), network(j+1))
                call swap_char32(plat_id(j), plat_id(j+1))
                call swap_char2(wmocode_id(j), wmocode_id(j+1))
                call swap_char2(fipscode_id(j), fipscode_id(j+1))
                call swap_int(pastwx(j), pastwx(j+1))
                call swap_int(preswx(j), preswx(j+1))
                call swap_int(wmoblk(j), wmoblk(j+1))
                call swap_char14(YYYYMMDDhhmmss(j), YYYYMMDDhhmmss(j+1))
                exchanges = .true.
             end if
          end do
          ipass = ipass + 1
       end do
    end if

    ! Report if any WMO obs are not legacy 5-digit SYNOP
    do j = 1, stncnt
       if (network(j) .ne. "WMO") cycle
       if (.not. len_trim(plat_id(j)) == 5) then
          write(LIS_logunit,*) '[INFO] Found non SYNOP WMO station ', &
               trim(plat_id(j))
       else if (verify(trim(plat_id(j)), "0123456789") .ne. 0) then
          write(LIS_logunit,*) '[INFO] Found non SYNOP WMO station ', &
               trim(plat_id(j))
       end if
    end do

    ! Copy into Gage_type object
    allocate(amts24(nsize_total)) ; amts24 = twfprc
    allocate(amts21(nsize_total)) ; amts21 = MISSING
    allocate(amts18(nsize_total)) ; amts18 = MISSING
    allocate(amts15(nsize_total)) ; amts15 = MISSING
    allocate(amts12(nsize_total)) ; amts12 = MISSING
    allocate(amts09(nsize_total)) ; amts09 = MISSING
    allocate(amts06(nsize_total)) ; amts06 = sixprc
    allocate(amts03(nsize_total)) ; amts03 = MISSING
    allocate(amts02(nsize_total)) ; amts02 = MISSING
    allocate(amts01(nsize_total)) ; amts01 = MISSING
    allocate(amts00(nsize_total)) ; amts00 = mscprc

    call obscur%new(date10, stncnt, YYYYMMDDhhmmss, network, plat_id, &
         wmocode_id, fipscode_id, wmoblk, bsn, ilat, ilon, &
         amts24, amts21, amts18, amts15, amts12, amts09, amts06, &
         amts03, amts02, amts01, amts00, duration, preswx, pastwx)

    ! Clean up the temporary arrays
    deallocate(twfprc)
    deallocate(duration)
    deallocate(sixprc)
    deallocate(mscprc)
    deallocate(ilat)
    deallocate(ilon)
    deallocate(bsn)
    deallocate(network)
    deallocate(plat_id)
    deallocate(wmocode_id)
    deallocate(fipscode_id)
    deallocate(pastwx)
    deallocate(preswx)
    deallocate(wmoblk)
    deallocate(YYYYMMDDhhmmss)
    deallocate(amts24)
    deallocate(amts21)
    deallocate(amts18)
    deallocate(amts15)
    deallocate(amts12)
    deallocate(amts09)
    deallocate(amts06)
    deallocate(amts03)
    deallocate(amts02)
    deallocate(amts01)
    deallocate(amts00)

    ! Check for unphysical values.
    call obscur%check_gross_errors()

    ! Reconcile values in same report (shorter durations should not be
    ! larger than longer durations)
    call obscur%reconcile_self()

    ! Fill gaps if possible
    call obscur%fill_gaps()

    ! Use the miscellaneous precip, with special rules for some
    ! countries.
    call obscur%use_misc_precip()

    ! Leverage present and past weather information to identify
    ! zero precip durations.
    !EMK...Disable until more accurate inspection of time in new
    !preobs files is implemented.
    !call obscur%use_preswx_pastwx()

    ! Try fetching presav files from earlier hours and reconcile.
    write(LIS_logunit,*)'[INFO] Will compare with earlier gage reports.'
    do deltahr = 3, 21, 3
       call esmf_timeintervalset(deltatime, yy=0, mm=0, d=0, &
            h=deltahr, m=0, s=0, rc=rc)
       prevtime = curtime - deltatime
       call esmf_timeget(prevtime, yy=prevyear, mm=prevmonth, &
            dd=prevday, h=prevhour, rc=rc)
       write(presav_filename,'(A,A,i4.4,i2.2,i2.2,i2.2)') &
            trim(presavdir), '/presav2.03hr.', &
            prevyear, prevmonth, prevday, prevhour
       inquire(file=presav_filename, exist=file_exists)
       if (file_exists) then
          write(LIS_logunit,*)'[INFO] Comparing against data in ', &
               trim(presav_filename)
          write(prevdate10,'(i4.4,i2.2,i2.2,i2.2)') &
               prevyear, prevmonth, prevday, prevhour
          call obsprev%read_data(presav_filename, prevdate10, &
               alert_number)
          if (deltahr == 3) then
             call obscur%reconcile_gages03(obsprev)
          else if (deltahr == 6) then
             call obscur%reconcile_gages06(obsprev)
          else if (deltahr == 9) then
             call obscur%reconcile_gages09(obsprev)
          else if (deltahr == 12) then
             call obscur%reconcile_gages12(obsprev)
          else if (deltahr == 15) then
             call obscur%reconcile_gages15(obsprev)
          else if (deltahr == 18) then
             call obscur%reconcile_gages18(obsprev)
          else if (deltahr == 21) then
             call obscur%reconcile_gages21(obsprev)
          end if
          call obsprev%delete()
       else
          write(LIS_logunit,*)'[WARN] Cannot find ', trim(presav_filename)
          write(LIS_logunit,*) &
               '[WARN] Will skip reconciling with obs from ', &
               abs(deltahr),' hours ago'
          message(1) = '[WARN] Program:  LIS'
          message(2) = '  Routine: USAF_read_preobs'
          message(3) = '  Cannot find earlier presav2 file ' // &
               trim(presav_filename)
          message(4) = ' Observation count will be reduced'
          if (LIS_masterproc) then
             alert_number = alert_number + 1
             call LIS_alert('LIS.USAF_read_preobs', &
                  alert_number, message)
          end if
       end if
    end do

    ! Correct for overnight reporting issues in South America
    write(LIS_logunit,*) &
        '[INFO] Correcting overnight reporting in South America'
    call obscur%correct_region3_12z()

    ! Have the master process write the data to file.
    if (LIS_masterproc) then
       write(presav_filename,'(A,A,i4.4,i2.2,i2.2,i2.2)') &
            trim(presavdir), '/presav2.03hr.', &
            year, month, day, hour
       write(LIS_logunit,*)'[INFO] Writing to ', trim(presav_filename)
       call obscur%write_data(presav_filename)
    end if
#if (defined SPMD)
    call MPI_Barrier(LIS_mpi_comm, ierr)
#endif

    ! Clean up
    call obscur%delete()

  end subroutine USAF_read_preobs

  ! Generate name of preobs file.  Based on code in AGRMET_getpcpobs.F90
  subroutine get_preobs_filename(filename, preobsdir, &
       use_timestamp, &
       ihemi, year, month, day, hour, use_expanded_station_ids)

    ! Defaults
    implicit none

    ! Arguments
    character(*), intent(out) :: filename
    character(*), intent(in) :: preobsdir
    integer, intent(in) :: use_timestamp
    integer, intent(in) :: ihemi
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    integer, intent(in) :: use_expanded_station_ids

    ! Locals
    character(8) :: ftime1
    character(10) :: ftime2
    character(2), parameter :: FHEMI(2) = (/'nh', 'sh'/)

    write(ftime2, '(i4.4, i2.2, i2.2, i2.2)') year, month, day, hour

    if (use_expanded_station_ids == 0) then
       !if (use_timestamp == 1) then
       !   write(ftime1, '(i4.4, i2.2, i2.2)') &
       !        year, month, day
       !   filename = ftime1 // '/' // trim(preobsdir) // &
       !        '/preobs_' // FHEMI(ihemi) // '.03hr.' // ftime2
       !else
          filename = trim(preobsdir) // &
               '/preobs_' // FHEMI(ihemi) // '.03hr.' // ftime2
       !endif
    else if (use_expanded_station_ids == 1) then
       !if (use_timestamp == 1) then
       !   write(ftime1, '(i4.4, i2.2, i2.2)') &
       !        year, month, day
       !   filename = ftime1 // '/' // trim(preobsdir) // &
       !        '/preobs_03hr_' // ftime2 // ".txt"
       !else
          filename = trim(preobsdir) // &
               '/preobs_03hr_' // ftime2 // ".txt"
       !endif
    end if
  end subroutine get_preobs_filename

  ! Set block number for a WMO station
  function set_bsn_wmo(use_expanded_station_id, plat_id, fipscode_id, &
       wmocode_id) result(bsn)

    ! Imports
    use USAF_GagesMod, only: JMOBS_SWEDEN_LOWLIMIT, &
         JMOBS_DENMARK_LOWLIMIT, JMOBS_RUSSIA_LOWLIMIT, &
         JMOBS_INDIA_LOWLIMIT, JMOBS_SRILANKA_LOWLIMIT, &
         JMOBS_S_AMER_LOWLIMIT

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: use_expanded_station_id
    character(32), intent(in) :: plat_id
    character(2), intent(in) :: fipscode_id
    character(2), intent(in) :: wmocode_id

    ! Return variable
    integer :: bsn

    bsn = 0

    ! Old preobs files don't have any country codes.  Just copy the
    ! 5-digit ID and use as is.
    if (use_expanded_station_id == 0) then
       read(plat_id, '(i5.5)') bsn
       return
    end if

    ! If this is a legacy 5-digit WMO ID, just copy it.
    if (verify(trim(plat_id), "0123456789") == 0 &
         .and. len_trim(plat_id) == 5) then
       read(plat_id, '(i5.5)') bsn
       return
    end if

    ! WIGOS supports legacy 5-digit WMO IDs with a special prefix
    if (len_trim(plat_id) == 15) then
       if (plat_id(1:10) == "0-20000-0-" .and. &
            verify(trim(plat_id(11:15)), "0123456789") == 0) then
          read(plat_id(11:15), '(i5.5)') bsn
          return
       end if
    end if

    ! Try using the FIPS country code.  This code seems to be more
    ! accurate than the WMO code.
    select case(fipscode_id)
    case('SW')
       bsn = JMOBS_SWEDEN_LOWLIMIT
       return
    case('DA')
       bsn = JMOBS_DENMARK_LOWLIMIT
       return
    case('AM', 'AJ', 'BO', 'EN', 'GG', 'LG', 'LH', 'KZ', &
         'KG', 'MD', 'RS', 'TI', 'TX', 'UP', 'UZ')
       ! AM is Armenia
       ! AJ is Azerbaijan
       ! BO is Belarus
       ! EN is Estonia
       ! GG is Georgia
       ! LG is Latvia
       ! LH is Lituania
       ! KZ is Kazakhstan
       ! KG is Kyrgyzstan
       ! MD is Moldova
       ! RS is Russia
       ! TI is Tajikistan
       ! TX is Turkmenistan
       ! UP is Ukraine
       ! UZ is Uzbekistan
       bsn = JMOBS_RUSSIA_LOWLIMIT
       return
    case ("IN")
       bsn = JMOBS_INDIA_LOWLIMIT
       return
    case ("CE")
       bsn = JMOBS_SRILANKA_LOWLIMIT
       return
    case ('AR', 'BL', 'BR', 'CI', 'CO', 'EC', 'FK', 'FG', 'GY', &
         'PA', 'PE', 'SX', 'NS', 'UY', 'VE')
       ! AR is Argentina
       ! BL is Bolivia
       ! BR is Brazil
       ! CI is Chile
       ! CO is Colombia
       ! EC is Ecuador
       ! FK is Falkland Islands
       ! FG is French Guiana
       ! GY is Guyana
       ! PA is Paraguay
       ! PE is Peru
       ! SX is South Georgia and the South Sandwich Islands
       ! NS is Suriname
       ! UY is Uruguay
       ! VE is Venezuela
       bsn = JMOBS_S_AMER_LOWLIMIT
       return
    case default
       bsn = 0
    end select

    ! Try using the WMO country code
    select case(wmocode_id)
    case('SE')
       bsn = JMOBS_SWEDEN_LOWLIMIT
       return
    case('DK')
       bsn = JMOBS_DENMARK_LOWLIMIT
       return
    case('AM', 'AZ', 'BY', 'EE', 'GE', 'LV', 'LT', 'KZ', 'KG', &
         'MD', 'RU', 'TJ', 'TM', 'UA', 'UZ')
       ! AM is Armenia
       ! AZ is Azerbaijan
       ! BY is Belarus
       ! EE is Estonia
       ! GE is Georgia
       ! LV is Latvia
       ! LT is Lithuania
       ! KZ is Kazakhstan
       ! KG is Kyrgyzstan
       ! MD is Moldova
       ! RU is Russia
       ! TJ is Tajikistan
       ! TM is Turkmenistan
       ! UA is Ukraine
       ! UZ is Uzbekistan
       bsn = JMOBS_RUSSIA_LOWLIMIT
       return
    case ("IN")
       bsn = JMOBS_INDIA_LOWLIMIT
       return
    case ("LK")
       bsn = JMOBS_SRILANKA_LOWLIMIT
       return
    case ('AR', 'BO', 'BR', 'CL', 'CO', 'EC', 'FK', 'GF', 'GY', &
         'PY', 'PE', 'GS', 'SR', 'UG', 'VE')
       ! AR is Argentina
       ! BO is Bolivia
       ! BR is Brazil
       ! CL is Chile
       ! CO is Colombia
       ! EC is Ecuador
       ! FK is Falkland Islands
       ! GF is French Guiana
       ! GY is Guyana
       ! PY is Paraguay
       ! PE is Peru
       ! GS is South Georgia and the South Sandwich Islands
       ! SR is Suriname
       ! UG is Uruguay
       ! VE is Venezuela
       bsn = JMOBS_S_AMER_LOWLIMIT
       return
    case default
       bsn = 0
    end select

  end function set_bsn_wmo

  subroutine swap_int(var1, var2)
    implicit none
    integer, intent(inout) :: var1
    integer, intent(inout) :: var2
    integer :: tmp
    tmp = var1
    var1 = var2
    var2 = tmp
  end subroutine swap_int

  subroutine swap_char2(var1, var2)
    implicit none
    character(2), intent(inout) :: var1
    character(2), intent(inout) :: var2
    character(2) :: tmp
    tmp = var1
    var1 = var2
    var2 = tmp
  end subroutine swap_char2

  subroutine swap_char10(var1, var2)
    implicit none
    character(10), intent(inout) :: var1
    character(10), intent(inout) :: var2
    character(10) :: tmp
    tmp = var1
    var1 = var2
    var2 = tmp
  end subroutine swap_char10

  subroutine swap_char14(var1, var2)
    implicit none
    character(14), intent(inout) :: var1
    character(14), intent(inout) :: var2
    character(14) :: tmp
    tmp = var1
    var1 = var2
    var2 = tmp
  end subroutine swap_char14

  subroutine swap_char32(var1, var2)
    implicit none
    character(32), intent(inout) :: var1
    character(32), intent(inout) :: var2
    character(32) :: tmp
    tmp = var1
    var1 = var2
    var2 = tmp
  end subroutine swap_char32

  ! Try to set FIPS country code based on WIGOS issuer number.  Only
  ! certain specific countries are checked for (for special accumulation
  ! rules).
  function get_fips_from_wigos_issuer(stn) result (fips)

    ! Defaults
    implicit none

    ! Arguments
    character(*),intent(in) :: stn

    ! Return variable
    character(2) :: fips

    ! Locals
    integer :: idx1, idx2
    integer :: id

    fips = '??' ! First guess

    ! See if this is a WIGOS station ID; if so, the issuer code is the
    ! integer string between the first two dashes.  If this is not
    ! WIGOS, just return.
    idx1 = scan(stn,'-')
    if (idx1 == 0) return
    idx2 = scan(stn(idx1+1:len(stn)),'-')
    if (idx2 == 0) return
    if (verify(trim(stn(idx1+1:idx2+1)), "0123456789") .ne. 0) return
    read(stn(idx1+1:idx2+1), '(i5.5)') id

    ! Select FIPS code for specific countries where we have special rules.
    select case(id)
    case(752)
       fips = 'SW' ! Sweden
    case(208)
       fips = 'DA' ! Denmark
    case(051) ! Post-Soviet countries start here
       fips = 'AM' ! Armenia
    case(031)
       fips = 'AJ' ! Azerbaijan
    case(112)
       fips = 'BO' ! Belarus
    case(233)
       fips = 'EN' ! Estonia
    case(268)
       fips = 'GG' ! Georgia
    case(428)
       fips = 'LG' ! Latvia
    case(440)
       fips = 'LH' ! Lithuania
    case(398)
       fips = 'KZ' ! Kazakhstan
    case(417)
       fips = 'KG' ! Kyrgyzstan
    case(498)
       fips = 'MD' ! Moldova
    case(643)
       fips = 'RS' ! Russia
    case(762)
       fips = 'TI' ! Tajikistan
    case(795)
       fips = 'TX' ! Turkmenistan
    case(804)
       fips = 'UP' ! Ukraine
    case(860)
       fips = 'UZ' ! Uzbekistan
    case(356) ! Post-Soviet countries end here
       fips = 'IN' ! India
    case(144)
       fips = 'CE' ! Sri Lanka
    case(032) ! South American countries start here
       fips = 'AR' ! Argentina
    case(068)
       fips = 'BL' ! Bolivia
    case(076)
       fips = 'BR' ! Brazil
    case(152)
       fips = 'CI' ! Chile
    case(170)
       fips = 'CO' ! Columbia
    case(218)
       fips = 'EC' ! Ecuador
    case(238)
       fips = 'FK' ! Falkland Islands
    case(254)
       fips = 'FG' ! French Guiana
    case(328)
       fips = 'GY' ! Guyana
    case(600)
       fips = 'PA' ! Paraguay
    case(604)
       fips = 'PE' ! Peru
    case(239)
       fips = 'SX' ! South Georgia and the South Sandwich Islands
    case(740)
       fips = 'NS' ! Suriname
    case(858)
       fips = 'UY' ! Uruguay
    case(862)
       fips = 'VE' ! Venezuela
    end select
  end function get_fips_from_wigos_issuer

end module USAF_PreobsReaderMod
