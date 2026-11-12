!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readSMOPSsmObs
! \label{readSMOPSsmObs}
!
! !REVISION HISTORY:
!  8 May 2013: Sujay Kumar, Initial Specification
!  22 Sep 2017: Mahdi Navari, filter added to screen the the ASCAT data
!               in the arid and semi-arid regions
!  9 Oct 2017: Mahdi Navari, added capability to read AMSR2 and SMAP
!  28 Oct 2026: Ehsan Jalilvand, read the 6-hourly SMOPS archive (three
!               naming conventions, as in the LIS SMOPS ASCAT reader) and
!               SMOPS v4 (ASCAT-B and ASCAT-C). The GRIB2 reading and the
!               interpolation are split so that the 6-hourly windows of a
!               time step are combined on the SMOPS grid and interpolated
!               once.
!  12 Nov 2026: Ehsan Jalilvand, SMOPS v5: any NPR-SMOPS-CMAP-6hrly
!               version is found (the newest is read), the version of
!               these files is taken from the name, and LDT stops on a
!               version whose GRIB2 parameter numbers it does not have.
!  12 Nov 2026: Ehsan Jalilvand, check the GRIB2 grid of each layer; flip
!               the rows only if they are stored north to south.
!
! !INTERFACE:
subroutine readSMOPSsmObs(n)
! !USES:
  use LDT_coreMod,      only : LDT_rc
  use LDT_timeMgrMod,   only : LDT_date2time, LDT_tick
  use LDT_logMod,       only : LDT_logunit
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use SMOPSsm_obsMod,   only : SMOPSsmobs

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!
! This subroutine provides the data reader for the SMOPS
! soil moisture retrieval product. The data has many layers
! and the reader provides the options to select the
! desired layer(s).
!
! Daily files (default):
!   <dir>/YYYY/NPR_SMOPS_CMAP_DYYYYMMDD.gr2, one per day.
!
! 6-hourly files ("SMOPS soil moisture use realtime data: 1"):
!   every 6-hourly window that starts within [current time, current
!   time + time step) is read. On the SMOPS grid, a later window replaces
!   an earlier one where it has a valid retrieval, so a 1-day time step
!   gives one value per grid cell per day (the latest), as in the NESDIS
!   daily product. The combined field is then interpolated once.
!
!EOP

  logical           :: file_exists, ok, found
  integer           :: c, r, iw, nwin, nnat
  integer           :: yr, mo, da, hr, mn, ss, doy
  real              :: gmt
  real*8            :: wtime
  character(len=LDT_CONST_PATH_LEN) :: fname
  character(len=16) :: fver
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real,      allocatable :: sm_nat(:), sm_win(:)
  logical*1, allocatable :: b_nat(:), b_win(:)

  SMOPSsmobs(n)%smobs = LDT_rc%udef
  smobs = LDT_rc%udef

  nnat = SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr
  allocate(sm_nat(nnat), b_nat(nnat))
  sm_nat = LDT_rc%udef
  b_nat = .false.
  found = .false.

  yr = LDT_rc%yr
  mo = LDT_rc%mo
  da = LDT_rc%da
  hr = LDT_rc%hr
  mn = 0
  ss = 0
  call LDT_date2time(wtime, doy, gmt, yr, mo, da, hr, mn, ss)

  if (SMOPSsmobs(n)%useRealtime == 1) then
     allocate(sm_win(nnat), b_win(nnat))
     nwin = nint(LDT_rc%ts / 21600.0)
     do iw = 1, nwin
        call create_SMOPSsm_6hrly_filename(n, SMOPSsmobs(n)%odir, &
             SMOPSsmobs(n)%conv, yr, mo, da, hr, fname, fver)
        if (trim(fname) /= 'NULL') then
           inquire(file=trim(fname), exist=file_exists)
           if (file_exists) then
              write(LDT_logunit,*) '[INFO] Reading ', trim(fname)
              call read_SMOPS_native(n, fname, fver, wtime, sm_win, b_win, ok)
              if (ok) then
                 where (b_win)
                    sm_nat = sm_win
                    b_nat = .true.
                 endwhere
                 found = .true.
              endif
           else
              write(LDT_logunit,*) '[WARN] Missing SMOPS file ', trim(fname)
           endif
        else
           write(LDT_logunit,'(a,i4.4,3i2.2)') &
                ' [WARN] No NPR-SMOPS-CMAP-6hrly file found for window starting ', &
                yr, mo, da, hr
        endif
        ! Next 6-hourly window
        call LDT_tick(wtime, doy, gmt, yr, mo, da, hr, mn, ss, 21600.0)
     enddo
     deallocate(sm_win, b_win)
  else
     call create_SMOPSsm_filename(SMOPSsmobs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
     inquire(file=trim(fname), exist=file_exists)
     if (file_exists) then
        write(LDT_logunit,*) '[INFO] Reading ', trim(fname)
        call read_SMOPS_native(n, fname, '', wtime, sm_nat, b_nat, found)
        write(LDT_logunit,*) '[INFO] Finished reading ', trim(fname)
     endif
  endif

  if (found) then
     call SMOPS_native_to_ldt(n, sm_nat, b_nat, smobs)
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if (smobs(c+(r-1)*LDT_rc%lnc(n)) .ne. -9999.0) then
              SMOPSsmobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  deallocate(sm_nat, b_nat)

  call LDT_logSingleDAobs(n, LDT_DAobsData(n)%soilmoist_obs, &
       SMOPSsmobs(n)%smobs, vlevel=1)

end subroutine readSMOPSsmObs


!BOP
!
! !ROUTINE: read_SMOPS_native
! \label{read_SMOPS_native}
!
! !INTERFACE:
subroutine read_SMOPS_native(n, fname, fversion, datatime, sm_out, b_out, ok)
!
! !USES:
#if (defined USE_GRIBAPI)
   use grib_api
#endif
   use LDT_coreMod,      only : LDT_rc
   use LDT_logMod,       only : LDT_logunit, LDT_endrun
   use SMOPSsm_obsMod,   only : SMOPSsmobs

   implicit none
!
! !ARGUMENTS:
!
   integer,           intent(in)  :: n
   character (len=*), intent(in)  :: fname
   character (len=*), intent(in)  :: fversion
   real*8,            intent(in)  :: datatime
   real,              intent(out) :: sm_out(SMOPSsmobs(n)%smopsnc*&
                                            SMOPSsmobs(n)%smopsnr)
   logical*1,         intent(out) :: b_out(SMOPSsmobs(n)%smopsnc*&
                                           SMOPSsmobs(n)%smopsnr)
   logical,           intent(out) :: ok
!
! !DESCRIPTION:
!  This subroutine reads one SMOPS grib2 file (daily or 6-hourly) and
!  applies the data quality flags. The result stays on the SMOPS 0.25 deg
!  grid, rows south to north.
!
!  The grid of every layer read is checked against the SMOPS grid LDT
!  uses: 1440 x 720 cells of 0.25 deg, row by row, the first column at
!  -179.875 (180.125 in GRIB2). Rows stored north to south are flipped;
!  rows stored south to north are kept. Any other grid stops LDT.
!
!  The GRIB2 parameter numbers depend on the SMOPS version, which is
!  taken from the file name for the NPR-SMOPS-CMAP-6hrly files (fversion,
!  e.g. 'v4r0'), and otherwise from "SMOPS soil moisture version:" or,
!  for "date-based", from the data time:
!
!  \begin{tabular}{lll}
!  version & ASCAT (QA)                       & dates \\
!  1.3     & A 213 (234), B 214 (235)         & before 2016-10-31T12 \\
!  2.0     & A 213 (234), B 214 (235)         & to 2017-08-24T12 \\
!  3.0     & A 213 (243), B 214 (244)         & to 2024-07-04T18 \\
!  4.0     & B 214 (244), C 250 (253)         & from 2024-07-04T18 \\
!  5.0     & not known yet; LDT stops         & from 2026-12-08 \\
!  \end{tabular}
!
!  LDT stops if a file's name gives a version other than v4r0, if it
!  differs from a version set in "SMOPS soil moisture version:", or if a
!  file without a version in its name is dated from 2026-12-08 and the
!  version is "date-based".
!
!  For ASCAT, the satellites are merged in the order listed; a later one
!  replaces an earlier one where it has a valid retrieval (as in LIS).
!
!  Quality flags are defined in:
!      NOAA NESDIS
!      CENTER FOR SATELLITE APPLICATIONS AND RESEARCH
!
!      SOIL MOISTURE OPERATIONAL PRODUCT SYSTEM (SMOPS)
!
!      ALGORITHM THEORETICAL BASIS DOCUMENT
!      Version 4.0
!
!  Found at http://www.ospo.noaa.gov/Products/land/smops/documents.html
!
!  The SMOPS QA flags are 16-bit (2-byte) integers, with the least
!  significant byte referred to as byte1 and the most significant byte
!  referred to as byte2.
!
!  bits: 15 | 14 | 13 | 12 | 11 | 10 | 9 | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0
!                   byte2                    |           byte1
!
!  The arguments are:
!  \begin{description}
!  \item[n]         index of the nest
!  \item[fname]     name of the SMOPS file
!  \item[fversion]  version from the file name ('v4r0', ...), or blank
!  \item[datatime]  LDT time of the data (selects the version if blank)
!  \item[sm\_out]   soil moisture on the SMOPS grid (udef where invalid)
!  \item[b\_out]    .true. where sm\_out is valid
!  \item[ok]        .false. if the file could not be read
! \end{description}
!
!EOP
   integer*1, parameter :: err_threshold = 5 ! in percent
   integer*2, parameter :: qa_undef = 9999
   ! EMK 21 Jul 2021...ISO Fortran restricts use of binary integer constants
   ! to DATA statements and as a principal argument to an intrinsic function.
   ! And gfortran 10 and 11 rejects their use in parameters.
   ! So we replace with arrays here.
   integer*1 :: AMSR2_accept(1)
   integer*2 :: SMOS_accept(8)
   data AMSR2_accept / b'00000001' /
   data SMOS_accept / b'0000000000000000', &
        b'0000000000000001', &
        b'0000000000001000', &
        b'0000000000001001', &
        b'0000000010000000', &
        b'0000000010000001', &
        b'0000000010001000', &
        b'0000000010001001' /

   ! GRIB2 parameter numbers; -1 = not in this SMOPS version
   integer   :: param_ASCAT(3), param_ASCAT_qa(3)
   integer   :: param_other, param_other_qa
   character(len=5) :: other_name

   real,      allocatable :: vals(:)
   real,      allocatable :: sm_ascat(:,:)
   integer*2, allocatable :: qa_ascat(:,:)
   real,      allocatable :: sm_other(:)
   integer*2, allocatable :: qa_other(:)

   integer        :: k, i, nnat
   integer        :: ftn, iret, igrib, param_num
   logical        :: north_first
   integer*2      :: qavalue
   integer*1      :: err
   character(len=3) :: vname

   ok = .false.
   sm_out = LDT_rc%udef
   b_out = .false.

#if (defined USE_GRIBAPI)
   nnat = SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr

   ! SMOPS version: from the file name (NPR-SMOPS-CMAP-6hrly files), else
   ! from "SMOPS soil moisture version:", else from the data time.
   if ( len_trim(fversion) > 0 ) then
      select case ( trim(fversion) )
      case ( 'v4r0' )
         vname = '4.0'
      case default
         ! v5r0 adds AMSR3 and may move the ASCAT layers; the numbers must
         ! come from a v5 file before v5 is read.
         write(LDT_logunit,*) '[ERR] ', trim(fname)
         write(LDT_logunit,*) '[ERR] is SMOPS ', trim(fversion), &
              '. LDT has the GRIB2 parameter numbers of SMOPS 1.3, 2.0, 3.0'
         write(LDT_logunit,*) '[ERR] and 4.0 (v4r0) only. Add those of ', &
              trim(fversion), ' to read_SMOPS_native.'
         call LDT_endrun()
      end select
      if ( SMOPSsmobs(n)%version /= 'date-based' .and. &
           SMOPSsmobs(n)%version /= 'NESDIS V3.0 REGEN' .and. &
           SMOPSsmobs(n)%version /= vname ) then
         write(LDT_logunit,*) '[ERR] ', trim(fname)
         write(LDT_logunit,*) '[ERR] is SMOPS ', vname, &
              ', but "SMOPS soil moisture version:" is ', &
              trim(SMOPSsmobs(n)%version)
         call LDT_endrun()
      endif
   elseif ( SMOPSsmobs(n)%version /= 'date-based' .and. &
            SMOPSsmobs(n)%version /= 'NESDIS V3.0 REGEN' ) then
      vname = SMOPSsmobs(n)%version(1:3)   ! checked in SMOPSsm_obsMod
   elseif ( datatime < SMOPSsmobs(n)%version2_time ) then
      vname = '1.3'
   elseif ( datatime < SMOPSsmobs(n)%version3_time ) then
      vname = '2.0'
   elseif ( datatime < SMOPSsmobs(n)%version4_time ) then
      vname = '3.0'
   elseif ( datatime < SMOPSsmobs(n)%version5_time ) then
      vname = '4.0'
   else
      write(LDT_logunit,*) '[ERR] ', trim(fname)
      write(LDT_logunit,*) '[ERR] has no version in its name and is dated ', &
           'from 2026-12-08, when SMOPS v5 starts.'
      write(LDT_logunit,*) '[ERR] LDT does not have the GRIB2 parameter ', &
           'numbers of SMOPS v5. If the file is v4,'
      write(LDT_logunit,*) '[ERR] set "SMOPS soil moisture version: 4.0".'
      call LDT_endrun()
   endif

   param_ASCAT = -1;  param_ASCAT_qa = -1
   param_other = -1;  param_other_qa = -1
   other_name = ''

   select case ( vname )
   case ( '1.3', '2.0' )
      param_ASCAT(1) = 213; param_ASCAT_qa(1) = 234
      param_ASCAT(2) = 214; param_ASCAT_qa(2) = 235
      if ( SMOPSsmobs(n)%useSMOS.eq.1 )  other_name = 'SMOS'
      if ( SMOPSsmobs(n)%useAMSR2.eq.1 ) other_name = 'AMSR2'
      if ( SMOPSsmobs(n)%useSMAP.eq.1 )  other_name = 'SMAP'
   case ( '3.0' )
      param_ASCAT(1) = 213; param_ASCAT_qa(1) = 243
      param_ASCAT(2) = 214; param_ASCAT_qa(2) = 244
      if ( SMOPSsmobs(n)%useSMOS.eq.1 ) then
         param_other = 212; param_other_qa = 242
      elseif ( SMOPSsmobs(n)%useAMSR2.eq.1 ) then
         param_other = 215; param_other_qa = 245
      elseif ( SMOPSsmobs(n)%useSMAP.eq.1 ) then
         param_other = 218; param_other_qa = 248
      endif
   case ( '4.0' )
      ! No ASCAT-A, no SMOS
      param_ASCAT(2) = 214; param_ASCAT_qa(2) = 244
      param_ASCAT(3) = 250; param_ASCAT_qa(3) = 253
      if ( SMOPSsmobs(n)%useSMOS.eq.1 ) then
         other_name = 'SMOS'
      elseif ( SMOPSsmobs(n)%useAMSR2.eq.1 ) then
         param_other = 215; param_other_qa = 245
      elseif ( SMOPSsmobs(n)%useSMAP.eq.1 ) then
         param_other = 218; param_other_qa = 248
      endif
   end select

   write(LDT_logunit,*) '[MSG] Reading SMOPS dataset as SMOPS version '// &
        vname

   if ( len_trim(other_name) > 0 ) then
      write(LDT_logunit,*) '[WARN] LDT does not process '//trim(other_name)// &
           ' from SMOPS version '//vname//'; file skipped.'
      return
   endif
   if ( SMOPSsmobs(n)%useASCAT.ne.1 ) param_ASCAT = -1

   allocate(vals(nnat))
   allocate(sm_ascat(nnat,3), qa_ascat(nnat,3))
   allocate(sm_other(nnat), qa_other(nnat))
   ! Fields missing from a file (e.g., ASCAT-A after Metop-A ended in
   ! Nov 2021) stay undefined.
   sm_ascat = LDT_rc%udef
   qa_ascat = qa_undef
   sm_other = LDT_rc%udef
   qa_other = qa_undef

   call grib_open_file(ftn, trim(fname), 'r', iret)
   if (iret .ne. 0) then
      write(LDT_logunit,*) '[WARN] Could not open file: ', trim(fname)
      deallocate(vals, sm_ascat, qa_ascat, sm_other, qa_other)
      return
   endif
   call grib_multi_support_on

   do
      call grib_new_from_file(ftn, igrib, iret)
      if ( iret == GRIB_END_OF_FILE ) exit
      if ( iret /= 0 ) then
         write(LDT_logunit,*) '[WARN] Problem reading GRIB message from ', &
              trim(fname)
         call grib_close_file(ftn)
         deallocate(vals, sm_ascat, qa_ascat, sm_other, qa_other)
         return
      endif

      call grib_get(igrib, 'parameterNumber', param_num, iret)
      if (iret /= 0) then
         write(LDT_logunit,*) '[WARN] Problem reading parameterNumber from ', &
              trim(fname)
         call grib_release(igrib, iret)
         call grib_close_file(ftn)
         deallocate(vals, sm_ascat, qa_ascat, sm_other, qa_other)
         return
      endif

      if ( param_num > 0 .and. ( any(param_ASCAT == param_num) .or. &
           any(param_ASCAT_qa == param_num) .or. &
           param_num == param_other .or. param_num == param_other_qa ) ) then

         call check_grid(igrib, param_num, north_first)
         call grib_get(igrib, 'values', vals, iret)
         if (iret /= 0) then
            write(LDT_logunit,*) '[WARN] Problem reading values of parameter ', &
                 param_num, ' from ', trim(fname)
            call grib_release(igrib, iret)
            call grib_close_file(ftn)
            deallocate(vals, sm_ascat, qa_ascat, sm_other, qa_other)
            return
         endif
         if ( north_first ) call flip_rows(vals)

         do k = 1, 3
            if ( param_num == param_ASCAT(k) )    sm_ascat(:,k) = vals
            if ( param_num == param_ASCAT_qa(k) ) qa_ascat(:,k) = int(vals, 2)
         enddo
         if ( param_num == param_other )    sm_other = vals
         if ( param_num == param_other_qa ) qa_other = int(vals, 2)
      endif

      call grib_release(igrib, iret)
   enddo

   call grib_close_file(ftn)

   if ( SMOPSsmobs(n)%useASCAT.eq.1 ) then
      ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
      ! (d) ASCAT Soil Moisture Product QA
      !
      ! Byte |  Description
      ! -----------------------------------------------------------------------
      ! 0    |  Estimated Error in Soil Moisture. (Integer. Scale factor: 0.01)
      ! 1    |  Soil Moisture Quality (Integer, Scale factor: 0.01)
      !
      ! The retrievals are rejected when the estimated error is above
      ! a predefined threshold (the recommeded value is 5%).
      !
      ! Technically speaking, err_threshold should be 0.05 and
      ! err should be scaled by 0.01.  But below is consistent with the NESDIS
      ! documentation.
      !
      ! Note that I am assuming that Byte 0 above refers to the least
      ! significant byte and Byte 1 above refers to the most significant byte.
      ! Meaning estimated error is get_byte1, and quality is get_byte2.
      do k = 1, 3
         if ( param_ASCAT(k) < 0 ) cycle
         do i = 1, nnat
            qavalue = qa_ascat(i,k)
            if ( qavalue == qa_undef ) then
               sm_ascat(i,k) = LDT_rc%udef
            else
               err = get_byte1(qavalue)
               if ( err >= err_threshold ) sm_ascat(i,k) = LDT_rc%udef
            endif
            if ( sm_ascat(i,k) < 0.001 ) sm_ascat(i,k) = LDT_rc%udef
         enddo
         ! A (1), then B (2), then C (3): later satellite wins
         where ( sm_ascat(:,k) /= LDT_rc%udef ) sm_out = sm_ascat(:,k)
      enddo

   elseif ( SMOPSsmobs(n)%useSMOS.eq.1 ) then
      ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
      ! (c) SMOS Soil Moisture Product QA
      !
      ! Byte 1:
      !
      ! Bit |  Description
      ! --------------------------------------------------------
      ! 0   |  Spare bit
      ! 1   |  1 = RFI for H pol above threshold, 0 = otherwise
      ! 2   |  1 = RFI for V pol above threshold, 0 = otherwise
      ! 3   |  Spare bit
      ! 4   |  1 = No products are generated, 0 = otherwise
      ! 5   |  1 = Retrieval values outside range, 0 = otherwise
      ! 6   |  1 = High retrieval DQX, 0 = otherwise
      ! 7   |  1 = Poor fit quality, 0 = otherwise
      !
      ! Byte 2:
      !
      ! Bit |  Description
      ! -------------------------------------------------------------
      ! 0   |  1 = Presence of other than nominal soil; 0 = otherwise
      ! 1   |  1 = Rocks; 0 = not rocks
      ! 2   |  1 = Moderate or strong topography; 0 = otherwise
      ! 3   |  1 = Open water; 0 = not open water
      ! 4   |  1 = Snow; 0 = not snow
      ! 5   |  1 = Forest; 0 = not forest
      ! 6   |  1 = Flood risk; 0 = no flood risk
      ! 7   |  1 = Urban area; 0 = not urban area
      !
      ! Accept when bytes 1 and 2 match one of SMOS_accept; otherwise reject.
      do i = 1, nnat
         qavalue = qa_other(i)
         if ( qavalue /= qa_undef .and. any(SMOS_accept == qavalue) ) then
            sm_out(i) = sm_other(i)
         endif
      enddo

   elseif ( SMOPSsmobs(n)%useAMSR2.eq.1 ) then
      ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
      ! (e) AMSR2 Soil Moisture Layer QA
      !
      ! Byte 1:
      !
      ! Bit |  Description
      ! -------------------------------------------------------------------
      ! 0   |  0 = overall quality is not good; 1 = overall quality is good
      ! 1   |  1 = retrieval attempted but quality is not good; 0 = otherwise
      ! 2   |  1 = retrieval attempted but unsuccessful due to input
      !     |      brightness temperature data quality; 0 = otherwise
      ! 3   |  1 = retrieval attempted but unsuccessful due to the quality
      !     |      of other input data; 0 = otherwise
      ! 4   |  1 = retrieval not attempted; 0 = retrieval attempted
      ! 5   |  0 = not cold desert; 1 = cold desert
      ! 6   |  0 = not snow or rain; 1 = snow or rain
      ! 7   |  0 = not frozen ground; 1 = frozen ground
      !
      ! Byte 2 (GVF classes, TB quality, NDVI source) is NOT considered.
      !
      ! Accept when byte 1 is b'00000001'; otherwise reject.
      do i = 1, nnat
         qavalue = qa_other(i)
         if ( qavalue /= qa_undef ) then
            if ( get_byte1(qavalue) == AMSR2_accept(1) ) sm_out(i) = sm_other(i)
         endif
      enddo

   elseif ( SMOPSsmobs(n)%useSMAP.eq.1 ) then
      ! SMAP: not yet documented.  Applying ASCAT QA logic.
      do i = 1, nnat
         qavalue = qa_other(i)
         if ( qavalue /= qa_undef ) then
            if ( get_byte1(qavalue) < err_threshold ) sm_out(i) = sm_other(i)
         endif
      enddo
   endif

   b_out = ( sm_out /= LDT_rc%udef )
   ok = .true.

   deallocate(vals, sm_ascat, qa_ascat, sm_other, qa_other)
#endif

contains

#if (defined USE_GRIBAPI)
subroutine check_grid(igrib, param_num, north_first)
   ! Stops LDT unless the message is on the SMOPS grid LDT uses (see
   ! above); north_first tells whether its rows run north to south.
   implicit none
   integer, intent(in)  :: igrib, param_num
   logical, intent(out) :: north_first
   integer :: ni, nj, jpos, ineg, jcons, st(7)
   real*8  :: lat1, lon1
   real*8, parameter :: tol = 1d-3

   call grib_get(igrib, 'Ni', ni, st(1))
   call grib_get(igrib, 'Nj', nj, st(2))
   call grib_get(igrib, 'jScansPositively', jpos, st(3))
   call grib_get(igrib, 'iScansNegatively', ineg, st(4))
   call grib_get(igrib, 'jPointsAreConsecutive', jcons, st(5))
   call grib_get(igrib, 'latitudeOfFirstGridPointInDegrees', lat1, st(6))
   call grib_get(igrib, 'longitudeOfFirstGridPointInDegrees', lon1, st(7))
   north_first = ( jpos == 0 )

   if ( any(st /= 0) ) then
      write(LDT_logunit,*) '[ERR] Cannot read the grid of parameter ', &
           param_num, ' in ', trim(fname)
      call LDT_endrun()
   endif
   if ( ni /= SMOPSsmobs(n)%smopsnc .or. nj /= SMOPSsmobs(n)%smopsnr .or. &
        ineg /= 0 .or. jcons /= 0 .or. &
        abs(modulo(lon1, 360d0) - 180.125d0) > tol .or. &
        abs(lat1 - merge(89.875d0, -89.875d0, north_first)) > tol ) then
      write(LDT_logunit,*) '[ERR] Parameter ', param_num, ' in ', trim(fname)
      write(LDT_logunit,*) '[ERR] is not on the SMOPS grid LDT expects.'
      write(LDT_logunit,'(a,2i6,a,3i3,a,2f10.4)') ' [ERR] Ni Nj', ni, nj, &
           '  jScansPositively iScansNegatively jPointsAreConsecutive', &
           jpos, ineg, jcons, '  first lat lon', lat1, lon1
      write(LDT_logunit,'(a,2i6,a)') ' [ERR] Expected: Ni Nj', &
           SMOPSsmobs(n)%smopsnc, SMOPSsmobs(n)%smopsnr, &
           ', iScansNegatively 0, jPointsAreConsecutive 0, first lon'// &
           ' 180.125 (= -179.875), first lat 89.875 (N to S) or -89.875'
      call LDT_endrun()
   endif
end subroutine check_grid
#endif

subroutine flip_rows(a)
   ! Rows north to south -> south to north.
   implicit none
   real, intent(inout) :: a(:)
   integer :: r, nc, nr
   real    :: row(SMOPSsmobs(n)%smopsnc)
   nc = SMOPSsmobs(n)%smopsnc
   nr = SMOPSsmobs(n)%smopsnr
   do r = 1, nr/2
      row = a((r-1)*nc+1:r*nc)
      a((r-1)*nc+1:r*nc) = a((nr-r)*nc+1:(nr-r+1)*nc)
      a((nr-r)*nc+1:(nr-r+1)*nc) = row
   enddo
end subroutine flip_rows

integer*1 function get_byte1(i)
   implicit none
   integer*2, intent(in) :: i
   integer*2 :: j
   ! This function expects a 16-bit integer as input.  It returns
   ! the least significant byte (as an 8-bit integer), referred
   ! to as byte1 in the NESDIS documention.
   !
   ! For example,
   ! i = b'0000001000000001' <--> 00000010|00000001
   !                         <--> byte2|byte1
   !                         <--> 0x0201
   ! Here byte1 is b'00000001'; byte2 is b'00000010'
   j = iand(i, int(z'00ff', 2))
   get_byte1 = int(j, 1)
end function get_byte1

end subroutine read_SMOPS_native


!BOP
!
! !ROUTINE: SMOPS_native_to_ldt
! \label{SMOPS_native_to_ldt}
!
! !INTERFACE:
subroutine SMOPS_native_to_ldt(n, sm_nat, b_nat, smobs_ip)
!
! !USES:
   use LDT_coreMod,      only : LDT_rc, LDT_domain
   use SMOPSsm_obsMod,   only : SMOPSsmobs
   ! MN : added to read the greenness, land cover and land mask
   use LDT_gfracMod,     only : LDT_gfrac_struc
   use LDT_LMLCMod
   use LDT_paramDataMod, only : LDT_LSMparam_struc

   implicit none
!
! !ARGUMENTS:
   integer,   intent(in)  :: n
   real,      intent(in)  :: sm_nat(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   logical*1, intent(in)  :: b_nat(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real,      intent(out) :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  Interpolates the SMOPS field to the LDT grid (nearest neighbour) and
!  screens it: for ASCAT, barren or sparsely vegetated cells (and, if
!  greenness is read, cells with maximum greenness <= 0.1); for all
!  sensors, cells within "SMOPS search radius for openwater proximity
!  detection" of a water cell.
!
!EOP
   real           :: sm_in(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   logical*1      :: b_in(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   logical*1      :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   integer        :: c, r, iret
   integer        :: search_rad, ix, jx, c_s, c_e, r_s, r_e

   sm_in = sm_nat
   b_in  = b_nat

   !--------------------------------------------------------------------------
   ! Interpolate to the LDT running domain
   !--------------------------------------------------------------------------
   call neighbor_interp(LDT_rc%gridDesc(n,:),&
        b_in, sm_in, smobs_b_ip, smobs_ip, &
        SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr, &
        LDT_rc%lnc(n)*LDT_rc%lnr(n), &
        LDT_domain(n)%lat, LDT_domain(n)%lon,&
        SMOPSsmobs(n)%n11, LDT_rc%udef, iret)

   if ( SMOPSsmobs(n)%useASCAT.eq.1 ) then
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            ! MN : add condition to check the GVF and land cover
            !      screen the ASCAT if GVF is less than 0.1 and
            !      land cover type is Barren or Sparsely Vegetated and
            !      the value of the pixel is larger than 0.5
            if ( allocated(LDT_gfrac_struc(n)%gfrac%value) ) then
               if ( maxval(LDT_gfrac_struc(n)%gfrac%value(c,r,:)) &
                    .le. 0.1 .and. &
                    LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%bareclass) &
                    .ge. 0.5  ) then
                  smobs_ip(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
               endif
            else
               if ( LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%bareclass) &
                    .ge. 0.5  ) then
                  smobs_ip(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
               endif
            endif
         enddo
      enddo
   endif

   !------------------------------------------------------------------------
   !  Remove pixel close to open water
   !------------------------------------------------------------------------
   search_rad = nint(SMOPSsmobs(n)%search_radius)
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         if (smobs_ip(c+(r-1)*LDT_rc%lnc(n)) .ne. LDT_rc%udef) then
            c_s = max(1,c-search_rad)
            c_e = min(LDT_rc%lnc(n),c+search_rad)
            r_s = max(1,r-search_rad)
            r_e = min(LDT_rc%lnr(n),r+search_rad)
            do ix=c_s,c_e
               do jx=r_s,r_e
                  if( LDT_LSMparam_struc(n)%landmask%value(ix,jx,1) .eq. 0 ) then
                     smobs_ip(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
                  end if
               enddo
            enddo
         endif
      enddo
   enddo

end subroutine SMOPS_native_to_ldt


!BOP
! !ROUTINE: create_SMOPSsm_filename
! \label{create_SMOPSsm_filename}
!
! !INTERFACE:
subroutine create_SMOPSsm_filename(ndir, yr, mo,da, filename)
! !USES:

  implicit none
! !ARGUMENTS:
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
!
! !DESCRIPTION:
!  This subroutine creates the daily SMOPS filename based on the date
!
!  The arguments are:
!  \begin{description}
!  \item[ndir] name of the SMOPS soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated RT SMOPS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
       //trim(fyr)//trim(fmo)//trim(fda)//'.gr2'

end subroutine create_SMOPSsm_filename


!BOP
! !ROUTINE: create_SMOPSsm_6hrly_filename
! \label{create_SMOPSsm_6hrly_filename}
!
! !INTERFACE:
subroutine create_SMOPSsm_6hrly_filename(n, ndir, conv, yr, mo, da, hr, &
     filename, fversion)
! !USES:
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
                               LDT_releaseUnitNumber, LDT_endrun
  use LDT_timeMgrMod,   only : LDT_date2time, LDT_tick
  use SMOPSsm_obsMod,   only : SMOPSsmobs

  implicit none
! !ARGUMENTS:
  integer,           intent(in)  :: n
  character (len=*), intent(in)  :: ndir
  character (len=*), intent(in)  :: conv
  integer,           intent(in)  :: yr, mo, da, hr
  character (len=*), intent(out) :: filename
  character (len=*), intent(out) :: fversion
!
! !DESCRIPTION:
!  Creates the name of the 6-hourly SMOPS file that holds the window
!  starting at yr-mo-da hr, following create\_SMOPS\_ASCATsm\_filename
!  in the LIS SMOPS ASCAT reader.
!
!  "AGRMET" naming convention (files directly in ndir):
!  \begin{description}
!  \item[before 2017-10-05] smops\_dYYYYMMDD\_sHH0000\_cness.gr2,
!     HH = window start
!  \item[2017-10-05 to 2024-07-04T18] NPR\_SMOPS\_CMAP\_DYYYYMMDDHH.gr2,
!     HH = window start + 6 hours (the file holds the 6 hours before HH)
!  \item[from 2024-07-04T18]
!     NPR-SMOPS-CMAP-6hrly\_vNrM\_blend\_sYYYYMMDDHH*\_e*\_c*.grib2
!     (v4r0 from 2024-07-04T18, v5r0 from 2026-12-08), HH = window start;
!     the version, end and creation times are found with glob. If several
!     files match, the newest version is used and, within it, the latest
!     creation time.
!  \end{description}
!
!  "LIS" naming convention: ndir/YYYY/NPR\_SMOPS\_CMAP\_DYYYYMMDDHH.gr2
!  for all dates, with the same 6-hour offset.
!
!  filename is 'NULL' if no NPR-SMOPS-CMAP-6hrly file is found. fversion
!  is the version in the name of an NPR-SMOPS-CMAP-6hrly file (e.g.
!  'v4r0'), and blank for the other names.
!EOP
  integer, external :: create_filelist ! C function (ldt/DAobs/WSF_sm)
  character (len=4) :: fyr
  character (len=2) :: fmo, fda, fhr
  character (len=*), parameter :: listfile = 'SMOPS_filelist_ldt.dat'
  character (len=len(filename)) :: line
  character (len=16) :: ver
  real*8            :: t, tfile
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy, rc, ftn, ierr
  integer           :: vnum, best
  real              :: gmt

  fversion = ''
  yr2 = yr; mo2 = mo; da2 = da; hr2 = hr; mn2 = 0; ss2 = 0
  call LDT_date2time(t, doy, gmt, yr2, mo2, da2, hr2, mn2, ss2)

  ! t (window start) picks the naming convention; yr2..hr2 go in the name.
  if ( t >= SMOPSsmobs(n)%naming3_time .and. &
       t <  SMOPSsmobs(n)%version4_time ) then
     ! NPR_SMOPS_CMAP_DYYYYMMDDHH.gr2 holds the 6 hours before HH
     call LDT_tick(tfile, doy, gmt, yr2, mo2, da2, hr2, mn2, ss2, 21600.0)
  endif

  write(unit=fyr, fmt='(i4.4)') yr2
  write(unit=fmo, fmt='(i2.2)') mo2
  write(unit=fda, fmt='(i2.2)') da2
  write(unit=fhr, fmt='(i2.2)') hr2

  if ( trim(conv) == 'LIS' ) then
     filename = trim(ndir)//'/'//fyr//'/NPR_SMOPS_CMAP_D' &
          //fyr//fmo//fda//fhr//'.gr2'

  elseif ( t < SMOPSsmobs(n)%naming3_time ) then
     filename = trim(ndir)//'/smops_d'//fyr//fmo//fda//'_s'//fhr// &
          '0000_cness.gr2'

  elseif ( t < SMOPSsmobs(n)%version4_time ) then
     filename = trim(ndir)//'/NPR_SMOPS_CMAP_D'//fyr//fmo//fda//fhr//'.gr2'

  else
     filename = 'NULL'
     rc = create_filelist(trim(ndir)//'/NPR-SMOPS-CMAP-6hrly_v*_blend_s' &
          //fyr//fmo//fda//fhr//'*.grib2'//char(0), listfile//char(0))
     if (rc /= 0) then
        write(LDT_logunit,*) '[WARN] Problem searching for ', &
             'NPR-SMOPS-CMAP-6hrly files in ', trim(ndir)
     endif
     best = -1
     ftn = LDT_getNextUnitNumber()
     open(ftn, file=listfile, status='old', iostat=ierr)
     do while (ierr == 0)
        read(ftn, '(a)', iostat=ierr) line
        if (ierr /= 0) exit
        if (len_trim(line) == 0) cycle
        call SMOPS_6hrly_version(line, ver, vnum)
        if (vnum < 0) then
           write(LDT_logunit,*) '[ERR] No SMOPS version (vNrM) in the name ', &
                trim(line)
           call LDT_endrun()
        endif
        ! glob lists the files by name, so within one version the last
        ! file has the latest creation time.
        if (vnum >= best) then
           best = vnum
           filename = line
           fversion = ver
        endif
     enddo
     close(ftn)
     call LDT_releaseUnitNumber(ftn)
  endif

end subroutine create_SMOPSsm_6hrly_filename

!BOP
! !ROUTINE: SMOPS_6hrly_version
! \label{SMOPS_6hrly_version}
!
! !INTERFACE:
subroutine SMOPS_6hrly_version(fname, ver, vnum)
  implicit none
! !ARGUMENTS:
  character (len=*), intent(in)  :: fname
  character (len=*), intent(out) :: ver
  integer,           intent(out) :: vnum
!
! !DESCRIPTION:
!  Reads the SMOPS version from the name of an NPR-SMOPS-CMAP-6hrly file,
!  e.g. NPR-SMOPS-CMAP-6hrly\_v5r0\_blend\_s...grib2 gives ver = 'v5r0'
!  and vnum = 5000 (1000 * N + M for vNrM, for comparing versions).
!  vnum is -1 and ver blank if the name has no vNrM.
!EOP
  character (len=*), parameter :: tag = '6hrly_v'
  character (len=*), parameter :: digits = '0123456789'
  integer :: i, j, k, major, minor, ios

  ver = ''
  vnum = -1
  i = index(fname, tag, back=.true.)
  if (i == 0) return
  i = i + len(tag) - 1                    ! the 'v'
  j = index(fname(i:), '_blend')
  if (j < 5) return                       ! at least 'vNrM'
  k = index(fname(i:i+j-2), 'r')          ! the 'r', relative to i
  if (k < 3 .or. k > j - 2) return
  if (verify(fname(i+1:i+k-2), digits) /= 0) return
  if (verify(fname(i+k:i+j-2), digits) /= 0) return
  read(fname(i+1:i+k-2), *, iostat=ios) major
  if (ios /= 0) return
  read(fname(i+k:i+j-2), *, iostat=ios) minor
  if (ios /= 0) return
  ver = fname(i:i+j-2)
  vnum = 1000*major + minor

end subroutine SMOPS_6hrly_version
