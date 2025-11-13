!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

! Main SNIP driver
subroutine SNIP_run(n)

  !************************************************************************
  !************************************************************************
  !**
  !**  NAME: SNOW DEPTH ANALYSIS MODEL
  !**
  !**  PURPOSE:  DRIVER ROUTINE FOR THE SNOW DEPTH ANALYSIS MODEL.
  !**
  !**   See USAFSI_arraysMod.F90 for common array descriptions
  !**   See USAFSI_paramsMod.F90 for common parameter descriptions
  !**
  !**  UPDATES
  !**  =======
  !**  26 APR 95  INITIAL UNISYS VERSION....................SSGT CONRY/SYSM
  !**  09 FEB 96  CHANGED ARGUMENTS PASSED BETWEEN GETOBS, GETSNO, GETSMI,
  !**             AND PROCES DUE TO CHANGES DESCRIBED IN THEIR PROLOGUES.
  !**             FORMATTED PRINT STATEMENTS.......DR KOPP, SSGT CONRY/SYSM
  !**  23 OCT 97  ADDED INCLUDE OF TUNES PROC AND CALL TO RDTUNE...........
  !**             ..........................................SSGT CONRY/DNXM
  !**  23 FEB 01  PORTED FROM UNISYS MAINFRAME.............SSGT MILLER/DNXM
  !**  30 JAN 03  ADDED ROUTINE BLACKLIST.....................MR GAYNO/DNXM
  !**  09 JUL 03  REMOVED CALL TO GETSST; NOW PART OF NT2ICE..MR GAYNO/DNXM
  !**  05 FEB 04  CHANGED GEOG RECORD LENGTH FROM 4-BYTE TO 2-BYTE WORD
  !**             (INCORPORATED 18 APR 03 AER CHANGES).....MR LEWISTON/DNXM
  !**  05 MAY 04  MERGED IN THE SSMIS NASA TEAM II SEA ICE CONCENTRATION
  !**             CAPABILITY.  CONVERTED TO FREE FORMAT FORTRAN 90 USING
  !**             MODULES RATHER THAN COMMON BLOCKS........MR EYLANDER/DNXM
  !**  21 OCT 04  PORTED TO CDFS II.  REMOVED GRIB SUBROUTINES SINCE
  !**             CDFS II WRITING C++ GRIB ROUTINES........MR EYLANDER/DNXM
  !**  21 MAR 05  REMOVED OBSOLETE ARRAY SMITMP............MR LEWISTON/DNXM
  !**  14 APR 05  SEPARATED SNOW AND ICE INTO SEPARATE ARRAYS.
  !**             ADDED CALL TO GETICE.....................MR LEWISTON/DNXM
  !**  07 JUN 05  REMOVED JULHR FROM GETOBS & GETSMI CALLS.................
  !**             .........................................MR LEWISTON/DNXM
  !**  08 MAY 09  ADDED AMSR-E SNOW. REMOVED SEA ICE...MR LEWISTON/2WXG/WEA
  !**  16 JUN 09  CONVERTED TO EQUIDISTANT CYLINDRICAL GRID...............
  !**             .....................................MR LEWISTON/2WXG/WEA
  !**  18 Dec 09  ADDED GRIB OUTPUT....................MR LEWISTON/16WS/WXE
  !**  05 JAN 10  CHANGED SSMIS FROM SDRS TO EDRS......MR LEWISTON/16WS/WXE
  !**  13 SEP 10  ADDED CALL TO GETSFC.................MR LEWISTON/16WS/WXE
  !**  29 APR 11  MOVED BLACKLIST AND VALIDATION TO DBPULL................
  !**             .....................................MR LEWISTON/16WS/WXE
  !**  08 NOV 11  REMOVED AMSR-E AND PORTED TO LINUX...MR LEWISTON/16WS/WXE
  !**  20 MAR 12  SWITCHED OBS SOURCE FROM CDMS TO JMOBS..................
  !**             .....................................MR LEWISTON/16WS/WXE
  !**  29 MAY 12  MOVED MODIFIED TEST FROM SCRIPT TO GETSNO...............
  !**             .....................................MR LEWISTON/16WS/WXE
  !**  07 NOV 12  ADDED FRACTIONAL SNOW................MR LEWISTON/16WS/WXE
  !**  10 OCT 13  CHANGED DAILY TO 4 TIMES/DAY. ADDED CYCLE_LOOP.
  !**             REPLACED NAMELIST WITH GETENV CALLS......................
  !**             .....................................MR LEWISTON/16WS/WXE
  !**  15 FEB 17  ADDED VIIRS DATA.......................MR PUSKAR/16WS/WXE
  !**  28 DEC 17  ADDED FILTER FOR MISSING ELEVATIONS.....................
  !**             .....................................MR LEWISTON/16WS/WXE
  !**  25 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
  !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
  !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
  !**  02 Nov 20  Removed blacklist logic...Eric Kemp, NASA GSFC/SSAI
  !**  06 Nov 20  Added GALWEM 2-m temperature....Eric Kemp, NASA GSFC/SSAI
  !**  12 Dec 20  Added GMI and AMSR2 capability..........................
  !**             ............................Yonghwan Kwon/NASA GSFC/ESSIC
  !**  26 Jan 21  Added revised 10-km snow climatology....................
  !**             ..............................Yeosang Yoon/NASA GSFC/SAIC
  !**  28 Jan 21  Updated messages for PMW snow retrievals and cleaned some
  !**             unused codes..................Yeosang Yoon/NASA GSFC/SAIC
  !**  13 Jan 22  Added support for FNMOC SST GRIB1 file..................
  !**             .................................Eric Kemp/NASA GSFC/SSAI
  !**  27 Jun 23  Removed LDT_endrun for normal termination, to avoid error
  !**             code 1.....................................Eric Kemp/SSAI
  !**  28 Jun 23  Extended station names to 31 characters....Eric Kemp/SSAI
  !**  24 Aug 23  Changed station names to 32 characters.....Eric Kemp/SSAI
  !**  19 Jul 24  Added ESPC-D support.......................Eric Kemp/SSAI
  !**  09 Jul 25  Migrated to SNIP...........................Eric Kemp/SSAI
  !**  11 Jul 25  Removed legacy satellite SD retrievals.  Added code
  !               to read SNIP or USAFSI prior analyses......Eric Kemp/SSAI
  !**  30 Jul 25  Add ASMR2 snow depth support...............Eric Kemp/SSAI
  !************************************************************************
  !************************************************************************

  ! Imports
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN
  use LDT_coreMod, only: LDT_masterproc, LDT_rc
  use LDT_logMod, only: LDT_logunit, LDT_endrun
  use LDT_pluginIndices
  use LDT_SNIPMod, only: SNIP_settings
  use map_utils
#if ( defined SPMD )
  use mpi
#endif
  use SNIP_analysisMod
  use SNIP_arraysMod, only: SNIP_arrays
  use SNIP_bratsethMod
  use SNIP_espcdMod
  use SNIP_galwemMod, only: SNIP_get_galwem_t2m
  use SNIP_gofsMod
  use SNIP_lisMod, only:  read_gr2_t2
  use SNIP_netcdfMod
  use SNIP_paramsMod
  use SNIP_utilMod

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: n

  ! Local variables
  character*10               ::  date10               ! DATE-TIME GROUP OF CYCLE
  character*90               ::  message    (msglns)  ! ERROR MESSAGE
  character*5,  allocatable  ::  netid      (:)       ! NETWORK ID OF AN OBSERVATION
  character(LDT_CONST_PATH_LEN) ::  snodep_modif         ! PATH TO MODIFIED DATA DIRECTORY
  character(LDT_CONST_PATH_LEN) ::  snodep_unmod         ! PATH TO UNMODIFIED DATA DIRECTORY
  character(LDT_CONST_PATH_LEN) ::  sfcobs               ! PATH TO DBPULL SNOW OBS DIRECTORY
  integer :: sfcobsfmt ! Format of sfcobs file
  character*32,  allocatable  ::  staid      (:)       ! STATION ID OF AN OBSERVATION
  character(LDT_CONST_PATH_LEN) ::  static               ! STATIC FILE DIRECTORY PATH
  character(LDT_CONST_PATH_LEN) ::  stmpdir              ! SFC TEMP DIRECTORY PATH
  character(LDT_CONST_PATH_LEN) :: sstdir ! EMK 20220113

  character(LDT_CONST_PATH_LEN) ::  viirsdir             ! PATH TO VIIRS DATA DIRECTORY
  integer                    ::  runcycle             ! CYCLE HOUR
  integer                    ::  hemi                 ! HEMISPHERE (1 = NH, 2 = SH, 3 = GLOBAL)
  integer                    ::  julhr                ! AFWA JULIAN HOUR BEING PROCESSED
  integer                    ::  julhr_beg            ! AFWA JULIAN HOUR OF BEGINNING CYCLE
  integer                    ::  julhr_end            ! AFWA JULIAN HOUR OF ENDING CYCLE
  integer                    ::  month                ! MONTH OF YEAR (1-12)
  integer                    ::  stacnt               ! TOTAL NUMBER OF OBSERVATIONS USED
  integer,    allocatable    ::  staelv     (:)       ! OBSERVATION ELEVATION (METERS)
  integer,    allocatable    ::  stalat     (:)       ! OBSERVATION LATITUDE
  integer,    allocatable    ::  stalon     (:)       ! OBSERVATION LONGITUDE
  logical                    ::  sfctmp_found         ! FLAG FOR SFC TEMP FILE FOUND
  real,       allocatable    ::  sfctmp (:, :)        ! GALWEM OR LIS SHELTER TEMPERATURE DATA
  real,       allocatable    ::  stadep     (:)       ! OBSERVATION SNOW DEPTH (METERS)
  character*20 :: routine_name
  type(SNIP_bratseth_t) :: bratseth
  character*10 :: network10
  character*32 :: platform32
  real :: rob, rlat, rlon, relev
  integer :: nc,nr
  real, allocatable :: landmask(:,:)
  real, allocatable :: elevations(:,:)
  real, allocatable :: landice(:,:)
  integer :: j
  integer :: c, r
  real :: arctlatr
  real, allocatable :: climo_tmp(:,:)
  integer :: maxsobs
  integer :: yyyy, mm, dd, hh, fh
  integer :: ierr
  logical :: found_navy_cice
  logical :: just_12z

  maxsobs = SNIP_settings%maxsobs

  ! Only the master process handles the file output
  if (LDT_masterproc) then
     routine_name = "SNIP"
     message = ' '
     hemi    = 3

     write (LDT_logunit,*) '[INFO] RUNNING SNIP MODEL'

     ! Check the LDT map projection.
     ! FIXME:  Support other projections in addition to LATLON
     if (trim(LDT_rc%lis_map_proj(1)) .ne. LDT_latlonId) then
        write(LDT_logunit,*)'[ERR] SNIP only supports lat/lon grid'
        call LDT_endrun()
     end if

     ! Get relevant fields from LDT parameter file
     call read_params(nc,nr,landmask,elevations, landice)

     ! Copy ldt.config files to local variables
     date10 = trim(SNIP_settings%date10)
     snodep_modif = trim(SNIP_settings%snodep_modif)
     snodep_unmod = trim(SNIP_settings%snodep_unmod)

     sfcobs = trim(SNIP_settings%sfcobs)
     sfcobsfmt = SNIP_settings%sfcobsfmt
     stmpdir = trim(SNIP_settings%stmpdir)
     sstdir = trim(SNIP_settings%sstdir)
     static = trim(SNIP_settings%static)
     viirsdir = trim(SNIP_settings%viirsdir)

     ! EXTRACT MONTH FROM DATE-TIME GROUP.
     read (date10(5:6), '(i2)', err=4200) month

     ! ALLOCATE GRIDDED DATA ARRAYS.
     allocate (SNIP_arrays%climo            (nc,     nr))
     allocate (SNIP_arrays%elevat           (nc,     nr))
     allocate (SNIP_arrays%iceage           (nc,     nr))
     allocate (SNIP_arrays%iceage12z        (nc,     nr))
     allocate (SNIP_arrays%icecon           (nc,     nr))
     allocate (SNIP_arrays%icemask          (nc,     nr))
     allocate (SNIP_arrays%oldcon           (nc,     nr))
     allocate (SNIP_arrays%olddep           (nc,     nr))
     allocate (SNIP_arrays%oldmask          (nc,     nr))
     allocate (SNIP_arrays%ptlat            (nc,     nr))
     allocate (SNIP_arrays%ptlon            (nc,     nr))
     ! This is the GALWEM or LIS data interpolated to the LDT grid
     allocate (sfctmp       (nc,     nr))
     allocate (SNIP_arrays%snoage           (nc,     nr))
     allocate (SNIP_arrays%snoage12z        (nc,     nr))
     allocate (SNIP_arrays%snoanl           (nc,     nr))
     allocate (SNIP_arrays%snofrac          (nc,     nr))
     allocate (SNIP_arrays%snow_poss        (nc,     nr))
     allocate (SNIP_arrays%sst              (nc,     nr))
     allocate (SNIP_arrays%viirsmap         (nc,     nr))
     allocate (SNIP_arrays%navy_icecon(nc,nr))
     allocate (SNIP_arrays%amsr2_snowdepth(nc,nr))

     ! RETRIEVE STATIC DATA SETS.
     write (LDT_logunit,*) '[INFO] CALLING GETGEO TO GET STATIC FIELDS'
     ! Pass LDT elevations to this routine to populate SNODEP elevat
     call getgeo (month, static, nc, nr, elevations)

     ! Retrieve the snow climatology for the month
     if (SNIP_settings%climo_option .eq. 1) then
        ! Mismatches can occur in landmask between 0.25 deg climo and
        ! LDT's grid.  In cases where LDT introduces "new" land, use a
        ! bogus value.
        arctlatr = float(arctlat) / 100.0
        allocate(climo_tmp(nc,nr))
        climo_tmp = SNIP_arrays%climo
        do r = 1, nr
           do c = 1, nc
              if (SNIP_arrays%ptlat(c,r) > -40.0 .and. &
                   SNIP_arrays%ptlat(c,r) < 40.0) cycle
              ! See if climo exists for glacier point.  If not, use a
              ! fill-in.
              if (landmask(c,r) > 0.5 .and. &
                   SNIP_arrays%climo(c,r) .le. 0) then
                 if (landice(c,r) > 0.5) then
                    climo_tmp(c,r) = SNIP_settings%fill_climo
                 else
                    climo_tmp(c,r) = 0
                 end if
              end if
           end do ! c
        end do ! r

        SNIP_arrays%climo(:,:) = climo_tmp(:,:)
        deallocate(climo_tmp)
     else if (SNIP_settings%climo_option .eq. 2) then
        ! newly produced 10-km snow climatology using ratio between
        ! SNODEP and USAFSI
        call getclimo (month, static)
     end if

     ! RETRIEVE THE PREVIOUS SNOW ANALYSIS.
     ! First, try reading SNIP in netCDF format.  If that doesn't work,
     ! fall back to reading USAFSI in netCDF format.
     ! fall back on the legacy SNODEP at 0.25 deg resolution.
     write (LDT_logunit,*) &
          '[INFO] CALLING GETSNO_SNIP_NC TO GET PREVIOUS SNIP DATA'
     call getsno_snip_nc(date10, julhr_beg, ierr)
     if (ierr .ne. 0) then
        write (LDT_logunit,*) &
             '[INFO] CALLING GETSNO_USAFSI_NC TO GET PREVIOUS USAFSI DATA'
        call getsno_usafsi_nc(date10, julhr_beg, ierr)
     end if
     if (ierr .ne. 0) then
        write (LDT_logunit,*) &
             '[INFO] CALLING GETSNO TO GET PREVIOUS SNODEP DATA'
        if (ierr == 1) then
           just_12z = .true.
        else
           just_12z = .false.
        end if
        call getsno (date10, snodep_modif, snodep_unmod, nc, nr, &
             landice, julhr_beg, just_12z)
     end if

     julhr = julhr_beg
     call date10_julhr (date10, julhr_end, program_name, routine_name)

     ! LOOP THROUGH UNCOMPLETED CYCLES.
     ! FIXME...Change this to handle only one cycle
     cycle_loop : do while (julhr .lt. julhr_end)

        ! EMK Create bratseth object
        call bratseth%new(maxsobs, SNIP_settings%back_err_var, &
             SNIP_settings%back_err_h_corr_len, &
             SNIP_settings%back_err_v_corr_len)

        julhr = julhr + 6
        call julhr_date10 (julhr, date10, program_name, routine_name)

        write(LDT_logunit,*) &
             '[INFO] ***PREPARING SNOW/ICE ANALYSIS VALID ', date10, &
             '***'

        write (ldt_logunit,6600) date10
        read (date10(9:10), '(i2)', err=4200) runcycle

        ! RETRIEVE GALWEM SURFACE TEMPERATURES.
        write (LDT_logunit,*) &
       '[INFO] CALLING SNIP_get_galwem_t2m TO GET GALWEM 2-M TEMPERATURES'
        call SNIP_get_galwem_t2m(n, julhr, nc, nr, sfctmp, ierr)
        if (ierr .ne. 0) then
           ! Fall back on LIS SURFACE temperatures
           write (LDT_logunit,*) &
              '[INFO] CALLING READ_GR2_T2 TO GET LIS SHELTER TEMPERATURES'
           call read_gr2_t2(date10, nc, nr, sfctmp, ierr)
           if (ierr .ne. 0) then
              write (LDT_logunit,*) &
                   '[INFO] CALLING GETSFC TO GET LIS SHELTER TEMPERATURES'
              call getsfc (date10, stmpdir, sfctmp_found, sfctmp)
           else
              sfctmp_found = .true.
           end if
        else
           sfctmp_found = .true.
        end if

        ! Retrieve sea surface temperature data.
        read (date10(1: 4), '(i4)', err=4200) yyyy
        read (date10(5: 6), '(i2)', err=4200) mm
        read (date10(7: 8), '(i2)', err=4200) dd
        read (date10(9:10), '(i2)', err=4200) hh
        fh = 0 ! Dummy value
        if (SNIP_settings%source_of_ocean_data == "ESPC-D") then
           write (LDT_logunit,*) &
        '[INFO] CALLING PROCESS_ESPCD_SST TO GET SEA SURFACE TEMPERATURES'
           call process_espcd_sst(SNIP_settings%espcd_sst_dir, &
                nc, nr, landmask, SNIP_arrays%sst, &
                yyyy, mm, dd, hh, ierr)
        else if (SNIP_settings%source_of_ocean_data == "GOFS") then
           write (LDT_logunit,*) &
         '[INFO] CALLING PROCESS_GOFS_SST TO GET SEA SURFACE TEMPERATURES'
           call process_gofs_sst(SNIP_settings%gofs_sst_dir, &
                nc, nr, landmask, SNIP_arrays%sst, &
                yyyy, mm, dd, hh, fh, ierr)
        end if
        if (ierr .ne. 0) then
           ! Fall back on legacy GETSST for 0.25 deg data.
           write (LDT_logunit,*) &
                '[INFO] CALLING GETSST TO GET SEA SURFACE TEMPERATURES'
           call getsst (date10, stmpdir, sstdir)
        end if

        ! RETRIEVE SURFACE OBSERVATIONS.
        write (LDT_logunit,*) &
             '[INFO] CALLING GETOBS TO PROCESS SURFACE OBSERVATIONS'
        ! ALLOCATE ARRAYS USED IN MAX # OF OBS RETURNED.
        allocate (netid            (maxsobs))
        allocate (stadep           (maxsobs))
        allocate (staelv           (maxsobs))
        allocate (staid            (maxsobs))
        allocate (stalat           (maxsobs))
        allocate (stalon           (maxsobs))

        if (sfcobsfmt == 1 .or. sfcobsfmt == 2) then
           call getobs (date10, month, sfcobs, netid, staid, stacnt, &
                stalat, stalon, staelv, stadep, sfcobsfmt)
        else
           write(LDT_logunit,*)'[ERR] Invalid sfcobs file format!'
           write(LDT_logunit,*)'[ERR] Expected 1 (old) or 2 (new)'
           write(LDT_logunit,*)'[ERR] Received ', sfcobsfmt
           call LDT_endrun()
        end if

        write(LDT_logunit,*) &
             '[INFO] TOTAL OBSERVATIONS RETURNED FROM GETOBS: ', &
             stacnt

        ! Copy observations to bratseth object
        do j = 1, stacnt
           network10 = trim(netid(j))
           platform32 = trim(staid(j))
           rob = stadep(j)
           rlat = real(stalat(j)) * 0.01
           rlon = real(stalon(j)) * 0.01
           relev = real(staelv(j))
           call bratseth%append_ob(network10, platform32, rob, &
                rlat, rlon,&
                relev, &
                SNIP_settings%ob_err_var, back=-1.)
        end do

        ! Clean up
        deallocate(netid)
        deallocate(staid)
        deallocate(stalat)
        deallocate(stalon)
        deallocate(staelv)
        deallocate(stadep)

        ! Choose between ESPC-D and GOFS.
        read (date10(1: 4), '(i4)', err=4200) yyyy
        read (date10(5: 6), '(i2)', err=4200) mm
        read (date10(7: 8), '(i2)', err=4200) dd
        read (date10(9:10), '(i2)', err=4200) hh
        fh = 0 ! Dummy value
        found_navy_cice = .false.
        if (SNIP_settings%source_of_ocean_data == "ESPC-D") then
           write(LDT_logunit,*) &
              '[INFO] CALLING PROCESS_ESPCD_CICE TO GET GOFS SEA ICE DATA'
           call process_espcd_cice(SNIP_settings%espcd_cice_dir, &
                nc, nr, landmask, SNIP_arrays%navy_icecon, &
                yyyy, mm, dd, hh, ierr)
           if (ierr == 0) found_navy_cice = .true.
        else if (SNIP_settings%source_of_ocean_data == "GOFS") then
           ! Try to get the GOFS sea ice data
           write(LDT_logunit,*) &
               '[INFO] CALLING PROCESS_GOFS_CICE TO GET GOFS SEA ICE DATA'
           call process_gofs_cice(SNIP_settings%gofs_cice_dir, &
                nc, nr, landmask, SNIP_arrays%navy_icecon, &
                yyyy, mm, dd, hh, fh, ierr)
           if (ierr == 0) found_navy_cice = .true.
        end if

        ! Read externally generated AMSR2 snow depth retrievals.
        call SNIP_read_netcdf_amsr2_sd(date10, ierr)

        ! RETRIEVE VIIRS DATA.
        if (SNIP_settings%useviirs) then
           write (LDT_logunit,*) &
                '[INFO] CALLING GETVIIRS TO GET VIIRS SNOW MAP'
           call getviirs(date10, viirsdir)
        end if

        ! PERFORM THE SNOW ANALYSIS.
        SNIP_arrays%snoanl(:,:) = -1
        SNIP_arrays%icecon(:,:) = -1
        SNIP_arrays%icemask(:,:) = -1
        write(LDT_logunit,*) &
             '[INFO] CALLING RUN_SNOW_ANALYSIS_NOGLACIER'
        call run_snow_analysis_noglacier(runcycle, nc, nr, landmask, &
             landice, elevations, sfctmp_found, sfctmp, bratseth)

        write(LDT_logunit,*) &
             '[INFO] CALLING RUN_SNOW_ANALYSIS_GLACIER'
        call run_snow_analysis_glacier(runcycle, nc, nr, landmask, &
             landice)

        ! Use ESPC-D or GOFS data.  SSMIS based analysis is retired.
        if (found_navy_cice) then
           write(LDT_logunit,*) &
                '[INFO] CALLING RUN_SEAICE_ANALYSIS_NAVY'
           call run_seaice_analysis_navy(month, runcycle, nc, nr, &
                landmask)
        else
           write(LDT_logunit,*)'[ERR] Navy sea ice not available!'
           call LDT_endrun()
        end if

        ! PRINT OUT TOTAL NUMBER OF STATIONS PROCESSED. PRINT ANY STATION
        ! WITH A SNOW DEPTH.  THIS WILL HELP IDENTIFY REGIONS MOST LIKELY
        ! TO REQUIRE MANUAL MODIFICATIONS.  ANY BAD OB SHOULD SHOW HERE.
        write(LDT_logunit,6800) bratseth%count_all_obs()
        call bratseth%sort_obs_by_id()
        call bratseth%print_snowdepths(SNIP_settings%minprt)

        ! Write out in netcdf format
        call SNIP_write_netcdf(date10)

        ! Clean up bratseth object
        call bratseth%delete()

     end do cycle_loop

     ! DEALLOCATE ALL ARRAYS.
     deallocate (SNIP_arrays%climo)
     deallocate (SNIP_arrays%elevat)
     deallocate (SNIP_arrays%iceage)
     deallocate (SNIP_arrays%iceage12z)
     deallocate (SNIP_arrays%icecon)
     deallocate (SNIP_arrays%icemask)
     deallocate (SNIP_arrays%oldcon)
     deallocate (SNIP_arrays%olddep)
     deallocate (SNIP_arrays%oldmask)
     deallocate (SNIP_arrays%ptlat)
     deallocate (SNIP_arrays%ptlon)
     deallocate (sfctmp)
     deallocate (SNIP_arrays%snoage)
     deallocate (SNIP_arrays%snoage12z)
     deallocate (SNIP_arrays%snoanl)
     deallocate (SNIP_arrays%snofrac)
     deallocate (SNIP_arrays%snow_poss)
     deallocate (SNIP_arrays%sst)
     deallocate (SNIP_arrays%viirsmap)
     deallocate (SNIP_arrays%navy_icecon)
     deallocate(landmask)
     deallocate(elevations)
     deallocate(landice)
  end if

#if (defined SPMD)
  call mpi_barrier(mpi_comm_world,ierr)
#endif

  write (LDT_logunit,*) '[INFO] NORMAL TERMINATION'
  return

  ! ERROR HANDLING SECTION.
4200 continue

  message(1) = ' [ERR] ERROR CONVERTING DATA FROM CHARACTER TO INTEGER'
  message(2) = ' [ERR] DATE10 = ' // date10
  call abort_message (program_name, routine_name, message)
  call LDT_endrun()

  ! FORMAT STATEMENTS.
6800 format (/, 1X, 55('-'),                                            &
       /, 3X, '[INFO] PROGRAM:  SNIP',                                  &
       /, 5X, '[INFO] TOTAL STATIONS PROCESSED = ', I5,                 &
       /, 5X, '[INFO] PRINTING DEPTH REPORTS BY NETWORK & STATION ID'   &
       /, 5X, '[INFO] ELEV OF "-1000" INDICATES ELEVATION NOT REPORTED' &
       /, 1X, 55('-'))
6600 format (/, '[INFO] CYCLE DTG = ', A10)

contains

  ! Internal subroutine for reading LDT parameter file and extracting
  ! several fields
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  subroutine read_params(nc, nr, landmask, elevation, landice)

    ! Imports
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_coreMod, only: LDT_rc
    use LDT_logMod, only: LDT_logunit, LDT_verify
    use LDT_paramDataMod, only: LDT_LSMparam_struc
    use netcdf

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(out) :: nc
    integer, intent(out) :: nr
    real, allocatable, intent(out) :: landmask(:,:)
    real, allocatable, intent(out) :: elevation(:,:)
    real, allocatable, intent(out) :: landice(:,:)

    ! Local variables
    character(LDT_CONST_PATH_LEN) :: filename
    integer :: ncid
    integer :: dimids(3)
    integer :: landmask_varid, elevation_varid, surfacetype_varid
    integer :: n
    logical :: file_exists
    integer :: nx, ny, nsfctypes
    real, allocatable :: surfacetype(:,:,:)
    integer :: glacierclass
    integer :: c, r

    n = 1 ! Just use outermost domain
    nc = LDT_rc%lnc(n)
    nr = LDT_rc%lnr(n)

    write(LDT_logunit,*) '[INFO] LDT dimensions: ',nc,' ',nr

    ! Construct filename
    filename = trim(LDT_LSMparam_struc(n)%param_filename)

    ! See if file exists
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
       write(LDT_logunit,*) &
            '[ERR] Cannot find ', trim(filename), ' for SNIP analysis'
       write(LDT_logunit,*)'LDT will stop'
       call LDT_endrun()
    end if

    ! Open the file
    write(ldt_logunit,*)'[INFO] Opening ', trim(filename)
    call LDT_verify(nf90_open(path=trim(filename), &
         mode=nf90_nowrite, ncid=ncid), &
         '[ERR] Cannot open '//trim(filename))

    ! Get the relevant dimensions
    call LDT_verify(nf90_inq_dimid(ncid, "east_west", dimids(1)), &
         '[ERR] nf90_inq_dimid failed')
    call LDT_verify(nf90_inquire_dimension(ncid, dimids(1), len=nx), &
         '[ERR] nf90_inquire_dimension failed')
    if (nx .ne. nc) then
       write(LDT_logunit,*)'[ERR] east_west dimension does not match!'
       write(LDT_logunit,*)'Expected ', nc, ', found ', nx
       call LDT_endrun()
    end if

    call LDT_verify(nf90_inq_dimid(ncid, "north_south", dimids(2)), &
         '[ERR] nf90_inq_dimid failed')
    call LDT_verify(nf90_inquire_dimension(ncid, dimids(2), len=ny), &
         '[ERR] nf90_inquire_dimension failed')
    if (ny .ne. nr) then
       write(LDT_logunit,*)'[ERR] north_south dimension does not match!'
       write(LDT_logunit,*)'Expected ', nr, ', found ', ny
       call LDT_endrun()
    end if

    call LDT_verify(nf90_inq_dimid(ncid, "sfctypes", dimids(3)), &
         '[ERR] nf90_inq_dimid failed')
    call LDT_verify(nf90_inquire_dimension(ncid, dimids(3), &
         len=nsfctypes), &
         '[ERR] nf90_inquire_dimension failed')

    ! Get the landmask
    call LDT_verify(nf90_inq_varid(ncid, 'LANDMASK', landmask_varid), &
         '[ERR] LANDMASK field not found in '//trim(filename))
    allocate(landmask(nc,nr))
    call LDT_verify(nf90_get_var(ncid, landmask_varid, landmask), &
         '[ERR] nf90_get_var failed')

    ! Get the elevation
    call LDT_verify(nf90_inq_varid(ncid, 'ELEVATION', elevation_varid), &
         '[ERR] ELEVATION field not found in '//trim(filename))
    allocate(elevation(nc,nr))
    call LDT_verify(nf90_get_var(ncid, elevation_varid, elevation), &
         '[ERR] nf90_get_var failed')

    ! Get the surfacetype
    call LDT_verify(nf90_inq_varid(ncid, 'SURFACETYPE', &
         surfacetype_varid), &
         '[ERR] SURFACETYPE field not found in '//trim(filename))
    allocate(surfacetype(nc,nr,nsfctypes))
    call LDT_verify(nf90_get_var(ncid, surfacetype_varid, surfacetype), &
         '[ERR] nf90_get_var failed')

    ! Fetch the glacier class number
    call LDT_verify(nf90_get_att(ncid, nf90_global, "GLACIERCLASS", &
         glacierclass), &
         '[ERR] GLACIERCLASS global attribute not found in '// &
         trim(filename))

    ! Create a landice mask
    allocate(landice(nc,nr))
    do r = 1,nr
       do c = 1,nc
          landice(c,r) = surfacetype(c, r, glacierclass)
       end do ! c
    end do ! r

    ! Clean up
    call LDT_verify(nf90_close(ncid))
    deallocate(surfacetype)

  end subroutine read_params

#else
  ! Dummy version
  subroutine read_params(nc, nr, landmask, elevations, landice)
    use LDT_logunit, only: LDT_logunit, LDT_endrun
    implicit none
    integer, intent(out) :: nc,nr
    real, allocatable, intent(out) :: landmask(:,:)
    real, allocatable, intent(out) :: elevations(:,:)
    real, allocatable, intent(out) :: landice(:,:)
    write(LDT_logunit,*) &
         '[ERR] LDT was not compiled with NETCDF support!'
    write(LDT_logunit,*) &
         '[ERR] Recompile with NETCDF support and try again!'
    call LDT_endrun()
  end subroutine read_params
#endif

end subroutine SNIP_run
