!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

! Main LDTSI driver
subroutine LDTSI_run()
   
   !*****************************************************************************************
   !*****************************************************************************************
   !**
   !**  NAME: SNOW DEPTH ANALYSIS MODEL
   !**
   !**  PURPOSE:  DRIVER ROUTINE FOR THE SNOW DEPTH ANALYSIS MODEL.
   !**
   !**   See LDTSI_arraysMod.F90 for common array descriptions
   !**   See LDTSI_paramsMod.F90 for common parameter descriptions
   !**
   !**  UPDATES
   !**  =======
   !**  26 APR 95  INITIAL UNISYS VERSION.....................................SSGT CONRY/SYSM
   !**  09 FEB 96  CHANGED ARGUMENTS PASSED BETWEEN GETOBS, GETSNO, GETSMI,
   !**             AND PROCES DUE TO CHANGES DESCRIBED IN THEIR PROLOGUES.
   !**             FORMATTED PRINT STATEMENTS........................DR KOPP, SSGT CONRY/SYSM
   !**  23 OCT 97  ADDED INCLUDE OF TUNES PROC AND CALL TO RDTUNE.............SSGT CONRY/DNXM
   !**  23 FEB 01  PORTED FROM UNISYS MAINFRAME..............................SSGT MILLER/DNXM
   !**  30 JAN 03  ADDED ROUTINE BLACKLIST......................................MR GAYNO/DNXM
   !**  09 JUL 03  REMOVED CALL TO GETSST; NOW PART OF NT2ICE...................MR GAYNO/DNXM
   !**  05 FEB 04  CHANGED GEOG RECORD LENGTH FROM 4-BYTE TO 2-BYTE WORD
   !**             (INCORPORATED 18 APR 03 AER CHANGES)......................MR LEWISTON/DNXM
   !**  05 MAY 04  MERGED IN THE SSMIS NASA TEAM II SEA ICE CONCENTRATION
   !**             CAPABILITY.  CONVERTED TO FREE FORMAT FORTRAN 90 USING
   !**             MODULES RATHER THAN COMMON BLOCKS.........................MR EYLANDER/DNXM
   !**  21 OCT 04  PORTED TO CDFS II.  REMOVED GRIB SUBROUTINES SINCE
   !**             CDFS II WRITING C++ GRIB ROUTINES.........................MR EYLANDER/DNXM
   !**  21 MAR 05  REMOVED OBSOLETE ARRAY SMITMP.............................MR LEWISTON/DNXM
   !**  14 APR 05  SEPARATED SNOW AND ICE INTO SEPARATE ARRAYS.
   !**             ADDED CALL TO GETICE......................................MR LEWISTON/DNXM
   !**  07 JUN 05  REMOVED JULHR FROM GETOBS & GETSMI CALLS..................MR LEWISTON/DNXM
   !**  08 MAY 09  ADDED AMSR-E SNOW. REMOVED SEA ICE....................MR LEWISTON/2WXG/WEA
   !**  16 JUN 09  CONVERTED TO EQUIDISTANT CYLINDRICAL GRID.............MR LEWISTON/2WXG/WEA
   !**  18 Dec 09  ADDED GRIB OUTPUT.....................................MR LEWISTON/16WS/WXE
   !**  05 JAN 10  CHANGED SSMIS FROM SDRS TO EDRS.......................MR LEWISTON/16WS/WXE
   !**  13 SEP 10  ADDED CALL TO GETSFC..................................MR LEWISTON/16WS/WXE
   !**  29 APR 11  MOVED BLACKLIST AND VALIDATION TO DBPULL..............MR LEWISTON/16WS/WXE
   !**  08 NOV 11  REMOVED AMSR-E AND PORTED TO LINUX....................MR LEWISTON/16WS/WXE
   !**  20 MAR 12  SWITCHED OBS SOURCE FROM CDMS TO JMOBS................MR LEWISTON/16WS/WXE
   !**  29 MAY 12  MOVED MODIFIED TEST FROM SCRIPT TO GETSNO.............MR LEWISTON/16WS/WXE
   !**  07 NOV 12  ADDED FRACTIONAL SNOW.................................MR LEWISTON/16WS/WXE
   !**  10 OCT 13  CHANGED DAILY TO 4 TIMES/DAY. ADDED CYCLE_LOOP.
   !**             REPLACED NAMELIST WITH GETENV CALLS...................MR LEWISTON/16WS/WXE
   !**  15 FEB 17  ADDED VIIRS DATA........................................MR PUSKAR/16WS/WXE
   !**  28 DEC 17  ADDED FILTER FOR MISSING ELEVATIONS...................MR LEWISTON/16WS/WXE
   !**  25 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
   !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
   !**
   !*****************************************************************************************
   !*****************************************************************************************
   
   ! Imports
   use LDT_bratsethMod 
   use LDT_coreMod, only: LDT_masterproc, LDT_rc
   use LDT_logMod, only: LDT_logunit, LDT_endrun
   use LDT_pluginIndices
   use LDT_ldtsiMod, only: ldtsi_settings 
   use map_utils 
#if ( defined SPMD )
   use mpi
#endif
   use LDTSI_analysisMod 
   use LDTSI_arraysMod, only: LDTSI_arrays
   use LDTSI_gofsMod
   use LDTSI_lisMod, only:  read_gr2_t2
   use LDTSI_netcdfMod 
   use LDTSI_paramsMod
   use LDTSI_ssmisMod, only: LDTSI_proc_ssmis
   use LDTSI_utilMod 

   ! Defaults
   implicit none

   ! Local variables
   character*10               ::  date10               ! DATE-TIME GROUP OF CYCLE
   character*100              ::  fracdir              ! FRACTIONAL SNOW DIRECTORY PATH
   character*90               ::  message    (msglns)  ! ERROR MESSAGE
   character*5,  allocatable  ::  netid      (:)       ! NETWORK ID OF AN OBSERVATION
   character*100              ::  modif                ! PATH TO MODIFIED DATA DIRECTORY
   character*100              ::  sfcobs               ! PATH TO DBPULL SNOW OBS DIRECTORY
   character*100              ::  ssmis                ! SSMIS FILE DIRECTORY PATH
   character*9,  allocatable  ::  staid      (:)       ! STATION ID OF AN OBSERVATION
   character*100              ::  static               ! STATIC FILE DIRECTORY PATH
   character*100              ::  stmpdir              ! SFC TEMP DIRECTORY PATH
   character*100              ::  unmod                ! PATH TO UNMODIFIED DATA DIRECTORY
   character*100              ::  viirsdir             ! PATH TO VIIRS DATA DIRECTORY
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
   real,       allocatable    ::  sfctmp_lis (:, :)    ! LIS SHELTER TEMPERATURE DATA
   real,       allocatable    ::  stadep     (:)       ! OBSERVATION SNOW DEPTH (METERS)      
   character*12 :: routine_name 
   type(LDT_bratseth_t) :: bratseth 
   character*10 :: network10, platform10
   real :: rob, rlat, rlon, relev
   integer :: nc,nr
   real, allocatable :: landmask(:,:)
   real, allocatable :: elevations(:,:)
   real, allocatable :: landice(:,:)
   integer :: j
   character*20, allocatable :: blacklist_stns(:)
   character*120 :: line
   integer :: num_blacklist_stns
   logical :: found_blacklist_file
   integer :: icount
   integer :: c, r
   real :: arctlatr
   real, allocatable :: climo_tmp(:,:)
   integer :: maxsobs
   integer :: yyyy, mm, dd, hh, fh
   integer :: ierr
   logical :: found_gofs_cice
   logical :: just_12z

   ! SSMIS snow depth, Yeosang Yoon
   character*100              ::  ssmis_raw_dir       ! SSMIS RAW FILE DIRECTORY PATH   
   integer                    ::  ssmis_option        ! option for snow depth retrieval algorithm

   maxsobs = ldtsi_settings%maxsobs

   ! Only the master process handles the file output
   if (LDT_masterproc) then
      routine_name = "LDTSI"
      message = ' '
      hemi    = 3

      write (LDT_logunit,*) '[INFO] RUNNING LDTSI MODEL'

      ! Check the LDT map projection.  
      ! FIXME:  Support other projections in addition to LATLON
      if (trim(LDT_rc%lis_map_proj) .ne. LDT_latlonId) then
         write(LDT_logunit,*)'[ERR] LDTSI only supports lat/lon grid'
         call LDT_endrun()
      end if

      ! Get relevant fields from LDT parameter file
      call read_params(nc,nr,landmask,elevations, landice)

      ! Copy ldt.config files to local variables
      date10 = trim(ldtsi_settings%date10)
      fracdir = trim(ldtsi_settings%fracdir)
      modif = trim(ldtsi_settings%modif)
      sfcobs = trim(ldtsi_settings%sfcobs)
      ssmis = trim(ldtsi_settings%ssmis)
      stmpdir = trim(ldtsi_settings%stmpdir)
      static = trim(ldtsi_settings%static)
      unmod = trim(ldtsi_settings%unmod)
      viirsdir = trim(ldtsi_settings%viirsdir)
      
      ! for SSMIS snow depth, Yeosang Yoon
      ssmis_raw_dir = trim(ldtsi_settings%ssmis_raw_dir)

      ! FIXME Specify blacklist filename in ldt.config
      ! EMK Read blacklist file.
      inquire(file='LDTSI_blacklist.cfg',exist=found_blacklist_file)
      num_blacklist_stns = 0
      if (found_blacklist_file) then
         open(2, file='LDTSI_blacklist.cfg')
         ! First, count the number of lines in the file
         icount = 0
         do
            read(2,*,end=10)
            icount = icount + 1
         end do
10       continue
         rewind(2)
         num_blacklist_stns = icount - 2 ! Skip first and last line
         if (num_blacklist_stns > 0) then
            allocate(blacklist_stns(num_blacklist_stns))
            ! Now, read each line, other than the first and last
            read(2,*)
            do j = 1, num_blacklist_stns
               read(2,*) line
               blacklist_stns(j) = line(1:len_trim(line))
            end do
         end if
         close(2)
      end if

      ! EXTRACT MONTH FROM DATE-TIME GROUP.      
      read (date10(5:6), '(i2)', err=4200) month

      ! ALLOCATE GRIDDED DATA ARRAYS.
      allocate (LDTSI_arrays%climo            (nc,     nr))
      allocate (LDTSI_arrays%elevat           (nc,     nr))
      allocate (LDTSI_arrays%iceage           (nc,     nr))
      allocate (LDTSI_arrays%iceage12z        (nc,     nr))
      allocate (LDTSI_arrays%icecon           (nc,     nr))
      allocate (LDTSI_arrays%icemask          (nc,     nr))
      allocate (LDTSI_arrays%oldcon           (nc,     nr))
      allocate (LDTSI_arrays%olddep           (nc,     nr))
      allocate (LDTSI_arrays%oldmask          (nc,     nr))
      allocate (LDTSI_arrays%ptlat            (nc,     nr))
      allocate (LDTSI_arrays%ptlon            (nc,     nr))
      ! This is the LIS data interpolated to the LDT grid
      allocate (sfctmp_lis       (nc,     nr))
      allocate (LDTSI_arrays%snoage           (nc,     nr))
      allocate (LDTSI_arrays%snoage12z        (nc,     nr))
      allocate (LDTSI_arrays%snoanl           (nc,     nr))
      allocate (LDTSI_arrays%snofrac          (nc,     nr))
      allocate (LDTSI_arrays%snow_poss        (nc,     nr))
      allocate (LDTSI_arrays%ssmis_depth      (nc,     nr))
      allocate (LDTSI_arrays%ssmis_icecon     (nc,     nr))
      allocate (LDTSI_arrays%sst              (nc,     nr))
      allocate (LDTSI_arrays%viirsmap         (nc,     nr))
      allocate (LDTSI_arrays%gofs_icecon(nc,nr))

      ! RETRIEVE STATIC DATA SETS.
      write (LDT_logunit,*) '[INFO] CALLING GETGEO TO GET STATIC FIELDS'
      ! EMK...Pass LDT elevations to this routine to populate SNODEP elevat
      call getgeo (month, static, nc, nr, elevations)

      ! EMK...Mismatches can occur in landmask between 0.25 deg climo and
      ! LDT's grid.  In cases where LDT introduces "new" land, use a bogus
      ! value.
      arctlatr = float(arctlat) / 100.0
      allocate(climo_tmp(nc,nr))
      climo_tmp(:,:) = LDTSI_arrays%climo(:,:)
      do r = 1, nr
         do c = 1, nc
            if (LDTSI_arrays%ptlat(c,r) > -40.0 .and. &
                 LDTSI_arrays%ptlat(c,r) < 40.0) cycle
            ! See if climo exists for glacier point.  If not, use a fill-in.
            !if (landmask(c,r) > 0.5 .and. &
            !     LDTSI_arrays%climo(c,r) .le. 0) then
            if (landmask(c,r) > 0.5 .and. &
                 LDTSI_arrays%climo(c,r) .le. 0) then
               if (landice(c,r) > 0.5) then
                  climo_tmp(c,r) = ldtsi_settings%fill_climo
               else
                  !climo_tmp(c,r) = ldtsi_settings%unkdep
                  climo_tmp(c,r) = 0 
               end if
            end if
         end do ! c
      end do ! r

      LDTSI_arrays%climo(:,:) = climo_tmp(:,:)
      deallocate(climo_tmp)

      ! RETRIEVE THE PREVIOUS SNOW ANALYSIS.
      ! First, try reading LDTSI in netCDF format.  If that doesn't work,
      ! fall back on the legacy SNODEP at 0.25 deg resolution.
      write (LDT_logunit,*) &
           '[INFO] CALLING GETSNO_NC TO GET PREVIOUS SNOW AND ICE DATA'
      call getsno_nc(date10, julhr_beg, ierr)
      if (ierr .ne. 0) then
         write (LDT_logunit,*) &
              '[INFO] CALLING GETSNO TO GET PREVIOUS SNOW AND ICE DATA'
         if (ierr == 1) then
            just_12z = .true.
         else
            just_12z = .false.
         end if
         call getsno (date10, modif, unmod, nc, nr, landice, julhr_beg, &
              just_12z)
      end if

      julhr = julhr_beg
      call date10_julhr (date10, julhr_end, program_name, routine_name)

      ! LOOP THROUGH UNCOMPLETED CYCLES.
      ! FIXME...Change this to handle only one cycle
      !YY cycle_loop : do while (julhr .lt. julhr_end)
      cycle_loop : do while (julhr .lt. julhr_end)

         ! EMK Create bratseth object
         call bratseth%new(maxsobs, ldtsi_settings%back_err_var, &
              ldtsi_settings%back_err_h_corr_len, &
              ldtsi_settings%back_err_v_corr_len)

         julhr = julhr + 6
         call julhr_date10 (julhr, date10, program_name, routine_name)

         write(LDT_logunit,*) &
              '[INFO] ***PREPARING SNOW/ICE ANALYSIS VALID ', date10, &
              '***'

         write (ldt_logunit,6600) date10
         read (date10(9:10), '(i2)', err=4200) runcycle

         ! RETRIEVE LIS SURFACE TEMPERATURES.
         ! FIXME...Try reading netCDF output from prior LIS history, and 
         ! pass back flag indicating success or failure.  If it fails, 
         ! fall back on GETSFC for 0.25deg data.         
         write (LDT_logunit,*) &
              '[INFO] CALLING READ_GR2_T2 TO GET LIS SHELTER TEMPERATURES'
         call read_gr2_t2(date10, nc, nr, sfctmp_lis, ierr)
         if (ierr .ne. 0) then
            write (LDT_logunit,*) &
                 '[INFO] CALLING GETSFC TO GET LIS SHELTER TEMPERATURES'
            call getsfc (date10, stmpdir, sfctmp_found, sfctmp_lis)
         else
            sfctmp_found = .true.
         end if

         ! RETRIEVE NAVY SEA SURFACE TEMPERATURE (SST) DATA.
         ! First try the US Navy 0.08 deg GOFS data
         read (date10(1: 4), '(i4)', err=4200) yyyy
         read (date10(5: 6), '(i2)', err=4200) mm
         read (date10(7: 8), '(i2)', err=4200) dd
         read (date10(9:10), '(i2)', err=4200) hh
         fh = 0 ! Dummy value
         write (LDT_logunit,*) &
              '[INFO] CALLING PROCESS_GOFS_SST TO GET SEA SURFACE TEMPERATURES'
         call process_gofs_sst(ldtsi_settings%gofs_sst_dir, &
              nc, nr, landmask, LDTSI_arrays%sst, &
              yyyy, mm, dd, hh, fh, ierr)
         if (ierr .ne. 0) then
            ! Fall back on legacy GETSST for 0.25 deg data.
            write (LDT_logunit,*) &
                 '[INFO] CALLING GETSST TO GET SEA SURFACE TEMPERATURES'
            call getsst (date10, stmpdir)
         end if

         ! RETRIEVE FRACTIONAL SNOW DATA.         
         if (ldtsi_settings%usefrac) then               
            write (LDT_logunit,*) &
                 '[INFO] CALLING GETFRAC TO GET FRACTIONAL SNOW DATA'
            call getfrac (date10, fracdir)               
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

         call getobs (date10, month,  sfcobs, netid, staid, stacnt, &
              stalat, stalon, staelv, stadep)         
         write(LDT_logunit,*) &
              '[INFO] TOTAL OBSERVATIONS RETURNED FROM GETOBS: ', &
              stacnt

         ! EMK Copy observations to bratseth object
         do j = 1, stacnt
            network10 = trim(netid(j))
            platform10 = trim(staid(j))
            rob = stadep(j)
            rlat = real(stalat(j)) * 0.01
            rlon = real(stalon(j)) * 0.01
            relev = real(staelv(j))
            call bratseth%append_ob(network10,platform10,rob,rlat,rlon,&
                 relev, &
                 ldtsi_settings%ob_err_var,back=-1.)
         end do

         ! Clean up
         deallocate(netid)
         deallocate(staid)
         deallocate(stalat)
         deallocate(stalon)
         deallocate(staelv)
         deallocate(stadep)

         ! Try to get the GOFS sea ice data
         write(LDT_logunit,*) &
              '[INFO] CALLING PROCESS_GOFS_CICE TO GET GOFS SEA ICE DATA'
         read (date10(1: 4), '(i4)', err=4200) yyyy
         read (date10(5: 6), '(i2)', err=4200) mm
         read (date10(7: 8), '(i2)', err=4200) dd
         read (date10(9:10), '(i2)', err=4200) hh
         fh = 0 ! Dummy value
         call process_gofs_cice(ldtsi_settings%gofs_cice_dir, &
              nc, nr, landmask, LDTSI_arrays%gofs_icecon, &
              yyyy, mm, dd, hh, fh, ierr)
         if (ierr == 0) then
            found_gofs_cice = .true.
         else
            found_gofs_cice = .false.
         end if
        
         ! Estimates SSMIS-based snow depth, Yeosang Yoon
         write (LDT_logunit,*) &
              '[INFO] CALLING LDTSI_PROC_SSMIS'
         call LDTSI_proc_ssmis(date10, ssmis_raw_dir, ssmis, &
              ldtsi_settings%ssmis_option) 
         
         ! RETRIEVE SSMIS DATA.         
         write (LDT_logunit,*) &
              '[INFO] CALLING GETSMI TO GET SSMIS EDR VALUES'
         call getsmi (date10, ssmis)

         ! RETRIEVE VIIRS DATA.
         if (ldtsi_settings%useviirs) then
            write (LDT_logunit,*) &
                 '[INFO] CALLING GETVIIRS TO GET VIIRS SNOW MAP'
            call getviirs(date10, viirsdir)
         end if

         ! PERFORM THE SNOW ANALYSIS.
         LDTSI_arrays%snoanl(:,:) = -1
         LDTSI_arrays%icecon(:,:) = -1
         LDTSI_arrays%icemask(:,:) = -1
         write(LDT_logunit,*) &
              '[INFO] CALLING RUN_SNOW_ANALYSIS_NOGLACIER'
         call run_snow_analysis_noglacier(runcycle,nc,nr,landmask, landice,  &
              elevations, sfctmp_found, sfctmp_lis, &
              num_blacklist_stns, blacklist_stns, bratseth)
         if (allocated(blacklist_stns)) deallocate(blacklist_stns)

         write(LDT_logunit,*) &
              '[INFO] CALLING RUN_SNOW_ANALYSIS_GLACIER'
         call run_snow_analysis_glacier(runcycle, nc, nr,landmask, landice)

         ! FIXME...Try using GOFS data first, and if unsuccessful, then run
         ! the old SSMIS analysis.
         if (found_gofs_cice) then
            write(LDT_logunit,*) &
                 '[INFO] CALLING RUN_SEAICE_ANALYSIS_GOFS'
            call run_seaice_analysis_gofs(month,runcycle,nc,nr,landmask)
         else
            write(LDT_logunit,*) &
                 '[INFO] CALLING RUN_SEAICE_ANALYSIS_SSMIS'
            call run_seaice_analysis_ssmis(month,runcycle,nc,nr,landmask)
         end if

         ! PRINT OUT TOTAL NUMBER OF STATIONS PROCESSED. PRINT ANY STATION
         ! WITH A SNOW DEPTH.  THIS WILL HELP IDENTIFY REGIONS MOST LIKELY
         ! TO REQUIRE MANUAL MODIFICATIONS.  ANY BAD OB SHOULD SHOW HERE.
         write(LDT_logunit,6800) bratseth%count_all_obs()
         call bratseth%sort_obs_by_id()
         call bratseth%print_snowdepths(ldtsi_settings%minprt)

         ! EMK Write out in netcdf format
         call LDTSI_write_netcdf(date10)

         ! Clean up bratseth object
         call bratseth%delete()

      !YY end do cycle_loop
      end do cycle_loop

      ! DEALLOCATE ALL ARRAYS.
      deallocate (ldtsi_arrays%climo)
      deallocate (ldtsi_arrays%elevat)
      deallocate (ldtsi_arrays%iceage)
      deallocate (ldtsi_arrays%iceage12z)
      deallocate (ldtsi_arrays%icecon)
      deallocate (ldtsi_arrays%icemask)
      deallocate (ldtsi_arrays%oldcon)
      deallocate (ldtsi_arrays%olddep)
      deallocate (ldtsi_arrays%oldmask)
      deallocate (ldtsi_arrays%ptlat)
      deallocate (ldtsi_arrays%ptlon)
      deallocate (sfctmp_lis)
      deallocate (ldtsi_arrays%snoage)
      deallocate (ldtsi_arrays%snoage12z)
      deallocate (ldtsi_arrays%snoanl)
      deallocate (ldtsi_arrays%snofrac)
      deallocate (ldtsi_arrays%snow_poss)
      deallocate (ldtsi_arrays%ssmis_depth)
      deallocate (ldtsi_arrays%ssmis_icecon)
      deallocate (ldtsi_arrays%sst)
      deallocate (ldtsi_arrays%viirsmap)
      deallocate (ldtsi_arrays%gofs_icecon)
      deallocate(landmask)
      deallocate(elevations)
      deallocate(landice)
   end if

#if (defined SPMD)
   call mpi_barrier(mpi_comm_world,ierr)
#endif

   write (LDT_logunit,*) '[INFO] NORMAL TERMINATION'      
   call LDT_endrun()

   ! ERROR HANDLING SECTION.   
4200 continue

   message(1) = ' [ERR] ERROR CONVERTING DATA FROM CHARACTER TO INTEGER'
   message(2) = ' [ERR] DATE10 = ' // date10
   call abort_message (program_name, program_name, message)
   call LDT_endrun()

   ! FORMAT STATEMENTS.   
6800 format (/, 1X, 55('-'),                                             &
        /, 3X, '[INFO] PROGRAM:  LDTSI',                                 &
        /, 5X, '[INFO] TOTAL STATIONS PROCESSED = ', I5,                 &
        /, 5X, '[INFO] PRINTING DEPTH REPORTS BY NETWORK & STATION ID'   &
        /, 5X, '[INFO] ELEV OF "-1000" INDICATES ELEVATION NOT REPORTED' &
        /, 1X, 55('-'))   
6600 format (/, '[INFO] CYCLE DTG = ', A10)

contains

   ! Internal subroutine for reading LDT parameter file and extracting
   ! several fields
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   subroutine read_params(nc,nr,landmask,elevation,landice)
      
      ! Imports
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
      character*100 :: filename
      integer :: ncid
      integer :: dimids(3)
      integer :: landmask_varid, elevation_varid, &
           surfacetype_varid
      integer :: n
      logical :: file_exists
      integer :: nx, ny, nsfctypes
      real, allocatable :: surfacetype(:,:,:)
      integer :: glacierclass
      integer :: c,r

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
              '[ERR] Cannot find ',trim(filename),' for LDTSI analysis'
         write(LDT_logunit,*)'LDT will stop'
         call LDT_endrun()
      end if
      
      ! Open the file
      call LDT_verify(nf90_open(path=trim(filename), &
           mode=nf90_nowrite, ncid=ncid), &
           '[ERR] Cannot open '//trim(filename))
      
      ! Get the relevant dimensions
      call LDT_verify(nf90_inq_dimid(ncid,"east_west",dimids(1)), &
           '[ERR] nf90_inq_dimid failed')
      call LDT_verify(nf90_inquire_dimension(ncid,dimids(1), len=nx), &
           '[ERR] nf90_inquire_dimension failed')
      if (nx .ne. nc) then
         write(LDT_logunit,*)'[ERR] east_west dimension does not match!'
         write(LDT_logunit,*)'Expected ',nc,', found ',nx
         call LDT_endrun()
      end if
      
      call LDT_verify(nf90_inq_dimid(ncid,"north_south",dimids(2)), &
           '[ERR] nf90_inq_dimid failed')
      call LDT_verify(nf90_inquire_dimension(ncid,dimids(2), len=ny), &
           '[ERR] nf90_inquire_dimension failed')
      if (ny .ne. nr) then
         write(LDT_logunit,*)'[ERR] north_south dimension does not match!'
         write(LDT_logunit,*)'Expected ',nr,', found ',ny
         call LDT_endrun()
      end if

      call LDT_verify(nf90_inq_dimid(ncid,"sfctypes",dimids(3)), &
           '[ERR] nf90_inq_dimid failed')
      call LDT_verify(nf90_inquire_dimension(ncid,dimids(3), len=nsfctypes), &
           '[ERR] nf90_inquire_dimension failed')
      
      ! Get the landmask
      call LDT_verify(nf90_inq_varid(ncid,'LANDMASK',landmask_varid), &
           '[ERR] LANDMASK field not found in '//trim(filename))
      allocate(landmask(nc,nr))
      call LDT_verify(nf90_get_var(ncid,landmask_varid,landmask), &
           '[ERR] nf90_get_var failed')
      
      ! Get the elevation
      call LDT_verify(nf90_inq_varid(ncid,'ELEVATION',elevation_varid), &
           '[ERR] ELEVATION field not found in '//trim(filename))
      allocate(elevation(nc,nr))
      call LDT_verify(nf90_get_var(ncid,elevation_varid,elevation), &
           '[ERR] nf90_get_var failed')

      ! Get the surfacetype
      call LDT_verify(nf90_inq_varid(ncid,'SURFACETYPE',surfacetype_varid), &
           '[ERR] SURFACETYPE field not found in '//trim(filename))
      allocate(surfacetype(nc,nr,nsfctypes))
      call LDT_verify(nf90_get_var(ncid,surfacetype_varid,surfacetype), &
           '[ERR] nf90_get_var failed')
      
      ! Fetch the glacier class number
      call LDT_verify(nf90_get_att(ncid,nf90_global,"GLACIERCLASS", &
           glacierclass), &
           '[ERR] GLACIERCLASS global attribute not found in '//trim(filename))

      ! Create a landice mask
      allocate(landice(nc,nr))
      do r = 1,nr
         do c = 1,nc
            landice(c,r) = surfacetype(c,r,glacierclass)
         end do ! c
      end do ! r

      ! Clean up
      call LDT_verify(nf90_close(ncid))
      deallocate(surfacetype)

   end subroutine read_params
   
#else
   ! Dummy version
   subroutine read_params(nc,nr,landmask,elevations,landice)
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
      
end subroutine LDTSI_run
