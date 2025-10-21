!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: TOOLSUBS_WSF
!
! DESCRIPTION: Module for reading WSF NetCDF data with ALL 17 channels
!              FIXED: Corrected Fortran dimension order (column-major)
!              READS: All channels, not just AMSR2 subset
!
!-------------------------------------------------------------------------

#include "LDT_misc.h"

MODULE TOOLSUBS_WSF
    USE LDT_logMod, only: LDT_logunit, LDT_endrun
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    USE netcdf
#endif
    IMPLICIT NONE

CONTAINS

    SUBROUTINE get_wsf_data_with_flags(filename, &
        tb_lowres, lat, lon, land_frac_low, quality_flag, &
        earth_inc_angle, snow, precip, &
        nscans, nfovs, nchans, chan_frequencies, chan_polarizations, ierr)
    
    ! =====================================================================
    ! CORRECTED DIMENSION ORDER FOR FORTRAN (column-major)
    ! ncdump shows: TbLowRes(nChan, nScanR, nFOVR)
    ! Fortran reads: TbLowRes(nFOVR, nScanR, nChan) - REVERSED!
    ! =====================================================================
    
    character(*), intent(in) :: filename
    real*4, allocatable, intent(out) :: tb_lowres(:,:,:)   ! (nFOVR, nScanR, nChan) - FORTRAN ORDER
    real*4, allocatable, intent(out) :: lat(:,:)           ! (nFOVR, nScanR) - FORTRAN ORDER
    real*4, allocatable, intent(out) :: lon(:,:)           ! (nFOVR, nScanR) - FORTRAN ORDER
    real*4, allocatable, intent(out) :: land_frac_low(:,:) ! (nFOVR, nScanR) - FORTRAN ORDER
    integer*1, allocatable, intent(out) :: quality_flag(:,:) ! (nFOVR, nScanR) - FORTRAN ORDER
    real*4, allocatable, intent(out) :: earth_inc_angle(:,:,:) ! (nFOVR, nScanR, 1) - dummy
    integer*4, allocatable, intent(out) :: snow(:,:)       ! (nFOVR, nScanR) - FORTRAN ORDER
    integer*4, allocatable, intent(out) :: precip(:,:)     ! (nFOVR, nScanR) - FORTRAN ORDER
    real*4, allocatable, intent(out) :: chan_frequencies(:)
    character*1, allocatable, intent(out) :: chan_polarizations(:)
    integer, intent(out) :: nscans, nfovs, nchans
    integer, intent(out) :: ierr
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    ! Local variables
    integer :: ncid, varid
    integer :: nscanr_dimid, nfovr_dimid, nchan_dimid, nband_dimid
    integer :: nband
    logical :: file_exists
    integer :: i, j, ichan
    real :: sil, tt18
    
    ! Temporary arrays for channel extraction (Fortran order)
    real*4, allocatable :: tb_18v(:,:), tb_18h(:,:)
    real*4, allocatable :: tb_23v(:,:), tb_36v(:,:), tb_89v(:,:)
    real*4, allocatable :: chan_freq(:)
    character*1, allocatable :: chan_pol(:)
    integer*1, allocatable :: qf_from_file(:,:,:)  ! (nFOVR, nScanR, nBand) - FORTRAN ORDER
    
    ierr = 0
    nscans = 0
    nfovs = 0
    nchans = 0
    nband = 0
    
    ! Check if file exists
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
        write(LDT_logunit,*)'[ERR] Cannot find file ', trim(filename)
        ierr = 1
        return
    end if
    
    ! Open the NetCDF file
    ierr = nf90_open(trim(filename), NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot open file ', trim(filename)
        write(LDT_logunit,*)nf90_strerror(ierr)
        ierr = 1
        return
    end if
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] WSF NetCDF Reader - ALL CHANNELS'
    write(LDT_logunit,*)'[INFO] File: ', trim(filename)
    
    ! Get dimensions
    ierr = nf90_inq_dimid(ncid, 'nScanR', nscanr_dimid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find dimension nScanR'
        ierr = nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inquire_dimension(ncid, nscanr_dimid, len=nscans)
    ierr = nf90_inq_dimid(ncid, 'nFOVR', nfovr_dimid)
    ierr = nf90_inquire_dimension(ncid, nfovr_dimid, len=nfovs)
    ierr = nf90_inq_dimid(ncid, 'nChan', nchan_dimid)
    ierr = nf90_inquire_dimension(ncid, nchan_dimid, len=nchans)
    ierr = nf90_inq_dimid(ncid, 'nBand', nband_dimid)
    ierr = nf90_inquire_dimension(ncid, nband_dimid, len=nband)
    
    write(LDT_logunit,*)'[INFO] NetCDF Dimensions (C order from ncdump):'
    write(LDT_logunit,*)'[INFO]   nScanR = ', nscans
    write(LDT_logunit,*)'[INFO]   nFOVR  = ', nfovs
    write(LDT_logunit,*)'[INFO]   nChan  = ', nchans, ' (ALL channels)'
    write(LDT_logunit,*)'[INFO]   nBand  = ', nband
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! =====================================================================
    ! ALLOCATE WITH FORTRAN DIMENSION ORDER (REVERSED from ncdump)
    ! =====================================================================
    write(LDT_logunit,*)'[INFO] Allocating arrays (Fortran column-major order):'
    
    ! Main data arrays - FORTRAN ORDER
    allocate(tb_lowres(nfovs, nscans, nchans))
    write(LDT_logunit,*)'[INFO]   tb_lowres(', nfovs, ',', nscans, ',', nchans, ')'
    
    allocate(lat(nfovs, nscans))
    write(LDT_logunit,*)'[INFO]   lat(', nfovs, ',', nscans, ')'
    
    allocate(lon(nfovs, nscans))
    write(LDT_logunit,*)'[INFO]   lon(', nfovs, ',', nscans, ')'
    
    allocate(land_frac_low(nfovs, nscans))
    write(LDT_logunit,*)'[INFO]   land_frac_low(', nfovs, ',', nscans, ')'
    
    allocate(qf_from_file(nfovs, nscans, nband))
    write(LDT_logunit,*)'[INFO]   quality_flag_file(', nfovs, ',', nscans, ',', nband, ')'
    
    allocate(earth_inc_angle(nfovs, nscans, 1))
    allocate(snow(nfovs, nscans))
    allocate(precip(nfovs, nscans))
    allocate(quality_flag(nfovs, nscans))
    
    ! Channel info arrays
    allocate(chan_frequencies(nchans))
    allocate(chan_polarizations(nchans))
    
    ! Temporary arrays for snow/precip detection
    allocate(tb_18v(nfovs, nscans))
    allocate(tb_18h(nfovs, nscans))
    allocate(tb_23v(nfovs, nscans))
    allocate(tb_36v(nfovs, nscans))
    allocate(tb_89v(nfovs, nscans))
    
    allocate(chan_freq(nchans))
    allocate(chan_pol(nchans))
    
    ! Initialize arrays
    tb_lowres = 0.0
    lat = 0.0
    lon = 0.0
    land_frac_low = 0.0
    qf_from_file = 0
    earth_inc_angle = 52.0  ! Default value (not used)
    snow = 0
    precip = 0
    quality_flag = 0
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Reading variables from NetCDF...'
    
    ! =====================================================================
    ! READ LATITUDE (nScanR, nFOVR) -> Fortran gets (nFOVR, nScanR)
    ! =====================================================================
    ierr = nf90_inq_varid(ncid, 'Latitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lat)
        write(LDT_logunit,*)'[INFO] ✓ Read Latitude'
    else
        write(LDT_logunit,*)'[WARN] Could not read Latitude'
    end if
    
    ! =====================================================================
    ! READ LONGITUDE (nScanR, nFOVR) -> Fortran gets (nFOVR, nScanR)
    ! =====================================================================
    ierr = nf90_inq_varid(ncid, 'Longitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lon)
        write(LDT_logunit,*)'[INFO] ✓ Read Longitude'
    else
        write(LDT_logunit,*)'[WARN] Could not read Longitude'
    end if
    
    ! =====================================================================
    ! READ TbLowRes - ALL 17 CHANNELS
    ! ncdump: TbLowRes(nChan, nScanR, nFOVR) = (17, 636, 286)
    ! Fortran: TbLowRes(nFOVR, nScanR, nChan) = (286, 636, 17)
    ! =====================================================================
    ierr = nf90_inq_varid(ncid, 'TbLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, tb_lowres)
        write(LDT_logunit,*)'[INFO] ✓ Read TbLowRes - ALL', nchans, 'channels'
        write(LDT_logunit,*)'[INFO]   Fortran array shape: (', nfovs, ',', nscans, ',', nchans, ')'
    else
        write(LDT_logunit,*)'[ERR] Could not read TbLowRes'
        ierr = nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ! =====================================================================
    ! READ LAND FRACTION (nScanR, nFOVR) -> Fortran gets (nFOVR, nScanR)
    ! =====================================================================
    ierr = nf90_inq_varid(ncid, 'LandFractionLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, land_frac_low)
        write(LDT_logunit,*)'[INFO] ✓ Read LandFractionLowRes'
    else
        write(LDT_logunit,*)'[WARN] Could not read LandFractionLowRes'
    end if
    
    ! =====================================================================
    ! READ QUALITY FLAG (nBand, nScanR, nFOVR) -> Fortran gets (nFOVR, nScanR, nBand)
    ! =====================================================================
    ierr = nf90_inq_varid(ncid, 'QualityFlag', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, qf_from_file)
        write(LDT_logunit,*)'[INFO] ✓ Read QualityFlag'
    else
        write(LDT_logunit,*)'[WARN] Could not read QualityFlag, assuming good quality'
        qf_from_file = 0
    end if
    
    ! =====================================================================
    ! READ CHANNEL FREQUENCIES (nChan) = (17)
    ! =====================================================================
    ierr = nf90_inq_varid(ncid, 'ChanFrequency', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, chan_freq)
        write(LDT_logunit,*)'[INFO] ✓ Read ChanFrequency'
        chan_frequencies(:) = chan_freq(:)
        write(LDT_logunit,*)'[INFO] Channel Frequencies (GHz):'
        do ichan = 1, min(nchans, 17)
            write(LDT_logunit,*)'[INFO]   Chan', ichan, ': ', chan_freq(ichan), ' GHz'
        end do
    else
        write(LDT_logunit,*)'[WARN] Could not read ChanFrequency, using defaults'
        ! Use WSF standard frequencies for first 10 channels
        if (nchans >= 10) then
            chan_frequencies(1:2) = 10.65
            chan_frequencies(3:4) = 18.7
            chan_frequencies(5:6) = 23.8
            chan_frequencies(7:8) = 36.5
            chan_frequencies(9:10) = 89.0
        end if
    end if
    
    ! =====================================================================
    ! READ CHANNEL POLARIZATIONS (nChan) = (17)
    ! =====================================================================
    ierr = nf90_inq_varid(ncid, 'ChanPolarization', varid)
    if (ierr == NF90_NOERR) then
        do i = 1, nchans
            ierr = nf90_get_var(ncid, varid, chan_pol(i:i), &
                start=(/i/), count=(/1/))
        end do
        write(LDT_logunit,*)'[INFO] ✓ Read ChanPolarization'
        chan_polarizations(:) = chan_pol(:)
        write(LDT_logunit,*)'[INFO] Channel Polarizations:'
        do ichan = 1, min(nchans, 17)
            write(LDT_logunit,*)'[INFO]   Chan', ichan, ': ', chan_pol(ichan)
        end do
    else
        write(LDT_logunit,*)'[WARN] Could not read ChanPolarization, using defaults'
        if (nchans >= 10) then
            chan_polarizations(1) = 'V'
            chan_polarizations(2) = 'H'
            chan_polarizations(3) = 'V'
            chan_polarizations(4) = 'H'
            chan_polarizations(5) = 'V'
            chan_polarizations(6) = 'H'
            chan_polarizations(7) = 'V'
            chan_polarizations(8) = 'H'
            chan_polarizations(9) = 'V'
            chan_polarizations(10) = 'H'
        end if
    end if
    
    ! Close file
    ierr = nf90_close(ncid)
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Extracting channels for snow/precip detection...'
    
    ! =====================================================================
    ! EXTRACT SPECIFIC CHANNELS FOR SNOW/PRECIP DETECTION
    ! Access pattern: tb_lowres(ifov, iscan, ichan) - FORTRAN ORDER
    ! =====================================================================
    tb_18v = 0.0
    tb_18h = 0.0
    tb_23v = 0.0
    tb_36v = 0.0
    tb_89v = 0.0
    
    do ichan = 1, nchans
        ! 18.7 GHz V
        if (abs(chan_freq(ichan) - 18.7) < 0.5 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_18v(:,:) = tb_lowres(:,:,ichan)  ! CORRECTED indexing
            write(LDT_logunit,*)'[INFO] ✓ Found TB_18V at channel ', ichan
        end if
        
        ! 18.7 GHz H
        if (abs(chan_freq(ichan) - 18.7) < 0.5 .and. &
            (chan_pol(ichan) == 'h' .or. chan_pol(ichan) == 'H')) then
            tb_18h(:,:) = tb_lowres(:,:,ichan)  ! CORRECTED indexing
            write(LDT_logunit,*)'[INFO] ✓ Found TB_18H at channel ', ichan
        end if
        
        ! 23.8 GHz V
        if (abs(chan_freq(ichan) - 23.8) < 0.5 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_23v(:,:) = tb_lowres(:,:,ichan)  ! CORRECTED indexing
            write(LDT_logunit,*)'[INFO] ✓ Found TB_23V at channel ', ichan
        end if
        
        ! 36.5 GHz V
        if (abs(chan_freq(ichan) - 36.5) < 1.0 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_36v(:,:) = tb_lowres(:,:,ichan)  ! CORRECTED indexing
            write(LDT_logunit,*)'[INFO] ✓ Found TB_36V at channel ', ichan
        end if
        
        ! 89.0 GHz V
        if (abs(chan_freq(ichan) - 89.0) < 2.0 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_89v(:,:) = tb_lowres(:,:,ichan)  ! CORRECTED indexing
            write(LDT_logunit,*)'[INFO] ✓ Found TB_89V at channel ', ichan
        end if
    end do
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Creating quality flags...'
    
    ! =====================================================================
    ! CREATE QUALITY FLAGS
    ! Access pattern: (ifov, iscan) - FORTRAN ORDER
    ! =====================================================================
    do j = 1, nfovs
        do i = 1, nscans
            quality_flag(j,i) = 0
            
            ! Bit 0: Ocean (land_frac < 0.2)
            if (land_frac_low(j,i) < 0.2) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 1)
            endif
            
            ! Bit 1: Precipitation detection
            if (tb_18v(j,i) > 0.0 .and. tb_18h(j,i) > 0.0 .and. &
                tb_23v(j,i) > 0.0) then
                sil = (tb_18v(j,i) - tb_18h(j,i)) / &
                      (tb_18v(j,i) + tb_18h(j,i))
                tt18 = (tb_18v(j,i) - tb_23v(j,i))
                if ((sil < 0.005) .and. (tt18 < -2.0)) then
                    precip(j,i) = 1
                    quality_flag(j,i) = IOR(quality_flag(j,i), 2)
                endif
            endif
            
            ! Bit 2: Snow detection
            if (tb_36v(j,i) > 0.0 .and. tb_89v(j,i) > 0.0) then
                if (tb_36v(j,i) - tb_89v(j,i) < -5.0) then
                    snow(j,i) = 1
                    quality_flag(j,i) = IOR(quality_flag(j,i), 4)
                endif
            endif
            
            ! Bit 3: Sensor quality (check bits 0-4 of band 0)
            ! nBand dimension in Fortran order: (nFOVR, nScanR, nBand)
            if (IAND(INT(qf_from_file(j,i,1)), 31) /= 0) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 8)
            endif
        end do
    end do
    
    write(LDT_logunit,*)'[INFO] ✓ Quality flags created'
    write(LDT_logunit,*)'[INFO]   Bit 0: Ocean (land_frac < 20%)'
    write(LDT_logunit,*)'[INFO]   Bit 1: Precipitation'
    write(LDT_logunit,*)'[INFO]   Bit 2: Snow'
    write(LDT_logunit,*)'[INFO]   Bit 3: Sensor quality (bits 0-4 check)'
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Cleanup temporary arrays
    deallocate(tb_18v, tb_18h, tb_23v, tb_36v, tb_89v)
    deallocate(chan_freq, chan_pol)
    deallocate(qf_from_file)
    
    write(LDT_logunit,*)'[INFO] ✓ Successfully read ALL WSF data with flags'
    write(LDT_logunit,*)'[INFO] ========================================='
    ierr = 0
    
#else
    ! Dummy version if LDT was compiled w/o NetCDF support
    write(LDT_logunit,*) '[ERR] get_wsf_data_with_flags called without NetCDF support!'
    write(LDT_logunit,*) '[ERR] Recompile LDT with NetCDF support and try again!'
    call LDT_endrun()
    ierr = 1
#endif
    
    END SUBROUTINE get_wsf_data_with_flags

END MODULE TOOLSUBS_WSF