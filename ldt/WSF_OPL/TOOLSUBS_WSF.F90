!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: TOOLSUBS_WSF
!
! DESCRIPTION: Module for reading WSF NetCDF data with ALL 17 channels
!              FIXED: Division by zero protection in quality flag creation
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
        tb_lowres, lat, lon, land_frac_low, quality_flag, &  ! This should be 3D output
        earth_inc_angle, snow, precip, &
        nscans, nfovs, nchans, & 
        chan_frequencies, chan_polarizations, ierr)
    
    character(*), intent(in) :: filename
    real*4, allocatable, intent(out) :: tb_lowres(:,:,:)
    real*4, allocatable, intent(out) :: lat(:,:)
    real*4, allocatable, intent(out) :: lon(:,:)
    real*4, allocatable, intent(out) :: land_frac_low(:,:)
    integer*1, allocatable, intent(out) :: quality_flag(:,:)
    real*4, allocatable, intent(out) :: earth_inc_angle(:,:,:)
    integer*4, allocatable, intent(out) :: snow(:,:)
    integer*4, allocatable, intent(out) :: precip(:,:)
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
    real :: sil, tt18, denom
    
    ! Temporary arrays for channel extraction (Fortran order)
    real*4, allocatable :: tb_18v(:,:), tb_18h(:,:)
    real*4, allocatable :: tb_23v(:,:), tb_36v(:,:), tb_89v(:,:)
    real*4, allocatable :: chan_freq(:)
    character*1, allocatable :: chan_pol(:)
    integer*1, allocatable :: qf_from_file(:,:,:)
    
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
    
    ! Allocate arrays
    write(LDT_logunit,*)'[INFO] Allocating arrays (Fortran column-major order):'
    
    allocate(tb_lowres(nfovs, nscans, nchans))
    allocate(lat(nfovs, nscans))
    allocate(lon(nfovs, nscans))
    allocate(land_frac_low(nfovs, nscans))
    allocate(qf_from_file(nfovs, nscans, nband))
    allocate(earth_inc_angle(nfovs, nscans, 1))
    allocate(snow(nfovs, nscans))
    allocate(precip(nfovs, nscans))
    allocate(quality_flag(nfovs, nscans))
    allocate(chan_frequencies(nchans))
    allocate(chan_polarizations(nchans))
    
    ! Temporary arrays
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
    earth_inc_angle = 52.0
    snow = 0
    precip = 0
    quality_flag = 0
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Reading variables from NetCDF...'
    
    ! Read Latitude
    ierr = nf90_inq_varid(ncid, 'Latitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lat)
        write(LDT_logunit,*)'[INFO] ✓ Read Latitude'
    else
        write(LDT_logunit,*)'[WARN] Could not read Latitude'
    end if
    
    ! Read Longitude
    ierr = nf90_inq_varid(ncid, 'Longitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lon)
        write(LDT_logunit,*)'[INFO] ✓ Read Longitude'
    else
        write(LDT_logunit,*)'[WARN] Could not read Longitude'
    end if
    
    ! Read TbLowRes
    ierr = nf90_inq_varid(ncid, 'TbLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, tb_lowres)
        
        ! Replace NaN and invalid values with -9999.0 using explicit loops
        do ichan = 1, nchans
            do i = 1, nscans
                do j = 1, nfovs
                    ! Check for NaN (value != itself)
                    if (tb_lowres(j,i,ichan) /= tb_lowres(j,i,ichan) .or. &
                        tb_lowres(j,i,ichan) < 0.0 .or. &
                        tb_lowres(j,i,ichan) > 400.0) then
                        tb_lowres(j,i,ichan) = -9999.0
                    endif
                end do
            end do
        end do
        
        write(LDT_logunit,*)'[INFO] ✓ Read TbLowRes - ALL', nchans, 'channels'
        write(LDT_logunit,*)'[INFO]   Replaced NaN/invalid values with -9999.0'
    else
        write(LDT_logunit,*)'[ERR] Could not read TbLowRes'
        ierr = nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ! Read Land Fraction
    ierr = nf90_inq_varid(ncid, 'LandFractionLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, land_frac_low)
        
        ! Replace NaN with -9999.0
        where (land_frac_low /= land_frac_low)
            land_frac_low = -9999.0
        end where
        
        write(LDT_logunit,*)'[INFO] ✓ Read LandFractionLowRes'
    else
        write(LDT_logunit,*)'[WARN] Could not read LandFractionLowRes'
    end if
    
    ! Read Quality Flag
    ierr = nf90_inq_varid(ncid, 'QualityFlag', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, qf_from_file)
        write(LDT_logunit,*)'[INFO] ✓ Read QualityFlag'
    else
        write(LDT_logunit,*)'[WARN] Could not read QualityFlag, assuming good quality'
        qf_from_file = 0
    end if
    
    ! Read Channel Frequencies
    ierr = nf90_inq_varid(ncid, 'ChanFrequency', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, chan_freq)
        write(LDT_logunit,*)'[INFO] ✓ Read ChanFrequency'
        chan_frequencies(:) = chan_freq(:)
    else
        write(LDT_logunit,*)'[WARN] Could not read ChanFrequency, using defaults'
        if (nchans >= 10) then
            chan_frequencies(1:2) = 10.65
            chan_frequencies(3:4) = 18.7
            chan_frequencies(5:6) = 23.8
            chan_frequencies(7:8) = 36.5
            chan_frequencies(9:10) = 89.0
        end if
    end if
    
    ! Read Channel Polarizations
    ierr = nf90_inq_varid(ncid, 'ChanPolarization', varid)
    if (ierr == NF90_NOERR) then
        do i = 1, nchans
            ierr = nf90_get_var(ncid, varid, chan_pol(i:i), &
                start=(/i/), count=(/1/))
        end do
        write(LDT_logunit,*)'[INFO] ✓ Read ChanPolarization'
        chan_polarizations(:) = chan_pol(:)
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
    
    ! Extract specific channels and handle fill values
    tb_18v = -9999.0
    tb_18h = -9999.0
    tb_23v = -9999.0
    tb_36v = -9999.0
    tb_89v = -9999.0
    
    do ichan = 1, nchans
        if (abs(chan_freq(ichan) - 18.7) < 0.5 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_18v(:,:) = tb_lowres(:,:,ichan)
            write(LDT_logunit,*)'[INFO] ✓ Found TB_18V at channel ', ichan
        end if
        
        if (abs(chan_freq(ichan) - 18.7) < 0.5 .and. &
            (chan_pol(ichan) == 'h' .or. chan_pol(ichan) == 'H')) then
            tb_18h(:,:) = tb_lowres(:,:,ichan)
            write(LDT_logunit,*)'[INFO] ✓ Found TB_18H at channel ', ichan
        end if
        
        if (abs(chan_freq(ichan) - 23.8) < 0.5 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_23v(:,:) = tb_lowres(:,:,ichan)
            write(LDT_logunit,*)'[INFO] ✓ Found TB_23V at channel ', ichan
        end if
        
        if (abs(chan_freq(ichan) - 36.5) < 1.0 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_36v(:,:) = tb_lowres(:,:,ichan)
            write(LDT_logunit,*)'[INFO] ✓ Found TB_36V at channel ', ichan
        end if
        
        if (abs(chan_freq(ichan) - 89.0) < 2.0 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_89v(:,:) = tb_lowres(:,:,ichan)
            write(LDT_logunit,*)'[INFO] ✓ Found TB_89V at channel ', ichan
        end if
    end do
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Creating quality flags...'
    
    ! =====================================================================
    ! CREATE QUALITY FLAGS - USE PROPER FILL VALUE CHECKS
    ! =====================================================================
    do j = 1, nfovs
        do i = 1, nscans
            quality_flag(j,i) = 0
            
            ! Bit 0: Ocean (land_frac < 0.5)
            if (land_frac_low(j,i) < 0.5) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 1)
            endif
            if (land_frac_low(j,i) > 0.5) then
                ! Bit 1: Precipitation detection
                ! Only compute if all required TBs are valid (> 0 means not fill value)
                if (tb_18v(j,i) > 0.0 .and. tb_18h(j,i) > 0.0 .and. &
                    tb_23v(j,i) > 0.0) then
                    
                    denom = tb_18v(j,i) + tb_18h(j,i)
                    sil = (tb_18v(j,i) - tb_18h(j,i)) / denom
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
                
                ! ================================================================
                ! BAND-SPECIFIC SENSOR QUALITY FLAGS (Bits 3-7)
                ! Check bits 0-4 of sensor quality flag for each band
                ! ================================================================
                
                ! Bit 3: Band 1 (10 GHz) sensor quality
                if (nband >= 1) then
                    if (IAND(INT(qf_from_file(j,i,1)), 7) /= 0) then
                        quality_flag(j,i) = IOR(quality_flag(j,i), 8)  ! 2^3 = 8
                    endif
                endif
                
                ! Bit 4: Band 2 (18 GHz) sensor quality
                if (nband >= 2) then
                    if (IAND(INT(qf_from_file(j,i,2)), 7) /= 0) then
                        quality_flag(j,i) = IOR(quality_flag(j,i), 16)  ! 2^4 = 16
                    endif
                endif
                
                ! Bit 5: Band 3 (23 GHz) sensor quality
                if (nband >= 3) then
                    if (IAND(INT(qf_from_file(j,i,3)), 7) /= 0) then
                        quality_flag(j,i) = IOR(quality_flag(j,i), 32)  ! 2^5 = 32
                    endif
                endif
                
                ! Bit 6: Band 4 (36 GHz) sensor quality
                if (nband >= 4) then
                    if (IAND(INT(qf_from_file(j,i,4)), 7) /= 0) then
                        quality_flag(j,i) = IOR(quality_flag(j,i), 64)  ! 2^6 = 64
                    endif
                endif
                
                ! Bit 7: Band 5 (89 GHz) sensor quality
                if (nband >= 5) then
                    if (IAND(INT(qf_from_file(j,i,5)), 7) /= 0) then
                        quality_flag(j,i) = IOR(quality_flag(j,i), 128)  ! 2^7 = 128
                    endif
                endif
            endif
        end do
    end do
    
    write(LDT_logunit,*)'[INFO] ✓ 8-bit combined quality flags created'
    write(LDT_logunit,*)'[INFO]   Bit 0: Ocean (land_frac < 50%)'
    write(LDT_logunit,*)'[INFO]   Bit 1: Precipitation'
    write(LDT_logunit,*)'[INFO]   Bit 2: Snow'
    write(LDT_logunit,*)'[INFO]   Bit 3: Sensor quality for 10GHz (Band 1)'
    write(LDT_logunit,*)'[INFO]   Bit 4: Sensor quality for 18GHz (Band 2)'
    write(LDT_logunit,*)'[INFO]   Bit 5: Sensor quality for 23GHz (Band 3)'
    write(LDT_logunit,*)'[INFO]   Bit 6: Sensor quality for 36GHz (Band 4)'
    write(LDT_logunit,*)'[INFO]   Bit 7: Sensor quality for 89GHz (Band 5)'
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[DEBUG] Ocean flags set: ', count(IBITS(quality_flag, 0, 1) == 1)
    write(LDT_logunit,*)'[DEBUG] Precip flags set: ', count(IBITS(quality_flag, 1, 1) == 1)
    write(LDT_logunit,*)'[DEBUG] Snow flags set: ', count(IBITS(quality_flag, 2, 1) == 1)
    write(LDT_logunit,*)'[DEBUG] Land frac range: ', minval(land_frac_low), maxval(land_frac_low)
    
    ! Cleanup temporary arrays
    deallocate(tb_18v, tb_18h, tb_23v, tb_36v, tb_89v)
    deallocate(chan_freq, chan_pol)
    deallocate(qf_from_file)
    
    write(LDT_logunit,*)'[INFO] ✓ Successfully read ALL WSF data with flags'
    ierr = 0
    
#else
    write(LDT_logunit,*) '[ERR] get_wsf_data_with_flags called without NetCDF support!'
    call LDT_endrun()
    ierr = 1
#endif
    
    END SUBROUTINE get_wsf_data_with_flags

END MODULE TOOLSUBS_WSF