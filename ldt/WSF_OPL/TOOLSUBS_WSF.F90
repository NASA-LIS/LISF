!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: TOOLSUBS_WSF
!
! DESCRIPTION: Module for reading WSF NetCDF data with proper band-specific
!              quality flag handling
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
    integer :: i, j, ichan, iband
    real :: sil, tt18, denom
    
    ! Temporary arrays for channel extraction (Fortran order)
    real*4, allocatable :: tb_18v(:,:), tb_18h(:,:)
    real*4, allocatable :: tb_23v(:,:), tb_36v(:,:), tb_89v(:,:)
    real*4, allocatable :: chan_freq(:)
    character*1, allocatable :: chan_pol(:)
    integer*1, allocatable :: qf_from_file(:,:,:)  ! (nfovs, nscans, nbands)
    integer*1, allocatable :: validation_flags(:,:,:)  ! ValidationConditionType
    integer*1, allocatable :: degradation_flags(:,:,:) ! DegradationConditionType
    integer*1, allocatable :: exclusion_flags(:,:,:)   ! ExclusionConditionType
    integer*1, allocatable :: band_quality_flags(:,:,:)  ! NEW: (nfovs, nscans, 6 bands)

    
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
    write(LDT_logunit,*)'[INFO] WSF NetCDF Reader - FIXED BAND QUALITY'
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
    
    write(LDT_logunit,*)'[INFO] NetCDF Dimensions:'
    write(LDT_logunit,*)'[INFO]   nScanR = ', nscans
    write(LDT_logunit,*)'[INFO]   nFOVR  = ', nfovs
    write(LDT_logunit,*)'[INFO]   nChan  = ', nchans
    write(LDT_logunit,*)'[INFO]   nBand  = ', nband
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Allocate arrays
    allocate(tb_lowres(nfovs, nscans, nchans))
    allocate(lat(nfovs, nscans))
    allocate(lon(nfovs, nscans))
    allocate(land_frac_low(nfovs, nscans))
    allocate(qf_from_file(nfovs, nscans, nband))
    allocate(validation_flags(nfovs, nscans, nband))
    allocate(degradation_flags(nfovs, nscans, nband))
    allocate(exclusion_flags(nfovs, nscans, nband))
    allocate(earth_inc_angle(nfovs, nscans, 1))
    allocate(snow(nfovs, nscans))
    allocate(precip(nfovs, nscans))
    allocate(quality_flag(nfovs, nscans))
    allocate(band_quality_flags(nfovs, nscans, 6))  ! 6 bands max
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
    validation_flags = 0
    degradation_flags = 0
    exclusion_flags = 0
    earth_inc_angle = 52.0
    snow = 0
    precip = 0
    quality_flag = 0
    band_quality_flags = 0
    
    ! Read variables
    write(LDT_logunit,*)'[INFO] Reading variables from NetCDF...'
    
    ! Read Latitude
    ierr = nf90_inq_varid(ncid, 'Latitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lat)
        write(LDT_logunit,*)'[INFO] ✓ Read Latitude'
    end if
    
    ! Read Longitude
    ierr = nf90_inq_varid(ncid, 'Longitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lon)
        write(LDT_logunit,*)'[INFO] ✓ Read Longitude'
    end if
    
    ! Read TbLowRes
    ierr = nf90_inq_varid(ncid, 'TbLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, tb_lowres)
        
        ! Replace NaN and invalid values with -9999.0
        do ichan = 1, nchans
            do i = 1, nscans
                do j = 1, nfovs
                    if (tb_lowres(j,i,ichan) /= tb_lowres(j,i,ichan) .or. &
                        tb_lowres(j,i,ichan) < 0.0 .or. &
                        tb_lowres(j,i,ichan) > 400.0) then
                        tb_lowres(j,i,ichan) = -9999.0
                    endif
                end do
            end do
        end do
        write(LDT_logunit,*)'[INFO] ✓ Read TbLowRes'
    end if
    
    ! Read Land Fraction
    ierr = nf90_inq_varid(ncid, 'LandFractionLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, land_frac_low)
        where (land_frac_low /= land_frac_low)
            land_frac_low = -9999.0
        end where
        write(LDT_logunit,*)'[INFO] ✓ Read LandFractionLowRes'
    end if
    
    ! Read QualityFlag (main quality flag per band)
    ierr = nf90_inq_varid(ncid, 'QualityFlag', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, qf_from_file)
        write(LDT_logunit,*)'[INFO] ✓ Read QualityFlag (3D array)'
        
        ! Debug: Show what's in the quality flags
        write(LDT_logunit,*)'[DEBUG] QualityFlag sample values (first pixel, each band):'
        do iband = 1, min(nband, 6)
            write(LDT_logunit,'(A,I1,A,I3)') '[DEBUG]   Band ', iband, ': ', &
                                             INT(qf_from_file(1,1,iband))
        end do
    else
        write(LDT_logunit,*)'[WARN] Could not read QualityFlag'
        qf_from_file = 0
    end if
    
    ! Read ValidationConditionType
    ierr = nf90_inq_varid(ncid, 'ValidationConditionType', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, validation_flags)
        write(LDT_logunit,*)'[INFO] ✓ Read ValidationConditionType'
    else
        write(LDT_logunit,*)'[INFO] ValidationConditionType not found'
    end if
    
    ! Read DegradationConditionType
    ierr = nf90_inq_varid(ncid, 'DegradationConditionType', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, degradation_flags)
        write(LDT_logunit,*)'[INFO] ✓ Read DegradationConditionType'
    else
        write(LDT_logunit,*)'[INFO] DegradationConditionType not found'
    end if
    
    ! Read ExclusionConditionType
    ierr = nf90_inq_varid(ncid, 'ExclusionConditionType', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, exclusion_flags)
        write(LDT_logunit,*)'[INFO] ✓ Read ExclusionConditionType'
    else
        write(LDT_logunit,*)'[INFO] ExclusionConditionType not found'
    end if
    
    ! Read Channel Frequencies
    ierr = nf90_inq_varid(ncid, 'ChanFrequency', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, chan_freq)
        chan_frequencies(:) = chan_freq(:)
        write(LDT_logunit,*)'[INFO] ✓ Read ChanFrequency'
    end if
    
    ! Read Channel Polarizations
    ierr = nf90_inq_varid(ncid, 'ChanPolarization', varid)
    if (ierr == NF90_NOERR) then
        do i = 1, nchans
            ierr = nf90_get_var(ncid, varid, chan_pol(i:i), &
                start=(/i/), count=(/1/))
        end do
        chan_polarizations(:) = chan_pol(:)
        write(LDT_logunit,*)'[INFO] ✓ Read ChanPolarization'
    end if
    
    ! Close file
    ierr = nf90_close(ncid)
    
    ! Extract specific channels for snow/precip detection
    tb_18v = -9999.0
    tb_18h = -9999.0
    tb_23v = -9999.0
    tb_36v = -9999.0
    tb_89v = -9999.0
    
    do ichan = 1, nchans
        if (abs(chan_freq(ichan) - 18.7) < 0.5) then
            if (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V') then
                tb_18v(:,:) = tb_lowres(:,:,ichan)
            else if (chan_pol(ichan) == 'h' .or. chan_pol(ichan) == 'H') then
                tb_18h(:,:) = tb_lowres(:,:,ichan)
            endif
        else if (abs(chan_freq(ichan) - 23.8) < 0.5) then
            if (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V') then
                tb_23v(:,:) = tb_lowres(:,:,ichan)
            endif
        else if (abs(chan_freq(ichan) - 36.5) < 1.0) then
            if (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V') then
                tb_36v(:,:) = tb_lowres(:,:,ichan)
            endif
        else if (abs(chan_freq(ichan) - 89.0) < 2.0) then
            if (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V') then
                tb_89v(:,:) = tb_lowres(:,:,ichan)
            endif
        endif
    end do
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Creating quality flags...'
    
    ! =====================================================================
    ! CREATE QUALITY FLAGS WITH PROPER BAND-SPECIFIC HANDLING
    ! =====================================================================
    do j = 1, nfovs
        do i = 1, nscans
            quality_flag(j,i) = 0
            
            ! Bit 0: Ocean (land_frac < 0.5)
            if (land_frac_low(j,i) >= 0.0) then
                if (land_frac_low(j,i) < 0.5) then
                    quality_flag(j,i) = IOR(quality_flag(j,i), 1)  ! Ocean
                endif
            endif
            
            ! Only check precip/snow over land
            if (land_frac_low(j,i) > 0.5) then
                ! Bit 1: Precipitation detection
                if (tb_18v(j,i) > 0.0 .and. tb_18h(j,i) > 0.0 .and. &
                    tb_23v(j,i) > 0.0) then
                    denom = tb_18v(j,i) + tb_18h(j,i)
                    if (denom > 0.0) then
                        sil = (tb_18v(j,i) - tb_18h(j,i)) / denom
                        tt18 = (tb_18v(j,i) - tb_23v(j,i))
                        if ((sil < 0.005) .and. (tt18 < -2.0)) then
                            precip(j,i) = 1
                            quality_flag(j,i) = IOR(quality_flag(j,i), 2)
                        endif
                    endif
                endif
                
                ! Bit 2: Snow detection
                if (tb_36v(j,i) > 0.0 .and. tb_89v(j,i) > 0.0) then
                    if (tb_36v(j,i) - tb_89v(j,i) < -5.0) then
                        snow(j,i) = 1
                        quality_flag(j,i) = IOR(quality_flag(j,i), 4)
                    endif
                endif
            endif
            
            ! ================================================================
            ! BAND-SPECIFIC QUALITY FLAGS (Bits 3-7)
            ! Check if each band has bad quality
            ! ================================================================
            
            ! Store per-band quality for later use
            do iband = 1, min(nband, 6)
                ! A band is considered bad if ANY of these conditions are true:
                ! - QualityFlag bit 0 (IsNotValid) is set
                ! - QualityFlag bit 1 (IsExclusionCondition) is set  
                ! - QualityFlag bit 2 (IsDegradationCondition) is set
                
                if (IBITS(INT(qf_from_file(j,i,iband)), 0, 1) == 1 .or. &  ! IsNotValid
                    IBITS(INT(qf_from_file(j,i,iband)), 1, 1) == 1 .or. &  ! IsExclusionCondition
                    !IBITS(INT(qf_from_file(j,i,iband)), 4, 1) == 1 .or. &   ! IsOverlap
                    !IBITS(INT(qf_from_file(j,i,iband)), 3, 1) == 1 .or. &   ! IsLimitedUtilityValidationCondition
                    IBITS(INT(qf_from_file(j,i,iband)), 2, 1) == 1) then   ! IsDegradationCondition
                    band_quality_flags(j,i,iband) = 1  ! Mark band as bad
                else
                    band_quality_flags(j,i,iband) = 0  ! Band is good
                endif
            end do
            
            ! Map band quality to the combined quality flag bits 3-7
            ! Band 1: 10.65 GHz -> bit 3
            if (nband >= 1 .and. band_quality_flags(j,i,1) == 1) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 8)   ! 2^3
            endif
            
            ! Band 2: 18.7 GHz -> bit 4
            if (nband >= 2 .and. band_quality_flags(j,i,2) == 1) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 16)  ! 2^4
            endif
            
            ! Band 3: 23.8 GHz -> bit 5
            if (nband >= 3 .and. band_quality_flags(j,i,3) == 1) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 32)  ! 2^5
            endif
            
            ! Band 4: 36.5 GHz -> bit 6
            if (nband >= 4 .and. band_quality_flags(j,i,4) == 1) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 64)  ! 2^6
            endif
            
            ! Band 5: 89.0 GHz -> bit 7
            if (nband >= 5 .and. band_quality_flags(j,i,5) == 1) then
                quality_flag(j,i) = IOR(quality_flag(j,i), 128) ! 2^7
            endif
        end do
    end do
    
    ! Debug output
    write(LDT_logunit,*)'[INFO] Quality flag statistics:'
    write(LDT_logunit,*)'[INFO]   Ocean pixels: ', count(IBITS(quality_flag, 0, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Precip pixels: ', count(IBITS(quality_flag, 1, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Snow pixels: ', count(IBITS(quality_flag, 2, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 10GHz: ', count(IBITS(quality_flag, 3, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 18GHz: ', count(IBITS(quality_flag, 4, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 23GHz: ', count(IBITS(quality_flag, 5, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 36GHz: ', count(IBITS(quality_flag, 6, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 89GHz: ', count(IBITS(quality_flag, 7, 1) == 1)
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Cleanup temporary arrays
    deallocate(tb_18v, tb_18h, tb_23v, tb_36v, tb_89v)
    deallocate(chan_freq, chan_pol)
    deallocate(qf_from_file)
    deallocate(validation_flags, degradation_flags, exclusion_flags)
    
    ierr = 0
    
#else
    write(LDT_logunit,*) '[ERR] get_wsf_data_with_flags called without NetCDF support!'
    call LDT_endrun()
    ierr = 1
#endif
    
    END SUBROUTINE get_wsf_data_with_flags

END MODULE TOOLSUBS_WSF
