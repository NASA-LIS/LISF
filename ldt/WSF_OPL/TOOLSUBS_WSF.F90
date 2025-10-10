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
! DESCRIPTION: Module for reading WSF NetCDF data with AMSR-style processing
!              CORRECTED: Fixed variable declarations and nf90_close usage
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
    
    ! Arguments  
    character(*), intent(in) :: filename
    real*4, allocatable, intent(out) :: tb_lowres(:,:,:)   ! (nchans, nscans, nfovs)
    real*4, allocatable, intent(out) :: lat(:,:)           ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: lon(:,:)           ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: land_frac_low(:,:) ! (nscans, nfovs)
    integer*1, allocatable, intent(out) :: quality_flag(:,:) ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: earth_inc_angle(:,:,:) ! (nchans, nscans, nfovs)
    integer*4, allocatable, intent(out) :: snow(:,:)       ! (nscans, nfovs)
    integer*4, allocatable, intent(out) :: precip(:,:)     ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: chan_frequencies(:)
    character*1, allocatable, intent(out) :: chan_polarizations(:)
    integer, intent(out) :: nscans, nfovs, nchans
    integer, intent(out) :: ierr
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    ! Local variables - ALL DECLARATIONS MUST BE HERE
    integer :: ncid, varid
    integer :: nscanr_dimid, nfovr_dimid, nchan_dimid
    logical :: file_exists
    integer :: i, j, ichan
    real :: sil, tt18
    real*4, allocatable :: tb_18v(:,:), tb_18h(:,:)
    real*4, allocatable :: tb_23v(:,:), tb_36v(:,:), tb_89v(:,:)
    real*4, allocatable :: chan_freq(:)      ! CORRECTED: Moved here
    character*1, allocatable :: chan_pol(:)  ! CORRECTED: Moved here
    
    ierr = 0
    nscans = 0
    nfovs = 0
    nchans = 0
    
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
    
    write(LDT_logunit,*)'[INFO] Successfully opened ', trim(filename)
    
    ! Get dimensions
    ierr = nf90_inq_dimid(ncid, 'nScanR', nscanr_dimid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find dimension nScanR'
        ierr = nf90_close(ncid)  ! CORRECTED: Function, not subroutine
        ierr = 1
        return
    end if
    
    ierr = nf90_inquire_dimension(ncid, nscanr_dimid, len=nscans)
    ierr = nf90_inq_dimid(ncid, 'nFOVR', nfovr_dimid)
    ierr = nf90_inquire_dimension(ncid, nfovr_dimid, len=nfovs)
    ierr = nf90_inq_dimid(ncid, 'nChan', nchan_dimid)
    ierr = nf90_inquire_dimension(ncid, nchan_dimid, len=nchans)
    
    write(LDT_logunit,*)'[INFO] Dimensions: nScanR=', nscans, &
        ' nFOVR=', nfovs, ' nChan=', nchans
    
    ! Allocate all arrays
    allocate(tb_lowres(nchans, nscans, nfovs))
    allocate(lat(nscans, nfovs))
    allocate(lon(nscans, nfovs))
    allocate(land_frac_low(nscans, nfovs))
    allocate(earth_inc_angle(nchans, nscans, nfovs))
    allocate(snow(nscans, nfovs))
    allocate(precip(nscans, nfovs))
    allocate(quality_flag(nscans, nfovs))
    allocate(chan_frequencies(nchans))
    allocate(chan_polarizations(nchans))
    
    ! Allocate temporary arrays for snow/precip detection
    allocate(tb_18v(nscans, nfovs))
    allocate(tb_18h(nscans, nfovs))
    allocate(tb_23v(nscans, nfovs))
    allocate(tb_36v(nscans, nfovs))
    allocate(tb_89v(nscans, nfovs))
    
    ! Allocate temporary arrays for channel info
    allocate(chan_freq(nchans))
    allocate(chan_pol(nchans))
    
    ! Initialize arrays
    tb_lowres = 0.0
    lat = 0.0
    lon = 0.0
    land_frac_low = 0.0
    earth_inc_angle = 0.0
    snow = 0
    precip = 0
    quality_flag = 0
    
    ! Read latitude
    ierr = nf90_inq_varid(ncid, 'Latitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lat)
        write(LDT_logunit,*)'[INFO] Read Latitude'
    else
        write(LDT_logunit,*)'[WARN] Could not read Latitude'
    end if
    
    ! Read longitude
    ierr = nf90_inq_varid(ncid, 'Longitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lon)
        write(LDT_logunit,*)'[INFO] Read Longitude'
    else
        write(LDT_logunit,*)'[WARN] Could not read Longitude'
    end if
    
    ! Read TbLowRes
    ierr = nf90_inq_varid(ncid, 'TbLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, tb_lowres)
        write(LDT_logunit,*)'[INFO] Read TbLowRes'
    else
        write(LDT_logunit,*)'[WARN] Could not read TbLowRes'
    end if
    
    ! Read LandFractionLowRes
    ierr = nf90_inq_varid(ncid, 'LandFractionLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, land_frac_low)
        write(LDT_logunit,*)'[INFO] Read LandFractionLowRes'
    else
        write(LDT_logunit,*)'[WARN] Could not read LandFractionLowRes'
    end if
    
    ! Read Earth Incidence Angle
    ierr = nf90_inq_varid(ncid, 'EarthIncAngle', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, earth_inc_angle)
        write(LDT_logunit,*)'[INFO] Read EarthIncAngle'
    else
        ! Try alternative name
        ierr = nf90_inq_varid(ncid, 'earth_incidence_angle', varid)
        if (ierr == NF90_NOERR) then
            ierr = nf90_get_var(ncid, varid, earth_inc_angle)
            write(LDT_logunit,*)'[INFO] Read earth_incidence_angle'
        else
            write(LDT_logunit,*)'[WARN] Could not read Earth Incidence Angle'
            earth_inc_angle = 55.0  ! Default value
        end if
    end if
    
    ! Read channel frequencies
    ierr = nf90_inq_varid(ncid, 'ChanFrequency', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, chan_freq)
        write(LDT_logunit,*)'[INFO] Read ChanFrequency'
        chan_frequencies(:) = chan_freq(:)
    else
        write(LDT_logunit,*)'[WARN] Could not read ChanFrequency, using defaults'
        ! Set default frequencies based on WSF standard channels
        if (nchans >= 10) then
            chan_frequencies(1:2) = 10.65
            chan_frequencies(3:4) = 18.7
            chan_frequencies(5:6) = 23.8
            chan_frequencies(7:8) = 36.5
            chan_frequencies(9:10) = 89.0
        end if
    end if
    
    ! Read channel polarizations
    ierr = nf90_inq_varid(ncid, 'ChanPolarization', varid)
    if (ierr == NF90_NOERR) then
        ! CORRECTED: Read as character array properly
        do i = 1, nchans
            ierr = nf90_get_var(ncid, varid, chan_pol(i:i), &
                start=(/i/), count=(/1/))
        end do
        write(LDT_logunit,*)'[INFO] Read ChanPolarization'
        chan_polarizations(:) = chan_pol(:)
    else
        write(LDT_logunit,*)'[WARN] Could not read ChanPolarization, using defaults'
        ! Set default polarizations (V,H pattern)
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
    ierr = nf90_close(ncid)  ! CORRECTED: Function, not subroutine
    
    ! ==================================================================
    ! EXTRACT SPECIFIC CHANNELS FOR SNOW/PRECIP DETECTION
    ! ==================================================================
    write(LDT_logunit,*)'[INFO] Extracting channels for snow/precip detection'
    
    ! Initialize
    tb_18v = 0.0
    tb_18h = 0.0
    tb_23v = 0.0
    tb_36v = 0.0
    tb_89v = 0.0
    
    ! Extract specific channels based on frequency and polarization
    do ichan = 1, nchans
        ! 18.7 GHz V
        if (abs(chan_freq(ichan) - 18.7) < 0.5 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_18v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_18V at channel ', ichan
        end if
        
        ! 18.7 GHz H
        if (abs(chan_freq(ichan) - 18.7) < 0.5 .and. &
            (chan_pol(ichan) == 'h' .or. chan_pol(ichan) == 'H')) then
            tb_18h(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_18H at channel ', ichan
        end if
        
        ! 23.8 GHz V
        if (abs(chan_freq(ichan) - 23.8) < 0.5 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_23v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_23V at channel ', ichan
        end if
        
        ! 36.5 GHz V
        if (abs(chan_freq(ichan) - 36.5) < 1.0 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_36v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_36V at channel ', ichan
        end if
        
        ! 89.0 GHz V
        if (abs(chan_freq(ichan) - 89.0) < 2.0 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_89v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_89V at channel ', ichan
        end if
    end do
    
    ! ==================================================================
    ! SNOW AND PRECIPITATION DETECTION
    ! Exact same algorithm as AMSR_OPL
    ! ==================================================================
    write(LDT_logunit,*)'[INFO] Computing snow and precipitation flags'
    
    snow = 0
    precip = 0
    
    do j = 1, nfovs
        do i = 1, nscans
            ! Only process over land (land_frac >= 50%)
            ! Note: WSF land_frac is 0-1, need to convert to percentage
            if (land_frac_low(i,j) >= 0.5) then
                ! Ensure all needed values are valid
                if (tb_18v(i,j) > 0 .and. tb_18h(i,j) > 0 .and. &
                    tb_23v(i,j) > 0 .and. tb_89v(i,j) > 0 .and. &
                    tb_36v(i,j) > 0) then
                    
                    ! Calculate intermediate variables (EXACT AMSR algorithm)
                    sil = 451.88 - 0.44*tb_18v(i,j) - 1.775*tb_23v(i,j) + &
                          0.00574*tb_23v(i,j)**2 - tb_89v(i,j)
                    tt18 = tb_18v(i,j) - tb_18h(i,j)
                    
                    if (sil > 10) then
                        if ((tb_23v(i,j) <= 264.0) .and. &
                            (tb_23v(i,j) <= (175.0 + 0.49*tb_89v(i,j)))) then
                            ! Snow branch
                            snow(i,j) = 1
                            if ((tt18 >= 18) .and. &
                                ((tb_18v(i,j) - tb_36v(i,j)) <= 10) .and. &
                                ((tb_36v(i,j) - tb_89v(i,j)) <= 10)) then
                                snow(i,j) = 0
                            endif
                            if ((tt18 >= 8) .and. &
                                ((tb_18v(i,j) - tb_36v(i,j)) <= 2) .and. &
                                ((tb_23v(i,j) - tb_89v(i,j)) <= 6)) then
                                snow(i,j) = 1
                            endif
                        else
                            ! Precipitation branch
                            snow(i,j) = 0
                            precip(i,j) = 1
                            if (tt18 > 20) then
                                precip(i,j) = 0
                            endif
                            if ((tb_89v(i,j) > 253) .and. (tt18 > 7)) then
                                precip(i,j) = 0
                            endif
                        endif
                    endif
                endif
            endif
        end do
    end do
    
    ! ==================================================================
    ! CREATE AMSR-STYLE QUALITY FLAG
    ! Bit 0: Ocean (1 if land_frac < 20%)
    ! Bit 1: Precipitation 
    ! Bit 2: Snow
    ! ==================================================================
    write(LDT_logunit,*)'[INFO] Packing quality flags at footprint level'
    
    quality_flag = 0
    
    do j = 1, nfovs
        do i = 1, nscans
            ! Bit 0: Ocean flag (land fraction < 20%)
            ! WSF land_frac is 0-1, so 0.2 = 20%
            if (land_frac_low(i,j) < 0.2) then
                quality_flag(i,j) = IOR(quality_flag(i,j), 1)
            endif
            
            ! Bit 1: Precipitation flag
            if (precip(i,j) == 1) then
                quality_flag(i,j) = IOR(quality_flag(i,j), 2)
            endif
            
            ! Bit 2: Snow flag
            if (snow(i,j) == 1) then
                quality_flag(i,j) = IOR(quality_flag(i,j), 4)
            endif
        end do
    end do
    
    write(LDT_logunit,*)'[INFO] Quality flags packed successfully'
    
    ! Cleanup temporary arrays
    deallocate(tb_18v, tb_18h, tb_23v, tb_36v, tb_89v)
    deallocate(chan_freq, chan_pol)
    
    write(LDT_logunit,*)'[INFO] Successfully read all WSF data with flags'
    ierr = 0
    
#else
    ! Dummy version if LDT was compiled w/o NetCDF support
    write(LDT_logunit,*) '[ERR] get_wsf_data_with_flags called without NetCDF support!'
    write(LDT_logunit,*) '[ERR] Recompile LDT with NetCDF support and try again!'
    call LDT_endrun()
    ierr = 1
#endif
    
    END SUBROUTINE get_wsf_data_with_flags

    ! ==================================================================
    ! Simplified version without snow/precip detection (for compatibility)
    ! ==================================================================
    SUBROUTINE get_wsf_data(filename, &
        tb_lowres, lat, lon, land_frac_low, &
        quality_flag, earth_inc_angle, &
        nscans, nfovs, nchans, ierr)
    
    character(*), intent(in) :: filename
    real*4, allocatable, intent(out) :: tb_lowres(:,:,:)
    real*4, allocatable, intent(out) :: lat(:,:)
    real*4, allocatable, intent(out) :: lon(:,:)
    real*4, allocatable, intent(out) :: land_frac_low(:,:)
    integer*1, allocatable, intent(out) :: quality_flag(:,:,:)
    real*4, allocatable, intent(out) :: earth_inc_angle(:,:,:)
    integer, intent(out) :: nscans, nfovs, nchans
    integer, intent(out) :: ierr
    
    ! Local variables
    integer*4, allocatable :: snow(:,:), precip(:,:)
    integer*1, allocatable :: quality_flag_2d(:,:)
    real*4, allocatable :: chan_freq(:)
    character*1, allocatable :: chan_pol(:)
    integer :: i, j, k
    
    ! Call the main routine
    call get_wsf_data_with_flags(filename, &
        tb_lowres, lat, lon, land_frac_low, quality_flag_2d, &
        earth_inc_angle, snow, precip, &
        nscans, nfovs, nchans, chan_freq, chan_pol, ierr)
    
    if (ierr == 0 .and. allocated(quality_flag_2d)) then
        ! Expand 2D quality flag to 3D for compatibility
        allocate(quality_flag(nchans, nscans, nfovs))
        do k = 1, nchans
            do j = 1, nfovs
                do i = 1, nscans
                    quality_flag(k,i,j) = quality_flag_2d(i,j)
                end do
            end do
        end do
        deallocate(quality_flag_2d)
    end if
    
    ! Cleanup
    if (allocated(snow)) deallocate(snow)
    if (allocated(precip)) deallocate(precip)
    if (allocated(chan_freq)) deallocate(chan_freq)
    if (allocated(chan_pol)) deallocate(chan_pol)
    
    END SUBROUTINE get_wsf_data

END MODULE TOOLSUBS_WSF