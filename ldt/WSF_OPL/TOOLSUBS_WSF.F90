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
        tb_lowres, lat, lon, land_frac_low, &
        snow, precip, quality_flag, &
        nscans, nfovs, nchans, ierr)
    
    ! Arguments  
    character(*), intent(in) :: filename
    real*4, allocatable, intent(out) :: tb_lowres(:,:,:)   ! (nchans, nscans, nfovs)
    real*4, allocatable, intent(out) :: lat(:,:)           ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: lon(:,:)           ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: land_frac_low(:,:) ! (nscans, nfovs)
    integer*4, allocatable, intent(out) :: snow(:,:)       ! (nscans, nfovs)
    integer*4, allocatable, intent(out) :: precip(:,:)     ! (nscans, nfovs)
    integer*1, allocatable, intent(out) :: quality_flag(:,:) ! (nscans, nfovs)
    integer, intent(out) :: nscans, nfovs, nchans
    integer, intent(out) :: ierr
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    ! Local variables
    integer :: ncid, varid
    integer :: nscanr_dimid, nfovr_dimid, nchan_dimid
    logical :: file_exists
    integer :: i, j, ichan
    real :: sil, tt18
    real*4, allocatable :: tb_18v(:,:), tb_18h(:,:)
    real*4, allocatable :: tb_23v(:,:), tb_36v(:,:), tb_89v(:,:)
    
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
        call nf90_close(ncid)
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
    
    ! Allocate arrays
    allocate(tb_lowres(nchans, nscans, nfovs))
    allocate(lat(nscans, nfovs))
    allocate(lon(nscans, nfovs))
    allocate(land_frac_low(nscans, nfovs))
    allocate(snow(nscans, nfovs))
    allocate(precip(nscans, nfovs))
    allocate(quality_flag(nscans, nfovs))
    
    ! Allocate temporary arrays for snow/precip detection
    allocate(tb_18v(nscans, nfovs))
    allocate(tb_18h(nscans, nfovs))
    allocate(tb_23v(nscans, nfovs))
    allocate(tb_36v(nscans, nfovs))
    allocate(tb_89v(nscans, nfovs))
    
    ! Read latitude
    ierr = nf90_inq_varid(ncid, 'Latitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lat)
        write(LDT_logunit,*)'[INFO] Read Latitude'
    end if
    
    ! Read longitude
    ierr = nf90_inq_varid(ncid, 'Longitude', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, lon)
        write(LDT_logunit,*)'[INFO] Read Longitude'
    end if
    
    ! Read TbLowRes
    ierr = nf90_inq_varid(ncid, 'TbLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, tb_lowres)
        write(LDT_logunit,*)'[INFO] Read TbLowRes'
    end if
    
    ! Read LandFractionLowRes
    ierr = nf90_inq_varid(ncid, 'LandFractionLowRes', varid)
    if (ierr == NF90_NOERR) then
        ierr = nf90_get_var(ncid, varid, land_frac_low)
        write(LDT_logunit,*)'[INFO] Read LandFractionLowRes'
    end if
    
    ! Close file for now
    ierr = nf90_close(ncid)
    
    ! ==================================================================
    ! EXTRACT SPECIFIC CHANNELS FOR SNOW/PRECIP DETECTION
    ! Match AMSR_OPL processing exactly
    ! ==================================================================
    write(LDT_logunit,*)'[INFO] Extracting channels for snow/precip detection'
    
    ! Initialize
    tb_18v = 0.0
    tb_18h = 0.0
    tb_23v = 0.0
    tb_36v = 0.0
    tb_89v = 0.0
    
    ! Extract channels based on frequency
    ! Need to read channel info from file
    ierr = nf90_open(trim(filename), NF90_NOWRITE, ncid)
    
    ! Read channel frequencies and polarizations to identify channels
    real*4, allocatable :: chan_freq(:)
    character*1, allocatable :: chan_pol(:)
    allocate(chan_freq(nchans))
    allocate(chan_pol(nchans))
    
    ierr = nf90_inq_varid(ncid, 'ChanFrequency', varid)
    if (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, varid, chan_freq)
    
    ierr = nf90_inq_varid(ncid, 'ChanPolarization', varid)
    if (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, varid, chan_pol)
    
    ierr = nf90_close(ncid)
    
    ! Extract specific channels needed for snow/precip algorithm
    do ichan = 1, nchans
        ! 18.7 GHz V
        if (abs(chan_freq(ichan) - 18.7) < 0.1 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_18v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_18V at channel ', ichan
        end if
        
        ! 18.7 GHz H
        if (abs(chan_freq(ichan) - 18.7) < 0.1 .and. &
            (chan_pol(ichan) == 'h' .or. chan_pol(ichan) == 'H')) then
            tb_18h(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_18H at channel ', ichan
        end if
        
        ! 23.8 GHz V
        if (abs(chan_freq(ichan) - 23.8) < 0.1 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_23v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_23V at channel ', ichan
        end if
        
        ! 36.5 GHz V
        if (abs(chan_freq(ichan) - 36.5) < 0.1 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_36v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_36V at channel ', ichan
        end if
        
        ! 89.0 GHz V
        if (abs(chan_freq(ichan) - 89.0) < 0.1 .and. &
            (chan_pol(ichan) == 'v' .or. chan_pol(ichan) == 'V')) then
            tb_89v(:,:) = tb_lowres(ichan,:,:)
            write(LDT_logunit,*)'[INFO] Found TB_89V at channel ', ichan
        end if
    end do
    
    deallocate(chan_freq, chan_pol)
    
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

END MODULE TOOLSUBS_WSF