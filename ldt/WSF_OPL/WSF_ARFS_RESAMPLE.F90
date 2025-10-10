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
! SUBROUTINE: WSF_ARFS_RESAMPLE
!
! DESCRIPTION: Resample WSF TB to Air Force Grid
!              Modified to output AMSR_OPL-style with quality and land filtering
!
!-------------------------------------------------------------------------

subroutine WSF_ARFS_RESAMPLE(wsf_filelist, nfiles, output_dir, n)

    USE VARIABLES
    USE DATADOMAIN
    USE FUNCTIONS
    USE TOOLSUBS_WSF
    USE invdist_wsf2arfs
    USE LDT_logMod
    USE LDT_wsf_oplMod
    USE LDT_coreMod, only: LDT_rc
    USE netcdf
    
    IMPLICIT NONE
    
    ! Arguments
    character(len=*), intent(in) :: wsf_filelist(:)
    integer, intent(in) :: nfiles
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: n
    
    ! Local variables
    integer :: ierr, nscans, nfovs, nchans, ifile
    real*4, allocatable :: tb_lowres(:,:,:)
    real*4, allocatable :: lat_in(:,:), lon_in(:,:)
    real*4, allocatable :: land_frac_low(:,:)
    integer*1, allocatable :: quality_flag_in(:,:,:)
    real*4, allocatable :: earth_inc_angle(:,:,:)
    
    real*8, allocatable :: ARFS_LAT(:), ARFS_LON(:)
    real*4, allocatable :: ARFS_TB(:,:,:)  ! Temporary 3D array
    real*4, allocatable :: ARFS_LAND_FRAC(:,:)
    integer*1, allocatable :: ARFS_QUALITY_FLAG(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_COUNT(:,:)
    
    ! Channel mapping arrays
    real*4, allocatable :: chan_frequencies(:)
    character*1, allocatable :: chan_polarizations(:)
    
    integer :: ichan, c, r, band_for_qc
    real :: freq
    character*1 :: pol
    integer*1 :: qc_bits
    real, parameter :: LAND_FRAC_THRESHOLD = 0.8
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Starting WSF ARFS Resampling'
    write(LDT_logunit,*)'[INFO] Processing ', nfiles, ' files'
    write(LDT_logunit,*)'[INFO] Using LOW resolution (30 km) data'
    write(LDT_logunit,*)'[INFO] Output: AMSR_OPL style (10 channels only)'
    write(LDT_logunit,*)'[INFO] Quality filter: bits 0-4 = 0'
    write(LDT_logunit,*)'[INFO] Land fraction filter: > ', LAND_FRAC_THRESHOLD
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Setup ARFS grid once
    write(LDT_logunit,*)'[INFO] Setting up ARFS grid...'
    call ARFS_GEO
    allocate(ARFS_LAT(arfs_nrow_lat))
    allocate(ARFS_LON(arfs_mcol_lon))
    ARFS_LAT = LAT(arfs_geo_lat_lo, arfs_geo_lat_up, -arfs_lat_space)
    ARFS_LON = LON(arfs_geo_lon_lf, arfs_geo_lon_rt, arfs_lon_space)
    
    write(LDT_logunit,*)'[INFO] ARFS grid dimensions: ', &
        arfs_mcol_lon, ' x ', arfs_nrow_lat
    
    ! Process each file
    do ifile = 1, nfiles
        write(LDT_logunit,*) ''
        write(LDT_logunit,*) '[INFO] ========================================'
        write(LDT_logunit,*) '[INFO] Processing file ', ifile, ' of ', nfiles
        write(LDT_logunit,*) '[INFO] File: ', trim(wsf_filelist(ifile))
        
        ! Read WSF data
        write(LDT_logunit,*)'[INFO] Reading WSF NetCDF file...'
        call get_wsf_data(wsf_filelist(ifile), &
            tb_lowres, lat_in, lon_in, land_frac_low, &
            quality_flag_in, earth_inc_angle, &
            nscans, nfovs, nchans, ierr)
        
        if (ierr /= 0) then
            write(LDT_logunit,*)'[ERR] Failed to read WSF data'
            write(LDT_logunit,*)'[ERR] Skipping file: ', trim(wsf_filelist(ifile))
            cycle
        end if
        
        write(LDT_logunit,*)'[INFO] Input dimensions: nscans=', nscans, &
            ' nfovs=', nfovs, ' nchans=', nchans
        
        ! Convert longitude from -180:180 to 0:360 if needed
        write(LDT_logunit,*)'[INFO] Converting longitude range...'
        do c = 1, nfovs
            do r = 1, nscans
                if (lon_in(r,c) < 0.0) then
                    lon_in(r,c) = lon_in(r,c) + 360.0
                end if
            end do
        end do
        
        ! Read channel frequencies and polarizations
        allocate(chan_frequencies(nchans))
        allocate(chan_polarizations(nchans))
        
        ierr = nf90_open(trim(wsf_filelist(ifile)), NF90_NOWRITE, ncid)
        if (ierr == NF90_NOERR) then
            ierr = nf90_inq_varid(ncid, 'ChanFrequency', varid)
            if (ierr == NF90_NOERR) then
                ierr = nf90_get_var(ncid, varid, chan_frequencies)
                write(LDT_logunit,*)'[INFO] Read channel frequencies'
            end if
            
            ierr = nf90_inq_varid(ncid, 'ChanPolarization', varid)
            if (ierr == NF90_NOERR) then
                ierr = nf90_get_var(ncid, varid, chan_polarizations)
                write(LDT_logunit,*)'[INFO] Read channel polarizations'
            end if
            
            ierr = nf90_close(ncid)
        end if
        
        ! Allocate temporary arrays
        allocate(ARFS_TB(nchans, arfs_mcol_lon, arfs_nrow_lat))
        allocate(ARFS_LAND_FRAC(arfs_mcol_lon, arfs_nrow_lat))
        allocate(ARFS_QUALITY_FLAG(arfs_mcol_lon, arfs_nrow_lat))
        allocate(ARFS_SAMPLE_COUNT(arfs_mcol_lon, arfs_nrow_lat))
        
        ! Perform inverse distance resampling with quality filtering
        write(LDT_logunit,*)'[INFO] Performing inverse distance resampling with QC...'
        call WSF2ARFS_INVDIS_WITH_QC(tb_lowres, lat_in, lon_in, &
            quality_flag_in, land_frac_low, &
            nscans, nfovs, nchans, &
            ARFS_LAT, ARFS_LON, &
            ARFS_TB, ARFS_QUALITY_FLAG, ARFS_LAND_FRAC, &
            ARFS_SAMPLE_COUNT, LAND_FRAC_THRESHOLD)
        
        ! Split 3D array into separate 2D arrays for V and H polarizations only
        write(LDT_logunit,*)'[INFO] Extracting V and H polarization channels...'
        
        do ichan = 1, nchans
            freq = chan_frequencies(ichan)
            pol = chan_polarizations(ichan)
            
            ! ONLY process V and H polarizations (no Stokes parameters)
            ! Map to AMSR_OPL structure
            
            ! 10.85 GHz
            if (abs(freq - 10.85) < 0.1) then
                if (pol == 'v' .or. pol == 'V') then
                    WSFopl%ARFS_TB_10V(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_10V'
                else if (pol == 'h' .or. pol == 'H') then
                    WSFopl%ARFS_TB_10H(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_10H'
                end if
            
            ! 18.7 GHz
            else if (abs(freq - 18.7) < 0.1) then
                if (pol == 'v' .or. pol == 'V') then
                    WSFopl%ARFS_TB_18V(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_18V'
                else if (pol == 'h' .or. pol == 'H') then
                    WSFopl%ARFS_TB_18H(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_18H'
                end if
            
            ! 23.8 GHz
            else if (abs(freq - 23.8) < 0.1) then
                if (pol == 'v' .or. pol == 'V') then
                    WSFopl%ARFS_TB_23V(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_23V'
                else if (pol == 'h' .or. pol == 'H') then
                    WSFopl%ARFS_TB_23H(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_23H'
                end if
            
            ! 36.5 GHz
            else if (abs(freq - 36.5) < 0.1) then
                if (pol == 'v' .or. pol == 'V') then
                    WSFopl%ARFS_TB_36V(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_36V'
                else if (pol == 'h' .or. pol == 'H') then
                    WSFopl%ARFS_TB_36H(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_36H'
                end if
            
            ! 89.0 GHz
            else if (abs(freq - 89.0) < 0.1) then
                if (pol == 'v' .or. pol == 'V') then
                    WSFopl%ARFS_TB_89V(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_89V'
                else if (pol == 'h' .or. pol == 'H') then
                    WSFopl%ARFS_TB_89H(:,:) = ARFS_TB(ichan,:,:)
                    write(LDT_logunit,*)'[INFO] Mapped channel ', ichan, ' to TB_89H'
                end if
            end if
            
            ! Note: Stokes parameters (pol='3' or '4') are intentionally skipped
        end do
        
        ! Apply final land fraction and quality flag masking
        write(LDT_logunit,*)'[INFO] Applying land fraction and quality masks...'
        do c = 1, arfs_mcol_lon
            do r = 1, arfs_nrow_lat
                ! Mask water pixels (land fraction <= threshold)
                if (ARFS_LAND_FRAC(c,r) > 0.0 .and. &
                    ARFS_LAND_FRAC(c,r) <= LAND_FRAC_THRESHOLD) then
                    WSFopl%ARFS_TB_10V(c,r) = -9999.0
                    WSFopl%ARFS_TB_10H(c,r) = -9999.0
                    WSFopl%ARFS_TB_18V(c,r) = -9999.0
                    WSFopl%ARFS_TB_18H(c,r) = -9999.0
                    WSFopl%ARFS_TB_23V(c,r) = -9999.0
                    WSFopl%ARFS_TB_23H(c,r) = -9999.0
                    WSFopl%ARFS_TB_36V(c,r) = -9999.0
                    WSFopl%ARFS_TB_36H(c,r) = -9999.0
                    WSFopl%ARFS_TB_89V(c,r) = -9999.0
                    WSFopl%ARFS_TB_89H(c,r) = -9999.0
                end if
                
                ! Check quality flag (bits 0-4 must all be 0)
                ! Bit 0 = IsNotValid
                ! Bit 1 = IsExclusionCondition
                ! Bit 2 = IsDegradationCondition
                ! Bit 3 = IsLimitedUtilityValidationCondition
                ! Bit 4 = IsOverlap
                qc_bits = IAND(ARFS_QUALITY_FLAG(c,r), 31)  ! 31 = 00011111 in binary
                if (qc_bits /= 0) then
                    WSFopl%ARFS_TB_10V(c,r) = -9999.0
                    WSFopl%ARFS_TB_10H(c,r) = -9999.0
                    WSFopl%ARFS_TB_18V(c,r) = -9999.0
                    WSFopl%ARFS_TB_18H(c,r) = -9999.0
                    WSFopl%ARFS_TB_23V(c,r) = -9999.0
                    WSFopl%ARFS_TB_23H(c,r) = -9999.0
                    WSFopl%ARFS_TB_36V(c,r) = -9999.0
                    WSFopl%ARFS_TB_36H(c,r) = -9999.0
                    WSFopl%ARFS_TB_89V(c,r) = -9999.0
                    WSFopl%ARFS_TB_89H(c,r) = -9999.0
                end if
            end do
        end do
        
        ! Copy auxiliary data
        WSFopl%ARFS_LAND_FRAC(:,:) = ARFS_LAND_FRAC(:,:)
        WSFopl%ARFS_QUALITY_FLAG(:,:) = ARFS_QUALITY_FLAG(:,:)
        
        write(LDT_logunit,*)'[INFO] Successfully processed file'
        
        ! Cleanup temporary arrays
        deallocate(tb_lowres)
        deallocate(lat_in, lon_in)
        deallocate(land_frac_low)
        deallocate(quality_flag_in, earth_inc_angle)
        deallocate(ARFS_TB, ARFS_LAND_FRAC)
        deallocate(ARFS_QUALITY_FLAG, ARFS_SAMPLE_COUNT)
        deallocate(chan_frequencies, chan_polarizations)
        
        write(LDT_logunit,*) '[INFO] Completed file ', ifile, ' of ', nfiles
        write(LDT_logunit,*) '[INFO] ========================================'
    end do
    
    ! Cleanup grid arrays
    deallocate(ARFS_LAT, ARFS_LON)
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] All WSF Files Processed'
    write(LDT_logunit,*)'[INFO] Output stored in WSFopl module structure'
    write(LDT_logunit,*)'[INFO] ========================================='

end subroutine WSF_ARFS_RESAMPLE