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
!              CORRECTED: Added snow/precip filtering to match AMSR_OPL
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
    integer*1, allocatable :: quality_flag_in(:,:)
    real*4, allocatable :: earth_inc_angle(:,:,:)
    
    ! CORRECTED: Added snow/precip arrays
    integer*4, allocatable :: snow_in(:,:), precip_in(:,:)
    
    ! Time array
    real*8, allocatable :: time_array(:)
    
    ! Output grid arrays
    real*8, allocatable :: ARFS_LAT(:), ARFS_LON(:)
    real*8, allocatable :: ARFS_TIME(:,:)
    real*4, allocatable :: ARFS_TB_10H(:,:), ARFS_TB_10V(:,:)
    real*4, allocatable :: ARFS_TB_18H(:,:), ARFS_TB_18V(:,:)
    real*4, allocatable :: ARFS_TB_23H(:,:), ARFS_TB_23V(:,:)
    real*4, allocatable :: ARFS_TB_36H(:,:), ARFS_TB_36V(:,:)
    real*4, allocatable :: ARFS_TB_89H(:,:), ARFS_TB_89V(:,:)
    real*4, allocatable :: ARFS_LAND_FRAC(:,:)
    integer*1, allocatable :: ARFS_QUALITY_FLAG(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_V(:,:), ARFS_SAMPLE_H(:,:)
    
    ! Individual channel arrays for extraction
    real*4, allocatable :: tb_10h(:,:), tb_10v(:,:)
    real*4, allocatable :: tb_18h(:,:), tb_18v(:,:)
    real*4, allocatable :: tb_23h(:,:), tb_23v(:,:)
    real*4, allocatable :: tb_36h(:,:), tb_36v(:,:)
    real*4, allocatable :: tb_89h(:,:), tb_89v(:,:)
    
    ! Channel mapping arrays
    real*4, allocatable :: chan_frequencies(:)
    character*1, allocatable :: chan_polarizations(:)
    
    integer :: ichan, c, r
    real :: freq
    character*1 :: pol
    integer*1 :: qc_bits
    real, parameter :: LAND_FRAC_THRESHOLD = 0.8
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Starting WSF ARFS Resampling (CORRECTED)'
    write(LDT_logunit,*)'[INFO] Processing ', nfiles, ' files'
    write(LDT_logunit,*)'[INFO] Using LOW resolution (30 km) data'
    write(LDT_logunit,*)'[INFO] Output: AMSR_OPL style (10 channels)'
    write(LDT_logunit,*)'[INFO] Snow/Precip filtering: ENABLED'
    write(LDT_logunit,*)'[INFO] Quality filter: bits 0-4 = 0'
    write(LDT_logunit,*)'[INFO] Land fraction filter: > ', LAND_FRAC_THRESHOLD
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Setup ARFS grid
    allocate(ARFS_LAT(1920))
    allocate(ARFS_LON(2560))
    
    ! Create lat/lon arrays for ARFS grid
    do r = 1, 1920
        ARFS_LAT(r) = 90.0 - (r - 0.5) * 0.09375
    end do
    
    do c = 1, 2560
        ARFS_LON(c) = -180.0 + (c - 0.5) * 0.140625
    end do
    
    ! Allocate output arrays
    allocate(ARFS_TIME(2560,1920))
    allocate(ARFS_TB_10H(2560,1920))
    allocate(ARFS_TB_10V(2560,1920))
    allocate(ARFS_TB_18H(2560,1920))
    allocate(ARFS_TB_18V(2560,1920))
    allocate(ARFS_TB_23H(2560,1920))
    allocate(ARFS_TB_23V(2560,1920))
    allocate(ARFS_TB_36H(2560,1920))
    allocate(ARFS_TB_36V(2560,1920))
    allocate(ARFS_TB_89H(2560,1920))
    allocate(ARFS_TB_89V(2560,1920))
    allocate(ARFS_LAND_FRAC(2560,1920))
    allocate(ARFS_QUALITY_FLAG(2560,1920))
    allocate(ARFS_SAMPLE_V(2560,1920))
    allocate(ARFS_SAMPLE_H(2560,1920))
    
    ! Process each WSF file
    do ifile = 1, nfiles
        write(LDT_logunit,*) '[INFO] ========================================'
        write(LDT_logunit,*) '[INFO] Processing file ', ifile, ' of ', nfiles
        write(LDT_logunit,*) '[INFO] File: ', trim(wsf_filelist(ifile))
        
        ! Read WSF data with quality flags
        call get_wsf_data_with_flags(wsf_filelist(ifile), &
            tb_lowres, lat_in, lon_in, land_frac_low, quality_flag_in, &
            earth_inc_angle, snow_in, precip_in, &
            nscans, nfovs, nchans, chan_frequencies, chan_polarizations, ierr)
        
        if (ierr /= 0) then
            write(LDT_logunit,*)'[WARN] Failed to read file, skipping'
            cycle
        end if
        
        write(LDT_logunit,*)'[INFO] File dimensions: ', nscans, 'x', nfovs, 'x', nchans
        write(LDT_logunit,*)'[INFO] Extracting channels for resampling'
        
        ! Allocate channel arrays
        allocate(tb_10h(nscans,nfovs))
        allocate(tb_10v(nscans,nfovs))
        allocate(tb_18h(nscans,nfovs))
        allocate(tb_18v(nscans,nfovs))
        allocate(tb_23h(nscans,nfovs))
        allocate(tb_23v(nscans,nfovs))
        allocate(tb_36h(nscans,nfovs))
        allocate(tb_36v(nscans,nfovs))
        allocate(tb_89h(nscans,nfovs))
        allocate(tb_89v(nscans,nfovs))
        
        ! Initialize to missing
        tb_10h = -9999.0
        tb_10v = -9999.0
        tb_18h = -9999.0
        tb_18v = -9999.0
        tb_23h = -9999.0
        tb_23v = -9999.0
        tb_36h = -9999.0
        tb_36v = -9999.0
        tb_89h = -9999.0
        tb_89v = -9999.0
        
        ! Extract channels based on frequency and polarization
        do ichan = 1, nchans
            freq = chan_frequencies(ichan)
            pol = chan_polarizations(ichan)
            
            ! Map to standard channels (approximate frequencies)
            if (abs(freq - 10.65) < 0.5) then
                if (pol == 'V') tb_10v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'H') tb_10h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 18.7) < 0.5) then
                if (pol == 'V') tb_18v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'H') tb_18h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 23.8) < 0.5) then
                if (pol == 'V') tb_23v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'H') tb_23h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 36.5) < 1.0) then
                if (pol == 'V') tb_36v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'H') tb_36h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 89.0) < 2.0) then
                if (pol == 'V') tb_89v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'H') tb_89h(:,:) = tb_lowres(:,:,ichan)
            end if
        end do
        
        ! Create time array (simplified - using dummy values)
        allocate(time_array(nscans))
        time_array = 0.0  ! Would be populated from actual file metadata
        
        write(LDT_logunit,*)'[INFO] Calling inverse distance resampling WITH snow/precip filtering'
        
        ! Call resampling routine WITH SNOW/PRECIP FILTERING
        call WSF2ARFS_INVDIS(time_array, &
            tb_10h, tb_10v, tb_18h, tb_18v, &
            tb_23h, tb_23v, tb_36h, tb_36v, &
            tb_89h, tb_89v, land_frac_low, &
            snow_in, precip_in, quality_flag_in, &  ! CORRECTED: Pass snow/precip flags
            lat_in, lon_in, nscans, nfovs, &
            ARFS_LAT, ARFS_LON, &
            ARFS_TIME, ARFS_LAND_FRAC, &
            ARFS_TB_10H, ARFS_TB_10V, ARFS_TB_18H, ARFS_TB_18V, &
            ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, &
            ARFS_TB_89H, ARFS_TB_89V, ARFS_QUALITY_FLAG, &
            ARFS_SAMPLE_V, ARFS_SAMPLE_H)
        
        write(LDT_logunit,*)'[INFO] Resampling complete, applying quality filters'
        
        ! Apply additional quality filtering (matching AMSR_OPL)
        do c = 1, 2560
            do r = 1, 1920
                ! Check land fraction threshold
                if (ARFS_LAND_FRAC(c,r) >= 0 .and. ARFS_LAND_FRAC(c,r) < LAND_FRAC_THRESHOLD) then
                    ARFS_TB_10V(c,r) = -9999.0
                    ARFS_TB_10H(c,r) = -9999.0
                    ARFS_TB_18V(c,r) = -9999.0
                    ARFS_TB_18H(c,r) = -9999.0
                    ARFS_TB_23V(c,r) = -9999.0
                    ARFS_TB_23H(c,r) = -9999.0
                    ARFS_TB_36V(c,r) = -9999.0
                    ARFS_TB_36H(c,r) = -9999.0
                    ARFS_TB_89V(c,r) = -9999.0
                    ARFS_TB_89H(c,r) = -9999.0
                end if
                
                ! Check quality flag (bits 0-4 must all be 0)
                qc_bits = IAND(ARFS_QUALITY_FLAG(c,r), 31)  ! 31 = 00011111 in binary
                if (qc_bits /= 0) then
                    ARFS_TB_10V(c,r) = -9999.0
                    ARFS_TB_10H(c,r) = -9999.0
                    ARFS_TB_18V(c,r) = -9999.0
                    ARFS_TB_18H(c,r) = -9999.0
                    ARFS_TB_23V(c,r) = -9999.0
                    ARFS_TB_23H(c,r) = -9999.0
                    ARFS_TB_36V(c,r) = -9999.0
                    ARFS_TB_36H(c,r) = -9999.0
                    ARFS_TB_89V(c,r) = -9999.0
                    ARFS_TB_89H(c,r) = -9999.0
                end if
            end do
        end do
        
        ! Copy to module arrays
        WSFopl%ARFS_TB_10V(:,:) = ARFS_TB_10V(:,:)
        WSFopl%ARFS_TB_10H(:,:) = ARFS_TB_10H(:,:)
        WSFopl%ARFS_TB_18V(:,:) = ARFS_TB_18V(:,:)
        WSFopl%ARFS_TB_18H(:,:) = ARFS_TB_18H(:,:)
        WSFopl%ARFS_TB_23V(:,:) = ARFS_TB_23V(:,:)
        WSFopl%ARFS_TB_23H(:,:) = ARFS_TB_23H(:,:)
        WSFopl%ARFS_TB_36V(:,:) = ARFS_TB_36V(:,:)
        WSFopl%ARFS_TB_36H(:,:) = ARFS_TB_36H(:,:)
        WSFopl%ARFS_TB_89V(:,:) = ARFS_TB_89V(:,:)
        WSFopl%ARFS_TB_89H(:,:) = ARFS_TB_89H(:,:)
        WSFopl%ARFS_LAND_FRAC(:,:) = ARFS_LAND_FRAC(:,:)
        WSFopl%ARFS_QUALITY_FLAG(:,:) = ARFS_QUALITY_FLAG(:,:)
        WSFopl%ARFS_SAMPLE_V(:,:) = ARFS_SAMPLE_V(:,:)
        WSFopl%ARFS_SAMPLE_H(:,:) = ARFS_SAMPLE_H(:,:)
        
        ! Store snow/precip flags if allocated in module
        if (allocated(WSFopl%ARFS_SNOW)) then
            ! Compute from quality flag bit pattern
            do c = 1, 2560
                do r = 1, 1920
                    if (IBITS(ARFS_QUALITY_FLAG(c,r), 2, 1) == 1) then
                        WSFopl%ARFS_SNOW(c,r) = 1
                    else
                        WSFopl%ARFS_SNOW(c,r) = 0
                    endif
                    
                    if (IBITS(ARFS_QUALITY_FLAG(c,r), 1, 1) == 1) then
                        WSFopl%ARFS_PRECIP(c,r) = 1
                    else
                        WSFopl%ARFS_PRECIP(c,r) = 0
                    endif
                end do
            end do
        endif
        
        write(LDT_logunit,*)'[INFO] Successfully processed file'
        
        ! Cleanup temporary arrays
        deallocate(tb_lowres)
        deallocate(lat_in, lon_in)
        deallocate(land_frac_low)
        deallocate(quality_flag_in, earth_inc_angle)
        deallocate(snow_in, precip_in)
        deallocate(time_array)
        deallocate(chan_frequencies, chan_polarizations)
        deallocate(tb_10h, tb_10v, tb_18h, tb_18v)
        deallocate(tb_23h, tb_23v, tb_36h, tb_36v)
        deallocate(tb_89h, tb_89v)
        
        write(LDT_logunit,*) '[INFO] Completed file ', ifile, ' of ', nfiles
        write(LDT_logunit,*) '[INFO] ========================================'
    end do
    
    ! Cleanup grid arrays
    deallocate(ARFS_LAT, ARFS_LON)
    deallocate(ARFS_TIME)
    deallocate(ARFS_TB_10H, ARFS_TB_10V)
    deallocate(ARFS_TB_18H, ARFS_TB_18V)
    deallocate(ARFS_TB_23H, ARFS_TB_23V)
    deallocate(ARFS_TB_36H, ARFS_TB_36V)
    deallocate(ARFS_TB_89H, ARFS_TB_89V)
    deallocate(ARFS_LAND_FRAC)
    deallocate(ARFS_QUALITY_FLAG)
    deallocate(ARFS_SAMPLE_V, ARFS_SAMPLE_H)
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] All WSF Files Processed'
    write(LDT_logunit,*)'[INFO] Snow/Precip filtering applied successfully'
    write(LDT_logunit,*)'[INFO] Output stored in WSFopl module structure'
    write(LDT_logunit,*)'[INFO] ========================================='

end subroutine WSF_ARFS_RESAMPLE