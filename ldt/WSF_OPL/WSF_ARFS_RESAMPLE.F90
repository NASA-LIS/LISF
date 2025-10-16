!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: WSF_ARFS_RESAMPLE
!
! DESCRIPTION: Resample WSF data to ARFS grid with ALL channels
!              FIXED: Corrected dimension order AND subroutine call
!              READS: All 17 channels from WSF, extracts 10 AMSR2-style
!
!-------------------------------------------------------------------------

subroutine WSF_ARFS_RESAMPLE(wsf_filename, output_dir, n)

    USE TOOLSUBS_WSF
    USE invdist_wsf2arfs
    USE LDT_logMod
    USE LDT_coreMod, only: LDT_rc
    USE LDT_WSF_ARFS_netcdfMod, only: LDT_WSF_ARFS_write_netcdf

    IMPLICIT NONE
    
    ! Arguments
    character(len=255), intent(in) :: wsf_filename
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: n
    
    ! WSF input data (from get_wsf_data_with_flags)
    ! NOTE: These are in FORTRAN ORDER (reversed from ncdump)
    real*4, allocatable :: tb_lowres(:,:,:)      ! (nFOVR, nScanR, nChan)
    real*4, allocatable :: lat_in(:,:)           ! (nFOVR, nScanR)
    real*4, allocatable :: lon_in(:,:)           ! (nFOVR, nScanR)
    real*4, allocatable :: land_frac_low(:,:)    ! (nFOVR, nScanR)
    integer*1, allocatable :: quality_flag_in(:,:) ! (nFOVR, nScanR)
    real*4, allocatable :: earth_inc_angle(:,:,:)  ! (nFOVR, nScanR, 1)
    integer*4, allocatable :: snow_in(:,:)       ! (nFOVR, nScanR)
    integer*4, allocatable :: precip_in(:,:)     ! (nFOVR, nScanR)
    
    ! Individual channel arrays (Fortran order: nFOVR, nScanR)
    real*4, allocatable :: tb_10h(:,:), tb_10v(:,:)
    real*4, allocatable :: tb_18h(:,:), tb_18v(:,:)
    real*4, allocatable :: tb_23h(:,:), tb_23v(:,:)
    real*4, allocatable :: tb_36h(:,:), tb_36v(:,:)
    real*4, allocatable :: tb_89h(:,:), tb_89v(:,:)
    
    ! Time array
    real*8, allocatable :: time_array(:)
    
    ! ARFS output grid arrays
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
    
    ! Channel mapping arrays
    real*4, allocatable :: chan_frequencies(:)
    character*1, allocatable :: chan_polarizations(:)
    
    integer :: nscans, nfovs, nchans
    integer :: ierr, ichan, c, r, ifov, iscan
    real :: freq
    character*1 :: pol
    character(len=255) :: output_filename
    character(len=255) :: basename
    integer :: filename_start_pos
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] WSF ARFS Resampling - ALL CHANNELS'
    write(LDT_logunit,*)'[INFO] File: ', trim(wsf_filename)
    write(LDT_logunit,*)'[INFO] Output: AMSR_OPL style (10 channels from 17 total)'
    write(LDT_logunit,*)'[INFO] Snow/Precip filtering: ENABLED'
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Setup ARFS grid (0.09375 deg lat x 0.140625 deg lon)
    allocate(ARFS_LAT(1920))
    allocate(ARFS_LON(2560))
    
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
    
    ! Initialize to missing values
    ARFS_TIME = 0.0
    ARFS_TB_10H = 0.0
    ARFS_TB_10V = 0.0
    ARFS_TB_18H = 0.0
    ARFS_TB_18V = 0.0
    ARFS_TB_23H = 0.0
    ARFS_TB_23V = 0.0
    ARFS_TB_36H = 0.0
    ARFS_TB_36V = 0.0
    ARFS_TB_89H = 0.0
    ARFS_TB_89V = 0.0
    ARFS_LAND_FRAC = 0.0
    ARFS_QUALITY_FLAG = 0
    ARFS_SAMPLE_V = 0
    ARFS_SAMPLE_H = 0
    
    ! =====================================================================
    ! READ WSF DATA - ALL 17 CHANNELS IN FORTRAN ORDER
    ! =====================================================================
    write(LDT_logunit,*)'[INFO] Reading WSF NetCDF data...'
    
    call get_wsf_data_with_flags(trim(wsf_filename), &
        tb_lowres, lat_in, lon_in, land_frac_low, quality_flag_in, &
        earth_inc_angle, snow_in, precip_in, &
        nscans, nfovs, nchans, chan_frequencies, chan_polarizations, ierr)
    
    if (ierr /= 0) then
        write(LDT_logunit,*)'[WARN] Failed to read file, skipping'
        return
    end if
    
    write(LDT_logunit,*)'[INFO] Dimensions (Fortran order):'
    write(LDT_logunit,*)'[INFO]   nFOVR  = ', nfovs
    write(LDT_logunit,*)'[INFO]   nScanR = ', nscans
    write(LDT_logunit,*)'[INFO]   nChan  = ', nchans
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Extracting 10 AMSR2-style channels...'
    
    ! =====================================================================
    ! ALLOCATE CHANNEL ARRAYS - FORTRAN ORDER (nFOVR, nScanR)
    ! =====================================================================
    allocate(tb_10h(nfovs, nscans))
    allocate(tb_10v(nfovs, nscans))
    allocate(tb_18h(nfovs, nscans))
    allocate(tb_18v(nfovs, nscans))
    allocate(tb_23h(nfovs, nscans))
    allocate(tb_23v(nfovs, nscans))
    allocate(tb_36h(nfovs, nscans))
    allocate(tb_36v(nfovs, nscans))
    allocate(tb_89h(nfovs, nscans))
    allocate(tb_89v(nfovs, nscans))
    
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
    
    ! =====================================================================
    ! EXTRACT CHANNELS FROM tb_lowres(nFOVR, nScanR, nChan)
    ! =====================================================================
    do ichan = 1, nchans
        freq = chan_frequencies(ichan)
        pol = chan_polarizations(ichan)
        
        ! 10.65/10.85 GHz channels
        if (abs(freq - 10.65) < 0.5 .or. abs(freq - 10.85) < 0.5) then
            if (pol == 'V' .or. pol == 'v') then
                tb_10v(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_10V from channel', ichan
            else if (pol == 'H' .or. pol == 'h') then
                tb_10h(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_10H from channel', ichan
            end if
            
        ! 18.7 GHz channels
        else if (abs(freq - 18.7) < 0.5) then
            if (pol == 'V' .or. pol == 'v') then
                tb_18v(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_18V from channel', ichan
            else if (pol == 'H' .or. pol == 'h') then
                tb_18h(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_18H from channel', ichan
            end if
            
        ! 23.8 GHz channels
        else if (abs(freq - 23.8) < 0.5) then
            if (pol == 'V' .or. pol == 'v') then
                tb_23v(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_23V from channel', ichan
            else if (pol == 'H' .or. pol == 'h') then
                tb_23h(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_23H from channel', ichan
            end if
            
        ! 36.5 GHz channels
        else if (abs(freq - 36.5) < 1.0) then
            if (pol == 'V' .or. pol == 'v') then
                tb_36v(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_36V from channel', ichan
            else if (pol == 'H' .or. pol == 'h') then
                tb_36h(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_36H from channel', ichan
            end if
            
        ! 89.0 GHz channels
        else if (abs(freq - 89.0) < 2.0) then
            if (pol == 'V' .or. pol == 'v') then
                tb_89v(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_89V from channel', ichan
            else if (pol == 'H' .or. pol == 'h') then
                tb_89h(:,:) = tb_lowres(:,:,ichan)
                write(LDT_logunit,*)'[INFO] ✓ TB_89H from channel', ichan
            end if
        end if
    end do
    
    ! Create time array
    allocate(time_array(nscans))
    time_array = 0.0
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Calling inverse distance resampling...'
    
    ! =====================================================================
    ! CALL RESAMPLING - MATCH ACTUAL SUBROUTINE SIGNATURE
    ! From invdist_wsf2arfs.F90:
    ! WSF2ARFS_INVDIS(tim, tb_10h, tb_10v, tb_18h, tb_18v,
    !                 tb_23h, tb_23v, tb_36h, tb_36v, tb_89h, tb_89v, 
    !                 land_water_frac, snow_flag, precip_flag, quality_flag,
    !                 lat_in, lon_in, nscans_in, nfovs_in,
    !                 ref_lat, ref_lon, arfs_time, arfs_land_water_frac,
    !                 arfs_tb_10h, arfs_tb_10v, arfs_tb_18h, arfs_tb_18v,
    !                 arfs_tb_23h, arfs_tb_23v, arfs_tb_36h, arfs_tb_36v,
    !                 arfs_tb_89h, arfs_tb_89v, arfs_quality_flag,
    !                 arfs_sample_v, arfs_sample_h)
    ! =====================================================================
    call WSF2ARFS_INVDIS(time_array, &
        tb_10h, tb_10v, tb_18h, tb_18v, &
        tb_23h, tb_23v, tb_36h, tb_36v, &
        tb_89h, tb_89v, land_frac_low, &
        snow_in, precip_in, quality_flag_in, &
        lat_in, lon_in, nfovs, nscans, &
        ARFS_LAT, ARFS_LON, &
        ARFS_TIME, ARFS_LAND_FRAC, &
        ARFS_TB_10H, ARFS_TB_10V, &
        ARFS_TB_18H, ARFS_TB_18V, &
        ARFS_TB_23H, ARFS_TB_23V, &
        ARFS_TB_36H, ARFS_TB_36V, &
        ARFS_TB_89H, ARFS_TB_89V, &
        ARFS_QUALITY_FLAG, &
        ARFS_SAMPLE_V, ARFS_SAMPLE_H)

    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] ✓ Resampling complete'
    write(LDT_logunit,*)'[INFO] Writing NetCDF output file...'
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! =====================================================================
    ! WRITE OUTPUT TO NETCDF FILE
    ! =====================================================================
    
    ! Extract basename from input filename to construct output filename
    filename_start_pos = index(wsf_filename, '/', back=.true.) + 1
    basename = wsf_filename(filename_start_pos:)
    
    ! Construct output filename: WSF_resampled_YYYYMMDD_tHHMMSS.nc
    ! Input format: WSFM_01_d20250601_t002900_e004059_gCOOK-B_r05899_c02_i0_v0522_res_sdr.nc
    ! Extract: d20250601 (positions 9-16) and t002900 (positions 18-24)
    output_filename = trim(output_dir)//'/WSF_resampled_'// &
                     basename(9:17)//'_'//basename(19:25)//'.nc'
    
    write(LDT_logunit,*)'[INFO] Calling NetCDF writer...'
    
    ! Call the NetCDF writing subroutine from the module
    call LDT_WSF_ARFS_write_netcdf(2560, 1920, &
        ARFS_TB_10H, ARFS_TB_10V, ARFS_TB_18H, ARFS_TB_18V, &
        ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, &
        ARFS_TB_89H, ARFS_TB_89V, ARFS_LAND_FRAC, &
        ARFS_QUALITY_FLAG, ARFS_SAMPLE_V, ARFS_SAMPLE_H, &
        ARFS_LAT, ARFS_LON, output_filename)

    
    ! Cleanup
    deallocate(ARFS_LAT, ARFS_LON)
    deallocate(ARFS_TIME)
    deallocate(ARFS_TB_10H, ARFS_TB_10V)
    deallocate(ARFS_TB_18H, ARFS_TB_18V)
    deallocate(ARFS_TB_23H, ARFS_TB_23V)
    deallocate(ARFS_TB_36H, ARFS_TB_36V)
    deallocate(ARFS_TB_89H, ARFS_TB_89V)
    deallocate(ARFS_LAND_FRAC, ARFS_QUALITY_FLAG)
    deallocate(ARFS_SAMPLE_V, ARFS_SAMPLE_H)
    
    deallocate(tb_lowres, lat_in, lon_in)
    deallocate(land_frac_low, quality_flag_in)
    deallocate(earth_inc_angle, snow_in, precip_in)
    deallocate(chan_frequencies, chan_polarizations)
    
    deallocate(tb_10h, tb_10v)
    deallocate(tb_18h, tb_18v)
    deallocate(tb_23h, tb_23v)
    deallocate(tb_36h, tb_36v)
    deallocate(tb_89h, tb_89v)
    deallocate(time_array)

end subroutine WSF_ARFS_RESAMPLE
