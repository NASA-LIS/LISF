!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: WSF_ARFS_RESAMPLE_HOURLY
!
! DESCRIPTION: Process multiple WSF files from the same hour
!              Stitches data together on ARFS grid
!              Uses mean for overlapping observations
!
!-------------------------------------------------------------------------

subroutine WSF_ARFS_RESAMPLE_HOURLY(hour_files, n_files, output_dir, &
                                    yyyymmdd, hour_str, n)

    USE TOOLSUBS_WSF
    USE invdist_wsf2arfs
    USE LDT_logMod
    USE LDT_coreMod, only: LDT_rc
    USE LDT_WSF_ARFS_netcdfMod, only: LDT_WSF_ARFS_write_netcdf_hourly
    USE LDT_wsf_oplMod, only: wsf_file_info

    IMPLICIT NONE
    
    ! Arguments
    integer, intent(in) :: n_files
    type(wsf_file_info), intent(in) :: hour_files(n_files)
    character(len=*), intent(in) :: output_dir
    character(len=8), intent(in) :: yyyymmdd
    character(len=2), intent(in) :: hour_str
    integer, intent(in) :: n
    
    ! ARFS grid arrays (accumulated across multiple files)
    real*8, allocatable :: ARFS_LAT(:), ARFS_LON(:)
    
    ! Accumulated data arrays
    real*8, allocatable :: ARFS_TIME_SUM(:,:)
    real*4, allocatable :: ARFS_TB_10H_SUM(:,:), ARFS_TB_10V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_18H_SUM(:,:), ARFS_TB_18V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_23H_SUM(:,:), ARFS_TB_23V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_36H_SUM(:,:), ARFS_TB_36V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_89H_SUM(:,:), ARFS_TB_89V_SUM(:,:)
    real*4, allocatable :: ARFS_LAND_FRAC_SUM(:,:)
    
    ! Count arrays for averaging
    integer*4, allocatable :: ARFS_COUNT_10H(:,:), ARFS_COUNT_10V(:,:)
    integer*4, allocatable :: ARFS_COUNT_18H(:,:), ARFS_COUNT_18V(:,:)
    integer*4, allocatable :: ARFS_COUNT_23H(:,:), ARFS_COUNT_23V(:,:)
    integer*4, allocatable :: ARFS_COUNT_36H(:,:), ARFS_COUNT_36V(:,:)
    integer*4, allocatable :: ARFS_COUNT_89H(:,:), ARFS_COUNT_89V(:,:)
    integer*4, allocatable :: ARFS_COUNT_LAND(:,:), ARFS_COUNT_TIME(:,:)
    
    ! Final averaged arrays
    real*8, allocatable :: ARFS_TIME(:,:)
    real*4, allocatable :: ARFS_TB_10H(:,:), ARFS_TB_10V(:,:)
    real*4, allocatable :: ARFS_TB_18H(:,:), ARFS_TB_18V(:,:)
    real*4, allocatable :: ARFS_TB_23H(:,:), ARFS_TB_23V(:,:)
    real*4, allocatable :: ARFS_TB_36H(:,:), ARFS_TB_36V(:,:)
    real*4, allocatable :: ARFS_TB_89H(:,:), ARFS_TB_89V(:,:)
    real*4, allocatable :: ARFS_LAND_FRAC(:,:)
    integer*1, allocatable :: ARFS_QUALITY_FLAG(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_V(:,:), ARFS_SAMPLE_H(:,:)
    
    ! Temporary arrays for single file processing
    real*8, allocatable :: TEMP_TIME(:,:)
    real*4, allocatable :: TEMP_TB_10H(:,:), TEMP_TB_10V(:,:)
    real*4, allocatable :: TEMP_TB_18H(:,:), TEMP_TB_18V(:,:)
    real*4, allocatable :: TEMP_TB_23H(:,:), TEMP_TB_23V(:,:)
    real*4, allocatable :: TEMP_TB_36H(:,:), TEMP_TB_36V(:,:)
    real*4, allocatable :: TEMP_TB_89H(:,:), TEMP_TB_89V(:,:)
    real*4, allocatable :: TEMP_LAND_FRAC(:,:)
    integer*1, allocatable :: TEMP_QUALITY_FLAG(:,:)
    integer*4, allocatable :: TEMP_SAMPLE_V(:,:), TEMP_SAMPLE_H(:,:)
    
    ! WSF input data
    real*4, allocatable :: tb_lowres(:,:,:)
    real*4, allocatable :: lat_in(:,:)
    real*4, allocatable :: lon_in(:,:)
    real*4, allocatable :: land_frac_low(:,:)
    integer*1, allocatable :: quality_flag_in(:,:)
    real*4, allocatable :: earth_inc_angle(:,:,:)
    integer*4, allocatable :: snow_in(:,:)
    integer*4, allocatable :: precip_in(:,:)
    real*4, allocatable :: tb_10h(:,:), tb_10v(:,:)
    real*4, allocatable :: tb_18h(:,:), tb_18v(:,:)
    real*4, allocatable :: tb_23h(:,:), tb_23v(:,:)
    real*4, allocatable :: tb_36h(:,:), tb_36v(:,:)
    real*4, allocatable :: tb_89h(:,:), tb_89v(:,:)
    real*8, allocatable :: time_array(:)
    
    real*4, allocatable :: chan_frequencies(:)
    character*1, allocatable :: chan_polarizations(:)
    
    integer :: i, j, r, c, ifile
    integer :: nscans, nfovs, nchans, ierr, ichan
    real :: freq
    character*1 :: pol
    character(len=255) :: output_filename
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] WSF HOURLY GROUP PROCESSING'
    write(LDT_logunit,*)'[INFO] Hour: ', hour_str, 'H'
    write(LDT_logunit,*)'[INFO] Number of files: ', n_files
    write(LDT_logunit,*)'[INFO] Date: ', yyyymmdd
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Setup ARFS grid
    allocate(ARFS_LAT(1920))
    allocate(ARFS_LON(2560))
    
    do r = 1, 1920
        ARFS_LAT(r) = 90.0 - (r - 0.5) * 0.09375
    end do
    
    do c = 1, 2560
        ARFS_LON(c) = -180.0 + (c - 0.5) * 0.140625
    end do
    
    ! Allocate accumulation arrays
    allocate(ARFS_TIME_SUM(2560,1920))
    allocate(ARFS_TB_10H_SUM(2560,1920))
    allocate(ARFS_TB_10V_SUM(2560,1920))
    allocate(ARFS_TB_18H_SUM(2560,1920))
    allocate(ARFS_TB_18V_SUM(2560,1920))
    allocate(ARFS_TB_23H_SUM(2560,1920))
    allocate(ARFS_TB_23V_SUM(2560,1920))
    allocate(ARFS_TB_36H_SUM(2560,1920))
    allocate(ARFS_TB_36V_SUM(2560,1920))
    allocate(ARFS_TB_89H_SUM(2560,1920))
    allocate(ARFS_TB_89V_SUM(2560,1920))
    allocate(ARFS_LAND_FRAC_SUM(2560,1920))
    
    ! Allocate count arrays
    allocate(ARFS_COUNT_TIME(2560,1920))
    allocate(ARFS_COUNT_10H(2560,1920))
    allocate(ARFS_COUNT_10V(2560,1920))
    allocate(ARFS_COUNT_18H(2560,1920))
    allocate(ARFS_COUNT_18V(2560,1920))
    allocate(ARFS_COUNT_23H(2560,1920))
    allocate(ARFS_COUNT_23V(2560,1920))
    allocate(ARFS_COUNT_36H(2560,1920))
    allocate(ARFS_COUNT_36V(2560,1920))
    allocate(ARFS_COUNT_89H(2560,1920))
    allocate(ARFS_COUNT_89V(2560,1920))
    allocate(ARFS_COUNT_LAND(2560,1920))
    
    ! Allocate temporary arrays
    allocate(TEMP_TIME(2560,1920))
    allocate(TEMP_TB_10H(2560,1920))
    allocate(TEMP_TB_10V(2560,1920))
    allocate(TEMP_TB_18H(2560,1920))
    allocate(TEMP_TB_18V(2560,1920))
    allocate(TEMP_TB_23H(2560,1920))
    allocate(TEMP_TB_23V(2560,1920))
    allocate(TEMP_TB_36H(2560,1920))
    allocate(TEMP_TB_36V(2560,1920))
    allocate(TEMP_TB_89H(2560,1920))
    allocate(TEMP_TB_89V(2560,1920))
    allocate(TEMP_LAND_FRAC(2560,1920))
    allocate(TEMP_QUALITY_FLAG(2560,1920))
    allocate(TEMP_SAMPLE_V(2560,1920))
    allocate(TEMP_SAMPLE_H(2560,1920))
    
    ! Initialize accumulation arrays to zero
    ARFS_TIME_SUM = 0.0
    ARFS_TB_10H_SUM = 0.0
    ARFS_TB_10V_SUM = 0.0
    ARFS_TB_18H_SUM = 0.0
    ARFS_TB_18V_SUM = 0.0
    ARFS_TB_23H_SUM = 0.0
    ARFS_TB_23V_SUM = 0.0
    ARFS_TB_36H_SUM = 0.0
    ARFS_TB_36V_SUM = 0.0
    ARFS_TB_89H_SUM = 0.0
    ARFS_TB_89V_SUM = 0.0
    ARFS_LAND_FRAC_SUM = 0.0
    
    ARFS_COUNT_TIME = 0
    ARFS_COUNT_10H = 0
    ARFS_COUNT_10V = 0
    ARFS_COUNT_18H = 0
    ARFS_COUNT_18V = 0
    ARFS_COUNT_23H = 0
    ARFS_COUNT_23V = 0
    ARFS_COUNT_36H = 0
    ARFS_COUNT_36V = 0
    ARFS_COUNT_89H = 0
    ARFS_COUNT_89V = 0
    ARFS_COUNT_LAND = 0
    
    ! =====================================================================
    ! PROCESS EACH FILE IN THE HOUR GROUP
    ! =====================================================================
    do ifile = 1, n_files
        write(LDT_logunit,*)'[INFO] Processing file ', ifile, '/', n_files, ':'
        write(LDT_logunit,*)'[INFO] ', trim(hour_files(ifile)%filename)
        
        ! Initialize temporary arrays
        TEMP_TIME = 0.0
        TEMP_TB_10H = 0.0
        TEMP_TB_10V = 0.0
        TEMP_TB_18H = 0.0
        TEMP_TB_18V = 0.0
        TEMP_TB_23H = 0.0
        TEMP_TB_23V = 0.0
        TEMP_TB_36H = 0.0
        TEMP_TB_36V = 0.0
        TEMP_TB_89H = 0.0
        TEMP_TB_89V = 0.0
        TEMP_LAND_FRAC = 0.0
        TEMP_QUALITY_FLAG = 0
        TEMP_SAMPLE_V = 0
        TEMP_SAMPLE_H = 0
        
        ! Read WSF data from file
        call get_wsf_data_with_flags(hour_files(ifile)%filename, &
            tb_lowres, lat_in, lon_in, land_frac_low, quality_flag_in, &
            earth_inc_angle, snow_in, precip_in, &
            chan_frequencies, chan_polarizations, &
            nscans, nfovs, nchans, ierr)
        
        if (ierr /= 0) then
            write(LDT_logunit,*)'[WARN] Failed to read file, skipping'
            cycle
        endif
        
        ! Extract channels
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
        allocate(time_array(nscans))
        
        ! Extract each channel
        do ichan = 1, nchans
            freq = chan_frequencies(ichan)
            pol = chan_polarizations(ichan)
            
            if (abs(freq - 10.65) < 1.0) then
                if (pol == 'v') tb_10v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'h') tb_10h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 18.7) < 1.0) then
                if (pol == 'v') tb_18v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'h') tb_18h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 23.8) < 1.0) then
                if (pol == 'v') tb_23v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'h') tb_23h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 36.5) < 1.0) then
                if (pol == 'v') tb_36v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'h') tb_36h(:,:) = tb_lowres(:,:,ichan)
            else if (abs(freq - 89.0) < 1.0) then
                if (pol == 'v') tb_89v(:,:) = tb_lowres(:,:,ichan)
                if (pol == 'h') tb_89h(:,:) = tb_lowres(:,:,ichan)
            end if
        end do
        
        time_array = 0.0
        
        ! Resample to ARFS grid
        call WSF2ARFS_INVDIS(time_array, &
            tb_10h, tb_10v, tb_18h, tb_18v, &
            tb_23h, tb_23v, tb_36h, tb_36v, &
            tb_89h, tb_89v, land_frac_low, &
            snow_in, precip_in, quality_flag_in, &
            lat_in, lon_in, nfovs, nscans, &
            ARFS_LAT, ARFS_LON, &
            TEMP_TIME, TEMP_LAND_FRAC, &
            TEMP_TB_10H, TEMP_TB_10V, &
            TEMP_TB_18H, TEMP_TB_18V, &
            TEMP_TB_23H, TEMP_TB_23V, &
            TEMP_TB_36H, TEMP_TB_36V, &
            TEMP_TB_89H, TEMP_TB_89V, &
            TEMP_QUALITY_FLAG, &
            TEMP_SAMPLE_V, TEMP_SAMPLE_H)
        
        ! ACCUMULATE DATA (only non-zero values)
        do r = 1, 1920
            do c = 1, 2560
                ! TB_10H
                if (TEMP_TB_10H(c,r) > 0.0) then
                    ARFS_TB_10H_SUM(c,r) = ARFS_TB_10H_SUM(c,r) + TEMP_TB_10H(c,r)
                    ARFS_COUNT_10H(c,r) = ARFS_COUNT_10H(c,r) + 1
                endif
                
                ! TB_10V
                if (TEMP_TB_10V(c,r) > 0.0) then
                    ARFS_TB_10V_SUM(c,r) = ARFS_TB_10V_SUM(c,r) + TEMP_TB_10V(c,r)
                    ARFS_COUNT_10V(c,r) = ARFS_COUNT_10V(c,r) + 1
                endif
                
                ! TB_18H
                if (TEMP_TB_18H(c,r) > 0.0) then
                    ARFS_TB_18H_SUM(c,r) = ARFS_TB_18H_SUM(c,r) + TEMP_TB_18H(c,r)
                    ARFS_COUNT_18H(c,r) = ARFS_COUNT_18H(c,r) + 1
                endif
                
                ! TB_18V
                if (TEMP_TB_18V(c,r) > 0.0) then
                    ARFS_TB_18V_SUM(c,r) = ARFS_TB_18V_SUM(c,r) + TEMP_TB_18V(c,r)
                    ARFS_COUNT_18V(c,r) = ARFS_COUNT_18V(c,r) + 1
                endif
                
                ! TB_23H
                if (TEMP_TB_23H(c,r) > 0.0) then
                    ARFS_TB_23H_SUM(c,r) = ARFS_TB_23H_SUM(c,r) + TEMP_TB_23H(c,r)
                    ARFS_COUNT_23H(c,r) = ARFS_COUNT_23H(c,r) + 1
                endif
                
                ! TB_23V
                if (TEMP_TB_23V(c,r) > 0.0) then
                    ARFS_TB_23V_SUM(c,r) = ARFS_TB_23V_SUM(c,r) + TEMP_TB_23V(c,r)
                    ARFS_COUNT_23V(c,r) = ARFS_COUNT_23V(c,r) + 1
                endif
                
                ! TB_36H
                if (TEMP_TB_36H(c,r) > 0.0) then
                    ARFS_TB_36H_SUM(c,r) = ARFS_TB_36H_SUM(c,r) + TEMP_TB_36H(c,r)
                    ARFS_COUNT_36H(c,r) = ARFS_COUNT_36H(c,r) + 1
                endif
                
                ! TB_36V
                if (TEMP_TB_36V(c,r) > 0.0) then
                    ARFS_TB_36V_SUM(c,r) = ARFS_TB_36V_SUM(c,r) + TEMP_TB_36V(c,r)
                    ARFS_COUNT_36V(c,r) = ARFS_COUNT_36V(c,r) + 1
                endif
                
                ! TB_89H
                if (TEMP_TB_89H(c,r) > 0.0) then
                    ARFS_TB_89H_SUM(c,r) = ARFS_TB_89H_SUM(c,r) + TEMP_TB_89H(c,r)
                    ARFS_COUNT_89H(c,r) = ARFS_COUNT_89H(c,r) + 1
                endif
                
                ! TB_89V
                if (TEMP_TB_89V(c,r) > 0.0) then
                    ARFS_TB_89V_SUM(c,r) = ARFS_TB_89V_SUM(c,r) + TEMP_TB_89V(c,r)
                    ARFS_COUNT_89V(c,r) = ARFS_COUNT_89V(c,r) + 1
                endif
                
                ! Land fraction
                if (TEMP_LAND_FRAC(c,r) > 0.0) then
                    ARFS_LAND_FRAC_SUM(c,r) = ARFS_LAND_FRAC_SUM(c,r) + TEMP_LAND_FRAC(c,r)
                    ARFS_COUNT_LAND(c,r) = ARFS_COUNT_LAND(c,r) + 1
                endif
            end do
        end do
        
        ! Cleanup for this file
        deallocate(tb_lowres, lat_in, lon_in)
        deallocate(land_frac_low, quality_flag_in)
        deallocate(earth_inc_angle, snow_in, precip_in)
        deallocate(chan_frequencies, chan_polarizations)
        deallocate(tb_10h, tb_10v, tb_18h, tb_18v)
        deallocate(tb_23h, tb_23v, tb_36h, tb_36v)
        deallocate(tb_89h, tb_89v, time_array)
        
    end do ! ifile loop
    
    ! =====================================================================
    ! CALCULATE MEANS FOR OVERLAPPING REGIONS
    ! =====================================================================
    
    ! Allocate final arrays
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
    
    ! Initialize to zero/missing
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
    
    ! Calculate means
    do r = 1, 1920
        do c = 1, 2560
            ! TB_10H
            if (ARFS_COUNT_10H(c,r) > 0) then
                ARFS_TB_10H(c,r) = ARFS_TB_10H_SUM(c,r) / real(ARFS_COUNT_10H(c,r))
                ARFS_SAMPLE_H(c,r) = ARFS_COUNT_10H(c,r)
            endif
            
            ! TB_10V
            if (ARFS_COUNT_10V(c,r) > 0) then
                ARFS_TB_10V(c,r) = ARFS_TB_10V_SUM(c,r) / real(ARFS_COUNT_10V(c,r))
                ARFS_SAMPLE_V(c,r) = ARFS_COUNT_10V(c,r)
            endif
            
            ! TB_18H
            if (ARFS_COUNT_18H(c,r) > 0) then
                ARFS_TB_18H(c,r) = ARFS_TB_18H_SUM(c,r) / real(ARFS_COUNT_18H(c,r))
            endif
            
            ! TB_18V
            if (ARFS_COUNT_18V(c,r) > 0) then
                ARFS_TB_18V(c,r) = ARFS_TB_18V_SUM(c,r) / real(ARFS_COUNT_18V(c,r))
            endif
            
            ! TB_23H
            if (ARFS_COUNT_23H(c,r) > 0) then
                ARFS_TB_23H(c,r) = ARFS_TB_23H_SUM(c,r) / real(ARFS_COUNT_23H(c,r))
            endif
            
            ! TB_23V
            if (ARFS_COUNT_23V(c,r) > 0) then
                ARFS_TB_23V(c,r) = ARFS_TB_23V_SUM(c,r) / real(ARFS_COUNT_23V(c,r))
            endif
            
            ! TB_36H
            if (ARFS_COUNT_36H(c,r) > 0) then
                ARFS_TB_36H(c,r) = ARFS_TB_36H_SUM(c,r) / real(ARFS_COUNT_36H(c,r))
            endif
            
            ! TB_36V
            if (ARFS_COUNT_36V(c,r) > 0) then
                ARFS_TB_36V(c,r) = ARFS_TB_36V_SUM(c,r) / real(ARFS_COUNT_36V(c,r))
            endif
            
            ! TB_89H
            if (ARFS_COUNT_89H(c,r) > 0) then
                ARFS_TB_89H(c,r) = ARFS_TB_89H_SUM(c,r) / real(ARFS_COUNT_89H(c,r))
            endif
            
            ! TB_89V
            if (ARFS_COUNT_89V(c,r) > 0) then
                ARFS_TB_89V(c,r) = ARFS_TB_89V_SUM(c,r) / real(ARFS_COUNT_89V(c,r))
            endif
            
            ! Land fraction
            if (ARFS_COUNT_LAND(c,r) > 0) then
                ARFS_LAND_FRAC(c,r) = ARFS_LAND_FRAC_SUM(c,r) / real(ARFS_COUNT_LAND(c,r))
            endif
        end do
    end do
    
    ! =====================================================================
    ! WRITE OUTPUT TO NETCDF FILE
    ! =====================================================================
    
    ! Construct output filename for hourly stitched data
    output_filename = trim(output_dir)//'/WSF_stitched_'//yyyymmdd//'_'//hour_str//'00.nc'
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Writing stitched hourly output'
    write(LDT_logunit,*)'[INFO] Output file: ', trim(output_filename)
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Call the NetCDF writing subroutine
    call LDT_WSF_ARFS_write_netcdf_hourly(2560, 1920, &
        ARFS_TB_10H, ARFS_TB_10V, ARFS_TB_18H, ARFS_TB_18V, &
        ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, &
        ARFS_TB_89H, ARFS_TB_89V, ARFS_LAND_FRAC, &
        ARFS_QUALITY_FLAG, ARFS_SAMPLE_V, ARFS_SAMPLE_H, &
        ARFS_LAT, ARFS_LON, output_filename, &
        yyyymmdd, hour_str, n_files)
    
    ! Report statistics
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Hourly stitching statistics:'
    write(LDT_logunit,*)'[INFO] Files processed: ', n_files
    write(LDT_logunit,*)'[INFO] Grid points with data:'
    write(LDT_logunit,*)'[INFO]   10H: ', count(ARFS_TB_10H > 0.0)
    write(LDT_logunit,*)'[INFO]   10V: ', count(ARFS_TB_10V > 0.0)
    write(LDT_logunit,*)'[INFO] Max samples per pixel:'
    write(LDT_logunit,*)'[INFO]   H-pol: ', maxval(ARFS_SAMPLE_H)
    write(LDT_logunit,*)'[INFO]   V-pol: ', maxval(ARFS_SAMPLE_V)
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Cleanup
    deallocate(ARFS_LAT, ARFS_LON)
    deallocate(ARFS_TIME_SUM, ARFS_TB_10H_SUM, ARFS_TB_10V_SUM)
    deallocate(ARFS_TB_18H_SUM, ARFS_TB_18V_SUM)
    deallocate(ARFS_TB_23H_SUM, ARFS_TB_23V_SUM)
    deallocate(ARFS_TB_36H_SUM, ARFS_TB_36V_SUM)
    deallocate(ARFS_TB_89H_SUM, ARFS_TB_89V_SUM)
    deallocate(ARFS_LAND_FRAC_SUM)
    
    deallocate(ARFS_COUNT_TIME, ARFS_COUNT_10H, ARFS_COUNT_10V)
    deallocate(ARFS_COUNT_18H, ARFS_COUNT_18V)
    deallocate(ARFS_COUNT_23H, ARFS_COUNT_23V)
    deallocate(ARFS_COUNT_36H, ARFS_COUNT_36V)
    deallocate(ARFS_COUNT_89H, ARFS_COUNT_89V)
    deallocate(ARFS_COUNT_LAND)
    
    deallocate(TEMP_TIME, TEMP_TB_10H, TEMP_TB_10V)
    deallocate(TEMP_TB_18H, TEMP_TB_18V)
    deallocate(TEMP_TB_23H, TEMP_TB_23V)
    deallocate(TEMP_TB_36H, TEMP_TB_36V)
    deallocate(TEMP_TB_89H, TEMP_TB_89V)
    deallocate(TEMP_LAND_FRAC, TEMP_QUALITY_FLAG)
    deallocate(TEMP_SAMPLE_V, TEMP_SAMPLE_H)
    
    deallocate(ARFS_TIME, ARFS_TB_10H, ARFS_TB_10V)
    deallocate(ARFS_TB_18H, ARFS_TB_18V)
    deallocate(ARFS_TB_23H, ARFS_TB_23V)
    deallocate(ARFS_TB_36H, ARFS_TB_36V)
    deallocate(ARFS_TB_89H, ARFS_TB_89V)
    deallocate(ARFS_LAND_FRAC, ARFS_QUALITY_FLAG)
    deallocate(ARFS_SAMPLE_V, ARFS_SAMPLE_H)

end subroutine WSF_ARFS_RESAMPLE_HOURLY