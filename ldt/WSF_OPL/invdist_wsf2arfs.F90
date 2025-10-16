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
! MODULE: invdist_wsf2arfs
!
! DESCRIPTION: Resample WSF TB onto Air Force Grid
!              CORRECTED: Added snow/precip filtering to match AMSR_OPL exactly
!
!-------------------------------------------------------------------------

MODULE invdist_wsf2arfs
    USE LDT_logMod, only: LDT_logunit
    IMPLICIT NONE

CONTAINS

    SUBROUTINE WSF2ARFS_INVDIS(tim, tb_10h, tb_10v, tb_18h, tb_18v, &
        tb_23h, tb_23v, tb_36h, tb_36v, tb_89h, tb_89v, land_water_frac, &
        snow_flag, precip_flag, quality_flag, &
        lat_in, lon_in, nscans_in, nfovs_in, &
        ref_lat, ref_lon, &
        arfs_time, arfs_land_water_frac, &
        arfs_tb_10h, arfs_tb_10v, arfs_tb_18h, arfs_tb_18v, &
        arfs_tb_23h, arfs_tb_23v, arfs_tb_36h, arfs_tb_36v, &
        arfs_tb_89h, arfs_tb_89v, arfs_quality_flag, &
        arfs_sample_v, arfs_sample_h)
    
    ! Arguments - EXACT match to AMSR_OPL signature
    integer, intent(in) :: nscans_in, nfovs_in
    real*8, intent(in) :: tim(nscans_in)
    real*4, intent(in) :: tb_10h(nscans_in, nfovs_in), tb_10v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_18h(nscans_in, nfovs_in), tb_18v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_23h(nscans_in, nfovs_in), tb_23v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_36h(nscans_in, nfovs_in), tb_36v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_89h(nscans_in, nfovs_in), tb_89v(nscans_in, nfovs_in)
    real*4, intent(in) :: land_water_frac(nscans_in, nfovs_in)
    integer*4, intent(in) :: snow_flag(nscans_in, nfovs_in), precip_flag(nscans_in, nfovs_in)
    integer*1, intent(in) :: quality_flag(nscans_in, nfovs_in)
    real*4, intent(in) :: lat_in(nscans_in, nfovs_in), lon_in(nscans_in, nfovs_in)
    real*8, intent(in) :: ref_lat(:), ref_lon(:)
    
    ! Outputs
    real*8, intent(out) :: arfs_time(2560,1920)
    real*4, intent(out) :: arfs_tb_10h(2560,1920), arfs_tb_10v(2560,1920)
    real*4, intent(out) :: arfs_tb_18h(2560,1920), arfs_tb_18v(2560,1920)
    real*4, intent(out) :: arfs_tb_23h(2560,1920), arfs_tb_23v(2560,1920)
    real*4, intent(out) :: arfs_tb_36h(2560,1920), arfs_tb_36v(2560,1920)
    real*4, intent(out) :: arfs_tb_89h(2560,1920), arfs_tb_89v(2560,1920)
    real*4, intent(out) :: arfs_land_water_frac(2560,1920)
    integer*1, intent(out) :: arfs_quality_flag(2560,1920)
    integer*4, intent(out) :: arfs_sample_v(2560,1920), arfs_sample_h(2560,1920)
    
    ! Local variables - EXACT match to AMSR_OPL
    integer :: ii, jj, r, c, rr, cc, rmin, rmax, cmin, cmax, i, j
    integer, allocatable :: zerodistflag(:,:)
    real*8, parameter :: RE_KM = 6371.228
    real*8, parameter :: search_radius = 20.0
    real*8, parameter :: PI = 3.141592653589793238
    real*8, parameter :: d2r = PI/180.0
    real*8 :: gcdist, lat1, lon1, lat2, lon2
    logical :: has_snow, has_precip
    
    ! Weight arrays
    real*8, allocatable :: arfs_wt_tim(:,:)
    real*4, allocatable :: arfs_wt_tb10h(:,:), arfs_wt_tb10v(:,:)
    real*4, allocatable :: arfs_wt_tb18h(:,:), arfs_wt_tb18v(:,:)
    real*4, allocatable :: arfs_wt_tb23h(:,:), arfs_wt_tb23v(:,:)
    real*4, allocatable :: arfs_wt_tb36h(:,:), arfs_wt_tb36v(:,:)
    real*4, allocatable :: arfs_wt_tb89h(:,:), arfs_wt_tb89v(:,:)
    real*4, allocatable :: arfs_wt_land_water_frac(:,:)
    
    integer, allocatable :: snow_count(:,:), precip_count(:,:)
    integer, allocatable :: ocean_count(:,:), total_count(:,:)
    integer, allocatable :: excluded_snow_count(:,:),excluded_precip_count(:,:)
    
    write(LDT_logunit,*)'[INFO] Starting WSF inverse distance resampling'
    ! Allocate weight arrays
    allocate(arfs_wt_tim(2560,1920))
    allocate(arfs_wt_tb10h(2560,1920))
    allocate(arfs_wt_tb10v(2560,1920))
    allocate(arfs_wt_tb18h(2560,1920))
    allocate(arfs_wt_tb18v(2560,1920))
    allocate(arfs_wt_tb23h(2560,1920))
    allocate(arfs_wt_tb23v(2560,1920))
    allocate(arfs_wt_tb36h(2560,1920))
    allocate(arfs_wt_tb36v(2560,1920))
    allocate(arfs_wt_tb89h(2560,1920))
    allocate(arfs_wt_tb89v(2560,1920))
    allocate(arfs_wt_land_water_frac(2560,1920))
    
    ! Allocate counter arrays
    allocate(snow_count(2560,1920))
    allocate(precip_count(2560,1920))
    allocate(ocean_count(2560,1920))
    allocate(total_count(2560,1920))
    allocate(excluded_snow_count(2560,1920))
    allocate(excluded_precip_count(2560,1920))
    
    write(LDT_logunit,*)'[INFO] WITH Snow/Precip filtering (CORRECTED to match AMSR_OPL)'
    write(LDT_logunit,*)'[INFO] Input dimensions: ', nscans_in, 'x', nfovs_in
    write(LDT_logunit,*)'[INFO] Output dimensions: ', size(ref_lon), 'x', size(ref_lat)
    
    ! Initialize arrays
    allocate(zerodistflag(size(ref_lon), size(ref_lat)))
    zerodistflag = 0
    
    ! Initialize output arrays
    arfs_time = 0.0
    arfs_tb_10h = 0.0
    arfs_tb_10v = 0.0
    arfs_tb_18h = 0.0
    arfs_tb_18v = 0.0
    arfs_tb_23h = 0.0
    arfs_tb_23v = 0.0
    arfs_tb_36h = 0.0
    arfs_tb_36v = 0.0
    arfs_tb_89h = 0.0
    arfs_tb_89v = 0.0
    arfs_land_water_frac = 0.0
    
    ! Initialize weight arrays
    arfs_wt_tim = 0.0
    arfs_wt_tb10h = 0.0
    arfs_wt_tb10v = 0.0
    arfs_wt_tb18h = 0.0
    arfs_wt_tb18v = 0.0
    arfs_wt_tb23h = 0.0
    arfs_wt_tb23v = 0.0
    arfs_wt_tb36h = 0.0
    arfs_wt_tb36v = 0.0
    arfs_wt_tb89h = 0.0
    arfs_wt_tb89v = 0.0
    arfs_wt_land_water_frac = 0.0
    
    ! Initialize quality flag tracking
    arfs_quality_flag = 0
    snow_count = 0
    precip_count = 0
    ocean_count = 0
    total_count = 0
    excluded_snow_count = 0
    excluded_precip_count = 0
    arfs_sample_v = 0
    arfs_sample_h = 0
    
    write(LDT_logunit,*)'[DEBUG] Array dimensions:'
    write(LDT_logunit,*)'   nscans_in, nfovs_in = ', nscans_in, nfovs_in
    write(LDT_logunit,*)'   size(ref_lat), size(ref_lon) = ', size(ref_lat), size(ref_lon)
    
    ! Main resampling loop - CORRECTED with snow/precip filtering
    do ii = 1, nfovs_in
        do jj = 1, nscans_in
            ! Skip invalid coordinates
            if (lat_in(jj,ii) < -90.0 .or. lat_in(jj,ii) > 90.0 .or. &
                lon_in(jj,ii) < -180.0 .or. lon_in(jj,ii) > 180.0) then
                cycle
            endif
            
            ! CRITICAL: Check for snow and precipitation flags
            has_snow = (snow_flag(jj,ii) == 1)
            has_precip = (precip_flag(jj,ii) == 1)
            
            ! Determine search bounds
            lat1 = lat_in(jj,ii)
            lon1 = lon_in(jj,ii)
            
            ! Find grid cell bounds
            c = MINLOC(ABS(lat_in(jj,ii) - ref_lat(:)), 1)  ! Lat direction
            r = MINLOC(ABS(lon_in(jj,ii) - ref_lon(:)), 1)  ! Lon direction
            
            rmin = r - 5
            IF (rmin < 1) rmin = 1
            rmax = r + 5  
            IF (rmax > size(ref_lon)) rmax = size(ref_lon)
            
            cmin = c - 5
            IF (cmin < 1) cmin = 1
            cmax = c + 5
            IF (cmax > size(ref_lat)) cmax = size(ref_lat)
            
            ! Process each grid cell in search radius
            do rr = rmin, rmax
                do cc = cmin, cmax
                    lat2 = ref_lat(cc)
                    lon2 = ref_lon(rr)
                    
                    ! Calculate great circle distance
                    gcdist = 2.0 * RE_KM * ASIN(SQRT((SIN(d2r*(lat1-lat2)/2.0))**2 + &
                             COS(d2r*lat1) * COS(d2r*lat2) * (SIN(d2r*(lon1-lon2)/2.0))**2))
                    
                    if (gcdist < search_radius) then
                        ! Always update counts for quality flag tracking
                        total_count(rr,cc) = total_count(rr,cc) + 1
                        
                        ! Update ocean/snow/precip counts based on quality flag bits
                        if (IBITS(quality_flag(jj,ii), 0, 1) == 1) then
                            ocean_count(rr,cc) = ocean_count(rr,cc) + 1
                        endif
                        if (has_precip) then
                            precip_count(rr,cc) = precip_count(rr,cc) + 1
                        endif
                        if (has_snow) then
                            snow_count(rr,cc) = snow_count(rr,cc) + 1
                        endif
                        
                        !! CRITICAL FILTERING: Skip resampling if snow or precip detected
                        !if (has_snow .OR. has_precip) then
                            !! Track excluded footprints
                            !if (has_snow) excluded_snow_count(rr,cc) = excluded_snow_count(rr,cc) + 1
                            !if (has_precip) excluded_precip_count(rr,cc) = excluded_precip_count(rr,cc) + 1
                            !! SKIP THIS FOOTPRINT - DO NOT RESAMPLE
                            !cycle
                        !endif
                        
                        ! No snow/precip - proceed with resampling
                        if (gcdist < 0.0001D0) then
                            ! Exact match
                            zerodistflag(rr,cc) = 1
                            arfs_quality_flag(rr,cc) = quality_flag(jj,ii)
                            
                            ! Time
                            if (jj <= size(tim)) then
                                if (ABS(tim(jj) - (-9999.0)) > 1.0E-6) then
                                    arfs_time(rr,cc) = tim(jj)
                                    arfs_wt_tim(rr,cc) = 1.0
                                endif
                            endif
                            
                            ! TB channels
                            if (ABS(tb_10h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_10h(rr,cc) = tb_10h(jj,ii)
                                arfs_wt_tb10h(rr,cc) = 1.0
                                arfs_sample_h(rr,cc) = arfs_sample_h(rr,cc) + 1
                            endif
                            if (ABS(tb_10v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_10v(rr,cc) = tb_10v(jj,ii)
                                arfs_wt_tb10v(rr,cc) = 1.0
                                arfs_sample_v(rr,cc) = arfs_sample_v(rr,cc) + 1
                            endif
                            if (ABS(tb_18h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_18h(rr,cc) = tb_18h(jj,ii)
                                arfs_wt_tb18h(rr,cc) = 1.0
                            endif
                            if (ABS(tb_18v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_18v(rr,cc) = tb_18v(jj,ii)
                                arfs_wt_tb18v(rr,cc) = 1.0
                            endif
                            if (ABS(tb_23h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_23h(rr,cc) = tb_23h(jj,ii)
                                arfs_wt_tb23h(rr,cc) = 1.0
                            endif
                            if (ABS(tb_23v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_23v(rr,cc) = tb_23v(jj,ii)
                                arfs_wt_tb23v(rr,cc) = 1.0
                            endif
                            if (ABS(tb_36h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_36h(rr,cc) = tb_36h(jj,ii)
                                arfs_wt_tb36h(rr,cc) = 1.0
                            endif
                            if (ABS(tb_36v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_36v(rr,cc) = tb_36v(jj,ii)
                                arfs_wt_tb36v(rr,cc) = 1.0
                            endif
                            if (ABS(tb_89h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_89h(rr,cc) = tb_89h(jj,ii)
                                arfs_wt_tb89h(rr,cc) = 1.0
                            endif
                            if (ABS(tb_89v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_89v(rr,cc) = tb_89v(jj,ii)
                                arfs_wt_tb89v(rr,cc) = 1.0
                            endif
                            
                            ! Land water fraction
                            if (ABS(land_water_frac(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_land_water_frac(rr,cc) = land_water_frac(jj,ii)
                                arfs_wt_land_water_frac(rr,cc) = 1.0
                            endif
                            
                        else if (zerodistflag(rr,cc) /= 1) then
                            ! Inverse distance weighting
                            
                            ! Time
                            if (jj <= size(tim)) then
                                if (ABS(tim(jj) - (-9999.0)) > 1.0E-6) then
                                    arfs_time(rr,cc) = arfs_time(rr,cc) + tim(jj) / gcdist
                                    arfs_wt_tim(rr,cc) = arfs_wt_tim(rr,cc) + 1.0 / gcdist
                                endif
                            endif
                            
                            ! TB channels
                            if (ABS(tb_10h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_10h(rr,cc) = arfs_tb_10h(rr,cc) + tb_10h(jj,ii) / gcdist
                                arfs_wt_tb10h(rr,cc) = arfs_wt_tb10h(rr,cc) + 1.0 / gcdist
                                arfs_sample_h(rr,cc) = arfs_sample_h(rr,cc) + 1
                            endif
                            if (ABS(tb_10v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_10v(rr,cc) = arfs_tb_10v(rr,cc) + tb_10v(jj,ii) / gcdist
                                arfs_wt_tb10v(rr,cc) = arfs_wt_tb10v(rr,cc) + 1.0 / gcdist
                                arfs_sample_v(rr,cc) = arfs_sample_v(rr,cc) + 1
                            endif
                            if (ABS(tb_18h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_18h(rr,cc) = arfs_tb_18h(rr,cc) + tb_18h(jj,ii) / gcdist
                                arfs_wt_tb18h(rr,cc) = arfs_wt_tb18h(rr,cc) + 1.0 / gcdist
                            endif
                            if (ABS(tb_18v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_18v(rr,cc) = arfs_tb_18v(rr,cc) + tb_18v(jj,ii) / gcdist
                                arfs_wt_tb18v(rr,cc) = arfs_wt_tb18v(rr,cc) + 1.0 / gcdist
                            endif
                            if (ABS(tb_23h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_23h(rr,cc) = arfs_tb_23h(rr,cc) + tb_23h(jj,ii) / gcdist
                                arfs_wt_tb23h(rr,cc) = arfs_wt_tb23h(rr,cc) + 1.0 / gcdist
                            endif
                            if (ABS(tb_23v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_23v(rr,cc) = arfs_tb_23v(rr,cc) + tb_23v(jj,ii) / gcdist
                                arfs_wt_tb23v(rr,cc) = arfs_wt_tb23v(rr,cc) + 1.0 / gcdist
                            endif
                            if (ABS(tb_36h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_36h(rr,cc) = arfs_tb_36h(rr,cc) + tb_36h(jj,ii) / gcdist
                                arfs_wt_tb36h(rr,cc) = arfs_wt_tb36h(rr,cc) + 1.0 / gcdist
                            endif
                            if (ABS(tb_36v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_36v(rr,cc) = arfs_tb_36v(rr,cc) + tb_36v(jj,ii) / gcdist
                                arfs_wt_tb36v(rr,cc) = arfs_wt_tb36v(rr,cc) + 1.0 / gcdist
                            endif
                            if (ABS(tb_89h(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_89h(rr,cc) = arfs_tb_89h(rr,cc) + tb_89h(jj,ii) / gcdist
                                arfs_wt_tb89h(rr,cc) = arfs_wt_tb89h(rr,cc) + 1.0 / gcdist
                            endif
                            if (ABS(tb_89v(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_tb_89v(rr,cc) = arfs_tb_89v(rr,cc) + tb_89v(jj,ii) / gcdist
                                arfs_wt_tb89v(rr,cc) = arfs_wt_tb89v(rr,cc) + 1.0 / gcdist
                            endif
                            
                            ! Land water fraction
                            if (ABS(land_water_frac(jj,ii) - (-9999.0)) > 1.0E-6) then
                                arfs_land_water_frac(rr,cc) = arfs_land_water_frac(rr,cc) + &
                                                                land_water_frac(jj,ii) / gcdist
                                arfs_wt_land_water_frac(rr,cc) = arfs_wt_land_water_frac(rr,cc) + 1.0 / gcdist
                            endif
                        endif
                    endif
                end do
            end do
        end do
    end do
    
    ! Apply weighting and set fill values
    WHERE(arfs_time /= 0.0 .AND. arfs_wt_tim /= 0.0)
        arfs_time = arfs_time / arfs_wt_tim
    ELSEWHERE
        arfs_time = -9999.0
    END WHERE
    
    WHERE(arfs_tb_10h /= 0.0 .AND. arfs_wt_tb10h /= 0.0)
        arfs_tb_10h = arfs_tb_10h / arfs_wt_tb10h
    ELSEWHERE
        arfs_tb_10h = -9999.0
    END WHERE
    
    WHERE(arfs_tb_10v /= 0.0 .AND. arfs_wt_tb10v /= 0.0)
        arfs_tb_10v = arfs_tb_10v / arfs_wt_tb10v
    ELSEWHERE
        arfs_tb_10v = -9999.0
    END WHERE
    
    WHERE(arfs_tb_18h /= 0.0 .AND. arfs_wt_tb18h /= 0.0)
        arfs_tb_18h = arfs_tb_18h / arfs_wt_tb18h
    ELSEWHERE
        arfs_tb_18h = -9999.0
    END WHERE
    
    WHERE(arfs_tb_18v /= 0.0 .AND. arfs_wt_tb18v /= 0.0)
        arfs_tb_18v = arfs_tb_18v / arfs_wt_tb18v
    ELSEWHERE
        arfs_tb_18v = -9999.0
    END WHERE
    
    WHERE(arfs_tb_23h /= 0.0 .AND. arfs_wt_tb23h /= 0.0)
        arfs_tb_23h = arfs_tb_23h / arfs_wt_tb23h
    ELSEWHERE
        arfs_tb_23h = -9999.0
    END WHERE
    
    WHERE(arfs_tb_23v /= 0.0 .AND. arfs_wt_tb23v /= 0.0)
        arfs_tb_23v = arfs_tb_23v / arfs_wt_tb23v
    ELSEWHERE
        arfs_tb_23v = -9999.0
    END WHERE
    
    WHERE(arfs_tb_36h /= 0.0 .AND. arfs_wt_tb36h /= 0.0)
        arfs_tb_36h = arfs_tb_36h / arfs_wt_tb36h
    ELSEWHERE
        arfs_tb_36h = -9999.0
    END WHERE
    
    WHERE(arfs_tb_36v /= 0.0 .AND. arfs_wt_tb36v /= 0.0)
        arfs_tb_36v = arfs_tb_36v / arfs_wt_tb36v
    ELSEWHERE
        arfs_tb_36v = -9999.0
    END WHERE
    
    WHERE(arfs_tb_89h /= 0.0 .AND. arfs_wt_tb89h /= 0.0)
        arfs_tb_89h = arfs_tb_89h / arfs_wt_tb89h
    ELSEWHERE
        arfs_tb_89h = -9999.0
    END WHERE
    
    WHERE(arfs_tb_89v /= 0.0 .AND. arfs_wt_tb89v /= 0.0)
        arfs_tb_89v = arfs_tb_89v / arfs_wt_tb89v
    ELSEWHERE
        arfs_tb_89v = -9999.0
    END WHERE
    
    WHERE(arfs_land_water_frac /= 0.0 .AND. arfs_wt_land_water_frac /= 0.0)
        arfs_land_water_frac = arfs_land_water_frac / arfs_wt_land_water_frac
    ELSEWHERE
        arfs_land_water_frac = -9999.0
    END WHERE
    
    ! Finalize quality flags using majority vote
    write(LDT_logunit,*)'[INFO] Computing majority vote for quality flags...'
    
    do i = 1, 2560
        do j = 1, 1920
            if (arfs_quality_flag(i,j) == 0) then
                if (total_count(i,j) > 0) then
                    ! Ocean flag (bit 0)
                    if (ocean_count(i,j) > total_count(i,j)/2) then
                        arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 1)
                    endif
                    
                    ! Precipitation flag (bit 1)
                    if (precip_count(i,j) > total_count(i,j)/2) then
                        arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 2)
                    endif
                    
                    ! Snow flag (bit 2)
                    if (snow_count(i,j) > total_count(i,j)/2) then
                        arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 4)
                    endif
                endif
            endif
        end do
    end do
    
    ! Report statistics
    write(LDT_logunit,*)'[INFO] Resampling statistics:'
    write(LDT_logunit,*)'[INFO]   Total excluded snow footprints: ', SUM(excluded_snow_count)
    write(LDT_logunit,*)'[INFO]   Total excluded precip footprints: ', SUM(excluded_precip_count)
    write(LDT_logunit,*)'[INFO]   Total V-pol samples: ', SUM(arfs_sample_v)
    write(LDT_logunit,*)'[INFO]   Total H-pol samples: ', SUM(arfs_sample_h)
    
    ! Cleanup
    deallocate(zerodistflag)
    deallocate(arfs_wt_tim)
    deallocate(arfs_wt_tb10h, arfs_wt_tb10v)
    deallocate(arfs_wt_tb18h, arfs_wt_tb18v)
    deallocate(arfs_wt_tb23h, arfs_wt_tb23v)
    deallocate(arfs_wt_tb36h, arfs_wt_tb36v)
    deallocate(arfs_wt_tb89h, arfs_wt_tb89v)
    deallocate(arfs_wt_land_water_frac)
    deallocate(snow_count, precip_count)
    deallocate(ocean_count, total_count)
    deallocate(excluded_snow_count, excluded_precip_count)
    
    write(LDT_logunit,*)'[INFO] Resampling complete with snow/precip filtering'
    
    END SUBROUTINE WSF2ARFS_INVDIS

END MODULE invdist_wsf2arfs