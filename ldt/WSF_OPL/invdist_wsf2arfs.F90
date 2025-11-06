!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: invdist_wsf2arfs
!
! DESCRIPTION: Inverse distance resampling with proper per-channel quality filtering
!
!-------------------------------------------------------------------------

MODULE invdist_wsf2arfs
    USE LDT_logMod, only: LDT_logunit
    IMPLICIT NONE

CONTAINS

    SUBROUTINE WSF2ARFS_INVDIS(tim, tb_10h, tb_10v, tb_18h, tb_18v, &
        tb_23h, tb_23v, tb_36h, tb_36v, tb_89h, tb_89v, land_water_frac, &
        snow_in, precip_in, quality_flag, &
        lat_in, lon_in, nscans_in, nfovs_in, &
        ref_lat, ref_lon, &
        arfs_time, arfs_land_water_frac, &
        arfs_tb_10h, arfs_tb_10v, arfs_tb_18h, arfs_tb_18v, &
        arfs_tb_23h, arfs_tb_23v, arfs_tb_36h, arfs_tb_36v, &
        arfs_tb_89h, arfs_tb_89v, arfs_quality_flag, &
        arfs_sample_v, arfs_sample_h)
    
    ! Arguments
    integer, intent(in) :: nscans_in, nfovs_in
    real*8, intent(in) :: tim(nscans_in)
    real*4, intent(in) :: tb_10h(nscans_in, nfovs_in), tb_10v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_18h(nscans_in, nfovs_in), tb_18v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_23h(nscans_in, nfovs_in), tb_23v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_36h(nscans_in, nfovs_in), tb_36v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_89h(nscans_in, nfovs_in), tb_89v(nscans_in, nfovs_in)
    real*4, intent(in) :: land_water_frac(nscans_in, nfovs_in)
    integer*4, intent(in) :: snow_in(nscans_in, nfovs_in), precip_in(nscans_in, nfovs_in)
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
    
    ! Local variables
    integer :: ii, jj, r, c, rr, cc, rmin, rmax, cmin, cmax
    real*8, parameter :: RE_KM = 6371.228
    real*8, parameter :: search_radius = 20.0
    real*8, parameter :: PI = 3.141592653589793238
    real*8, parameter :: d2r = PI/180.0
    real*8 :: gcdist, lat1, lon1, lat2, lon2, weight
    logical :: has_snow, has_precip, has_water
    logical :: band1_bad, band2_bad, band3_bad, band4_bad, band5_bad
    
    ! Weight arrays for accumulation
    real*4, allocatable :: arfs_wt_tb10h(:,:), arfs_wt_tb10v(:,:)
    real*4, allocatable :: arfs_wt_tb18h(:,:), arfs_wt_tb18v(:,:)
    real*4, allocatable :: arfs_wt_tb23h(:,:), arfs_wt_tb23v(:,:)
    real*4, allocatable :: arfs_wt_tb36h(:,:), arfs_wt_tb36v(:,:)
    real*4, allocatable :: arfs_wt_tb89h(:,:), arfs_wt_tb89v(:,:)
    real*4, allocatable :: arfs_wt_land_water_frac(:,:)
    
    ! Quality flag counters
    integer, allocatable :: ocean_count(:,:), precip_count(:,:), snow_count(:,:)
    integer, allocatable :: total_count(:,:)
    integer, allocatable :: band1_bad_count(:,:), band2_bad_count(:,:)
    integer, allocatable :: band3_bad_count(:,:), band4_bad_count(:,:)
    integer, allocatable :: band5_bad_count(:,:)
    
    ! Statistics counters
    integer :: excluded_snow, excluded_precip, excluded_water
    integer :: excluded_band1, excluded_band2, excluded_band3
    integer :: excluded_band4, excluded_band5
    
    write(LDT_logunit,*)'[INFO] Starting WSF inverse distance resampling'
    write(LDT_logunit,*)'[INFO] WITH per-channel quality filtering'
    
    ! Allocate weight and counter arrays
    allocate(arfs_wt_tb10h(2560,1920), arfs_wt_tb10v(2560,1920))
    allocate(arfs_wt_tb18h(2560,1920), arfs_wt_tb18v(2560,1920))
    allocate(arfs_wt_tb23h(2560,1920), arfs_wt_tb23v(2560,1920))
    allocate(arfs_wt_tb36h(2560,1920), arfs_wt_tb36v(2560,1920))
    allocate(arfs_wt_tb89h(2560,1920), arfs_wt_tb89v(2560,1920))
    allocate(arfs_wt_land_water_frac(2560,1920))
    
    allocate(ocean_count(2560,1920), precip_count(2560,1920), snow_count(2560,1920))
    allocate(total_count(2560,1920))
    allocate(band1_bad_count(2560,1920), band2_bad_count(2560,1920))
    allocate(band3_bad_count(2560,1920), band4_bad_count(2560,1920))
    allocate(band5_bad_count(2560,1920))
    
    ! Initialize arrays
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
    arfs_quality_flag = 0
    arfs_sample_v = 0
    arfs_sample_h = 0
    
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
    
    ocean_count = 0
    precip_count = 0
    snow_count = 0
    total_count = 0
    band1_bad_count = 0
    band2_bad_count = 0
    band3_bad_count = 0
    band4_bad_count = 0
    band5_bad_count = 0
    
    excluded_snow = 0
    excluded_precip = 0
    excluded_water = 0
    excluded_band1 = 0
    excluded_band2 = 0
    excluded_band3 = 0
    excluded_band4 = 0
    excluded_band5 = 0
    
    ! =====================================================================
    ! MAIN RESAMPLING LOOP WITH PER-CHANNEL QUALITY FILTERING
    ! =====================================================================
    
    do ii = 1, nfovs_in
        do jj = 1, nscans_in
            ! Skip invalid coordinates
            if (lat_in(jj,ii) < -90.0 .or. lat_in(jj,ii) > 90.0 .or. &
                lon_in(jj,ii) < -180.0 .or. lon_in(jj,ii) > 360.0) cycle
            
            ! Extract quality conditions for this footprint
            has_water = (IBITS(quality_flag(jj,ii), 0, 1) == 1)
            has_precip = (IBITS(quality_flag(jj,ii), 1, 1) == 1)
            has_snow = (IBITS(quality_flag(jj,ii), 2, 1) == 1)
            
            ! Extract band-specific quality (1 = bad, 0 = good)
            band1_bad = (IBITS(quality_flag(jj,ii), 3, 1) == 1)
            band2_bad = (IBITS(quality_flag(jj,ii), 4, 1) == 1)
            band3_bad = (IBITS(quality_flag(jj,ii), 5, 1) == 1)
            band4_bad = (IBITS(quality_flag(jj,ii), 6, 1) == 1)
            band5_bad = (IBITS(quality_flag(jj,ii), 7, 1) == 1)
            
            ! Find grid cell bounds
            lat1 = lat_in(jj,ii)
            lon1 = lon_in(jj,ii)

            ! Convert longitude from 0-360 to -180 to +180 convention
            if (lon1 > 180.0) then
                lon1 = lon1 - 360.0
            endif
            
            c = MINLOC(ABS(lat1 - ref_lat(:)), 1)
            r = MINLOC(ABS(lon1 - ref_lon(:)), 1)
            
            rmin = max(1, r - 5)
            rmax = min(size(ref_lon), r + 5)
            cmin = max(1, c - 5)
            cmax = min(size(ref_lat), c + 5)
            
            ! Process each grid cell in search radius
            do rr = rmin, rmax
                do cc = cmin, cmax
                    lat2 = ref_lat(cc)
                    lon2 = ref_lon(rr)
                    
                    ! Calculate great circle distance
                    gcdist = 2.0 * RE_KM * ASIN(SQRT((SIN(d2r*(lat1-lat2)/2.0))**2 + &
                             COS(d2r*lat1) * COS(d2r*lat2) * (SIN(d2r*(lon1-lon2)/2.0))**2))
                    
                    if (gcdist < search_radius) then
                        ! Update quality flag counters
                        total_count(rr,cc) = total_count(rr,cc) + 1
                        
                        if (has_water) ocean_count(rr,cc) = ocean_count(rr,cc) + 1
                        if (has_precip) precip_count(rr,cc) = precip_count(rr,cc) + 1
                        if (has_snow) snow_count(rr,cc) = snow_count(rr,cc) + 1
                        if (band1_bad .and. .not. has_water) band1_bad_count(rr,cc) = band1_bad_count(rr,cc) + 1
                        if (band2_bad .and. .not. has_water) band2_bad_count(rr,cc) = band2_bad_count(rr,cc) + 1
                        if (band3_bad .and. .not. has_water) band3_bad_count(rr,cc) = band3_bad_count(rr,cc) + 1
                        if (band4_bad .and. .not. has_water) band4_bad_count(rr,cc) = band4_bad_count(rr,cc) + 1
                        if (band5_bad .and. .not. has_water) band5_bad_count(rr,cc) = band5_bad_count(rr,cc) + 1
                        
                        ! Calculate weight (inverse distance or 1.0 for exact match)
                        if (gcdist < 0.0001D0) then
                            weight = 1.0D0
                        else
                            weight = 1.0D0 / gcdist
                        endif
                        
                        ! ============================================================
                        ! RESAMPLE EACH CHANNEL BASED ON ITS QUALITY
                        ! ============================================================
                        
                        ! 10 GHz channels (Band 1)
                        if (.not. band1_bad .and. .not. has_snow .and. .not. has_precip) then
                            if (ABS(tb_10h(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_10h(jj,ii) > 0.0) then
                                arfs_tb_10h(rr,cc) = arfs_tb_10h(rr,cc) + tb_10h(jj,ii) * weight
                                arfs_wt_tb10h(rr,cc) = arfs_wt_tb10h(rr,cc) + weight
                                arfs_sample_h(rr,cc) = arfs_sample_h(rr,cc) + 1
                            endif
                            if (ABS(tb_10v(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_10v(jj,ii) > 0.0) then
                                arfs_tb_10v(rr,cc) = arfs_tb_10v(rr,cc) + tb_10v(jj,ii) * weight
                                arfs_wt_tb10v(rr,cc) = arfs_wt_tb10v(rr,cc) + weight
                                arfs_sample_v(rr,cc) = arfs_sample_v(rr,cc) + 1
                            endif
                        else if (band1_bad) then
                            excluded_band1 = excluded_band1 + 1
                        endif
                        
                        ! 18 GHz channels (Band 2)
                        if (.not. band2_bad .and. .not. has_snow .and. .not. has_precip) then
                            if (ABS(tb_18h(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_18h(jj,ii) > 0.0) then
                                arfs_tb_18h(rr,cc) = arfs_tb_18h(rr,cc) + tb_18h(jj,ii) * weight
                                arfs_wt_tb18h(rr,cc) = arfs_wt_tb18h(rr,cc) + weight
                            endif
                            if (ABS(tb_18v(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_18v(jj,ii) > 0.0) then
                                arfs_tb_18v(rr,cc) = arfs_tb_18v(rr,cc) + tb_18v(jj,ii) * weight
                                arfs_wt_tb18v(rr,cc) = arfs_wt_tb18v(rr,cc) + weight
                            endif
                        else if (band2_bad) then
                            excluded_band2 = excluded_band2 + 1
                        endif
                        
                        ! 23 GHz channels (Band 3)
                        if (.not. band3_bad .and. .not. has_snow .and. .not. has_precip) then
                            if (ABS(tb_23h(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_23h(jj,ii) > 0.0) then
                                arfs_tb_23h(rr,cc) = arfs_tb_23h(rr,cc) + tb_23h(jj,ii) * weight
                                arfs_wt_tb23h(rr,cc) = arfs_wt_tb23h(rr,cc) + weight
                            endif
                            if (ABS(tb_23v(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_23v(jj,ii) > 0.0) then
                                arfs_tb_23v(rr,cc) = arfs_tb_23v(rr,cc) + tb_23v(jj,ii) * weight
                                arfs_wt_tb23v(rr,cc) = arfs_wt_tb23v(rr,cc) + weight
                            endif
                        else if (band3_bad) then
                            excluded_band3 = excluded_band3 + 1
                        endif
                        
                        ! 36 GHz channels (Band 4)
                        if (.not. band4_bad .and. .not. has_snow .and. .not. has_precip) then
                            if (ABS(tb_36h(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_36h(jj,ii) > 0.0) then
                                arfs_tb_36h(rr,cc) = arfs_tb_36h(rr,cc) + tb_36h(jj,ii) * weight
                                arfs_wt_tb36h(rr,cc) = arfs_wt_tb36h(rr,cc) + weight
                            endif
                            if (ABS(tb_36v(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_36v(jj,ii) > 0.0) then
                                arfs_tb_36v(rr,cc) = arfs_tb_36v(rr,cc) + tb_36v(jj,ii) * weight
                                arfs_wt_tb36v(rr,cc) = arfs_wt_tb36v(rr,cc) + weight
                            endif
                        else if (band4_bad) then
                            excluded_band4 = excluded_band4 + 1
                        endif
                        
                        ! 89 GHz channels (Band 5)
                        if (.not. band5_bad .and. .not. has_snow .and. .not. has_precip) then
                            if (ABS(tb_89h(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_89h(jj,ii) > 0.0) then
                                arfs_tb_89h(rr,cc) = arfs_tb_89h(rr,cc) + tb_89h(jj,ii) * weight
                                arfs_wt_tb89h(rr,cc) = arfs_wt_tb89h(rr,cc) + weight
                            endif
                            if (ABS(tb_89v(jj,ii) - (-9999.0)) > 1.0E-6 .and. tb_89v(jj,ii) > 0.0) then
                                arfs_tb_89v(rr,cc) = arfs_tb_89v(rr,cc) + tb_89v(jj,ii) * weight
                                arfs_wt_tb89v(rr,cc) = arfs_wt_tb89v(rr,cc) + weight
                            endif
                        else if (band5_bad) then
                            excluded_band5 = excluded_band5 + 1
                        endif
                        
                        ! Land fraction (always resample if valid)
                        if (ABS(land_water_frac(jj,ii) - (-9999.0)) > 1.0E-6 .and. &
                            land_water_frac(jj,ii) >= 0.0) then
                            arfs_land_water_frac(rr,cc) = arfs_land_water_frac(rr,cc) + &
                                                          land_water_frac(jj,ii) * weight
                            arfs_wt_land_water_frac(rr,cc) = arfs_wt_land_water_frac(rr,cc) + weight
                        endif
                        
                        ! Track exclusions
                        if (has_snow) excluded_snow = excluded_snow + 1
                        if (has_precip) excluded_precip = excluded_precip + 1
                        if (has_water) excluded_water = excluded_water + 1
                        
                    endif ! within search radius
                end do
            end do
        end do
    end do
    
    ! =====================================================================
    ! NORMALIZE BY WEIGHTS
    ! =====================================================================
    
    WHERE(arfs_wt_tb10h > 0.0)
        arfs_tb_10h = arfs_tb_10h / arfs_wt_tb10h
    ELSEWHERE
        arfs_tb_10h = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb10v > 0.0)
        arfs_tb_10v = arfs_tb_10v / arfs_wt_tb10v
    ELSEWHERE
        arfs_tb_10v = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb18h > 0.0)
        arfs_tb_18h = arfs_tb_18h / arfs_wt_tb18h
    ELSEWHERE
        arfs_tb_18h = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb18v > 0.0)
        arfs_tb_18v = arfs_tb_18v / arfs_wt_tb18v
    ELSEWHERE
        arfs_tb_18v = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb23h > 0.0)
        arfs_tb_23h = arfs_tb_23h / arfs_wt_tb23h
    ELSEWHERE
        arfs_tb_23h = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb23v > 0.0)
        arfs_tb_23v = arfs_tb_23v / arfs_wt_tb23v
    ELSEWHERE
        arfs_tb_23v = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb36h > 0.0)
        arfs_tb_36h = arfs_tb_36h / arfs_wt_tb36h
    ELSEWHERE
        arfs_tb_36h = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb36v > 0.0)
        arfs_tb_36v = arfs_tb_36v / arfs_wt_tb36v
    ELSEWHERE
        arfs_tb_36v = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb89h > 0.0)
        arfs_tb_89h = arfs_tb_89h / arfs_wt_tb89h
    ELSEWHERE
        arfs_tb_89h = -9999.0
    END WHERE
    
    WHERE(arfs_wt_tb89v > 0.0)
        arfs_tb_89v = arfs_tb_89v / arfs_wt_tb89v
    ELSEWHERE
        arfs_tb_89v = -9999.0
    END WHERE
    
    WHERE(arfs_wt_land_water_frac > 0.0)
        arfs_land_water_frac = arfs_land_water_frac / arfs_wt_land_water_frac
    ELSEWHERE
        arfs_land_water_frac = -9999.0
    END WHERE
    
    ! =====================================================================
    ! CREATE QUALITY FLAGS USING MAJORITY VOTE
    ! =====================================================================
    
    do rr = 1, 2560
        do cc = 1, 1920
            if (total_count(rr,cc) > 0) then
                arfs_quality_flag(rr,cc) = 0
                
                ! Bit 0: Ocean
                if (ocean_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 1)
                endif
                
                ! Bit 1: Precipitation
                if (precip_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 2)
                endif
                
                ! Bit 2: Snow
                if (snow_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 4)
                endif
                
                ! Bit 3: Band 1 (10 GHz) quality
                if (band1_bad_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 8)
                endif
                
                ! Bit 4: Band 2 (18 GHz) quality
                if (band2_bad_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 16)
                endif
                
                ! Bit 5: Band 3 (23 GHz) quality
                if (band3_bad_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 32)
                endif
                
                ! Bit 6: Band 4 (36 GHz) quality
                if (band4_bad_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 64)
                endif
                
                ! Bit 7: Band 5 (89 GHz) quality
                if (band5_bad_count(rr,cc) > total_count(rr,cc)/2) then
                    arfs_quality_flag(rr,cc) = IOR(arfs_quality_flag(rr,cc), 128)
                endif
            else
                arfs_quality_flag(rr,cc) = -1  ! No data
            endif
        end do
    end do
    
    ! Report statistics
    write(LDT_logunit,*)'[INFO] ========================================'
    write(LDT_logunit,*)'[INFO] Resampling statistics:'
    write(LDT_logunit,*)'[INFO]   Excluded snow footprints: ', excluded_snow
    write(LDT_logunit,*)'[INFO]   Excluded precip footprints: ', excluded_precip
    write(LDT_logunit,*)'[INFO]   Excluded water footprints: ', excluded_water
    write(LDT_logunit,*)'[INFO]   Excluded bad band 1 (10GHz): ', excluded_band1
    write(LDT_logunit,*)'[INFO]   Excluded bad band 2 (18GHz): ', excluded_band2
    write(LDT_logunit,*)'[INFO]   Excluded bad band 3 (23GHz): ', excluded_band3
    write(LDT_logunit,*)'[INFO]   Excluded bad band 4 (36GHz): ', excluded_band4
    write(LDT_logunit,*)'[INFO]   Excluded bad band 5 (89GHz): ', excluded_band5
    write(LDT_logunit,*)'[INFO]   Total V-pol samples: ', SUM(arfs_sample_v)
    write(LDT_logunit,*)'[INFO]   Total H-pol samples: ', SUM(arfs_sample_h)
    write(LDT_logunit,*)'[INFO] Quality flag bits set in output:'
    write(LDT_logunit,*)'[INFO]   Ocean (bit 0): ', count(IBITS(arfs_quality_flag, 0, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Precip (bit 1): ', count(IBITS(arfs_quality_flag, 1, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Snow (bit 2): ', count(IBITS(arfs_quality_flag, 2, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 10GHz (bit 3): ', count(IBITS(arfs_quality_flag, 3, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 18GHz (bit 4): ', count(IBITS(arfs_quality_flag, 4, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 23GHz (bit 5): ', count(IBITS(arfs_quality_flag, 5, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 36GHz (bit 6): ', count(IBITS(arfs_quality_flag, 6, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 89GHz (bit 7): ', count(IBITS(arfs_quality_flag, 7, 1) == 1)
    write(LDT_logunit,*)'[INFO] ========================================'
    
    ! Cleanup
    deallocate(arfs_wt_tb10h, arfs_wt_tb10v)
    deallocate(arfs_wt_tb18h, arfs_wt_tb18v)
    deallocate(arfs_wt_tb23h, arfs_wt_tb23v)
    deallocate(arfs_wt_tb36h, arfs_wt_tb36v)
    deallocate(arfs_wt_tb89h, arfs_wt_tb89v)
    deallocate(arfs_wt_land_water_frac)
    deallocate(ocean_count, precip_count, snow_count, total_count)
    deallocate(band1_bad_count, band2_bad_count, band3_bad_count)
    deallocate(band4_bad_count, band5_bad_count)
    
    END SUBROUTINE WSF2ARFS_INVDIS

END MODULE invdist_wsf2arfs
