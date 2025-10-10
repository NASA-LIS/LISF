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
!              EXACT match to AMSR_OPL processing
!
!-------------------------------------------------------------------------

MODULE invdist_wsf2arfs
    USE LDT_logMod, only: LDT_logunit
    IMPLICIT NONE

CONTAINS

    SUBROUTINE WSF2ARFS_INVDIS_EXACT(tim, tb_10h, tb_10v, tb_18h, tb_18v, &
        tb_23h, tb_23v, tb_36h, tb_36v, tb_89h, tb_89v, land_water_frac, &
        snow, precip, quality_flag, &
        lat_in, lon_in, nscans_in, nfovs_in, &
        ref_lat, ref_lon, &
        arfs_time, arfs_land_water_frac, &
        arfs_tb_10h, arfs_tb_10v, arfs_tb_18h, arfs_tb_18v, &
        arfs_tb_23h, arfs_tb_23v, arfs_tb_36h, arfs_tb_36v, &
        arfs_tb_89h, arfs_tb_89v, arfs_quality_flag, &
        arfs_sample_v, arfs_sample_h)
    
    ! Arguments - EXACT match to AMSR_OPL signature
    integer, intent(in) :: nscans_in, nfovs_in
    real*8, intent(in) :: tim(:)
    real*4, intent(in) :: tb_10h(nscans_in, nfovs_in), tb_10v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_18h(nscans_in, nfovs_in), tb_18v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_23h(nscans_in, nfovs_in), tb_23v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_36h(nscans_in, nfovs_in), tb_36v(nscans_in, nfovs_in)
    real*4, intent(in) :: tb_89h(nscans_in, nfovs_in), tb_89v(nscans_in, nfovs_in)
    real*4, intent(in) :: land_water_frac(nscans_in, nfovs_in)
    integer*4, intent(in) :: snow(nscans_in, nfovs_in), precip(nscans_in, nfovs_in)
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
    real*4 :: arfs_wt_tim(2560,1920)
    real*4 :: arfs_wt_tb10h(2560,1920), arfs_wt_tb10v(2560,1920)
    real*4 :: arfs_wt_tb18h(2560,1920), arfs_wt_tb18v(2560,1920)
    real*4 :: arfs_wt_tb23h(2560,1920), arfs_wt_tb23v(2560,1920)
    real*4 :: arfs_wt_tb36h(2560,1920), arfs_wt_tb36v(2560,1920)
    real*4 :: arfs_wt_tb89h(2560,1920), arfs_wt_tb89v(2560,1920)
    real*4 :: arfs_wt_land_water_frac(2560,1920)
    
    ! Quality flag tracking - EXACT match to AMSR_OPL
    integer :: snow_count(2560,1920), precip_count(2560,1920)
    integer :: ocean_count(2560,1920), total_count(2560,1920)
    integer :: excluded_snow_count(2560,1920), excluded_precip_count(2560,1920)
    
    write(LDT_logunit,*)'[INFO] Starting WSF inverse distance resampling'
    write(LDT_logunit,*)'[INFO] Using EXACT AMSR_OPL processing'
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
    
    ! Loop through all input footprints
    do jj = 1, nfovs_in
        if (mod(jj, 50) == 0) then
            write(LDT_logunit,*)'[INFO] Processing FOV ', jj, ' of ', nfovs_in
        end if
        
        do ii = 1, nscans_in
            ! Skip invalid coordinates
            if (lat_in(ii,jj) < -90.0 .or. lat_in(ii,jj) > 90.0 .or. &
                lon_in(ii,jj) < 0.0 .or. lon_in(ii,jj) > 360.0) then
                cycle
            endif
            
            ! Check if footprint has snow or precip
            has_snow = (IBITS(quality_flag(ii,jj), 2, 1) == 1)
            has_precip = (IBITS(quality_flag(ii,jj), 1, 1) == 1)
            
            ! Find nearest ARFS grid point
            c = MINLOC(ABS(lat_in(ii,jj)-ref_lat(:)), 1)
            r = MINLOC(ABS(lon_in(ii,jj)-ref_lon(:)), 1)
            
            ! Ensure valid indices
            if (r < 1 .or. r > size(ref_lon) .or. &
                c < 1 .or. c > size(ref_lat)) then
                cycle
            endif
            
            ! Define search window
            rmin = max(1, r-5)
            rmax = min(size(ref_lon), r+5)
            cmin = max(1, c-5)
            cmax = min(size(ref_lat), c+5)
            
            ! Bounds check
            if (rmin < 1 .or. rmax > size(ref_lon) .or. &
                cmin < 1 .or. cmax > size(ref_lat)) then
                cycle
            endif
            
            ! Loop through search window
            do cc = cmin, cmax
                lat2 = ref_lat(cc) * d2r
                
                do rr = rmin, rmax
                    lon2 = ref_lon(rr) * d2r
                    lat1 = lat_in(ii,jj) * d2r
                    lon1 = lon_in(ii,jj) * d2r
                    
                    ! Calculate great circle distance
                    if (lat1 == lat2 .and. lon1 == lon2) then
                        gcdist = 0.0
                    else
                        gcdist = RE_KM * DACOS(DSIN(lat1) * DSIN(lat2) + &
                            DCOS(lat1) * DCOS(lat2) * DCOS(lon1-lon2))
                    endif
                    
                    ! Process if within search radius
                    if (gcdist < search_radius) then
                        ! ALWAYS update counts for quality flag tracking (ALL footprints)
                        total_count(rr,cc) = total_count(rr,cc) + 1
                        if (IBITS(quality_flag(ii,jj), 0, 1) == 1) &
                            ocean_count(rr,cc) = ocean_count(rr,cc) + 1
                        if (IBITS(quality_flag(ii,jj), 1, 1) == 1) &
                            precip_count(rr,cc) = precip_count(rr,cc) + 1
                        if (IBITS(quality_flag(ii,jj), 2, 1) == 1) &
                            snow_count(rr,cc) = snow_count(rr,cc) + 1
                        
                        ! Check if this footprint should be excluded from resampling
                        if (has_snow .or. has_precip) then
                            ! Track excluded footprints
                            if (has_snow) excluded_snow_count(rr,cc) = excluded_snow_count(rr,cc) + 1
                            if (has_precip) excluded_precip_count(rr,cc) = excluded_precip_count(rr,cc) + 1
                        else
                            ! No snow/precip - proceed with actual resampling
                            if (gcdist < 0.0001D0) then
                                ! Exact match
                                zerodistflag(rr,cc) = 1
                                arfs_quality_flag(rr,cc) = quality_flag(ii,jj)
                                
                                ! Time
                                if (ii <= size(tim)) then
                                    if (tim(ii) >= 0) then
                                        arfs_time(rr,cc) = tim(ii)
                                        arfs_wt_tim(rr,cc) = 1.0
                                    endif
                                endif
                                
                                ! TB channels
                                if (tb_10h(ii,jj) > 0 .and. tb_10h(ii,jj) < 400) then
                                    arfs_tb_10h(rr,cc) = tb_10h(ii,jj)
                                    arfs_wt_tb10h(rr,cc) = 1.0
                                endif
                                if (tb_10v(ii,jj) > 0 .and. tb_10v(ii,jj) < 400) then
                                    arfs_tb_10v(rr,cc) = tb_10v(ii,jj)
                                    arfs_wt_tb10v(rr,cc) = 1.0
                                endif
                                if (tb_18h(ii,jj) > 0 .and. tb_18h(ii,jj) < 400) then
                                    arfs_tb_18h(rr,cc) = tb_18h(ii,jj)
                                    arfs_wt_tb18h(rr,cc) = 1.0
                                endif
                                if (tb_18v(ii,jj) > 0 .and. tb_18v(ii,jj) < 400) then
                                    arfs_tb_18v(rr,cc) = tb_18v(ii,jj)
                                    arfs_wt_tb18v(rr,cc) = 1.0
                                endif
                                if (tb_23h(ii,jj) > 0 .and. tb_23h(ii,jj) < 400) then
                                    arfs_tb_23h(rr,cc) = tb_23h(ii,jj)
                                    arfs_wt_tb23h(rr,cc) = 1.0
                                endif
                                if (tb_23v(ii,jj) > 0 .and. tb_23v(ii,jj) < 400) then
                                    arfs_tb_23v(rr,cc) = tb_23v(ii,jj)
                                    arfs_wt_tb23v(rr,cc) = 1.0
                                endif
                                if (tb_36h(ii,jj) > 0 .and. tb_36h(ii,jj) < 400) then
                                    arfs_tb_36h(rr,cc) = tb_36h(ii,jj)
                                    arfs_wt_tb36h(rr,cc) = 1.0
                                endif
                                if (tb_36v(ii,jj) > 0 .and. tb_36v(ii,jj) < 400) then
                                    arfs_tb_36v(rr,cc) = tb_36v(ii,jj)
                                    arfs_wt_tb36v(rr,cc) = 1.0
                                endif
                                if (tb_89h(ii,jj) > 0 .and. tb_89h(ii,jj) < 400) then
                                    arfs_tb_89h(rr,cc) = tb_89h(ii,jj)
                                    arfs_wt_tb89h(rr,cc) = 1.0
                                endif
                                if (tb_89v(ii,jj) > 0 .and. tb_89v(ii,jj) < 400) then
                                    arfs_tb_89v(rr,cc) = tb_89v(ii,jj)
                                    arfs_wt_tb89v(rr,cc) = 1.0
                                endif
                                
                                ! Land water fraction
                                if (land_water_frac(ii,jj) >= 0) then
                                    arfs_land_water_frac(rr,cc) = land_water_frac(ii,jj)
                                    arfs_wt_land_water_frac(rr,cc) = 1.0
                                endif
                                
                            else
                                ! Weighted by inverse distance
                                if (zerodistflag(rr,cc) /= 1) then
                                    ! Time
                                    if (ii <= size(tim)) then
                                        if (tim(ii) >= 0) then
                                            arfs_time(rr,cc) = arfs_time(rr,cc) + tim(ii) / gcdist
                                            arfs_wt_tim(rr,cc) = arfs_wt_tim(rr,cc) + 1.0 / gcdist
                                        endif
                                    endif
                                    
                                    ! TB channels
                                    if (tb_10h(ii,jj) > 0 .and. tb_10h(ii,jj) < 400) then
                                        arfs_tb_10h(rr,cc) = arfs_tb_10h(rr,cc) + tb_10h(ii,jj) / gcdist
                                        arfs_wt_tb10h(rr,cc) = arfs_wt_tb10h(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_10v(ii,jj) > 0 .and. tb_10v(ii,jj) < 400) then
                                        arfs_tb_10v(rr,cc) = arfs_tb_10v(rr,cc) + tb_10v(ii,jj) / gcdist
                                        arfs_wt_tb10v(rr,cc) = arfs_wt_tb10v(rr,cc) + 1.0 / gcdist
                                    endif
                                    ! ... (similar for all other channels)
                                    if (tb_18h(ii,jj) > 0 .and. tb_18h(ii,jj) < 400) then
                                        arfs_tb_18h(rr,cc) = arfs_tb_18h(rr,cc) + tb_18h(ii,jj) / gcdist
                                        arfs_wt_tb18h(rr,cc) = arfs_wt_tb18h(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_18v(ii,jj) > 0 .and. tb_18v(ii,jj) < 400) then
                                        arfs_tb_18v(rr,cc) = arfs_tb_18v(rr,cc) + tb_18v(ii,jj) / gcdist
                                        arfs_wt_tb18v(rr,cc) = arfs_wt_tb18v(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_23h(ii,jj) > 0 .and. tb_23h(ii,jj) < 400) then
                                        arfs_tb_23h(rr,cc) = arfs_tb_23h(rr,cc) + tb_23h(ii,jj) / gcdist
                                        arfs_wt_tb23h(rr,cc) = arfs_wt_tb23h(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_23v(ii,jj) > 0 .and. tb_23v(ii,jj) < 400) then
                                        arfs_tb_23v(rr,cc) = arfs_tb_23v(rr,cc) + tb_23v(ii,jj) / gcdist
                                        arfs_wt_tb23v(rr,cc) = arfs_wt_tb23v(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_36h(ii,jj) > 0 .and. tb_36h(ii,jj) < 400) then
                                        arfs_tb_36h(rr,cc) = arfs_tb_36h(rr,cc) + tb_36h(ii,jj) / gcdist
                                        arfs_wt_tb36h(rr,cc) = arfs_wt_tb36h(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_36v(ii,jj) > 0 .and. tb_36v(ii,jj) < 400) then
                                        arfs_tb_36v(rr,cc) = arfs_tb_36v(rr,cc) + tb_36v(ii,jj) / gcdist
                                        arfs_wt_tb36v(rr,cc) = arfs_wt_tb36v(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_89h(ii,jj) > 0 .and. tb_89h(ii,jj) < 400) then
                                        arfs_tb_89h(rr,cc) = arfs_tb_89h(rr,cc) + tb_89h(ii,jj) / gcdist
                                        arfs_wt_tb89h(rr,cc) = arfs_wt_tb89h(rr,cc) + 1.0 / gcdist
                                    endif
                                    if (tb_89v(ii,jj) > 0 .and. tb_89v(ii,jj) < 400) then
                                        arfs_tb_89v(rr,cc) = arfs_tb_89v(rr,cc) + tb_89v(ii,jj) / gcdist
                                        arfs_wt_tb89v(rr,cc) = arfs_wt_tb89v(rr,cc) + 1.0 / gcdist
                                    endif
                                    
                                    ! Land water fraction
                                    if (land_water_frac(ii,jj) >= 0) then
                                        arfs_land_water_frac(rr,cc) = arfs_land_water_frac(rr,cc) + &
                                            land_water_frac(ii,jj) / gcdist
                                        arfs_wt_land_water_frac(rr,cc) = arfs_wt_land_water_frac(rr,cc) + &
                                            1.0 / gcdist
                                    endif
                                endif ! zerodistflag check
                            endif ! gcdist < 0.0001
                        endif ! has_snow .or. has_precip
                    endif ! gcdist < search_radius
                end do !rr
            end do !cc
        end do !ii
    end do !jj
    
    ! Apply weighting and set fill values
    write(LDT_logunit,*)'[INFO] Applying weights and computing averages...'
    
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
    
    ! ... (similar for all channels)
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
                    ! Ocean flag
                    if (ocean_count(i,j) > total_count(i,j)/2) then
                        arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 1)
                    endif
                    
                    ! Precipitation flag
                    if (precip_count(i,j) > total_count(i,j)/2) then
                        arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 2)
                    endif
                    
                    ! Snow flag
                    if (snow_count(i,j) > total_count(i,j)/2) then
                        arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 4)
                    endif
                endif
            endif
        end do
    end do
    
    ! Cleanup
    deallocate(zerodistflag)
    
    write(LDT_logunit,*)'[INFO] Resampling complete'
    
    END SUBROUTINE WSF2ARFS_INVDIS_EXACT

END MODULE invdist_wsf2arfs