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
! DESCRIPTION: Subroutine to resample WSF TB onto Air Force Grid
!
!-------------------------------------------------------------------------

MODULE invdist_wsf2arfs
    USE LDT_logMod, only: LDT_logunit
    IMPLICIT NONE

CONTAINS

    SUBROUTINE WSF2ARFS_INVDIS(tb_in, lat_in, lon_in, &
        quality_flag_in, land_frac_in, &
        nscans_in, nfovs_in, nchans_in, &
        ref_lat, ref_lon, &
        tb_out, quality_flag_out, land_frac_out, &
        sample_count)
    
    ! Arguments
    integer, intent(in) :: nscans_in, nfovs_in, nchans_in
    real*4, intent(in) :: tb_in(nchans_in, nscans_in, nfovs_in)
    real*4, intent(in) :: lat_in(nscans_in, nfovs_in)
    real*4, intent(in) :: lon_in(nscans_in, nfovs_in)
    integer*1, intent(in) :: quality_flag_in(:,:,:)  ! (nbands, nscans, nfovs)
    real*4, intent(in) :: land_frac_in(nscans_in, nfovs_in)
    real*8, intent(in) :: ref_lat(:)
    real*8, intent(in) :: ref_lon(:)
    real*4, intent(out) :: tb_out(nchans_in, size(ref_lon), size(ref_lat))
    integer*1, intent(out) :: quality_flag_out(size(ref_lon), size(ref_lat))
    real*4, intent(out) :: land_frac_out(size(ref_lon), size(ref_lat))
    integer*4, intent(out) :: sample_count(size(ref_lon), size(ref_lat))
    
    ! Local variables
    integer :: ii, jj, ichan, r, c, rr, cc, rmin, rmax, cmin, cmax
    integer :: zerodistflag(size(ref_lon), size(ref_lat))
    real*8, parameter :: RE_KM = 6371.228
    real*8, parameter :: SEARCH_RADIUS = 20.0
    real*8, parameter :: PI = 3.141592653589793238
    real*8, parameter :: D2R = PI/180.0
    real*8 :: gcdist, lat1, lon1, lat2, lon2
    real*4 :: wt_tb(nchans_in, size(ref_lon), size(ref_lat))
    real*4 :: wt_land(size(ref_lon), size(ref_lat))
    integer :: quality_sum(size(ref_lon), size(ref_lat))
    integer :: quality_count(size(ref_lon), size(ref_lat))
    logical :: is_valid
    integer :: band_idx
    
    write(LDT_logunit,*)'[INFO] Starting inverse distance resampling'
    write(LDT_logunit,*)'[INFO] Input dimensions: ', nscans_in, 'x', nfovs_in
    write(LDT_logunit,*)'[INFO] Output dimensions: ', size(ref_lon), 'x', size(ref_lat)
    write(LDT_logunit,*)'[INFO] Number of channels: ', nchans_in
    
    ! Initialize output arrays
    zerodistflag = 0
    tb_out = 0.0
    wt_tb = 0.0
    land_frac_out = 0.0
    wt_land = 0.0
    sample_count = 0
    quality_flag_out = 0
    quality_sum = 0
    quality_count = 0
    
    ! Loop through all input points
    do jj = 1, nfovs_in
        if (mod(jj, 50) == 0) then
            write(LDT_logunit,*)'[INFO] Processing FOV ', jj, ' of ', nfovs_in
        end if
        
        do ii = 1, nscans_in
            ! Skip invalid coordinates
            if (lat_in(ii,jj) < -90.0 .or. lat_in(ii,jj) > 90.0 .or. &
                lon_in(ii,jj) < 0.0 .or. lon_in(ii,jj) > 360.0) then
                cycle
            end if
            
            ! Check data quality - use first band as representative
            is_valid = .true.
            if (size(quality_flag_in, 1) > 0) then
                band_idx = 1  ! Use first band for quality check
                ! Bit 0 = IsNotValid
                if (BTEST(quality_flag_in(band_idx,ii,jj), 0)) then
                    is_valid = .false.
                end if
            end if
            
            ! Skip if invalid
            if (.not. is_valid) cycle
            
            ! Find nearest ARFS grid point
            c = MINLOC(ABS(lat_in(ii,jj)-ref_lat(:)), 1)
            r = MINLOC(ABS(lon_in(ii,jj)-ref_lon(:)), 1)
            
            ! Ensure valid indices
            if (r < 1 .or. r > size(ref_lon) .or. &
                c < 1 .or. c > size(ref_lat)) then
                cycle
            end if
            
            ! Define search window
            rmin = max(1, r-5)
            rmax = min(size(ref_lon), r+5)
            cmin = max(1, c-5)
            cmax = min(size(ref_lat), c+5)
            
            ! Loop through search window
            do cc = cmin, cmax
                lat2 = ref_lat(cc) * D2R
                
                do rr = rmin, rmax
                    lon2 = ref_lon(rr) * D2R
                    lat1 = lat_in(ii,jj) * D2R
                    
                    ! Convert longitude to -180 to 180 range for distance calculation
                    lon1 = lon_in(ii,jj)
                    if (lon1 > 180.0) lon1 = lon1 - 360.0
                    lon1 = lon1 * D2R
                    
                    ! Calculate great circle distance
                    if (lat1 == lat2 .and. lon1 == lon2) then
                        gcdist = 0.0
                    else
                        gcdist = RE_KM * DACOS(DSIN(lat1) * DSIN(lat2) + &
                            DCOS(lat1) * DCOS(lat2) * DCOS(lon1-lon2))
                    end if
                    
                    ! Process if within search radius
                    if (gcdist < SEARCH_RADIUS) then
                        sample_count(rr,cc) = sample_count(rr,cc) + 1
                        
                        if (gcdist < 0.0001) then
                            ! Exact match
                            zerodistflag(rr,cc) = 1
                            
                            do ichan = 1, nchans_in
                                if (tb_in(ichan,ii,jj) > 0.0 .and. &
                                    tb_in(ichan,ii,jj) < 400.0) then
                                    tb_out(ichan,rr,cc) = tb_in(ichan,ii,jj)
                                    wt_tb(ichan,rr,cc) = 1.0
                                end if
                            end do
                            
                            if (land_frac_in(ii,jj) >= 0.0) then
                                land_frac_out(rr,cc) = land_frac_in(ii,jj)
                                wt_land(rr,cc) = 1.0
                            end if
                            
                            ! Copy quality flag for exact match
                            if (size(quality_flag_in, 1) > 0) then
                                quality_flag_out(rr,cc) = quality_flag_in(1,ii,jj)
                            end if
                            
                        else
                            ! Weighted by inverse distance
                            if (zerodistflag(rr,cc) /= 1) then
                                do ichan = 1, nchans_in
                                    if (tb_in(ichan,ii,jj) > 0.0 .and. &
                                        tb_in(ichan,ii,jj) < 400.0) then
                                        tb_out(ichan,rr,cc) = tb_out(ichan,rr,cc) + &
                                            tb_in(ichan,ii,jj) / gcdist
                                        wt_tb(ichan,rr,cc) = wt_tb(ichan,rr,cc) + &
                                            1.0 / gcdist
                                    end if
                                end do
                                
                                if (land_frac_in(ii,jj) >= 0.0) then
                                    land_frac_out(rr,cc) = land_frac_out(rr,cc) + &
                                        land_frac_in(ii,jj) / gcdist
                                    wt_land(rr,cc) = wt_land(rr,cc) + 1.0 / gcdist
                                end if
                                
                                ! Accumulate quality flags for majority vote
                                if (size(quality_flag_in, 1) > 0) then
                                    quality_sum(rr,cc) = quality_sum(rr,cc) + &
                                        quality_flag_in(1,ii,jj)
                                    quality_count(rr,cc) = quality_count(rr,cc) + 1
                                end if
                            end if
                        end if
                    end if
                end do  ! rr
            end do  ! cc
        end do  ! ii
    end do  ! jj
    
    write(LDT_logunit,*)'[INFO] Applying weights to resampled data'
    
    ! Apply weights
    do cc = 1, size(ref_lat)
        do rr = 1, size(ref_lon)
            if (zerodistflag(rr,cc) == 0) then
                ! Not an exact match - apply weights
                do ichan = 1, nchans_in
                    if (wt_tb(ichan,rr,cc) > 0.0) then
                        tb_out(ichan,rr,cc) = tb_out(ichan,rr,cc) / wt_tb(ichan,rr,cc)
                    end if
                end do
                
                if (wt_land(rr,cc) > 0.0) then
                    land_frac_out(rr,cc) = land_frac_out(rr,cc) / wt_land(rr,cc)
                end if
                
                ! Average quality flag
                if (quality_count(rr,cc) > 0) then
                    quality_flag_out(rr,cc) = int(quality_sum(rr,cc) / &
                        quality_count(rr,cc), kind=1)
                end if
            end if
        end do
    end do
    
    ! Count valid output points
    write(LDT_logunit,*)'[INFO] Resampling complete'
    write(LDT_logunit,*)'[INFO] Valid output points: ', &
        count(sample_count > 0)
    write(LDT_logunit,*)'[INFO] Max samples per point: ', &
        maxval(sample_count)
    
    END SUBROUTINE WSF2ARFS_INVDIS

END MODULE invdist_wsf2arfs