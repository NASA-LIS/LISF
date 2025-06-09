!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_galwemge
! \label{get_galwemge}
!
!
! !REVISION HISTORY:
! 15 Mar 2022: Yeosang Yoon, initial code
! 04 Apr 2023: Yeosang Yoon, Update code to fit new format
! 01 Jun 2025: Yeosang Yoon, update codes for new precip. bias-correction
!
! !INTERFACE:
subroutine get_galwemge(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use galwemge_forcingMod
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
!
! !DESCRIPTION:
!  Opens, reads, and interpolates GALWEM-GE forecast forcing.
!
!  At the beginning of a simulation, the code reads the most recent
!  past data (nearest 3-hour interval), and the nearest future data.
!  These two datasets are used to temporally interpolate the data to
!  the current model timestep.

!EOP
  integer           :: order, ferror, m, t
  character(len=LIS_CONST_PATH_LEN) :: fname
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2
  real              :: gmt1, gmt2
  real              :: ts1, ts2

  integer           :: valid_hour
  integer           :: fcsthr_intv
  integer           :: openfile

  ! precipitation bias correction
  real              :: pcp_tmp(10)
  integer           :: lead_time

  external :: get_galwemge_filename
  external :: read_galwemge
  external :: apply_cdf_correction_hybrid

  ! GALWEM-GE cycles every 12 hours; ecch cycle provide up to 384 hours
  ! (16 days) forecast; <=192 (every 3-hour); > 192 (every 6-hour)

  if (LIS_rc%ts.gt.10800) then
     write(LIS_logunit,*) &
          '[ERR] The model timestep is > forcing data timestep ...'
     write(LIS_logunit,*) &
          '[ERR] LIS does not support this mode currently.'
     call LIS_endrun()
  endif

  openfile=0

  if (LIS_rc%tscount(n).eq.1 .or. LIS_rc%rstflag(n).eq.1) then  !beginning of run
     LIS_rc%rstflag(n) = 0
  endif

  ! First timestep of run
  if (LIS_rc%tscount(n).eq.1 .or. LIS_rc%rstflag(n).eq.1) then
    ! Bookend-time record 1
     yr1 = LIS_rc%yr
     mo1=LIS_rc%mo
     da1=LIS_rc%da
     hr1=LIS_rc%hr
     mn1=0
     ss1=0
     ts1=0
     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

     ! Bookend-time record 2
     yr2=LIS_rc%yr    !next hour
     mo2=LIS_rc%mo
     da2=LIS_rc%da
     hr2=3            !3 hour in first 192 hours
     mn2=0
     ss2=0
     ts2=0
     call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
     openfile=1
  endif

  ! Determine valid times when forecasts are available to be read in:
  if (galwemge_struc(n)%fcst_hour < 192) then
     fcsthr_intv = 3
     valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
  elseif (galwemge_struc(n)%fcst_hour >= 192) then
     fcsthr_intv = 6
     valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
  endif

  if (galwemge_struc(n)%fcst_hour == 384) then
     !write(LIS_logunit,*) 'here: ', galwemge_struc(n)%fcst_hour
     galwemge_struc(n)%fcst_hour = 378 ! for the last forecast run
  endif

  if ((valid_hour==LIS_rc%hr .and. LIS_rc%mn==0) .or. openfile == 1)  then
     ! Forecast hour condition within each file:
     galwemge_struc(n)%fcst_hour = &
          galwemge_struc(n)%fcst_hour + fcsthr_intv

     ! Check if local forecast hour exceeds max grib file forecast hour:
     if (galwemge_struc(n)%fcst_hour > 384 ) then
        write(LIS_logunit,*) &
             "[ERR] GALWEM-GE Forecast hour has exceeded the grib file's final"
        write(LIS_logunit,*) &
             '[ERR] Forecast hour: ', galwemge_struc(n)%fcst_hour, &
             ' Run will end here for now ... '
        call LIS_endrun
     endif

     ! Update bookend-time record 2:
     if (LIS_rc%tscount(n).ne.1) then
        galwemge_struc(n)%fcsttime1 = galwemge_struc(n)%fcsttime2
        galwemge_struc(n)%metdata1 = galwemge_struc(n)%metdata2

        yr2=LIS_rc%yr
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=valid_hour
        mn2=fcsthr_intv*60    ! Backward looking
        ss2=0
        ts2=0
        call LIS_tick(time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, &
             ts2)
     endif

     do m=1, galwemge_struc(n)%max_ens_members

        ! Read in file contents:
        if (LIS_rc%tscount(n) == 1) then  ! Read in first two book-ends
           ferror=0
           order=1
           ! Obtain filenames
           call get_galwemge_filename(galwemge_struc(n)%odir, &
                galwemge_struc(n)%init_yr, &
                galwemge_struc(n)%init_mo, &
                galwemge_struc(n)%init_da, &
                galwemge_struc(n)%init_hr, &
                0, m, fname)

           write(LIS_logunit,*) &
                '[INFO] Getting GALWEM-GE forecast file1 ... ', &
                trim(fname)
           call read_galwemge(n, m, findex, order, fname, ferror)
           if (ferror.ge.1) galwemge_struc(n)%fcsttime1=time1

           ferror=0
           order=2
           call get_galwemge_filename(galwemge_struc(n)%odir, &
                galwemge_struc(n)%init_yr, &
                galwemge_struc(n)%init_mo, &
                galwemge_struc(n)%init_da, &
                galwemge_struc(n)%init_hr, &
                galwemge_struc(n)%fcst_hour, m, fname)

           write(LIS_logunit,*) &
                '[INFO] Getting GALWEM-GE forecast file2 ... ', &
                trim(fname)
           call read_galwemge(n, m, findex, order, fname, ferror)
           if (ferror.ge.1) galwemge_struc(n)%fcsttime2=time2
        else
           ferror=0
           order=2
           call get_galwemge_filename(galwemge_struc(n)%odir, &
                galwemge_struc(n)%init_yr, &
                galwemge_struc(n)%init_mo, &
                galwemge_struc(n)%init_da, &
                galwemge_struc(n)%init_hr, &
                galwemge_struc(n)%fcst_hour, m, fname)

           write(LIS_logunit,*) &
                '[INFO] Getting GALWEM-GE forecast file2 ... ', &
                trim(fname)
           call read_galwemge(n, m, findex, order, fname, ferror)
           if (ferror.ge.1) galwemge_struc(n)%fcsttime2=time2
        endif
     enddo

     ! apply precipitation bias correction (cdf from difference bewteen
     ! NAPFA and GALWEM-GE using MOGREPS-G parameters as a proxy)
     if (galwemge_struc(n)%bc == 1) then
        lead_time = &
             min(8, floor(real(galwemge_struc(n)%fcst_hour) / 24.0) + 1)

        do t = 1, LIS_rc%ngrid(n)
           if (galwemge_struc(n)%fcst_hour==3) then
              pcp_tmp = galwemge_struc(n)%metdata2(8,:,t)  !Inital time: no precp.
           else
              pcp_tmp = &
                   (galwemge_struc(n)%metdata2(8,:,t)-galwemge_struc(n)%metdata1(8,:,t))
           endif

           ! Apply CDF correction only over land
           if (galwemge_struc(n)%landmask(t) == 1.0) then
              call apply_cdf_correction_hybrid(pcp_tmp, &
                   galwemge_struc(n)%max_ens_members, &
                   galwemge_struc(n)%model_cdf(t,:,lead_time), &
                   galwemge_struc(n)%ref_cdf(t,:,lead_time), &
                   galwemge_struc(n)%percentiles, 101, &
                   galwemge_struc(n)%pcp_bc(:,t))
           else  ! For non-land points
              galwemge_struc(n)%pcp_bc(:,t) = pcp_tmp
           endif
        enddo
     endif
  endif
  openfile = 0

end subroutine get_galwemge

!BOP
!
! !ROUTINE: get_galwemge_filename
! \label{get_galwemge_filename}
!
! !INTERFACE:
subroutine get_galwemge_filename(rootdir, yr, mo, da, hr, fc_hr, ens_id, &
     filename)

  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_logMod, only: LIS_endrun

  implicit none
! !ARGUMENTS:
  character(len=*), intent(in)  :: rootdir
  integer,          intent(in)  :: yr,mo,da,hr
  integer,          intent(in)  :: fc_hr
  integer,          intent(in)  :: ens_id
  character(len=*), intent(out) :: filename
!
! !DESCRIPTION:
!  This subroutine puts together GALWEM-GE file name for
!   operational products
!EOP
  character(8) :: ftime
  character(2) :: chr
  character(3) :: fchr
  character(3) :: ens

  character(len=LIS_CONST_PATH_LEN) :: fname

  write (UNIT=chr, FMT='(i2.2)') hr
  write (UNIT=fchr, FMT='(i3.3)') fc_hr
  write (UNIT=ftime, FMT='(i4, i2.2, i2.2)') yr, mo, da

  fname = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GE-MEMB'

  !TODO: need to check, memebers 00,28-36
  if (ens_id == 1) then
     write (UNIT=ens, FMT='(i3.3)') ens_id-1  ! start 00
  else
     write (UNIT=ens, FMT='(i3.3)') ens_id+26 ! start 28-36
  endif

  filename = trim(rootdir)//'/'//ftime//'/'//'member'//ens//'/'//&
             trim(fname)//ens//'_GR.C20KM_AR.GLOBAL_DD.'//   &
             ftime//'_CY.'//chr//'_FH.'//fchr//'_DF.GR2'
end subroutine get_galwemge_filename

subroutine apply_cdf_correction_hybrid(val_ens, nens, model_cdf, &
     ref_cdf, percentiles, npercentile, corrected_ens)
  implicit none

  ! Inputs
  integer, intent(in) :: nens, npercentile
  real, intent(in)    :: val_ens(nens)
  real, intent(in)    :: model_cdf(npercentile), ref_cdf(npercentile), &
       percentiles(npercentile)

  ! Output
  real, intent(out)   :: corrected_ens(nens)

  ! Internal
  real :: model_cdf_smooth(npercentile), ref_cdf_smooth(npercentile)
  real :: val_mean, inv_cdf, model_val_at_inv_cdf, ref_val_at_inv_cdf
  real :: corr_mean, spread_corr(nens)
  real :: sigma, scaling, bias_corr
  real :: model_std, ref_std, ratio
  integer :: i
  real :: upper_clip, lower_clip

  external :: gaussian_smooth
  external :: interp1_linear

  ! 0. Check for degenerate CDFs (all values same or not strictly increasing)
  if (maxval(model_cdf) - minval(model_cdf) < 1e-6 .or. maxval(ref_cdf) - minval(ref_cdf) < 1e-6) then
    corrected_ens = val_ens  ! skip correction
    return
  end if

  ! 1. Gussian smoothing
  sigma = 1.0
  call gaussian_smooth(model_cdf, npercentile, model_cdf_smooth, sigma)
  call gaussian_smooth(ref_cdf, npercentile, ref_cdf_smooth, sigma)

  ! 2. Validate monotonicity (safety check)
  do i = 1, npercentile - 1
     if (model_cdf_smooth(i+1) <= model_cdf_smooth(i) .or. &
          ref_cdf_smooth(i+1) <= ref_cdf_smooth(i)) then
      corrected_ens = val_ens  ! skip correction
      return
    end if
  end do

  ! 3. Compute val_mean
  val_mean = sum(val_ens) / real(nens)

  ! 4. Invert model CDF (val_mean → inv_cdf percentile)
  call interp1_linear(model_cdf_smooth, percentiles, npercentile, &
       val_mean, inv_cdf)

  !!5. Clamp inv_cdf to valid percentile range (e.g., 0–100)
  inv_cdf = max(min(inv_cdf, 100.0), 0.0)

  ! 6. Map inv_cdf to ref/model values
  call interp1_linear(percentiles, ref_cdf_smooth, npercentile, &
       inv_cdf, ref_val_at_inv_cdf)
  call interp1_linear(percentiles, model_cdf_smooth, npercentile, &
       inv_cdf, model_val_at_inv_cdf)

  ! 7. Bias correction of mean
  corr_mean = val_mean + (ref_val_at_inv_cdf - model_val_at_inv_cdf)

  ! 8. Estimate spread from model and reference
  ! --- Compute spread (half-range approximation) ---
  model_std = 0.5 * (model_cdf_smooth(npercentile) - model_cdf_smooth(1))
  ref_std   = 0.5 * (ref_cdf_smooth(npercentile) - ref_cdf_smooth(1))

  ! --- Adjust spread scaling: allow mild expansion ---
  if (model_std > 1e-6) then
    !ratio = ref_std / model_std
    ratio = min(ref_std / model_std, 2.0)
    scaling = 1.0 + 0.5 * tanh(ratio - 1.0)  ! more modest scaling

    ! Allow range: [0.5, 1.5]
    scaling = min(max(scaling, 0.5), 1.5)
  else
    scaling = 1.0
  end if

  ! 9. Apply scaling around val_mean and reconstruct corrected_ens
  do i = 1, nens
    spread_corr(i) = (val_ens(i) - val_mean) * scaling
    corrected_ens(i) = corr_mean + spread_corr(i)
  end do

  ! 10. Re-center corrected ensemble
  bias_corr = sum(corrected_ens) / real(nens) - corr_mean
  do i = 1, nens
    corrected_ens(i) = corrected_ens(i) - bias_corr
  end do

  ! 11. Cap values using reference CDF range (e.g., P0–P100)
  lower_clip = ref_cdf_smooth(1)
  upper_clip = ref_cdf_smooth(npercentile)
  do i = 1, nens
    corrected_ens(i) = max(min(corrected_ens(i), upper_clip), lower_clip)
  end do

end subroutine apply_cdf_correction_hybrid

subroutine interp1_linear(x, y, n, xq, yq)
  use LIS_logMod

  implicit none
  integer, intent(in) :: n
  real, intent(in)    :: x(n), y(n), xq
  real, intent(out)   :: yq
  integer :: i
  real :: h, t

  ! Check monotonicity
  do i = 1, n - 1
    if (x(i+1) <= x(i)) then
       write(LIS_logunit,*) &
            '[ERR] in interp1_linear: x must be strictly increasing.'
       write(LIS_logunit,*) 'x(', i, ') = ', x(i), ', x(', i+1, ') = ', &
            x(i+1)
       call LIS_endrun
    end if
  end do

  ! Handle left extrapolation
  if (xq <= x(1)) then
     h = x(2) - x(1)
     yq = y(1) + (y(2) - y(1)) / h * (xq - x(1))
     return
  end if

  ! Handle right extrapolation
  if (xq >= x(n)) then
     h = x(n) - x(n-1)
     yq = y(n-1) + (y(n) - y(n-1)) / h * (xq - x(n-1))
     return
  end if

  ! Locate the interval
  do i = 1, n - 1
     if (xq >= x(i) .and. xq <= x(i+1)) then
        h = x(i+1) - x(i)
        t = (xq - x(i)) / h
        yq = (1.0 - t) * y(i) + t * y(i+1)
        return
     end if
  end do

  ! Fallback (should not occur)
  yq = y(1)
end subroutine interp1_linear

subroutine gaussian_smooth(input, n, output, sigma)
  implicit none
  integer, intent(in) :: n
  real, intent(in)    :: input(n), sigma
  real, intent(out)   :: output(n)

  ! Locals
  integer :: i, j, w
  real :: weight_sum, dist, wgt

  w = 3  ! smoothing window: [-3, 3]

  do i = 1, n
     output(i) = 0.0
     weight_sum = 0.0
     do j = max(1, i - w), min(n, i + w)
        dist = real(j - i)
        wgt = exp(-0.5 * (dist / sigma) ** 2)
        output(i) = output(i) + input(j) * wgt
        weight_sum = weight_sum + wgt
     end do
     output(i) = output(i) / weight_sum
  end do
end subroutine gaussian_smooth

