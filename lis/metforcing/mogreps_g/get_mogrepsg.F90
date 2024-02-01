!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_mogrepsg
! \label{get_mogrepsg}
!
!
! !REVISION HISTORY:
! 26 Jan 2023: Yeosang Yoon, initial code
! 13 Mar 2023: Yeosang Yoon, update codes to fit new format
! 01 Jan 2024; Yeosang Yoon; update codes for precpi. bias-correction
!
! !INTERFACE:
subroutine get_mogrepsg(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use mogrepsg_forcingMod
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
!
! !DESCRIPTION:
!  Opens, reads, and interpolates MOGREPS-G forecast forcing.
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
  real              :: pcp1, pcp2
  integer           :: lead_time

  external :: get_mogrepsg_filename
  external :: read_mogrepsg

  ! MOGREPS-G cycles every 6 hours; ecch cycle provide up to 192 hours (8 days; 3-hour interval) forecast
  if (LIS_rc%ts .gt. 10800) then
     write(LIS_logunit,*) '[ERR] The model timestep is > forcing data timestep ...'
     write(LIS_logunit,*) '[ERR] LIS does not support this mode currently.'
     call LIS_endrun()
  endif

  openfile = 0

  if (LIS_rc%tscount(n) .eq. 1 .or. LIS_rc%rstflag(n) .eq. 1) then  !beginning of run
     LIS_rc%rstflag(n) = 0
  endif

  ! First timestep of run
  if (LIS_rc%tscount(n) .eq. 1 .or. LIS_rc%rstflag(n) .eq. 1) then
     ! Bookend-time record 1
     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = 0
     ss1 = 0
     ts1 = 0
     call LIS_tick(time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1)

     ! Bookend-time record 2
     yr2 = LIS_rc%yr    !next hour
     mo2 = LIS_rc%mo
     da2 = LIS_rc%da
     hr2 = 3
     mn2 = 0
     ss2 = 0
     ts2 = 0
     call LIS_tick(time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2)
     openfile=1
  endif

  ! 3 hourly interval
  fcsthr_intv = 3
  valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)

  if ((valid_hour == LIS_rc%hr .and. LIS_rc%mn == 0) .or. &
      openfile == 1)  then

     ! Forecast hour condition within each file:
     mogrepsg_struc(n)%fcst_hour = &
          mogrepsg_struc(n)%fcst_hour + fcsthr_intv

     ! Check if local forecast hour exceeds max grib file forecast hour:
     if (mogrepsg_struc(n)%fcst_hour > 195 ) then
        write(LIS_logunit,*) &
             "[ERR] MOGREPS-G Forecast hour has exceeded the grib file's final"
        write(LIS_logunit,*) &
              '  forecast hour (record). Run will end here for now ... '
        call LIS_endrun
     endif

     ! Update bookend-time record 2:
     if (LIS_rc%tscount(n) .ne. 1) then
        mogrepsg_struc(n)%fcsttime1 = mogrepsg_struc(n)%fcsttime2
        mogrepsg_struc(n)%metdata1(:,:,:) = &
             mogrepsg_struc(n)%metdata2(:,:,:)

        yr2 = LIS_rc%yr
        mo2 = LIS_rc%mo
        da2 = LIS_rc%da
        hr2 = valid_hour
        mn2 = fcsthr_intv * 60    ! Backward looking
        ss2 = 0
        ts2 = 0
        call LIS_tick(time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, &
             ts2)
     endif

     do m = 1, mogrepsg_struc(n)%max_ens_members

        ! Read in file contents:
        if (LIS_rc%tscount(n) == 1) then  ! Read in first two book-ends
           ferror = 0
           order = 1
           call get_mogrepsg_filename(mogrepsg_struc(n)%odir, &
                mogrepsg_struc(n)%init_yr, &
                mogrepsg_struc(n)%init_mo, &
                mogrepsg_struc(n)%init_da, &
                mogrepsg_struc(n)%init_hr, &
                0, m, fname)

           write(LIS_logunit,*)&
                '[INFO] Getting MOGREPS-G forecast file1 ... ', &
                trim(fname)
           call read_mogrepsg(n, m, findex, order, fname, ferror)
           if (ferror .ge. 1) mogrepsg_struc(n)%fcsttime1 = time1

           ferror = 0
           order = 2
           call get_mogrepsg_filename(mogrepsg_struc(n)%odir, &
                mogrepsg_struc(n)%init_yr, &
                mogrepsg_struc(n)%init_mo, &
                mogrepsg_struc(n)%init_da, &
                mogrepsg_struc(n)%init_hr, &
                mogrepsg_struc(n)%fcst_hour, m, fname)

           write(LIS_logunit,*) &
                '[INFO] Getting MOGREPS-G forecast file2 ... ', &
                trim(fname)
           call read_mogrepsg(n, m, findex, order, fname, ferror)
           if (ferror .ge. 1) mogrepsg_struc(n)%fcsttime2 = time2

           !only for T+0 due to mssing LW variable
           mogrepsg_struc(n)%metdata1(4,m,:) = &
                mogrepsg_struc(n)%metdata2(4,m,:)
        else
           ferror = 0
           order = 2
           ! met forcings except for pcp
           call get_mogrepsg_filename(mogrepsg_struc(n)%odir, &
                mogrepsg_struc(n)%init_yr, &
                mogrepsg_struc(n)%init_mo, &
                mogrepsg_struc(n)%init_da, &
                mogrepsg_struc(n)%init_hr, &
                mogrepsg_struc(n)%fcst_hour, m, fname)

           write(LIS_logunit,*) &
                '[INFO] Getting MOGREPS-G forecast file2 ... ', &
                trim(fname)
           call read_mogrepsg(n, m, findex, order, fname, ferror)
           if (ferror .ge. 1) mogrepsg_struc(n)%fcsttime2 = time2

           !only for T+141 due to mssing LW varaible
           if (mogrepsg_struc(n)%fcst_hour == 141) then
              mogrepsg_struc(n)%metdata2(4,m,:) = &
                   mogrepsg_struc(n)%metdata1(4,m,:)
           endif
        endif

        ! apply precipitation bias correction (cdf from difference bewteen NAPFA and MOGREPS-G)
        if (mogrepsg_struc(n)%bc == 1) then
           lead_time=floor((float(mogrepsg_struc(n)%fcst_hour))/24) + 1

           if (lead_time > 8) then
              lead_time = 8
           endif

           do t = 1, LIS_rc%ngrid(n)
              if (mogrepsg_struc(n)%metdata2(8,m,t) .ne. LIS_rc%udef) then
                 ! only for land pixels
                 if (mogrepsg_struc(n)%bc_param_a(t,lead_time) .ne. LIS_rc%udef) then
                    ! perform centering and scaling
                    pcp1= &
                         mogrepsg_struc(n)%metdata2(8,m,t) - &
                         mogrepsg_struc(n)%metdata1(8,m,t)
                    if (mogrepsg_struc(n)%bc_std(t,lead_time) .ne. 0) then
                       pcp2=(pcp1-mogrepsg_struc(n)%bc_mean(t,lead_time))/&
                          mogrepsg_struc(n)%bc_std(t,lead_time)
                    else
                       pcp2=pcp1
                    endif

                    ! apply cdf params
                    pcp2 = pcp2 * &
                         mogrepsg_struc(n)%bc_param_a(t,lead_time)+mogrepsg_struc(n)%bc_param_b(t,lead_time)
                    ! check for negative precipitation; if the corrected value has negative, keep the original value.
                    if (pcp2 >= 0) then
                       mogrepsg_struc(n)%pcp_bc(m,t) = pcp2
                    else
                       mogrepsg_struc(n)%pcp_bc(m,t) = pcp1
                    endif
                    ! additionally, avoid bias correction for values that are too samll to reduce abnormal noise.
                    if(pcp1 < 0.01) then
                       mogrepsg_struc(n)%pcp_bc(m,t) = pcp1
                    endif
                 else ! for water pixels
                    mogrepsg_struc(n)%pcp_bc(m,t) = &
                         mogrepsg_struc(n)%metdata2(8,m,t) - &
                         mogrepsg_struc(n)%metdata1(8,m,t)
                 endif
              endif
           enddo
        endif
     enddo
  endif
  openfile = 0

end subroutine get_mogrepsg

!BOP
!
! !ROUTINE: get_mogrepsg_filename
! \label{get_mogrepsg_filename}
!
! !INTERFACE:
subroutine get_mogrepsg_filename(rootdir, yr, mo, da, hr, fc_hr, &
     ens_id, filename)

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
!  This subroutine puts together MOGREPS-G file name for
!   operational products
!EOP
  character(8) :: ftime
  character(2) :: chr
  character(3) :: fchr
  character(2) :: ens
  character(len=36) :: fname

  write (UNIT=chr, FMT='(i2.2)') hr        ! cycle 00/06/12/18
  write (UNIT=fchr, FMT='(i3.3)') fc_hr    ! forecast time
  write (UNIT=ftime, FMT='(i4, i2.2, i2.2)') yr, mo, da

  fname = 'prods_op_mogreps-g_'

  !00/12z cycle memebers 00,01-17
  !06/18z cycle memebers 00,18-34
  if ((hr == 0) .or. (hr == 12)) then
       write (UNIT=ens, FMT='(i2.2)') ens_id-1   ! start 00, 01 - 17
  else
     if (ens_id == 1) then
        write (UNIT=ens, FMT='(i2.2)') ens_id-1  ! start 00
     else
        write (UNIT=ens, FMT='(i2.2)') ens_id+16 ! start 18-34
     endif
  endif

  filename = trim(rootdir)//'/'//ftime//chr//'/'//trim(fname) &
             //ftime//'_'//chr//'_'//ens//'_'//fchr//'.grib2'
end subroutine get_mogrepsg_filename

