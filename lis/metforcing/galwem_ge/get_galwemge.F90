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
! !ROUTINE: get_galwemge
! \label{get_galwemge}
!
!
! !REVISION HISTORY:
! 15 Mar 2022: Yeosang Yoon, initial code
! 04 Apr 2023: Yeosang Yoon, Update code to fit new format
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
  integer           :: order, ferror, m
  character(len=LIS_CONST_PATH_LEN) :: fname
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2
  real              :: gmt1, gmt2
  real              :: ts1, ts2
  integer           :: fc_hr

  integer           :: hr_int1, hr_int2
  integer           :: valid_hour
  integer           :: fcsthr_intv
  integer           :: fcst_hour
  integer           :: openfile

  ! GALWEM-GE cycles every 12 hours; ecch cycle provide up to 384 hours (16 days) forecast; 
  ! <=192 (every 3-hour); > 192 (every 6-hour)

  if(LIS_rc%ts.gt.10800) then
     write(LIS_logunit,*) '[WARN] The model timestep is > forcing data timestep ...'
     write(LIS_logunit,*) '[WARN] LIS does not support this mode currently.'
     call LIS_endrun()
  endif

  openfile=0

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1) then  !beginning of run
     LIS_rc%rstflag(n) = 0
  endif

  ! First timestep of run
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1) then
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
  if(galwemge_struc(n)%fcst_hour < 192) then
     fcsthr_intv = 3
     valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
  elseif(galwemge_struc(n)%fcst_hour >= 192) then
     fcsthr_intv = 6
     valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
  endif

  if((valid_hour==LIS_rc%hr .and. LIS_rc%mn==0) .or. &
      openfile == 1)  then

     ! Forecast hour condition within each file:
      galwemge_struc(n)%fcst_hour = galwemge_struc(n)%fcst_hour + fcsthr_intv
 
     ! Check if local forecast hour exceeds max grib file forecast hour:
     if(galwemge_struc(n)%fcst_hour > 384 ) then
        write(LIS_logunit,*) &
              "[INFO] GALWEM-GE Forecast hour has exceeded the grib file's final"
        write(LIS_logunit,*) &
              '  forecast hour (record). Run will end here for now ... '
        call LIS_endrun
     endif
  
     ! Update bookend-time record 2:
     if(LIS_rc%tscount(n).ne.1) then
        galwemge_struc(n)%fcsttime1=galwemge_struc(n)%fcsttime2
        galwemge_struc(n)%metdata1=galwemge_struc(n)%metdata2

        yr2=LIS_rc%yr
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=valid_hour
        mn2=fcsthr_intv*60    ! Backward looking
        ss2=0
        ts2=0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
     endif

     do m=1,galwemge_struc(n)%max_ens_members    

        ! Read in file contents:
        if(LIS_rc%tscount(n) == 1) then  ! Read in first two book-ends 
           ferror=0
           order=1
          ! Obtain filenames   
           call get_galwemge_filename(galwemge_struc(n)%odir,galwemge_struc(n)%init_yr,&
                galwemge_struc(n)%init_mo,galwemge_struc(n)%init_da,galwemge_struc(n)%init_hr,&
                0,m,fname)

           write(LIS_logunit,*)'[INFO] Getting GALWEM-GE forecast file1 ... ',trim(fname)
           call read_galwemge(n, m, findex, order, fname, ferror)
           if(ferror.ge.1) galwemge_struc(n)%fcsttime1=time1
                 
           ferror=0
           order=2
           call get_galwemge_filename(galwemge_struc(n)%odir,galwemge_struc(n)%init_yr,&
                galwemge_struc(n)%init_mo,galwemge_struc(n)%init_da,galwemge_struc(n)%init_hr,&
                galwemge_struc(n)%fcst_hour,m,fname)

           write(LIS_logunit,*)'[INFO] Getting GALWEM-GE forecast file2 ... ',trim(fname)
           call read_galwemge(n, m, findex, order, fname, ferror)
           if(ferror.ge.1) galwemge_struc(n)%fcsttime2=time2
        else
           ferror=0
           order=2
           call get_galwemge_filename(galwemge_struc(n)%odir,galwemge_struc(n)%init_yr,&
                galwemge_struc(n)%init_mo,galwemge_struc(n)%init_da,galwemge_struc(n)%init_hr,&
                galwemge_struc(n)%fcst_hour,m,fname)

           write(LIS_logunit,*)'[INFO] Getting GALWEM-GE forecast file2 ... ',trim(fname)
           call read_galwemge(n, m, findex, order, fname, ferror)
           if(ferror.ge.1) galwemge_struc(n)%fcsttime2=time2
        endif
     enddo
  endif
  openfile=0

end subroutine get_galwemge

!BOP
!
! !ROUTINE: get_galwemge_filename
! \label{get_galwemge_filename}
!
! !INTERFACE:
subroutine get_galwemge_filename(rootdir,yr,mo,da,hr,fc_hr,ens_id,filename)

  use LIS_logMod, only: LIS_logunit, LIS_endrun
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

  character(len=37) :: fname

  write (UNIT=chr, FMT='(i2.2)') hr
  write (UNIT=fchr, FMT='(i3.3)') fc_hr
  write (UNIT=ftime, FMT='(i4, i2.2, i2.2)') yr, mo, da

  fname = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GE-MEMB'

  !TODO: need to check, 12z cycle memebers 00,28-44
  if (hr == 0) then
       write (UNIT=ens, FMT='(i3.3)') ens_id-1   ! start 00, 01 - 17
  else
     if (ens_id == 1) then
        write (UNIT=ens, FMT='(i3.3)') ens_id-1  ! start 00
     else
        write (UNIT=ens, FMT='(i3.3)') ens_id+26 ! start 28-44
     endif
  endif

  filename = trim(rootdir)//'/'//ftime//'T'//chr//'00Z'//'/'// &
             trim(fname)//ens//'_GR.C0P5DEG_AR.GLOBAL_DD.'//   &
             ftime//'_CY.'//chr//'_FH.'//fchr//'_DF.GR2'
end subroutine get_galwemge_filename


