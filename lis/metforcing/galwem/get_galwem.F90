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
! !ROUTINE: get_galwem
! \label{get_galwem}
!
!
! !REVISION HISTORY:
! 15 Mar 2022: Yeosang Yoon, initial code
! 08 Sep 2022; Yeosang Yoon, Add codes to read GALWEM 25 DEG dataset
! 25 May 2023: Eric Kemp, extended GALWEM 0.25 deg to 240 hours.
! 01 Jun 2023: Eric Kemp, fixes to time code to support runs other than
!              00Z, to correct time intervals between GALWEM-GD files,
!              and to stop searching for (nonexistent) GALWEM-GD file
!              when LIS reaches it's end time.
! !INTERFACE:
subroutine get_galwem(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use galwem_forcingMod
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates GALWEM forecast forcing.

!EOP
  integer           :: order,ferror
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

  ! GALWEM cycles every 6 hours; each cycle provide up to 168 hours (7 days) forecast for GALWEM-17km;
  ! each cycle provide up to 240 hours (10 days) forecast for GALWEM-25deg;
  ! <=42 (every 1-hour); > 42 (every 3-hour)

  if(LIS_rc%ts.gt.10800) then
     write(LIS_logunit,*) '[WARN] The model timestep is > forcing data timestep ...'
     write(LIS_logunit,*) '[WARN] LIS does not support this mode currently.'
     call LIS_endrun()
  endif

  ! EMK...Return if LIS has reached the end time (meaning a new GALWEM
  ! file does not need to be read).
  if (LIS_rc%yr == LIS_rc%eyr .and. &
       LIS_rc%mo == LIS_rc%emo .and. &
       LIS_rc%da == LIS_rc%eda .and. &
       LIS_rc%hr == LIS_rc%ehr .and. &
       LIS_rc%mn == LIS_rc%emn .and. &
       LIS_rc%ss == LIS_rc%ess) then
     write(LIS_logunit,*) &
          '[INFO] LIS run has reached end time, will not read more GALWEM'
     return
  end if
  
  openfile=0

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1) then  !beginning of run
     LIS_rc%rstflag(n) = 0
  endif

  ! First timestep of run
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1) then
    ! Bookend-time record 1
     yr1 = LIS_rc%syr
     mo1=LIS_rc%smo
     da1=LIS_rc%sda
     hr1=LIS_rc%shr
     mn1=0
     ss1=0
     ts1=0
     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

     ! Bookend-time record 2
     yr2=LIS_rc%syr
     mo2=LIS_rc%smo
     da2=LIS_rc%sda
     hr2=LIS_rc%shr
     mn2=0
     ss2=0
     ts2=3600 ! Advance 1 hour
     call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

     !movetime = 1
     openfile=1

     !write(LIS_logunit,*)'EMK: time1,yr1,mo1,da1,hr1,ts1 = ', &
     !     time1,yr1,mo1,da1,hr1,ts1
     !write(LIS_logunit,*)'EMK: time2,yr2,mo2,da2,hr2,ts2 = ', &
     !     time2,yr2,mo2,da2,hr2,ts2

  endif

  ! Determine valid times when forecasts are available to be read in:
  ! EMK...Revised based on sample forecast runs for 20230523.
  if(galwem_struc(n)%fcst_hour < 42) then
     fcsthr_intv = 1
     valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
  elseif(galwem_struc(n)%fcst_hour >= 42 .and. &
       galwem_struc(n)%fcst_hour < 168) then
     fcsthr_intv = 3
     valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
  else
     fcsthr_intv = 6
     valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
  endif

  if((valid_hour==LIS_rc%hr .and. LIS_rc%mn==0) .or. &
      openfile == 1)  then

     ! Forecast hour condition within each file:
     if(galwem_struc(n)%fcst_hour < 42) then
        galwem_struc(n)%fcst_hour = galwem_struc(n)%fcst_hour + fcsthr_intv
     elseif(galwem_struc(n)%fcst_hour >= 42) then
        galwem_struc(n)%fcst_hour = galwem_struc(n)%fcst_hour + fcsthr_intv
     endif

     ! Check if local forecast hour exceeds max grib file forecast hour (GALWEM-17km):
     if(galwem_struc(n)%resol == 17) then
        if(galwem_struc(n)%fcst_hour > 168 ) then
           write(LIS_logunit,*) &
                 "[INFO] GALWEM Forecast hour has exceeded the grib file's final"
           write(LIS_logunit,*) &
                 "  forecast hour (record). Run will end here for now ... "
           call LIS_endrun
        endif
     endif
     ! Check if local forecast hour exceeds max grib file forecast hour (GALWEM-25deg):
     if(galwem_struc(n)%resol == 25) then
        if(galwem_struc(n)%fcst_hour > 240 ) then
           write(LIS_logunit,*) &
                 "[INFO] GALWEM Forecast hour has exceeded the grib file's final"
           write(LIS_logunit,*) &
                 "  forecast hour (record). Run will end here for now ... "
           call LIS_endrun
        endif
     endif
   
     ! Update bookend-time record 2:
     if(LIS_rc%tscount(n).ne.1) then
        galwem_struc(n)%fcsttime1=galwem_struc(n)%fcsttime2
        galwem_struc(n)%metdata1=galwem_struc(n)%metdata2

        yr2=LIS_rc%syr
        mo2=LIS_rc%smo
        da2=LIS_rc%sda
        hr2=LIS_rc%shr
        mn2=0
        ss2=0
        ts2=3600 * galwem_struc(n)%fcst_hour
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
      endif

      ! Read in file contents:
      if(LIS_rc%tscount(n) == 1) then  ! Read in first two book-ends 

         !write(LIS_logunit,*)'EMK: galwem_struc(n)%init_yr = ', &
         !     galwem_struc(n)%init_yr
         !write(LIS_logunit,*)'EMK: galwem_struc(n)%init_mo = ', &
         !     galwem_struc(n)%init_mo
         !write(LIS_logunit,*)'EMK: galwem_struc(n)%init_da = ', &
         !     galwem_struc(n)%init_da
         !write(LIS_logunit,*)'EMK: galwem_struc(n)%init_hr = ', &
         !     galwem_struc(n)%init_hr
         !write(LIS_logunit,*)'EMK: galwem_struc(n)%fcst_hour = ', &
         ! galwem_struc(n)%fcst_hour

         ferror=0
         order=1   
         call getGALWEMfilename(n,galwem_struc(n)%odir,galwem_struc(n)%init_yr,&
              galwem_struc(n)%init_mo,galwem_struc(n)%init_da,galwem_struc(n)%init_hr,&
              0,fname)

         write(LIS_logunit,*)'[INFO] Getting GALWEM forecast file1 ... ',trim(fname)
         call read_galwem(n, findex, order, fname, ferror)
         if(ferror.ge.1) galwem_struc(n)%fcsttime1=time1
                 
         ferror=0
         order=2
         call getGALWEMfilename(n,galwem_struc(n)%odir,galwem_struc(n)%init_yr,&
              galwem_struc(n)%init_mo,galwem_struc(n)%init_da,galwem_struc(n)%init_hr,&
              galwem_struc(n)%fcst_hour,fname)

         write(LIS_logunit,*)'[INFO] Getting GALWEM forecast file2 ... ',trim(fname)
         call read_galwem(n, findex, order, fname, ferror)
         if(ferror.ge.1) galwem_struc(n)%fcsttime2=time2
      else
         ferror=0
         order=2
         call getGALWEMfilename(n,galwem_struc(n)%odir,galwem_struc(n)%init_yr,&
              galwem_struc(n)%init_mo,galwem_struc(n)%init_da,galwem_struc(n)%init_hr,&
              galwem_struc(n)%fcst_hour,fname)

         write(LIS_logunit,*)'[INFO] Getting GALWEM forecast file2 ... ',trim(fname)
         call read_galwem(n, findex, order, fname, ferror)
         if(ferror.ge.1) galwem_struc(n)%fcsttime2=time2
      endif
  endif
  openfile=0

  !write(LIS_logunit,*)"EMK: galwem_struc(n)%fcsttime1 = ", &
  !     galwem_struc(n)%fcsttime1
  !write(LIS_logunit,*)"EMK: galwem_struc(n)%fcsttime2 = ", &
  !     galwem_struc(n)%fcsttime2

end subroutine get_galwem

!BOP
!
! !ROUTINE: getGALWEMfilename
! \label{getGALWEMfilename}
!
! !INTERFACE:
subroutine getGALWEMfilename(n,rootdir,yr,mo,da,hr,fc_hr,filename)

  use LIS_logMod, only: LIS_logunit, LIS_endrun
  use galwem_forcingMod

  implicit none
! !ARGUMENTS:
  integer,          intent(in)  :: n
  character(len=*), intent(in)  :: rootdir
  integer,          intent(in)  :: yr,mo,da,hr
  integer,          intent(in)  :: fc_hr
  character(len=*), intent(out) :: filename

!EOP
  character(8) :: ftime
  character(2) :: chr
  character(3) :: fchr

  character(len=100) :: fname

  write (UNIT=chr, FMT='(i2.2)') hr
  write (UNIT=fchr, FMT='(i3.3)') fc_hr
 
  ! GALWEM-17km 
  if(galwem_struc(n)%resol == 17) then   
     fname = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GD_GR.C17KM_AR.GLOBAL_DD.'
  endif

  ! GALWEM-25deg
  if(galwem_struc(n)%resol == 25) then
     fname = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GD_GR.C0P25DEG_AR.GLOBAL_DD.'
  endif

  write (UNIT=ftime, FMT='(i4, i2.2, i2.2)') yr, mo, da

  filename = trim(rootdir) // '/' // ftime // '/' // &
             trim(fname) // ftime // '_CY.' // chr // '_FH.' // fchr // '_DF.GR2'
end subroutine getGALWEMfilename


