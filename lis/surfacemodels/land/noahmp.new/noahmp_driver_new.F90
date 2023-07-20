!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Oct 15 2018 Shugong Wang started for implementing Noah-MP 4.0.1 based on the version of 3.6  
! Oct 15 2018 Zhuo Wang modifed for implementing Noah-MP 4.0.1 
! May 2023: Cenlin He modified for refactored Noah-MP v5 and later

#undef LIS_NoahMP_TEST 
! !INTERFACE
subroutine noahmp_driver_new(n, NoahmpIO, LISparam)   

  use LIS_coreMod, only    : LIS_rc
  use LIS_logMod,  only    : LIS_logunit, LIS_endrun
  use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
  use NoahmpIOVarType
  use LisNoahmpParamType 
  use NoahmpDriverMainMod

  implicit none 
  integer,                   intent(in)    :: n          ! nest id 
  type(NoahmpIO_type),       intent(inout) :: NoahmpIO   ! noahmp IO data type
  type(LisNoahmpParam_type), intent(in)    :: LISparam   ! lis/noahmp parameter

!--------------------------------------------------------------------------------
  ! external function
  real, external      :: month_d_new

  ! local variables 
  character(len=12)   :: nowdate
  integer :: k 

  ! Added by David Mocko on 11/19/2018
  logical :: Bondvillecheck
  integer :: i,local_hour
  integer :: locyr,locmo,locda,lochr,locmn,locss,locdoy
  real*8  :: loctime
  real    :: locgmt,change

  ! ---------------------------------------
   NoahmpIO%xland(1,1) = 1
   NoahmpIO%xice(1,1)  = 0
   ! Fraction of grid determining seaice (from WRF and HRLDAS)
   NoahmpIO%xice_threshold = 0.5

#ifndef WRF_HYDRO
   NoahmpIO%sfcheadrt(1,1) = 0.0
#endif
! SR = frozen precip fraction.  For offline, it is set to zero.
! If running coupled to WRF, then SR is set by the WRF model.
   if (NoahmpIO%IOPT_SNF .ne. 4) then
      NoahmpIO%SR(1,1) = 0.0
   else
      write(LIS_logunit,*) "[ERR] SR should be set by the WRF model."
      write(LIS_logunit,*) "[ERR] Code needs fixing.  Stopping LIS."
      call LIS_endrun
   endif

   ! dx is horizontal grid spacing. dx is not used past this point,
   ! but is used earlier when run_opt=5 (new groundwater scheme).
   NoahmpIO%dx = 0.0

  !!!!! print all the options not supported in offline mode 
  if (NoahmpIO%IOPT_SFC > 2) then
     stop "(opt_sfc == 3) and (opt_sfc == 4) are not for offline use"
  endif


  ! cosz, yearlen, and julian are calculated in subroutine calc_declin.
  ! Be careful here!!!, LIS uses GMT; the date for calc_declin_401 should
  ! be local time.  Longitude is need to convert GMT into local time!!!.
  ! This should only be done when using Bondville single_point forcing,
  ! as the forcing file uses local time instead of GMT.
  Bondvillecheck = .false.
  do i=1,LIS_rc%nmetforc
     if (trim(LIS_rc%metforc(i)).eq."Bondville") Bondvillecheck = .true.
  enddo

! For a true benchmark against the HRLDAS Noah-MP testcase
! from NCAR, set "change = -21600.0".  This line allows the code
! to _incorrectly_ run in the same way as the HRLDAS testcase,
! which runs on local time instead of on UTC time.  The changes
! is 6 hours * 3600.0 seconds per hour, to go from GMT to the
! Bondville location in Illinois, which uses Central time.
  change = 0.0
  if (Bondvillecheck) change = -21600.0
  locyr = LIS_rc%yr
  locmo = LIS_rc%mo
  locda = LIS_rc%da
  lochr = LIS_rc%hr
  locmn = LIS_rc%mn
  locss = LIS_rc%ss

  call LIS_date2time(loctime,locdoy,locgmt,locyr,locmo,locda,lochr,locmn,locss)
  call LIS_tick(loctime,locdoy,locgmt,locyr,locmo,locda,lochr,locmn,locss,change)

  write(nowdate,'(I4.4,4I2.2)') locyr, locmo, locda, lochr, locmn
  call calc_declin(nowdate(1:4)//"-"//nowdate(5:6)//"-"//nowdate(7:8)//"_"//nowdate(9:10)//":"//nowdate(11:12)//":00", &
       NoahmpIO%xlat(1,1), NoahmpIO%xlon(1,1), NoahmpIO%coszen(1,1), NoahmpIO%yearlen, NoahmpIO%julian)

  read(nowdate(1:4), '(I4)') NoahmpIO%YR

  if ((NoahmpIO%IOPT_DVEG .eq. 1).or.(NoahmpIO%IOPT_DVEG .eq. 6).or.(NoahmpIO%IOPT_DVEG .eq. 7)) then
    ! with dveg_opt==1/6/7, shdfac is fed directly to fveg
    NoahmpIO%vegfra(1,1) = month_d_new(NoahmpIO%shdfac_monthly(1,:,1), nowdate)
  else
    ! with dveg_opt==2/3/8, fveg is computed from lai and sai, and shdfac is unused
    ! with dveg_opt==4/5/9, fveg is set to the maximum shdfac, and shdfac is unused
    NoahmpIO%vegfra(1,1) = -1.E36   ! To match with HRLDAS initialization
  endif
  NoahmpIO%gvfmax(1,1) = maxval(NoahmpIO%shdfac_monthly(1,:,1))

  ! assign forcing variables 
  NoahmpIO%dz8w(1,2,1)     = NoahmpIO%dz8w(1,1,1)
  NoahmpIO%T_PHY(1,2,1)    = NoahmpIO%T_PHY(1,1,1)  
  NoahmpIO%P8W(1,2,1)      = NoahmpIO%P8W(1,1,1)
  NoahmpIO%QV_CURR(1,1,1)  = NoahmpIO%QV_CURR(1,1,1)/&
                             (1.0 - NoahmpIO%QV_CURR(1,1,1)) ! Convert specific humidity to water vapor mixing ratio
  NoahmpIO%QV_CURR(1,2,1)  = NoahmpIO%QV_CURR(1,1,1)
  NoahmpIO%U_PHY(1,2,1)    = NoahmpIO%U_PHY(1,1,1)
  NoahmpIO%V_PHY(1,2,1)    = NoahmpIO%V_PHY(1,1,1)
  NoahmpIO%QSFC(1,1)       = NoahmpIO%QV_CURR(1,1,1)
  NoahmpIO%SNOWBL(1,1)     = 0.0
  NoahmpIO%SR(1,1)         = 0.0                              ! Will only use component if opt_snf=4
  NoahmpIO%RAINCV(1,1)     = 0.0
  NoahmpIO%RAINNCV(1,1)    = NoahmpIO%RAINBL(1,1)
  NoahmpIO%RAINSHV(1,1)    = 0.0
  NoahmpIO%SNOWNCV(1,1)    = NoahmpIO%SNOWBL(1,1)
  NoahmpIO%GRAUPELNCV(1,1) = 0.0
  NoahmpIO%HAILNCV(1,1)    = 0.0

  ! If coupled to WRF, set these variables to realistic values,
  ! and then pass back to WRF after the call to noahmplsm_401.
  NoahmpIO%z0(1,1)     = 0.0
  NoahmpIO%znt(1,1)    = 0.0
  ! z0 and znt should be passed to WRF, if coupled. - dmm

  ! Code from module_NoahMP_hrldas_driver.F.  Initial guess only.
  if ((trim(LIS_rc%startcode).eq."coldstart").and.(NoahmpIO%itimestep.eq.1)) then
     NoahmpIO%eahxy(1,1) = NoahmpIO%P8W(1,1,1) * (NoahmpIO%QV_CURR(1,1,1)/(0.622+NoahmpIO%QV_CURR(1,1,1)))
     NoahmpIO%tahxy(1,1) = NoahmpIO%T_PHY(1,1,1)
     NoahmpIO%cmxy(1,1)  = 0.1
     NoahmpIO%chxy(1,1)  = 0.1
  endif
  
  ! main NoahMP driver physics
  call NoahmpDriverMain(NoahmpIO,LISparam)

end subroutine noahmp_driver_new

real function month_d_new(a12, nowdate) result (nowval)
  !
  ! Given a set of 12 values, taken to be valid on the fifteenth of each month (Jan through Dec)
  ! and a date in the form <YYYYMMDD[HHmmss]> ....
  ! 
  ! Return a value valid for the day given in <nowdate>, as an interpolation from the 12
  ! monthly values.
  !
  use kwm_date_utilities_new
  implicit none
  real, dimension(12), intent(in) :: a12 ! 12 monthly values, taken to be valid on the 15th of
  !                                      ! the month
  character(len=12), intent(in) :: nowdate ! Date, in the form <YYYYMMDD[HHmmss]>
  integer :: nowy, nowm, nowd
  integer :: prevm, postm
  real    :: factor
  integer, dimension(12) :: ndays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  !
  ! Handle leap year by setting the number of days in February for the year in question.
  !
  read(nowdate(1:8),'(I4,I2,I2)') nowy, nowm, nowd
  ndays(2) = nfeb(nowy)

  !
  ! Do interpolation between the fifteenth of two successive months.
  !
  if (nowd == 15) then
     nowval = a12(nowm)
     return
  else if (nowd < 15) then
     postm = nowm
     prevm = nowm - 1
     if (prevm == 0) prevm = 12
     factor = real(ndays(prevm)-15+nowd)/real(ndays(prevm))
  else if (nowd > 15) then
     prevm = nowm
     postm = nowm + 1
     if (postm == 13) postm = 1
     factor = real(nowd-15)/real(ndays(prevm))
  endif

  nowval = a12(prevm)*(1.0-factor) + a12(postm)*factor

end function month_d_new 

SUBROUTINE calc_declin ( nowdate, latitude, longitude, cosz, yearlen, julian)
  use kwm_date_utilities_new 
!---------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
  character(len=19), intent(in)  :: nowdate    ! YYYY-MM-DD_HH:mm:ss
  real,              intent(in)  :: latitude
  real,              intent(in)  :: longitude
  real,              intent(out) :: cosz
  integer,           intent(out) :: yearlen
  real,              intent(out) :: JULIAN

  REAL                           :: hrang
  real                           :: DECLIN
  real                           :: GMT
  real                           :: tloctim
  REAL                           :: OBECL
  REAL                           :: SINOB
  REAL                           :: SXLONG
  REAL                           :: ARG
  integer                        :: iyear
  integer                        :: iday
  integer                        :: ihour
  integer                        :: iminute
  integer                        :: isecond

  REAL, PARAMETER :: DEGRAD = 3.14159265/180.
  REAL, PARAMETER :: DPD    = 360./365.

  !
  ! Determine the number of days in the year
  !

  read(nowdate(1:4), '(I4)') iyear
  yearlen = 365
  if (mod(iyear,4) == 0) then
     yearlen = 366
     if (mod(iyear,100) == 0) then
        yearlen = 365
        if (mod(iyear,400) == 0) then
           yearlen = 366
           if (mod(iyear,3600) == 0) then
              yearlen = 365
           endif
        endif
     endif
  endif

  !
  ! Determine the Julian time (floating-point day of year).
  !

  call geth_idts(nowdate(1:10), nowdate(1:4)//"-01-01", iday)
  read(nowdate(12:13), *) ihour
  read(nowdate(15:16), *) iminute
  read(nowdate(18:19), *) isecond
  GMT = REAL(IHOUR) + IMINUTE/60.0 + ISECOND/3600.0
  JULIAN = REAL(IDAY) + GMT/24.

! for short wave radiation

  DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

  OBECL=23.5*DEGRAD
  SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

  IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
  IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
  ARG=SINOB*SIN(SXLONG)
  DECLIN=ASIN(ARG)

  TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! Local time in hours
  tloctim = AMOD(tloctim+24.0, 24.0)
  HRANG=15.*(TLOCTIM-12.)*DEGRAD
  COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)
  COSZ=MIN(COSZ,1.0);   !Added by kwH 3/1/16 to address floating point roundoff errors 
  COSZ=MAX(COSZ,-1.0);  !

!KWM   write(wrf_err_message,10)DECDEG/DEGRAD
!KWM10 FORMAT(1X,'*** SOLAR DECLINATION ANGLE = ',F6.2,' DEGREES.',' ***')
!KWM   CALL wrf_debug (50, wrf_err_message)

END SUBROUTINE calc_declin

