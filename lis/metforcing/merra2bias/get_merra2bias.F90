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
! !ROUTINE: get_merra2bias
! \label{get_merra2bias}
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it)
! 07 Jan 2026: Kristen Whitney, initial code for using dynamic lapse rate
! 17 Jun 2026: Kristen Whitney, simplified metadata handling
!              and updated dynamic lapse-rate filename construction.
! 29 Jun 2026: Kristen Whitney, initial code (based on geos-itbias)
!
! !INTERFACE:
subroutine get_merra2bias(n,findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use merra2bias_forcingMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN, LIS_CONST_LAPSE_RATE

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly MERRA2bias forcing.
!
!  The MERRA2bias forcing data are organized into hourly files.
!  The data are considered valid at the mid-point of the hourly interval.
!
!  In general, metforcing readers read the forcing data before the
!  current time, referred to as bookend1, and after the current time,
!  referred to as bookend2.  Then the readers temporally interpolate
!  between bookend1 and bookend2.
!
!  Below are some examples to illustrate the timing logic of the
!  MERRA2bias reader.
!
!  \begin{verbatim}
!          ---*---|---*---|---*---|---*---|---*---|---*---|---*---|---*---
!  hour          21      22      23       0       1       2       3
!  hr_int         <---22--X--23---X--24---X---1---X---2---X---3--->
!
!  where:
!  hour is the hour UTC
!  hr_int is the hour-interval
!  * marks the valid point for the interval <--- hr_int --->
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing dataset index
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    call to advance or retract time
!  \item[merra2biasfiles](\ref{merra2biasfiles}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_merra2bias](\ref{read_merra2bias}) \newline
!    call to read the MERRA2bias data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: order
  integer           :: ferror
  character(len=LIS_CONST_PATH_LEN) :: merra2name
  integer           :: c,r,kk
  integer           :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer           :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8            :: time1,time2,timenow
  real              :: gmt1,gmt2
  real              :: ts1,ts2
  integer           :: gid
  integer           :: movetime ! Flag to move bookend2 files to bookend1
  character(len=LIS_CONST_PATH_LEN) :: lapseratefname
  character(len=13) :: fdate

  external :: merra2biasfiles
  external :: read_merra2bias
! _________________________________________________________
!

  if (LIS_rc%nts(n).gt.3600) then ! > 1-hr timestep
     write(LIS_logunit,*) '[ERR] When running LIS with MERRA2bias,'
     write(LIS_logunit,*) '[ERR] the clock should run with a time'
     write(LIS_logunit,*) '[ERR] step less than or equal to 1 hr.'
     call LIS_endrun()
  endif

  merra2bias_struc(n)%findtime1 = 0
  merra2bias_struc(n)%findtime2 = 0
  movetime = 0

!----------------------------------------------------------
! Determine current time
!----------------------------------------------------------
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  if (LIS_rc%ts.gt.3600.0) then
     write(LIS_logunit,*)                                          &
          '[ERR] The model timestep is > forcing data timestep'
     write(LIS_logunit,*)                                          &
          '[ERR] LIS does not support this mode currently'
     write(LIS_logunit,*) '[ERR] Program stopping ...'
     call LIS_endrun()
  endif

  if (mod(nint(LIS_rc%ts),3600).eq.0) then
     if (timenow.ge.merra2bias_struc(n)%merra2biastime2) then
        yr1 = LIS_rc%yr
        mo1 = LIS_rc%mo
        da1 = LIS_rc%da
        hr1 = LIS_rc%hr
        mn1 = 0
        ss1 = 0
        ts1 = -60*60
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2 = LIS_rc%yr     !next hour
        mo2 = LIS_rc%mo
        da2 = LIS_rc%da
        hr2 = LIS_rc%hr
        mn2 = 0
        ss2 = 0
        ts2 = 0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

        merra2bias_struc(n)%findtime1 = 1
        merra2bias_struc(n)%merra2biastime2 = time2
     endif
  else
     if (timenow.ge.merra2bias_struc(n)%merra2biastime2) then
        yr1 = LIS_rc%yr
        mo1 = LIS_rc%mo
        da1 = LIS_rc%da
        hr1 = LIS_rc%hr
        mn1 = 0
        ss1 = 0
        ts1 = 0
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        merra2bias_struc(n)%findtime1 = 1

        yr2 = LIS_rc%yr     !next hour
        mo2 = LIS_rc%mo
        da2 = LIS_rc%da
        hr2 = LIS_rc%hr
        mn2 = 0
        ss2 = 0
        ts2 = 60*60
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        merra2bias_struc(n)%merra2biastime2 = time2
     endif
  endif

! Beginning of the run
  if ((LIS_rc%tscount(n).eq.1).or.(LIS_rc%rstflag(n).eq.1)) then
     merra2bias_struc(n)%findtime1 = 1
     merra2bias_struc(n)%merra2biastime2 = time2
     LIS_rc%rstflag(n) = 0
  endif

! Read MERRA2bias - Bookend 1 files:
  if (merra2bias_struc(n)%findtime1.eq.1) then
     order = 1
     do kk = merra2bias_struc(n)%st_iterid,merra2bias_struc(n)%en_iterid
        call merra2biasfiles(n,kk,findex,          &
             merra2bias_struc(n)%merra2biasdir,    &
             doy1,yr1,mo1,da1,hr1,                 &
             merra2name)

        ! Get lapse rate filename
        write(fdate, fmt='(i4.4,"-",i2.2,"-",i2.2,"-",i2.2)') &
             yr1,mo1,da1,hr1
        lapseratefname = trim(merra2bias_struc(n)%dynlapseratedir)// &
             trim(merra2bias_struc(n)%dynlapseratepfx)// &
             fdate//trim(merra2bias_struc(n)%dynlapseratesfx)
        call read_merra2bias(n,order,mo1,findex,                       &
             merra2name, lapseratefname,                                 &
             merra2bias_struc(n)%merra2biasforc1(kk,:,:),ferror)
     enddo
     merra2bias_struc(n)%merra2biastime1 = time1
  endif

! Assign MERRA2bias forcing fields to two LIS time-interp placeholders:
  do r = 1,LIS_rc%lnr(n)
     do c = 1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r).ne.-1) then
           merra2bias_struc(n)%metdata(:,:,LIS_domain(n)%gindex(c,r)) = &
                merra2bias_struc(n)%merra2biasforc1(:,:,(c+(r-1)*LIS_rc%lnc(n)))
        endif
     enddo
  enddo

  ! Assign lapse rate values
  LIS_forc(n,findex)%lapseRate(:) = LIS_CONST_LAPSE_RATE

  if ((merra2bias_struc(n)%usedynlapserate.eq.1).and.                  &
       ((LIS_rc%met_ecor(findex).eq."lapse-rate").or.                  &
       (LIS_rc%met_ecor(findex).eq.                                    &
       "lapse-rate and slope-aspect").or.                              &
       (LIS_rc%met_ecor(findex).eq."micromet"))) then
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              gid = LIS_domain(n)%gindex(c,r)
              LIS_forc(n,findex)%lapseRate(gid) = &
                   merra2bias_struc(n)%lapserate1(gid)

              if(merra2bias_struc(n)%applydynlapseratecutoff.eq.1) then
                 if(LIS_forc(n,findex)%lapseRate(gid).gt. &
                      merra2bias_struc(n)%dynlapseratemaxcutoff) then
                    LIS_forc(n,findex)%lapseRate(gid) = &
                         merra2bias_struc(n)%dynlapseratemaxcutoff
                 elseif(LIS_forc(n,findex)%lapseRate(gid).lt. &
                      merra2bias_struc(n)%dynlapseratemincutoff) then
                    LIS_forc(n,findex)%lapseRate(gid) = &
                         merra2bias_struc(n)%dynlapseratemincutoff
                 endif
              endif
           endif
        enddo
     enddo
  endif
end subroutine get_merra2bias

!BOP
! !ROUTINE: merra2biasfiles
! \label{merra2biasfiles}
!
! !INTERFACE:
subroutine merra2biasfiles(n,kk,findex,merra2biasdir,doy,yr,mo,da,hr,    &
     merra2name)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod

  implicit none
! !ARGUMENTS:
  integer                       :: n
  integer                       :: kk
  integer                       :: findex
  character(len=*), intent(in)  :: merra2biasdir
  integer, intent(in)           :: doy,yr,mo,da,hr
  character(len=*), intent(out) :: merra2name

! !DESCRIPTION:
!   This subroutine puts together MERRA2bias file names for
!   hourly netCDF files.
!
!  The arguments are:
!  \begin{description}
!  \item[merra2biasdir]
!    Name of the MERRA2bias directory
!  \item[doy]
!    day of year
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour
!  \item[slvname]
!   name of the timestamped single level file
!  \item[flxname]
!   name of the timestamped flux file
!  \item[lfoname]
!   name of the timestamped land surface forcings file
!  \item[radname]
!   name of the timestamped radiation forcings file
!  \end{description}
!
!EOP

  character*4  :: cyear
  character*2  :: cmo,cdy,chr
  character*9 :: prefix

  write(unit=cyear,fmt='(i4.4)') yr
  write(unit=cmo,  fmt='(i2.2)') mo
  write(unit=cdy,  fmt='(i2.2)') da
  write(unit=chr,  fmt='(i2.2)') hr

  prefix = 'merra2bias_'
  if (yr.lt.1998) then
     write(LIS_logunit,*) '[ERR] MERRA2bias data starts 1 Jan 1998.'
     call LIS_endrun()
  endif

! Single level fields:
  merra2name = trim(merra2biasdir)//'/'//cyear//   &
       '/'//prefix//cyear//'-'//                 &
       cmo//'-'//cdy//'-'//chr//'.nc'
end subroutine merra2biasfiles

