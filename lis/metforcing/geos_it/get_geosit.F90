!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_geosit
! \label{get_geosit}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 08 Dec 2015: James Geiger, update timing logic
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine get_geosit(n,findex)
! !USES:
      use LIS_coreMod
      use LIS_timeMgrMod
      use LIS_logMod
      use LIS_metforcingMod
      use geosit_forcingMod
      use LIS_constantsMod, only : LIS_CONST_PATH_LEN

      implicit none

! !ARGUMENTS:
      integer, intent(in) :: n
      integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly GEOS-IT forcing.
!
!  The GEOS-IT forcing data are organized into hourly files.
!  The data are considered valid at the mid-point of the hourly interval.
!
!  In general, metforcing readers read the forcing data before the
!  current time, referred to as bookend1, and after the current time,
!  referred to as bookend2.  Then the readers temporally interpolate
!  between bookend1 and bookend2.
!
!  Below are some examples to illustrate the timing logic of the
!  GEOS-IT reader.
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
!  \item[geositfiles](\ref{geositfiles}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_geosit](\ref{read_geosit}) \newline
!    call to read the GEOS-IT data and perform spatial interpolation
!  \end{description}
!EOP
      integer           :: order
      integer           :: ferror
      character(len=LIS_CONST_PATH_LEN) :: slvname,flxname,lfoname,radname
      integer           :: c,f,r,kk
      integer           :: yr1,mo1,da1,hr1,mn1,ss1,doy1
      integer           :: yr2,mo2,da2,hr2,mn2,ss2,doy2
      real*8            :: time1,time2,timenow
      real              :: gmt1,gmt2
      real              :: ts1,ts2
      integer           :: hr_int1,hr_int2
      integer           :: movetime ! Flag to move bookend2 files to bookend1

! _________________________________________________________
!

      if (LIS_rc%nts(n).gt.3600) then ! > 1-hr timestep
         write(LIS_logunit,*) '[ERR] When running LIS with GEOS-IT,'
         write(LIS_logunit,*) '[ERR] the clock should run with a time'
         write(LIS_logunit,*) '[ERR] step less than or equal to 1 hr.'
         call LIS_endrun()
      endif

      geosit_struc(n)%findtime1 = 0
      geosit_struc(n)%findtime2 = 0
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
         if (timenow.ge.geosit_struc(n)%geosittime2) then 
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
         endif
      else
         if (timenow.ge.geosit_struc(n)%geosittime2) then 
            yr1 = LIS_rc%yr
            mo1 = LIS_rc%mo
            da1 = LIS_rc%da
            hr1 = LIS_rc%hr
            mn1 = 0
            ss1 = 0
            ts1 = 0
            call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
            geosit_struc(n)%findtime1 = 1

            yr2 = LIS_rc%yr     !next hour
            mo2 = LIS_rc%mo
            da2 = LIS_rc%da
            hr2 = LIS_rc%hr
            mn2 = 0
            ss2 = 0
            ts2 = 60*60
            call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
            geosit_struc(n)%geosittime2 = time2
         endif
      endif

! Beginning of the run
      if ((LIS_rc%tscount(n).eq.1).or.(LIS_rc%rstflag(n).eq.1)) then
         geosit_struc(n)%findtime1 = 1
         geosit_struc(n)%geosittime2 = time2
         LIS_rc%rstflag(n) = 0
      endif

! Read GEOS-IT - Bookend 1 files:
      if (geosit_struc(n)%findtime1.eq.1) then
         order = 1
         do kk = geosit_struc(n)%st_iterid,geosit_struc(n)%en_iterid
            call geositfiles(n,kk,findex,geosit_struc(n)%geositdir,    &
                           doy1,yr1,mo1,da1,hr1,                       &
                           slvname,flxname,lfoname,radname)
            call read_geosit(n,order,mo1,findex,                       &
                           slvname,flxname,lfoname,radname,            &
                           geosit_struc(n)%geositforc1(kk,:,:),ferror)
         enddo
         geosit_struc(n)%geosittime1 = time1
      endif

! Assign GEOS-IT forcing fields to two LIS time-interp placeholders:
      do r = 1,LIS_rc%lnr(n)
         do c = 1,LIS_rc%lnc(n)
            if (LIS_domain(n)%gindex(c,r).ne.-1) then
               geosit_struc(n)%metdata1(:,:,LIS_domain(n)%gindex(c,r)) = &
                  geosit_struc(n)%geositforc1(:,:,(c+(r-1)*LIS_rc%lnc(n)))
               geosit_struc(n)%metdata2(:,:,LIS_domain(n)%gindex(c,r)) = &
                  geosit_struc(n)%geositforc2(:,:,(c+(r-1)*LIS_rc%lnc(n)))
            endif
         enddo
      enddo

      end subroutine get_geosit

!BOP
! !ROUTINE: geositfiles
! \label{geositfiles}
!
! !INTERFACE:
      subroutine geositfiles(n,kk,findex,geositdir,doy,yr,mo,da,hr,    &
                             slvname,flxname,lfoname,radname)

! !USES:
      use LIS_coreMod
      use LIS_logMod
      use LIS_timeMgrMod

      implicit none
! !ARGUMENTS:
      integer                       :: n
      integer                       :: kk
      integer                       :: findex
      character(len=*), intent(in)  :: geositdir
      integer, intent(in)           :: doy,yr,mo,da,hr
      character(len=*), intent(out) :: slvname
      character(len=*), intent(out) :: flxname
      character(len=*), intent(out) :: lfoname
      character(len=*), intent(out) :: radname

! !DESCRIPTION:
!   This subroutine puts together GEOS-IT file names for
!   hourly netCDF files.
!
!  The arguments are:
!  \begin{description}
!  \item[geositdir]
!    Name of the GEOS-IT directory
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
      character*2  :: cyr,cmo,cdy,chr
      character*16 :: prefix
      character*31 :: slv_spec,flx_spec,lfo_spec,rad_spec

      write(unit=cyear,fmt='(i4.4)') yr
      write(unit=cmo,  fmt='(i2.2)') mo
      write(unit=chr,  fmt='(i2.2)') hr
      write(unit=cdy,  fmt='(i2.2)') da

      prefix = 'd5294_geosit_jan'
      if (yr.lt.1998) then
         write(LIS_logunit,*) '[ERR] GEOS-IT data starts 1 Jan 1998.'
         call LIS_endrun()
      else
         cyr = '98'
      endif
      if (yr.ge.2008) then
         cyr = '08'
      endif
      if (yr.ge.2018) then
         cyr = '18'
      endif

      slv_spec = '.slv_tavg_1hr_glo_L576x361_slv.'
      flx_spec = '.flx_tavg_1hr_glo_L576x361_slv.'
      lfo_spec = '.lfo_tavg_1hr_glo_L576x361_slv.'
      rad_spec = '.rad_tavg_1hr_glo_L576x361_slv.'

! Single level fields:
      slvname = trim(geositdir)//'/'//prefix//cyr//'/diag/Y'//cyear//  &
                '/M'//cmo//'/'//prefix//cyr//slv_spec//cyear//'-'//    &
                cmo//'-'//cdy//'T'//chr//'30Z.nc4'

! Flux fields:
      flxname = trim(geositdir)//'/'//prefix//cyr//'/diag/Y'//cyear//  &
                '/M'//cmo//'/'//prefix//cyr//flx_spec//cyear//'-'//    &
                cmo//'-'//cdy//'T'//chr//'30Z.nc4'

! Land surface forcing level:
      lfoname = trim(geositdir)//'/'//prefix//cyr//'/diag/Y'//cyear//  &
                '/M'//cmo//'/'//prefix//cyr//lfo_spec//cyear//'-'//    &
                cmo//'-'//cdy//'T'//chr//'30Z.nc4'

! Radiation fields:
      radname = trim(geositdir)//'/'//prefix//cyr//'/diag/Y'//cyear//  &
                '/M'//cmo//'/'//prefix//cyr//rad_spec//cyear//'-'//    &
                cmo//'-'//cdy//'T'//chr//'30Z.nc4'

      end subroutine geositfiles

