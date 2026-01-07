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
! !ROUTINE: get_geositbias
! \label{get_geositbias}
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it) 
! 07 Jan 2026: Kristen Whitney, initial code for using dynamic lapse rate
!
! !INTERFACE:
      subroutine get_geositbias(n,findex)
! !USES:
      use LIS_coreMod
      use LIS_timeMgrMod
      use LIS_logMod
      use LIS_metforcingMod
      use geositbias_forcingMod
      use LIS_constantsMod, only : LIS_CONST_PATH_LEN, LIS_CONST_LAPSE_RATE

      implicit none

! !ARGUMENTS:
      integer, intent(in) :: n
      integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly GEOS-ITbias forcing.
!
!  The GEOS-ITbias forcing data are organized into hourly files.
!  The data are considered valid at the mid-point of the hourly interval.
!
!  In general, metforcing readers read the forcing data before the
!  current time, referred to as bookend1, and after the current time,
!  referred to as bookend2.  Then the readers temporally interpolate
!  between bookend1 and bookend2.
!
!  Below are some examples to illustrate the timing logic of the
!  GEOS-ITbias reader.
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
!  \item[geositbiasfiles](\ref{geositbiasfiles}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_geositbias](\ref{read_geositbias}) \newline
!    call to read the GEOS-ITbias data and perform spatial interpolation
!  \end{description}
!EOP
      integer           :: order
      integer           :: ferror
      character(len=LIS_CONST_PATH_LEN) :: geosname
      integer           :: c,f,r,kk
      integer           :: yr1,mo1,da1,hr1,mn1,ss1,doy1
      integer           :: yr2,mo2,da2,hr2,mn2,ss2,doy2
      real*8            :: time1,time2,timenow
      real              :: gmt1,gmt2
      real              :: ts1,ts2
      integer           :: gid
      integer           :: hr_int1,hr_int2
      integer           :: movetime ! Flag to move bookend2 files to bookend1
      character(len=LIS_CONST_PATH_LEN) :: lapseratefname
      character*20                      :: fdate
! _________________________________________________________
!

      if (LIS_rc%nts(n).gt.3600) then ! > 1-hr timestep
         write(LIS_logunit,*) '[ERR] When running LIS with GEOS-ITbias,'
         write(LIS_logunit,*) '[ERR] the clock should run with a time'
         write(LIS_logunit,*) '[ERR] step less than or equal to 1 hr.'
         call LIS_endrun()
      endif

      geositbias_struc(n)%findtime1 = 0
      geositbias_struc(n)%findtime2 = 0
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
         if (timenow.ge.geositbias_struc(n)%geositbiastime2) then 
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
         if (timenow.ge.geositbias_struc(n)%geositbiastime2) then 
            yr1 = LIS_rc%yr
            mo1 = LIS_rc%mo
            da1 = LIS_rc%da
            hr1 = LIS_rc%hr
            mn1 = 0
            ss1 = 0
            ts1 = 0
            call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
            geositbias_struc(n)%findtime1 = 1

            yr2 = LIS_rc%yr     !next hour
            mo2 = LIS_rc%mo
            da2 = LIS_rc%da
            hr2 = LIS_rc%hr
            mn2 = 0
            ss2 = 0
            ts2 = 60*60
            call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
            geositbias_struc(n)%geositbiastime2 = time2
         endif
      endif

! Beginning of the run
      if ((LIS_rc%tscount(n).eq.1).or.(LIS_rc%rstflag(n).eq.1)) then
         geositbias_struc(n)%findtime1 = 1
         geositbias_struc(n)%geositbiastime2 = time2
         LIS_rc%rstflag(n) = 0
      endif

! Read GEOS-ITbias - Bookend 1 files:
      if (geositbias_struc(n)%findtime1.eq.1) then
         order = 1
         do kk = geositbias_struc(n)%st_iterid,geositbias_struc(n)%en_iterid
            call geositbiasfiles(n,kk,findex,geositbias_struc(n)%geositbiasdir,    &
                           doy1,yr1,mo1,da1,hr1,                       &
                           geosname)
            
                   ! Get lapse rate filename
            write(fdate, fmt='(i4.4,"-",i2.2,"-",i2.2,"-",i2.2)') yr1,mo1,da1,hr1
            lapseratefname = trim(geositbias_struc(n)%dynlapseratedir)//&
                    trim(geositbias_struc(n)%dynlapseratepfx)//&
                    trim(fdate)//trim(geositbias_struc(n)%dynlapseratesfx)
            call read_geositbias(n,order,mo1,findex,                       &
                           geosname, lapseratefname,           &
                           geositbias_struc(n)%geositbiasforc1(kk,:,:),ferror)
         enddo
         geositbias_struc(n)%geositbiastime1 = time1
      endif

! Assign GEOS-ITbias forcing fields to two LIS time-interp placeholders:
      do r = 1,LIS_rc%lnr(n)
         do c = 1,LIS_rc%lnc(n)
            if (LIS_domain(n)%gindex(c,r).ne.-1) then
               geositbias_struc(n)%metdata1(:,:,LIS_domain(n)%gindex(c,r)) = &
                  geositbias_struc(n)%geositbiasforc1(:,:,(c+(r-1)*LIS_rc%lnc(n)))
               geositbias_struc(n)%metdata2(:,:,LIS_domain(n)%gindex(c,r)) = &
                  geositbias_struc(n)%geositbiasforc2(:,:,(c+(r-1)*LIS_rc%lnc(n)))
            endif
         enddo
      enddo

      ! Assign lapse rate values
      LIS_forc(n,findex)%lapseRate(:) = LIS_CONST_LAPSE_RATE

      if ((geositbias_struc(n)%usedynlapserate.eq.1).and.                      &
         ((LIS_rc%met_ecor(findex).eq."lapse-rate").or.                    &
          (LIS_rc%met_ecor(findex).eq."lapse-rate and slope-aspect").or.   &
          (LIS_rc%met_ecor(findex).eq."micromet"))) then
          do r=1,LIS_rc%lnr(n)
             do c=1,LIS_rc%lnc(n)
                if(LIS_domain(n)%gindex(c,r).ne.-1) then
                   gid = LIS_domain(n)%gindex(c,r)
                   LIS_forc(n,findex)%lapseRate(gid) = &
                      geositbias_struc(n)%lapserate1(gid)
                   
                   if(geositbias_struc(n)%applydynlapseratecutoff.eq.1) then
                      if(LIS_forc(n,findex)%lapseRate(gid).gt.geositbias_struc(n)%dynlapseratemaxcutoff) then
                         LIS_forc(n,findex)%lapseRate(gid) = geositbias_struc(n)%dynlapseratemaxcutoff
                      elseif(LIS_forc(n,findex)%lapseRate(gid).lt.geositbias_struc(n)%dynlapseratemincutoff) then
                         LIS_forc(n,findex)%lapseRate(gid) = geositbias_struc(n)%dynlapseratemincutoff
                      endif
                   endif
                endif
             enddo
          enddo
      endif
      end subroutine get_geositbias

!BOP
! !ROUTINE: geositbiasfiles
! \label{geositbiasfiles}
!
! !INTERFACE:
      subroutine geositbiasfiles(n,kk,findex,geositbiasdir,doy,yr,mo,da,hr,    &
                             geosname)

! !USES:
      use LIS_coreMod
      use LIS_logMod
      use LIS_timeMgrMod

      implicit none
! !ARGUMENTS:
      integer                       :: n
      integer                       :: kk
      integer                       :: findex
      character(len=*), intent(in)  :: geositbiasdir
      integer, intent(in)           :: doy,yr,mo,da,hr
      character(len=*), intent(out) :: geosname

! !DESCRIPTION:
!   This subroutine puts together GEOS-ITbias file names for
!   hourly netCDF files.
!
!  The arguments are:
!  \begin{description}
!  \item[geositbiasdir]
!    Name of the GEOS-ITbias directory
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

      prefix = 'geosbias_'
      if (yr.lt.1998) then
         write(LIS_logunit,*) '[ERR] GEOS-ITbias data starts 1 Jan 1998.'
         call LIS_endrun()
      endif


! Single level fields:
      geosname = trim(geositbiasdir)//'/'//cyear//   &
                '/'//prefix//cyear//'-'//       &
                cmo//'-'//cdy//'-'//chr//'.nc'
      end subroutine geositbiasfiles

