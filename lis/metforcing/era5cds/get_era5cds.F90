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
! !ROUTINE: get_era5cds
! \label{get_era5cds}
!
!
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code
! 04 Mar 2025: Hiroko Beudoing, adopted ERA5 routines for the public CDS
!                               data format
!
! !INTERFACE:
subroutine get_era5cds(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use era5cds_forcingMod
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly ERA5 forcing from
!  the Climate Data Store.
!
!  The ERA5 forcing data are organized into monthly files, where each
!  file contains 24 one-hourly records over days in a month per instantaneous 
!  or accumulation forcing fields. Entire month of data is read in at the
!  beginning of month.   
!
!  In general, metforcing readers read the forcing data before the current
!  time, referred to as bookend1, and after the current time, referred to as
!  bookend2.  Then the readers temporally interpolate between bookend1 and
!  bookend2.  
!
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
!  \item[era5cdsfiles](\ref{era5cdsfiles}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_era5cds](\ref{read_era5cds}) \newline
!    call to read the ERA5 data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: order
  integer           :: ferror
  character(len=LIS_CONST_PATH_LEN) :: instfilename, lmlfilename
  character(len=LIS_CONST_PATH_LEN) :: avgfilename, prevavgfilename
  integer           :: c, r,kk,f,try
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2, timenow
  real*8            :: dtime1, dtime2
  real              :: gmt1, gmt2
  real              :: ts1, ts2

  integer           :: movetime  ! Flag to move bookend2 files to bookend1
  logical           :: retrieve_file

! _________________________________________________________

  if( LIS_rc%nts(n).gt.3600 ) then   ! > 1-hr timestep
     write(LIS_logunit,*) '[ERR] When running LIS with ERA5CDS, the clock '
     write(LIS_logunit,*) '[ERR] should run with a timestep less than or '
     write(LIS_logunit,*) '[ERR] equal to one hour.',LIS_rc%nts(n)
     call LIS_endrun()
  endif
  if(LIS_rc%ts.gt.3600) then 
     write(LIS_logunit,*) '[ERR] The model timestep is > forcing data timestep'
     write(LIS_logunit,*) '[ERR] LIS does not support this mode currently'
     write(LIS_logunit,*) '[ERR] Program stopping ...',LIS_rc%ts
     call LIS_endrun()
  endif

  era5cds_struc(n)%findtime1 = 0
  era5cds_struc(n)%findtime2 = 0
  movetime = 0
  retrieve_file = .false.

!=== Determine Required ERA5 Data Times (The previous hour and the future hour)
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 
  if( timenow < era5cds_struc(n)%validstart ) then
     write(LIS_logunit,*) '[ERR] ERA5CDS forecast begins at 7z on 1940/01/01'
     write(LIS_logunit,*) '[ERR] Starting hour should be set to 7z '
     write(LIS_logunit,*) '[ERR] in lis.config file. Stopping ... '
     call LIS_endrun()
  endif

  yr1 = LIS_rc%yr  !previous hour
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=LIS_rc%yr    !next hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=LIS_rc%hr
  mn2=0
  ss2=0
  ts2=60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  ! First time step, read in data and set metdata1
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then  
     era5cds_struc(n)%findtime1=1
     era5cds_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
     era5cds_struc(n)%mon = mo1   ! current month
  endif

  ! Flag date/time when it is the switch of a month:
  if( era5cds_struc(n)%mon .ne. mo2 ) then
     era5cds_struc(n)%findtime1 = 1
     era5cds_struc(n)%mon = mo2
     ! advance filename to next hour
     yr1 = yr2
     mo1 = mo2
     da1 = da2
     hr1 = hr2
  endif

  ! Flag date/time for when to open and read next ERA5 file:
  if( LIS_get_nstep(LIS_rc,n) == 1 &
     .or. (era5cds_struc(n)%findtime1 == 1 ) ) then
     retrieve_file = .true.
     era5cds_struc(n)%findtime1=0
  endif

  if(timenow.ge.era5cds_struc(n)%era5cdstime2) then 
     movetime = 1
     era5cds_struc(n)%findtime2 = 1
  endif

  if(movetime.eq.1) then
     era5cds_struc(n)%era5cdstime1=era5cds_struc(n)%era5cdstime2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           era5cds_struc(n)%metdata1(:,f,c)=era5cds_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if( retrieve_file ) then
     !- Obtaining ERA5 File:
        do kk= era5cds_struc(n)%st_iterid, era5cds_struc(n)%en_iterid
           if (LIS_get_nstep(LIS_rc,n) == 1) then
            order = 1
           else
            order = 2
            era5cds_struc(n)%findtime2 = 0
           endif
           call era5cdsfiles(n,kk,findex,era5cds_struc(n)%era5cdsdir, &
                yr1, mo1, da1, hr1, instfilename, avgfilename, lmlfilename, &
                prevavgfilename)
           write(LIS_logunit,*) '[INFO] opening Bookend1 ',LIS_get_nstep(LIS_rc,n),hr1
           call read_era5cds(n, kk, order, yr1, mo1, da1, hr1, retrieve_file,&
                findex, instfilename, avgfilename, lmlfilename, &
                prevavgfilename, ferror)
        enddo

        if(ferror.ge.1) then !successfully retrieved forcing data
          era5cds_struc(n)%era5cdstime1=time1
          retrieve_file = .false.
          if(era5cds_struc(n)%findtime2==0) era5cds_struc(n)%era5cdstime2=time2
        else  !ferror still=0
           write(LIS_logunit,*)'[ERR] ERA5CDS data missing file 1'
           call LISrun()
        endif
  endif   

  if(era5cds_struc(n)%findtime2.eq.1) then
     ! just need to assign metdata2
        do kk= era5cds_struc(n)%st_iterid, era5cds_struc(n)%en_iterid
           order = 2
           call era5cdsfiles(n,kk,findex,era5cds_struc(n)%era5cdsdir, &
                yr2, mo2, da2, hr2, instfilename, avgfilename, lmlfilename, &
                prevavgfilename)
           write(LIS_logunit,*) '[INFO] using Bookend2 ',LIS_get_nstep(LIS_rc,n),hr2
           call read_era5cds(n, kk, order, yr2, mo2, da2, hr2, retrieve_file, &
               findex, instfilename, avgfilename, lmlfilename, &
               prevavgfilename, ferror)
        end do

        if(ferror.ge.1) then !successfully retrieved forcing data
           era5cds_struc(n)%era5cdstime2=time2
        else
           write(LIS_logunit,*)'[ERR] ERA5CDS data missing file 2'
           call LIS_endrun()
        endif
  endif 

end subroutine get_era5cds


!BOP
! !ROUTINE: era5cdsfiles
! \label{era5cdsfiles}
!
! !INTERFACE:
subroutine era5cdsfiles(n, kk, findex, era5cdsdir, yr, mo, da, hr, &
  instfilename, avgfilename, lmlfilename, prevavgfilename)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_forecastMod
  use LIS_timeMgrMod

  implicit none
! !ARGUMENTS:
  integer                       :: n 
  integer                       :: kk
  integer                       :: findex
  character(len=*), intent(in)  :: era5cdsdir
  integer, intent(in)           :: yr,mo,da,hr
  character(len=*), intent(out) :: instfilename, avgfilename, lmlfilename
  character(len=*), intent(out) :: prevavgfilename

! !DESCRIPTION:
!   This subroutine puts together ERA5 file names for
!   daily netcdf files
!
!  The arguments are:
!  \begin{description}
!  \item[era5cdsdir]
!    Name of the ERA5 directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[instfilename]
!   name of the timestamped ERA5 instantaneous file
!  \item[avgfilename]
!   name of the timestamped ERA5 accumulation file
!  \item[lmlfilename]
!   name of the timestamped ERA5 lowest model level file
!  \item[prevavgfilename]
!   name of the previous month timestamped ERA5 accumulation file
!  \end{description}
!
!EOP

  character*4  :: cyear
  character*2  :: cmonth
  character*2  :: cdate
  character*2  :: chour
  integer      :: mn, ss
  real*8       :: time
  integer      :: doy
  real         :: gmt
  integer      :: iyr,imo,ida,ihr,imn,iss,ts,idoy
  real         :: igmt
  real*8       :: itime
  character    :: ciyr*4, cimo*2, cida*2

  mn = 0 
  ss = 0 

  if(LIS_rc%forecastMode.eq.0) then !hindcast run
     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     write(unit=cdate,fmt='(i2.2)') da
     write(unit=chour,fmt='(i2.2)') hr
     
     instfilename = trim(era5cdsdir)//'/'//cyear//'/era5_forcing_'//cyear//cmonth//'_instant.nc'
     lmlfilename = trim(era5cdsdir)//'/'//cyear//'/era5_forcing_'//cyear//cmonth//'_ml_instant.nc'
     avgfilename = trim(era5cdsdir)//'/'//cyear//'/era5_forcing_'//cyear//cmonth//'_accum.nc'
     ! accum file starts from 7z on the 1st of current month and ends on 6z on 
     ! the 1st of next month, and need previous months
     !== roll back for last day of last month
     iyr=yr;  imo=mo;  ida=da
     ihr=hr;  imn=0;   iss=0
     ts = -24*60*60*da
     call LIS_tick(itime,idoy,igmt,iyr,imo,ida,ihr,imn,iss,real(ts))
     write(ciyr, '(i4.4)') iyr
     write(cimo, '(i2.2)') imo
     write(cida, '(i2.2)') ida
     if ( imo /= mo ) then
      prevavgfilename = trim(era5cdsdir)//'/'//ciyr//'/era5_forcing_'//ciyr//cimo//'_accum.nc'
     else
      prevavgfilename = "none"
     endif

  else !forecast mode
     !sample yr, mo, da

     call LIS_sample_forecastDate(n, kk, findex, yr,mo,da)

     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     write(unit=cdate,fmt='(i2.2)') da
     
     instfilename = trim(era5cdsdir)//'/'//cyear//'/era5_forcing_'//cyear//cmonth//'_instant.nc'
     avgfilename = trim(era5cdsdir)//'/'//cyear//'/era5_forcing_'//cyear//cmonth//'_accum.nc'
     lmlfilename = trim(era5cdsdir)//'/'//cyear//'/era5_forcing_'//cyear//cmonth//'_ml_instant.nc'
     ! accum file starts from 7z on the 1st of current month and ends on 6z on 
     ! the 1st of next month, and need previous months
     !== roll back for last day of last month
     iyr=yr;  imo=mo;  ida=da
     ihr=hr;  imn=0;   iss=0
     ts = -24*60*60*da
     call LIS_tick(itime,idoy,igmt,iyr,imo,ida,ihr,imn,iss,real(ts))
     write(ciyr, '(i4.4)') iyr
     write(cimo, '(i2.2)') imo
     write(cida, '(i2.2)') ida
     if ( imo /= mo ) then
      prevavgfilename = trim(era5cdsdir)//'/'//ciyr//'/era5_forcing_'//ciyr//cimo//'_accum.nc'
     else
      prevavgfilename = "none"
     endif
     
  endif
end subroutine era5cdsfiles

