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
! !ROUTINE: get_hmaens
! \label{get_hmaens}
!
!
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code
!
! !INTERFACE:
subroutine get_hmaens(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use hmaens_forcingMod

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly HMAENS forcing.
!
!  The HMAENS forcing data are organized into monthly files, where each
!  file contains 24 one-hourly records of forcing fields.  
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
!  \item[hmaensfiles](\ref{hmaensfiles}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_hmaens](\ref{read_hmaens}) \newline
!    call to read the HMAENS data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: order
  integer           :: ferror
  character*100     :: fname
  integer           :: c, r,kk,f,try
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2, timenow
  real*8            :: dtime1, dtime2
  real              :: gmt1, gmt2
  real              :: ts1, ts2

  integer           :: hr_int1, hr_int2
  integer           :: movetime  ! Flag to move bookend2 files to bookend1

! _________________________________________________________

  if( LIS_rc%nts(n).gt.3600 ) then   ! > 1-hr timestep
     write(LIS_logunit,*) '[ERR] When running LIS with HMAENS, the clock '
     write(LIS_logunit,*) '[ERR] should run with a timestep less than or '
     write(LIS_logunit,*) '[ERR] equal to one hour.'
     call LIS_endrun()
  endif

  hmaens_struc(n)%findtime1 = 0
  hmaens_struc(n)%findtime2 = 0
  movetime = 0

!=== Determine Required HMAENS Data Times (The previous hour and the future hour)
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 
  if(LIS_rc%ts.gt.3600) then 
     write(LIS_logunit,*) '[ERR] The model timestep is > forcing data timestep'
     write(LIS_logunit,*) '[ERR] LIS does not support this mode currently'
     write(LIS_logunit,*) '[ERR] Program stopping ...'
     call LIS_endrun()
  endif

  if(mod(nint(LIS_rc%ts),3600).eq.0) then 
     if(timenow.ge.hmaens_struc(n)%hmaenstime2) then 
        yr1 = LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=0
        ss1=0
        ts1=-60*60
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        
        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=0
        ss2=0
        ts2=0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        hmaens_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.hmaens_struc(n)%hmaenstime2) then 
        yr1 = LIS_rc%yr
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

        movetime = 1
        hmaens_struc(n)%findtime2 = 1
     endif
  endif
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then  
     hmaens_struc(n)%findtime1=1
     hmaens_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  
  if(movetime.eq.1) then
     hmaens_struc(n)%hmaenstime1=hmaens_struc(n)%hmaenstime2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           hmaens_struc(n)%metdata1(:,f,c)=hmaens_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(hmaens_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        do kk= hmaens_struc(n)%st_iterid, hmaens_struc(n)%en_iterid
           order = 1
           call hmaensfiles(n,kk,findex,hmaens_struc(n)%hmaensdir, yr1, mo1, da1, &
                fname)
           call read_hmaens(n, kk,order, yr1,mo1, da1, hr1, &
                findex, fname, ferror)
        enddo

        if(ferror.ge.1) hmaens_struc(n)%hmaenstime1=time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] HMAENS data gap exceeds 10 days on file 1'
           call LISrun()
        endif
     enddo
!=== end of data search
  endif   

  if(hmaens_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1

     !- Obtaining HMAENS File:
        do kk= hmaens_struc(n)%st_iterid, hmaens_struc(n)%en_iterid
          order = 2
          call hmaensfiles(n,kk,findex,hmaens_struc(n)%hmaensdir, yr2, mo2, da2, &
               fname)
          call read_hmaens(n, kk,order, yr2,mo2, da2, hr2, &
               findex, fname, ferror)
        end do

        if(ferror.ge.1) then
           hmaens_struc(n)%hmaenstime2=time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] HMAENS data gap exceeds 10 days on file 2'
           call LIS_endrun()
        endif
     enddo
  endif 

end subroutine get_hmaens


!BOP
! !ROUTINE: hmaensfiles
! \label{hmaensfiles}
!
! !INTERFACE:
subroutine hmaensfiles(n, kk, findex, hmaensdir, yr, mo, da, hr, fname)


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
  character(len=*), intent(in)  :: hmaensdir
  integer, intent(in)           :: yr,mo,da
  character(len=*), intent(out) :: fname

! !DESCRIPTION:
!   This subroutine puts together HMAENS file names for
!   daily netcdf files
!
!  The arguments are:
!  \begin{description}
!  \item[hmaensdir]
!    Name of the HMAENS directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[fname]
!   name of the timestamped HMAENS file
!  \end{description}
!
!EOP

  character*4  :: cyear
  character*2  :: cmonth
  character*8  :: cdate
  character*4  :: cday
  character*4  :: chour
  character*20 :: dir
  integer      :: hr, mn, ss
  real*8       :: time
  integer      :: doy
  real         :: gmt

  hr = 0 
  mn = 0 
  ss = 0 

  if(LIS_rc%forecastMode.eq.0) then !hindcast run
     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     write(unit=cday,fmt='(i2.2)') da
     write(unit=chour,fmt='(i2.2)') hr

     fname = trim(hmaensdir)//'/ERA5_LPM_5km_'//trim(cyear)//trim(cmonth)//trim(cday)//trim(chour)//'.nc'
    

  else !forecast mode
     !sample yr, mo, da

     call LIS_sample_forecastDate(n, kk, findex, yr,mo,da)

     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     write(unit=cday,fmt='(i2.2)') da
     write(unit=chour,fmt='(i2.2)') hr
     
     fname = trim(hmaensdir)//'/ERA5_LPM_5km_'//trim(cyear)//trim(cmonth)//trim(cday)//trim(chour)//'.nc'     
     
  endif
end subroutine hmaensfiles
