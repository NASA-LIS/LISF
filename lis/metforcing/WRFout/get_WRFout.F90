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
! !ROUTINE: get_WRFout
! \label{get_WRFout}
!
!
! !REVISION HISTORY:
! 14 Mar  2013 : Sujay Kumar ; Initial Version in LIS
!
! !INTERFACE:
subroutine get_WRFout(n, findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_timeMgrMod,    only : LIS_tick
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
  use WRFout_forcingMod, only : WRFout_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens and reads 1-hourly WRF output forcing. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the WRF output data times
!  \item[WRFoutfile](\ref{WRFoutfile}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_WRFout](\ref{read_WRFout}) \newline
!      Reads WRF output data to LIS grid
!  \end{description}
!EOP

  integer      :: ferror,try
  real*8       :: time1,time2,timenow
  integer      :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer      :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real         :: ts1, ts2
  character(len=LIS_CONST_PATH_LEN) :: fname
  real         :: gmt1,gmt2
  integer      :: movetime     ! 1=move time 2 data into time 1

!=== End Variable Definition =============================================
  try=-999
  
!====Assumption will be not to find or move any data
  WRFout_struc(n)%findtime1=0
  WRFout_struc(n)%findtime2=0
  movetime=0

!=== Determine Required NCEP Data Times (The previous hour and the future hour)
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

  if(mod(LIS_rc%ts,3600.0).eq.0) then 
     if(timenow.ge.WRFout_struc(n)%WRFouttime2) then 
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
        WRFout_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.WRFout_struc(n)%WRFouttime2) then 
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
        WRFout_struc(n)%findtime2 = 1
     endif
  endif
    
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     WRFout_struc(n)%findtime1=1
     WRFout_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  if(movetime.eq.1) then
     WRFout_struc(n)%WRFouttime1=WRFout_struc(n)%WRFouttime2
     WRFout_struc(n)%metdata1 = WRFout_struc(n)%metdata2
  endif    !end of movetime=1
  
! ---- File 1 ----

  if(WRFout_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1

        call WRFoutfile(fname,WRFout_struc(n)%WRFoutdir,&
             WRFout_struc(n)%nest_id, yr1,mo1,da1,hr1,mn1,ss1)

        write(unit=LIS_logunit,fmt=*)'[INFO] getting file1.. ',trim(fname)
        call read_WRFout(n,findex,1,fname,ferror)

        if(ferror.ge.1) WRFout_struc(n)%WRFouttime1=time1
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: WRFout data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
!=== end of data search
  endif   !end of LIS_rc%findtime=1	   	

! ---- File 2 ----

  if(WRFout_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        
        call WRFoutfile(fname,WRFout_struc(n)%WRFoutdir,&
             WRFout_struc(n)%nest_id, yr2,mo2,da2,hr2,mn2,ss2)

        write(unit=LIS_logunit,fmt=*)'[INFO] getting file2.. ',trim(fname)
        call read_WRFout(n,findex,2,fname,ferror)

        if(ferror.ge.1) then
           WRFout_struc(n)%WRFouttime2=time2
        endif
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: WRFout data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1

end subroutine get_WRFout

!BOP
! !ROUTINE: WRFoutfile
! \label{WRFoutfile}
!
! !REVISION HISTORY:
! 14 Mar 2013: Sujay Kumar, initial specification
!
! !INTERFACE:
 subroutine WRFoutfile(filename,wrfdir,nest,yr,mo,da,hr,mn,ss)

   implicit none
! !ARGUMENTS: 
   character(len=*), intent(out) :: filename
   character(len=*), intent(in)  :: wrfdir
   integer, intent(in)       :: nest
   integer, intent(in)       :: yr,mo,da,hr,mn,ss

! !DESCRIPTION:
!
!EOP

   character*10     :: ftime1
   character*2      :: fnest
   character*4      :: fyr
   character*2      :: fmo
   character*2      :: fda
   character*2      :: fhr
   character*2      :: fmn
   character*2      :: fss

   write(unit=ftime1, fmt='(i4.4,i2.2,i2.2,i2.2)') yr,mo,da,hr
   write(unit=fnest,fmt='(i2.2)') nest
   write(unit=fyr,fmt='(i4.4)') yr
   write(unit=fmo,fmt='(i2.2)') mo
   write(unit=fda,fmt='(i2.2)') da
   write(unit=fhr,fmt='(i2.2)') hr
   write(unit=fmn,fmt='(i2.2)') mn
   write(unit=fss,fmt='(i2.2)') ss

   if(nest.ne.1) then 
      filename = trim(wrfdir)//'/'//& 
           'wrfout_d'//trim(fnest)//'_'//trim(fyr)//'-'//trim(fmo)//'-'//&
           trim(fda)//'_'//trim(fhr)//':'//trim(fmn)//':'//trim(fss)
   else

! EMK...No special rules. Users can set up symbolic links 
! if their files aren't written at exact times due to time step constraints.

      filename = trim(wrfdir)//'/'//& 
           'wrfout_d'//trim(fnest)//'_'//trim(fyr)//'-'//trim(fmo)//'-'//&
           trim(fda)//'_'//trim(fhr)//':'//trim(fmn)//':'//trim(fss)

!       if(hr.eq.0.or.hr.eq.3.or.hr.eq.6.or.hr.eq.9.or.hr.eq.12.or.hr.eq.15.or.&
!            hr.eq.18.or.hr.eq.21.or.hr.eq.24) then 
!          filename = trim(wrfdir)//'/'//& 
!               'wrfout_d'//trim(fnest)//'_'//trim(fyr)//'-'//trim(fmo)//'-'//&
!               trim(fda)//'_'//trim(fhr)//':'//trim(fmn)//':'//trim(fss)
!       elseif(hr.eq.1.or.hr.eq.4.or.hr.eq.7.or.hr.eq.10.or.hr.eq.13.or.&
!            hr.eq.16.or.hr.eq.19.or.hr.eq.22) then 
!          filename = trim(wrfdir)//'/'//& 
!               'wrfout_d'//trim(fnest)//'_'//trim(fyr)//'-'//trim(fmo)//'-'//&
!               trim(fda)//'_'//trim(fhr)//':'//trim(fmn)//':18'
!       elseif(hr.eq.2.or.hr.eq.5.or.hr.eq.8.or.hr.eq.11.or.hr.eq.14.or.&
!            hr.eq.17.or.hr.eq.20.or.hr.eq.23) then 
!          filename = trim(wrfdir)//'/'//& 
!               'wrfout_d'//trim(fnest)//'_'//trim(fyr)//'-'//trim(fmo)//'-'//&
!               trim(fda)//'_'//trim(fhr)//':'//trim(fmn)//':09'
!       endif
         
   endif
      
 end subroutine WRFoutfile

