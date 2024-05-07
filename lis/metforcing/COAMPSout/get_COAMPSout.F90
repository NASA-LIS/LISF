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
! !ROUTINE: get_COAMPSout
! \label{get_COAMPSout}
!
!
! !REVISION HISTORY:
! 14 Mar  2013 : Sujay Kumar ; Initial Version in LIS
!
! !INTERFACE:
subroutine get_COAMPSout(n, findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_timeMgrMod,    only : LIS_tick
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
  use COAMPSout_forcingMod, only : COAMPSout_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens and reads 1-hourly COAMPS output forcing. 
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
!    determines the COAMPS output data times
!  \item[COAMPSoutfile](\ref{COAMPSoutfile}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_COAMPSout](\ref{read_COAMPSout}) \newline
!      Reads COAMPS output data to LIS grid
!  \end{description}
!EOP

  integer      :: ferror,try
  real*8       :: time1,time2,timenow
  integer      :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer      :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real         :: ts1, ts2
  integer      :: fcsthr
  character(len=LIS_CONST_PATH_LEN) :: fname
  real         :: gmt1,gmt2
  integer      :: movetime     ! 1=move time 2 data into time 1

!=== End Variable Definition =============================================
  try=-999
  
!====Assumption will be not to find or move any data
  COAMPSout_struc(n)%findtime1=0
  COAMPSout_struc(n)%findtime2=0
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
     if(timenow.ge.COAMPSout_struc(n)%COAMPSouttime2) then 
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
        COAMPSout_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.COAMPSout_struc(n)%COAMPSouttime2) then 
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
        COAMPSout_struc(n)%findtime2 = 1
     endif
  endif
    
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     COAMPSout_struc(n)%findtime1=1
     COAMPSout_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  if(movetime.eq.1) then
     COAMPSout_struc(n)%COAMPSouttime1=COAMPSout_struc(n)%COAMPSouttime2
     COAMPSout_struc(n)%metdata1 = COAMPSout_struc(n)%metdata2
  endif    !end of movetime=1
  
! ---- File 1 ----

  if(COAMPSout_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        
        fcsthr = (hr1/6)*6
        call COAMPSoutfile(fname,COAMPSout_struc(n)%COAMPSoutdir,&
             COAMPSout_struc(n)%nest_id, fcsthr, yr1,mo1,da1,hr1,mn1,ss1)

        write(unit=LIS_logunit,fmt=*)'[INFO] getting file1.. ',trim(fname)
        call read_COAMPSout(n,findex,1,fname,ferror)

        if(ferror.ge.1) COAMPSout_struc(n)%COAMPSouttime1=time1
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] COAMPSout data gap exceeds 10 days on file 1'
           write(LIS_logunit,*)'[ERR] Program stopping'
           call LIS_endrun()
        endif
     enddo
!=== end of data search
  endif   !end of LIS_rc%findtime=1	   	

! ---- File 2 ----

  if(COAMPSout_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1

        fcsthr = (hr2/6)*6
        call COAMPSoutfile(fname,COAMPSout_struc(n)%COAMPSoutdir,&
             COAMPSout_struc(n)%nest_id, fcsthr, yr2,mo2,da2,hr2,mn2,ss2)

        write(unit=LIS_logunit,fmt=*)'[INFO] getting file2.. ',trim(fname)
        call read_COAMPSout(n,findex,2,fname,ferror)

        if(ferror.ge.1) then
           COAMPSout_struc(n)%COAMPSouttime2=time2
        endif
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] COAMPSout data gap exceeds 10 days on file 2'
           write(LIS_logunit,*)'[ERR] Program stopping...'
           call LIS_endrun()
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1

end subroutine get_COAMPSout

!BOP
! !ROUTINE: COAMPSoutfile
! \label{COAMPSoutfile}
!
! !REVISION HISTORY:
! 14 Mar 2013: Sujay Kumar, initial specification
!
! !INTERFACE:
 subroutine COAMPSoutfile(filename,coampsdir,nest,fcsthr,&
      yr,mo,da,hr,mn,ss)

   implicit none
! !ARGUMENTS: 
   character(len=*), intent(out) :: filename
   character(len=*), intent(in)  :: coampsdir
   integer, intent(in)       :: nest
   integer, intent(in)       :: fcsthr
   integer, intent(in)       :: yr,mo,da,hr,mn,ss

! !DESCRIPTION:
!
!EOP

   integer          :: hr1
   character*10     :: ftime1
   character*1      :: fnest
   character*4      :: fyr
   character*2      :: fmo
   character*2      :: fda
   character*2      :: fhr,fhr1
   character*2      :: fmn
   character*2      :: fss

   hr1 = hr - fcsthr
   write(unit=ftime1, fmt='(i4.4,i2.2,i2.2,i2.2)') yr,mo,da,hr
   write(unit=fnest,fmt='(i1.1)') nest
   write(unit=fyr,fmt='(i4.4)') yr
   write(unit=fmo,fmt='(i2.2)') mo
   write(unit=fda,fmt='(i2.2)') da
   write(unit=fhr,fmt='(i2.2)') hr1
   write(unit=fhr1,fmt='(i2.2)') fcsthr
   write(unit=fmn,fmt='(i2.2)') mn
   write(unit=fss,fmt='(i2.2)') ss
   
   filename = trim(coampsdir)//'/'//& 
        trim(fyr)//trim(fmo)//trim(fda)//trim(fhr1)//&
        '/coamps_'//trim(fnest)//'_'//trim(fyr)//trim(fmo)//&
        trim(fda)//trim(fhr1)//'_00'//trim(fhr)//'0000.nc'
      
 end subroutine COAMPSoutfile
