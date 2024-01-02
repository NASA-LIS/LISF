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
! !ROUTINE: get_narr
!  \label{get_narr}
!
! !REVISION HISTORY:
!  30 APR 2009; Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine get_narr(n, findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_metforcingMod, only : LIS_forc
  use LIS_timeMgrMod,     only : LIS_get_nstep, LIS_tick
  use LIS_logMod,         only : LIS_logunit
  use narr_forcingMod,  only : narr_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!EOP
  
  integer             :: yr1,mo1,da1,hr1,mn1,ss1
  integer             :: yr2,mo2,da2,hr2,mn2,ss2
  real*8              :: time1, time2, timenow
  integer             :: movetime
  integer             :: doy1, doy2
  real                :: gmt1, gmt2,ts1,ts2
  integer             :: order
  integer             :: nstep
  character(len=LIS_CONST_PATH_LEN) :: narrfile

  narr_struc(n)%findtime1 = 0 
  narr_struc(n)%findtime2 = 0 
  movetime = 0 

  nstep=LIS_get_nstep(LIS_rc,n)
  
  if(nstep.eq.0.or.nstep.eq.1.or.LIS_rc%rstflag(n).eq.1) then
     narr_struc(n)%findtime1=1
     narr_struc(n)%findtime2=1
     movetime=0        ! movetime is not properly set at time-step = 1
     LIS_rc%rstflag(n) = 0
  endif

!-----------------------------------------------------------------
! Determine required NARR data times 
! (previous assimilation, current & future assimilation hours) 
!-----------------------------------------------------------------
! Current time
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0 
  ts1 = 0 
  call LIS_tick(timenow, doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

!previous hour (nearest 3 hour interval)
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = 3*(int(real(LIS_rc%hr)/6.0))
  mn1 = 0 
  ss1 = 0 
  ts1 = 0 

  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
!Next hour (nearest 6 hour interval)
  yr2 = LIS_rc%yr
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 3*(int(real(LIS_rc%hr)/3.0))
  mn2 = 0 
  ss2 = 0 
  ts2 = 3*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  if(nstep.eq.0) then 
     narr_struc(n)%narrtime1 = time1
     narr_struc(n)%narrtime2 = time2
  endif

  if(timenow > narr_struc(n)%narrtime2) then 
     movetime = 1
     narr_struc(n)%findtime2 = 1
  endif

!-----------------------------------------------------------------
! Reading first bookmark
!-----------------------------------------------------------------

  if(narr_struc(n)%findtime1.eq.1) then 
    
     call narr_filename(narr_struc(n)%narrdir,narrfile,&
          yr1,mo1,da1,hr1)

     write(LIS_logunit,*) 'NARR file1 ',trim(narrfile)
     order = 1
     call read_narr(n, findex, order, narrfile)
     
     narr_struc(n)%narrtime1 = time1
  endif

  if(movetime.eq.1) then 
     narr_struc(n)%narrtime1 = narr_struc(n)%narrtime2
     narr_struc(n)%findtime2 = 1
     
! need to transfer data...

  endif
  
  if(narr_struc(n)%findtime2.eq.1) then 

     call narr_filename(narr_struc(n)%narrdir,narrfile,&
          yr2,mo2,da2,hr2)

     write(LIS_logunit,*) 'NARR file2 ',trim(narrfile)
     order = 2
     call read_narr(n, findex, order, narrfile)
     narr_struc(n)%narrtime2 = time2
  endif

end subroutine get_narr

 
subroutine narr_filename(odir,filename,yr,mo,da,hr)
  
  implicit none

  character (len=*)    :: odir
  character (len=*)    :: filename
  integer              :: yr
  integer              :: mo
  integer              :: da
  integer              :: hr

  character(len=8)     :: ftime1
  character(len=4)     :: fyr
  character(len=2)     :: fmo
  character(len=2)     :: fhr

  write(unit=ftime1,fmt='(i4,i2.2,i2.2)') yr, mo, da
  write(unit=fyr,fmt='(i4.4)') yr
  write(unit=fmo,fmt='(i2.2)') mo
  write(unit=fhr,fmt='(i2.2)') hr
  
  filename = trim(odir)//'/'//trim(fyr)//'/'//trim(fyr)//trim(fmo)//'/narr-a_221_'//trim(ftime1)//'_'//&
       trim(fhr)//'00_000.grb'

  
end subroutine narr_filename


