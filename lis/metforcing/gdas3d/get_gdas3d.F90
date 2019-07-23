!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_gdas3d
!  \label{get_gdas3d}
!
! !REVISION HISTORY:

!
! !INTERFACE:
subroutine get_gdas3d(n, findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_timeMgrMod,     only : LIS_get_nstep, LIS_tick
  use LIS_logMod,         only : LIS_logunit
  use gdas3d_forcingMod,  only : gdas3d_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!EOP
  
  integer             :: order
  integer             :: yr1,mo1,da1,hr1,mn1,ss1
  integer             :: yr2,mo2,da2,hr2,mn2,ss2
  real*8              :: time1, time2, timenow
  integer             :: movetime
  integer             :: doy1, doy2
  real                :: gmt1, gmt2,ts1,ts2
  integer             :: nstep
  character*100       :: sanlfile

  gdas3d_struc(n)%findtime1 = 0 
  gdas3d_struc(n)%findtime2 = 0 
  movetime = 0 

  nstep=LIS_get_nstep(LIS_rc,n)
  
  if(nstep.eq.0.or.nstep.eq.1.or.LIS_rc%rstflag(n).eq.1) then
     gdas3d_struc(n)%findtime1=1
     gdas3d_struc(n)%findtime2=1
     movetime=0        ! movetime is not properly set at time-step = 1
     LIS_rc%rstflag(n) = 0
  endif

!-----------------------------------------------------------------
! Determine required GDAS data times 
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

!previous hour (nearest 6 hour interval)
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = 6*(int(real(LIS_rc%hr)/6.0))
  mn1 = 0 
  ss1 = 0 
  ts1 = 0 

  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
!Next hour (nearest 6 hour interval)
  yr2 = LIS_rc%yr
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 6*(int(real(LIS_rc%hr)/6.0))
  mn2 = 0 
  ss2 = 0 
  ts2 = 6*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  if(nstep.eq.0) then 
     gdas3d_struc(n)%gdastime1 = time1
     gdas3d_struc(n)%gdastime2 = time2
  endif

  if(timenow > gdas3d_struc(n)%gdastime2) then 
     movetime = 1
     gdas3d_struc(n)%findtime2 = 1
  endif

!-----------------------------------------------------------------
! Reading first bookmark
!-----------------------------------------------------------------

  if(gdas3d_struc(n)%findtime1.eq.1) then 
    
     call gdas3d_sanlfile(gdas3d_struc(n)%gdasdir,sanlfile,&
          yr1,mo1,da1,hr1)

     write(LIS_logunit,*) 'sanl file1 ',trim(sanlfile)
     order = 1
     call read_gdas3d(n, findex, order, sanlfile)
     
     gdas3d_struc(n)%gdastime1 = time1
  endif

  if(movetime.eq.1) then 
     gdas3d_struc(n)%gdastime1 = gdas3d_struc(n)%gdastime2
     gdas3d_struc(n)%findtime2 = 1
     
! need to transfer data...

  endif
  
  if(gdas3d_struc(n)%findtime2.eq.1) then 

     call gdas3d_sanlfile(gdas3d_struc(n)%gdasdir,sanlfile,&
          yr2,mo2,da2,hr2)

     write(LIS_logunit,*) 'sanl file2 ',trim(sanlfile)
     order = 2
     call read_gdas3d(n, findex, order, sanlfile)
     gdas3d_struc(n)%gdastime2 = time2
  endif

end subroutine get_gdas3d

 
subroutine gdas3d_sanlfile(odir,filename,yr,mo,da,hr)
  
  implicit none

  character (len=*)    :: odir
  character (len=*)    :: filename
  integer              :: yr
  integer              :: mo
  integer              :: da
  integer              :: hr

  character(len=8)     :: ftime1
  character(len=2)     :: fyr

  write(unit=ftime1,fmt='(i4,i2.2,i2.2)') yr, mo, da
  write(unit=fyr,fmt='(i2.2)') hr
  
  filename = trim(odir)//trim(ftime1)//'/gdas.t'//trim(fyr)//'z.'//&
       trim(ftime1)//'.sanl.bin.be'

  
end subroutine gdas3d_sanlfile


