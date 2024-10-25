!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: get_gddp
! \label{get_gddp}
!  
! !REVISION HISTORY:
!
! 03 Feb 2022; Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine get_gddp(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_metforcingMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_constantsMod
  use gddp_forcingMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates the daily GDDP forcing. 
!.
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the GDDP data times
!  \item[read\_gddp](\ref{read_gddp}) \newline
!      Interpolates GDDP data to LIS grid
!  \end{description}
  !EOP
  
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer            :: c,f,order,ferror,try
  integer            :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer            :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  integer            :: ryr1,rmo1,rda1,rhr1,rmn1,rss1,rdoy1
  integer            :: ryr2,rmo2,rda2,rhr2,rmn2,rss2,rdoy2
  real*8             :: time1,time2,rtime1,rtime2
  real*8             :: dtime1, dtime2
  real*8             :: timenow
  real               :: gmt1,gmt2,rgmt1,rgmt2,ts1,ts2,ts3,rts

  character(len=LIS_CONST_PATH_LEN) :: names(gddp_struc(n)%nmodels,7)
  character(len=LIS_CONST_PATH_LEN) :: ref_dailyclimofile
  character(len=LIS_CONST_PATH_LEN) :: ref_hourlyclimofile(7)
    
  integer :: movetime      ! 1=move time 2 data into time 1
  integer :: nforce     ! GDDP forcing file time, # forcing variables
  integer :: nstep

  nstep = LIS_get_nstep(LIS_rc,n)
!-------------------------------------------------------------------
! Determine the correct number of forcing variables
!-------------------------------------------------------------------
  nforce = LIS_rc%met_nf(findex)
  gddp_struc(n)%findtime1=0
  gddp_struc(n)%findtime2=0
  LIS_rc%shortflag = 2
  LIS_rc%longflag=2             !Time averaged LW 
  movetime=0
!-------------------------------------------------------------------
! Determine Required GDDP Data Times 
! (The previous hour & the future hour)
!-------------------------------------------------------------------
  yr1=LIS_rc%yr    !Time now
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0        
  
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  if(mod(nint(LIS_rc%ts),3600).eq.0) then

     ryr1=LIS_rc%yr
     rmo1=LIS_rc%mo
     rda1=LIS_rc%da
     rhr1=0 
     rmn1=0
     rss1=0
     rts=-24*60*60
     call LIS_tick(rtime1,rdoy1,rgmt1,ryr1,rmo1,rda1,rhr1,rmn1,rss1,rts)

     ryr2=LIS_rc%yr
     rmo2=LIS_rc%mo
     rda2=LIS_rc%da
     rhr2=0 
     rmn2=0
     rss2=0     
     rts=0
     call LIS_tick(rtime2,rdoy2,rgmt2,ryr2,rmo2,rda2,rhr2,rmn2,rss2,rts)
     
     if(timenow.ge.gddp_struc(n)%gddptime2) then
        yr1=LIS_rc%yr
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
        gddp_struc(n)%findtime2 = 1
     endif
  else

     if(timenow.ge.gddp_struc(n)%gddptime2) then
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
        gddp_struc(n)%findtime2 = 1
     endif
  endif
  
 if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then   
     gddp_struc(n)%findtime1=1
     gddp_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif

  if(movetime.eq.1) then
     gddp_struc(n)%gddptime1=gddp_struc(n)%gddptime2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           gddp_struc(n)%metdata1(:,f,c)=gddp_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1

  if(gddp_struc(n)%findtime1.eq.1) then

     ferror=0
     try=0
     ts1=-60*60*24
     do
        if ( ferror /= 0 ) exit
        try=try+1
        
        call create_gddpfilenames(n,       &
             gddp_struc(n)%odir,           &
             gddp_struc(n)%ref_dclimodir,  &
             gddp_struc(n)%ref_hclimodir,  &
             names,                        &
             ref_dailyclimofile,           &
             ref_hourlyclimofile,          & 
             ryr1,rmo1,rda1,hr1)
        
        order = 1
        call read_gddp(n,findex,order,ryr1,rdoy1, &
             names,                         &
             ref_dailyclimofile,            &
             ref_hourlyclimofile,           &
             ferror)

        if(ferror.ge.1) gddp_struc(n)%gddptime1=time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] GDDP data gap exceeds 10 days on file 1'
           write(LIS_logunit,*)'[ERR] Program stopping...'
           call LIS_endrun()
        endif
     enddo
!=== end of data search
  endif   !end of LIS_rc%findtime=1  


  if(gddp_struc(n)%findtime2.eq.1) then
     ferror=0
     try=0
     ts2=-60*60*24
     do
        if ( ferror /= 0 ) exit
        try=try+1
        
        call create_gddpfilenames(n,       &
             gddp_struc(n)%odir,           &
             gddp_struc(n)%ref_dclimodir,  &
             gddp_struc(n)%ref_hclimodir,  &
             names,                        &
             ref_dailyclimofile,           &
             ref_hourlyclimofile,          & 
             yr2,mo2,da2,hr2)
       
        order = 2
        call read_gddp(n,findex,order,yr2,rdoy2,&
             names,                         &
             ref_dailyclimofile,            &
             ref_hourlyclimofile,           &
             ferror)
        
        if(ferror.ge.1) then
           gddp_struc(n)%gddptime2=time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] GDDP data gap exceeds 10 days on file 2'
           write(LIS_logunit,*)'[ERR] Program stopping...'
           call LIS_endrun()
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1
  
end subroutine get_gddp

!BOP
!
! !ROUTINE: create_gddpfilename
! \label{create_gddpfilename}
!  
!
! !INTERFACE:
subroutine create_gddpfilenames(n,&
     gddpdir,                &
     refdclimodir,           &
     refhclimodir,           & 
     fnames,                 &
     ref_dailyclimofile,     &
     ref_hourlyclimofile,    & 
     yr, mo, da, hr)
! !USES:
  use gddp_forcingMod,    only : gddp_struc
  use LIS_constantsMod
!
! !DESCRIPTION:
! This routine creates the filenames for the GDDP data. 25 ensemble members
! are employed in the current setup
!
!EOP
  implicit none

! !ARGUMENTS:
  integer                         :: n

  character(len=*), intent(out)   :: fnames(gddp_struc(n)%nmodels,7)
  character(len=*), intent(in)    :: gddpdir
  character(len=*), intent(in)    :: refdclimodir
  character(len=*), intent(in)    :: refhclimodir
  character(len=*)                :: ref_dailyclimofile
  character(len=*)                :: ref_hourlyclimofile(7)
  integer, intent(in)             :: yr,mo,da,hr

  integer                         :: doy
  integer                         :: t,k
  
  character*4                     :: fyr
  character*2                     :: fmo
  character*2                     :: fda
  character*2                     :: fhr
  character*3                     :: fdoy

  character(len=LIS_CONST_PATH_LEN) :: gddpdir2
  character(len=LIS_CONST_PATH_LEN) :: dname(gddp_struc(n)%nmodels)
  character(len=LIS_CONST_PATH_LEN) :: mname(gddp_struc(n)%nmodels)
  character(len=LIS_CONST_PATH_LEN) :: extn1(gddp_struc(n)%nmodels)
  character(len=LIS_CONST_PATH_LEN) :: extn2(gddp_struc(n)%nmodels)

  
  dname(1) = 'ACCESS-CM2'
  dname(2) = 'ACCESS-ESM1-5'
  dname(3) = 'CESM2'
  dname(4) = 'CESM2-WACCM'
  dname(5) = 'CMCC-CM2-SR5'
  dname(6) = 'CMCC-ESM2'
  dname(7) = 'CNRM-CM6-1'
  dname(8) = 'CNRM-ESM2-1'
  dname(9) = 'EC-Earth3'
  dname(10) = 'FGOALS-g3'
  dname(11) = 'GFDL-CM4'
  dname(12) = 'GFDL-CM4_gr2'
  dname(13) = 'GFDL-ESM4'
  dname(14) = 'GISS-E2-1-G'
  dname(15) = 'IITM-ESM'
  dname(16) = 'INM-CM4-8'
  dname(17) = 'INM-CM5-0'
  dname(18) = 'KACE-1-0-G'
  dname(19) = 'MIROC-ES2L'
  dname(20) = 'MPI-ESM1-2-HR'
  dname(21) = 'MPI-ESM1-2-LR'
  dname(22) = 'MRI-ESM2-0'
  dname(23) = 'NorESM2-LM'
  dname(24) = 'NorESM2-MM'
  dname(25) = 'TaiESM1'

  mname = dname
  mname(12) = 'GFDL-CM4'
 
  extn2(1) =  'gn'
  extn2(2) =  'gn'
  extn2(3) =  'gn'
  extn2(4) =  'gn'
  extn2(5) =  'gn'
  extn2(6) =  'gn'
  extn2(7) =  'gr'
  extn2(8) =  'gr'
  extn2(9) =  'gr'
  extn2(10) = 'gn'
  extn2(11) = 'gr1'
  extn2(12) = 'gr2'
  extn2(13) = 'gr1'
  extn2(14) = 'gn'
  extn2(15) = 'gn'
  extn2(16) = 'gr1'
  extn2(17) = 'gr1'
  extn2(18) = 'gr'
  extn2(19) = 'gn'
  extn2(20) = 'gn'
  extn2(21) = 'gn'
  extn2(22) = 'gn'
  extn2(23) = 'gn'
  extn2(24) = 'gn'
  extn2(25) = 'gn'

  extn1(1) = 'r1i1p1f1'
  extn1(2) = 'r1i1p1f1'
  extn1(3) = 'r4i1p1f1'
  extn1(4) = 'r3i1p1f1'
  extn1(5) = 'r1i1p1f1'
  extn1(6) = 'r1i1p1f1'
  extn1(7) = 'r1i1p1f2'
  extn1(8) = 'r1i1p1f2'
  extn1(9) = 'r1i1p1f1'
  extn1(10) = 'r3i1p1f1'
  extn1(11) = 'r1i1p1f1'
  extn1(12) = 'r1i1p1f1'
  extn1(13) = 'r1i1p1f1'  
  extn1(14) = 'r1i1p1f2'
  extn1(15) = 'r1i1p1f1'
  extn1(16) = 'r1i1p1f1'
  extn1(17) = 'r1i1p1f1'
  extn1(18) = 'r1i1p1f1'
  extn1(19) = 'r1i1p1f2'
  extn1(20) = 'r1i1p1f1'  
  extn1(21) = 'r1i1p1f1'
  extn1(22) = 'r1i1p1f1'
  extn1(23) = 'r1i1p1f1'
  extn1(24) = 'r1i1p1f1'
  extn1(25) = 'r1i1p1f1'


  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  !because the files use a 1 to 24 convention
  write(unit=fhr, fmt='(i2.2)') hr+1

  !
  ! compute DOY consistent with the climatology calculation where
  ! the number of days in a year are from 1 to 366.
  !

  call compute_doy(yr,mo,da,doy)
  write(unit=fdoy, fmt='(i3.3)') doy


  ref_dailyclimofile =trim(refdclimodir)//&
       '/MERRA2_DAILY_CLIMO'//trim(fdoy)//'.nc' 

  ref_hourlyclimofile(1) = trim(refhclimodir)//&
       '/TAIR/'//'MERRA2_HOURLY_CLIMO_TAIR_'//trim(fdoy)//&
       '_'//trim(fhr)//'.nc'
  ref_hourlyclimofile(2) = trim(refhclimodir)//&
       '/QAIR/'//'MERRA2_HOURLY_CLIMO_QAIR_'//trim(fdoy)//&
       '_'//trim(fhr)//'.nc'
  ref_hourlyclimofile(3) = trim(refhclimodir)//&
       '/SWDOWN/'//'MERRA2_HOURLY_CLIMO_SWDOWN_'//trim(fdoy)//&
       '_'//trim(fhr)//'.nc'
  ref_hourlyclimofile(4) = trim(refhclimodir)//&
       '/LWDOWN/'//'MERRA2_HOURLY_CLIMO_LWDOWN_'//trim(fdoy)//&
       '_'//trim(fhr)//'.nc'
  ref_hourlyclimofile(5) = trim(refhclimodir)//&
       '/WIND/'//'MERRA2_HOURLY_CLIMO_WIND_'//trim(fdoy)//&
       '_'//trim(fhr)//'.nc'
  ref_hourlyclimofile(6) = trim(refhclimodir)//&
       '/PSURF/'//'MERRA2_HOURLY_CLIMO_PSURF_'//trim(fdoy)//&
       '_'//trim(fhr)//'.nc'
  ref_hourlyclimofile(7) = trim(refhclimodir)//&
       '/PRCP/'//'MERRA2_HOURLY_CLIMO_PRCP_'//trim(fdoy)//&
       '_'//trim(fhr)//'.nc' 
  
  do t=1,gddp_struc(n)%nmodels
     
     gddpdir2 = trim(gddpdir)//'/'//trim(dname(t))//'/'//&
          trim(gddp_struc(n)%scenario)//'/'//&
          trim(extn1(t))//'/'
     
     do k=1,7        
        fnames(t,1) = trim(gddpdir2)//'/tas/tas_day_'//trim(mname(t))//&
             '_'//trim(gddp_struc(n)%scenario)//'_'//&
             trim(extn1(t))//'_'//trim(extn2(t))//'_'//&
             trim(fyr)//'.nc'
        fnames(t,2) = trim(gddpdir2)//'/huss/huss_day_'//trim(mname(t))//&
             '_'//trim(gddp_struc(n)%scenario)//'_'//&
             trim(extn1(t))//'_'//trim(extn2(t))//'_'//&
             trim(fyr)//'.nc'
        fnames(t,3) = trim(gddpdir2)//'/rsds/rsds_day_'//trim(mname(t))//&
             '_'//trim(gddp_struc(n)%scenario)//'_'//&
             trim(extn1(t))//'_'//trim(extn2(t))//'_'//&
             trim(fyr)//'.nc'
        fnames(t,4) = trim(gddpdir2)//'/rlds/rlds_day_'//trim(mname(t))//&
             '_'//trim(gddp_struc(n)%scenario)//'_'//&
             trim(extn1(t))//'_'//trim(extn2(t))//'_'//&
             trim(fyr)//'.nc'
        fnames(t,5) = trim(gddpdir2)//'/sfcWind/sfcWind_day_'//trim(mname(t))//&
             '_'//trim(gddp_struc(n)%scenario)//'_'//&
             trim(extn1(t))//'_'//trim(extn2(t))//'_'//&
          trim(fyr)//'.nc'
        fnames(t,6) = trim(gddpdir2)//'/pr/pr_day_'//trim(mname(t))//&
             '_'//trim(gddp_struc(n)%scenario)//'_'//&
             trim(extn1(t))//'_'//trim(extn2(t))//'_'//&
             trim(fyr)//'.nc'
        fnames(t,7) = trim(gddpdir2)//'/hurs/hurs_day_'//trim(mname(t))//&
             '_'//trim(gddp_struc(n)%scenario)//'_'//&
             trim(extn1(t))//'_'//trim(extn2(t))//'_'//&
             trim(fyr)//'.nc'

     enddo
  enddo
  
end subroutine create_gddpfilenames

subroutine compute_doy(yr,mo,da,doy)

  implicit none
! !ARGUMENTS:
  integer :: yr,mo,da,doy

! !DESCRIPTION:
!
!  determines the time, time in GMT, and the day of the year
!  based on the value of year, month, day of month, hour of
!  the day, minute and second. This method is the inverse of
!  time2date
!
!  The arguments are:
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[time]
!   lvt time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
  integer :: yrdays,days(13),k
  data days /31,28,31,30,31,30,31,31,30,31,30,31,30/

  if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
       .or.(mod(yr,400).eq.0))then             !correct for y2k
     yrdays=366
     days(2) = 29
  else
     yrdays=365
     days(2) = 29
  endif
  yrdays = 366

  doy=0
  do k=1,(mo-1)
     doy=doy+days(k)
  enddo
  doy=doy+da
  
end subroutine compute_doy
