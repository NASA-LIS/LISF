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
! !ROUTINE: get_gswp1
! \label{get_gswp1}
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine get_gswp1(n, findex)
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use LIS_timeMgrMod,      only : LIS_get_nstep, LIS_tick
  use LIS_logMod,          only : LIS_logunit, LIS_endrun
  use gswp1_forcingMod,    only : gswp1_struc
  use LIS_constantsMod,    only : LIS_CONST_PATH_LEN
  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n 
  integer, intent(in)  :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly, GSWP1 forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest hourly interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day
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
!    determines the GSWP data times
!  \item[gswp1file](\ref{gswp1file}) \newline
!    determines the GSWP1 data file
!  \item[read\_gswp1](\ref{read_gswp1}) \newline
!     reads GSWP1 data to LIS grid
!  \end{description}
!EOP
  integer, parameter :: ndays = 10 ! # days to look back for forcing data
  integer :: try, ferror
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8 :: time1,time2,dumbtime1,dumbtime2
  real*8 :: timenow
  character(len=LIS_CONST_PATH_LEN) :: name
  real :: gmt1,gmt2,ts1,ts2
  integer :: movetime       ! 1=move time 2 data into time 1
  integer :: nforce         ! GSWP-1 forcing file time, # forcing variables
  integer :: nstep

  nstep = LIS_get_nstep(LIS_rc,n)

!-----------------------------------------------------------------------
! Determine the correct number of forcing variables
!-----------------------------------------------------------------------
  nforce = LIS_rc%met_nf(findex)

  gswp1_struc(n)%findtime1=0
  gswp1_struc(n)%findtime2=0
  movetime=0
!-----------------------------------------------------------------------
! Determine Required GSWP-1 Data Times
! (The previous hour & the future hour)
!-----------------------------------------------------------------------
  yr1=LIS_rc%yr              !Time now
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr1=LIS_rc%yr              !Previous Hour
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=1*((LIS_rc%hr)/1)
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr2=LIS_rc%yr              !Next Hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=1*((LIS_rc%hr)/1)
  mn2=0
  ss2=0
  ts2=1*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
  
  if (timenow.ge.gswp1_struc(n)%gswp1time2) then
     movetime=1
     gswp1_struc(n)%findtime2=1
  endif
  
  if ((nstep.eq.0).or.(nstep.eq.1).or.(LIS_rc%rstflag(n).eq.1)) then
     gswp1_struc(n)%findtime1=1
     gswp1_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag = 0
  endif
  LIS_rc%shortflag=2         !Time averaged SW
  LIS_rc%longflag=2          !Time averaged LW

!-----------------------------------------------------------------------
! Establish gswp1time1
!-----------------------------------------------------------------------
  if (gswp1_struc(n)%findtime1.eq.1) then
     order = 1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if (ferror.ne.0) then
           exit
        endif
        try = try+1
        call gswp1file(name,gswp1_struc(n)%gswp1dir,yr1,mo1,da1,hr1,     &
             gswp1_struc(n)%ncold)
        call read_gswp1(order,n,findex,name,LIS_rc%tscount(n),ferror)
        if (ferror.eq.1) then
!-----------------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------------
           gswp1_struc(n)%gswp1time1=time1
        else
!-----------------------------------------------------------------------
! ferror still=0, so roll back one day
!-----------------------------------------------------------------------
           call LIS_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        endif
        if (try.gt.ndays) then
           print *,'ERR: GSWP-1 data gap exceeds 10 days on file 1'
           call LIS_endrun
        endif
     enddo
  endif
  if (movetime.eq.1) then
     gswp1_struc(n)%gswp1time1=gswp1_struc(n)%gswp1time2
     gswp1_struc(n)%findtime2=1
     do f=1,nforce
        do c=1,LIS_rc%ngrid(n)
           gswp1_struc(n)%metdata1(f,c)=gswp1_struc(n)%metdata2(f,c)
        enddo
     enddo
  endif
  
!-----------------------------------------------------------------------
! Establish gswp1time2
!-----------------------------------------------------------------------
  if (gswp1_struc(n)%findtime2.eq.1) then
     order = 2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if (ferror.ne.0) exit
        try = try+1
        call gswp1file(name,gswp1_struc(n)%gswp1dir,yr2,mo2,da2,hr2,     &
             gswp1_struc(n)%ncold)
        call read_gswp1(order,n,findex,name,LIS_rc%tscount(n),ferror)
        if (ferror.eq.1) then
!-----------------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------------
           gswp1_struc(n)%gswp1time2=time2
        else
!-----------------------------------------------------------------------
! ferror still=0, so roll back one day
!-----------------------------------------------------------------------
           call LIS_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        endif
        if (try.gt.ndays) then
           print *,'ERR: GSWP-1 data gap exceeds 10 days on file 2'
           call LIS_endrun
        endif
     enddo
  endif
  
84 format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
  
  return
end subroutine get_gswp1

!BOP
! !ROUTINE: gswp1file
! \label{gswp1file}
!
! !INTERFACE:
subroutine gswp1file(name,gswp1dir,yr,mo,da,hr,ncold)

  implicit none
      
! !INPUT PARAMETERS:
  character(len=*), intent(in) :: gswp1dir
  integer, intent(in)      :: yr,mo,da,hr,ncold
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: name
!
! !DESCRIPTION:
!  This subroutine puts together GSWP-1 file name
!
!EOP
  integer uyr,umo,uda,uhr,i,c,ii,jj
  character(len=2) :: initcode
  character(len=6) :: fdir, fsubs
  character(len=10) :: ftime
  character*2 hrstr(24)
  data hrstr /'01','02','03','04','05','06','07','08','09','10',   &
       '11','12','13','14','15','16','17','18','19','20',   &
       '21','22','23','00'/
  
  ii = ncold
  jj = 150
!-----------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!-----------------------------------------------------------------------
  uyr=yr
  umo=mo
  uda=da
  uhr = 1*(hr/1)  !hour needs to be a multiple of 1 hour
  if (uhr.eq.0) uhr = 24
!-----------------------------------------------------------------------
!  Determine initcode for the hour of the forecast file
!  If the time is 12 or later the file is time stamped
!  with the next day.  So check for that first
!-----------------------------------------------------------------------
  initcode = hrstr(uhr)

  write(UNIT=fdir, FMT='(i4, i2.2)') uyr, umo

  write(UNIT=ftime, FMT='(i4, i2.2, i2.2, a2)') uyr, umo, uda, initcode

  fsubs = '.GSWP1'

  name = trim(gswp1dir) // '/' // fdir // '/' // ftime // fsubs
  
  return

end subroutine gswp1file

