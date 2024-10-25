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
! !ROUTINE: get_scan
! \label{get_scan}
!
! !REVISION HISTORY:
! 13 Apr 2007: Bailing Li; Initial Specification
!
! !INTERFACE:
subroutine get_scan(n,findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_logMod, only : LIS_logunit
  use LIS_metforcingMod,  only : LIS_forc
  use scan_forcingMod, only : scan_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates the SCAN station data. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data, and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the SCAN data times
!  \item[read\_scan](\ref{read_scan}) \newline
!      Interpolates the appropriate SCAN station data to LIS grid
!  \end{description}
!EOP

  integer :: order
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real :: gmt1,gmt2,ts1,ts2
  real*8 :: timenow,time1,time2
  integer :: movetime      ! if 1=move time 2 data into time 1
  integer :: t

  movetime = 0 

  yr1 = LIS_rc%yr  !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   

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

  scan_struc(n)%findtime1 = 0
  scan_struc(n)%findtime2 = 0
  if(timenow .ge. scan_struc(n)%starttime .and. &
       .not.scan_struc(n)%startRead) then 
     scan_struc(n)%findtime1 = 1
     scan_struc(n)%findtime2 = 1
     scan_struc(n)%startRead = .true.
     movetime = 0
  endif
  if(scan_struc(n)%startRead) then 
     if(timenow.ge.scan_struc(n)%scantime2) then
        movetime = 1
        scan_struc(n)%findtime2 = 1
     endif

!Time to open file and start reading..
!keep on reading until the obstime is reached.     
     if(scan_struc(n)%findtime1.eq.1) then 
        write(LIS_logunit,*) 'reading time1 data...'
        order = 1
        call read_scan(n,222,findex,order)
        scan_struc(n)%scantime1 = time1
     endif
     if(movetime .eq. 1) then 
        scan_struc(n)%scantime1 = scan_struc(n)%scantime2
        do t=1,LIS_rc%ngrid(n)
           scan_struc(n)%metdata1(1,t) = scan_struc(n)%metdata2(1,t)
        enddo
     endif

     if(scan_struc(n)%findtime2.eq.1) then 
        write(LIS_logunit,*) 'reading time2 data...'
        order = 2
        call read_scan(n,222,findex,order)
        scan_struc(n)%scantime2 = time2
     endif

  endif
end subroutine get_scan


