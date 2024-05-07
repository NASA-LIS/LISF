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
! !ROUTINE: get_gswp2
! \label{get_gswp2}
!  
! !REVISION HISTORY:
!
! 20Feb2004; Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine get_gswp2(n, findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_metforcingMod, only : LIS_forc
  use LIS_timeMgrMod,     only : LIS_get_nstep, LIS_tick
  use LIS_logMod,         only : LIS_logunit
  use gswp2_forcingMod,    only : gswp2_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, GSWP2 forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
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
!    determines the GSWP2 data times
!  \item[read\_gswp2](\ref{read_gswp2}) \newline
!      Interpolates GSWP2 data to LIS grid
!  \end{description}
!EOP
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  integer :: yr3,mo3,da3,hr3,mn3,ss3,doy3
  real*8 :: time1,time2,time3
  real*8 :: timenow
  real :: gmt1,gmt2,gmt3,ts1,ts2,ts3
  integer :: movetime      ! 1=move time 2 data into time 1
  integer :: nforce     ! GSWP2 forcing file time, # forcing variables
  integer :: nstep

  nstep = LIS_get_nstep(LIS_rc,n)
!-------------------------------------------------------------------
! Determine the correct number of forcing variables
!-------------------------------------------------------------------
  nforce = LIS_rc%met_nf(findex)
  gswp2_struc(n)%findtime1=0
  gswp2_struc(n)%findtime2=0
  LIS_rc%shortflag = 2
  LIS_rc%longflag=2             !Time averaged LW 
  movetime=0
!-------------------------------------------------------------------
! Determine Required GSWP2 Data Times 
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
      
  yr1=LIS_rc%yr    !Previous Hour
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=3*((LIS_rc%hr)/3)
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr2=LIS_rc%yr    !Next Hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=3*((LIS_rc%hr)/3)
  mn2=0
  ss2=0
  ts2=3*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  yr3=LIS_rc%yr    !Uber-Next Hour
  mo3=LIS_rc%mo
  da3=LIS_rc%da
  hr3=3*((LIS_rc%hr)/3)
  mn3=0
  ss3=0
  ts3=2*3*60*60
  call LIS_tick(time3,doy3,gmt3,yr3,mo3,da3,hr3,mn3,ss3,ts3)

  ! For GSWP2, we must read past (previous hour), current (next hour), and
  ! future (uber-next hour) forcing time-steps.
  !
  ! So, for GSWP2, LIS_rc%findtime1 means update past and current forcing data.
  !
  ! LIS_rc%findtime2 means update the future forcing data.
  ! 
  ! gswp2_struc(n)%gswp2time1 is the update time for the past forcing time-step,
  ! which is only needed at initialization.
  !
  ! gswp2_struc(n)%gswp2time2 is the update time for the current forcing time-step.
  !
  ! When it is time to update the current forcing time-step, we actually shift
  ! the old current to past, the old future to current, and we read in
  ! the new future.
  if ( timenow .gt.gswp2_struc(n)%gswp2time2 ) then
     movetime        = 1
     gswp2_struc(n)%findtime2 = 1
  endif
  
  if ( nstep == 0 .or. nstep == 1 .or. LIS_rc%rstflag(n) == 1 ) then 
     gswp2_struc(n)%findtime1 = 1
     gswp2_struc(n)%findtime2 = 1
     movetime        = 0
     LIS_rc%rstflag(n)   = 0
  endif

  if ( gswp2_struc(n)%findtime1 == 1 ) then 
     write(LIS_logunit,*)'MSG: get_gswp2 -- reading time1 data'
     write(LIS_logunit,*)'get_gswp2',yr1,mo1,da1,hr1,mn1,ss1
     order = 1
     call read_gswp2(order,n, findex, yr1,mo1,da1,hr1,mn1,ss1)
     gswp2_struc(n)%gswp2time1 = time1

     if(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
        write(LIS_logunit,*)'MSG: get_gswp2 -- reading time2 data'
        write(LIS_logunit,*)'get_gswp2',yr2,mo2,da2,hr2,mn2,ss2
        order = 2   
        call read_gswp2(order,n, findex, yr2,mo2,da2,hr2,mn2,ss2)
        gswp2_struc(n)%gswp2time2 = time2
     endif
  endif

  if ( movetime == 1 ) then
     gswp2_struc(n)%gswp2time1=gswp2_struc(n)%gswp2time2
     gswp2_struc(n)%findtime2=1 
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           gswp2_struc(n)%metdata1(f,c)=gswp2_struc(n)%metdata2(f,c)
           if(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              gswp2_struc(n)%metdata2(f,c)=gswp2_struc(n)%metdata3(f,c)
           endif
        enddo
     enddo
  endif 

  if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
     if ( gswp2_struc(n)%findtime2 == 1 ) then 
        write(LIS_logunit,*)'MSG: get_gswp2 -- reading time2 data'
        write(LIS_logunit,*)'get_gswp2',yr2,mo2,da2,hr2,mn2,ss2
        order = 2 
        call read_gswp2(order,n, findex, yr2,mo2,da2,hr2,mn2,ss2)
        gswp2_struc(n)%gswp2time2 = time2
     endif
  elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
     if ( gswp2_struc(n)%findtime2 == 1 ) then 
        write(LIS_logunit,*)'MSG: get_gswp2 -- reading time3 data'
        write(LIS_logunit,*)'get_gswp2',yr3,mo3,da3,hr3,mn3,ss3
        order = 3   
        call read_gswp2(order,n,findex, yr3,mo3,da3,hr3,mn3,ss3)
        gswp2_struc(n)%gswp2time2 = time2
     endif
  endif

end subroutine get_gswp2

