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
! !ROUTINE:  get_WRFoutv2
!  \label{get_WRFoutv2}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Kato; Initial Specification adopted from LIS
!  20 Nov 2020: K.R. Arsenault; Updated for different WRF output files
! 
! !INTERFACE:
subroutine get_WRFoutv2(n,findex)

! !USES:
  use LIS_coreMod,          only : LIS_rc 
  use LIS_timeMgrMod,       only : LIS_get_nstep, LIS_tick
  use LIS_logMod,           only : LIS_logunit, LIS_endrun
  use WRFoutv2_forcingMod,  only : WRFoutv2_struc

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly, 4-km Lambert
!  WRFout forcing. At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 1-hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to go backwards up to 10 days 
!  to get forcing at the same time of day.
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
!    call to advance or retract time
!  \item[read\_WRFoutv2](\ref{read_WRFoutv2}) \newline
!    call to read the WRFout data and perform spatial interpolation 
!   \item[read\_WRFoutv2\_elev](\ref{read_WRFoutv2_elev}) \newline
!    reads the native elevation of the WRFout data to be used
!    for topographic adjustments to the forcing
!  \end{description}
!
!EOP
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer :: try, ferror
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8  :: time1,time2,dumbtime1,dumbtime2
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime      ! 1=move time 2 data into time 1
  !===

  WRFoutv2_struc(n)%findtime1=0
  WRFoutv2_struc(n)%findtime2=0
  movetime=0

  !=== Determine Required Data Times (The previous hour & the future hour)

  !=== Check if time interval boundary was crossed
  ! Input files are hourly
 if(mod(nint(LIS_rc%ts),3600).eq.0) then
  if( LIS_rc%time.ge.WRFoutv2_struc(n)%WRFouttime2 ) then
    yr1=LIS_rc%yr    ! Previous Hour
    mo1=LIS_rc%mo
    da1=LIS_rc%da
    hr1=LIS_rc%hr   
    mn1=0
    ss1=0
    ts1=-60*60
    call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

    yr2=LIS_rc%yr    ! Next Hour
    mo2=LIS_rc%mo
    da2=LIS_rc%da
    hr2=LIS_rc%hr    
    mn2=0
    ss2=0
    ts2=0
    call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
    movetime=1
    WRFoutv2_struc(n)%findtime2=1
  endif

 else
  if( LIS_rc%time.ge.WRFoutv2_struc(n)%WRFouttime2 ) then
    yr1=LIS_rc%yr    ! Previous Hour
    mo1=LIS_rc%mo
    da1=LIS_rc%da
    hr1=LIS_rc%hr  
    mn1=0
    ss1=0
    ts1=0
    call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

    yr2=LIS_rc%yr    ! Next Hour
    mo2=LIS_rc%mo
    da2=LIS_rc%da
    hr2=LIS_rc%hr    
    mn2=0
    ss2=0
    ts2=60*60
    call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
    movetime=1
    WRFoutv2_struc(n)%findtime2=1
  endif
 endif

  !=== Beginning of the run
  if( LIS_rc%tscount(n).eq.1 .or. LIS_rc%rstflag(n).eq.1 ) then
     WRFoutv2_struc(n)%findtime1=1
     WRFoutv2_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif

  !=== Find new time2 value and tranfer time2 data to time1 
  if(movetime.eq.1) then
     WRFoutv2_struc(n)%WRFouttime1 = WRFoutv2_struc(n)%WRFouttime2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           WRFoutv2_struc(n)%metdata1(:,f,c)=WRFoutv2_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif  ! if movetime=1
    
  !=== Look or "roll" back 10 days to see if previous days exist

  if( WRFoutv2_struc(n)%findtime1==1 ) then  ! Need to get new time1 from the past
     order = 1   ! Get data for glbdata1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if( ferror /= 0 ) then
          exit
        end if
        try = try+1
        ! Obtain file and variable records:
        call read_WRFoutv2(order,n,findex,yr1,mo1,da1,hr1,ferror)

        if ( ferror == 1 ) then !successfully retrieved forcing data
           WRFoutv2_struc(n)%WRFouttime1=time1
        else  ! ferror still=0, so roll back one day
           call LIS_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        end if
        if ( try > ndays ) then 
           write(LIS_logunit,*) '[ERR] WRFout data file 1 gap exceeds 10 days'
           write(LIS_logunit,*) '[ERR] LIS run stopping ...'
           call LIS_endrun
        end if
     end do
  endif

  if( WRFoutv2_struc(n)%findtime2.eq.1 ) then ! need new time2 data
     order = 2     ! Get data for glbdata2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if ( ferror /= 0 ) then
          exit
        endif
        try = try+1
        call read_WRFoutv2(order,n,findex,yr2,mo2,da2,hr2,ferror)

        if ( ferror == 1 ) then !successfully retrieved forcing data
           write(LIS_logunit,*) '[INFO] reset WRFouttime2 to time2'
           WRFoutv2_struc(n)%WRFouttime2=time2
        else  ! ferror still=0, so roll back one day
           call LIS_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        end if
        if ( try > ndays ) then
           write(LIS_logunit,*)'[ERR] WRFout data file 2 gap exceeds 10 days'
           write(LIS_logunit,*)'[ERR] Stopping LIS run.'
           call LIS_endrun
        end if
     end do
  endif
  
end subroutine get_WRFoutv2
   
