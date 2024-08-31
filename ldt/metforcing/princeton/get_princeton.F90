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
! !ROUTINE:  get_princeton
!  \label{get_princeton}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Kato; Initial Specification adopted from LDT
! 
! !INTERFACE:
subroutine get_princeton(n,findex)
! !USES:
  use LDT_coreMod,          only : LDT_rc 
  use LDT_metforcingMod,    only : LDT_forc
  use LDT_timeMgrMod,       only : LDT_get_nstep, LDT_tick
  use LDT_logMod,           only : LDT_logunit
  use princeton_forcingMod, only : princeton_struc

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, 1 degree 
!  Princeton forcing. At the beginning of a simulation, the code 
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
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    call to advance or retract time
!  \item[read\_princeton](\ref{read_princeton}) \newline
!    call to read the Princeton data and perform spatial interpolation 
!   \item[read\_princeton\_elev](\ref{read_princeton_elev}) \newline
!    reads the native elevation of the Princeton data to be used
!    for topographic adjustments to the forcing
!  \end{description}
!EOP
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer, parameter :: nforce = 9  ! # forcing variables
  integer :: try, ferror
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8  :: time1,time2,dumbtime1,dumbtime2
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime      ! 1=move time 2 data into time 1

  princeton_struc(n)%findtime1=0
  princeton_struc(n)%findtime2=0
  movetime=0

  !=== Determine Required Data Times (The previous hour & the future hour)
  yr1=LDT_rc%yr    !Previous Hour
  mo1=LDT_rc%mo
  da1=LDT_rc%da
  hr1=3*((LDT_rc%hr)/3)
  mn1=0
  ss1=0
  ts1=0
  call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=LDT_rc%yr    !Next Hour
  mo2=LDT_rc%mo
  da2=LDT_rc%da
  hr2=3*((LDT_rc%hr)/3)
  mn2=0
  ss2=0
  ts2=3*60*60
  call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2) 

  !=== Check if time interval boundary was crossed
  if(LDT_rc%time.gt.princeton_struc(n)%princetontime2) then
     movetime=1
     princeton_struc(n)%findtime2=1
  endif

  if (LDT_get_nstep(LDT_rc,n) .eq. 1.or.LDT_rc%rstflag(n).eq.1) then ! beginning of the run
     princeton_struc(n)%findtime1=1
     princeton_struc(n)%findtime2=1
     movetime=0
     LDT_rc%rstflag(n) = 0
  endif
    
  !=== Establish fmodeltime1
  if (princeton_struc(n)%findtime1==1) then  !need to get new time1 from the past
     order=1   !Get data for glbdata1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if ( ferror /= 0 ) then
           exit
        end if
        try = try+1
        call read_princeton(order,n, findex, yr1,mo1,da1,hr1,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           princeton_struc(n)%princetontime1=time1
        else  !ferror still=0, so roll back one day
           call LDT_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        end if
        if ( try > ndays ) then 
           write(LDT_logunit,*) 'ERROR: Princeton data gap exceeds 10 days'
           STOP
        end if
     end do
  endif

  !=== Find new time2 value and tranfer time2 data to time1 
  if(movetime.eq.1) then
     princeton_struc(n)%princetontime1=princeton_struc(n)%princetontime2
     princeton_struc(n)%findtime2=1 !to ensure getting new time2 data
     do f=1,nforce
        do c=1,LDT_rc%ngrid(n)
           LDT_forc(n,findex)%metdata1(f,c)=LDT_forc(n,findex)%metdata2(f,c)
        enddo
     enddo
  endif  ! if movetime=1
  
  if(princeton_struc(n)%findtime2.eq.1) then ! need new time2 data
     order=2   !Get data for glbdata2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call read_princeton(order,n,findex,yr2,mo2,da2,hr2,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           write(LDT_logunit,*) 'reset princetontime2 to time2'
           princeton_struc(n)%princetontime2=time2
        else  !ferror still=0, so roll back one day
           call LDT_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        end if
        if ( try > ndays ) then
           write(LDT_logunit,*)'ERROR: Princeton data gap exceeds 10 days'
           STOP
        end if
     end do
  endif
  
end subroutine get_princeton
   






