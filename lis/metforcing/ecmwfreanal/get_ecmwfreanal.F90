!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_ecmwfreanal
!  \label{get_ecmwfreanal}
!
! !REVISION HISTORY:
!  11 Apr 2002: Urszula Jambor; original code based on getgeos.f
!  22 Oct 2002: Urszula Jambor; Limited SW forcing processing to 
!               land-only grid points
!  24 Nov 2003: Sujay Kumar; Included the scheme in LIS
! !INTERFACE:
subroutine get_ecmwfreanal(n,findex)
! !USES:
  use LIS_coreMod,            only : LIS_rc
  use LIS_logMod,             only : LIS_logunit
  use LIS_timeMgrMod,         only : LIS_get_nstep, LIS_tick
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 6-hrly, 1/2 degree Reanalysis 
!  ECMWF forcing. At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 6 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
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
!  \item[read\_ecmwfreanal](\ref{read_ecmwfreanal}) \newline
!    call to read the ECMWF Reanalysis data and perform spatial interpolation 
!  \end{description}
!EOP
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer, parameter :: nforce=9  ! # forcing variables
  integer :: try, ferror
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8  :: time1,time2,dumbtime1,dumbtime2
  real*8  :: timenow
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime      ! 1=move time 2 data into time 1
 
  !=== Assumption will be not to find or move any data
  ecmwfreanal_struc(n)%findtime1=0
  ecmwfreanal_struc(n)%findtime2=0
  movetime=0
  !=== Determine Required Data Times (The previous hour & the future hour)

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
  hr1=6*((LIS_rc%hr)/6)
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=LIS_rc%yr    !Next Hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=6*((LIS_rc%hr)/6)
  mn2=0
  ss2=0
  ts2=6*60*60

  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2) 

  if(timenow.gt.ecmwfreanal_struc(n)%fmodeltime2) then
     movetime=1
     ecmwfreanal_struc(n)%findtime2=1
  endif

  if ( LIS_get_nstep(LIS_rc,n) == 1 .or. &
       LIS_rc%rstflag(n) == 1 ) then !beginning of the run
     ecmwfreanal_struc(n)%findtime1=1
     ecmwfreanal_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
      
  !=== Establish fmodeltime1
  if (ecmwfreanal_struc(n)%findtime1==1) then  !need to get new time1 from the past
     order=1   !Get data for glbdata1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if ( ferror /= 0 ) then
           exit
        end if
        try = try+1
        call read_ecmwfreanal(order,n, findex, yr1,mo1,da1,hr1,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           ecmwfreanal_struc(n)%fmodeltime1=time1
        else  !ferror still=0, so roll back one day
           call LIS_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        end if
        if ( try > ndays ) then 
           write(LIS_logunit,*) 'ERROR: ECMWF Reanalysis data gap exceeds 10 days on file 1'
           STOP
        end if
     end do
  endif
  
  !  Repeat for time 2
  
  if(movetime.eq.1) then !transfer time2 data to time1
     ecmwfreanal_struc(n)%fmodeltime1=ecmwfreanal_struc(n)%fmodeltime2
     ecmwfreanal_struc(n)%findtime2=1 !include to ensure getting new time2 data
     
     do f=1,nforce
        do c=1,LIS_rc%ngrid(n)
           ecmwfreanal_struc(n)%metdata1(f,c)=ecmwfreanal_struc(n)%metdata2(f,c)
        enddo
     enddo
     
  endif  ! if movetime=1
  
  if(ecmwfreanal_struc(n)%findtime2.eq.1) then ! need new time2 data
     order=2   !Get data for glbdata2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call read_ecmwfreanal(order,n, findex, yr2,mo2,da2,hr2,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           ecmwfreanal_struc(n)%fmodeltime2=time2
        else  !ferror still=0, so roll back one day
           call LIS_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        end if
        if ( try > ndays ) then 
           write(LIS_logunit,*) 'ERROR: ECMWF Reanalysis data gap exceeds 10 days on file 2'
           STOP
        end if
     end do
  endif

end subroutine get_ecmwfreanal
   






