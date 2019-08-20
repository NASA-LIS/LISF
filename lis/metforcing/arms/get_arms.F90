!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_arms
! \label{get_arms}
!
! !REVISION HISTORY:
! 08 Dec 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine get_arms(n, findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_logMod, only : LIS_logunit
  use arms_forcingMod, only : arms_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates the ARMS station data. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data, and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
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
!    determines the ARMS data times
!  \item[read\_arms](\ref{read_arms}) \newline
!      Interpolates the appropriate ARMS station data to LIS grid
!  \end{description}
!EOP

  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real :: gmt1,gmt2,ts1,ts2
  real*8 :: timenow,time1,time2
  integer :: movetime      ! if 1=move time 2 data into time 1
  integer :: f,t
  integer :: order

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

  arms_struc(n)%findtime1=0
  arms_struc(n)%findtime2=0
  movetime = 0

  if(timenow .ge. arms_struc(n)%armstime2) then 
     arms_struc(n)%findtime2=1
     movetime = 1
  endif

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     arms_struc(n)%findtime1=1
     arms_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif

  if(movetime.eq.1) then 
     arms_struc(n)%armstime1 = arms_struc(n)%armstime2
     do f=1,8
        do t=1,LIS_rc%ngrid(n)
           arms_struc(n)%metdata1(f,t) = arms_struc(n)%metdata2(f,t)
        enddo
     enddo
  endif

  if(arms_struc(n)%findtime1.eq.1) then
     write(LIS_logunit,*) 'reading ARMS time1 data...'
     arms_struc(n)%armstime1 = time1
     order = 1
     call read_arms(n,arms_struc(n)%armstime1, findex,order)
  endif

  if(arms_struc(n)%findtime2.eq.1) then
     write(LIS_logunit,*) 'reading ARMS time2 data...'
     arms_struc(n)%armstime2 = time2
     order = 2
     call read_arms(n,arms_struc(n)%armstime2,findex,order)
  endif

end subroutine get_arms


