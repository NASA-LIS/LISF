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
! !ROUTINE: get_agrradps
! \label{get_agrradps}
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
! 16 Feb 2007: Chuck Alonge; Changed file name creation to use internal 
!               files instead of external file (avoids failure in parallel runs)
! !INTERFACE:
subroutine get_agrradps(n, findex)
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use agrradps_forcingMod, only : agrradps_struc
  use LIS_timeMgrMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 6-hrly, CMAP forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 6 hour interval), and
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
!    determines the data times
!  \item[read\_agrradps](\ref{read_agrradps}) \newline
!    routines that reads the AGRMET radiation (polar stereographic) and interpolates
!    it to the LIS projection. 
!  \end{description}
!EOP
   
!==== Local Variables=======================
  integer                  :: c,f
  integer                  :: order
  integer                  :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer                  :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real                     :: ts1, ts2
  real                     :: gmt1,gmt2
  real*8                   :: time1,time2,timenow
  integer                  :: movetime
  real :: gridDesci(50)
!=== End Variable Definition =======================

  agrradps_struc(n)%findtime1 = 0
  agrradps_struc(n)%findtime2 = 0 
  movetime = 0 
!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end time)
! Model current time
!------------------------------------------------------------------------
  yr1 = LIS_rc%yr  !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   
!------------------------------------------------------------------------ 
! AGRRAD product start and end time
!------------------------------------------------------------------------
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call LIS_tick(time1, doy1, gmt1, yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2 = LIS_rc%yr
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = LIS_rc%hr
  mn2 = 0
  ss2 = 0
  ts2 = 60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  if(timenow.ge.agrradps_struc(n)%agrtime2) then
     movetime = 1
     agrradps_struc(n)%findtime2 = 1
  endif
  if(LIS_rc%tscount(n).eq.1.or.LIS_rc%rstflag(n).eq.1) then
     agrradps_struc(n)%findtime1 = 1
     agrradps_struc(n)%findtime2 = 1
     movetime = 0
     LIS_rc%rstflag(n) = 0
  endif
  if(movetime.eq.1) then
     agrradps_struc(n)%agrtime1 = agrradps_struc(n)%agrtime2
     do f=1,2
        do c=1,LIS_rc%ngrid(n)
          agrradps_struc(n)%metdata1(f,c) = agrradps_struc(n)%metdata2(f,c)
        enddo
     enddo
  endif
  if(agrradps_struc(n)%findtime1.eq.1) then
     order = 1
     call read_agrradps(n,findex,order,yr1,mo1,da1,hr1)
     agrradps_struc(n)%agrtime1 = time1
  endif

  if(agrradps_struc(n)%findtime2.eq.1) then
     order =2
     call read_agrradps(n,findex,order,yr2,mo2,da2,hr2)
     agrradps_struc(n)%agrtime2 = time2
  endif
 return
end subroutine get_agrradps

