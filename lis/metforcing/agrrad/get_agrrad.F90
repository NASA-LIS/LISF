!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_agrrad
! \label{get_agrrad}
!
! !REVISION HISTORY:
! 29 Jul 2005; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine get_agrrad(n, findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use LIS_timeMgrMod,    only : LIS_tick
  use agrrad_forcingMod, only : agrrad_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  This is the entry point for calling various AGRMET radiation analyses. 
!  At the beginning of a simulation, the code invokes call to 
!  read the most recent past data (nearest hourly interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[read\_agrrad](\ref{read_agrrad}) \newline
!    routines that reads the AGRMET radiation and interpolates
!    it to the LIS projection. 
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    computes the AGRMET read/processing times. 
!  \end{description}
!EOP
  integer                  :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer                  :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real                     :: ts1,ts2
  real                     :: gmt1,gmt2
  real*8                   :: time1,time2,timenow
  integer                  :: movetime

  agrrad_struc(n)%findtime1 = 0
  agrrad_struc(n)%findtime2 = 0 
  movetime = 0 

  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo 
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0 
  ts1 = 0 
  call LIS_tick(timenow, doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = 3*(int(real(LIS_rc%hr)/3.0))
  mn1 = 0 
  ss1 = 0
  ts1 = 0 
  call LIS_tick(time1, doy1, gmt1, yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr2 = LIS_rc%yr
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 3*(int(real(LIS_rc%hr)/3.0))
  mn2 = 0 
  ss2 = 0
  ts2 = 3*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
  
  if(timenow.ge.agrrad_struc(n)%agrtime2) then 
     movetime = 1
     agrrad_struc(n)%findtime2 = 1
  endif
  
  if(LIS_rc%tscount(n).eq.0 .or. LIS_rc%tscount(n).eq.1 .or. &
     LIS_rc%rstflag(n).eq.1) then 
     agrrad_struc(n)%findtime1 = 1
     agrrad_struc(n)%findtime2 = 1
     movetime = 0 
     LIS_rc%rstflag(n) = 0 
  endif
  if(movetime.eq.1) then 
     agrrad_struc(n)%agrtime1 = agrrad_struc(n)%agrtime2
     agrrad_struc(n)%metdata1 = agrrad_struc(n)%metdata2
  endif
  if(agrrad_struc(n)%findtime1.eq.1) then 
     call read_agrrad(n,findex,1,yr1,mo1,da1,hr1)
     agrrad_struc(n)%agrtime1 = time1
  endif

  if(agrrad_struc(n)%findtime2.eq.1) then 
     call read_agrrad(n,findex,2,yr2,mo2,da2,hr2)
     agrrad_struc(n)%agrtime2 = time2
  endif
end subroutine get_agrrad
