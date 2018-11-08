!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_ceop
! \label{get_ceop}
!
! !REVISION HISTORY:
! 08 Dec 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine get_ceop(n, findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_logMod, only : LIS_logunit
  use ceop_forcingMod, only : ceop_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates the CEOP station data. 
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
!    determines the CEOP data times
!  \item[read\_ceop](\ref{read_ceop}) \newline
!      Interpolates the appropriate CEOP station data to LIS grid
!  \end{description}
!EOP

  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real :: gmt1,gmt2,ts1,ts2
  integer :: order
  real*8 :: timenow,time1,time2
  integer :: movetime      ! if 1=move time 2 data into time 1
  integer :: f,t

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

  ceop_struc(n)%findtime1 = 0
  ceop_struc(n)%findtime2 = 0
  if(timenow .ge. ceop_struc(n)%starttime .and. &
       .not.ceop_struc(n)%startRead) then 
     ceop_struc(n)%findtime1 = 1
     ceop_struc(n)%findtime2 = 1
     ceop_struc(n)%startRead = .true.
     movetime = 0
  endif
  if(ceop_struc(n)%startRead) then 
     if(timenow.ge.ceop_struc(n)%ceoptime2) then
        movetime = 1
        ceop_struc(n)%findtime2 = 1
     endif

!Time to open file and start reading..
!keep on reading until the obstime is reached.     
     if(ceop_struc(n)%findtime1.eq.1) then 
        write(LIS_logunit,*) 'reading time1 data...'
        order = 1
        call read_ceop(n,findex,order)
        ceop_struc(n)%ceoptime1 = time1
     endif
     if(movetime .eq. 1) then 
        ceop_struc(n)%ceoptime1 = ceop_struc(n)%ceoptime2
        do f=1,8
           do t=1,LIS_rc%ngrid(n)
              ceop_struc(n)%metdata1(f,t) = ceop_struc(n)%metdata2(f,t)
           enddo
        enddo
     endif

     if(ceop_struc(n)%findtime2.eq.1) then 
        write(LIS_logunit,*) 'reading time2 data...'
        order = 2
        call read_ceop(n,findex,order)
        ceop_struc(n)%ceoptime2 = time2
     endif

  endif
end subroutine get_ceop


