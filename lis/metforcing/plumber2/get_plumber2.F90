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
! !ROUTINE: get_plumber2
! \label{get_plumber2}
!
! !REVISION HISTORY:
! 15 Sep 2021: Mark Beauharnois, Derived from Bondville and PUMET readers
!
! !INTERFACE:
subroutine get_plumber2(n,findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_get_nstep, LIS_tick
  use LIS_logMod, only : LIS_logunit
  use plumber2_forcingMod, only : plumber2_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Opens, reads, and interpolates the PLUMBER2 station data.
!  At the beginning of a simulation, the code reads the most
!  recent past data, and the nearest future data.  These two
!  datasets are used to temporally interpolate the data to
!  the current model timestep.
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
!    determines the PLUMBER2 data times
!  \item[read\_plumber2](\ref{read_plumber2}) \newline
!    Interpolates the appropriate PLUMBER2 station data to LIS grid
!  \end{description}
!EOP

  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  integer :: ferror
  real :: ts1, ts2
  real :: gmt1,gmt2
  real*8 :: timenow,time1,time2
  integer :: movetime       ! if 1=move time 2 data into time 1
  integer :: nstep
  integer :: f,t
  integer :: order

  ! write(LIS_logunit,*) '[DIAG] Starting get_plumber2'

  plumber2_struc(n)%findtime1=0
  plumber2_struc(n)%findtime2=0
  movetime = 0


  if(LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1 .or. &
     plumber2_struc(n)%reset_flag) then

     plumber2_struc(n)%findtime1 = 1
     plumber2_struc(n)%findtime2 = 1

     LIS_rc%rstflag(n) = 0
     plumber2_struc(n)%reset_flag = .false.

     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = LIS_rc%mn
     ss1 = 0
     ts1 = plumber2_struc(n)%ts

     call LIS_tick(plumber2_struc(n)%ringtime,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
     !write(LIS_logunit,*) '[DIAGM] plumber2_struc ringtime: ',&
     !       plumber2_struc(n)%ringtime
  endif

  ! Current time:
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0

  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,0.0)

  !write(LIS_logunit,*) '[DIAGM] timenow: ',timenow
  
  !! If current time >= time for when to move bookend, move bookend 2 to
  !! bookend 1

   if(timenow >= plumber2_struc(n)%ringtime) then
     plumber2_struc(n)%findtime2 = 1
     if(plumber2_struc(n)%findtime1 == 0) then
       movetime = 1
     endif

     !! reset ringtime to next PLUMBER2 ts increment
     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = LIS_rc%mn
     ss1 = 0
     ts1 = plumber2_struc(n)%ts

     call LIS_tick(plumber2_struc(n)%ringtime,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
   endif

   if(plumber2_struc(n)%findtime1 == 1) then
     !! determine PLUMBER2 forcing 1 time, ONLY first time step
     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = LIS_rc%mn
     ss1 = 0
     ts1 = 0
     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
   endif

   if(plumber2_struc(n)%findtime2 == 1) then
     !! determine PLUMBER2 forcing 2 time, first time step AND each time
     !! ringtime is exceeded
     yr2 = LIS_rc%yr
     mo2 = LIS_rc%mo
     da2 = LIS_rc%da
     hr2 = LIS_rc%hr
     mn2 = LIS_rc%mn
     ss2 = 0
     ts2 = plumber2_struc(n)%ts   !! 'ts' seconds ahead
     call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
   endif

   !Time to read...
   if (plumber2_struc(n)%findtime1.eq.1) then
        write(LIS_logunit,*) '[INFO] Reading PLUMBER2 time1 data...'
        write(LIS_logunit,*) '[INFO] time1: ',time1
        write(LIS_logunit,*) '[INFO] plumber2_struc%ts: ',plumber2_struc(n)%ts
        order = 1
        plumber2_struc(n)%read_index = plumber2_struc(n)%read_index + 1
        call read_plumber2(n, order, findex, &
                plumber2_struc(n)%plumber2file, &
                plumber2_struc(n)%metdata1, ferror)
        plumber2_struc(n)%plumber2time1 = time1
   endif

   if (movetime.eq.1) then 
      plumber2_struc(n)%metdata1 = plumber2_struc(n)%metdata2
      plumber2_struc(n)%metdata2 = LIS_rc%udef

      plumber2_struc(n)%plumber2time1 = &
            plumber2_struc(n)%plumber2time2
   endif

   if (plumber2_struc(n)%findtime2.eq.1) then
      write(LIS_logunit,*) '[INFO] Reading PLUMBER2 time2 data...'
      write(LIS_logunit,*) '[INFO] time2: ',time2
      write(LIS_logunit,*) '[INFO] plumber2_struc%ts: ',plumber2_struc(n)%ts
      order = 2
      plumber2_struc(n)%read_index = plumber2_struc(n)%read_index + 1
      call read_plumber2(n, order, findex, &
              plumber2_struc(n)%plumber2file, &
              plumber2_struc(n)%metdata2, ferror)
      plumber2_struc(n)%plumber2time2 = time2
  endif

end subroutine get_plumber2
