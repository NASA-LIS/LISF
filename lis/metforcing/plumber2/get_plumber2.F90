!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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

  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2
  integer :: ferror
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

  yr1 = LIS_rc%yr           !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0

  !      write(LIS_logunit,*) [DIAG] yr1,mo1,da1,hr1,mn1,ss1
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,real(ts1))

  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr

!!!! MCB Note. 'mn1' here most likely needs to be something like:
!!!! mn1 = plumber2_struc(n)%ts / 60 since PLUMBER2 data could be
!!!!  30min or 1hr data
  mn1 = plumber2_struc(n)%ts / 60
  ss1 = 0
  ts1 = 0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,real(ts1))

  yr2 = LIS_rc%yr           !next hour
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = LIS_rc%hr
  mn2 = LIS_rc%mn
  ss2 = 0

!!!! MCB Note. 'ts2' here most likely needs to be something like:
!!!! ts2 = 60 * (plumber2_struc(n)%ts / 60)
  ts2 = 60 * (plumber2_struc(n)%ts / 60)
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,real(ts2))

  plumber2_struc(n)%findtime1 = 0
  plumber2_struc(n)%findtime2 = 0
  write(LIS_logunit,*) '[DIAG] PLUMBER2 get_plumber2:'
  write(LIS_logunit,*) '[DIAG] timenow: ',timenow
  write(LIS_logunit,*) '[DIAG] time1: ',time1
  write(LIS_logunit,*) '[DIAG] time2: ',time2

  if ((timenow.ge.plumber2_struc(n)%starttime).and.              &
       (.not.plumber2_struc(n)%startRead)) then
     plumber2_struc(n)%findtime1 = 1
     plumber2_struc(n)%findtime2 = 1
     plumber2_struc(n)%startRead = .true.
     movetime = 0
  endif

  if (plumber2_struc(n)%startRead) then
     if (timenow.ge.plumber2_struc(n)%plumber2time2) then
        movetime = 1
        plumber2_struc(n)%findtime2 = 1
     endif

     !Time to open file and start reading.
     !keep on reading until the obstime is reached.
     if (plumber2_struc(n)%findtime1.eq.1) then
        write(LIS_logunit,*) '[INFO] Reading PLUMBER2 time1 data...'
        order = 1
        !!call read_plumber2(n,222,findex,order,1)
        !! MCB Note: New PLUMBER2 reader based on PUMET's reader
        call read_plumber2(n, order, findex, &
                plumber2_struc(n)%plumber2file, &
                plumber2_struc(n)%metdata1, ferror)
        plumber2_struc(n)%plumber2time1 = time1
     endif
     if (movetime.eq.1) then 
        plumber2_struc(n)%plumber2time1 =                      &
             plumber2_struc(n)%plumber2time2

!       MCB Note: LIS_rc%met_nf(findex) is set in plumber2_forcingMod.F90 'init_plumber2'
!       It's used here so that the number of PLUMBER2 met variables used is not hardcoded
        do f = 1,LIS_rc%met_nf(findex)
           do t = 1,LIS_rc%ngrid(n)
              plumber2_struc(n)%metdata1(f,t) = plumber2_struc(n)%metdata2(f,t)
           enddo
        enddo
     endif

     if (plumber2_struc(n)%findtime2.eq.1) then
        write(LIS_logunit,*) '[INFO] Reading PLUMBER2 time2 data...'
        order = 2
        !! call read_plumber2(n,222,findex,order,2)
        !! MCB Note: New PLUMBER2 reader based on PUMET's reader
        call read_plumber2(n, order, findex, &
                plumber2_struc(n)%plumber2file, &
                plumber2_struc(n)%metdata2, ferror)
        plumber2_struc(n)%plumber2time2 = time2
     endif
  endif

end subroutine get_plumber2

