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
! !ROUTINE: get_Loobos
! \label{get_Loobos}
!
! !REVISION HISTORY:
! 07 Oct 2010: David Mocko, Updated for Loobos test case
!
! !INTERFACE:
subroutine get_Loobos(n,findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_logMod, only : LIS_logunit
  use Loobos_forcingMod, only : Loobos_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Opens, reads, and interpolates the Loobos station data.
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
!    determines the Loobos data times
!  \item[read\_Loobos](\ref{read_Loobos}) \newline
!    Interpolates the appropriate Loobos station data to LIS grid
!  \end{description}
!EOP

  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2
  real :: gmt1,gmt2
  real*8 :: timenow,time1,time2
  integer :: movetime       ! if 1=move time 2 data into time 1
  integer :: f,t
  integer :: order

  !      write(LIS_logunit,*) 'starting get_Loobos'
  movetime = 0

  yr1 = LIS_rc%yr           !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  !      write(LIS_logunit,*) yr1,mo1,da1,hr1,mn1,ss1
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,real(ts1))

  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,real(ts1))

  yr2 = LIS_rc%yr           !next hour
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = LIS_rc%hr
  mn2 = LIS_rc%mn
  ss2 = 0
  ts2 = 60*30
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,real(ts2))

  Loobos_struc(n)%findtime1 = 0
  Loobos_struc(n)%findtime2 = 0
  !write(LIS_logunit,*) 'timenow: ',timenow
  !write(LIS_logunit,*) 'time1: ',time1
  !write(LIS_logunit,*) 'time2: ',time2
  !write(LIS_logunit,*) 'starttime: ',Loobos_struc(n)%starttime
  if ((timenow.ge.Loobos_struc(n)%starttime).and.              &
       (.not.Loobos_struc(n)%startRead)) then
     Loobos_struc(n)%findtime1 = 1
     Loobos_struc(n)%findtime2 = 1
     Loobos_struc(n)%startRead = .true.
     movetime = 0
  endif
  if (Loobos_struc(n)%startRead) then
     if (timenow.ge.Loobos_struc(n)%Loobostime2) then
        movetime = 1
        Loobos_struc(n)%findtime2 = 1
     endif

     !Time to open file and start reading.
     !keep on reading until the obstime is reached.
     if (Loobos_struc(n)%findtime1.eq.1) then
        write(LIS_logunit,*) 'reading time1 data...'
        order = 1
        call read_Loobos(n,222,findex,order,1)
        Loobos_struc(n)%Loobostime1 = time1
     endif
     if (movetime.eq.1) then 
        Loobos_struc(n)%Loobostime1 =                      &
             Loobos_struc(n)%Loobostime2
        do f = 1,8
           do t = 1,LIS_rc%ngrid(n)
              Loobos_struc(n)%metdata1(f,t) = Loobos_struc(n)%metdata2(f,t)
           enddo
        enddo
     endif

     if (Loobos_struc(n)%findtime2.eq.1) then
        write(LIS_logunit,*) 'reading time2 data...'
        order = 2
        call read_Loobos(n,222,findex,order,2)
        Loobos_struc(n)%Loobostime2 = time2
     endif
  endif

end subroutine get_Loobos

