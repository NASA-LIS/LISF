!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_FASSTsingle
! \label{get_FASSTsingle}
!
! !REVISION HISTORY:
! 13 Apr 2007: Bailing Li, Initial Specification
! 21 Oct 2010: David Mocko, Updated for FASST single-point test case
!
! !INTERFACE:
subroutine get_FASSTsingle(n,findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_logMod, only : LIS_logunit
  use FASSTsingle_forcingMod, only : FASSTsingle_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
  real(kind=8) :: metdata3(LIS_rc%met_nf(findex),LIS_rc%ngrid(n))
!
! !DESCRIPTION:
!  Opens, reads, and interpolates the FASSTsingle station data.
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
!    determines the FASSTsingle data times
!  \item[read\_FASSTsingle](\ref{read_FASSTsingle}) \newline
!    Interpolates the appropriate FASSTsingle station data to LIS grid
!  \end{description}
!EOP

  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2
  real :: gmt1,gmt2
  real*8 :: timenow,time1,time2
  integer :: movetime       ! if 1=move time 2 data into time 1
  integer :: f,t

  !      write(LIS_logunit,*) 'starting get_FASSTsingle'
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
  mn1 = 30
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

  FASSTsingle_struc(n)%findtime1 = 0
  FASSTsingle_struc(n)%findtime2 = 0
  !      write(LIS_logunit,*) 'timenow: ',timenow
  !      write(LIS_logunit,*) 'time1: ',time1
  !      write(LIS_logunit,*) 'time2: ',time2
  !      write(LIS_logunit,*) 'starttime: ',FASSTsingle_struc(n)%starttime
  if ((timenow.ge.FASSTsingle_struc(n)%starttime).and.             &
       (.not.FASSTsingle_struc(n)%startRead)) then
     FASSTsingle_struc(n)%findtime1 = 1
     FASSTsingle_struc(n)%findtime2 = 1
     FASSTsingle_struc(n)%startRead = .true.
     movetime = 0
  endif
  if (FASSTsingle_struc(n)%startRead) then
     if (timenow.ge.FASSTsingle_struc(n)%FASSTsingletime2) then
        movetime = 1
        FASSTsingle_struc(n)%findtime2 = 1
     endif

     !Time to open file and start reading.
     !keep on reading until the obstime is reached.
     if (FASSTsingle_struc(n)%findtime1.eq.1) then
        write(LIS_logunit,*) 'reading time1 data...'
        call read_FASSTsingle(n,222,findex,metdata3,1)
        FASSTsingle_struc(n)%FASSTsingletime1 = time1
        do f = 1,35
           do t = 1,LIS_rc%ngrid(n)
              FASSTsingle_struc(n)%metm_back(f) = metdata3(f,t)
           enddo
        enddo
     endif
     if (movetime.eq.1) then 
        FASSTsingle_struc(n)%FASSTsingletime1 =                    &
             FASSTsingle_struc(n)%FASSTsingletime2
        do f = 1,35
           do t = 1,LIS_rc%ngrid(n)
              FASSTsingle_struc(n)%metm_back(f) =                  &
                   FASSTsingle_struc(n)%metdata4(f)
           enddo
        enddo
     endif

     if (FASSTsingle_struc(n)%findtime2.eq.1) then
        write(LIS_logunit,*) 'reading time2 data...'
        call read_FASSTsingle(n,222,findex,metdata3,2)
        FASSTsingle_struc(n)%FASSTsingletime2 = time2
        do f = 1,35
           do t = 1,LIS_rc%ngrid(n)
              FASSTsingle_struc(n)%metdata4(f) = metdata3(f,t)
           enddo
        enddo
     endif
  endif

end subroutine get_FASSTsingle

