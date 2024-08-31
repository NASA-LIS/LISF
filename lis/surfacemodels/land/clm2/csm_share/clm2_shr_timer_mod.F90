!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module clm2_shr_timer_mod

   !----------------------------------------------------------------------------
   !
   ! routines that support multiple CPU timers via F90 intrisics
   !
   ! Note: 
   ! o if   an operation is requested on an invalid timer number n
   !   then nothing is done in a routine
   ! o if   more than max_timers are requested, 
   !   then timer n=max_timers is "overloaded" and becomes invalid/undefined
   !----------------------------------------------------------------------------

   use clm2_shr_kind_mod

   implicit none

   private  ! resticted access
   public  :: clm2_shr_timer_init , clm2_shr_timer_get      , &
   &          clm2_shr_timer_start, clm2_shr_timer_stop     , &
   &          clm2_shr_timer_print, clm2_shr_timer_print_all, &
   &          clm2_shr_timer_check, clm2_shr_timer_check_all, &
   &          clm2_shr_timer_zero , clm2_shr_timer_zero_all , &
   &          clm2_shr_timer_free , clm2_shr_timer_free_all , &
   &          clm2_shr_timer_sleep

   integer(clm2_shr_kind_in),parameter :: stat_free    = 0  ! timer status constants
   integer(clm2_shr_kind_in),parameter :: stat_inuse   = 1
   integer(clm2_shr_kind_in),parameter :: stat_started = 2
   integer(clm2_shr_kind_in),parameter :: stat_stopped = 3
   integer(clm2_shr_kind_in),parameter :: max_timers   = 100 ! max number of timers

   integer(clm2_shr_kind_in) :: status (max_timers) ! status of each timer
   integer(clm2_shr_kind_in) :: cycles1(max_timers) ! cycle number at timer start 
   integer(clm2_shr_kind_in) :: cycles2(max_timers) ! cycle number at timer stop  
   character   (len=80) :: name   (max_timers) ! name assigned to each timer
   real   (clm2_shr_kind_r8) :: dt     (max_timers) ! accumulated time
   integer(clm2_shr_kind_in) :: calls  (max_timers) ! # of samples in accumulation
   integer(clm2_shr_kind_in) :: cycles_max = -1     ! max cycles before wrapping
   real   (clm2_shr_kind_r8) :: clock_rate          ! clock_rate: seconds per cycle

   save

!===============================================================================

   contains

!===============================================================================

subroutine clm2_shr_timer_init

   !----- local -----
   integer(clm2_shr_kind_in) :: cycles ! count rate return by system clock

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_init) ',a,i5)"

!-------------------------------------------------------------------------------
!
! This routine initializes:
! 1) values in all timer array locations
! 2) machine parameters necessary for computing cpu time from F90 intrinsics.
!    F90 intrinsic: system_clock(count_rate=cycles, count_max=cycles_max)
!-------------------------------------------------------------------------------

   call clm2_shr_timer_free_all

   call system_clock(count_rate=cycles, count_max=cycles_max)

   if (cycles /= 0) then
     clock_rate = 1.0/real(cycles)
   else
     clock_rate = 0
     write(6,F00) 'ERROR: no system clock available' 
   endif

end subroutine clm2_shr_timer_init

!===============================================================================

subroutine clm2_shr_timer_get(n, str)

   !----- arguments -----
   integer(clm2_shr_kind_in),intent(out) :: n    ! timer number 
   character (*)       ,intent( in) :: str  ! text string with timer name

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_get) ',a,i5)"

!-----------------------------------------------------------------------
!
!  search for next free timer
!
!-----------------------------------------------------------------------

   do n=1,max_timers
     if (status(n) == stat_free) then
       status(n) = stat_inuse
       name  (n) = str
       calls (n) = 0
       return
     endif
   end do

   n=max_timers
   name  (n) = "<invalid - undefined - overloaded>"
   write(6,F00) 'ERROR: exceeded maximum number of timers'

end subroutine clm2_shr_timer_get

!===============================================================================

subroutine clm2_shr_timer_start(n)

   !----- arguments -----
   integer(clm2_shr_kind_in), intent(in) :: n      ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_start) ',a,i5)"

!-----------------------------------------------------------------------
!
!  This routine starts a given timer.
!
!-----------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) call clm2_shr_timer_stop(n)

     status(n) = stat_started
     call system_clock(count=cycles1(n))
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine clm2_shr_timer_start
 
!===============================================================================

subroutine clm2_shr_timer_stop(n)

   !----- arguments -----
   integer(clm2_shr_kind_in), intent(in) :: n  ! timer number

   !----- local -----
!   real (clm2_shr_kind_r8) :: elapse      ! elapsed time returned by system counter

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_stop) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine stops a given timer, checks for cycle wrapping, computes the 
!  elpased time, and accumulates the elpased time in the dt(n) array
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if ( status(n) == stat_started) then
       call system_clock(count=cycles2(n))
       if (cycles2(n) >= cycles1(n)) then
         dt(n) = dt(n) + clock_rate*(cycles2(n) - cycles1(n))
       else
         dt(n) = dt(n) + clock_rate*(cycles_max + cycles2(n) - cycles1(n))
       endif
       calls (n) = calls (n) + 1
       status(n) = stat_stopped
     end if
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine clm2_shr_timer_stop
 
!===============================================================================

subroutine clm2_shr_timer_print(n)

   !----- arguments -----
   integer(clm2_shr_kind_in), intent(in) :: n     ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_print) ',a,i5)"
   character(len=*),parameter :: F01 = "('(clm2_shr_timer_print) timer',i3,&
   &                                     ':',i8,' calls,',f10.3,'s, id: ',a)"
!-------------------------------------------------------------------------------
!
!  prints the accumulated time for a given timer
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) then
       call clm2_shr_timer_stop(n)
       write (6,F01) n,dt(n),trim(name(n))
       call clm2_shr_timer_start(n)
     else
       write (6,F01) n,calls(n),dt(n),trim(name(n))
     endif
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine clm2_shr_timer_print

!===============================================================================

subroutine clm2_shr_timer_print_all

   !----- local -----
   integer(clm2_shr_kind_in) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_print_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  prints accumulated time for all timers in use
!
!-------------------------------------------------------------------------------

   write(6,F00) 'print all timing info:'

   do n=1,max_timers
     if (status(n) /= stat_free) call clm2_shr_timer_print(n)
   end do

end subroutine clm2_shr_timer_print_all

!===============================================================================

subroutine clm2_shr_timer_zero(n)

   !----- arguments -----
   integer(clm2_shr_kind_in), intent(in) :: n       ! timer number
   
   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_zero) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine resets a given timer.
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     dt(n) = 0.0
     calls(n) = 0
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine clm2_shr_timer_zero

!===============================================================================

subroutine clm2_shr_timer_zero_all

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_zero_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine resets all timers.
!
!-------------------------------------------------------------------------------

   dt = 0.0
   calls = 0

end subroutine clm2_shr_timer_zero_all

!===============================================================================

subroutine clm2_shr_timer_check(n)

   !----- arguments -----
   integer(clm2_shr_kind_in), intent(in) ::  n   ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_check) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine checks a given timer.  This is primarily used to
!  periodically accumulate time in the timer to prevent timer cycles
!  from wrapping around max_cycles.
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) then
       call clm2_shr_timer_stop (n)
       call clm2_shr_timer_start(n)
     endif
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine clm2_shr_timer_check

!===============================================================================

subroutine clm2_shr_timer_check_all

   !----- local -----
   integer(clm2_shr_kind_in) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_check_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  Call clm2_shr_timer_check for all timers in use
!
!-------------------------------------------------------------------------------

   do n=1,max_timers
     if (status(n) == stat_started) then
       call clm2_shr_timer_stop (n)
       call clm2_shr_timer_start(n)
     endif
   end do

end subroutine clm2_shr_timer_check_all

!===============================================================================

subroutine clm2_shr_timer_free(n)

   !----- arguments -----
   integer(clm2_shr_kind_in),intent(in) :: n    ! timer number 

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_free) ',a,i5)"

!-----------------------------------------------------------------------
!
!  initialize/free all timer array values
!
!-----------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     status (n) = stat_free
     name   (n) = "<invalid - undefined>"
     dt     (n) = 0.0
     cycles1(n) = 0
     cycles2(n) = 0
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine clm2_shr_timer_free

!===============================================================================

subroutine clm2_shr_timer_free_all

   !----- local -----
   integer(clm2_shr_kind_in) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_timer_free_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  initialize/free all timer array values
!
!-------------------------------------------------------------------------------

   do n=1,max_timers
     call clm2_shr_timer_free(n)
   end do

end subroutine clm2_shr_timer_free_all

!===============================================================================

subroutine clm2_shr_timer_sleep(sec)

   use clm2_shr_sys_mod     ! share system calls (namely, shr_sys_sleep)

   !----- local -----
   real   (clm2_shr_kind_r8),intent(in) :: sec  ! number of seconds to sleep

!-------------------------------------------------------------------------------
! Sleep for approximately sec seconds
!
! Note: sleep is typically a system call, hence it is implemented in 
!       clm2_shr_sys_mod, although it probably would only be used in a timing 
!       context, which is why there is a clm2_shr_timer_* wrapper provided here.
!-------------------------------------------------------------------------------

   call clm2_shr_sys_sleep(sec)

end subroutine clm2_shr_timer_sleep

!===============================================================================
end module clm2_shr_timer_mod
!===============================================================================
