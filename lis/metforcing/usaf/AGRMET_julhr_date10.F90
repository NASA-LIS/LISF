!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: AGRMET_julhr_date10
!  \label{AGRMET_julhr_date10}
!
! !REVISION HISTORY:
!     15 oct 1998  initial version.........................mr moore/dnxm
!     10 aug 1999  ported to ibm sp2.  added intent attributes to
!                  arguments...............................mr gayno/dnxm
!     29 oct 2005  Sujay Kumar, Adopted in LIS
!     11 Mar 2010  Changed program names in messages to LIS..............
!                  ............................Chris Franks/16WS/WXE/SEMS
!
! !INTERFACE:    
subroutine AGRMET_julhr_date10( julhr, date10)
! !USES: 
  use LIS_timeMgrMod, only : LIS_tmjul4
  use LIS_logMod,     only : LIS_abort, LIS_endrun

  implicit none 
! !ARGUMENTS: 
  character(len=*), intent(out)  :: date10
  integer,          intent(in)   :: julhr   

! !DESCRIPTION:
!     to convert from a julian hour to a 10 digit date/time group
!     (yyyymmddhh)
!  
!     \textbf{Method}
!    
!     - call utility routine LIS\_tmjul4 to convert from julian hours
!       to year, month, day, hour. \newline
!     - perform several checks to ensure date information passed
!       back from LIS\_tmjul4 is valid. \newline
!     - if date is good, convert from integer data to a 10 digit
!       date/time group and pass back to calling routine. \newline
!
! The arguments are variables are: 
! \begin{description}
!   \item[date10]          output 10-digit character date time group
!                     (yyyymmddhh)
!   \item[dd]              day of the month
!   \item[hh]              time of day in hours
!   \item[j]               loop counter
!   \item[julhr]           input julian hour
!   \item[mm]              month of the year
!   \item[yyyy]            four digit year
!   \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tmjul4](\ref{LIS_tmjul4}) \newline
!   convert julian hour to hour,day,month and year
!  \end{description}
!EOP
  character*255                  :: message ( 20 )
  integer                        :: dd     
  integer                        :: hh     
  integer                        :: j      
  integer                        :: mm       
  integer                        :: yyyy    
!     ------------------------------------------------------------------
!     executable code begins here... use LIS_tmjul4 to convert julhr to 
!     hour, day, month and year
!     ------------------------------------------------------------------

  call LIS_tmjul4( hh, dd, mm, yyyy, julhr )
 
!     ------------------------------------------------------------------
!     check for valid hour, day, month, and year
!     ------------------------------------------------------------------ 

  if( (  hh .lt.    0 .or.   hh .gt.   23) .or. &
       (  dd .lt.    1 .or.   dd .gt.   31) .or. &
       (  mm .lt.    1 .or.   mm .gt.   12) .or. &
       (yyyy .lt. 1968 .or. yyyy .gt. 9999) )then
     
     message(1) = 'program: LIS'
     message(2) = '  routine AGRMET_julhr_date10'
     message(3) = '  invalid julhr to date/time conversion'
     
     call lis_abort( message )
     call LIS_endrun()
     return
  endif

!     ------------------------------------------------------------------ 
!     output errors:
!     flag invalid day/month combinations, for example...
!     if february check for correct number of days in month
!     ------------------------------------------------------------------ 
    
  if( mm .eq. 2 )then
     
     if( (mod((yyyy - 68), 4) .ne. 0) .and. (dd .gt. 28) )then
        
        message(1) = 'program: LIS'
        message(2) = '  routine AGRMET_julhr_date10'
        message(3) = '  oops, created 29 feb in non-leap year'
        
        call lis_abort( message )
        call LIS_endrun        
     endif
     
  elseif (  (dd .gt. 30) .and. &
       ( (mm .eq. 4)  .or. &
       (mm .eq. 6)  .or. &
       (mm .eq. 9)  .or. &
       (mm .eq. 11) ) ) then  
     
     message(1) = 'program: LIS'
     message(2) = '  routine AGRMET_julhr_date10'
     message(3) = '  oops, created 31st day for month'
     message(4) = '  with only 30 days'
     
     call lis_abort( message )
     call LIS_endrun
     
  endif
  
!    ------------------------------------------------------------------
!     convert integer values into a 10 character date/time string
!     ------------------------------------------------------------------

  write(date10(1:4) ,'(i4)', err=9300) yyyy
  write(date10(5:6) ,'(i2)', err=9300) mm
  write(date10(7:8) ,'(i2)', err=9300) dd
  write(date10(9:10),'(i2)', err=9300) hh
  
  do j = 1, 10
     if( date10(j:j) .eq. ' ' )then
        date10(j:j) = '0'
     endif
  enddo
  
  return
      
!     ------------------------------------------------------------------
!     error handling section.
!     ------------------------------------------------------------------

9300 continue
  
  message(1) = 'program: LIS'
  message(2) = '  routine AGRMET_julhr_date10'
  message(3) = '  error converting integer date/times to character'
  
  call lis_abort( message )
  call LIS_endrun
  
end subroutine AGRMET_julhr_date10
