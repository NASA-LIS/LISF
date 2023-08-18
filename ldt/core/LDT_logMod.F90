!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_logMod
!BOP
! !MODULE: LDT_logMod
!
! !DESCRIPTION:
!  The code in this file provides routines for diagnostic and error
!  logging in LDT. A file unit number is reserved for the diagnostic 
!  log file. When a multiprocessor run is conducted, a diagnostic file
!  from each processor is written to a separate file. An entry to the 
!  diagnostic log can be made by using: 
!  \begin{verbatim}
!   write(LDT_logunit,*) msg
!  \end{verbatim}
! 
! !REVISION HISTORY:
!  02 Oct 2008: Sujay Kumar; Initial version
!
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_log_init     ! initializes LDT logging 
  public :: LDT_log_msg      ! writes a time-processor stamped message 
                             !  to the diagnostic log
  public :: LDT_abort        ! generates a standard abort message
  public :: LDT_alert        ! generates a standard alert message
  public :: LDT_endrun       ! ends the program
  public :: LDT_verify       ! checks a return a code for success
  public :: LDT_warning      ! checks a return code and issue a warning
  public :: LDT_getNextUnitNumber   ! get the next available unit number
  public :: LDT_releaseUnitNumber   ! release the unit number
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LDT_logunit ! file unit number used for diagnostic logging
!EOP  

  integer :: LDT_logunit 
  logical :: IOU(99)  = .false. ! I/O file unit numbers
  
!BOP 
! 
! !ROUTINE: LDT_verify
! \label{LDT_verify}
! 
! !INTERFACE:
  interface LDT_verify
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure check_error
     module procedure verify
  end interface
!EOP

contains
!BOP
! !ROUTINE: LDT_log_init
! \label{LDT_log_init}
!
! !INTERFACE: 
  subroutine LDT_log_init(iunit)
    integer, intent(IN) :: iunit
! 
! !DESCRIPTION:
! Initializes the lis log management. Reserves
! the assigned file unit number for lis diagnostic
! logging. 
!
!EOP    
    LDT_logunit = iunit

  end subroutine LDT_log_init


!BOP
! 
! !ROUTINE: LDT_log_msg
! \label{LDT_log_msg}
! 
! !INTERFACE:
  subroutine LDT_log_msg(msg)      
    implicit none
      
! !ARGUMENTS: 
    character(len=*), intent(in) :: msg
! 
! !DESCRIPTION: 
! 
!  This routine formats a given message by prepending a
!  time stamp and appending the process id number.  This newly formatted
!  message is then written to standard out.
!
!  The arguments are: 
!  \begin{description}
!   \item [msg]
!     message to be written to the logfile
!  \end{description}
!
!EOP
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    integer, dimension(8) :: values
    
    call date_and_time(date,time,zone,values)
    write(unit=LDT_logunit,fmt=*)date(1:4),'-',date(5:6),'-',date(7:8),'T',  &
         time(1:2),':',time(3:4),':',time(5:10),' ', &
         trim(msg)
    
  end subroutine LDT_log_msg
  

!BOP
!
! !ROUTINE: LDT_abort
! \label{LDT_abort}
!       
! !REVISION HISTORY:
!     03 mar 1999 initial version .........................mr moore/dnxm
!     09 aug 1999 ported to IBM SP-2. increase size of abort message
!                 to 100 characters. removed stop 255 with a call
!                 to IBM utility "abort"...................mr gayno/dnxm     
!     29 oct 2005  Sujay Kumar, Adopted in LDT
!
! !INTERFACE:    
  subroutine LDT_abort( abort_message )

    implicit none
    
    character*100  :: abort_message(20)
    
! !DESCRIPTION:
!
!  to generate standard abort message and abort calling program.
!
!  \begin{itemize}
!   \item{
!   open abort file, write abort message and then call
!   the system abort routine. this will cause program termination
!   and cause script calling this routine to abort.}
!   \item{if there is an error opening or writing abort file, send
!    message to unit 6 and call the system abort routine}
!  \end{itemize}
!
!   The arguments are: 
!   \begin{description}
!   \item [abort\_message]
!     abort message
!  \end{description}
!
!EOP
    character*7                :: iofunc
    character*13               :: message_file
    integer                    :: i
    integer                    :: ftn 
    integer                    :: istat
    
!    ------------------------------------------------------------------
!     executable code begins here ...
!     open and write message to message file
!     ------------------------------------------------------------------

    message_file = 'abort_message'
    
    iofunc = 'opening'
    
    ftn = LDT_getNextUnitNumber()
    open (unit     = ftn, &
         file     = message_file,&
         err      = 9000, &
         iostat   = istat)
    iofunc = 'writing'
    
    do i = 1, 20
       write (ftn, 6000, err = 9000) trim(abort_message(i))
    enddo
    
    call LDT_releaseUnitNumber(ftn)
    
    call LDT_endrun
      
!     ----------------------------------------------------------------
!     format statements.
!     ----------------------------------------------------------------
                      
6000 format (a)   
    
8000 format ( /, 1x, 62('-'),    &
          /, 1x, 'program: unknown', &
          /, 3x, 'subroutine: agr_abort', &
          /, 5x, 'error while ', a7, ' file -> abort_message',&
          /, 5x, 'istat = ', i6, &
          /, 5x, 'unable to abort normally.',&
          /, 1x, 62('-'))
  
!     ----------------------------------------------------------------
!     error handling section.
!     ----------------------------------------------------------------

9000 write (6, 8000) iofunc, istat
      
    call LDT_endrun
  
  end subroutine LDT_abort

!BOP
!
! !ROUTINE: LDT_alert
! \label{LDT_alert}
!
! !REVISION HISTORY:
!
!     03 mar 1999 initial version .........................mr moore/dnxm
!     09 aug 1999 ported to IBM SP-2.  increase size of message to
!                 100 characters.  added intent attributes to 
!                 arguments................................mr gayno/dnxm     
!     29 oct 2005  Sujay Kumar, Adopted in LDT
!
! !INTERFACE:    
  subroutine LDT_alert( program_name, alert_number, message )
!EOP    
    implicit none
    
    character(len=*),  intent(in)    :: program_name  
    integer,           intent(in)    :: alert_number
    character(len=*), intent(inout)  :: message(20)
    
!
! !DESCRIPTION:
!
!     to generate alert messages when model is degraded but
!     not to the point where we would want to abort.
!
!     \begin{itemize}
!     \item{write alert message to file.}
!     \item{if there is an error opening or writing alert message file, 
!       call abort routine}
!
!   The arguments are: 
!   \begin{description}
!   \item [program\_name]
!     name of the alerting program 
!   \item [alert\_number]
!     number of the current alert
!   \item [message]
!     alert message
!  \end{description}
!
!EOP
    character*3                   :: calert_number
    character*7                   :: iofunc
    character*37                  :: message_file
    integer                       :: i
    integer                       :: istat
    integer                       :: ftn
    
!     ------------------------------------------------------------------
!     executable code starts here ... set output message file name
!     ------------------------------------------------------------------

    write (calert_number, '(i3.3)', iostat = istat) alert_number
    
    if (istat .ne. 0) then
       calert_number = 'xx'
    endif
    
    message_file = 'alert.' // trim(program_name) // &
         '.' // calert_number
      
!     ------------------------------------------------------------------
!     open file
!     ------------------------------------------------------------------

    iofunc = 'opening'
    ftn = LDT_getNextUnitNumber()
    open (unit   = ftn, &
         file   = trim(message_file), &
         iostat = istat)
    
!     ------------------------------------------------------------------
!     write to file and close
!     ------------------------------------------------------------------
    
    iofunc = 'writing'
    
    do i = 1, 20
       write (ftn, 8000, iostat=istat) trim(message(i))
    enddo
    
    call LDT_releaseUnitNumber(ftn)
    
    return
  
!     ----------------------------------------------------------------
!     format statements
!     ----------------------------------------------------------------
    
8000 format (a)   
      
  end subroutine LDT_alert

!BOP
! !ROUTINE: check_error
! 
! !INTERFACE:
  subroutine check_error(ierr,msg)
! 
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ierr
    character(len=*), intent(in) :: msg  
! 
! !DESCRIPTION:
! This is an error check routine. Program exits in case of error
! with an associated error message written to the 'abort message'
! file. 
! 
!EOP

    if ( ierr /= 0 ) then
       write(LDT_logunit,*) "[ERR] ",msg," Stopping.",ierr
       call LDT_endrun
    endif

  end subroutine check_error

!BOP
! !ROUTINE: verify
! \label{verify}
! 
! !INTERFACE:
  subroutine verify(ierr)
! 
! !ARGUMENTS:

    implicit none
    integer, intent(in)  :: ierr
! 
! !DESCRIPTION:
! This is an error check routine. Program exits in case of error
! with an associated error message written to the 'abort message'
! file. 
! 
!EOP

    if( ierr /= 0 ) then
       write(LDT_logunit,*) "[ERR] Error from an LDT Subroutine Call, Program Stopping ...",ierr
       call LDT_endrun
    endif

  end subroutine verify

!BOP
! !ROUTINE: LDT_warning
! 
! !INTERFACE:
  subroutine LDT_warning(ierr,msg)
! 
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ierr
    character(len=*), intent(in) :: msg
! 
! !DESCRIPTION:
! This is an error check routine. Program issues a warning in case of error
! with an associated error message.
! 
!EOP
    if ( ierr /= 0 ) then
!       write(LDT_logunit,*) '****************WARNING*********************'
       write(LDT_logunit,*) msg
!       write(LDT_logunit,*) '****************WARNING*********************'
    endif

  end subroutine LDT_warning

!BOP
!
! !ROUTINE: LDT_getNextUnitNumber
! \label{LDT_getNextUnitNumber}
!    
! !INTERFACE: 
    integer function LDT_getNextUnitNumber()
      
! !DESCRIPTION: 
! get next available Fortran unit number
!
! Method: 
! Get next available Fortran unit number itst. This routine is
! modified from the CCSM codes
!EOP
      integer :: i 
! starting with 13 because 12 is reserved for grib
      do i = 13, 99
         if(.not.IOU(i)) then 
            LDT_getNextUnitNumber = i
            IOU(i) = .true.
            RETURN
         endif
      enddo
      write(*,*) 'LDT_fileUnitsMod: Ran out of Unit numbers'
    end function LDT_getNextUnitNumber
    
!BOP
! 
! !ROUTINE: LDT_releaseUnitNumber
! \label{LDT_releaseUnitNumber}
! 
! !INTERFACE: 
    subroutine LDT_releaseUnitNumber(iunit)
!
! !DESCRIPTION: 
! 
! close and release Fortran unit no longer in use
!
!EOP
      integer, intent(in)    :: iunit
      
      if(.not.IOU(iunit)) then 
         write(*,*) 'file unit release error: ',iunit,'is not flaggged as in use'
      endif
      close(iunit)
      IOU(iunit) = .false.
      return
    end subroutine LDT_releaseUnitNumber

!BOP
!
! !ROUTINE: LDT_endrun
! \label{LDT_endrun}
!
! !REVISION HISTORY: 
!  14 Nov 2002    Sujay Kumar  Initial Specification
!
! !INTERFACE:
  subroutine LDT_endrun

! !USES: 
#if (defined SPMD)
    use LDT_mpiMod
#endif
    implicit none

    integer :: ierr
!  
!  !DESCRIPTION: 
!  Routine to be called to terminate the program. This routine
!  flushes the output streams and aborts the mpi processes.
!
!EOP
    write(LDT_logunit,*)'endrun is being called'
    flush( LDT_logunit )   ! Flush all output to standard output

#if (defined SPMD) 
!    call mpi_abort (MPI_COMM_WORLD, 1)         ! LDT
    call mpi_abort (MPI_COMM_WORLD, 1, ierr)   ! LIS-7
#else
    error stop 1
#endif   
   end subroutine LDT_endrun

  end module LDT_logMod
