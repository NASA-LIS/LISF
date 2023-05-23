!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: LVT_logMod
! \label(LVT_logMod)
!
! !INTERFACE:
module LVT_logMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  The code in this file provides routines for diagnostic and error
!  logging in LVT. A file unit number is reserved for the diagnostic 
!  log file. When a multiprocessor run is conducted, a diagnostic file
!  from each processor is written to a separate file. An entry to the 
!  diagnostic log can be made by using: 
!  \begin{verbatim}
!   write(LVT_logunit,*) msg
!  \end{verbatim}
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  02 Oct 2008: Sujay Kumar; Initial version
! 
!EOP
!BOP
!
! 
!
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_log_init ! initializes LVT logging. 
  public :: LVT_log_msg ! writes a time-processor stamped message 
                        ! to the diagnostic log
  public :: LVT_abort  ! generates a standard abort message
  public :: LVT_alert  ! generates a standard alert message
  public :: LVT_getNextUnitNumber   !get the next available unit number
  public :: LVT_releaseUnitNumber   !release the unit number!
  public :: LVT_endrun
  public :: LVT_warning
  public :: LVT_verify
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LVT_logunit ! file unit number used for diagnostic logging
!EOP  

  integer, parameter :: stacksize= 100000
  integer :: LVT_logunit 
  logical :: IOU(stacksize)  = .false. ! I/O file unit numbers
  
!BOP
! 
! !ROUTINE: LVT_verify
! \label{LVT_verify}
!
! !INTERFACE:
  interface LVT_verify
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP 
! 
! 
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure check_error
     module procedure verify
  end interface
!EOP

contains
!BOP
! 
! !ROUTINE: LVT_log_init
! \label{LVT_log_init}
!
! !INTERFACE: 
  subroutine LVT_log_init(iunit)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    integer, intent(IN) :: iunit
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! Initializes the lis log management. Reserves
! the assigned file unit number for lis diagnostic
! logging. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    LVT_logunit = iunit

  end subroutine LVT_log_init


!BOP
! 
! !ROUTINE: LVT_log_msg
! \label{LVT_log_msg}
!
! !INTERFACE:
  subroutine LVT_log_msg(msg)      
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
    character(len=*), intent(in) :: msg
! 
! !OUTPUT PARAMETERS:
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
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    integer, dimension(8) :: values
    
    call date_and_time(date,time,zone,values)
    write(unit=LVT_logunit,fmt=*)date(1:4),'-',date(5:6),'-',date(7:8),'T',  &
         time(1:2),':',time(3:4),':',time(5:10),' ', &
         trim(msg)
    
  end subroutine LVT_log_msg
  

!BOP
! 
! !ROUTINE: LVT_abort
! \label{LVT_abort}
!
! !INTERFACE:    
  subroutine LVT_abort( abort_message )
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
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
! !FILES USED:
!
! !REVISION HISTORY:
!     03 mar 1999 initial version .........................mr moore/dnxm
!     09 aug 1999 ported to IBM SP-2. increase size of abort message
!                 to 100 characters. removed stop 255 with a call
!                 to IBM utility "abort"...................mr gayno/dnxm     
!     29 oct 2005  Sujay Kumar, Adopted in LVT
! 
!EOP
!BOP
    character*100              :: abort_message(20)
    
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

    ! EMK...Write messages to standard out before writing to abort_message file
    do i = 1,20
       write(LVT_logunit,6000,err=9000) trim(abort_message(i))
    end do

    message_file = 'abort_message'
    
    iofunc = 'opening'
    
    ftn = LVT_getNextUnitNumber()
    open (unit     = ftn, &
         file     = message_file,&
         err      = 9000, &
         iostat   = istat)
    iofunc = 'writing'
    
    do i = 1, 20
       write (ftn, 6000, err = 9000) trim(abort_message(i))
    enddo
    
    call LVT_releaseUnitNumber(ftn)
    
    call LVT_endrun
      
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
      
    call LVT_endrun
  
  end subroutine LVT_abort

!BOP
! 
! !ROUTINE: LVT_alert
! \label{LVT_alert}
!
! !INTERFACE:    
  subroutine LVT_alert( program_name, alert_number, message )
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
    character(len=*),  intent(in)     :: program_name  
    integer,           intent(in)     :: alert_number
    character(len=*), intent(inout)  :: message(20)
! 
! !OUTPUT PARAMETERS:
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
! !FILES USED:
!
! !REVISION HISTORY:
!     03 mar 1999 initial version .........................mr moore/dnxm
!     09 aug 1999 ported to IBM SP-2.  increase size of message to
!                 100 characters.  added intent attributes to 
!                 arguments................................mr gayno/dnxm     
!     29 oct 2005  Sujay Kumar, Adopted in LVT
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
    ftn = LVT_getNextUnitNumber()
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
    
    call LVT_releaseUnitNumber(ftn)
    
    return
  
!     ----------------------------------------------------------------
!     format statements
!     ----------------------------------------------------------------
    
8000 format (a)   
      
  end subroutine LVT_alert
 
!BOP
! 
! !ROUTINE: check_error
! \label(check_error)
!
! !INTERFACE:
  subroutine check_error(ierr,msg)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)          :: ierr
    character(len=*), intent(in) :: msg  
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This is an error check routine. Program exits in case of error
! with an associated error message written to the 'abort message'
! file. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  
    if ( ierr /= 0 ) then
       write(LVT_logunit,*) '[ERR] --------------------------------------------'
       write(LVT_logunit,*) '[ERR] ',msg
       write(LVT_logunit,*) '[ERR] program stopping ...'
       write(LVT_logunit,*) '[ERR] --------------------------------------------'
       call LVT_endrun
    endif

  end subroutine check_error

!BOP
! 
! !ROUTINE: verify
! \label{verify}
!
! !INTERFACE:
  subroutine verify(ierr)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)          :: ierr
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This is an error check routine. Program exits in case of error
! with an associated error message written to the 'abort message'
! file. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  
    if ( ierr /= 0) then
       write(LVT_logunit,*) '[ERR] --------------------------------------------'
       print*,'[ERR] return error in LVT call, Prgm Stopping...'
       write(LVT_logunit,*) '[ERR] --------------------------------------------'
       call LVT_endrun
    endif

  end subroutine verify

!BOP
! 
! !ROUTINE: LVT_getNextUnitNumber
! \label{LVT_getNextUnitNumber}
!
! !INTERFACE: 
    integer function LVT_getNextUnitNumber()
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! get next available Fortran unit number
!
! Method: 
! Get next available Fortran unit number itst. This routine is
! modified from the CCSM codes
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
      integer :: i 
! starting with 13 because 12 is reserved for grib
      do i = 13, stacksize
         if(.not.IOU(i)) then 
            LVT_getNextUnitNumber = i
            IOU(i) = .true.
            RETURN
         endif
      enddo
      write(*,*) 'LVT_fileUnitsMod: Ran out of Unit numbers'
    end function LVT_getNextUnitNumber
    
!BOP
! 
! !ROUTINE: LVT_releaseUnitNumber
! \label{LVT_releaseUnitNumber}
!
! !INTERFACE: 
    subroutine LVT_releaseUnitNumber(iunit)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
      integer, intent(in)    :: iunit
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! close and release Fortran unit no longer in use
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
      
      if(.not.IOU(iunit)) then 
         write(*,*) 'file unit release error: ',iunit,'is not flaggged as in use'
      endif
      close(iunit)
      IOU(iunit) = .false.
      return
    end subroutine LVT_releaseUnitNumber

!BOP
! 
! !ROUTINE: LVT_endrun
! \label{LVT_endrun}
!
! !INTERFACE:
  subroutine LVT_endrun
! 
! !USES: 
#if (defined SPMD)
    use LVT_mpiMod
#endif
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! Routine to be called to terminate the program. This routines 
! flushes the output streams and aborts the mpi processes.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  14 Nov 2002    Sujay Kumar  Initial Specification
! 
!EOP

    write(LVT_logunit,*)'[ERR] endrun is being called'
#if (defined SPMD) 
    call mpi_abort (MPI_COMM_WORLD, 1)  
#else
    error stop 1
#endif   
  end subroutine LVT_endrun

!BOP
! !ROUTINE: LVT_warning
! 
! !INTERFACE:
  subroutine LVT_warning(ierr,msg)
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
       write(LVT_logunit,*) '[WARN] ',msg
    endif

  end subroutine LVT_warning

  end module LVT_logMod
