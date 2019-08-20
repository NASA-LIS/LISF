!!! Shugong Wang, 07/12/2019 LIS version of logging_mod with write_to_log disabled
#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE logging_mod

USE mpi

USE iso_fortran_env, ONLY: output_unit, error_unit

IMPLICIT NONE

! Log levels - these can be combined using bitwise operators to indicate
! any combination of log levels
INTEGER, PARAMETER ::                                                         &
  log_level_info  = 1,                                                        &
  log_level_debug = 2,                                                        &
  log_level_warn  = 4,                                                        &
  log_level_error = 8,                                                        &
  log_level_fatal = 16

! Determines what log levels are printed to log_unit - this is a bitwise
! combination of values from above.
! The default (31) is to print everything
INTEGER, PARAMETER :: log_print_level = 31

! Determines what log levels cause a program to stop - this is a bitwise
! combination of values from above.
! Fatal errors will always cause the program to stop, by definition.
! The default (0) is to stop only for fatal errors
! Setting this to 15 (i.e. stop for everything, even info) or 14 (i.e. stop
! for everything except info) are useful options for debugging
INTEGER, PARAMETER :: log_stop_level = 0

! The maximum line length for log entries. Any message larger than this is
! truncated
INTEGER, PARAMETER :: log_max_line_len = 500

! The calculated number of tasks id and task prefix
! These are cached the first time they are calculated to avoid calling the
! MPI routines every time a log message is written
INTEGER :: ntasks = 0
CHARACTER(LEN=20) :: task_prefix = ""


! Visibilities
PRIVATE
PUBLIC :: log_info, log_debug, log_warn, log_error, log_fatal

CONTAINS

SUBROUTINE write_to_log(log_level, proc_name, message)

#if defined(INTEL_FORTRAN)
USE ifcore, ONLY: tracebackqq
#endif

USE string_utils_mod, ONLY: to_string

!-----------------------------------------------------------------------------
! Description:
!   Logs the given message at the given level and performs any necessary
!   action
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: log_level
  ! The level at which to log the given message
CHARACTER(LEN=*), INTENT(IN) :: proc_name
  ! The name of the originating routine/function
CHARACTER(LEN=*), INTENT(IN) :: message
  ! The message to log

! Work variables
CHARACTER(LEN=log_max_line_len) :: full_message
INTEGER :: task_id, error

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Construct the full log message
!-----------------------------------------------------------------------------
! Combine the name of the originating procedure with the message
full_message = TRIM(proc_name) // ': ' // message

! Prepend a fragment to the message depending on the log level
SELECT CASE ( log_level )
CASE ( log_level_info )
  full_message = "[INFO] " // full_message

CASE ( log_level_debug )
  full_message = "[DEBUG] " // full_message

CASE ( log_level_warn )
  full_message = "[WARNING] " // full_message

CASE ( log_level_error )
  full_message = "[ERROR] " // full_message

CASE ( log_level_fatal )
  full_message = "[FATAL ERROR] " // full_message

CASE DEFAULT
  ! This should never happen since the only public access to write_to_log is
  ! through the log_* routines defined below
  CALL log_fatal("write_to_log", "Unknown log level")
END SELECT

!!! disabled by Shugong Wang 
#if 0
! Calculate the number of tasks and task prefix if we have not already
IF ( ntasks <= 0 ) THEN
  CALL mpi_comm_size(mpi_comm_world, ntasks, error)
  CALL mpi_comm_rank(mpi_comm_world, task_id, error)
  IF ( ntasks > 1 ) THEN
    task_prefix = "{MPI Task " // TRIM(to_string(task_id)) // "}"
  END IF
END IF

! Prepend the task prefix
full_message = ADJUSTL(TRIM(task_prefix) // " " // full_message)

!-----------------------------------------------------------------------------
! Use a bitwise and to check if we want to print log messages for the
! given log level
!-----------------------------------------------------------------------------
IF ( IAND(log_print_level, log_level) > 0 ) THEN
  IF ( log_level <= log_level_debug ) THEN
    ! Print info and debug to stdout
    WRITE(output_unit, "(A)") TRIM(full_message)
  ELSE
    ! Print errors to stdout
    WRITE(error_unit, "(A)") TRIM(full_message)
  END IF
END IF

!-----------------------------------------------------------------------------
! Check if we need to stop the program
!-----------------------------------------------------------------------------
IF ( IAND(log_stop_level, log_level) > 0 .OR.                                 &
     log_level == log_level_fatal ) THEN
  ! If stopping abnormally, attempt to print a stack trace
#if defined(GNU_FORTRAN)
  ! Currently, the default Met Office version of gfortran is too old to support
  ! BACKTRACE, as is the default gfortran version in some common OS package
  ! systems
  ! We'll add this back in when gfortran 4.8 is fully available
  !      CALL BACKTRACE()
#elif defined(INTEL_FORTRAN)
  CALL tracebackqq(user_exit_code = -1)
#endif

  ! Abort MPI with a non-zero exit code
  CALL mpi_abort(mpi_comm_world, 1, error)
END IF
#endif
END SUBROUTINE write_to_log


SUBROUTINE log_info(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level info
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: proc_name
CHARACTER(LEN=*), INTENT(IN) :: message

CALL write_to_log(log_level_info, proc_name, message)

END SUBROUTINE log_info


SUBROUTINE log_debug(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level debug
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: proc_name
CHARACTER(LEN=*), INTENT(IN) :: message

CALL write_to_log(log_level_debug, proc_name, message)

END SUBROUTINE log_debug


SUBROUTINE log_warn(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level warn
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: proc_name
CHARACTER(LEN=*), INTENT(IN) :: message

CALL write_to_log(log_level_warn, proc_name, message)

END SUBROUTINE log_warn


SUBROUTINE log_error(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level error
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: proc_name
CHARACTER(LEN=*), INTENT(IN) :: message

CALL write_to_log(log_level_error, proc_name, message)

END SUBROUTINE log_error


SUBROUTINE log_fatal(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level fatal
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: proc_name
CHARACTER(LEN=*), INTENT(IN) :: message

CALL write_to_log(log_level_fatal, proc_name, message)

END SUBROUTINE log_fatal

END MODULE logging_mod
#endif
