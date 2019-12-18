!!! Shugong Wang, 07/12/2019 LIS version of logging_mod with write_to_log disabled
!!! 
#if !defined(UM_JULES)

MODULE logging_mod
USE mpi
USE iso_fortran_env, ONLY: output_unit, error_unit
IMPLICIT NONE

INTEGER, PARAMETER ::                                                         &
  log_level_info  = 1,                                                        &
  log_level_debug = 2,                                                        &
  log_level_warn  = 4,                                                        &
  log_level_error = 8,                                                        &
  log_level_fatal = 16
INTEGER, PARAMETER :: log_print_level = 31
INTEGER, PARAMETER :: log_stop_level = 0
INTEGER, PARAMETER :: log_max_line_len = 500
INTEGER :: ntasks = 0
CHARACTER(LEN=20) :: task_prefix = ""
PRIVATE
PUBLIC :: log_info, log_debug, log_warn, log_error, log_fatal

CONTAINS
SUBROUTINE write_to_log(log_level, proc_name, message)
  USE string_utils_mod, ONLY: to_string
  INTEGER, INTENT(IN) :: log_level
  CHARACTER(LEN=*), INTENT(IN) :: proc_name
  CHARACTER(LEN=*), INTENT(IN) :: message
  CHARACTER(LEN=log_max_line_len) :: full_message
  INTEGER :: task_id, error
  full_message = TRIM(proc_name) // ': ' // message
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
    CALL log_fatal("write_to_log", "Unknown log level")
  END SELECT
END SUBROUTINE write_to_log

SUBROUTINE log_info(proc_name, message)
  CHARACTER(LEN=*), INTENT(IN) :: proc_name
  CHARACTER(LEN=*), INTENT(IN) :: message
  CALL write_to_log(log_level_info, proc_name, message)
END SUBROUTINE log_info
SUBROUTINE log_debug(proc_name, message)
  CHARACTER(LEN=*), INTENT(IN) :: proc_name
  CHARACTER(LEN=*), INTENT(IN) :: message
  CALL write_to_log(log_level_debug, proc_name, message)
END SUBROUTINE log_debug
SUBROUTINE log_warn(proc_name, message)
  CHARACTER(LEN=*), INTENT(IN) :: proc_name
  CHARACTER(LEN=*), INTENT(IN) :: message
  CALL write_to_log(log_level_warn, proc_name, message)
END SUBROUTINE log_warn
SUBROUTINE log_error(proc_name, message)
  CHARACTER(LEN=*), INTENT(IN) :: proc_name
  CHARACTER(LEN=*), INTENT(IN) :: message
  CALL write_to_log(log_level_error, proc_name, message)
END SUBROUTINE log_error
SUBROUTINE log_fatal(proc_name, message)
  CHARACTER(LEN=*), INTENT(IN) :: proc_name
  CHARACTER(LEN=*), INTENT(IN) :: message
  CALL write_to_log(log_level_fatal, proc_name, message)
END SUBROUTINE log_fatal
END MODULE logging_mod
#endif
