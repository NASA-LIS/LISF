subroutine wrf_debug(unit_number, message)
  use LIS_logMod, only     : LIS_logunit
  implicit none
  integer :: unit_number 
  character(len=*) :: message 
  write(LIS_logunit, *) "Noah-MP.4.0.1 (wrf_debug): ", message 
end subroutine wrf_debug
