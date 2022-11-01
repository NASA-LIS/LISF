!*******************************************************************************
!Subroutine - rapid_close_Qfor_file 
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_close_Qfor_file

!Purpose:
!Close Qfor_file from Fortran.
!Author: 
!Cedric H. David, 2013-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
use rapid_var, only :                                                          &
                   rank
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************


!*******************************************************************************
!Close file
!*******************************************************************************
if (rank==0) close(34)

!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_close_Qfor_file

#else

! Dummy version
subroutine rapid_close_Qfor_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_close_Qfor_file

#endif
