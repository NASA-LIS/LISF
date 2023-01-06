!*******************************************************************************
!Subroutine - rapid_close_Qhum_file 
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_close_Qhum_file

!Purpose:
!Close Qhum_file from Fortran.
!Author: 
!Cedric H. David, 2014-2020.


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
if (rank==0) close(36)

!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_close_Qhum_file

#else

! Dummy version
subroutine rapid_close_Qhum_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_close_Qhum_file

#endif
