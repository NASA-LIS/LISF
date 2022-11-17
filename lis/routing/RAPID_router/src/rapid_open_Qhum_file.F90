!*******************************************************************************
!Subroutine - rapid_open_Qhum
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_open_Qhum_file(Qhum_file) 

!Purpose:
!Open Qhum_file from Fortran.
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
character(len=*), intent(in):: Qhum_file


!*******************************************************************************
!Open file
!*******************************************************************************
if (rank==0) open(36,file=Qhum_file,status='old')


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_open_Qhum_file

#else

! Dummy version
subroutine rapid_open_Qhum_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_open_Qhum_file

#endif
