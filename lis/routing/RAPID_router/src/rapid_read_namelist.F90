!*******************************************************************************
!Subroutine - rapid_read_namelist
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_read_namelist

!Purpose:
!This subroutine allows to read the RAPID namelist and hence to run the model
!multiple times without ever have to recompile.  Some information on the options
!used is also printed in the stdout.
!Author: 
!Cedric H. David, 2011-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
use rapid_var, only :                                                          &
                     NL_namelist,namelist_file
implicit none


!*******************************************************************************
!Read namelist file 
!*******************************************************************************
open(88,file=namelist_file,status='old',form='formatted')
read(88, NL_namelist)
close(88)


!*******************************************************************************
!Optional prints what was read 
!*******************************************************************************
!print *, namelist_file


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_read_namelist

#else

! Dummy version
subroutine rapid_read_namelist
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_read_namelist

#endif
