!*******************************************************************************
!Subroutine - rapid_read_Qfor_file
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_read_Qfor_file

!Purpose:
!Read Qfor_file from Fortran.
!Author: 
!Cedric H. David, 2013-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petscvec.h>
use petscvec
use rapid_var, only :                                                          &
                   rank,ierr,ZV_read_for_tot,                                  &
                   ZV_Qfor,IS_for_bas,IV_for_loc2,IV_for_index,ZV_read_for_tot
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************


!*******************************************************************************
!Read file
!*******************************************************************************
if (rank==0) read(34,*) ZV_read_for_tot


!*******************************************************************************
!Set values in PETSc vector
!*******************************************************************************
if (rank==0) then
call VecSetValues(ZV_Qfor,IS_for_bas,IV_for_loc2,                              &
                  ZV_read_for_tot(IV_for_index),INSERT_VALUES,ierr)
                  !here we only look at the forcing within the basin studied 
end if

!*******************************************************************************
!Assemble PETSc vector
!*******************************************************************************
call VecAssemblyBegin(ZV_Qfor,ierr)
call VecAssemblyEnd(ZV_Qfor,ierr)


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_read_Qfor_file

#else

! Dummy version
subroutine rapid_read_Qfor_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_read_Qfor_file

#endif
