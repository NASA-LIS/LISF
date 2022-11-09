!*******************************************************************************
!Subroutine - rapid_read_Qhum_file
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_read_Qhum_file

!Purpose:
!Read Qhum_file from Fortran.
!Author: 
!Cedric H. David, 2014-2020.


!!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petscvec.h>
use petscvec
use rapid_var, only :                                                          &
                   rank,ierr,ZV_read_hum_tot,                                  &
                   ZV_Qhum,IS_hum_bas,IV_hum_loc1,IV_hum_index,ZV_read_hum_tot
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************


!*******************************************************************************
!Read file
!*******************************************************************************
if (rank==0) read(36,*) ZV_read_hum_tot


!*******************************************************************************
!Set values in PETSc vector
!*******************************************************************************
if (rank==0) then
call VecSetValues(ZV_Qhum,IS_hum_bas,IV_hum_loc1,                              &
                  ZV_read_hum_tot(IV_hum_index),INSERT_VALUES,ierr)
                  !here we only look at the human-induced flows within the basin 
                  !studied 
end if

!*******************************************************************************
!Assemble PETSc vector
!*******************************************************************************
call VecAssemblyBegin(ZV_Qhum,ierr)
call VecAssemblyEnd(ZV_Qhum,ierr)


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_read_Qhum_file

#else

! Dummy version
subroutine rapid_read_Qhum_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_read_Qhum_file

#endif
