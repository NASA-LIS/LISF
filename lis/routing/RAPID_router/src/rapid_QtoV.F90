!*******************************************************************************
!Subroutine - rapid_QtoV
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_QtoV(ZV_k,ZV_x,ZV_QoutbarR,ZV_Qext,ZV_VbarR) 

!Purpose:
!Computes the volume of water in each river reach from the flows based on the 
!Muskingum method (McCarthy 1938). 
!Author: 
!Cedric H. David, 2015-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************

#include <petsc/finclude/petscmat.h>
use petscmat

use rapid_var, only :                                                          &
                   ZM_Net,ierr
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************
Vec, intent(in)    :: ZV_QoutbarR,ZV_Qext,                                     &
                      ZV_k,ZV_x 
Vec, intent(inout) :: ZV_VbarR

PetscInt :: IS_localsize,JS_localsize

PetscScalar, pointer :: ZV_QoutbarR_p(:),ZV_Qext_p(:),                         &
                        ZV_k_p(:),ZV_x_p(:),                                   &
                        ZV_VbarR_p(:)


!*******************************************************************************
!Get local sizes for vectors
!*******************************************************************************
call VecGetLocalSize(ZV_QoutbarR,IS_localsize,ierr)


!*******************************************************************************
!Calculation of Volume 
!*******************************************************************************
call MatMult(ZM_Net,ZV_QoutbarR,ZV_VbarR,ierr)                
!Vbar=Net*QoutbarR

call VecGetArrayF90(ZV_QoutbarR,ZV_QoutbarR_p,ierr)
call VecGetArrayF90(ZV_Qext,ZV_Qext_p,ierr)
call VecGetArrayF90(ZV_k,ZV_k_p,ierr)
call VecGetArrayF90(ZV_x,ZV_x_p,ierr)
call VecGetArrayF90(ZV_VbarR,ZV_VbarR_p,ierr)

do JS_localsize=1,IS_localsize
     ZV_VbarR_p(JS_localsize)=ZV_VbarR_p(JS_localsize)+ZV_Qext_p(JS_localsize) 
     !VbarR=VbarR+Qext
     ZV_VbarR_p(JS_localsize)=ZV_VbarR_p(JS_localsize)*ZV_x_p(JS_localsize)
     !VbarR=VbarR*x
     ZV_VbarR_p(JS_localsize)=ZV_VbarR_p(JS_localsize)                         &
                             +(1-ZV_x_p(JS_localsize))                         &
                             *ZV_QoutbarR_p(JS_localsize) 
     !VbarR=VbarR+(1-x)*QoutbarR
     ZV_VbarR_p(JS_localsize)=ZV_VbarR_p(JS_localsize)*ZV_k_p(JS_localsize)
     !VbarR=VbarR*x
end do

call VecRestoreArrayF90(ZV_QoutbarR,ZV_QoutbarR_p,ierr)
call VecRestoreArrayF90(ZV_Qext,ZV_Qext_p,ierr)
call VecRestoreArrayF90(ZV_k,ZV_k_p,ierr)
call VecRestoreArrayF90(ZV_x,ZV_x_p,ierr)
call VecRestoreArrayF90(ZV_VbarR,ZV_VbarR_p,ierr)


!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_QtoV

#else

! Dummy version
subroutine rapid_QtoV
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_QtoV

#endif
