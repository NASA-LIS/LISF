!*******************************************************************************
!Subroutine - rapid_uq
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_uq

!Purpose:
!This subroutine generates river flow uncertainty estimates based on the bias,
!variance, and average covariances of the combined surface/subsurface flow into
!the river network.  
!Author: 
!Cedric H. David, 2016-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petscksp.h>
use petscksp
use rapid_var, only :                                                          &
                   IS_riv_tot,IS_riv_bas,                                      &
                   IV_riv_index,IV_riv_loc1,                                   &
                   ZS_TauR,                                                    &
                   ZV_bQlat,ZV_vQlat,ZV_caQlat,ZV_bQout,ZV_sQout,ZV_rQout,     &
                   ZM_Net,ZM_A,ZV_C1,ZV_one,                                   &
                   ZV_riv_tot_bQlat,ZV_riv_tot_vQlat,ZV_riv_tot_caQlat,        &
                   ZV_riv_bas_bQout,ZV_riv_bas_sQout,ZV_riv_bas_rQout,         &
                   ZV_nbuptot,ZV_one,ZV_temp1,ZV_temp2,                        &
                   ZV_SeqZero,ZV_pointer,ZS_one,                               &
                   ierr,rank,vecscat,ksp
implicit none


!*******************************************************************************
!Initialize variables
!*******************************************************************************
call VecSet(ZV_bQlat,0*ZS_one,ierr)
!Make sure the bias of lateral inflow is initialized to zero.

call VecSet(ZV_vQlat,0*ZS_one,ierr)
!Make sure the error variance of lateral inflow is initialized to zero.

call VecSet(ZV_caQlat,0*ZS_one,ierr)
!Make sure the average error covariance of lateral inflow is initialized to zero.

call VecSet(ZV_bQout,0*ZS_one,ierr)
!Make sure the bias of outflow is initialized to zero.

call VecSet(ZV_sQout,0*ZS_one,ierr)
!Make sure the standard error of outflow is initialized to zero.

call VecSet(ZV_rQout,0*ZS_one,ierr)
!Make sure the RMSE of outflow is initialized to zero.


!*******************************************************************************
!Modifying the linear system matrix for UQ and setting it in the solver
!*******************************************************************************
call MatCopy(ZM_Net,ZM_A,DIFFERENT_NONZERO_PATTERN,ierr)   !A=Net
call MatScale(ZM_A,-ZS_one,ierr)                           !A=-A
call MatShift(ZM_A,ZS_one,ierr)                            !A=A+1*I
!Result:A=I-Net

call KSPSetOperators(ksp,ZM_A,ZM_A,ierr)
!Set KSP to use matrix ZM_A


!*******************************************************************************
!Check that error variance provided is always positive
!*******************************************************************************
if (minval(ZV_riv_tot_vQlat) < 0) then
     print *, 'ERROR - The standard error provided includes negative values'
     stop 99
end if


!*******************************************************************************
!Apply the lateral inflow error estimates to the sub-basin of interest
!*******************************************************************************
if (rank==0) then
     call VecSetValues(ZV_bQlat,IS_riv_bas,IV_riv_loc1,                        &
                       ZV_riv_tot_bQlat(IV_riv_index),INSERT_VALUES,ierr)
end if
call VecAssemblyBegin(ZV_bQlat,ierr)
call VecAssemblyEnd(ZV_bQlat,ierr)

if (rank==0) then
     call VecSetValues(ZV_vQlat,IS_riv_bas,IV_riv_loc1,                        &
                       ZV_riv_tot_vQlat(IV_riv_index),INSERT_VALUES,ierr)
end if
call VecAssemblyBegin(ZV_vQlat,ierr)
call VecAssemblyEnd(ZV_vQlat,ierr)

if (rank==0) then
     call VecSetValues(ZV_caQlat,IS_riv_bas,IV_riv_loc1,                       &
                       ZV_riv_tot_caQlat(IV_riv_index),INSERT_VALUES,ierr)
end if
call VecAssemblyBegin(ZV_caQlat,ierr)
call VecAssemblyEnd(ZV_caQlat,ierr)


!*******************************************************************************
!Compute the bias of discharge
!*******************************************************************************
call KSPSolve(ksp,ZV_bQlat,ZV_bQout,ierr)
!solves A*bQout=bQlat


!*******************************************************************************
!Compute the standard error of discharge
!*******************************************************************************

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Compute the accumulation of error variances in lateral inflows
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call KSPSolve(ksp,ZV_vQlat,ZV_sQout,ierr)
!solves A*sQout=vQlat

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Compute the accumulation of error covariances in lateral inflows
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call KSPSolve(ksp,ZV_one,ZV_nbuptot,ierr)
!solves A*nbuptot=1

call VecShift(ZV_nbuptot,-ZS_one,ierr)
!nbuptot=nbuptot-1

call KSPSolve(ksp,ZV_caQlat,ZV_temp1,ierr)
!solves A*temp1=caQlat

call VecPointwiseMult(ZV_temp2,ZV_nbuptot,ZV_temp1,ierr)
!temp2=nbuptot.*temp1 (pointwise multiplication)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Compute the standard error in discharge from the variance in discharge
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecAXPY(ZV_sQout,ZS_one,ZV_temp2,ierr)
!sQout=sQout+temp2

call VecSqrtAbs(ZV_sQout,ierr)
!sQout=sqrt(sQout)


!*******************************************************************************
!Compute the RMSE of discharge
!*******************************************************************************
call VecCopy(ZV_bQout,ZV_rQout,ierr)
!rQout=bQout

call VecPow(ZV_rQout,ZS_one*2,ierr)
!rQout=rQout^2

call VecCopy(ZV_sQout,ZV_temp1,ierr)
!temp1=sQout

call VecPow(ZV_temp1,ZS_one*2,ierr)
!temp1=temp1^2

call VecAXPY(ZV_rQout,ZS_one,ZV_temp1,ierr)
!rQout=rQout+temp1

call VecSqrtAbs(ZV_rQout,ierr)
!rQout=sqrt(rQout)


!*******************************************************************************
!Gather PETSc vector on processor zero and get the values of the array
!*******************************************************************************

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Bias
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecScatterBegin(vecscat,ZV_bQout,ZV_SeqZero,                              &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat,ZV_bQout,ZV_SeqZero,                                &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)
!Gather PETSc vector

if (rank==0) call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)
!Get array from PETSc vector

if  (rank==0) ZV_riv_bas_bQout=ZV_pointer
!Copy values into the variable that will be written in netCDF file

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Standard error
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecScatterBegin(vecscat,ZV_sQout,ZV_SeqZero,                              &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat,ZV_sQout,ZV_SeqZero,                                &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)
!Gather PETSc vector

if (rank==0) call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)
!Get array from PETSc vector

if  (rank==0) ZV_riv_bas_sQout=ZV_pointer
!Copy values into the variable that will be written in netCDF file

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!RMSE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecScatterBegin(vecscat,ZV_rQout,ZV_SeqZero,                              &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat,ZV_rQout,ZV_SeqZero,                                &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)
!Gather PETSc vector

if (rank==0) call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)
!Get array from PETSc vector

if  (rank==0) ZV_riv_bas_rQout=ZV_pointer
!Copy values into the variable that will be written in netCDF file


!*******************************************************************************
!Recomputing the linear system matrix for RAPID and setting it in the solver
!*******************************************************************************
call MatCopy(ZM_Net,ZM_A,DIFFERENT_NONZERO_PATTERN,ierr)   !A=Net
call MatDiagonalScale(ZM_A,ZV_C1,ZV_one,ierr)              !A=diag(C1)*A
call MatScale(ZM_A,-ZS_one,ierr)                           !A=-A
call MatShift(ZM_A,ZS_one,ierr)                            !A=A+1*I
!Result:A=I-diag(C1)*Net

call KSPSetOperators(ksp,ZM_A,ZM_A,ierr)
!Set KSP to use matrix ZM_A. This is because this UQ subroutine had modified the
!matrix that was computed in rapid_init and that is later needed for routing. 


!*******************************************************************************
!End subroutine
!*******************************************************************************
call PetscPrintf(PETSC_COMM_WORLD,'Uncertainty quantification completed'       &
                 //char(10),ierr)
call PetscPrintf(PETSC_COMM_WORLD,'--------------------------'//char(10),ierr)

end subroutine rapid_uq

#else

! Dummy version
subroutine rapid_uq
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_uq

#endif
