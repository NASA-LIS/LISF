!*******************************************************************************
!Subroutine - rapid_kf_obs_mat
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_kf_obs_mat

!Purpose:
!Using the set of available observation gage,
!Compute the observation operator for data assimilation (ZM_H)
!This operator is a submatrix of runoff-to-discharge operator (ZM_L)
!As it contains only the rows of ZM_L that corresponds to potentially 
!observed reach
!Authors: 
!Charlotte M. Emery, and Cedric H. David, 2018-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************

#include <petsc/finclude/petscmat.h>
use petscmat

use rapid_var, only :                                                          &
                IS_riv_bas,                                                    &
                IV_obs_loc1,IS_obs_bas,JS_obs_bas,                             &
                ierr,rank,                                                     &
                IS_one,ZS_one,                                                 &
                ZM_L,ZM_S,ZM_H
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************

PetscInt :: IS_val
PetscInt :: ISobs_ownfirst, ISobs_ownlast

PetscInt, dimension(:), allocatable :: IV_sorted_obs_loc1, IV_sort
PetscInt, dimension(:), allocatable :: IVobs_nz, IVobs_dnz, IVobs_onz

!*******************************************************************************
!Prepare for selection operator preallocation
!*******************************************************************************

!-------------------------------------------------------------------------------
!Allocate/Initialize local (temporary variables)
!-------------------------------------------------------------------------------

allocate(IV_sorted_obs_loc1(IS_obs_bas))
allocate(IV_sort(IS_obs_bas))

allocate(IVobs_nz(IS_obs_bas))
allocate(IVobs_dnz(IS_obs_bas))
allocate(IVobs_onz(IS_obs_bas))

IVobs_nz(:) = 1
IVobs_dnz(:) = 1
IVobs_onz(:) = 0

!-------------------------------------------------------------------------------
!Sort observation gages from upstream to downstream
!-------------------------------------------------------------------------------

do JS_obs_bas = 1,IS_obs_bas
    IV_sort(JS_obs_bas) = IV_obs_loc1(JS_obs_bas)
end do

do JS_obs_bas = 1,IS_obs_bas
    IS_val = MINLOC(IV_sort,DIM=1)
    IV_sorted_obs_loc1(JS_obs_bas) = IV_sort(IS_val)
    IV_sort(IS_val) = MAXVAL(IV_sort)+1   
end do

!-------------------------------------------------------------------------------
!Create observation operator matrix
!-------------------------------------------------------------------------------

call MatCreate(PETSC_COMM_WORLD,ZM_S,ierr)
call MatSetSizes(ZM_S,PETSC_DECIDE,PETSC_DECIDE,IS_obs_bas,IS_riv_bas,ierr)
call MatSetFromOptions(ZM_S,ierr)
call MatSetUp(ZM_S,ierr)

!-------------------------------------------------------------------------------
!Count non-zeros in ZM_S
!-------------------------------------------------------------------------------

call MatGetOwnershipRange(ZM_S,ISobs_ownfirst,ISobs_ownlast,ierr) 

do JS_obs_bas = 1,IS_obs_bas

    IS_val = IV_sorted_obs_loc1(JS_obs_bas)

    if (((JS_obs_bas.ge.ISobs_ownfirst+1).and.(JS_obs_bas.lt.ISobs_ownlast+1)).and.   &
        ((IS_val.ge.ISobs_ownfirst).and.(IS_val.lt.ISobs_ownlast))) then
        IVobs_dnz(JS_obs_bas) = 1
    end if

    if (((JS_obs_bas.ge.ISobs_ownfirst+1).and.(JS_obs_bas.lt.ISobs_ownlast+1)).and.   &
        ((IS_val.lt.ISobs_ownfirst).or.(IS_val.ge.ISobs_ownlast))) then
        IVobs_onz(JS_obs_bas) = 1
    end if
end do


!*******************************************************************************
!Selection operator preallocation
!*******************************************************************************

call MatSeqAIJSetPreallocation(ZM_S,PETSC_DEFAULT_INTEGER,IVobs_nz,ierr)
call MatMPIAIJSetPreallocation(ZM_S,                                           &
                               PETSC_DEFAULT_INTEGER,                          &
                               IVobs_dnz(ISobs_ownfirst+1:ISobs_ownlast),      &
                               PETSC_DEFAULT_INTEGER,                          &
                               IVobs_onz(ISobs_ownfirst+1:ISobs_ownlast),ierr)

!*******************************************************************************
!Compute operator preallocation
!*******************************************************************************

if (rank.eq.0) then

do JS_obs_bas = 1,IS_obs_bas
    IS_val = IV_sorted_obs_loc1(JS_obs_bas)

    call MatSetValues(ZM_S,                         &
                      IS_one,                       &
                      JS_obs_bas-1,                 &
                      IS_one,                       &
                      IS_val,                       &
                      ZS_one,                       &
                      INSERT_VALUES,ierr)
                      
end do

end if

call MatAssemblyBegin(ZM_S,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(ZM_S,MAT_FINAL_ASSEMBLY,ierr)

!*******************************************************************************
!Compute observation operator ZM_H = ZM_S*ZM_L
!*******************************************************************************

call MatMatMult(ZM_S,ZM_L,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,ZM_H,ierr)

call PetscPrintf(PETSC_COMM_WORLD,'Observation operator ZM_H created'//char(10),ierr)

!*******************************************************************************
!Free up memory used by local (temporary) variables
!*******************************************************************************

deallocate(IV_sorted_obs_loc1)
deallocate(IV_sort)

call MatDestroy(ZM_L,ierr)

!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_kf_obs_mat

#else

! Dummy version
subroutine rapid_kf_obs_mat
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_kf_obs_mat

#endif
