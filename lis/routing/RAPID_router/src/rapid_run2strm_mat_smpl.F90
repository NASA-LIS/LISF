!*******************************************************************************
!Subroutine - rapid_run2strm_mat_smpl
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_run2strm_mat_smpl

!Purpose:
!Compute the full-network observation operator - simplified/optimized algorithm
!This operator turns daily-averaged runoff into daily-averaged discharge
!This operator is used within the data assimilation as part of the
!observations operator H in the Kalman gain
!Authors: 
!Charlotte M. Emery, and Cedric H. David, 2018-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************

#include <petsc/finclude/petscmat.h>
use petscmat

use rapid_var, only :                                                          &
                IS_riv_bas,JS_riv_bas,JS_riv_bas2,JS_up,                       &
                IM_index_up,IV_riv_index,IV_nbup,                              &
                ZV_C3,ZV_k,                                                    &
                ierr,rank,                                                     &
                IS_one,ZS_val,                                                 &
                IV_nz,IV_dnz,IV_onz,                                           &
                IS_ownfirst,IS_ownlast,                                        &
                IS_R,IS_RpM,                                                   &
                ZM_L
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************

PetscInt :: IS_n
PetscInt, dimension(:), allocatable :: IV_index_down

PetscScalar :: ZS_val2, ZS_val3

VecScatter :: vecscat2
Vec :: ZV_C3_SeqAll, ZV_k_SeqAll

!*******************************************************************************
!Prepare for preallocation of ZM_L and associated temporary matrices
!*******************************************************************************

!-------------------------------------------------------------------------------
!Allocate and initialize temporary variables
!-------------------------------------------------------------------------------

IS_n = IS_RpM*IS_R

allocate(IV_index_down(IS_riv_bas))
IV_index_down(:) = 0
!Indexes of downstream reaches

IV_nz(:)=1
IV_dnz(:)=1
IV_onz(:)=0
!The number of non-zero elements per row.

call MatGetOwnershipRange(ZM_L,IS_ownfirst,IS_ownlast,ierr)

call VecScatterCreateToAll(ZV_k,vecscat2,ZV_C3_SeqAll,ierr)
call VecScatterCreateToAll(ZV_k,vecscat2,ZV_k_SeqAll,ierr)

!-------------------------------------------------------------------------------
!Populate temporary variables
!-------------------------------------------------------------------------------

do JS_riv_bas2=1,IS_riv_bas 
    do JS_up=1,IV_nbup(IV_riv_index(JS_riv_bas2))
        if (IM_index_up(JS_riv_bas2,JS_up)/=0) then

                JS_riv_bas=IM_index_up(JS_riv_bas2,JS_up)
                IV_index_down(JS_riv_bas)=JS_riv_bas2

        end if
    end do
end do
!IV_index_down gives, for each reach, the index of the unique downstream reach

!-------------------------------------------------------------------------------
!Create a temporary VecScatter object
!-------------------------------------------------------------------------------

call VecScatterBegin(vecscat2,ZV_C3,ZV_C3_SeqAll,                              &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat2,ZV_C3,ZV_C3_SeqAll,                                &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)

call VecScatterBegin(vecscat2,ZV_k,ZV_k_SeqAll,                                &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat2,ZV_k,ZV_k_SeqAll,                                  &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)

!-------------------------------------------------------------------------------
!Count number of non-zero elements of ZM_L
!-------------------------------------------------------------------------------

do JS_riv_bas=1,IS_riv_bas  !loop over columns

        ! diagonal term (already accounted for at initialization)
        call VecGetValues(ZV_C3_SeqAll,                                        &
                          IS_one,                                              &
                          JS_riv_bas-1,                                        &
                           ZS_val3,ierr)
        
        ZS_val = 1.0-(ZS_val3/real(IS_n))*                                     &
                     ((1.0-ZS_val3**(IS_n))/                                   &
                      (1.0-ZS_val3))

        ! next downstream reach
        JS_riv_bas2 = IV_index_down(JS_riv_bas)

        if (JS_riv_bas2 .ne. 0) then

                ! next downstream term
                call VecGetValues(ZV_k_SeqAll,                                 &
                                  IS_one,                                      &
                                  JS_riv_bas2-1,                               &
                                  ZS_val3,ierr)

                ZS_val2 = ZS_val - ZS_val3/86400.0

                do while ((JS_riv_bas2 .ne. 0).and.(ZS_val2 .ge. 0.0))

                        IV_nz(JS_riv_bas2) = IV_nz(JS_riv_bas2) + 1

                        ! diagonal block
                        if (((JS_riv_bas2.ge.IS_ownfirst+1).and.               &
                              (JS_riv_bas2.lt.IS_ownlast+1)).and.               &
                            ((JS_riv_bas.ge.IS_ownfirst+1).and.                &
                              (JS_riv_bas.lt.IS_ownlast+1))) then
   
                                IV_dnz(JS_riv_bas2) = IV_dnz(JS_riv_bas2)+1

                            end if

                        ! off-diagonal block
                        if (((JS_riv_bas2.ge.IS_ownfirst+1).and.               &
                              (JS_riv_bas2.lt.IS_ownlast+1)).and.               &
                            ((JS_riv_bas.lt.IS_ownfirst+1).or.                 &
                              (JS_riv_bas.ge.IS_ownlast+1))) then
   
                                IV_onz(JS_riv_bas2) = IV_onz(JS_riv_bas2)+1

                            end if

                        ZS_val = ZS_val2
                        JS_riv_bas2 = IV_index_down(JS_riv_bas2)
                        if (JS_riv_bas2 .ne. 0) then

                                call VecGetValues(ZV_k_SeqAll,                &
                                                    IS_one,                     &
                                                  JS_riv_bas2-1,              &
                                                  ZS_val3,ierr)

                                ZS_val2 = ZS_val - ZS_val3/86400.0

                        endif

                end do

        end if
end do

!*******************************************************************************
!Matrix preallocation (ZM_L)
!*******************************************************************************

call MatSeqAIJSetPreallocation(ZM_L,PETSC_DEFAULT_INTEGER,IV_nz,ierr)
call MatMPIAIJSetPreallocation(ZM_L,                                           &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_dnz(IS_ownfirst+1:IS_ownlast),               &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_onz(IS_ownfirst+1:IS_ownlast),ierr)

call PetscSynchronizedPrintf(PETSC_COMM_WORLD,'ZM_L pre-allocated'//char(10),ierr)


!*******************************************************************************
!Populate ZM_L
!*******************************************************************************

!-------------------------------------------------------------------------------
!Poupulate ZM_L
!-------------------------------------------------------------------------------

if (rank==0) then

do JS_riv_bas=1,IS_riv_bas  !loop over columns

        ! diagonal term (already accounted for at initialization)
        call VecGetValues(ZV_C3_SeqAll,                                        &
                          IS_one,                                              &
                          JS_riv_bas-1,                                        &
                           ZS_val3,ierr)

        ZS_val = 1.0-(ZS_val3/real(IS_n))*                                     &
                     ((1.0-ZS_val3**(IS_n))/                                   &
                      (1.0-ZS_val3))

        call MatSetValues(ZM_L,                                                &
                          IS_one,                                              &
                          JS_riv_bas-1,                                        &
                          IS_one,                                              &
                          JS_riv_bas-1,                                        &
                          ZS_val,                                              &
                          INSERT_VALUES,                                       &
                          ierr)

        ! next downstream reach
        JS_riv_bas2 = IV_index_down(JS_riv_bas)

        if (JS_riv_bas2 .ne. 0) then

                ! next downstream term
                call VecGetValues(ZV_k_SeqAll,                                 &
                                  IS_one,                                      &
                                  JS_riv_bas2-1,                               &
                                  ZS_val3,ierr)

                ZS_val2 = ZS_val - ZS_val3/86400.0

                do while ((JS_riv_bas2 .ne. 0).and.(ZS_val2 .ge. 0.0))

                        call MatSetValues(ZM_L,                                &
                                          IS_one,                              &
                                            JS_riv_bas2-1,                       &
                                            IS_one,                              &
                                            JS_riv_bas-1,                        &
                                            ZS_val2,                             &
                                            INSERT_VALUES,                       &
                                            ierr)

                        ZS_val = ZS_val2
                        JS_riv_bas2 = IV_index_down(JS_riv_bas2)
                        if (JS_riv_bas2 .ne. 0) then

                                call VecGetValues(ZV_k_SeqAll,                 &
                                                  IS_one,                      &
                                                  JS_riv_bas2-1,               &
                                                  ZS_val3,ierr)

                                ZS_val2 = ZS_val - ZS_val3/86400.0

                        end if

                end do
        end if

end do

end if

!-------------------------------------------------------------------------------
!Matrices assembly
!-------------------------------------------------------------------------------

call MatAssemblyBegin(ZM_L,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(ZM_L,MAT_FINAL_ASSEMBLY,ierr)

call PetscPrintf(PETSC_COMM_WORLD,'Inverse runoff-discharge operator created'    &
                                   //char(10),ierr)

!*******************************************************************************
!Free up memory used by local (temporary) variables
!*******************************************************************************

deallocate(IV_index_down)

call VecDestroy(ZV_C3_SeqAll,ierr)
call VecDestroy(ZV_k_SeqAll,ierr)
call VecScatterDestroy(vecscat2,ierr)

!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_run2strm_mat_smpl

#else

! Dummy version
subroutine rapid_run2strm_mat_smpl
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_run2strm_mat_smpl

#endif
