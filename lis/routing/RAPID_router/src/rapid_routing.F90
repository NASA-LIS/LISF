!*******************************************************************************
!Subroutine - rapid_routing
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_routing(ZV_C1,ZV_C2,ZV_C3,ZV_Qext,                            &
                         ZV_QoutinitR,                                         &
                         ZV_QoutR,ZV_QoutbarR)

!Purpose:
!Performs flow calculation in each reach of a river network using the Muskingum
!method (McCarthy 1938).  
!Author: 
!Cedric H. David, 2008-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petscksp.h>
use petscksp
use netcdf
use rapid_var, only :                                                          &
                   ZS_dtR,IS_R,JS_R,                                           &
                   ZM_Net,ZM_TC1,ZM_M,                                         &
                   ZV_b,ZV_babsmax,ZV_bhat,                                    &
                   ZV_QoutprevR,ZV_QoutRabsmin,ZV_QoutRabsmax,                 &
                   ZV_QoutRhat,                                                &
                   ierr,ksp,                                                   &
                   ZS_one,IS_ksp_iter,IS_ksp_iter_max,                         &
                   IS_riv_bas,JS_riv_bas,IM_index_up,                          &
                   IS_opt_routing,IV_nbup,IV_riv_index,                        &
                   BS_opt_influence
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************
Vec, intent(in)    :: ZV_C1,ZV_C2,ZV_C3,ZV_Qext,                               &
                      ZV_QoutinitR 
Vec, intent(inout) :: ZV_QoutR,ZV_QoutbarR

PetscInt :: IS_localsize,JS_localsize
PetscScalar, pointer :: ZV_QoutR_p(:),ZV_QoutprevR_p(:),                       &
                        ZV_Qext_p(:),ZV_C1_p(:),ZV_C2_p(:),                    &
                        ZV_C3_p(:),ZV_b_p(:),                                  &
                        ZV_babsmax_p(:),ZV_QoutRabsmin_p(:),ZV_QoutRabsmax_p(:)


!*******************************************************************************
!Get local sizes for vectors
!*******************************************************************************
call VecGetLocalSize(ZV_QoutR,IS_localsize,ierr)


!*******************************************************************************
!Set mean values to zero initialize QoutprevR with QoutinitR
!*******************************************************************************
call VecSet(ZV_QoutbarR,0*ZS_one,ierr)                     !Qoutbar=0 
!set the means to zero at beginning of iterations over routing time step

call VecCopy(ZV_QoutinitR,ZV_QoutprevR,ierr)               !QoutprevR=QoutinitR
!set the previous value to the initial value given as input to subroutine


!*******************************************************************************
!Temporal loop 
!*******************************************************************************
call VecGetArrayF90(ZV_C1,ZV_C1_p,ierr)
call VecGetArrayF90(ZV_C2,ZV_C2_p,ierr)
call VecGetArrayF90(ZV_C3,ZV_C3_p,ierr)
call VecGetArrayF90(ZV_Qext,ZV_Qext_p,ierr)

do JS_R=1,IS_R
!-------------------------------------------------------------------------------
!Update mean
!-------------------------------------------------------------------------------
call VecAXPY(ZV_QoutbarR,ZS_one/IS_R,ZV_QoutprevR,ierr) 
!Qoutbar=Qoutbar+Qoutprev/IS_R

!-------------------------------------------------------------------------------
!Calculation of the right hand size, b
!-------------------------------------------------------------------------------
call MatMult(ZM_Net,ZV_QoutprevR,ZV_b,ierr)                !b2=Net*Qoutprev

call VecGetArrayF90(ZV_b,ZV_b_p,ierr)
call VecGetArrayF90(ZV_QoutprevR,ZV_QoutprevR_p,ierr)

do JS_localsize=1,IS_localsize
     ZV_b_p(JS_localsize)=ZV_b_p(JS_localsize)*ZV_C2_p(JS_localsize)           &
                         +(ZV_C1_p(JS_localsize)+ZV_C2_p(JS_localsize))        &
                         *ZV_Qext_p(JS_localsize)                              &
                         +ZV_C3_p(JS_localsize)*ZV_QoutprevR_p(JS_localsize)
end do

call VecRestoreArrayF90(ZV_QoutprevR,ZV_QoutprevR_p,ierr)
call VecRestoreArrayF90(ZV_b,ZV_b_p,ierr)

!-------------------------------------------------------------------------------
!Routing with PETSc using a matrix method
!-------------------------------------------------------------------------------
if (IS_opt_routing==1) then

call KSPSolve(ksp,ZV_b,ZV_QoutR,ierr)                      !solves A*Qout=b
call KSPGetIterationNumber(ksp,IS_ksp_iter,ierr)
if (IS_ksp_iter>IS_ksp_iter_max) IS_ksp_iter_max=IS_ksp_iter

end if

!-------------------------------------------------------------------------------
!Routing with Fortran using the traditional Muskingum method
!-------------------------------------------------------------------------------
if (IS_opt_routing==2) then

call VecGetArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
call VecGetArrayF90(ZV_QoutprevR,ZV_QoutprevR_p,ierr)
call VecGetArrayF90(ZV_b,ZV_b_p,ierr)

do JS_riv_bas=1,IS_riv_bas
     ZV_QoutR_p(JS_riv_bas)=ZV_b_p(JS_riv_bas)                                 &
                            +sum(ZV_C1_p(JS_riv_bas)                           &
                                  *ZV_QoutR_p(IM_index_up(JS_riv_bas,1:        &
                                   IV_nbup(IV_riv_index(JS_riv_bas))))) 
end do
!Taking into account the knowledge of how many upstream locations exist.
!Similar to exact preallocation of network matrix

call VecRestoreArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
call VecRestoreArrayF90(ZV_QoutprevR,ZV_QoutprevR_p,ierr)
call VecRestoreArrayF90(ZV_b,ZV_b_p,ierr)
end if

!-------------------------------------------------------------------------------
!Routing with PETSc using a matrix method with transboundary matrix
!-------------------------------------------------------------------------------
if (IS_opt_routing==3) then

call KSPSetType(ksp,KSPPREONLY,ierr)                         !preonly enforced

call KSPSolve(ksp,ZV_b,ZV_QoutRhat,ierr)                     !solves A*Qouthat=b
call KSPGetIterationNumber(ksp,IS_ksp_iter,ierr)
if (IS_ksp_iter>IS_ksp_iter_max) IS_ksp_iter_max=IS_ksp_iter

call MatMult(ZM_TC1,ZV_QoutRhat,ZV_bhat,ierr)
call VecAYPX(ZV_bhat,ZS_one,ZV_b,ierr)

call KSPSolve(ksp,ZV_bhat,ZV_QoutR,ierr)                     !solves A*Qout=bhat
call KSPGetIterationNumber(ksp,IS_ksp_iter,ierr)
if (IS_ksp_iter>IS_ksp_iter_max) IS_ksp_iter_max=IS_ksp_iter

end if

!-------------------------------------------------------------------------------
!Routing with PETSc using a matrix method with Muskingum operator
!-------------------------------------------------------------------------------
if (IS_opt_routing==4) then

call MatMult(ZM_M,ZV_b,ZV_QoutR,ierr)

end if

!-------------------------------------------------------------------------------
!Calculation of babsmax, QoutRabsmin and QoutRabsmax
!-------------------------------------------------------------------------------
if (BS_opt_influence) then

call VecGetArrayF90(ZV_b,ZV_b_p,ierr)
call VecGetArrayF90(ZV_babsmax,ZV_babsmax_p,ierr)
do JS_localsize=1,IS_localsize
     if (ZV_babsmax_p(JS_localsize)<=abs(ZV_b_p(JS_localsize))) then
         ZV_babsmax_p(JS_localsize) =abs(ZV_b_p(JS_localsize))
     end if
end do
call VecRestoreArrayF90(ZV_b,ZV_b_p,ierr)
call VecRestoreArrayF90(ZV_babsmax,ZV_babsmax_p,ierr)

call VecGetArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
call VecGetArrayF90(ZV_QoutRabsmin,ZV_QoutRabsmin_p,ierr)
call VecGetArrayF90(ZV_QoutRabsmax,ZV_QoutRabsmax_p,ierr)
do JS_localsize=1,IS_localsize
     if (ZV_QoutRabsmin_p(JS_localsize)>=abs(ZV_QoutR_p(JS_localsize))) then
         ZV_QoutRabsmin_p(JS_localsize) =abs(ZV_QoutR_p(JS_localsize))
     end if
     if (ZV_QoutRabsmax_p(JS_localsize)<=abs(ZV_QoutR_p(JS_localsize))) then
         ZV_QoutRabsmax_p(JS_localsize) =abs(ZV_QoutR_p(JS_localsize))
     end if
end do
call VecRestoreArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
call VecRestoreArrayF90(ZV_QoutRabsmin,ZV_QoutRabsmin_p,ierr)
call VecRestoreArrayF90(ZV_QoutRabsmax,ZV_QoutRabsmax_p,ierr)

end if

!-------------------------------------------------------------------------------
!Reset previous
!-------------------------------------------------------------------------------
call VecCopy(ZV_QoutR,ZV_QoutprevR,ierr)              !Qoutprev=Qout
!reset previous 


!-------------------------------------------------------------------------------
!optional write outputs
!-------------------------------------------------------------------------------
!call VecScatterBegin(vecscat,ZV_QoutR,ZV_SeqZero,                              &
!                     INSERT_VALUES,SCATTER_FORWARD,ierr)
!call VecScatterEnd(vecscat,ZV_QoutR,ZV_SeqZero,                                &
!                        INSERT_VALUES,SCATTER_FORWARD,ierr)
!call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)
!!if (rank==0) write (99,'(10e10.3)') ZV_pointer
!if (rank==0) IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qout,IS_nc_id_var_Qout,    &
!                                       ZV_pointer,                             &
!                     [IV_nc_start(1),(IV_nc_start(2)-1)*IS_R+JS_R],IV_nc_count2)
!call VecRestoreArrayF90(ZV_SeqZero,ZV_pointer,ierr)


!-------------------------------------------------------------------------------
!End temporal loop
!-------------------------------------------------------------------------------
end do

call VecRestoreArrayF90(ZV_C1,ZV_C1_p,ierr)
call VecRestoreArrayF90(ZV_C2,ZV_C2_p,ierr)
call VecRestoreArrayF90(ZV_C3,ZV_C3_p,ierr)
call VecRestoreArrayF90(ZV_Qext,ZV_Qext_p,ierr)


!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_routing

#else

! Dummy version
subroutine rapid_routing
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_routing

#endif
