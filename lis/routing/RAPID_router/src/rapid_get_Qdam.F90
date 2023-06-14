!*******************************************************************************
!Subroutine - rapid_get_Qdam
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_get_Qdam

!Purpose:
!Communicate with a dam subroutine to exchange inflows and outflows.
!Authors:
!Cedric H. David and Ahmad A. Tavakoly, 2013-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************

#include <petsc/finclude/petscmat.h>
use petscmat

use rapid_var, only :                                                          &
                   rank,ierr,vecscat,ZV_pointer,ZV_SeqZero,ZS_one,             &
                   ZM_Net,ZV_Qext,ZV_Qdam,ZV_QoutbarR,ZV_QinbarR,              &
                   IS_dam_bas,IV_dam_index,IV_dam_loc2,IV_dam_pos,             &
                   JS_dam_tot,IS_dam_tot,                                      &
                   ZV_Qout_dam,ZV_Qin_dam_prev,ZV_Qout_dam_prev,               &
                   ZV_k_dam,ZV_p_dam,ZV_S_dam,ZV_Smax_dam,ZV_Smin_dam,         &
                   ZS_TauR
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************


!*******************************************************************************
!Compute previous inflow from river network and outside of river network to dams
!*******************************************************************************
!-------------------------------------------------------------------------------
!Compute inflow into dams from previous river flow
!-------------------------------------------------------------------------------
call MatMult(ZM_Net,ZV_QoutbarR,ZV_QinbarR,ierr)           
call VecAXPY(ZV_QinbarR,ZS_one,ZV_Qext,ierr)
!QinbarR=Net*QoutbarR+Qext

!-------------------------------------------------------------------------------
!Set values from PETSc vector into Fortran vector 
!-------------------------------------------------------------------------------
if (rank==0) ZV_Qin_dam_prev=0 
call VecScatterBegin(vecscat,ZV_QinbarR,ZV_SeqZero,                            &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat,ZV_QinbarR,ZV_SeqZero,                              &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)
if (rank==0) ZV_Qin_dam_prev=ZV_pointer(IV_dam_pos) 
call VecRestoreArrayF90(ZV_SeqZero,ZV_pointer,ierr)
!Get values from ZV_QinbarR (PETSc) into ZV_Qin_dam_prev (Fortran)


!*******************************************************************************
!Compute outflow from dams
!*******************************************************************************
!-------------------------------------------------------------------------------
!If dam module does not exist, outflow is computed from this subroutine
!-------------------------------------------------------------------------------
!if (rank==0) then 
!     ZV_Qout_dam=ZV_Qin_dam_prev
!end if

!-------------------------------------------------------------------------------
!If dam module does exist, use it
!-------------------------------------------------------------------------------
!if (rank==0) then 
!     call dam_linear(ZV_Qin_dam_prev,ZV_Qout_dam_prev,ZV_Qout_dam)
!end if

!-------------------------------------------------------------------------------
!If using Equation (6) in Doll et al. (2003, JH)
!-------------------------------------------------------------------------------
if (rank==0) then
do JS_dam_tot=1,IS_dam_tot

ZV_Qout_dam(JS_dam_tot)=(ZV_k_dam(JS_dam_tot)/ZS_TauR)                         &
                       *(ZV_S_dam(JS_dam_tot)-ZV_Smin_dam(JS_dam_tot))         &
                       *(                                                      &
                          (ZV_S_dam(JS_dam_tot)-ZV_Smin_dam(JS_dam_tot))       &
                         /(ZV_Smax_dam(JS_dam_tot)-ZV_Smin_dam(JS_dam_tot))    &
                                                         )**ZV_p_dam(JS_dam_tot)

!Equation (6) in Doll et al (2003)

ZV_S_dam(JS_dam_tot)=ZV_S_dam(JS_dam_tot)                                      &
                    +(ZV_Qin_dam_prev(JS_dam_tot)-ZV_Qout_dam(JS_dam_tot))     &
                    *ZS_TauR
!Update the storage value

if (ZV_S_dam(JS_dam_tot) <= ZV_Smin_dam(JS_dam_tot)) then
   ZV_S_dam(JS_dam_tot)=ZV_Smin_dam(JS_dam_tot)
   ZV_Qout_dam(JS_dam_tot)=ZV_Qin_dam_prev(JS_dam_tot)                         &
                          +(ZV_S_dam(JS_dam_tot)-ZV_Smin_dam(JS_dam_tot))      &
                          /ZS_TauR
end if

!Check storage value is not below minimum allowable

end do
end if

!*******************************************************************************
!Optional - Write information in stdout 
!*******************************************************************************
!if (rank==0) print *, 'ZS_TauR =',     ',', ZS_TauR,                           &
!                      'max storage =', ',', ZV_Smax_dam(1),                    &
!                      'min storage =', ',', ZV_Smin_dam(1),                    &
!                      'doll power =',  ',', ZV_p_dam(1)
!if (rank==0) print *, 'Qin_dam  =',    ',', ZV_Qin_dam(1),                     &
!                      'Qout_dam  =',   ',', ZV_Qout_dam(1)
!if (rank==0) print *, 'Qin_dam_prev  =', ',', ZV_Qin_dam_prev
!if (rank==0) print *, 'Qin_dam_prev  =', ',', ZV_Qin_dam_prev(1)
!if (rank==0) print *, 'Qout_dam_prev =', ',', ZV_Qout_dam_prev
!if (rank==0) print *, 'Qout_dam_prev =', ',', ZV_Qout_dam_prev(1)
!if (rank==0) print *, ZV_Qin_dam_prev(1), ',', ZV_Qout_dam_prev(1)
!call VecView(ZV_Qdam,PETSC_VIEWER_STDOUT_WORLD,ierr)


!*******************************************************************************
!Set values from Fortran vector into PETSc vector 
!*******************************************************************************
if (rank==0) then
     call VecSetValues(ZV_Qdam,IS_dam_bas,IV_dam_loc2,                         &
                       ZV_Qout_dam(IV_dam_index),INSERT_VALUES,ierr)
end if

call VecAssemblyBegin(ZV_Qdam,ierr)
call VecAssemblyEnd(ZV_Qdam,ierr)           


!*******************************************************************************
!Update ZV_Qout_dam_prev - After calling dam_linear to not override init. values 
!*******************************************************************************
if (rank==0) then 
     ZV_Qout_dam_prev=ZV_Qout_dam
end if


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_get_Qdam

#else

! Dummy version
subroutine rapid_get_Qdam
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_get_Qdam

#endif
