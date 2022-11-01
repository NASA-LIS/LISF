!*******************************************************************************
!Subroutine - rapid_set_Qext0
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_set_Qext0

!Purpose:
!This subroutine is only useful if a dam model is used and its goal is to 
!properly initialize the flow of water into the dams.
!The inflow of water ZV_Qin_dam_prev from the river network and from outside of 
!the river network into the dams is computed based on ZV_QoutbarR and ZV_Qext
!in the subroutine rapid_get_Qdam.F90. 
!Therefore, one has to inject the initial value of ZV_Qin_dam_prev (ZV_Qin_dam0) 
!into either ZV_QoutbarR or ZV_Qext otherwise the initial value will be 
!overwritten in rapid_get_Qdam.F90. The latter is used here (through ZV_Qdam)
!since the modifications made on the network matrix make it difficult to use
!ZV_Qin_dam_prev without creating a new variable.
!Author: 
!Cedric H. David, 2013-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petscvec.h>
use petscvec
use rapid_var, only:                                                           &
                   rank,ierr,IS_one,ZS_one,                                    &
                   ZV_Qdam,ZV_Qext,                                            &
                   IS_dam_tot,JS_dam_tot,IV_dam_pos,                           & 
                   ZV_Qin_dam0,                                                &
                   ZV_Qout_dam_prev,ZV_Qout_dam0
implicit none


!*******************************************************************************
!Set Qdam to zero, because this is called at the beginning of every phiroutine
!*******************************************************************************
call VecSet(ZV_Qdam,0*ZS_one,ierr)                            !Qdam=0
!call VecView(ZV_Qdam,PETSC_VIEWER_STDOUT_WORLD,ierr)


!*******************************************************************************
!Set values of Qin_dam0 into Qdam to allow proper initialization
!*******************************************************************************
if (rank==0) then
     do JS_dam_tot=1,IS_dam_tot

if (IV_dam_pos(JS_dam_tot)/=0) then
     call VecSetValues(ZV_Qdam,IS_one,IV_dam_pos(JS_dam_tot)-1,                &
                       ZV_Qin_dam0(JS_dam_tot),INSERT_VALUES,ierr)
     !print *, IV_dam_pos(JS_dam_tot)-1, ZV_Qin_dam0(JS_dam_tot)
end if

     end do
end if

call VecAssemblyBegin(ZV_Qdam,ierr)
call VecAssemblyEnd(ZV_Qdam,ierr)      
!call VecView(ZV_Qdam,PETSC_VIEWER_STDOUT_WORLD,ierr)
!the values of Qindam0 are set here where the dams are, not downstream of them


!*******************************************************************************
!Copy Qdam into Qext and reset Qdam to zero
!*******************************************************************************
call VecCopy(ZV_Qdam,ZV_Qext,ierr)                            !Qext=Qdam
call VecSet(ZV_Qdam,0*ZS_one,ierr)                            !Qdam=0
!call VecView(ZV_Qext,PETSC_VIEWER_STDOUT_WORLD,ierr)


!*******************************************************************************
!Initialize Qout_dam_prev again or its values differ with each phiroutine call
!*******************************************************************************
if (rank==0) then
     ZV_Qout_dam_prev=ZV_Qout_dam0
end if


!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_set_Qext0

#else

! Dummy version
subroutine rapid_set_Qext0
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_set_Qext0

#endif
