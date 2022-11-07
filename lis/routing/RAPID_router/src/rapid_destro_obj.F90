!*******************************************************************************
!Subroutine - rapid_destro_obj
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_destro_obj 

!Purpose:
!All PETSc and TAO objects need be destroyed (requirement of both mathematical 
!libraries).  PETSc and TAO also need be finalized.  This is what's done here
!Note: only finilized here, need to add destroy of vectors.
!Author: 
!Cedric H. David, 2008-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petsctao.h>
use petsctao
use rapid_var, only :                                                          &
                   IS_riv_bas,                                                 &
                   ZM_hsh_tot,ZM_hsh_bas,                                      &
                   ZM_Net,ZM_A,ZM_T,ZM_TC1,ZM_M,                               &
                   ZM_Obs,ZV_Qobs,ZV_temp1,ZV_temp2,ZV_kfac,                   &
                   ZV_k,ZV_x,ZV_p,ZV_pnorm,ZV_pfac,                            &
                   ZV_C1,ZV_C2,ZV_C3,ZV_Cdenom,                                &
                   ZV_b,ZV_babsmax,ZV_bhat,                                    &
                   ZV_Qext,ZV_Qfor,ZV_Qlat,ZV_Qhum,ZV_Qdam,                    &
                   ZV_Vext,ZV_Vfor,ZV_Vlat,                                    &
                   ZV_VinitM,ZV_QoutinitM,ZV_QoutinitO,ZV_QoutbarO,            &
                   ZV_QoutR,ZV_QoutinitR,ZV_QoutprevR,ZV_QoutbarR,ZV_QinbarR,  &
                   ZV_QoutRabsmin,ZV_QoutRabsmax,ZV_QoutRhat,                  &
                   ZV_VR,ZV_VinitR,ZV_VprevR,ZV_VbarR,ZV_VoutR,                &
                   ZV_Qobsbarrec,                                              &
                   ZV_bQlat,ZV_vQlat,ZV_caQlat,ZV_bQout,ZV_sQout,ZV_rQout,     &
                   ZV_nbuptot,                                                 &
                   ierr,ksp,vecscat,ZV_SeqZero,ZV_one,                         &
                   tao,ZV_1stIndex,ZV_2ndIndex,                                &
                   ZM_Pb,ZM_L,ZM_H,ZM_S,ZM_HPbt,ZM_HPbHt,                      &
                   ZV_Qbmean,ZV_dQeb,ksp2,ZV_QoutinitR_save,                   &
                   IS_opt_run
implicit none

!*******************************************************************************
!Destroy all objects and finalize PETSc and TAO
!*******************************************************************************
call VecDestroy(ZV_1stIndex,ierr)
call VecDestroy(ZV_2ndIndex,ierr)
call TaoDestroy(tao,ierr)

call KSPDestroy(ksp,ierr)
call KSPDestroy(ksp2,ierr)

call MatDestroy(ZM_hsh_tot,ierr)
call MatDestroy(ZM_hsh_bas,ierr)

call MatDestroy(ZM_A,ierr)
call MatDestroy(ZM_Net,ierr)
call MatDestroy(ZM_T,ierr)
call MatDestroy(ZM_TC1,ierr)
call MatDestroy(ZM_M,ierr)
call MatDestroy(ZM_Obs,ierr)

call MatDestroy(ZM_Pb,ierr)
if (IS_opt_run.ne.3) call MatDestroy(ZM_L,ierr)
if (IS_opt_run.eq.3) then
    call MatDestroy(ZM_H,ierr)
    call MatDestroy(ZM_S,ierr)
    call MatDestroy(ZM_HPbt,ierr)
    call MatDestroy(ZM_HPbHt,ierr)
end if

call VecDestroy(ZV_k,ierr)
call VecDestroy(ZV_x,ierr)
call VecDestroy(ZV_C1,ierr)
call VecDestroy(ZV_C2,ierr)
call VecDestroy(ZV_C3,ierr)
call VecDestroy(ZV_Cdenom,ierr)

call VecDestroy(ZV_b,ierr)
call VecDestroy(ZV_babsmax,ierr)
call VecDestroy(ZV_bhat,ierr)

call VecDestroy(ZV_Qext,ierr)
call VecDestroy(ZV_Qfor,ierr)
call VecDestroy(ZV_Qlat,ierr)
call VecDestroy(ZV_Qhum,ierr)
call VecDestroy(ZV_Qdam,ierr)
call VecDestroy(ZV_Vext,ierr)
call VecDestroy(ZV_Vfor,ierr)
call VecDestroy(ZV_Vlat,ierr)

call VecDestroy(ZV_QoutinitM,ierr)
call VecDestroy(ZV_QoutinitO,ierr)
call VecDestroy(ZV_QoutbarO,ierr)

call VecDestroy(ZV_QoutR,ierr)
call VecDestroy(ZV_QoutinitR,ierr)
call VecDestroy(ZV_QoutprevR,ierr)
call VecDestroy(ZV_QoutbarR,ierr)
call VecDestroy(ZV_QinbarR,ierr)
call VecDestroy(ZV_QoutRabsmin,ierr)
call VecDestroy(ZV_QoutRabsmax,ierr)
call VecDestroy(ZV_QoutRhat,ierr)

call VecDestroy(ZV_QoutinitR_save,ierr)

call VecDestroy(ZV_VinitM,ierr)

call VecDestroy(ZV_VR,ierr)
call VecDestroy(ZV_VinitR,ierr)
call VecDestroy(ZV_VprevR,ierr)
call VecDestroy(ZV_VbarR,ierr)
call VecDestroy(ZV_VoutR,ierr)

call VecDestroy(ZV_temp1,ierr)
call VecDestroy(ZV_temp2,ierr)
call VecDestroy(ZV_Qobs,ierr)
call VecDestroy(ZV_kfac,ierr)
call VecDestroy(ZV_Qobsbarrec,ierr)

call VecDestroy(ZV_one,ierr)

call VecDestroy(ZV_p,ierr)
call VecDestroy(ZV_pnorm,ierr)
call VecDestroy(ZV_pfac,ierr)

call VecDestroy(ZV_nbuptot,ierr)
call VecDestroy(ZV_bQlat,ierr)
call VecDestroy(ZV_vQlat,ierr)
call VecDestroy(ZV_caQlat,ierr)
call VecDestroy(ZV_bQout,ierr)
call VecDestroy(ZV_sQout,ierr)
call VecDestroy(ZV_rQout,ierr)

call VecDestroy(ZV_Qbmean,ierr)
call VecDestroy(ZV_dQeb,ierr)

call VecDestroy(ZV_SeqZero,ierr)
call VecScatterDestroy(vecscat,ierr)
!Need to be destroyed separately even though created together

call PetscFinalize(ierr)

!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_destro_obj 

#else

! Dummy version
subroutine rapid_destro_obj
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_destro_obj

#endif
