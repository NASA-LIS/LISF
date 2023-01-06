!*******************************************************************************
!Subroutine - rapid_write_Qfinal_file
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_write_Qfinal_file

!Purpose:
!Write into Qfinal_file from Fortran/netCDF.
!Author: 
!Cedric H. David, 2017-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petscvec.h>
use petscvec
use netcdf
use rapid_var, only :                                                          &
                   rank,ierr,vecscat,ZV_SeqZero,ZV_pointer,                    &
                   IS_nc_status,IS_nc_id_fil_Qfinal,IS_nc_id_var_Qfinal,       &
                   IS_nc_id_var_time,                                          &
                   IV_time,IS_time,ZS_TauR,                                    &
                   IS_riv_tot,ZV_read_riv_tot,IV_riv_index,                    &
                   ZV_QoutR
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************


!*******************************************************************************
!Gather PETSc vector on processor zero
!*******************************************************************************
call VecScatterBegin(vecscat,ZV_QoutR,ZV_SeqZero,                              &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat,ZV_QoutR,ZV_SeqZero,                                &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)


!*******************************************************************************
!Get array from PETSc vector
!*******************************************************************************
if (rank==0) call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)


!*******************************************************************************
!Transform the array from size of IS_riv_bas to size of IS_riv_tot
!*******************************************************************************
if (rank==0) then
ZV_read_riv_tot=0
ZV_read_riv_tot(IV_riv_index)=ZV_pointer
end if


!*******************************************************************************
!Write data
!*******************************************************************************
if (rank==0) IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qfinal,IS_nc_id_var_Qfinal,&
                                       ZV_read_riv_tot,(/1,1/),(/IS_riv_tot,1/))

if (rank==0 .and. IV_time(1)/=-9999) then 
     !The default value for 'no data' in rapid_init.F90 is -9999 for time
     IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qfinal,IS_nc_id_var_time,          &
                               IV_time(IS_time)+int(ZS_TauR),(/1/))
     !IV_time(IS_time) is the number of seconds since epoch at the beginning of 
     !the last time step within Vlat_file. We need to shift forward by the
     !number of seconds in a time step to obtain the time corresponding to 
     !Qfinal.
end if


!*******************************************************************************
!Restore array to PETSc vector
!*******************************************************************************
if (rank==0) call VecRestoreArrayF90(ZV_SeqZero,ZV_pointer,ierr)


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_write_Qfinal_file

#else

! Dummy version
subroutine rapid_write_Qfinal_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_write_Qfinal_file

#endif
