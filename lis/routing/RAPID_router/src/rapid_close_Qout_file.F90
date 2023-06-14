!*******************************************************************************
!Subroutine - rapid_close_Qout_file 
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_close_Qout_file

!Purpose:
!Close Qout_file from Fortran/netCDF.
!Author: 
!Cedric H. David, 2013-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
use netcdf
use rapid_var, only :                                                          &
                   rank,IS_nc_status,IS_nc_id_fil_Qout
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************


!*******************************************************************************
!Close file
!*******************************************************************************
if (rank==0) IS_nc_status=NF90_CLOSE(IS_nc_id_fil_Qout)


!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_close_Qout_file

#else

! Dummy version
subroutine rapid_close_Qout_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_close_Qout_file

#endif
