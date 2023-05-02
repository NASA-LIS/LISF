!*******************************************************************************
!Subroutine - rapid_close_V_file 
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_close_V_file

!Purpose:
!Close V_file from Fortran/netCDF.
!Author: 
!Cedric H. David, 2015-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
use netcdf
use rapid_var, only :                                                          &
                   rank,IS_nc_status,IS_nc_id_fil_V
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************


!*******************************************************************************
!Close file
!*******************************************************************************
if (rank==0) IS_nc_status=NF90_CLOSE(IS_nc_id_fil_V)


!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_close_V_file

#else

subroutine rapid_close_V_file
end subroutine rapid_close_V_file

#endif
