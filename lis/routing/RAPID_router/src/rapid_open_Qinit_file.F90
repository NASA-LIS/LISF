!*******************************************************************************
!Subroutine - rapid_open_Qinit_file
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_open_Qinit_file(Qinit_file) 

!Purpose:
!Open Qinit_file from Fortran/netCDF.
!Author: 
!Cedric H. David, 2017-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
use netcdf
use rapid_var, only :                                                          &
                   rank,IS_nc_status,IS_nc_id_fil_Qinit,IS_nc_id_var_Qinit
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************
character(len=*), intent(in):: Qinit_file


!*******************************************************************************
!Open file
!*******************************************************************************
if (rank==0) then 
     open(30,file=Qinit_file,status='old')
     close(30)
     IS_nc_status=NF90_OPEN(Qinit_file,NF90_WRITE,IS_nc_id_fil_Qinit)
     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Qinit,'Qout',IS_nc_id_var_Qinit)
end if


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_open_Qinit_file

#else

! Dummy version
subroutine rapid_open_Qinit_file
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_open_Qinit_file

#endif
