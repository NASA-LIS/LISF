!*******************************************************************************
!Subroutine - rapid_Vlat
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_Vlat(nc,nr,runsf,runsb)
!PURPOSE
!This coupler allows to convert runoff information from a land surface model
!to a volume of water entering RAPID river reaches.

!Dec 17, 2020: Yeosang Yoon, Initial Implementation

#include <petsc/finclude/petscvec.h>
use petscvec
use netcdf

use rapid_var, only :                                              &
                   rank,ierr,IS_time,JS_time,IS_riv_bas,           &
                   ZS_TauR,ZS_crs_sma,ZS_crs_iflat,                &
                   IV_time,IM_time_bnds,IV_riv_loc1,IV_riv_index,  &
                   ZV_read_riv_tot,ZV_Vlat,                        &
                   weight_table_file,n_weight_table,               &
                   rivid,npt,idx_i,idx_j,area_sqm,lat,lon,         &
                   ZV_riv_tot_lat,ZV_riv_tot_lon

implicit none

!-------------------------------------------------------------------------------
!Declaration of variables
!-------------------------------------------------------------------------------
PetscInt                               :: ncid, var_runsf, var_runsb ! variables for netcdf

! Runoff data are in kg/m2 accumulated over a time step
PetscInt, intent(in)                   :: nc, nr
real                                   :: runsf(nc,nr)   ! surface runoff
real                                   :: runsb(nc,nr)   ! subsurface runoff

PetscInt                               :: nreach_new
PetscInt                               :: col, row                !
PetscScalar                            :: conversion_factor=0.001 !convert from kg/m^2 (i.e. mm) to m

!PetscInt,    dimension(:), allocatable :: rivid_new
PetscScalar, dimension(:), allocatable :: m3_riv                  ! inflow data to RAPID river reaches are in m3 accumulated over a time step
PetscScalar                            :: m3_riv_np

PetscInt            :: eof, status
PetscInt            :: i,j,k

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Obtain a new subset of data & Calculate water inflows
!-------------------------------------------------------------------------------
! initialization

if (rank==0) then
j=1
k=1
allocate(m3_riv(IS_riv_bas))
!allocate(rivid_new(IS_riv_bas))

do i=1,n_weight_table
     m3_riv(k)=0;

     ! check if all npoints points correspond to the same streamID
     if (npt(i) > 1) then
        ! EMK Avoid array bounds read error
        !if (i > 1 .AND. (rivid(i-1) == rivid(i)))  CYCLE
        if (i > 1) then
           if (rivid(i-1) == rivid(i)) cycle
        end if

        do j=1,npt(i)
           col=idx_i(i+j-1)+1
           row=idx_j(i+j-1)+1

           !Set negative values to zero including fill values (i.e., -9999)
           if (runsf(col,row)<0) runsf(col,row)=0
           if (runsb(col,row)<0) runsb(col,row)=0

           ! combine data
           m3_riv_np=(runsf(col,row)                        &
                     +runsb(col,row))*ZS_TauR               & !kg m-2 s-1 -> kg m-2
                     *area_sqm(i+j-1)*conversion_factor       !kg m-2 (mm) -> m
           m3_riv(k)=m3_riv(k)+m3_riv_np

           ! for lat/lon
           ZV_riv_tot_lat(k)=lat(i+j-1)
           ZV_riv_tot_lon(k)=lon(i+j-1)

        end do
!        rivid_new(k)=rivid(i)
        k=k+1
        
     else
        col=idx_i(i)+1
        row=idx_j(i)+1

        !Set negative values to zero including fill values (i.e., -9999)
        if (runsf(col,row)<0) runsf(col,row)=0
        if (runsb(col,row)<0) runsb(col,row)=0

        m3_riv(k)=(runsf(col,row)                       &
                  +runsb(col,row))*ZS_TauR              & !kg m-2 s-1 -> kg m-2
                  *area_sqm(i)*conversion_factor          !kg m-2 (mm) -> m

        ! for lat/lon
        ZV_riv_tot_lat(k)=lat(i)
        ZV_riv_tot_lon(k)=lon(i)

 !       rivid_new(k)=rivid(i)
        k=k+1

     end if
end do

nreach_new=k-1

end if

!*******************************************************************************
! open Vlat_file (from rapid_open_Vlat_file.F90)
!*******************************************************************************

!if (rank==0) then
!     open(99,file=Vlat_file,status='old')
!     close(99)
!     IS_nc_status=NF90_OPEN(Vlat_file,NF90_NOWRITE,IS_nc_id_fil_Vlat)
!     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Vlat,'m3_riv',IS_nc_id_var_Vlat)
!     if (IS_nc_status<0) IS_nc_id_var_Vlat=-9999
!     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Vlat,'time',IS_nc_id_var_time)
!     if (IS_nc_status<0) IS_nc_id_var_time=-9999
!     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Vlat,'time_bnds',                \
!                                 IS_nc_id_var_time_bnds)
!     if (IS_nc_status<0) IS_nc_id_var_time_bnds=-9999
!     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Vlat,'lon',IS_nc_id_var_lon)
!     if (IS_nc_status<0) IS_nc_id_var_lon=-9999
!     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Vlat,'lat',IS_nc_id_var_lat)
!     if (IS_nc_status<0) IS_nc_id_var_lat=-9999
!     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Vlat,'crs',IS_nc_id_var_crs)
!     if (IS_nc_status<0) IS_nc_id_var_crs=-9999
!     IS_nc_status=NF90_INQ_VARID(IS_nc_id_fil_Vlat,'m3_riv_err',               \
!                                 IS_nc_id_var_Vlat_err)
!     if (IS_nc_status<0) IS_nc_id_var_Vlat_err=-9999
     !A negative value for IS_nc_id_var_* is used if the variable doesn't exist,
     !this is because the default value of "1" might match another existing
     !variable.
!end if

!*******************************************************************************
! set metadata for Vlat_file (from rapid_meta_Vlat_file.F90)
!*******************************************************************************
!Read global attributes
!if (rank==0) then
!     IS_nc_status=NF90_GET_ATT(IS_nc_id_fil_Vlat,NF90_GLOBAL,                  \
!                  "title", YV_title)
!     IS_nc_status=NF90_GET_ATT(IS_nc_id_fil_Vlat,NF90_GLOBAL,                  \
!                  "institution", YV_institution)
!     IS_nc_status=NF90_GET_ATT(IS_nc_id_fil_Vlat,NF90_GLOBAL,                  \
!                  "comment", YV_comment)
!end if

!Read variable attributes
if (rank==0) then
!     IS_nc_status=NF90_GET_ATT(IS_nc_id_fil_Vlat,IS_nc_id_var_time,            \
!                  "units", YV_time_units)
!     IS_nc_status=NF90_GET_ATT(IS_nc_id_fil_Vlat,IS_nc_id_var_crs,             \
!                  "semi_major_axis", ZS_crs_sma)
!     IS_nc_status=NF90_GET_ATT(IS_nc_id_fil_Vlat,IS_nc_id_var_crs,             \
!                  "inverse_flattening", ZS_crs_iflat)

     ZS_crs_sma=6378137.0        ! semi major axis of the spheroid
     ZS_crs_iflat=298.257223563  ! inverse flattening of the spheroid
end if

!Read space and time variable values
! TODO: re-review lon/lat, time, time_bnds variables
!if (rank==0) then
!     if (IS_nc_id_var_lon>=0) then
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_lon,             \
!                               ZV_riv_tot_lon,(/1/),(/IS_riv_tot/))
!     end if
!     if (IS_nc_id_var_lat>=0) then
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_lat,             \
!                               ZV_riv_tot_lat,(/1/),(/IS_riv_tot/))
!     end if
!     if (IS_nc_id_var_time>=0) then
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_time,            \
!                               IV_time,(/1/),(/IS_time/))
!     end if
!     if (IS_nc_id_var_time_bnds>=0) then
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_time_bnds,       \
!                               IM_time_bnds,(/1,1/),(/2,IS_time/))
!     end if
!end if

!Read uncertainty quantification inputs, convert from volume to flow
! TODO: re-review later
!if (rank==0) then
!     if (IS_nc_id_var_Vlat_err>=0) then
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_Vlat_err,        \
!                               ZV_riv_tot_bQlat,(/1,1/),(/IS_riv_tot,1/))
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_Vlat_err,        \
!                               ZV_riv_tot_vQlat,(/1,2/),(/IS_riv_tot,1/))
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_Vlat_err,        \
!                               ZV_riv_tot_caQlat,(/1,3/),(/IS_riv_tot,1/))
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_Vlat_err,        \
!                               ZV_riv_tot_cdownQlat,(/1,4/),(/IS_riv_tot,IS_radius/))

!     ZV_riv_tot_bQlat=sign(sqrt(abs(ZV_riv_tot_bQlat)),ZV_riv_tot_bQlat)/ZS_dtUQ
!     ZV_riv_tot_vQlat=ZV_riv_tot_vQlat/(ZS_dtUQ**2)
!     ZV_riv_tot_caQlat=ZV_riv_tot_caQlat/(ZS_dtUQ**2)
!     ZV_riv_tot_cdownQlat=ZV_riv_tot_cdownQlat/(ZS_dtUQ**2)
!     end if
!end if

if (rank==0) then
!Check temporal consistency if metadata present
if (IV_time(1)/=-9999) then
     do JS_time=1,IS_time-1
     if (IV_time(JS_time+1)-IV_time(JS_time)/=int(ZS_TauR)) then
     !Checking that interval between values of the time variable is ZS_TauR
          print '(a53)','Inconsistent time intervals in namelist and Vlat_file'
          stop 99
     end if
     end do
end if

if (IM_time_bnds(1,1)/=-9999) then
     do JS_time=1,IS_time-1
     if (IM_time_bnds(1,JS_time+1)-IM_time_bnds(1,JS_time)/=int(ZS_TauR)) then
     !Checking that interval between values of the time_bnd variable is ZS_TauR
          print '(a53)','Inconsistent time intervals in namelist and Vlat_file'
          stop 99
     end if
     end do
     do JS_time=1,IS_time
     if (IM_time_bnds(2,JS_time)-IM_time_bnds(1,JS_time)/=int(ZS_TauR)) then
     !Checking that interval for each value of the time_bnd variable is ZS_TauR
          print '(a53)','Inconsistent time intervals in namelist and Vlat_file'
          stop 99
     end if
     end do
end if

end if
!*******************************************************************************
! Read Vlat_file (from rapid_read_Vlat_file.F90)
!*******************************************************************************

! read file
if (rank==0) then
!     IS_nc_status=NF90_GET_VAR(IS_nc_id_fil_Vlat,IS_nc_id_var_Vlat,            &
!                               ZV_read_riv_tot,IV_nc_start,IV_nc_count)
      ZV_read_riv_tot=m3_riv(1:nreach_new)
end if

! Set values in PETSc vector
if (rank==0) then
     call VecSetValues(ZV_Vlat,IS_riv_bas,IV_riv_loc1,                         &
                       ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
end if

! Assemble PETSc vector
call VecAssemblyBegin(ZV_Vlat,ierr)
call VecAssemblyEnd(ZV_Vlat,ierr)

end subroutine rapid_Vlat

#else

! Dummy version
subroutine rapid_Vlat
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_Vlat

#endif
