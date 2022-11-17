!*******************************************************************************
!Subroutine - rapid_create_Qout_file
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_create_Qout_file(Qout_file) 

!Purpose:
!Create Qout_file from Fortran/netCDF.
!Author: 
!Cedric H. David, 2013-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
use netcdf
use rapid_var, only :                                                          &
                   rank,                                                       &
                   IS_nc_status,IS_nc_id_fil_Qout,                             &
                   IS_nc_id_dim_time,IS_nc_id_dim_rivid,IV_nc_id_dim,          &
                   IS_nc_id_dim_nv,IS_nc_id_dim_nerr,                          &
                   IS_nc_id_var_Qout,IS_nc_id_var_rivid,                       &
                   IS_nc_id_var_time,IS_nc_id_var_lon,IS_nc_id_var_lat,        &
                   IS_nc_id_var_time_bnds,IS_nc_id_var_crs,                    &
                   IS_nc_id_var_Qout_err,                                      &
                   IV_riv_bas_id,IS_riv_bas,                                   &
                   YV_now,YV_version,                                          &
                   Vlat_file,                                                  &
                   YV_title,YV_institution,YV_comment,                         &
                   YV_time_units,ZS_crs_sma,ZS_crs_iflat,                      &
                   ZV_riv_tot_lon,ZV_riv_tot_lat,IV_riv_index,                 &
                   ZV_riv_bas_bQout,ZV_riv_bas_sQout,ZV_riv_bas_rQout
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************
character(len=*), intent(in):: Qout_file


!*******************************************************************************
!Create file
!*******************************************************************************
if (rank==0) then 

!-------------------------------------------------------------------------------
!Create file
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_CREATE(Qout_file,NF90_CLOBBER,IS_nc_id_fil_Qout)

!-------------------------------------------------------------------------------
!Define dimensions
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_DEF_DIM(IS_nc_id_fil_Qout,'time',NF90_UNLIMITED,        &
                               IS_nc_id_dim_time)
     IS_nc_status=NF90_DEF_DIM(IS_nc_id_fil_Qout,'rivid',IS_riv_bas,           &
                               IS_nc_id_dim_rivid)
     IS_nc_status=NF90_DEF_DIM(IS_nc_id_fil_Qout,'nv',2,                       &
                               IS_nc_id_dim_nv)
     IS_nc_status=NF90_DEF_DIM(IS_nc_id_fil_Qout,'nerr',3,                     &
                               IS_nc_id_dim_nerr)
     IV_nc_id_dim(1)=IS_nc_id_dim_rivid
     IV_nc_id_dim(2)=IS_nc_id_dim_time

!-------------------------------------------------------------------------------
!Define variables
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'Qout',NF90_REAL,             &
                               IV_nc_id_dim,IS_nc_id_var_Qout)
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'Qout_err',NF90_REAL,         &
                               (/IS_nc_id_dim_rivid,IS_nc_id_dim_nerr/),       &
                               IS_nc_id_var_Qout_err)
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'rivid',NF90_INT,             &
                               IS_nc_id_dim_rivid,IS_nc_id_var_rivid)
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'time',NF90_INT,              &
                               IS_nc_id_dim_time,IS_nc_id_var_time)
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'time_bnds',NF90_INT,         &
                               (/IS_nc_id_dim_nv,IS_nc_id_dim_time/),          &
                               IS_nc_id_var_time_bnds)
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'lon',NF90_DOUBLE,            &
                               IS_nc_id_dim_rivid,IS_nc_id_var_lon)
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'lat',NF90_DOUBLE,            &
                               IS_nc_id_dim_rivid,IS_nc_id_var_lat)
     IS_nc_status=NF90_DEF_VAR(IS_nc_id_fil_Qout,'crs',NF90_INT,               &
                               IS_nc_id_var_crs)

!-------------------------------------------------------------------------------
!Define variable attributes
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout,            &
                               'long_name','average river water discharge '    &
                               // 'downstream of each river reach')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout,            &
                               'units','m3 s-1')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout,            &
                               'coordinates','lon lat')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout,            &
                               'grid_mapping','crs')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout,            &
                               'cell_methods','time: mean')

     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               'long_name','average river water discharge '    &
                               // 'uncertainty downstream of each river reach')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               'units','m3 s-1')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               'coordinates','lon lat')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               'grid_mapping','crs')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               'cell_methods','time: mean')

     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_rivid,           &
                               'long_name','unique identifier for each river ' &
                               // 'reach')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_rivid,           &
                               'units','1')  
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_rivid,           &
                               'cf_role','timeseries_id')  

     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_time,            &
                               'standard_name','time')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_time,            &
                               'long_name','time')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_time,            &
                               'units',YV_time_units)
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_time,            &
                               'calendar','gregorian')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_time,            &
                               'axis','T')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_time,            &
                               'bounds','time_bnds')

     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lon,             &
                               'standard_name','longitude')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lon,             &
                               'long_name','longitude of a point related '     &
                               // 'to each river reach')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lon,             &
                               'units','degrees_east')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lon,             &
                               'axis','X')

     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lat,             &
                               'standard_name','latitude')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lat,             &
                               'long_name','latitude of a point related '      &
                               // 'to each river reach')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lat,             &
                               'units','degrees_north')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_lat,             &
                               'axis','Y')

     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_crs,             &
                               'grid_mapping_name','latitude_longitude')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_crs,             &
                               'semi_major_axis',ZS_crs_sma)
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,IS_nc_id_var_crs,             &
                               'inverse_flattening',ZS_crs_iflat)

!-------------------------------------------------------------------------------
!Define global attributes
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'Conventions','CF-1.6')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'title',YV_title)
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'institution',YV_institution)
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'source','RAPID: '//TRIM(YV_version)//          &
                               ', water inflow: '//Vlat_file)
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'history','date_created: '//YV_now)
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'references','https://github.com/c-h-david/ra'//&
                                'pid/, http://dx.doi.org/10.1175/2011JHM1345.1')
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'comment',YV_comment)
     IS_nc_status=NF90_PUT_ATT(IS_nc_id_fil_Qout,NF90_GLOBAL,                  &
                               'featureType','timeSeries')

!-------------------------------------------------------------------------------
!End definition and close file
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_ENDDEF(IS_nc_id_fil_Qout)
     IS_nc_status=NF90_CLOSE(IS_nc_id_fil_Qout)

!-------------------------------------------------------------------------------
!End create file
!-------------------------------------------------------------------------------
end if


!*******************************************************************************
!Populate variables
!*******************************************************************************
if (rank==0) then 

!-------------------------------------------------------------------------------
!Open file
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_OPEN(Qout_file,NF90_WRITE,IS_nc_id_fil_Qout)

!-------------------------------------------------------------------------------
!Populate variables
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qout,IS_nc_id_var_rivid,           &
                               IV_riv_bas_id)
     if (ZV_riv_tot_lon(1)/=-9999) then 
     !The default value for 'no data' in rapid_init.F90 is -9999 for longitude
     IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qout,IS_nc_id_var_lon,             &
                               ZV_riv_tot_lon(IV_riv_index))
     end if
     if (ZV_riv_tot_lat(1)/=-9999) then 
     !The default value for 'no data' in rapid_init.F90 is -9999 for latitude
     IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qout,IS_nc_id_var_lat,             &
                               ZV_riv_tot_lat(IV_riv_index))
     end if
     if (maxval(ZV_riv_bas_bQout)/=0) then 
     !The default value for 'no data' in rapid_init.F90 is 0 for bias of Qout
     IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               ZV_riv_bas_bQout,(/1,1/),(/IS_riv_bas,1/))
     end if
     if (maxval(ZV_riv_bas_sQout)/=0) then 
     !The default value for 'no data' in rapid_init.F90 is 0 for stderr of Qout
     IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               ZV_riv_bas_sQout,(/1,2/),(/IS_riv_bas,1/))
     end if
     if (maxval(ZV_riv_bas_rQout)/=0) then 
     !The default value for 'no data' in rapid_init.F90 is 0 for RMSE of Qout
     IS_nc_status=NF90_PUT_VAR(IS_nc_id_fil_Qout,IS_nc_id_var_Qout_err,        &
                               ZV_riv_bas_rQout,(/1,3/),(/IS_riv_bas,1/))
     end if

!-------------------------------------------------------------------------------
!Close file
!-------------------------------------------------------------------------------
     IS_nc_status=NF90_CLOSE(IS_nc_id_fil_Qout)

!-------------------------------------------------------------------------------
!End populate variables
!-------------------------------------------------------------------------------
end if


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_create_Qout_file

#else

subroutine rapid_create_Qout_file
end subroutine rapid_create_Qout_file

#endif
