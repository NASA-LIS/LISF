!*******************************************************************************
!Subroutine - rapid_init 
!*******************************************************************************

#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_init

!Purpose:
!This subroutine allows to initialize RAPID for both regular runs and 
!optimization runs, by performing slightly different tasks depending on what 
!option is chosen.  
!Initialization tasks common to all RAPID options:
!     -Read namelist file (sizes of domain, duration, file names, options, etc.)  
!     -Compute number of time steps based on durations
!     -Allocate Fortran arrays
!     -Create all PETSc and TAO objects 
!     -Print information and warnings
!     -Determine IDs for various computing cores
!     -Compute helpful arrays 
!     -Compute the network matrix
!     -Initialize values of flow and volume for main procedure
!Initialization tasks specific to Option 1
!     -Copy main initial flow and vol to routing initial flow and vol
!     -Read k and x 
!     -Compute linear system matrix
!Initialization tasks specific to Option 2
!     -Copy main initial flow to optimization initial flow
!     -Compute the observation matrix
!     -Read kfac and Qobsbarrec
!     -Set initial values for the vector pnorm
!Author: 
!Cedric H. David, 2012-2020.

! !REVISION HISTORY:
! 07 May 2021: Yeosang Yoon: Update log message
! 07 Jul 2021: Yeosang Yoon: Add memory allocation for weight table

!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************
#include <petsc/finclude/petscksp.h>
use petscksp
use rapid_var, only :                                                          &
                   IS_riv_tot,IS_riv_bas,                                      &
                   IV_riv_bas_id,IV_riv_index,IV_riv_loc1,IV_riv_tot_id,       &
                   IV_down,IV_nbup,IM_up,IM_index_up,IS_max_up,                &
                   IV_nz,IV_dnz,IV_onz,                                        &
                   BS_opt_Qinit,BS_opt_Qfinal,BS_opt_V,BS_opt_influence,       & 
                   BS_opt_dam,BS_opt_for,BS_opt_hum,BS_opt_uq,                 &
                   IS_opt_run,IS_opt_routing,IS_opt_phi,                       &
                   ZV_read_riv_tot,ZV_read_obs_tot,ZV_read_hum_tot,            &
                   ZV_read_for_tot,ZV_read_dam_tot,                            &
                   ZS_TauM,ZS_TauO,ZS_TauR,ZS_dtO,ZS_dtR,ZS_dtM,ZS_dtF,ZS_dtH, &
                   IS_obs_tot,IS_obs_use,IS_obs_bas,                           &
                   IV_obs_tot_id,IV_obs_use_id,                                &
                   IV_obs_index,IV_obs_loc1,                                   &
                   IS_hum_tot,IS_hum_use,                                      &
                   IV_hum_tot_id,IV_hum_use_id,                                &
                   IS_for_tot,IS_for_use,                                      &
                   IV_for_tot_id,IV_for_use_id,                                &
                   IS_dam_tot,JS_dam_tot,IS_dam_use,                           &
                   IV_dam_tot_id,IV_dam_use_id,                                &
                   ZV_Qin_dam,ZV_Qout_dam,ZV_Qin_dam_prev,ZV_Qout_dam_prev,    &
                   ZV_Qin_dam0,ZV_Qout_dam0,                                   &
                   ZV_S_dam,ZV_Smax_dam,ZV_Smin_dam,                           &
                   ZV_k_dam,ZV_p_dam,                                          &
                   ZV_riv_tot_bQlat,ZV_riv_tot_vQlat,ZV_riv_tot_caQlat,        &
                   ZV_riv_bas_bQout,ZV_riv_bas_sQout,ZV_riv_bas_rQout,         &
                   ZV_riv_bas_bV,ZV_riv_bas_sV,ZV_riv_bas_rV,                  &
                   ZV_QoutinitM,ZV_QoutinitO,ZV_QoutinitR,                     &
                   ZV_VinitM,ZV_VinitR,                                        &
                   ZV_babsmax,ZV_QoutRabsmin,ZV_QoutRabsmax,                   &
                   IS_M,IS_O,IS_R,IS_RpO,IS_RpM,IS_RpF,IS_RpH,IS_time,         &
                   kfac_file,x_file,k_file,Vlat_file,Qinit_file,               &
                   Qobsbarrec_file,                                            &
                   ZS_Qout0,ZS_V0,                                             &
                   ZV_Qobsbarrec,dam_file,                                     &
                   ZV_k,ZV_x,ZV_kfac,ZV_pnorm,                                 &
                   ZS_knorm_init,ZS_xnorm_init,ZS_kfac,ZS_xfac,                &
                   ZV_C1,ZV_C2,ZV_C3,ZM_A,                                     &
                   IV_now,YV_now,YV_version,                                   &
                   ZV_riv_tot_lon,ZV_riv_tot_lat,IV_time,IM_time_bnds,         &
                   ierr,ksp,rank,IS_one,ZS_one,                                &
                   IS_radius,ZV_riv_tot_cdownQlat,                             &
                   IV_nbrows,IV_lastrow,                                       &
                   n_weight_table,rivid,npt,idx_i,idx_j,area_sqm,lat,lon

use LIS_logMod

implicit none


!*******************************************************************************
!Initialization procedure common to all options
!*******************************************************************************

!-------------------------------------------------------------------------------
!Get current date-time and format using ISO 8601 international standard 
!-------------------------------------------------------------------------------
call date_and_time(VALUES=IV_now)

write(YV_now,'(i4,a1,i0.2,a1,i0.2,a1,i0.2,a1,i0.2,a1,i0.2,a1,i0.2,a1,i0.2)')   &
              IV_now(1), '-', IV_now(2), '-', IV_now(3), 'T',                  &
              IV_now(5), ':', IV_now(6), ':', IV_now(7),                       &
              '*', abs(IV_now(4)/60), ':', mod(abs(IV_now(4)),60) 
if (IV_now(4)>=0) then
     YV_now=YV_now(1:19)//'+'//YV_now(21:25)
else
     YV_now=YV_now(1:19)//'-'//YV_now(21:25)
end if
!Using the ISO 8601 international standard: 2016-01-31T16:45:46+00:00

!-------------------------------------------------------------------------------
!Get the version of RAPID determined during build
!-------------------------------------------------------------------------------
#ifdef RAPID_VERSION
     YV_version=RAPID_VERSION
#else
     YV_version='v1.8.0'    !Yeosang Yoon
#endif
!Compilation examples: -D RAPID_VERSION="'v1.4.0'"  
!                      -D RAPID_VERSION="'20131114'" 

!-------------------------------------------------------------------------------
!Read name list
!-------------------------------------------------------------------------------
call rapid_read_namelist

!-------------------------------------------------------------------------------
!Compute number of time steps
!-------------------------------------------------------------------------------
IS_M=int(ZS_TauM/ZS_dtM)
IS_O=int(ZS_TauO/ZS_dtO)
IS_R=int(ZS_TauR/ZS_dtR)
IS_RpO=int(ZS_dtO/ZS_TauR)
IS_RpM=int(ZS_dtM/ZS_TauR)
IS_RpF=int(ZS_dtF/ZS_TauR)
IS_RpH=int(ZS_dtH/ZS_TauR)
IS_time=IS_M*IS_RpM

!-------------------------------------------------------------------------------
!Allocate Fortran arrays
!-------------------------------------------------------------------------------
allocate(IV_riv_bas_id(IS_riv_bas))
allocate(IV_riv_index(IS_riv_bas))
allocate(IV_riv_loc1(IS_riv_bas))

allocate(IV_riv_tot_id(IS_riv_tot))
allocate(IV_down(IS_riv_tot))
allocate(IV_nbup(IS_riv_tot))
allocate(IM_up(IS_riv_tot,IS_max_up))
allocate(IM_index_up(IS_riv_tot,IS_max_up))

allocate(IV_nz(IS_riv_bas))
allocate(IV_dnz(IS_riv_bas))
allocate(IV_onz(IS_riv_bas))

allocate(ZV_read_riv_tot(IS_riv_tot))

allocate(ZV_riv_tot_lon(IS_riv_tot))
allocate(ZV_riv_tot_lat(IS_riv_tot))

allocate(IV_time(IS_time))
allocate(IM_time_bnds(2,IS_time))

if ((IS_opt_run==2).or.(IS_opt_run==3).or.(IS_opt_run==4)) then
     allocate(IV_obs_tot_id(IS_obs_tot))
     allocate(IV_obs_use_id(IS_obs_use))
     allocate(ZV_read_obs_tot(IS_obs_tot))
end if

if (BS_opt_hum) then
     allocate(IV_hum_tot_id(IS_hum_tot))
     allocate(IV_hum_use_id(IS_hum_use))
     allocate(ZV_read_hum_tot(IS_hum_tot))
end if

if (BS_opt_for) then
     allocate(IV_for_tot_id(IS_for_tot))
     allocate(IV_for_use_id(IS_for_use))
     allocate(ZV_read_for_tot(IS_for_tot))
end if

if (BS_opt_dam) then
     allocate(IV_dam_tot_id(IS_dam_tot))
     allocate(IV_dam_use_id(IS_dam_use))
     allocate(ZV_read_dam_tot(IS_dam_tot))
     allocate(ZV_Qin_dam(IS_dam_tot))
     allocate(ZV_Qin_dam_prev(IS_dam_tot))
     allocate(ZV_Qout_dam(IS_dam_tot))
     allocate(ZV_Qout_dam_prev(IS_dam_tot))
     allocate(ZV_Qin_dam0(IS_dam_tot))
     allocate(ZV_Qout_dam0(IS_dam_tot))
     allocate(ZV_k_dam(IS_dam_tot))
     allocate(ZV_p_dam(IS_dam_tot))
     allocate(ZV_S_dam(IS_dam_tot))
     allocate(ZV_Smin_dam(IS_dam_tot))
     allocate(ZV_Smax_dam(IS_dam_tot))
end if

allocate(ZV_riv_tot_bQlat(IS_riv_tot))
allocate(ZV_riv_tot_vQlat(IS_riv_tot))
allocate(ZV_riv_tot_caQlat(IS_riv_tot))
!Used in rapid_meta_Vlat_file regardless of BS_opt_uq

allocate(ZV_riv_bas_bQout(IS_riv_bas))
allocate(ZV_riv_bas_sQout(IS_riv_bas))
allocate(ZV_riv_bas_rQout(IS_riv_bas))
!Used in rapid_create_Qout_file regardless of BS_opt_uq

allocate(ZV_riv_bas_bV(IS_riv_bas))
allocate(ZV_riv_bas_sV(IS_riv_bas))
allocate(ZV_riv_bas_rV(IS_riv_bas))
!Used in rapid_create_V_file regardless of BS_opt_uq

allocate(ZV_riv_tot_cdownQlat(IS_riv_tot,IS_radius))
!Used in rapid_meta_Vlat_file and rapid_cov_mat for data assimilation

allocate(IV_nbrows(IS_riv_bas))
allocate(IV_lastrow(IS_riv_bas))
!Used in rapid_mus_mat and rapid_runoff2streamflow_mat for data assimilation

! weight table (Yeosang Yoon; 08 Jul 2021)
allocate(rivid(n_weight_table))
allocate(area_sqm(n_weight_table))
allocate(idx_i(n_weight_table))
allocate(idx_j(n_weight_table))
allocate(npt(n_weight_table))
allocate(lat(n_weight_table))
allocate(lon(n_weight_table))

!-------------------------------------------------------------------------------
!Make sure some Fortran arrays are initialized to zero
!-------------------------------------------------------------------------------
ZV_riv_tot_lon=-9999
ZV_riv_tot_lat=-9999
IV_time=-9999
IM_time_bnds=-9999
!The value of '-9999' is used here as 'No Data' in case not present in Vlat_file

if (BS_opt_dam) then
     ZV_Qin_dam0 =0
     ZV_Qout_dam0=0
     ZV_S_dam=0
end if
!These are not populated anywhere before being used and hold meaningless values

ZV_riv_tot_bQlat=0
ZV_riv_tot_vQlat=0
ZV_riv_tot_caQlat=0
!Used in rapid_meta_Vlat_file regardless of BS_opt_uq

ZV_riv_bas_bQout=0
ZV_riv_bas_sQout=0
ZV_riv_bas_rQout=0
!Used in rapid_create_Qout_file regardless of BS_opt_uq

ZV_riv_bas_bV=0
ZV_riv_bas_sV=0
ZV_riv_bas_rV=0
!Used in rapid_create_V_file regardless of BS_opt_uq

ZV_riv_tot_cdownQlat=0
!Used in rapid_meta_Vlat_file and rapid_cov_mat for data assimilation

!-------------------------------------------------------------------------------
!Initialize libraries and create objects common to all options
!-------------------------------------------------------------------------------
call rapid_create_obj
!Initialize libraries and create PETSc and TAO objects (Mat,Vec,taoapp...)

!-------------------------------------------------------------------------------
!Prints information about current model run based on info from namelist
!-------------------------------------------------------------------------------
!call PetscPrintf(PETSC_COMM_WORLD,'--------------------------'//char(10),ierr)
!if (rank==0)                                               print '(a70)',      &
!       'RAPID: ' // YV_version //       '                                      ' 
!if (rank==0)                                               print '(a70)',      &
!       'Current ISO 8601 time: ' // YV_now //          '                       ' 
!if (rank==0 .and. .not. BS_opt_Qinit)                      print '(a70)',      &
!       'Not reading initial flows from a file                                  '
!if (rank==0 .and. BS_opt_Qinit)                            print '(a70)',      &
!       'Reading initial flows from a file                                      '
!if (rank==0 .and. .not. BS_opt_Qfinal .and. IS_opt_run==1) print '(a70)',      &
!       'Not writing final flows into a file                                    '
!if (rank==0 .and. BS_opt_Qfinal .and. IS_opt_run==1)       print '(a70)',      &
!       'Writing final flows into a file                                        '
!if (rank==0 .and. .not. BS_opt_V .and. IS_opt_run==1)      print '(a70)',      &
!       'Not computing water volumes in river reaches                           '
!if (rank==0 .and. BS_opt_V .and. IS_opt_run==1)            print '(a70)',      &
!       'Computing water volumes in river reaches                               '
!if (rank==0 .and. .not. BS_opt_for)                        print '(a70)',      &
!       'Not using forcing                                                      '
!if (rank==0 .and. BS_opt_for)                              print '(a70)',      &
!       'Using forcing                                                          '
!if (rank==0 .and. .not. BS_opt_hum)                        print '(a70)',      &
!       'Not using human-induced flows                                          '
!if (rank==0 .and. BS_opt_hum)                              print '(a70)',      &
!       'Using human-induced flows                                              '
!if (rank==0 .and. .not. BS_opt_uq .and. IS_opt_run==1)     print '(a70)',      &
!       'Not quantifying uncertainty                                            '
!if (rank==0 .and. BS_opt_uq .and. IS_opt_run==1)           print '(a70)',      &
!       'Quantifying uncertainty                                                '
!if (rank==0 .and. IS_opt_routing==1)                       print '(a70)',      &
!       'Routing with matrix-based Muskingum method                             '
!if (rank==0 .and. IS_opt_routing==2)                       print '(a70)',      &
!       'Routing with traditional Muskingum method                              '
!if (rank==0 .and. IS_opt_routing==3)                       print '(a70)',      &
!       'Routing with matrix-based Muskingum method using transboundary matrix  '
!if (rank==0 .and. IS_opt_routing==4)                       print '(a70)',      &
!       'Routing with matrix-based Muskingum method using Muskingum operator    '
!if (rank==0 .and. IS_opt_run==1)                           print '(a70)',      &
!       'RAPID mode: computing flowrates                                        '
!if (rank==0 .and. IS_opt_run==2 .and. IS_opt_phi==1)       print '(a70)',      &
!       'RAPID mode: optimizing parameters, using phi1                          ' 
!if (rank==0 .and. IS_opt_run==2 .and. IS_opt_phi==2)       print '(a70)',      &
!       'RAPID mode: optimizing parameters, using phi2                          ' 
!if (rank==0 .and. IS_opt_run==3)                           print '(a70)',      &
!       'RAPID mode: data assimilation                                          '
!if (rank==0 .and. IS_opt_run==3)                           print '(a70)',      &
!       'RAPID mode: data assimilation with simplified observation operator     '
!!if (rank==0)                                               print '(a10,a60)',  &
!!       'Using    :', Vlat_file 
!if (rank==0 .and. IS_opt_run==1)                           print '(a10,a60)',  &
!       'Using    :',k_file 
!if (rank==0 .and. IS_opt_run==1)                           print '(a10,a60)',  &
!       'Using    :',x_file 
!if (rank==0 .and. IS_opt_run==2)                           print '(a10,a60)',  &
!       'Using    :',kfac_file 
!call PetscPrintf(PETSC_COMM_WORLD,'--------------------------'//char(10),ierr)

! Yeosang Yoon, update log message to fit LIS 
if (rank==0)                                                                                            &
    write(LIS_logunit,*) '[INFO] RAPID: Initializing RAPID runs.........................................'
if (rank==0)                                                                                            &
    write(LIS_logunit,*) '[INFO] RAPID: ' // trim(YV_version)
if (rank==0)                                                                                            &
    write(LIS_logunit,*) '[INFO] Current ISO 8601 time: ' // YV_now //          '                       '
if (rank==0 .and. .not. BS_opt_Qinit)                                                                   &
    write(LIS_logunit,*) '[INFO] Not reading initial flows from a file                                  '
if (rank==0 .and. BS_opt_Qinit)                                                                         &
    write(LIS_logunit,*) '[INFO] Reading initial flows from a file                                      '
if (rank==0 .and. .not. BS_opt_Qfinal .and. IS_opt_run==1)                                              &
    write(LIS_logunit,*) '[INFO] Not writing final flows into a file                                    '
if (rank==0 .and. BS_opt_Qfinal .and. IS_opt_run==1)                                                    &
    write(LIS_logunit,*) '[INFO] Writing final flows into a file                                        '
if (rank==0 .and. .not. BS_opt_V .and. IS_opt_run==1)                                                   &
    write(LIS_logunit,*) '[INFO] Not computing water volumes in river reaches                           '
if (rank==0 .and. BS_opt_V .and. IS_opt_run==1)                                                         &
    write(LIS_logunit,*) '[INFO] Computing water volumes in river reaches                               '
if (rank==0 .and. .not. BS_opt_for)                                                                     &
    write(LIS_logunit,*) '[INFO] Not using forcing                                                      '
if (rank==0 .and. BS_opt_for)                                                                           &
    write(LIS_logunit,*) '[INFO] Using forcing                                                          '
if (rank==0 .and. .not. BS_opt_hum)                                                                     & 
    write(LIS_logunit,*) '[INFO] Not using human-induced flows                                          '
if (rank==0 .and. BS_opt_hum)                                                                           &
    write(LIS_logunit,*) '[INFO] Using human-induced flows                                              '
if (rank==0 .and. .not. BS_opt_uq .and. IS_opt_run==1)                                                  &
    write(LIS_logunit,*) '[INFO] Not quantifying uncertainty                                            '
if (rank==0 .and. BS_opt_uq .and. IS_opt_run==1)                                                        &
    write(LIS_logunit,*) '[INFO] Quantifying uncertainty                                                '
if (rank==0 .and. IS_opt_routing==1)                                                                    &
    write(LIS_logunit,*) '[INFO] Routing with matrix-based Muskingum method                             '
if (rank==0 .and. IS_opt_routing==2)                                                                    &
    write(LIS_logunit,*) '[INFO] Routing with traditional Muskingum method                              '
if (rank==0 .and. IS_opt_routing==3)                                                                    &
    write(LIS_logunit,*) '[INFO] Routing with matrix-based Muskingum method using transboundary matrix  '
if (rank==0 .and. IS_opt_routing==4)                                                                    &
    write(LIS_logunit,*) '[INFO] Routing with matrix-based Muskingum method using Muskingum operator    '
if (rank==0 .and. IS_opt_run==1)                                                                        &
    write(LIS_logunit,*) '[INFO] RAPID mode: computing flowrates                                        '
if (rank==0 .and. IS_opt_run==2 .and. IS_opt_phi==1)                                                    &
    write(LIS_logunit,*) '[INFO] RAPID mode: optimizing parameters, using phi1                          '
if (rank==0 .and. IS_opt_run==2 .and. IS_opt_phi==2)                                                    &
    write(LIS_logunit,*) '[INFO] RAPID mode: optimizing parameters, using phi2                          '
if (rank==0 .and. IS_opt_run==3)                                                                        &
    write(LIS_logunit,*) '[INFO] RAPID mode: data assimilation                                          '
if (rank==0 .and. IS_opt_run==3)                                                                        &
    write(LIS_logunit,*) '[INFO] RAPID mode: data assimilation with simplified observation operator     '
if (rank==0 .and. IS_opt_run==1)                                                                        &
    write(LIS_logunit,*) '[INFO] Using    :',trim(k_file)
if (rank==0 .and. IS_opt_run==1)                                                                        &
    write(LIS_logunit,*) '[INFO] Using    :',trim(x_file)
if (rank==0 .and. IS_opt_run==2)                                                                        &
    write(LIS_logunit,*) '[INFO] Using    :',trim(kfac_file)
!-------------------------------------------------------------------------------
!Calculate helpful arrays
!-------------------------------------------------------------------------------
call rapid_arrays

!-------------------------------------------------------------------------------
!Calculate Network matrix
!-------------------------------------------------------------------------------
call rapid_net_mat

!-------------------------------------------------------------------------------
!Breaks connections in Network matrix
!-------------------------------------------------------------------------------
if (BS_opt_for .or. BS_opt_dam) call rapid_net_mat_brk

!-------------------------------------------------------------------------------
!calculates or set initial flows and volumes
!-------------------------------------------------------------------------------
if (.not. BS_opt_Qinit) then
call VecSet(ZV_QoutinitM,ZS_Qout0,ierr)
end if

if (BS_opt_Qinit) then
call rapid_open_Qinit_file(Qinit_file)
call rapid_read_Qinit_file
call rapid_close_Qinit_file
end if

call VecSet(ZV_VinitM,ZS_V0,ierr)
!Set initial volumes for Main procedure

!-------------------------------------------------------------------------------
!Initialize default values for ZV_QoutRabsmin, ZV_QoutRabsmax and ZV_babsmax
!-------------------------------------------------------------------------------
if (BS_opt_influence) then
call VecSet(ZV_babsmax    ,ZS_one*0        ,ierr)
call VecSet(ZV_QoutRabsmin,ZS_one*999999999,ierr)
call VecSet(ZV_QoutRabsmax,ZS_one*0        ,ierr)
end if

!*******************************************************************************
!Initialization procedure for OPTION 1
!*******************************************************************************
if (IS_opt_run==1) then

!-------------------------------------------------------------------------------
!copy main initial values into routing initial values 
!-------------------------------------------------------------------------------
call VecCopy(ZV_QoutinitM,ZV_QoutinitR,ierr)
call VecCopy(ZV_VinitM,ZV_VinitR,ierr)

!-------------------------------------------------------------------------------
!Read/set k and x
!-------------------------------------------------------------------------------
open(20,file=k_file,status='old')
read(20,*) ZV_read_riv_tot
call VecSetValues(ZV_k,IS_riv_bas,IV_riv_loc1,                                 &
                  ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_k,ierr)
call VecAssemblyEnd(ZV_k,ierr)
close(20)
!get values for k in a file and create the corresponding ZV_k vector

open(21,file=x_file,status='old')
read(21,*) ZV_read_riv_tot
call VecSetValues(ZV_x,IS_riv_bas,IV_riv_loc1,                                 &
                  ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_x,ierr)
call VecAssemblyEnd(ZV_x,ierr)
close(21)
!get values for x in a file and create the corresponding ZV_x vector

!-------------------------------------------------------------------------------
!Read dam_file with storage information and parameters
!-------------------------------------------------------------------------------
if (BS_opt_dam) then
open(24,file=dam_file,status='old')
do JS_dam_tot=1,IS_dam_tot
     read(24,*) ZV_Smax_dam(JS_dam_tot),ZV_Smin_dam(JS_dam_tot),               &
                ZV_p_dam(JS_dam_tot),ZV_k_dam(JS_dam_tot)
end do
close(24)

ZV_S_dam=ZV_Smin_dam
!Initialize all dam storage to the minimum storage
end if

!-------------------------------------------------------------------------------
!Compute routing parameters and linear system matrix
!-------------------------------------------------------------------------------
call rapid_routing_param(ZV_k,ZV_x,ZV_C1,ZV_C2,ZV_C3,ZM_A)
!calculate Muskingum parameters and matrix ZM_A

call KSPSetOperators(ksp,ZM_A,ZM_A,ierr)
!Set KSP to use matrix ZM_A

!-------------------------------------------------------------------------------
!Calculate Muskingum matrix
!-------------------------------------------------------------------------------
if (IS_opt_routing==4) call rapid_mus_mat

!-------------------------------------------------------------------------------
!End of initialization procedure for OPTION 1
!-------------------------------------------------------------------------------
end if


!*******************************************************************************
!Initialization procedure for OPTION 2
!*******************************************************************************
if (IS_opt_run==2) then

!-------------------------------------------------------------------------------
!Create observation matrix
!-------------------------------------------------------------------------------
call rapid_obs_mat
!Create observation matrix

!-------------------------------------------------------------------------------
!copy main initial values into optimization initial values 
!-------------------------------------------------------------------------------
call VecCopy(ZV_QoutinitM,ZV_QoutinitO,ierr)
!copy initial main variables into initial optimization variables

!-------------------------------------------------------------------------------
!Read/set kfac, xfac and Qobsbarrec
!-------------------------------------------------------------------------------
open(22,file=kfac_file,status='old')
read(22,*) ZV_read_riv_tot
close(22)
call VecSetValues(ZV_kfac,IS_riv_bas,IV_riv_loc1,                              &
                  ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
                  !only looking at basin, doesn't have to be whole domain here 
call VecAssemblyBegin(ZV_kfac,ierr)
call VecAssemblyEnd(ZV_kfac,ierr)  
!reads kfac and assigns to ZV_kfac

if (IS_opt_phi==2) then
open(35,file=Qobsbarrec_file,status='old')
read(35,*) ZV_read_obs_tot
close(35)
call VecSetValues(ZV_Qobsbarrec,IS_obs_bas,IV_obs_loc1,                        &
                  ZV_read_obs_tot(IV_obs_index),INSERT_VALUES,ierr)
                  !here we only look at the observations within the basin
                  !studied
call VecAssemblyBegin(ZV_Qobsbarrec,ierr)
call VecAssemblyEnd(ZV_Qobsbarrec,ierr)  
!reads Qobsbarrec and assigns to ZV_Qobsbarrec
end if

!-------------------------------------------------------------------------------
!Set pnorm
!-------------------------------------------------------------------------------
call VecSetValues(ZV_pnorm,IS_one,IS_one-1,ZS_knorm_init,INSERT_VALUES,ierr)
call VecSetValues(ZV_pnorm,IS_one,IS_one,ZS_xnorm_init,INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_pnorm,ierr)
call VecAssemblyEnd(ZV_pnorm,ierr)
!set pnorm to pnorm=(knorm,xnorm)

!-------------------------------------------------------------------------------
!End of OPTION 2
!-------------------------------------------------------------------------------
end if

!*******************************************************************************
!Initialization procedure for OPTION 3
!*******************************************************************************
if (IS_opt_run==3) then

!-------------------------------------------------------------------------------
!copy main initial values into routing initial values 
!-------------------------------------------------------------------------------
call VecCopy(ZV_QoutinitM,ZV_QoutinitR,ierr)
call VecCopy(ZV_VinitM,ZV_VinitR,ierr)

!-------------------------------------------------------------------------------
!Read/set k and x
!-------------------------------------------------------------------------------
open(20,file=k_file,status='old')
read(20,*) ZV_read_riv_tot
call VecSetValues(ZV_k,IS_riv_bas,IV_riv_loc1,                                 &
                  ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_k,ierr)
call VecAssemblyEnd(ZV_k,ierr)
close(20)
!get values for k in a file and create the corresponding ZV_k vector

open(21,file=x_file,status='old')
read(21,*) ZV_read_riv_tot
call VecSetValues(ZV_x,IS_riv_bas,IV_riv_loc1,                                 &
                  ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_x,ierr)
call VecAssemblyEnd(ZV_x,ierr)
close(21)
!get values for x in a file and create the corresponding ZV_x vector

!-------------------------------------------------------------------------------
!Compute routing parameters and linear system matrix
!-------------------------------------------------------------------------------
call rapid_routing_param(ZV_k,ZV_x,ZV_C1,ZV_C2,ZV_C3,ZM_A)
!calculate Muskingum parameters and matrix ZM_A

call KSPSetOperators(ksp,ZM_A,ZM_A,ierr)
!Set KSP to use matrix ZM_A

!-------------------------------------------------------------------------------
!Calculate Muskingum matrix
!-------------------------------------------------------------------------------
call rapid_mus_mat

!-------------------------------------------------------------------------------
!Calculate observation operator
!-------------------------------------------------------------------------------
call rapid_runoff2streamflow_mat
!Create runoff-to-discharge operator (ZM_L)
call rapid_kf_obs_mat
!Create selection operator ZM_S
!Extract "observed" rows of ZM_L to build ZM_H = ZM_S*ZM_L
!Destroy ZM_L to free memory


!-------------------------------------------------------------------------------
!End of OPTION 3
!-------------------------------------------------------------------------------
end if

!*******************************************************************************
!Initialization procedure for OPTION 4
!*******************************************************************************
if (IS_opt_run==4) then

!-------------------------------------------------------------------------------
!copy main initial values into routing initial values 
!-------------------------------------------------------------------------------
call VecCopy(ZV_QoutinitM,ZV_QoutinitR,ierr)
call VecCopy(ZV_VinitM,ZV_VinitR,ierr)

!-------------------------------------------------------------------------------
!Read/set k and x
!-------------------------------------------------------------------------------
open(20,file=k_file,status='old')
read(20,*) ZV_read_riv_tot
call VecSetValues(ZV_k,IS_riv_bas,IV_riv_loc1,                                 &
                  ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_k,ierr)
call VecAssemblyEnd(ZV_k,ierr)
close(20)
!get values for k in a file and create the corresponding ZV_k vector

open(21,file=x_file,status='old')
read(21,*) ZV_read_riv_tot
call VecSetValues(ZV_x,IS_riv_bas,IV_riv_loc1,                                 &
                  ZV_read_riv_tot(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_x,ierr)
call VecAssemblyEnd(ZV_x,ierr)
close(21)
!get values for x in a file and create the corresponding ZV_x vector

!-------------------------------------------------------------------------------
!Compute routing parameters and linear system matrix
!-------------------------------------------------------------------------------
call rapid_routing_param(ZV_k,ZV_x,ZV_C1,ZV_C2,ZV_C3,ZM_A)
!calculate Muskingum parameters and matrix ZM_A

call KSPSetOperators(ksp,ZM_A,ZM_A,ierr)
!Set KSP to use matrix ZM_A

!-------------------------------------------------------------------------------
!Calculate observation operator
!-------------------------------------------------------------------------------
call rapid_run2strm_mat_smpl
!Create simplified runoff-to-discharge operator (ZM_L)

call rapid_kf_obs_mat
!Create selection operator ZM_S
!Extract "observed" rows of ZM_L to build ZM_H = ZM_S*ZM_L
!Destroy ZM_L to free memory

!-------------------------------------------------------------------------------
!End of OPTION 4
!-------------------------------------------------------------------------------
end if

!*******************************************************************************
!End subroutine
!*******************************************************************************
end subroutine rapid_init

#else

! Dummy version
subroutine rapid_init
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_init

#endif
