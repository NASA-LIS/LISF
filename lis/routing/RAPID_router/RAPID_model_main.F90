!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!BOP
!
! !DESCRIPTION:
!
!  Reference:
!
! !REVISION HISTORY:
! 26 Mar 2021: Yeosang Yoon: Initial implementation in LIS based on the
!                            RAPID offline routing code (rapid_main.F90). 


!*******************************************************************************
! Subroutine- RAPID_model_main (rapid_main)
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine RAPID_model_main (n,bQinit,bQfinal,bV,bhum,bfor,   &
                             bdam,binfluence,buq,             &
                             run_opt,routing_opt,phi_opt,     &
                             connectfile,max_reach,n_riv_tot, &
                             weightfile,n_wei_table,          &
                             basinIDfile,n_riv_bas,           &
                             kfile,xfile,                     &
                             nmlfile,qfile,                   &
                             nc,nr,runsf,runsb,initCheck,     &
                             dt,routingInterval)

!Purpose:
!Allows to route water through a river network, and to estimate optimal 
!parameters using the inverse method 
!Author: 
!Cedric H. David, 2008-2020.

!Fortran includes, modules, and implicity
#include <petsc/finclude/petsctao.h>
use petsctao
use rapid_var, only :                                                          &    !Yeosang Yoon
                   namelist_file,                                              &
                   Vlat_file,Qfor_file,Qhum_file,                              &
                   Qout_file,V_file,                                           &
                   IS_M,JS_M,JS_RpM,IS_RpM,IS_RpF,IS_RpH,                      &
                   ZS_TauM,ZS_dtM,ZS_TauR,ZS_dtR,                              &
                   ZV_pnorm,                                                   &
                   ZV_k,ZV_x,ZV_C1,ZV_C2,ZV_C3,                                &
                   ZV_Qext,ZV_Qfor,ZV_Qlat,ZV_Qhum,ZV_Qdam,                    &
                   ZV_Vlat,                                                    &
                   ZV_QoutR,ZV_QoutinitR,ZV_QoutbarR,                          &
                   ZV_VR,ZV_VinitR,ZV_VbarR,                                   &
                   ierr,rank,stage,temp_char,temp_char2,                       &
                   ZS_one,                                                     &
                   IS_riv_tot,IS_riv_bas,IS_for_bas,IV_riv_bas_id,             &
                   IS_dam_bas,IS_hum_bas,                                      &
                   ZS_time1,ZS_time2,ZS_time3,ZS_time4,                        &
                   IV_nc_start,IV_nc_count,IV_nc_count2,                       &
                   BS_opt_Qinit,BS_opt_Qfinal,BS_opt_V,BS_opt_hum,             &
                   BS_opt_for,BS_opt_dam,BS_opt_influence,BS_opt_uq,           &
                   IS_opt_run,IS_opt_routing,IS_opt_phi,                       &
                   tao,                                                        &
                   Qobs_file,ZV_Qobs,                                          &
                   ZV_Qbmean,ZV_dQeb,ZS_val,                                   &
                   ZV_QoutinitR_save,                                          &
                   rapid_connect_file,IS_max_up,weight_table_file,             &
                   n_weight_table,riv_bas_id_file,k_file,x_file,               &
                   vecscat,ZV_SeqZero,ZV_pointer,IV_riv_loc1,IV_riv_index,     &
                   ZV_riv_tot_lat,ZV_riv_tot_lon       

use LIS_coreMod, only: LIS_rc
use LIS_logMod
use LIS_timeMgrMod
use RAPID_routingMod, only : RAPID_routing_struc

implicit none
external rapid_phiroutine
!because the subroutine is called by a function

!Yeosang Yoon
!Arguments
integer,       intent(in)     :: n
logical,       intent(in)     :: bQinit       ! initial flow
logical,       intent(in)     :: bQfinal      ! write final flow
logical,       intent(in)     :: bV           ! compute volume
logical,       intent(in)     :: bhum         ! human-induced flows
logical,       intent(in)     :: bfor         ! forcing
logical,       intent(in)     :: bdam         ! dam model used
logical,       intent(in)     :: binfluence   ! output influence
logical,       intent(in)     :: buq          ! uncertainty quantif.

integer,       intent(in)     :: run_opt      ! run option
integer,       intent(in)     :: routing_opt  ! routing option
integer,       intent(in)     :: phi_opt      ! phi option
character*200, intent(in)     :: connectfile  ! river connectivity file
integer,       intent(in)     :: n_riv_tot    ! number of river connectivity
integer,       intent(in)     :: max_reach    ! max number of upstream reaches
character*200, intent(in)     :: weightfile   ! river weight table
integer,       intent(in)     :: n_wei_table  ! number of reach in weight table file
character*200, intent(in)     :: basinIDfile  ! river basin ID file
integer,       intent(in)     :: n_riv_bas    ! number of river basins
character*200, intent(in)     :: kfile        ! Muskingum parameter k file
character*200, intent(in)     :: xfile        ! Muskingum parameter x file

character*200, intent(in)     :: nmlfile
character*200, intent(in)     :: qfile

!Runoff data are in kg/m2 accumulated over a time step
integer,       intent(in)     :: nc, nr
real,          intent(in)     :: runsf(nc,nr)        ! surface runoff
real,          intent(in)     :: runsb(nc,nr)        ! subsurface runoff

logical,       intent(inout)  :: initCheck
real,          intent(in)     :: dt                  ! internal time step (in seconds) 
real,          intent(in)     :: routingInterval     ! routing time step (in seconds)
logical                       :: alarmCheck
PetscScalar,   allocatable    :: Qinit(:)
 
!*******************************************************************************
!Initialize
!*******************************************************************************
!Yeosang Yoon, for one-time initialization
if (initCheck .eqv. .true.) then
   BS_opt_Qinit=bQinit
   BS_opt_Qfinal=bQfinal
   BS_opt_V=bV
   BS_opt_hum=bhum
   BS_opt_for=bfor
   BS_opt_dam=bdam
   BS_opt_influence=binfluence
   BS_opt_uq=buq
   
   IS_opt_run=run_opt
   IS_opt_routing=routing_opt
   IS_opt_phi=phi_opt
   rapid_connect_file=trim(connectfile)
   IS_max_up=max_reach
   weight_table_file=trim(weightfile)
   riv_bas_id_file=trim(basinIDfile)
   k_file=kfile
   x_file=xfile

   namelist_file=nmlfile

   ! size of domain
   IS_riv_tot=n_riv_tot
   IS_riv_bas=n_riv_bas
   n_weight_table=n_wei_table
   
   ! temporal parameters
   ZS_TauM=int(dt)              ! duration of main procedure; same as ZS_dtM for LIS-RAPID       
   ZS_dtM=int(dt)               ! time step of main procedure
   ZS_TauR=int(dt)              ! duration of river routing procedure; same as ZS_dtM
   ZS_dtR=int(routingInterval)  ! time step of river routing procedure   
   
   call rapid_init

   ! for RAPID output
   if (rank==0) then 
       RAPID_routing_struc(n)%riv_bas_id=IV_riv_bas_id
   endif

   initCheck = .false.

   ! for RAPID restart
   if(RAPID_routing_struc(n)%startmode.eq."restart") then
      if (rank==0) then
          allocate(Qinit(IS_riv_bas))
          Qinit=RAPID_routing_struc(n)%Qout
          
          call VecSetValues(ZV_QoutinitR,IS_riv_bas,IV_riv_loc1, &
                            Qinit(IV_riv_index),INSERT_VALUES,ierr)
          deallocate(Qinit)
      endif
      call VecAssemblyBegin(ZV_QoutinitR,ierr)
      call VecAssemblyEnd(ZV_QoutinitR,ierr)
   endif
endif

!Qout_file=trim(qfile)   ! LIS-RAPID output filename
!alarmCheck = LIS_isAlarmRinging(LIS_rc,"RAPID router output alarm")

!*******************************************************************************
!OPTION 1 - use to calculate flows and volumes and generate output data
!*******************************************************************************
if (IS_opt_run==1) then

!-------------------------------------------------------------------------------
!Open Vlat file
!-------------------------------------------------------------------------------
!call rapid_open_Vlat_file(Vlat_file)    ! Yeosang Yoon, not use these subroutines within LIS
!call rapid_meta_Vlat_file(Vlat_file)

!-------------------------------------------------------------------------------
!Quantify uncertainty
!-------------------------------------------------------------------------------
!if (BS_opt_uq) call rapid_uq

!-------------------------------------------------------------------------------
!Create and open Qout file; Yeosang Yoon; not use
!-------------------------------------------------------------------------------
!if(alarmCheck) call rapid_create_Qout_file(Qout_file)
!if(alarmCheck) call rapid_open_Qout_file(Qout_file)

!-------------------------------------------------------------------------------
!Create and open V_file
!-------------------------------------------------------------------------------
!if(alarmCheck) then
if (BS_opt_V) call rapid_create_V_file(V_file)
if (BS_opt_V) call rapid_open_V_file(V_file)
!!endif

!-------------------------------------------------------------------------------
!Open remaining files
!-------------------------------------------------------------------------------
if (BS_opt_for) call rapid_open_Qfor_file(Qfor_file)
if (BS_opt_hum) call rapid_open_Qhum_file(Qhum_file)

!-------------------------------------------------------------------------------
!Make sure the vectors potentially used for inflow to dams are initially null
!-------------------------------------------------------------------------------
call VecSet(ZV_Qext,0*ZS_one,ierr)                         !Qext=0
call VecSet(ZV_QoutbarR,0*ZS_one,ierr)                     !QoutbarR=0
!This should be done by PETSc but just to be safe

!-------------------------------------------------------------------------------
!Set initial value of Qext from Qout_dam0
!-------------------------------------------------------------------------------
if (BS_opt_dam .and. IS_dam_bas>0) then
     call rapid_set_Qext0                                  !Qext from Qout_dam0
     !call VecView(ZV_Qext,PETSC_VIEWER_STDOUT_WORLD,ierr)
end if

!-------------------------------------------------------------------------------
!Read, compute and write
!-------------------------------------------------------------------------------
!call PetscLogStageRegister('Read Comp Write',stage,ierr)
!call PetscLogStagePush(stage,ierr)
ZS_time3=0

IV_nc_start=(/1,1/)
IV_nc_count=(/IS_riv_tot,1/)
IV_nc_count2=(/IS_riv_bas,1/)

do JS_M=1,IS_M                   ! Yeosang Yoon, both IS_M and IS_RpM are 1 for LIS-RAPID
do JS_RpM=1,IS_RpM

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Read/set surface and subsurface volumes
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!call rapid_read_Vlat_file                    ! Yeosang Yoon, not use this subroutine within LIS
call rapid_Vlat(nc,nr,runsf,runsb)            ! Yeosang Yoon, read LSM runoff and transfer to boundary inflows

call VecCopy(ZV_Vlat,ZV_Qlat,ierr)            !Qlat=Vlat
call VecScale(ZV_Qlat,1/ZS_TauR,ierr)         !Qlat=Qlat/TauR

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Read/set upstream forcing
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (BS_opt_for .and. IS_for_bas>0                                              &
                   .and. mod((JS_M-1)*IS_RpM+JS_RpM,IS_RpF)==1) then

call rapid_read_Qfor_file

end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Run dam model based on previous values of QoutbarR and Qext to get Qdam
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (BS_opt_dam .and. IS_dam_bas>0) then

call rapid_get_Qdam

end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Read/set human induced flows
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (BS_opt_hum .and. IS_hum_bas>0                                              &
                   .and. mod((JS_M-1)*IS_RpM+JS_RpM,IS_RpH)==1) then

call rapid_read_Qhum_file

end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!calculation of Qext
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecCopy(ZV_Qlat,ZV_Qext,ierr)                            !Qext=Qlat
if (BS_opt_for) call VecAXPY(ZV_Qext,ZS_one,ZV_Qfor,ierr)     !Qext=Qext+1*Qfor
if (BS_opt_dam) call VecAXPY(ZV_Qext,ZS_one,ZV_Qdam,ierr)     !Qext=Qext+1*Qdam
if (BS_opt_hum) call VecAXPY(ZV_Qext,ZS_one,ZV_Qhum,ierr)     !Qext=Qext+1*Qhum

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Routing procedure
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call PetscTime(ZS_time1,ierr)

call rapid_routing(ZV_C1,ZV_C2,ZV_C3,ZV_Qext,                                  &
                   ZV_QoutinitR,                                               &
                   ZV_QoutR,ZV_QoutbarR)

if (BS_opt_V) call rapid_QtoV(ZV_k,ZV_x,ZV_QoutbarR,ZV_Qext,ZV_VbarR)

call PetscTime(ZS_time2,ierr)
ZS_time3=ZS_time3+ZS_time2-ZS_time1

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Update variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecCopy(ZV_QoutR,ZV_QoutinitR,ierr)
call VecCopy(ZV_VR,ZV_VinitR,ierr)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!write outputs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!if(alarmCheck) then
!   write(LIS_logunit,*) '[INFO] Writing routing model output to:'
!   write(LIS_logunit,*) trim(Qout_file)
!
!   call rapid_write_Qout_file
if (BS_opt_V) call rapid_write_V_file
!endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Yeosang Yoon, for RAPID output/restart file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecScatterBegin(vecscat,ZV_QoutR,ZV_SeqZero,                              &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat,ZV_QoutR,ZV_SeqZero,                                &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
if(rank==0) then
   RAPID_routing_struc(n)%riv_tot_lat=ZV_riv_tot_lat
   RAPID_routing_struc(n)%riv_tot_lon=ZV_riv_tot_lon 

   call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)
   RAPID_routing_struc(n)%Qout=ZV_pointer
   call VecRestoreArrayF90(ZV_SeqZero,ZV_pointer,ierr)
endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Update netCDF location
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!if (rank==0) IV_nc_start(2)=IV_nc_start(2)+1
!do not comment out if writing directly from the routing subroutine


end do
end do

!-------------------------------------------------------------------------------
!Performance statistics
!-------------------------------------------------------------------------------
!call PetscPrintf(PETSC_COMM_WORLD,'Cumulative time for routing only'           &
!                                  //char(10),ierr)
!write(temp_char ,'(i10)')   rank
!write(temp_char2,'(f10.2)') ZS_time3
!call PetscSynchronizedPrintf(PETSC_COMM_WORLD,'Rank     :'//temp_char //', '// &
!                                              'Time     :'//temp_char2//       &
!                                               char(10),ierr)
!call PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_NULL_INTEGER,ierr)
!
!call PetscLogStagePop(ierr)
!call PetscPrintf(PETSC_COMM_WORLD,'Output data created'//char(10),ierr)


!-------------------------------------------------------------------------------
!Close files
!-------------------------------------------------------------------------------
!call rapid_close_Qout_file
!call rapid_close_Vlat_file
if (BS_opt_for) call rapid_close_Qfor_file(Qfor_file)
if (BS_opt_hum) call rapid_close_Qhum_file(Qhum_file)
if (BS_opt_V) call rapid_close_V_file(V_file)


!-------------------------------------------------------------------------------
!End of OPTION 1
!-------------------------------------------------------------------------------
end if

#if 0
!*******************************************************************************
!OPTION 2 - Optimization
!*******************************************************************************
if (IS_opt_run==2) then

!-------------------------------------------------------------------------------
!Only one computation of phi - For testing purposes only
!-------------------------------------------------------------------------------
!call PetscLogStageRegister('One comp of phi',stage,ierr)
!call PetscLogStagePush(stage,ierr)
!!do JS_M=1,5
!call rapid_phiroutine(tao,ZV_pnorm,ZS_phi,PETSC_NULL,ierr)
!!enddo
!call PetscLogStagePop(ierr)

!-------------------------------------------------------------------------------
!Optimization procedure
!-------------------------------------------------------------------------------
call PetscLogStageRegister('Optimization   ',stage,ierr)
call PetscLogStagePush(stage,ierr)
call TaoSetObjectiveRoutine(tao,rapid_phiroutine,PETSC_NULL_INTEGER,ierr)
call TaoSetInitialVector(tao,ZV_pnorm,ierr)
call TaoSolve(tao,ierr)

call TaoView(tao,PETSC_VIEWER_STDOUT_WORLD,ierr)
call PetscPrintf(PETSC_COMM_WORLD,'final normalized p=(k,x)'//char(10),ierr)
call VecView(ZV_pnorm,PETSC_VIEWER_STDOUT_WORLD,ierr)
call PetscLogStagePop(ierr)

!-------------------------------------------------------------------------------
!End of OPTION 2
!-------------------------------------------------------------------------------
end if


!*******************************************************************************
!OPTION 3/4 - data assimilation (of discharge to correct runoff input)
!*******************************************************************************
if ((IS_opt_run==3).or.(IS_opt_run==4)) then

!-------------------------------------------------------------------------------
!Open Vlat file
!-------------------------------------------------------------------------------
call rapid_open_Vlat_file(Vlat_file)
call rapid_meta_Vlat_file(Vlat_file)

!-------------------------------------------------------------------------------
!Compute runoff error covariance matrix
!-------------------------------------------------------------------------------
call rapid_cov_mat
call rapid_kf_cov_mat

!-------------------------------------------------------------------------------
!Open observation file
!-------------------------------------------------------------------------------
call rapid_open_Qobs_file(Qobs_file)

!-------------------------------------------------------------------------------
!Create and open Qout file
!-------------------------------------------------------------------------------
call rapid_create_Qout_file(Qout_file)
call rapid_open_Qout_file(Qout_file)

!-------------------------------------------------------------------------------
!Create and open V_file
!-------------------------------------------------------------------------------
if (BS_opt_V) call rapid_create_V_file(V_file)
if (BS_opt_V) call rapid_open_V_file(V_file)

!-------------------------------------------------------------------------------
!Make sure the vectors potentially used for inflow to dams are initially null
!-------------------------------------------------------------------------------
call VecSet(ZV_Qext,0*ZS_one,ierr)                         !Qext=0
call VecSet(ZV_QoutbarR,0*ZS_one,ierr)                     !QoutbarR=0
!This should be done by PETSc but just to be safe

!-------------------------------------------------------------------------------
!Read, compute, Kalman, and write
!-------------------------------------------------------------------------------
call PetscLogStageRegister('Read Comp KF Write',stage,ierr)
call PetscLogStagePush(stage,ierr)
ZS_time3=0
ZS_time4=0

IV_nc_start=(/1,1/)
IV_nc_count=(/IS_riv_tot,1/)
IV_nc_count2=(/IS_riv_bas,1/)

do JS_M=1,IS_M

!write(temp_char,'(i10)') IS_M
!write(temp_char2,'(i10)') JS_M
!call PetscPrintf(PETSC_COMM_WORLD,'Assimilation day '                          &
!                                  //temp_char2//'/'//temp_char//               &
!                                  char(10),ierr)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Initialize ZV_Qbmean and ZV_dQeb for assimilation and save initial condition
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call VecSet(ZV_Qbmean,0*ZS_one,ierr)    !Daily-averaged simulated discharge
call VecSet(ZV_dQeb,0*ZS_one,ierr)      !Kalman filter correction
call VecSet(ZV_Qobs,0*ZS_one,ierr)      !Observation vector (size IS_riv_bas)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Run RAPID forecast
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Save initial condition for later (analysis) run
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call VecCopy(ZV_QoutinitR,ZV_QoutinitR_save,ierr)

do JS_RpM=1,IS_RpM

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Read/set surface and subsurface volumes
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call rapid_read_Vlat_file

call VecCopy(ZV_Vlat,ZV_Qlat,ierr)            !Qlat=Vlat
call VecScale(ZV_Qlat,1/ZS_TauR,ierr)         !Qlat=Qlat/TauR

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!calculation of Qext
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call VecCopy(ZV_Qlat,ZV_Qext,ierr)                            !Qext=Qlat

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Routing procedure with background runoff
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call PetscTime(ZS_time1,ierr)

call rapid_routing(ZV_C1,ZV_C2,ZV_C3,ZV_Qext,                                  &
                   ZV_QoutinitR,                                               &
                   ZV_QoutR,ZV_QoutbarR)

if (BS_opt_V) call rapid_QtoV(ZV_k,ZV_x,ZV_QoutbarR,ZV_Qext,ZV_VbarR)

call PetscTime(ZS_time2,ierr)
ZS_time3=ZS_time3+ZS_time2-ZS_time1

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Update ZV_Qbmean
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
ZS_val= ZS_one/(REAL(IS_RpM))
call VecAXPY(ZV_Qbmean,ZS_val,ZV_QoutbarR,ierr)

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Update variables
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call VecCopy(ZV_QoutR,ZV_QoutinitR,ierr)
call VecCopy(ZV_VR,ZV_VinitR,ierr)

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Update netCDF location
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
if (rank==0) IV_nc_start(2)=IV_nc_start(2)+1
!do not comment out if writing directly from the routing subroutine

end do

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Kalman filtering
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call PetscTime(ZS_time1,ierr)

call rapid_read_Qobs_file
!Build observation vector

call rapid_kf_update
!Kalman filter update

call PetscTime(ZS_time2,ierr)
ZS_time4=ZS_time4+ZS_time2-ZS_time1

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!Run RAPID analysis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Re-set initial condition and netcdf location
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call VecCopy(ZV_QoutinitR_save,ZV_QoutinitR,ierr)
if (rank==0) IV_nc_start(2)=IV_nc_start(2)-IS_RpM

do JS_RpM=1,IS_RpM

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Read/set surface and subsurface volumes
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call rapid_read_Vlat_file

call VecCopy(ZV_Vlat,ZV_Qlat,ierr)            !Qlat=Vlat
call VecScale(ZV_Qlat,1/ZS_TauR,ierr)         !Qlat=Qlat/TauR

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Update surface and subsurface flow with Kalman Filter correction
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call VecAXPY(ZV_Qlat,ZS_one,ZV_dQeb,ierr)         !Qlat=Qlat+dQeb

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!calculation of Qext
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call VecCopy(ZV_Qlat,ZV_Qext,ierr)                            !Qext=Qlat+dQeb

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Routing procedure with analysis runoff
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call PetscTime(ZS_time1,ierr)

call rapid_routing(ZV_C1,ZV_C2,ZV_C3,ZV_Qext,                                  &
                   ZV_QoutinitR,                                               &
                   ZV_QoutR,ZV_QoutbarR)

if (BS_opt_V) call rapid_QtoV(ZV_k,ZV_x,ZV_QoutbarR,ZV_Qext,ZV_VbarR)

call PetscTime(ZS_time2,ierr)
ZS_time3=ZS_time3+ZS_time2-ZS_time1

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Update variables
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call VecCopy(ZV_QoutR,ZV_QoutinitR,ierr)
call VecCopy(ZV_VR,ZV_VinitR,ierr)

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!write outputs
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
call rapid_write_Qout_file
if (BS_opt_V) call rapid_write_V_file

!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!Update netCDF location
!  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
if (rank==0) IV_nc_start(2)=IV_nc_start(2)+1
!do not comment out if writing directly from the routing subroutine


end do

end do

!-------------------------------------------------------------------------------
!Performance statistics
!-------------------------------------------------------------------------------
call PetscPrintf(PETSC_COMM_WORLD,'Cumulative time for routing (background '// &
                                  'and analysis)'//char(10),ierr)
write(temp_char ,'(i10)')   rank
write(temp_char2,'(f10.2)') ZS_time3
call PetscSynchronizedPrintf(PETSC_COMM_WORLD,'Rank     :'//temp_char //', '// &
                                              'Time     :'//temp_char2//       &
                                               char(10),ierr)
call PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_NULL_INTEGER,ierr)

call PetscPrintf(PETSC_COMM_WORLD,'Cumulative time for Kalman filtering'       &
                                  //char(10),ierr)
write(temp_char ,'(i10)')   rank
write(temp_char2,'(f10.2)') ZS_time4
call PetscSynchronizedPrintf(PETSC_COMM_WORLD,'Rank     :'//temp_char //', '// &
                                              'Time     :'//temp_char2//       &
                                               char(10),ierr)
call PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_NULL_INTEGER,ierr)

call PetscLogStagePop(ierr)
call PetscPrintf(PETSC_COMM_WORLD,'Output data created'//char(10),ierr)


!-------------------------------------------------------------------------------
!Close files
!-------------------------------------------------------------------------------
call rapid_close_Qout_file
call rapid_close_Vlat_file
call rapid_close_Qobs_file
if (BS_opt_V) call rapid_close_V_file(V_file)

!-------------------------------------------------------------------------------
!End of OPTION 3/4
!-------------------------------------------------------------------------------
end if

#endif

!*******************************************************************************
!Finalize
!*******************************************************************************
!call rapid_clean_var
if (LIS_rc%endtime==1) then
    call rapid_final
endif

end subroutine RAPID_model_main

#else

! Dummy version
subroutine RAPID_model_main
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine RAPID_model_main

#endif
