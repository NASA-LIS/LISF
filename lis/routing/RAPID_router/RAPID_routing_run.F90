!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! 
! !DESCRIPTION: 
! 
!  Reference: 
!
! !REVISION HISTORY: 
! 18 Mar 2021: Yeosang Yoon: Initial implementation in LIS based on the 
!                            RAPID offline routing code. 
! 
! !USES: 
subroutine RAPID_routing_run(n)

  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_routingMod
  use LIS_logMod
  use LIS_historyMod
  use LIS_histDataMod
  use LIS_constantsMod
  use RAPID_routingMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  implicit none
  
  integer, intent(in)   :: n  
  integer               :: m
  type(ESMF_Field)      :: sf_runoff_field
  type(ESMF_Field)      :: baseflow_field
  real,   allocatable   :: surface_runoff(:,:)
  real,   pointer       :: surface_runoff_t(:)
  real,   pointer       :: baseflow_t(:)
  real,   allocatable   :: qs(:,:)
  real,   allocatable   :: qsb(:,:)
  real,   allocatable   :: baseflow(:,:)
  real,   allocatable   :: streamflow_lvec(:)
  integer               :: status
  logical               :: alarmCheck
  integer               :: c,r,t
  logical               :: dummy
  integer               :: ios, nid,qsid,qsbid
! _______________________________________________

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "RAPID router model alarm")
  if(alarmCheck) then 
     allocate(streamflow_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))

     streamflow_lvec = LIS_rc%udef

     if(LIS_masterproc) then     
        call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",sf_runoff_field,&
             rc=status)
        call LIS_verify(status, "ESMF_StateGet failed for Surface Runoff")
              
        call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
             baseflow_field, rc=status)
        call LIS_verify(status, "ESMF_StateGet failed for Subsurface Runoff")
              
        call ESMF_FieldGet(sf_runoff_field,localDE=0,farrayPtr=surface_runoff_t,&
             rc=status)
        call LIS_verify(status, "ESMF_FieldGet failed for Surface Runoff")
              
        call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow_t,&
             rc=status)
        call LIS_verify(status, "ESMF_FieldGet failed for Subsurface Runoff")
                         
        allocate(surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        allocate(baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n)))
              
        call LIS_tile2grid(n,surface_runoff,surface_runoff_t,1)
        call LIS_tile2grid(n,baseflow,baseflow_t,1)
        
        !TODO               
        ! run RAPID
        call RAPID_model_main (RAPID_routing_struc(n)%run_opt,RAPID_routing_struc(n)%routing_opt,RAPID_routing_struc(n)%phi_opt,&
                               RAPID_routing_struc(n)%connectfile,RAPID_routing_struc(n)%max_reach,RAPID_routing_struc(n)%weightfile,&
                               RAPID_routing_struc(n)%basinIDfile,RAPID_routing_struc(n)%kfile,RAPID_routing_struc(n)%xfile,&
                               RAPID_routing_struc(n)%nmlfile)
        
        deallocate(surface_runoff)
        deallocate(baseflow)
     endif

     !TODO: need to change streamflow format
     call LIS_grid2patch(n,LIS_rc%lsm_index, RAPID_routing_struc(n)%streamflow,&
             streamflow_lvec)
     !call LIS_grid2tile
     !call LIS_grid2patch(n,LIS_rc%lsm_index,rnfsto_mm(:,:,1),rnfsto_lvec,dummy)
      

     do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_STREAMFLOW,&
             value=streamflow_lvec(t),vlevel=1,unit="m3 s-1",&
             direction="-")
       !call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
       !     value=rnfsto_lvec(t),vlevel=1,unit="mm",&
       !     direction="-")
          
     enddo
                
     deallocate(streamflow_lvec)
     !deallocate(rnfsto_lvec)

     !if(LIS_masterproc)then
       !if(allocated(rnfsto_mm))deallocate(rnfsto_mm)
     !endif
     
  endif
end subroutine RAPID_routing_run
