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
! 18 Nov 2021: Yeosang Yoon: Initial implementation in LIS based on the 
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
  real,   allocatable   :: rnfsto_lvec(:)
  integer               :: status
  logical               :: alarmCheck
  integer               :: c,r,t
  logical               :: dummy
  integer               :: ios, nid,qsid,qsbid

  !real,   allocatable   :: rnfsto_mm(:,:,:),bsfsto_mm(:,:,:)
! _______________________________________________

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "RAPID router model alarm")
  if(alarmCheck) then 
     
    ! if(RAPID_routing_struc(n)%useens.eq.0) then 
    !   allocate(rnfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),1))
    !   allocate(bsfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),1))
    ! elseif(RAPID_routing_struc(n)%useens.eq.1) then
    !   allocate(rnfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n)))
    !   allocate(bsfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n)))
    ! endif

     if(RAPID_routing_struc(n)%useens.eq.1) then 
        if(LIS_masterproc) then     
           !run the routing model at 1 hour output interval
           
           call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",&
                sf_runoff_field,&
                rc=status)
           call LIS_verify(status, "ESMF_StateGet failed for Surface Runoff")
           
           call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
                baseflow_field, rc=status)
           call LIS_verify(status, "ESMF_StateGet failed for Subsurface Runoff")
           
           call ESMF_FieldGet(sf_runoff_field,localDE=0,&
                farrayPtr=surface_runoff_t,&
                rc=status)
           call LIS_verify(status, "ESMF_FieldGet failed for Surface Runoff")
           
           call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow_t,&
                rc=status)
           call LIS_verify(status, "ESMF_FieldGet failed for Subsurface Runoff")
           
           allocate(surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n)))
           allocate(baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n)))           

           do m=1,LIS_rc%nensem(n)
              
              surface_runoff = 0.0
              baseflow = 0.0

              call LIS_tile2grid(n,m,surface_runoff,surface_runoff_t,1)
              call LIS_tile2grid(n,m,baseflow,baseflow_t,1)
             
              ! call rapid ... 
              
           enddo
           deallocate(surface_runoff)
           deallocate(baseflow)
        endif

           
          !call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
          !     value=rnfsto_lvec(t),vlevel=1,unit="mm",&
          !     direction="-")
           
     else
        if(LIS_masterproc) then     
           if(LIS_rc%lsm.ne."none") then 
              
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
             
              print *, LIS_rc%gnc(n), LIS_rc%gnr(n) 
              print *, surface_runoff
              stop
         
           else !read from previous output. 
              
              allocate(surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n)))
              allocate(baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n)))

              !call readrunoffdata(trim(LIS_rc%runoffdatasource)//char(0),&
              !     n,surface_runoff, baseflow)
              
           endif
           

           !call model(LIS_rc%udef,&

        
           deallocate(surface_runoff)
           deallocate(baseflow)
        endif

        !call LIS_grid2patch(n,LIS_rc%lsm_index,rnfsto_mm(:,:,1),rnfsto_lvec,dummy)
      

        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
          !call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
          !     value=rnfsto_lvec(t),vlevel=1,unit="mm",&
          !     direction="-")
          
        enddo
                
     endif
     !deallocate(rnfsto_lvec)

     !if(LIS_masterproc)then
       !if(allocated(rnfsto_mm))deallocate(rnfsto_mm)
     !endif
     
  endif
end subroutine RAPID_routing_run
