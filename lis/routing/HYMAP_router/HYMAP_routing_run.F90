!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! 
! !DESCRIPTION: 
! 
!  Reference: Augusto C.V. Getirana, A. Boone, D. Yamazaki, B. Decharme,
!             F. Papa, and N. Mognard, 2012: The Hydrological Modeling
!             and Analysis Platform (HyMAP): Evaluation in the Amazon
!             Basin.  Journal of Hydrometeorology, 13, 1641-1665.
!             doi: http://dx.doi.org/10.1175/JHM-D-12-021.1
!
! !REVISION HISTORY: 
! 08 Nov 2011: Augusto Getirana, Initial implementation in LIS based on the 
!                                HYMAP offline routing code. 
! 
! !USES: 
subroutine HYMAP_routing_run(n)

  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_routingMod
  use LIS_logMod
  use LIS_historyMod
  use LIS_histDataMod
  use LIS_constantsMod
  use HYMAP_routingMod
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
  real,   allocatable   :: rivsto_lvec(:)
  real,   allocatable   :: rivdph_lvec(:)
  real,   allocatable   :: rivvel_lvec(:)
  real,   allocatable   :: rivout_lvec(:)
  real,   allocatable   :: evpout_lvec(:)
  real,   allocatable   :: fldout_lvec(:)
  real,   allocatable   :: fldsto_lvec(:)
  real,   allocatable   :: flddph_lvec(:)
  real,   allocatable   :: fldvel_lvec(:)
  real,   allocatable   :: fldfrc_lvec(:)
  real,   allocatable   :: fldare_lvec(:)
  real,   allocatable   :: sfcelv_lvec(:)
  real,   allocatable   :: rnfsto_lvec(:)
  real,   allocatable   :: bsfsto_lvec(:)
  integer               :: status
  logical               :: alarmCheck
  integer               :: c,r,t
  logical               :: dummy
  integer               :: ios, nid,qsid,qsbid

  real,   allocatable   :: rnfsto_mm(:,:,:),bsfsto_mm(:,:,:)
! _______________________________________________

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "HYMAP router model alarm")
  if(alarmCheck) then 

     allocate(rivsto_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(rivdph_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(rivvel_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(rivout_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(evpout_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(fldout_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(fldsto_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(flddph_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(fldvel_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(fldfrc_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(fldare_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(sfcelv_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(rnfsto_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(bsfsto_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     
     if(HYMAP_routing_struc(n)%useens.eq.0) then 
       allocate(rnfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),1))
       allocate(bsfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),1))
     elseif(HYMAP_routing_struc(n)%useens.eq.1) then
       allocate(rnfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n)))
       allocate(bsfsto_mm(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n)))
     endif

     if(HYMAP_routing_struc(n)%useens.eq.1) then 
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

              call LIS_tile2grid(n,m,surface_runoff,surface_runoff_t)
              call LIS_tile2grid(n,m,baseflow,baseflow_t)
              
              call model(LIS_rc%udef,&
                   LIS_rc%gnc(n),&
                   LIS_rc%gnr(n),&
                   HYMAP_routing_struc(n)%inz,&
                   nint(HYMAP_routing_struc(n)%dt),&
                   HYMAP_routing_struc(n)%flowtype,&
                   HYMAP_routing_struc(n)%linresflag,&
                   HYMAP_routing_struc(n)%evapflag,&
                   HYMAP_routing_struc(n)%seqx,&
                   HYMAP_routing_struc(n)%seqy,&
                   HYMAP_routing_struc(n)%nseqriv,&
                   HYMAP_routing_struc(n)%nseqall,&
                   HYMAP_routing_struc(n)%nextx,&
                   HYMAP_routing_struc(n)%nexty,&
                   HYMAP_routing_struc(n)%elevtn,&
                   HYMAP_routing_struc(n)%nxtdst,&
                   HYMAP_routing_struc(n)%grarea,&
                   HYMAP_routing_struc(n)%fldgrd,&
                   HYMAP_routing_struc(n)%fldman,&
                   HYMAP_routing_struc(n)%fldstomax,&
                   HYMAP_routing_struc(n)%rivman,&
                   HYMAP_routing_struc(n)%rivelv,&
                   HYMAP_routing_struc(n)%rivstomax,&
                   HYMAP_routing_struc(n)%rivlen,&
                   HYMAP_routing_struc(n)%rivwth,&
                   HYMAP_routing_struc(n)%rivhgt,&
                   HYMAP_routing_struc(n)%rivare,&
                   HYMAP_routing_struc(n)%rslpmin,&
                   HYMAP_routing_struc(n)%trnoff,&
                   HYMAP_routing_struc(n)%tbsflw,&
                   HYMAP_routing_struc(n)%cntime,&
                   HYMAP_routing_struc(n)%mask,&
                   surface_runoff,&
                   baseflow,&
                   
                   !==== TEMPORARILY COMMENTED ====
                   !input evap        
                   !HYMAP_routing_struc(n)%tsfc,&
                   !HYMAP_routing_struc(n)%psur,&
                   !HYMAP_routing_struc(n)%pres,&
                   !HYMAP_routing_struc(n)%tair,&
                   !HYMAP_routing_struc(n)%relh,&
                   !HYMAP_routing_struc(n)%wind,&
                   !HYMAP_routing_struc(n)%zref,&
                   
                   !outputs
                   HYMAP_routing_struc(n)%rivsto(:,:,m),&
                   HYMAP_routing_struc(n)%rivdph(:,:,m),&
                   HYMAP_routing_struc(n)%rivvel(:,:,m),&
                   HYMAP_routing_struc(n)%rivout(:,:,m),&
                   HYMAP_routing_struc(n)%evpout(:,:,m),&
                   HYMAP_routing_struc(n)%fldout(:,:,m),&
                   HYMAP_routing_struc(n)%fldsto(:,:,m),&
                   HYMAP_routing_struc(n)%flddph(:,:,m),&
                   HYMAP_routing_struc(n)%fldvel(:,:,m),&
                   HYMAP_routing_struc(n)%fldfrc(:,:,m),&
                   HYMAP_routing_struc(n)%fldare(:,:,m),&
                   HYMAP_routing_struc(n)%sfcelv(:,:,m),&
                   HYMAP_routing_struc(n)%rnfsto(:,:,m),&
                   HYMAP_routing_struc(n)%bsfsto(:,:,m))

!========================================================================			
! the following is done to distribute the global arrays to the individual 
! processors, so that they can use the generic LIS interfaces for writing
! output -- which works based on the assumption of multiprocessors. 
! When the HYMAP routine itself is parallelized, this can (should) be taken
! out. 
             !ag (26Oct2017) - converting surface runoff and baseflow storage units from m3 to mm
             do c=1,LIS_rc%gnc(n)
               do r=1,LIS_rc%gnr(n)
                 if(HYMAP_routing_struc(n)%rnfsto(c,r,m)/=LIS_rc%udef.and.&
                    HYMAP_routing_struc(n)%bsfsto(c,r,m)/=LIS_rc%udef.and.&
                    HYMAP_routing_struc(n)%grarea(c,r)/=LIS_rc%udef)then
                   rnfsto_mm(c,r,m)=1000*HYMAP_routing_struc(n)%rnfsto(c,r,m)/HYMAP_routing_struc(n)%grarea(c,r)
                   bsfsto_mm(c,r,m)=1000*HYMAP_routing_struc(n)%bsfsto(c,r,m)/HYMAP_routing_struc(n)%grarea(c,r)
                 else
                   rnfsto_mm(c,r,m)=LIS_rc%udef
                   bsfsto_mm(c,r,m)=LIS_rc%udef             
                 endif
               enddo
             enddo
             
              !where(HYMAP_routing_struc(n)%rnfsto(:,:,m)/=LIS_rc%udef.and.&
              !      HYMAP_routing_struc(n)%bsfsto(:,:,m)/=LIS_rc%udef.and.&
              !      HYMAP_routing_struc(n)%grarea/=LIS_rc%udef)
              !  rnfsto_mm(:,:,m)=1000*HYMAP_routing_struc(n)%rnfsto(:,:,m)/HYMAP_routing_struc(n)%grarea
              !  bsfsto_mm(:,:,m)=1000*HYMAP_routing_struc(n)%bsfsto(:,:,m)/HYMAP_routing_struc(n)%grarea
              !end where
              
           enddo
           deallocate(surface_runoff)
           deallocate(baseflow)
        endif

        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivsto,&
             rivsto_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivdph,&
             rivdph_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivvel,&
             rivvel_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivout,&
             rivout_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%evpout,&
             evpout_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldout,&
             fldout_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldsto,&
             fldsto_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%flddph,&
             flddph_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldvel,&
             fldvel_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldfrc,&
             fldfrc_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldare,&
             fldare_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%sfcelv,&
             sfcelv_lvec,1)
        
        !ag (26Oct2017) - converting surface runoff and baseflow storage units from m3 to mm
        !call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rnfsto,&
        !     rnfsto_lvec,1)
        !
        !call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%bsfsto,&
        !     bsfsto_lvec,1)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,rnfsto_mm,rnfsto_lvec,1)
      
        call LIS_grid2patch(n,LIS_rc%lsm_index,bsfsto_mm,bsfsto_lvec,1)

        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RIVSTO,&
                value=rivsto_lvec(t),vlevel=1,unit="m3",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RIVDPH,&
                value=rivdph_lvec(t),vlevel=1,unit="m",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RIVVEL,&
                value=rivvel_lvec(t),vlevel=1,unit="m s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_STREAMFLOW,&
                value=rivout_lvec(t),vlevel=1,unit="m3 s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDEVAP,&
                value=evpout_lvec(t),vlevel=1,unit="m3",&   
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDOUT,&
                value=fldout_lvec(t),vlevel=1,unit="m3 s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDSTO,&
                value=fldsto_lvec(t),vlevel=1,unit="m3",&
                direction="-")
              
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDDPH,&
                value=flddph_lvec(t),vlevel=1,unit="m",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDVEL,&
                value=fldvel_lvec(t),vlevel=1,unit="m s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDFRC,&
                value=fldfrc_lvec(t),vlevel=1,unit="-",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDARE,&
                value=fldare_lvec(t),vlevel=1,unit="m2",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_SFCELV,&
                value=sfcelv_lvec(t),vlevel=1,unit="m",&
                direction="-")
           
           !ag (26Oct2017) - convert surface runoff and baseflow storage units from m3 to mm
           !call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
           !     value=rnfsto_lvec(t),vlevel=1,unit="m3",&
           !     direction="-")
           !
           !call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFSTO,&
           !     value=bsfsto_lvec(t),vlevel=1,unit="m3",&
           !     direction="-")

          call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
               value=rnfsto_lvec(t),vlevel=1,unit="mm",&
               direction="-")
          
          call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFSTO,&
               value=bsfsto_lvec(t),vlevel=1,unit="mm",&
               direction="-")  
           
        enddo
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
              
           else !read from previous output. 
              
              allocate(surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n)))
              allocate(baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n)))

              call readrunoffdata(trim(LIS_rc%runoffdatasource)//char(0),&
                   n,surface_runoff, baseflow)
              
           endif
           

           call model(LIS_rc%udef,&
                LIS_rc%gnc(n),&
                LIS_rc%gnr(n),&
                HYMAP_routing_struc(n)%inz,&
                nint(HYMAP_routing_struc(n)%dt),&
                HYMAP_routing_struc(n)%flowtype,&
                HYMAP_routing_struc(n)%linresflag,&
                HYMAP_routing_struc(n)%evapflag,&
                HYMAP_routing_struc(n)%seqx,&
                HYMAP_routing_struc(n)%seqy,&
                HYMAP_routing_struc(n)%nseqriv,&
                HYMAP_routing_struc(n)%nseqall,&
                HYMAP_routing_struc(n)%nextx,&
                HYMAP_routing_struc(n)%nexty,&
                HYMAP_routing_struc(n)%elevtn,&
                HYMAP_routing_struc(n)%nxtdst,&
                HYMAP_routing_struc(n)%grarea,&
                HYMAP_routing_struc(n)%fldgrd,&
                HYMAP_routing_struc(n)%fldman,&
                HYMAP_routing_struc(n)%fldstomax,&
                HYMAP_routing_struc(n)%rivman,&
                HYMAP_routing_struc(n)%rivelv,&
                HYMAP_routing_struc(n)%rivstomax,&
                HYMAP_routing_struc(n)%rivlen,&
                HYMAP_routing_struc(n)%rivwth,&
                HYMAP_routing_struc(n)%rivhgt,&
                HYMAP_routing_struc(n)%rivare,&
                HYMAP_routing_struc(n)%rslpmin,&
                HYMAP_routing_struc(n)%trnoff,&
                HYMAP_routing_struc(n)%tbsflw,&
                HYMAP_routing_struc(n)%cntime,&
                HYMAP_routing_struc(n)%mask,&
                surface_runoff,&
                baseflow,&
                   
                   !==== TEMPORARILY COMMENTED ====
                   !input evap        
                   !HYMAP_routing_struc(n)%tsfc,&
                   !HYMAP_routing_struc(n)%psur,&
                   !HYMAP_routing_struc(n)%pres,&
                   !HYMAP_routing_struc(n)%tair,&
                   !HYMAP_routing_struc(n)%relh,&
                   !HYMAP_routing_struc(n)%wind,&
                   !HYMAP_routing_struc(n)%zref,&
                   
                   !outputs
                HYMAP_routing_struc(n)%rivsto,&
                HYMAP_routing_struc(n)%rivdph,&
                HYMAP_routing_struc(n)%rivvel,&
                HYMAP_routing_struc(n)%rivout,&
                HYMAP_routing_struc(n)%evpout,&
                HYMAP_routing_struc(n)%fldout,&
                HYMAP_routing_struc(n)%fldsto,&
                HYMAP_routing_struc(n)%flddph,&
                HYMAP_routing_struc(n)%fldvel,&
                HYMAP_routing_struc(n)%fldfrc,&
                HYMAP_routing_struc(n)%fldare,&
                HYMAP_routing_struc(n)%sfcelv,&
                HYMAP_routing_struc(n)%rnfsto,&
                HYMAP_routing_struc(n)%bsfsto)

!========================================================================			
! the following is done to distribute the global arrays to the individual 
! processors, so that they can use the generic LIS interfaces for writing
! output -- which works based on the assumption of multiprocessors. 
! When the HYMAP routine itself is parallelized, this can (should) be taken
! out. 

           !ag (26Oct2017) - converting surface runoff and baseflow storage units from m3 to mm
           do c=1,LIS_rc%gnc(n)
             do r=1,LIS_rc%gnr(n)
               if(HYMAP_routing_struc(n)%rnfsto(c,r,1)/=LIS_rc%udef.and.&
                 HYMAP_routing_struc(n)%bsfsto(c,r,1)/=LIS_rc%udef.and.&
                 HYMAP_routing_struc(n)%grarea(c,r)/=LIS_rc%udef)then
                   rnfsto_mm(c,r,1)=1000*HYMAP_routing_struc(n)%rnfsto(c,r,1)/HYMAP_routing_struc(n)%grarea(c,r)
                   bsfsto_mm(c,r,1)=1000*HYMAP_routing_struc(n)%bsfsto(c,r,1)/HYMAP_routing_struc(n)%grarea(c,r)
               else
                 rnfsto_mm(c,r,1)=LIS_rc%udef
                 bsfsto_mm(c,r,1)=LIS_rc%udef             
               endif
             enddo
           enddo
                   
           !where(HYMAP_routing_struc(n)%rnfsto(:,:,m)/=LIS_rc%udef.and.&
           !      HYMAP_routing_struc(n)%bsfsto(:,:,m)/=LIS_rc%udef.and.&
           !      HYMAP_routing_struc(n)%grarea/=LIS_rc%udef)
           !  rnfsto_mm(:,:,1)=1000*HYMAP_routing_struc(n)%rnfsto(:,:,1)/HYMAP_routing_struc(n)%grarea
           !  bsfsto_mm(:,:,1)=1000*HYMAP_routing_struc(n)%bsfsto(:,:,1)/HYMAP_routing_struc(n)%grarea
           !end where
        
           deallocate(surface_runoff)
           deallocate(baseflow)
        endif
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivsto(:,:,1),&
             rivsto_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivdph(:,:,1),&
             rivdph_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivvel(:,:,1),&
             rivvel_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rivout(:,:,1),&
             rivout_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%evpout(:,:,1),&
             evpout_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldout(:,:,1),&
             fldout_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldsto(:,:,1),&
             fldsto_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%flddph(:,:,1),&
             flddph_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldvel(:,:,1),&
             fldvel_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldfrc(:,:,1),&
             fldfrc_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%fldare(:,:,1),&
             fldare_lvec,dummy)
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%sfcelv(:,:,1),&
             sfcelv_lvec,dummy)
        
        !ag (26Oct2017) - converting surface runoff and baseflow storage units from m3 to mm
        !call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%rnfsto(:,:,1),&
        !     rnfsto_lvec,dummy)
        !
        !call LIS_grid2patch(n,LIS_rc%lsm_index,HYMAP_routing_struc(n)%bsfsto(:,:,1),&
        !     bsfsto_lvec,dummy)

        call LIS_grid2patch(n,LIS_rc%lsm_index,rnfsto_mm(:,:,1),rnfsto_lvec,dummy)
      
        call LIS_grid2patch(n,LIS_rc%lsm_index,bsfsto_mm(:,:,1),bsfsto_lvec,dummy)

        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RIVSTO,&
                value=rivsto_lvec(t),vlevel=1,unit="m3",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RIVDPH,&
                value=rivdph_lvec(t),vlevel=1,unit="m",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RIVVEL,&
                value=rivvel_lvec(t),vlevel=1,unit="m s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_STREAMFLOW,&
                value=rivout_lvec(t),vlevel=1,unit="m3 s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDEVAP,&
                value=evpout_lvec(t),vlevel=1,unit="m3",&   
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDOUT,&
                value=fldout_lvec(t),vlevel=1,unit="m3 s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDSTO,&
                value=fldsto_lvec(t),vlevel=1,unit="m3",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDDPH,&
                value=flddph_lvec(t),vlevel=1,unit="m",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDVEL,&
                value=fldvel_lvec(t),vlevel=1,unit="m s-1",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDFRC,&
                value=fldfrc_lvec(t),vlevel=1,unit="-",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_FLDARE,&
                value=fldare_lvec(t),vlevel=1,unit="m2",&
                direction="-")
           
           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_SFCELV,&
                value=sfcelv_lvec(t),vlevel=1,unit="m",&
                direction="-")
           
           !ag (26Oct2017) - convert surface runoff and baseflow storage units from m3 to mm
           !call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
           !     value=rnfsto_lvec(t),vlevel=1,unit="m3",&
           !     direction="-")
           !
           !call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFSTO,&
           !     value=bsfsto_lvec(t),vlevel=1,unit="m3",&
           !     direction="-")

          call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
               value=rnfsto_lvec(t),vlevel=1,unit="mm",&
               direction="-")
          
          call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFSTO,&
               value=bsfsto_lvec(t),vlevel=1,unit="mm",&
               direction="-")  

        enddo
                
     endif
     deallocate(rivsto_lvec)
     deallocate(rivdph_lvec)
     deallocate(rivvel_lvec)
     deallocate(rivout_lvec)
     deallocate(evpout_lvec)
     deallocate(fldout_lvec)
     deallocate(fldsto_lvec)
     deallocate(flddph_lvec)
     deallocate(fldvel_lvec)
     deallocate(fldfrc_lvec)
     deallocate(fldare_lvec)
     deallocate(sfcelv_lvec)
     deallocate(rnfsto_lvec)
     deallocate(bsfsto_lvec)

     if(LIS_masterproc)then
       if(allocated(rnfsto_mm))deallocate(rnfsto_mm)
       if(allocated(bsfsto_mm))deallocate(bsfsto_mm)
     endif
     
  endif
end subroutine HYMAP_routing_run
