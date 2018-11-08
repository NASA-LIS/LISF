!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
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
!             doi:10.1175/JHM-D-12-021.1
!
!             Getirana, A., Peters-Lidard, C., Rodell, M., Bates, P.D., 2017. 
!             Trade-off between cost and accuracy in large-scale surface water 
!             dynamic modeling. Water Resources Research. doi: 10.1002/2017WR020519
!
! !REVISION HISTORY: 
! 08 Nov 2011: Augusto Getirana, Initial implementation in LIS based on the 
!                                HYMAP offline routing code. 
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia formulation, 
!                                adaptive time step and reservoir operation. 
! 13 Apr 2016: Augusto Getirana, Inclusion of option for hybrid runs with a 
!                                river flow map. 
!
! !USES: 
subroutine HYMAP2_routing_run(n)

  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_routingMod
  use LIS_logMod
  use LIS_historyMod
  use LIS_histDataMod
  use LIS_constantsMod
  use HYMAP2_routingMod
  use HYMAP2_evapMod
  use HYMAP2_initMod, only : HYMAP2_grid2vector,HYMAP2_vector2grid
  
  
  !use clsmf25_lsmMod, only : clsmf25_struc
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod 

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  implicit none
  
  integer, intent(in)   :: n  
  integer               :: m
  type(ESMF_Field)      :: sf_runoff_field
  type(ESMF_Field)      :: baseflow_field  
  real,   pointer       :: surface_runoff_t(:)
  real,   pointer       :: baseflow_t(:)
  real,   allocatable   :: surface_runoff(:)
  real,   allocatable   :: baseflow(:)
  real,   allocatable   :: tmpr(:,:),tmpb(:,:)
  
  real,   allocatable   :: evap(:)
  real,   allocatable   :: tair(:)
  real,   allocatable   :: qair(:)
  real,   allocatable   :: wind(:)
  real,   allocatable   :: pres(:)
  real,   allocatable   :: qnet(:)
  real,   allocatable   :: tmp_tmp(:,:),tmp_q2(:,:)
  real,   allocatable   :: tmp_wind(:,:),tmp_psurf(:,:),tmp_qnet(:,:)
  real,   allocatable   :: gvar1(:),gvar2(:),gvar3(:),gvar4(:),gvar5(:),gvar6(:),gvar7(:)
  
  real,   allocatable   :: rnfsto_mm(:,:),bsfsto_mm(:,:)

  real,   allocatable   :: tmp_nensem(:,:,:)

  type(ESMF_Field)      :: evapotranspiration_field
  type(ESMF_Field)      :: potential_evaporation_field
  real,   pointer       :: evapotranspiration_t(:)
  real,   pointer       :: potential_evaporation_t(:)
  !real,   allocatable   :: dif_evap(:)
  !real,   allocatable   :: tmppe(:,:),tmpet(:,:),tmpde(:,:)
  real,   allocatable   ::  tmpet(:,:)
 
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
  real,   allocatable   :: rnfdwi_lvec(:)
  real,   allocatable   :: bsfdwi_lvec(:)
  real,   allocatable   :: surfws_lvec(:)

  real,   allocatable   :: ewat_lvec(:)
  real,   allocatable   :: edif_lvec(:)

  integer               :: status
  logical               :: alarmCheck
  integer               :: c,r,t
  integer               :: ios, nid,qsid,qsbid

  integer            :: tid
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),swd(:),lwd(:)
  real,pointer       :: psurf(:)
  
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "HYMAP2 router model alarm")
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
     allocate(rnfdwi_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(bsfdwi_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(surfws_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))

     allocate(ewat_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(edif_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))

     rivsto_lvec = LIS_rc%udef
     rivdph_lvec = LIS_rc%udef
     rivvel_lvec = LIS_rc%udef
     rivout_lvec = LIS_rc%udef
     evpout_lvec = LIS_rc%udef
     fldout_lvec = LIS_rc%udef
     fldsto_lvec = LIS_rc%udef
     flddph_lvec = LIS_rc%udef
     fldvel_lvec = LIS_rc%udef
     fldfrc_lvec = LIS_rc%udef
     fldare_lvec = LIS_rc%udef
     sfcelv_lvec = LIS_rc%udef
     rnfsto_lvec = LIS_rc%udef
     bsfsto_lvec = LIS_rc%udef
     rnfdwi_lvec = LIS_rc%udef
     bsfdwi_lvec = LIS_rc%udef
     surfws_lvec = LIS_rc%udef

     ewat_lvec = LIS_rc%udef
     edif_lvec = LIS_rc%udef
     
     if(LIS_masterproc)then
        allocate(surface_runoff(HYMAP2_routing_struc(n)%nseqall))     
        allocate(baseflow(HYMAP2_routing_struc(n)%nseqall))
        allocate(tmpr(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        allocate(tmpb(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        if(HYMAP2_routing_struc(n)%evapflag.ne.0)then
           allocate(tmp_tmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
           allocate(tmp_q2(LIS_rc%gnc(n),LIS_rc%gnr(n)))
           allocate(tmp_psurf(LIS_rc%gnc(n),LIS_rc%gnr(n)))
           allocate(tmp_wind(LIS_rc%gnc(n),LIS_rc%gnr(n)))
           allocate(tmp_qnet(LIS_rc%gnc(n),LIS_rc%gnr(n)))

           allocate(evap(HYMAP2_routing_struc(n)%nseqall))     
           allocate(tair(HYMAP2_routing_struc(n)%nseqall))     
           allocate(qair(HYMAP2_routing_struc(n)%nseqall))
           allocate(pres(HYMAP2_routing_struc(n)%nseqall))     
           allocate(wind(HYMAP2_routing_struc(n)%nseqall))     
           allocate(qnet(HYMAP2_routing_struc(n)%nseqall))
           allocate(tmpet(LIS_rc%gnc(n),LIS_rc%gnr(n)))
           
        endif

        if(HYMAP2_routing_struc(n)%useens.eq.0) then 
           allocate(rnfsto_mm(HYMAP2_routing_struc(n)%nseqall,1))
           allocate(bsfsto_mm(HYMAP2_routing_struc(n)%nseqall,1))
        elseif(HYMAP2_routing_struc(n)%useens.eq.1) then
           allocate(rnfsto_mm(HYMAP2_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
           allocate(bsfsto_mm(HYMAP2_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
        endif
     endif

     if(HYMAP2_routing_struc(n)%useens.eq.1) then 
        if(LIS_masterproc) then     
           allocate(tmp_nensem(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n)))
           !run the routing model at 1 hour output interval   
           call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",&
                sf_runoff_field,rc=status)
           call LIS_verify(status, "ESMF_StateGet failed for Surface Runoff")

           call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
                baseflow_field, rc=status)
           call LIS_verify(status, "ESMF_StateGet failed for Subsurface Runoff")

           call ESMF_FieldGet(sf_runoff_field,localDE=0,&
                farrayPtr=surface_runoff_t,rc=status)
           call LIS_verify(status, "ESMF_FieldGet failed for Surface Runoff")

           call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow_t,&
                rc=status)
           call LIS_verify(status, "ESMF_FieldGet failed for Subsurface Runoff")


           do m=1,LIS_rc%nensem(n)

              surface_runoff = 0.0
              baseflow = 0.0

              !temporary solution  
              tmpr=0.
              call LIS_tile2grid(n,m,tmpr,surface_runoff_t)
              call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmpr,surface_runoff)

              !temporary solution  
              tmpb=0.
              call LIS_tile2grid(n,m,tmpb,baseflow_t)
              call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmpb,baseflow)

              call HYMAP2_model(n,real(HYMAP2_routing_struc(n)%imis),&
                   LIS_rc%gnc(n),&
                   LIS_rc%gnr(n),&
                   LIS_rc%yr,&
                   LIS_rc%mo,&
                   LIS_rc%da,&
                   LIS_rc%hr,&
                   LIS_rc%mn,&
                   LIS_rc%ss,&				   
                   HYMAP2_routing_struc(n)%nseqall,&
                   HYMAP2_routing_struc(n)%nz,&
                   HYMAP2_routing_struc(n)%dt,&
                   HYMAP2_routing_struc(n)%flowmap,&
                   HYMAP2_routing_struc(n)%linresflag,&
                   HYMAP2_routing_struc(n)%evapflag,&
                                !ag (19Jan2016)
                   HYMAP2_routing_struc(n)%rivout_pre(:,m),&
                   HYMAP2_routing_struc(n)%rivdph_pre(:,m),&
                   HYMAP2_routing_struc(n)%fldout_pre(:,m),&
                   HYMAP2_routing_struc(n)%flddph_pre(:,m),&
                   HYMAP2_routing_struc(n)%fldelv1(:,m),&
                   HYMAP2_routing_struc(n)%grv,&
                   HYMAP2_routing_struc(n)%cadp,&
                   HYMAP2_routing_struc(n)%steptype,&
                   HYMAP2_routing_struc(n)%resopflag,&                   
                   HYMAP2_routing_struc(n)%floodflag,&                   
                   HYMAP2_routing_struc(n)%outlet,&
                   HYMAP2_routing_struc(n)%next,&
                   HYMAP2_routing_struc(n)%elevtn,&
                   HYMAP2_routing_struc(n)%nxtdst,&
                   HYMAP2_routing_struc(n)%grarea,&	   
                   HYMAP2_routing_struc(n)%fldgrd,&
                   HYMAP2_routing_struc(n)%fldman,&
                   HYMAP2_routing_struc(n)%fldhgt,&
                   HYMAP2_routing_struc(n)%fldstomax,&
                   HYMAP2_routing_struc(n)%rivman,&
                   HYMAP2_routing_struc(n)%rivelv,&
                   HYMAP2_routing_struc(n)%rivstomax,&
                   HYMAP2_routing_struc(n)%rivlen,&
                   HYMAP2_routing_struc(n)%rivwth,&
                   HYMAP2_routing_struc(n)%rivhgt,&
                   HYMAP2_routing_struc(n)%rivare,&
                   HYMAP2_routing_struc(n)%rslpmin,&
                   HYMAP2_routing_struc(n)%trnoff,&
                   HYMAP2_routing_struc(n)%tbsflw,&
                   HYMAP2_routing_struc(n)%cntime,&
                   HYMAP2_routing_struc(n)%dwiflag,&
                   HYMAP2_routing_struc(n)%rnfdwi_ratio,&
                   HYMAP2_routing_struc(n)%bsfdwi_ratio,&
                   surface_runoff,&
                   baseflow,&
                   HYMAP2_routing_struc(n)%edif(:,m),&
                   HYMAP2_routing_struc(n)%rivsto(:,m),&
                   HYMAP2_routing_struc(n)%rivdph(:,m),&
                   HYMAP2_routing_struc(n)%rivvel(:,m),&
                   HYMAP2_routing_struc(n)%rivout(:,m),&
                   HYMAP2_routing_struc(n)%evpout(:,m),&
                   HYMAP2_routing_struc(n)%fldout(:,m),&
                   HYMAP2_routing_struc(n)%fldsto(:,m),&
                   HYMAP2_routing_struc(n)%flddph(:,m),&
                   HYMAP2_routing_struc(n)%fldvel(:,m),&
                   HYMAP2_routing_struc(n)%fldfrc(:,m),&
                   HYMAP2_routing_struc(n)%fldare(:,m),&
                   HYMAP2_routing_struc(n)%sfcelv(:,m),&
                   HYMAP2_routing_struc(n)%rnfsto(:,m),&
                   HYMAP2_routing_struc(n)%bsfsto(:,m),&
                   HYMAP2_routing_struc(n)%rnfdwi(:,m),&
                   HYMAP2_routing_struc(n)%bsfdwi(:,m),&
                   HYMAP2_routing_struc(n)%surfws(:,m),&
                   HYMAP2_routing_struc(n)%dtaout(:,m) )
              !========================================================================			
              ! the following is done to distribute the global arrays to the individual 
              ! processors, so that they can use the generic LIS interfaces for writing
              ! output -- which works based on the assumption of multiprocessors. 
              ! When the HYMAP routine itself is parallelized, this can (should) be taken
              ! out. 
              rnfsto_mm(:,m)=1000*HYMAP2_routing_struc(n)%rnfsto(:,m)/HYMAP2_routing_struc(n)%grarea
              bsfsto_mm(:,m)=1000*HYMAP2_routing_struc(n)%bsfsto(:,m)/HYMAP2_routing_struc(n)%grarea
           enddo

        endif

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,&
                tmp_nensem,HYMAP2_routing_struc(n)%rivsto)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             rivsto_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%rivdph)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             rivdph_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%rivvel)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             rivvel_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%rivout)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             rivout_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%evpout)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             evpout_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%fldout)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             fldout_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%fldsto)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             fldsto_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%flddph)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             flddph_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%fldvel)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             fldvel_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%fldfrc)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             fldfrc_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%fldare)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             fldare_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%sfcelv)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             sfcelv_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,rnfsto_mm)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             rnfsto_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,bsfsto_mm)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             bsfsto_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%rnfdwi)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             rnfdwi_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%bsfdwi)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             bsfdwi_lvec,1)

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n),HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_nensem,HYMAP2_routing_struc(n)%surfws)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem,&
             surfws_lvec,1)

        !do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        do t=1, LIS_rc%ntiles(n)        
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
                value=evpout_lvec(t),vlevel=1,unit="kg m-2 s-1",&   
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

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
                value=rnfsto_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFSTO,&
                value=bsfsto_lvec(t),vlevel=1,unit="mm",&
                direction="-")        

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFDWI,&
                value=rnfdwi_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFDWI,&
                value=bsfdwi_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_SURFWS,&
                value=surfws_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_ewat,&
                value=ewat_lvec(t),vlevel=1,unit="kg m-2 s-1",&   
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_edif,&
                value=edif_lvec(t),vlevel=1,unit="kg m-2 s-1",&   
                direction="-")
        enddo
     else    
        !ag (03May2017)
        !import evaporation from open water       
        if(HYMAP2_routing_struc(n)%evapflag.ne.0)then !"compute" 
           if(LIS_rc%lsm.ne."none") then !from current lsm run

              if(LIS_masterproc)then     
                 call ESMF_StateGet(LIS_runoff_state(n),"Total Evapotranspiration",&
                      evapotranspiration_field, rc=status)
                 call LIS_verify(status, "HYMAP2_routing_run: ESMF_StateGet failed for Total Evapotranspiration")

                 call ESMF_FieldGet(evapotranspiration_field,localDE=0,&
                      farrayPtr=evapotranspiration_t,rc=status)
                 call LIS_verify(status, "HYMAP2_routing_run: ESMF_FieldGet failed for Total Evapotranspiration")
          
!print*,minval(gvar6),maxval(gvar6),minval(lwd),maxval(lwd)
!print*,minval(gvar7),maxval(gvar7),minval(swd),maxval(swd)

                 call LIS_tile2grid(n,tmpet,evapotranspiration_t,1)
                 call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                      HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                      HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmpet,evap)
              endif

              call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,rc=status)
              call LIS_verify(status,'HYMAP2_routing_run: ESMF_FieldGet failed for Tair')

              call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,rc=status)
              call LIS_verify(status,'HYMAP2_routing_run: ESMF_FieldGet failed for Qair')

              call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,rc=status)
              call LIS_verify(status,'HYMAP2_routing_run: ESMF_FieldGet failed for SWdown')

              call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,rc=status)
              call LIS_verify(status,'HYMAP2_routing_run: ESMF_FieldGet failed for LWdown')

              call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,rc=status)
              call LIS_verify(status,'HYMAP2_routing_run: ESMF_FieldGet failed for Wind_E')

              call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,rc=status)
              call LIS_verify(status,'HYMAP2_routing_run: ESMF_FieldGet failed for Wind_N')

              call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,rc=status)
              call LIS_verify(status, 'HYMAP2_routing_run: ESMF_FieldGet failed for PSurf')


              call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp,rc=status)
              call LIS_verify(status)

              call ESMF_FieldGet(q2Field, localDE=0, farrayPtr=q2,rc=status)
              call LIS_verify(status)

              call ESMF_FieldGet(uField, localDE=0, farrayPtr=uwind,rc=status)
              call LIS_verify(status)

              call ESMF_FieldGet(vField, localDE=0, farrayPtr=vwind,rc=status)
              call LIS_verify(status)

              call ESMF_FieldGet(psurfField, localDE=0, farrayPtr=psurf,rc=status)
              call LIS_verify(status)

              call ESMF_FieldGet(lwdField, localDE=0, farrayPtr=lwd,rc=status)
              call LIS_verify(status)

              call ESMF_FieldGet(swdField, localDE=0, farrayPtr=swd,rc=status)
              call LIS_verify(status)


              call LIS_gather_tiled_vector_output(n, gvar1, tmp)
              call LIS_gather_tiled_vector_output(n, gvar2, q2)
              call LIS_gather_tiled_vector_output(n, gvar3, psurf)
              call LIS_gather_tiled_vector_output(n, gvar4, uwind)
              call LIS_gather_tiled_vector_output(n, gvar5, vwind)
              call LIS_gather_tiled_vector_output(n, gvar6, lwd)
              call LIS_gather_tiled_vector_output(n, gvar7, swd)


              if(LIS_masterproc)then     
                 call LIS_tile2grid(n,tmp_tmp,gvar1-273.15,1)
                 call LIS_tile2grid(n,tmp_q2,gvar2,1)
                 call LIS_tile2grid(n,tmp_psurf,gvar3/1e3,1)
                 call LIS_tile2grid(n,tmp_wind,sqrt(gvar4**2+gvar5**2),1)
                 call LIS_tile2grid(n,tmp_qnet,(gvar6+gvar7)/1e6,1)

                 deallocate(gvar1)
                 deallocate(gvar2)
                 deallocate(gvar3)
                 deallocate(gvar4)
                 deallocate(gvar5)
                 deallocate(gvar6)
                 deallocate(gvar7)

                 call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                      HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                      HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_tmp,tair)
                 call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                      HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                      HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_q2,qair)
                 call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                      HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                      HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_psurf,pres)
                 call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                      HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                      HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_qnet,qnet)
                 call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                      HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                      HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmp_wind,wind)
              endif

           else
              print*,'wrong value for evapflag or LSM id'
              stop
           endif

           if(LIS_masterproc)then     
              call HYMAP2_evap_main(HYMAP2_routing_struc(n)%evapflag,n,&
                   HYMAP2_routing_struc(n)%nseqall,real(HYMAP2_routing_struc(n)%imis),&
                   HYMAP2_routing_struc(n)%outlet,pres,tair,qair,&
                   wind,qnet,evap,HYMAP2_routing_struc(n)%ewat,HYMAP2_routing_struc(n)%edif)
           endif
        endif

        if(LIS_masterproc)then     
           allocate(tmp_nensem(LIS_rc%gnc(n),LIS_rc%gnr(n),1))
           tmpr=0.
           tmpb=0.

           !import surface runoff and baseflow
           if(LIS_rc%lsm.ne."none") then !from current lsm run
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
              !call LIS_tile2grid(n,surface_runoff,surface_runoff_t,1)
              !call LIS_tile2grid(n,baseflow,baseflow_t,1)

              !temporary solution  
              call LIS_tile2grid(n,tmpr,surface_runoff_t,1)
              call LIS_tile2grid(n,tmpb,baseflow_t,1)
           else !read from previous output. 
              call readrunoffdata(trim(LIS_rc%runoffdatasource)//char(0),&
                   n,tmpr, tmpb)
           endif
           call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmpr,surface_runoff)
           call HYMAP2_grid2vector(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,HYMAP2_routing_struc(n)%seqy,tmpb,baseflow)


           call HYMAP2_model(n,real(HYMAP2_routing_struc(n)%imis),&
                LIS_rc%gnc(n),&
                LIS_rc%gnr(n),&
                LIS_rc%yr,&
                LIS_rc%mo,&
                LIS_rc%da,&
                LIS_rc%hr,&
                LIS_rc%mn,&
                LIS_rc%ss,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%nz,&
                HYMAP2_routing_struc(n)%dt,&
                HYMAP2_routing_struc(n)%flowmap,&
                HYMAP2_routing_struc(n)%linresflag,&
                HYMAP2_routing_struc(n)%evapflag,&
                                !ag (19Jan2016)
                HYMAP2_routing_struc(n)%rivout_pre,&
                HYMAP2_routing_struc(n)%rivdph_pre,&
                HYMAP2_routing_struc(n)%fldout_pre,&
                HYMAP2_routing_struc(n)%flddph_pre,&
                HYMAP2_routing_struc(n)%fldelv1,&
                HYMAP2_routing_struc(n)%grv,&
                HYMAP2_routing_struc(n)%cadp,&
                HYMAP2_routing_struc(n)%steptype,&
                HYMAP2_routing_struc(n)%resopflag,&                
                HYMAP2_routing_struc(n)%floodflag,&                   
                HYMAP2_routing_struc(n)%outlet,&
                HYMAP2_routing_struc(n)%next,&
                HYMAP2_routing_struc(n)%elevtn,&
                HYMAP2_routing_struc(n)%nxtdst,&
                HYMAP2_routing_struc(n)%grarea,&	   
                HYMAP2_routing_struc(n)%fldgrd,&
                HYMAP2_routing_struc(n)%fldman,&
                HYMAP2_routing_struc(n)%fldhgt,&
                HYMAP2_routing_struc(n)%fldstomax,&
                HYMAP2_routing_struc(n)%rivman,&
                HYMAP2_routing_struc(n)%rivelv,&
                HYMAP2_routing_struc(n)%rivstomax,&
                HYMAP2_routing_struc(n)%rivlen,&
                HYMAP2_routing_struc(n)%rivwth,&
                HYMAP2_routing_struc(n)%rivhgt,&
                HYMAP2_routing_struc(n)%rivare,&
                HYMAP2_routing_struc(n)%rslpmin,&
                HYMAP2_routing_struc(n)%trnoff,&
                HYMAP2_routing_struc(n)%tbsflw,&
                HYMAP2_routing_struc(n)%cntime,&
                HYMAP2_routing_struc(n)%dwiflag,&
                HYMAP2_routing_struc(n)%rnfdwi_ratio,&
                HYMAP2_routing_struc(n)%bsfdwi_ratio,&
                surface_runoff,&
                baseflow,&
                HYMAP2_routing_struc(n)%edif,&
                HYMAP2_routing_struc(n)%rivsto,&
                HYMAP2_routing_struc(n)%rivdph,&
                HYMAP2_routing_struc(n)%rivvel,&
                HYMAP2_routing_struc(n)%rivout,&
                HYMAP2_routing_struc(n)%evpout,&
                HYMAP2_routing_struc(n)%fldout,&
                HYMAP2_routing_struc(n)%fldsto,&
                HYMAP2_routing_struc(n)%flddph,&
                HYMAP2_routing_struc(n)%fldvel,&
                HYMAP2_routing_struc(n)%fldfrc,&
                HYMAP2_routing_struc(n)%fldare,&
                HYMAP2_routing_struc(n)%sfcelv,&
                HYMAP2_routing_struc(n)%rnfsto,&
                HYMAP2_routing_struc(n)%bsfsto,&
                HYMAP2_routing_struc(n)%rnfdwi,&
                HYMAP2_routing_struc(n)%bsfdwi,&
                HYMAP2_routing_struc(n)%surfws,&
                HYMAP2_routing_struc(n)%dtaout)
           !========================================================================			
           ! the following is done to distribute the global arrays to the individual 
           ! processors, so that they can use the generic LIS interfaces for writing
           ! output -- which works based on the assumption of multiprocessors. 
           ! When the HYMAP routine itself is parallelized, this can (should) be taken
           ! out. 

           rnfsto_mm(:,1)=1000*HYMAP2_routing_struc(n)%rnfsto(:,1)/HYMAP2_routing_struc(n)%grarea
           bsfsto_mm(:,1)=1000*HYMAP2_routing_struc(n)%bsfsto(:,1)/HYMAP2_routing_struc(n)%grarea
        else
           allocate(tmp_nensem(1,1,1))
        endif

        if(LIS_masterproc)then
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%rivsto)    
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             rivsto_lvec)

        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%rivdph)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             rivdph_lvec)
        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%rivvel)             
        endif

        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             rivvel_lvec)
        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%rivout)  
        endif
        
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             rivout_lvec)
        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%evpout)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             evpout_lvec)
        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%fldout)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             fldout_lvec)
        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%fldsto)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             fldsto_lvec)

        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%flddph)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             flddph_lvec)
        
        if(LIS_masterproc) then 
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%fldvel)  
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             fldvel_lvec)
        if(LIS_masterproc) then            
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%fldfrc)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             fldfrc_lvec)
        if(LIS_masterproc) then            

           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),HYMAP2_routing_struc(n)%fldare)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             fldare_lvec)
        if(LIS_masterproc) then            
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%sfcelv)    
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             sfcelv_lvec)

        if(LIS_masterproc) then            
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                rnfsto_mm)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
                   rnfsto_lvec)
        
        if(LIS_masterproc) then            
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),bsfsto_mm)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             bsfsto_lvec)
        if(LIS_masterproc) then            
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%rnfdwi) 
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             rnfdwi_lvec)
        if(LIS_masterproc) then             
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),HYMAP2_routing_struc(n)%bsfdwi)   
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             bsfdwi_lvec)
        if(LIS_masterproc) then           
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),&
                HYMAP2_routing_struc(n)%surfws)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             surfws_lvec)
        if(LIS_masterproc) then           
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),HYMAP2_routing_struc(n)%ewat)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             ewat_lvec)
        if(LIS_masterproc) then           
           call HYMAP2_vector2grid(LIS_rc%gnc(n),LIS_rc%gnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1),HYMAP2_routing_struc(n)%edif)             
        endif
        call LIS_grid2patch(n,LIS_rc%lsm_index,tmp_nensem(:,:,1),&
             edif_lvec)

        deallocate(tmp_nensem)

        !do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        do t=1, LIS_rc%ntiles(n)  
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
                value=evpout_lvec(t),vlevel=1,unit="kg m-2 s-1",&   
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

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFSTO,&
                value=rnfsto_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFSTO,&
                value=bsfsto_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_RNFDWI,&
                value=rnfdwi_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_BSFDWI,&
                value=bsfdwi_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_SURFWS,&
                value=surfws_lvec(t),vlevel=1,unit="mm",&
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_ewat,&
                value=ewat_lvec(t),vlevel=1,unit="kg m-2 s-1",&   
                direction="-")

           call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_edif,&
                value=edif_lvec(t),vlevel=1,unit="kg m-2 s-1",&   
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
     deallocate(rnfdwi_lvec)
     deallocate(bsfdwi_lvec)
     deallocate(surfws_lvec)
     deallocate(ewat_lvec)
     deallocate(edif_lvec)

     if(LIS_masterproc)then
        deallocate(surface_runoff)     
        deallocate(baseflow)
        deallocate(tmpr)
        deallocate(tmpb)
        !deallocate(tmp_nensem)

        deallocate(rnfsto_mm)
        deallocate(bsfsto_mm)

        !ag (22Sep2016)
        if(HYMAP2_routing_struc(n)%evapflag.ne.0)then
           deallocate(tmp_tmp)
           deallocate(tmp_q2)
           deallocate(tmp_qnet)
           deallocate(tmp_wind)
           deallocate(tmp_psurf)
           deallocate(tair)     
           deallocate(qair)
           deallocate(pres)     
           deallocate(wind)     
           deallocate(qnet)
           deallocate(tmpet)
        endif
     endif

  endif

end subroutine HYMAP2_routing_run
