!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
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
! 27 Apr 2020: Augusto Getirana,  Added support for urban drainage 
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
  use HYMAP2_initMod
  
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod 

  use LIS_mpiMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  implicit none
  
  integer, intent(in)   :: n  
  integer               :: m
  type(ESMF_Field)      :: sf_runoff_field
  type(ESMF_Field)      :: baseflow_field  
  type(ESMF_Field)      :: rivsto_field  
  type(ESMF_Field)      :: fldsto_field
  type(ESMF_Field)      :: fldfrc_field
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
  real,   pointer       :: rivstotmp_lvec(:)
  real,   pointer       :: fldstotmp_lvec(:)
  real,   pointer       :: fldfrctmp_lvec(:)
  real                  :: fldfrctmp1_lvec(LIS_rc%ntiles(n))
  
  !ag(27Apr2020)
  real,   allocatable   :: drsto_lvec(:)
  real,   allocatable   :: drout_lvec(:)

  integer               :: status
  logical               :: alarmCheck
  integer               :: c,r,t
  integer               :: ios, nid,qsid,qsbid

  integer            :: tid
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),swd(:),lwd(:)
  real,pointer       :: psurf(:)

  integer            :: i,ix, iy, ix1, iy1
!TBD:SVK - need to redo the code to work with ntiles rather than npatches. 
!
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "HYMAP2 router model alarm")
  if(alarmCheck) then

     allocate(rivsto_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(rivdph_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(rivvel_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(rivout_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(evpout_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(fldout_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(fldsto_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(flddph_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(fldvel_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(fldfrc_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(fldare_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(sfcelv_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(rnfsto_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(bsfsto_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(rnfdwi_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(bsfdwi_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(surfws_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     
     allocate(ewat_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(edif_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     
     !ag(27Apr2020)
     allocate(drsto_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))
     allocate(drout_lvec(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)))

     
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

     drsto_lvec = LIS_rc%udef
     drout_lvec = LIS_rc%udef

     allocate(surface_runoff(HYMAP2_routing_struc(n)%nseqall))     
     allocate(baseflow(HYMAP2_routing_struc(n)%nseqall))
     allocate(tmpr(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     allocate(tmpb(LIS_rc%lnc(n),LIS_rc%lnr(n)))

     if(HYMAP2_routing_struc(n)%evapflag.ne.0)then
        allocate(tmp_tmp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        allocate(tmp_q2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        allocate(tmp_psurf(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        allocate(tmp_wind(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        allocate(tmp_qnet(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        
        allocate(evap(HYMAP2_routing_struc(n)%nseqall))     
        allocate(tair(HYMAP2_routing_struc(n)%nseqall))     
        allocate(qair(HYMAP2_routing_struc(n)%nseqall))
        allocate(pres(HYMAP2_routing_struc(n)%nseqall))     
        allocate(wind(HYMAP2_routing_struc(n)%nseqall))     
        allocate(qnet(HYMAP2_routing_struc(n)%nseqall))
        allocate(tmpet(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        
     endif
     
     if(HYMAP2_routing_struc(n)%useens.eq.0) then 
        allocate(rnfsto_mm(HYMAP2_routing_struc(n)%nseqall,1))
        allocate(bsfsto_mm(HYMAP2_routing_struc(n)%nseqall,1))
     elseif(HYMAP2_routing_struc(n)%useens.eq.1) then
        allocate(rnfsto_mm(HYMAP2_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
        allocate(bsfsto_mm(HYMAP2_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
     endif
  
     !ag (03May2017)
     !import evaporation from open water     
     !TBD: SVK - block that needs update for parallelism  
     if(HYMAP2_routing_struc(n)%evapflag.ne.0)then !"compute" 
        call ESMF_StateGet(LIS_runoff_state(n),"Total Evapotranspiration",&
             evapotranspiration_field, rc=status)
        call LIS_verify(status, "HYMAP2_routing_run: ESMF_StateGet failed for Total Evapotranspiration")
        
        call ESMF_FieldGet(evapotranspiration_field,localDE=0,&
             farrayPtr=evapotranspiration_t,rc=status)
        call LIS_verify(status, "HYMAP2_routing_run: ESMF_FieldGet failed for Total Evapotranspiration")
        
        call LIS_tile2grid(n,tmpet,evapotranspiration_t)

        call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
             HYMAP2_routing_struc(n)%nseqall,&
             HYMAP2_routing_struc(n)%imis,&
             HYMAP2_routing_struc(n)%seqx,&
             HYMAP2_routing_struc(n)%seqy,tmpet,evap)

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
        
        
        call LIS_tile2grid(n,tmp_tmp,tmp-273.15)
        call LIS_tile2grid(n,tmp_q2,q2)
        call LIS_tile2grid(n,tmp_psurf,psurf/1e3)
        call LIS_tile2grid(n,tmp_wind,sqrt(uwind**2+vwind**2))
        call LIS_tile2grid(n,tmp_qnet,(lwd+swd)/1e6)
        
        
        call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
             HYMAP2_routing_struc(n)%nseqall,&
             HYMAP2_routing_struc(n)%imis,&
             HYMAP2_routing_struc(n)%seqx,&
             HYMAP2_routing_struc(n)%seqy,tmp_tmp,tair)

        call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
             HYMAP2_routing_struc(n)%nseqall,&
             HYMAP2_routing_struc(n)%imis,&
             HYMAP2_routing_struc(n)%seqx,&
             HYMAP2_routing_struc(n)%seqy,tmp_q2,qair)

        call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
             HYMAP2_routing_struc(n)%nseqall,&
             HYMAP2_routing_struc(n)%imis,&
             HYMAP2_routing_struc(n)%seqx,&
             HYMAP2_routing_struc(n)%seqy,tmp_psurf,pres)

        call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
             HYMAP2_routing_struc(n)%nseqall,&
             HYMAP2_routing_struc(n)%imis,&
             HYMAP2_routing_struc(n)%seqx,&
             HYMAP2_routing_struc(n)%seqy,tmp_qnet,qnet)

        call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
             HYMAP2_routing_struc(n)%nseqall,&
             HYMAP2_routing_struc(n)%imis,&
             HYMAP2_routing_struc(n)%seqx,&
             HYMAP2_routing_struc(n)%seqy,tmp_wind,wind)
        
        call HYMAP2_evap_main(HYMAP2_routing_struc(n)%evapflag,n,&
             HYMAP2_routing_struc(n)%nseqall,&
             real(HYMAP2_routing_struc(n)%imis),&
             HYMAP2_routing_struc(n)%outlet,&
             pres,tair,qair,wind,qnet,evap,&
             HYMAP2_routing_struc(n)%ewat,&
             HYMAP2_routing_struc(n)%edif)

     endif

     if(HYMAP2_routing_struc(n)%useens.eq.1) then !ensemble mode

        tmpr=0.
        tmpb=0.
        
        !import surface runoff and baseflow
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


        do m=1,LIS_rc%nensem(n)

           call LIS_tile2grid(n,m,tmpr,surface_runoff_t)
           call LIS_tile2grid(n,m,tmpb,baseflow_t)
        
           call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmpr,surface_runoff)
           call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmpb,baseflow)

           call HYMAP2_model(n,real(HYMAP2_routing_struc(n)%imis),&
                LIS_rc%lnc(n),&
                LIS_rc%lnr(n),&
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
                HYMAP2_routing_struc(n)%fldgrd(:,:,m),&
                HYMAP2_routing_struc(n)%fldman,&
                HYMAP2_routing_struc(n)%fldhgt,&
                HYMAP2_routing_struc(n)%fldstomax(:,:,m),&
                HYMAP2_routing_struc(n)%rivman(:,m),&
                HYMAP2_routing_struc(n)%rivelv,&
                HYMAP2_routing_struc(n)%rivstomax(:,m),&
                HYMAP2_routing_struc(n)%rivlen,&
                HYMAP2_routing_struc(n)%rivwth(:,m),&
                HYMAP2_routing_struc(n)%rivhgt(:,m),&
                HYMAP2_routing_struc(n)%rivare(:,m),&
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
                HYMAP2_routing_struc(n)%dtaout(:,m),&
                !ag (27Apr2020)
                !urban drainage variables/parameters
                HYMAP2_routing_struc(n)%flowtype,&
                HYMAP2_routing_struc(n)%drvel,&
                HYMAP2_routing_struc(n)%drtotwth,&
                HYMAP2_routing_struc(n)%drrad,&
                HYMAP2_routing_struc(n)%drman,&
                HYMAP2_routing_struc(n)%drslp,&
                HYMAP2_routing_struc(n)%drtotlgh,&
                HYMAP2_routing_struc(n)%drnoutlet,&
                !HYMAP2_routing_struc(n)%droutlet,&
                HYMAP2_routing_struc(n)%drstomax,&
                HYMAP2_routing_struc(n)%drout(:,m),&
                HYMAP2_routing_struc(n)%drsto(:,m)) 
        
           rnfsto_mm(:,m)=1000*HYMAP2_routing_struc(n)%rnfsto(:,m)/&
                HYMAP2_routing_struc(n)%grarea
           bsfsto_mm(:,m)=1000*HYMAP2_routing_struc(n)%bsfsto(:,m)/&
                HYMAP2_routing_struc(n)%grarea
           
! The following calls are done because the vector size in LIS and HYMAP are different 
! (and because the LIS I/O system is being leveraged for output

           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%rivsto(:,m),&
                rivsto_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%rivdph(:,m),&
                rivdph_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%rivvel(:,m),&
                rivvel_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%rivout(:,m),&
                rivout_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%evpout(:,m),&
                evpout_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%fldout(:,m),&
                fldout_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%fldsto(:,m),&
                fldsto_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%flddph(:,m),&
                flddph_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%fldvel(:,m),&
                fldvel_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%fldfrc(:,m),&
                fldfrc_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%fldare(:,m),&
                fldare_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%sfcelv(:,m),&
                sfcelv_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%rnfsto(:,m),&
                rnfsto_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%bsfsto(:,m),&
                bsfsto_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%rnfdwi(:,m),&
                rnfdwi_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%bsfdwi(:,m),&
                bsfdwi_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%surfws(:,m),&
                surfws_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%ewat(:,m),&
                ewat_lvec)
           call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%edif(:,m),&
                edif_lvec)
           !ag(27Apr2020)
           if(HYMAP2_routing_struc(n)%flowtype==4)then
             call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%drsto(:,m),&
                  drsto_lvec)
             call HYMAP2_grid2tile(n,m,HYMAP2_routing_struc(n)%drout(:,m),&
                  drout_lvec)
          endif
        enddo

     else !single member run

        if(LIS_rc%lsm.ne."none") then 
           tmpr=0.
           tmpb=0.
           
           !import surface runoff and baseflow
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
     
           !temporary solution  
           call LIS_tile2grid(n,tmpr,surface_runoff_t)
           call LIS_tile2grid(n,tmpb,baseflow_t)
           
           call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmpr,surface_runoff)
           call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmpb,baseflow)
        else

           call readrunoffdata(trim(LIS_rc%runoffdatasource)//char(0),&
                n,tmpr,tmpb)

           call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmpr,surface_runoff)
           call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
                HYMAP2_routing_struc(n)%nseqall,&
                HYMAP2_routing_struc(n)%imis,&
                HYMAP2_routing_struc(n)%seqx,&
                HYMAP2_routing_struc(n)%seqy,tmpb,baseflow)           
           
        endif
        
        call HYMAP2_model(n,real(HYMAP2_routing_struc(n)%imis),&
             LIS_rc%lnc(n),&
             LIS_rc%lnr(n),&
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
             HYMAP2_routing_struc(n)%rivout_pre(:,1),&
             HYMAP2_routing_struc(n)%rivdph_pre(:,1),&
             HYMAP2_routing_struc(n)%fldout_pre(:,1),&
             HYMAP2_routing_struc(n)%flddph_pre(:,1),&
             HYMAP2_routing_struc(n)%fldelv1(:,1),&
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
             HYMAP2_routing_struc(n)%fldgrd(:,:,1),&
             HYMAP2_routing_struc(n)%fldman,&
             HYMAP2_routing_struc(n)%fldhgt,&
             HYMAP2_routing_struc(n)%fldstomax(:,:,1),&
             HYMAP2_routing_struc(n)%rivman(:,1),&
             HYMAP2_routing_struc(n)%rivelv,&
             HYMAP2_routing_struc(n)%rivstomax,&
             HYMAP2_routing_struc(n)%rivlen,&
             HYMAP2_routing_struc(n)%rivwth(:,1),&
             HYMAP2_routing_struc(n)%rivhgt(:,1),&
             HYMAP2_routing_struc(n)%rivare(:,1),&
             HYMAP2_routing_struc(n)%rslpmin,&
             HYMAP2_routing_struc(n)%trnoff,&
             HYMAP2_routing_struc(n)%tbsflw,&
             HYMAP2_routing_struc(n)%cntime,&
             HYMAP2_routing_struc(n)%dwiflag,&
             HYMAP2_routing_struc(n)%rnfdwi_ratio,&
             HYMAP2_routing_struc(n)%bsfdwi_ratio,&
             surface_runoff,&
             baseflow,&
             HYMAP2_routing_struc(n)%edif(:,1),&
             HYMAP2_routing_struc(n)%rivsto(:,1),&
             HYMAP2_routing_struc(n)%rivdph(:,1),&
             HYMAP2_routing_struc(n)%rivvel(:,1),&
             HYMAP2_routing_struc(n)%rivout(:,1),&
             HYMAP2_routing_struc(n)%evpout(:,1),&
             HYMAP2_routing_struc(n)%fldout(:,1),&
             HYMAP2_routing_struc(n)%fldsto(:,1),&
             HYMAP2_routing_struc(n)%flddph(:,1),&
             HYMAP2_routing_struc(n)%fldvel(:,1),&
             HYMAP2_routing_struc(n)%fldfrc(:,1),&
             HYMAP2_routing_struc(n)%fldare(:,1),&
             HYMAP2_routing_struc(n)%sfcelv(:,1),&
             HYMAP2_routing_struc(n)%rnfsto(:,1),&
             HYMAP2_routing_struc(n)%bsfsto(:,1),&
             HYMAP2_routing_struc(n)%rnfdwi(:,1),&
             HYMAP2_routing_struc(n)%bsfdwi(:,1),&
             HYMAP2_routing_struc(n)%surfws(:,1),&
             HYMAP2_routing_struc(n)%dtaout(:,1),&
             !ag (27Apr2020)
             !urban drainage variables/parameters
             HYMAP2_routing_struc(n)%flowtype,&
             HYMAP2_routing_struc(n)%drvel,&
             HYMAP2_routing_struc(n)%drtotwth,&
             HYMAP2_routing_struc(n)%drrad,&
             HYMAP2_routing_struc(n)%drman,&
             HYMAP2_routing_struc(n)%drslp,&
             HYMAP2_routing_struc(n)%drtotlgh,&
             HYMAP2_routing_struc(n)%drnoutlet,&
             !HYMAP2_routing_struc(n)%droutlet,&
             HYMAP2_routing_struc(n)%drstomax,&
             HYMAP2_routing_struc(n)%drout,&
             HYMAP2_routing_struc(n)%drsto)              
              
        rnfsto_mm(:,1)=1000*HYMAP2_routing_struc(n)%rnfsto(:,1)/&
             HYMAP2_routing_struc(n)%grarea
        bsfsto_mm(:,1)=1000*HYMAP2_routing_struc(n)%bsfsto(:,1)/&
             HYMAP2_routing_struc(n)%grarea

        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%rivsto(:,1),&
             rivsto_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%rivdph(:,1),&
                rivdph_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%rivvel(:,1),&
             rivvel_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%rivout(:,1),&
             rivout_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%evpout(:,1),&
             evpout_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%fldout(:,1),&
             fldout_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%fldsto(:,1),&
             fldsto_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%flddph(:,1),&
             flddph_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%fldvel(:,1),&
             fldvel_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%fldfrc(:,1),&
             fldfrc_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%fldare(:,1),&
             fldare_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%sfcelv(:,1),&
             sfcelv_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%rnfsto(:,1),&
             rnfsto_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%bsfsto(:,1),&
             bsfsto_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%rnfdwi(:,1),&
             rnfdwi_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%bsfdwi(:,1),&
             bsfdwi_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%surfws(:,1),&
             surfws_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%ewat(:,1),&
             ewat_lvec)
        call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%edif(:,1),&
             edif_lvec)
        !ag(27Apr2020)
        if(HYMAP2_routing_struc(n)%flowtype==4)then
           call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%drsto(:,1),&
                drsto_lvec)
           call HYMAP2_grid2tile_noens(n,HYMAP2_routing_struc(n)%drout(:,1),&
                drout_lvec)
        endif

        !ag (12Sep2019)
        if (HYMAP2_routing_struc(n)%enable2waycpl==1) then
           allocate(tmp_nensem(LIS_rc%lnc(n),LIS_rc%lnr(n),1))
          ! River Storage
          call ESMF_StateGet(LIS_runoff_state(n),"River Storage",rivsto_field, rc=status)
          call LIS_verify(status, "HYMAP2_routing_run: ESMF_StateGet failed for River Storage")

          call ESMF_FieldGet(rivsto_field,localDE=0,farrayPtr=rivstotmp_lvec,rc=status)
          call LIS_verify(status)

          ! convert from m3 to m 
          do i=1,HYMAP2_routing_struc(n)%nseqall
             HYMAP2_routing_struc(n)%rivstotmp(i,:)=HYMAP2_routing_struc(n)%rivsto(i,:)/&
                  HYMAP2_routing_struc(n)%grarea(i)
          enddo
          
          !HYMAP2_routing_struc(n)%rivstotmp=HYMAP2_routing_struc(n)%rivsto
          
          call HYMAP2_vector2grid(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
               HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1), &
               HYMAP2_routing_struc(n)%rivstotmp)

          call LIS_grid2tile(n,tmp_nensem(:,:,1),rivstotmp_lvec)

          ! Flood Storage
          call ESMF_StateGet(LIS_runoff_state(n),"Flood Storage",fldsto_field, rc=status)
          call LIS_verify(status, "HYMAP2_routing_run: ESMF_StateGet failed for Flood Storage")

          call ESMF_FieldGet(fldsto_field,localDE=0,farrayPtr=fldstotmp_lvec,rc=status)
          call LIS_verify(status)
 
          ! convert from m3 to m  
          do i=1,HYMAP2_routing_struc(n)%nseqall
             HYMAP2_routing_struc(n)%fldstotmp(i,:)=HYMAP2_routing_struc(n)%fldsto(i,:)/&
                  HYMAP2_routing_struc(n)%grarea(i)
          enddo


          call HYMAP2_vector2grid(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
               HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1), &
               HYMAP2_routing_struc(n)%fldstotmp)

          call LIS_grid2tile(n,tmp_nensem(:,:,1),fldstotmp_lvec)

          ! Flooded Fraction
          call ESMF_StateGet(LIS_runoff_state(n),"Flooded Fraction",fldfrc_field, rc=status)
          call LIS_verify(status, "HYMAP2_routing_run: ESMF_StateGet failed for Flooded Fraction")

          call ESMF_FieldGet(fldfrc_field,localDE=0,farrayPtr=fldfrctmp_lvec,rc=status)
          call LIS_verify(status)
          
          !create flooded fraction flags
          call HYMAP2_vector2grid(LIS_rc%lnc(n),LIS_rc%lnr(n),1,&
               HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,tmp_nensem(:,:,1), &
               HYMAP2_routing_struc(n)%fldfrc)

          call LIS_grid2tile(n,tmp_nensem(:,:,1),fldfrctmp1_lvec)

          do t=1,LIS_rc%ntiles(n)
             if(fldfrctmp1_lvec(t) >= HYMAP2_routing_struc(n)%fldfrc2waycpl) then 
                fldfrctmp_lvec(t)=1.
             else
                fldfrctmp_lvec(t)=0.
             endif
          enddo
          deallocate(tmp_nensem)
       endif

    endif

    do t=1, LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)
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
       
       call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_DRSTO,&
            value=drsto_lvec(t),vlevel=1,unit="m3",&  
            direction="-")
       
       call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_DROUT,&
            value=drout_lvec(t),vlevel=1,unit="m3 s-1",&
            direction="-")
       
    enddo
    
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
    
    !ag(27Apr2020)
    deallocate(drsto_lvec)
    deallocate(drout_lvec)
    
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
       deallocate(tmp_psurf)
       deallocate(tmp_wind)
       deallocate(tmp_qnet)

       deallocate(evap)
       deallocate(tair)
       deallocate(qair)
       deallocate(pres)
       deallocate(wind)
       deallocate(qnet)
       deallocate(tmpet)
    endif
    
 endif
end subroutine HYMAP2_routing_run

!BOP
! !ROUTINE: HYMAP2_grid2tile
! \label{HYMAP2_grid2tile}
!
! !INTERFACE:
!
! !INTERFACE:
  subroutine HYMAP2_grid2tile(n,m,tvar,tvar_ens)
! !USES:
    use LIS_coreMod

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    integer, intent(in) :: m
    real                :: tvar(LIS_rc%nroutinggrid(n))
    real                :: tvar_ens(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n))

! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: i,t,c,r

    do i=1,LIS_rc%nroutinggrid(n)
       t = m + (i-1)*LIS_rc%nensem(n)
       tvar_ens(t) = tvar(i)
    enddo

  end subroutine HYMAP2_grid2tile

  !BOP
! !ROUTINE: HYMAP2_grid2tile_noens
! \label{HYMAP2_grid2tile_noens}
!
! !INTERFACE:
!
! !INTERFACE:
  subroutine HYMAP2_grid2tile_noens(n,tvar,tvar_ens)
! !USES:
    use LIS_coreMod

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    real                :: tvar(LIS_rc%nroutinggrid(n))
    real                :: tvar_ens(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n))

! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: i,t,c,r,m

    do i=1,LIS_rc%nroutinggrid(n)
       do m=1,LIS_rc%nensem(n)
          t = m + (i-1)*LIS_rc%nensem(n)
          tvar_ens(t) = tvar(i)
       enddo
    enddo
    
  end subroutine HYMAP2_grid2tile_noens
