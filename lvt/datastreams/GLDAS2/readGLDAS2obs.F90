!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGLDAS2Obs
! \label{readGLDAS2Obs}
!
! !INTERFACE: 
subroutine readGLDAS2Obs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use GLDAS2obsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This plugin processes the Global Land Data Assimilation System (GLDAS)
!   version 2 data available from NASA GES-DISC.
!   
!  NOTES: 
!  Currently the NOAH model-based monthly outputs in NetCDF format is 
!  is supported. The data can be downloaded from: 
!  http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP

  integer                :: c,r, k,nc,nr,tindex
  integer                :: nDims, nVars, nAtts, unlimId
  character(len=100)     :: var_name
  character(len=10)      :: var_suffix
  integer                :: flag
  integer                :: ftn
  character*100          :: fname
  logical                :: file_exists
  integer                :: qsid, qsbid, canopintid, lwdownid, swdownid
  integer                :: rainfid, snowfid, qgid, avgsurftid, qhid, qleid
  integer                :: lwnetid, swnetid, qsmid, smid, tsoilid, sweid, evapid
  integer                :: qairid, windid, tairid, psurfid
  real                   :: qs(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: qsb(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: canopint(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: lwdown(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: swdown(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: rainf(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: snowf(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: qg(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: avgsurft(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: qh(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: qle(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: lwnet(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: swnet(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: qsm(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: sm(gldas2obs(source)%nc, gldas2obs(source)%nr,4)
  real                   :: tsoil(gldas2obs(source)%nc, gldas2obs(source)%nr,4)
  real                   :: swe(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: evap(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: qair(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: wind(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: tair(gldas2obs(source)%nc, gldas2obs(source)%nr)
  real                   :: psurf(gldas2obs(source)%nc, gldas2obs(source)%nr)

  real                   :: qs_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: qsb_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: canopint_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: lwdown_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: swdown_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: rainf_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: snowf_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: qg_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: avgsurft_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: qh_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: qle_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: lwnet_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: swnet_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: qsm_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: sm_ip(LVT_rc%lnc,LVT_rc%lnr,4)
  real                   :: tsoil_ip(LVT_rc%lnc,LVT_rc%lnr,4)
  real                   :: swe_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: evap_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: qair_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: wind_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: tair_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: psurf_ip(LVT_rc%lnc,LVT_rc%lnr)

  integer                :: iret

  qs_ip       =  LVT_rc%udef
  qsb_ip      =  LVT_rc%udef
  canopint_ip =  LVT_rc%udef
  lwdown_ip   =  LVT_rc%udef
  swdown_ip   =  LVT_rc%udef
  rainf_ip    =  LVT_rc%udef
  snowf_ip    =  LVT_rc%udef
  qg_ip       =  LVT_rc%udef
  avgsurft_ip =  LVT_rc%udef
  qh_ip       =  LVT_rc%udef
  qle_ip      =  LVT_rc%udef
  lwnet_ip    =  LVT_rc%udef
  swnet_ip    =  LVT_rc%udef
  qsm_ip      =  LVT_rc%udef
  sm_ip       =  LVT_rc%udef
  tsoil_ip    =  LVT_rc%udef
  swe_ip      =  LVT_rc%udef
  evap_ip     =  LVT_rc%udef
  qair_ip     =  LVT_rc%udef
  wind_ip     =  LVT_rc%udef
  tair_ip     =  LVT_rc%udef
  psurf_ip    =  LVT_rc%udef

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  nc = GLDAS2obs(source)%nc
  nr = GLDAS2obs(source)%nr

  if((GLDAS2obs(source)%mo.ne.LVT_rc%d_nmo(source)).or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 

     if(GLDAS2obs(source)%startFlag) then 
        GLDAS2obs(source)%startFlag = .false. 
     endif

     GLDAS2obs(source)%yr = LVT_rc%d_nyr(source)
     GLDAS2obs(source)%mo = LVT_rc%d_nmo(source)

     call create_GLDAS2_filename(GLDAS2obs(source)%odir,&
          GLDAS2obs(source)%model_name, &
          LVT_rc%dyr(source),&
          LVT_rc%dmo(source),&
          fname)
     
     inquire(file=trim(fname),exist=file_exists) 

     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading GLDAS2 file ',trim(fname)
        
        iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
             ncid = ftn)
        if(iret.eq.0) then 

           iret = nf90_inquire(ftn, nDims, nVars, nAtts, unlimId)
           iret = nf90_inquire_variable(ftn, nVars, var_name)
           var_suffix = var_name(index(var_name, 'ave')+3:len(var_name))
           
           call LVT_verify(nf90_inq_varid(ftn,"Qs_GDS0_SFC_ave"//trim(var_suffix),qsid),&
                'nf90_inq_varid failed for Qs_GDS0_SFC_ave'//trim(var_suffix))
           call LVT_verify(nf90_inq_varid(ftn,"Qsb_GDS0_SFC_ave"//trim(var_suffix),qsbid),&
                'nf90_inq_varid failed for Qsb_GDS0_SFC_ave'//trim(var_suffix))
           call LVT_verify(nf90_inq_varid(ftn,"Canopint_GDS0_SFC_ave"//trim(var_suffix),canopintid),&
                 'nf90_inq_varid failed for Canopint_GDS0_SFC_ave'//trim(var_suffix))
           call LVT_verify(nf90_inq_varid(ftn,"LWdown_GDS0_SFC_ave"//trim(var_suffix),lwdownid),&
                 'nf90_inq_varid failed for LWdown_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"SWdown_GDS0_SFC_ave"//trim(var_suffix),swdownid),&
                 'nf90_inq_varid failed for SWdown_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Rainf_GDS0_SFC_ave"//trim(var_suffix),rainfid),&
                 'nf90_inq_varid failed for Rainf_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Snowf_GDS0_SFC_ave"//trim(var_suffix),snowfid),&
                 'nf90_inq_varid failed for Snowf_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Qg_GDS0_SFC_ave"//trim(var_suffix),qgid),&
                 'nf90_inq_varid failed for Qg_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"AvgSurfT_GDS0_SFC_ave"//trim(var_suffix),avgsurftid),&
                 'nf90_inq_varid failed for AvgSurfT_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Qh_GDS0_SFC_ave"//trim(var_suffix),qhid),&
                 'nf90_inq_varid failed for Qh_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Qle_GDS0_SFC_ave"//trim(var_suffix),qleid),&
                 'nf90_inq_varid failed for Qle_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"LWnet_GDS0_SFC_ave"//trim(var_suffix),lwnetid),&
                 'nf90_inq_varid failed for LWnet_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"SWnet_GDS0_SFC_ave"//trim(var_suffix),swnetid),&
                 'nf90_inq_varid failed for SWnet_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Qsm_GDS0_SFC_ave"//trim(var_suffix),qsmid),&
                 'nf90_inq_varid failed for Qsm_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"SoilM_GDS0_DBLY_ave"//trim(var_suffix),smid),&
                 'nf90_inq_varid failed for SoilM_GDS0_DBLY_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"TSoil_GDS0_DBLY_ave"//trim(var_suffix),tsoilid),&
                 'nf90_inq_varid failed for TSoil_GDS0_DBLY_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"SWE_GDS0_SFC_ave"//trim(var_suffix),sweid),&
                 'nf90_inq_varid failed for SWE_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Evap_GDS0_SFC_ave"//trim(var_suffix),evapid),&
                 'nf90_inq_varid failed for Evap_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Qair_GDS0_SFC_ave"//trim(var_suffix),qairid),&
                 'nf90_inq_varid failed for Qair_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Wind_GDS0_SFC_ave"//trim(var_suffix),windid),&
                 'nf90_inq_varid failed for Wind_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"Tair_GDS0_SFC_ave"//trim(var_suffix),tairid),&
                 'nf90_inq_varid failed for Tair_GDS0_SFC_ave'//trim(var_suffix))                
           call LVT_verify(nf90_inq_varid(ftn,"PSurf_GDS0_SFC_ave"//trim(var_suffix),psurfid),&
                 'nf90_inq_varid failed for PSurf_GDS0_SFC_ave'//trim(var_suffix))                

           call LVT_verify(nf90_get_var(ftn,qsid,qs,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Qs_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,qsbid,qsb,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Qsb_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,canopintid,canopint,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Canopint_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,lwdownid,lwdown,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var LWdown_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,swdownid,swdown,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var SWdown_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,rainfid,rainf,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Rainf_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,snowfid,snowf,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Snowf_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,qgid,qg,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Qg_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,avgsurftid,avgsurft,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var AvgSurfT_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,qhid,qh,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Qh_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,qleid,qle,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Qle_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,lwnetid,lwnet,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var LWnet_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,swnetid,swnet,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var SWnet_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,qsmid,qsm,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Qsm_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,smid,sm,&
                start=(/1,1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr,4/)),&
                'Error in nf90_get_var SoilM_GDS0_DBLY_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,tsoilid,tsoil,&
                start=(/1,1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr,4/)),&
                'Error in nf90_get_var TSoil_GDS0_DBLY_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,sweid,swe,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var SWE_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,evapid,evap,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Evap_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,qairid,qair,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Qair_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,windid,wind,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Wind_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,tairid,tair,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var Tair_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_get_var(ftn,psurfid,psurf,&
                start=(/1,1/),&
                count=(/gldas2obs(source)%nc,gldas2obs(source)%nr/)),&
                'Error in nf90_get_var PSurf_GDS0_SFC_ave'//trim(var_suffix))

           call LVT_verify(nf90_close(ftn))

           call interp_gldas2var2d(source,qs,qs_ip)
           call interp_gldas2var2d(source,qsb,qsb_ip)
           call interp_gldas2var2d(source,canopint,canopint_ip)
           call interp_gldas2var2d(source,lwdown,lwdown_ip)
           call interp_gldas2var2d(source,swdown,swdown_ip)
           call interp_gldas2var2d(source,rainf,rainf_ip)
           call interp_gldas2var2d(source,snowf,snowf_ip)
           call interp_gldas2var2d(source,qg,qg_ip)
           call interp_gldas2var2d(source,avgsurft,avgsurft_ip)
           call interp_gldas2var2d(source,qh,qh_ip)
           call interp_gldas2var2d(source,qle,qle_ip)
           call interp_gldas2var2d(source,lwnet,lwnet_ip)
           call interp_gldas2var2d(source,swnet,swnet_ip)
           call interp_gldas2var2d(source,qsm,qsm_ip)
           call interp_gldas2var3d(source,4,sm,sm_ip)
           call interp_gldas2var3d(source,4,tsoil,tsoil_ip)
           call interp_gldas2var2d(source,swe,swe_ip)
           call interp_gldas2var2d(source,evap,evap_ip)
           call interp_gldas2var2d(source,qair,qair_ip)
           call interp_gldas2var2d(source,wind,wind_ip)
           call interp_gldas2var2d(source,tair,tair_ip)
           call interp_gldas2var2d(source,psurf,psurf_ip)

        endif
        
     endif
     
  endif
#endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_QS,source,qs_ip,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_QSB,source,qsb_ip,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_CANOPINT,source,canopint_ip,&
       vlevel=1,units="kg/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_LWDOWNFORC,source,lwdown_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_SWDOWNFORC,source,swdown_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_RAINFFORC,source,rainf_ip,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWFFORC,source,snowf_ip,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_QG,source,qg_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,avgsurft_ip,&
       vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,qle_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_QH,source,qh_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_LWNET,source,lwnet_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_SWNET,source,swnet_ip,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_QSM,source,qsm_ip,&
       vlevel=1,units="kg/m2s")
  do k=1,4
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(sm_ip(c,r,k).gt.0) then 
              if(k.eq.1) then 
                 sm_ip(c,r,k) = sm_ip(c,r,k)/100.0
              elseif(k.eq.2) then 
                 sm_ip(c,r,k) = sm_ip(c,r,k)/300.0
              elseif(k.eq.3) then 
                 sm_ip(c,r,k) = sm_ip(c,r,k)/600.0
              elseif(k.eq.4) then 
                 sm_ip(c,r,k) = sm_ip(c,r,k)/1000.0
              endif
           endif
        enddo
     enddo

     call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,sm_ip(:,:,k),&
          vlevel=k,units="m3/m3")
     call LVT_logSingleDataStreamVar(LVT_MOC_SOILTEMP,source,tsoil_ip(:,:,k),&
          vlevel=k,units="K")
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,swe_ip,&
       vlevel=1,units="kg/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_EVAP,source,evap_ip,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_QAIRFORC,source,qair_ip,&
       vlevel=1,units="kg/kg")
  call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,source,tair_ip,&
       vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_PSURFFORC,source,psurf_ip,&
       vlevel=1,units="Pa")
end subroutine readGLDAS2Obs


!BOP
!
! !ROUTINE: interp_gldas2var2d
! \label{interp_gldas2var2d}
!
! !INTERFACE: 
subroutine interp_gldas2var2d(source, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use GLDAS2obsMod
! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(gldas2obs(source)%nc,gldas2obs(source)%nr)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the GLDAS fields to the 
!  target LVT domain
!
!EOP

  real              :: var_inp_1d(gldas2obs(source)%nc*gldas2obs(source)%nr)
  logical*1         :: input_bitmap(gldas2obs(source)%nc*gldas2obs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r
  integer           :: iret

  nc = gldas2obs(source)%nc
  nr = gldas2obs(source)%nr
  
  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        if(var_inp(c,r).ne.1e20) then 
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
           input_bitmap(c+(r-1)*nc) = .true. 
        else
           var_inp(c,r) = LVT_rc%udef
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
        endif
     enddo
  enddo
  
  if(LVT_isAtAfinerResolution(gldas2obs(source)%datares)) then
     call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d, &
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          gldas2obs(source)%rlat, & 
          gldas2obs(source)%rlon, &
          gldas2obs(source)%n11, &
          LVT_rc%udef, iret)
     
  else
     call upscaleByAveraging(&
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          gldas2obs(source)%n11, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d)
     
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
           var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
        endif
     enddo
  enddo

end subroutine interp_gldas2var2d

!BOP
!
! !ROUTINE: interp_gldas2var3d
! \label{interp_gldas2var3d}
!
! !INTERFACE: 
subroutine interp_gldas2var3d(source, nlevs, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use GLDAS2obsMod
! !ARGUMENTS: 
  integer           :: source
  integer           :: nlevs
  real              :: var_inp(gldas2obs(source)%nc,gldas2obs(source)%nr,nlevs)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr,nlevs)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the GLDAS fields to the 
!  target LVT domain
!
!EOP
  integer           :: k
  real              :: var_inp_1d(gldas2obs(source)%nc*gldas2obs(source)%nr)
  logical*1         :: input_bitmap(gldas2obs(source)%nc*gldas2obs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r
  integer           :: iret

  do k=1,nlevs

     nc = gldas2obs(source)%nc
     nr = gldas2obs(source)%nr
     
     input_bitmap = .false. 
     do r=1,nr
        do c=1,nc
           if(var_inp(c,r,k).ne.1e20) then 
              var_inp_1d(c+(r-1)*nc) = var_inp(c,r,k)
              input_bitmap(c+(r-1)*nc) = .true. 
           else
              var_inp(c,r,k) = LVT_rc%udef
              var_inp_1d(c+(r-1)*nc) = var_inp(c,r,k)
           endif
        enddo
     enddo
     
     if(LVT_isAtAfinerResolution(gldas2obs(source)%datares)) then
        call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
             var_inp_1d, output_bitmap, var_out_1d, &
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, &
             gldas2obs(source)%rlat, & 
             gldas2obs(source)%rlon, &
             gldas2obs(source)%n11, &
             LVT_rc%udef, iret)
        
     else
        call upscaleByAveraging(&
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
             gldas2obs(source)%n11, input_bitmap, &
             var_inp_1d, output_bitmap, var_out_1d)
        
     endif
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
              var_out(c,r,k) = var_out_1d(c+(r-1)*LVT_rc%lnc)
           endif
        enddo
     enddo
  enddo
     
end subroutine interp_gldas2var3d

!BOP
! 
! !ROUTINE: create_GLDAS2_filename
! \label{create_GLDAS2_filename}
!
! !INTERFACE: 
subroutine create_GLDAS2_filename(odir,model_name, yr,mo,filename)
! 
! !USES:   
  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  character(len=*)             :: model_name
  integer                      :: yr
  integer                      :: mo
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the GLDAS2 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            GLDAS2 base directory
!   \item[model\_name]     name of the model used in the GLDAS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the GLDAS2 file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*2             :: fmo

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  
  filename = trim(odir)//'/GLDAS_'//trim(model_name)//'025_M.A'//&
       trim(fyr)//trim(fmo)//'.020.nc'
  
end subroutine create_GLDAS2_filename


