!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_DAobs_pluginMod
!BOP
!
! !MODULE: LIS_DAobs_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions that are
!   used to read observation data for data assimiliation.  
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  27 Feb 2005;   Sujay Kumar  Initial Specification
!  11 Aug 2016:   Mahdi Navari, PILDAS added 
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_DAobs_plugin
  PUBLIC :: LIS_isDAinstanceValid
  
  type DAobsFuncEntry

     character*100            :: name
     character*100            :: class
     type(DAobsFuncEntry), pointer :: next

  end type DAobsFuncEntry
  
  type, public :: daobsfunc_meta

     type(daobsfuncentry), pointer :: head_daobsfunc_list

  end type daobsfunc_meta
     
  type(daobsfunc_meta) :: LIS_DAobsFuncEntry     
contains
!BOP
! !ROUTINE: LIS_DAobs_plugin
!  \label{LIS_DAobs_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing routines to handle the 
!  observation data for assimilation.As explaind in 
!  the {\tt dataassim\_module}, there are three
!  different abstractions associated with the 
!  data assimilation implementation. The implementations defined in this
!  subroutine complete the "wirings" required to complete the 
!  observation data-related abstractions. Other required
!  interfaces defined in {\tt lsmda\_pluginMod} and {\tt dataassim\_pluginMod}
!  should be completed for a successful implementation. 
! 
!  The following interfaces should be implemented in this routine. 
!  \begin{description}
!  \item[Define DAobs class]
!      Define the "class" of DAobservation (LSM/Routing), which helps to
!      setup the appropriate surface model for data assimilation
!  \item[Setup]
!      Initialization of data and memory structures
!      (to be registered using {\tt registerreaddaobssetup} and later called 
!       through {\tt daobssetup})
!  \item[read observations]
!      Routines to read the observation data and perform any 
!      spatial transformation. 
!      (to be registered using {\tt registerreaddaobs} and later called
!       through {\tt readobservations})
!  \item[get number of selected observations]
!      routines to retrieve the number of selected observations for 
!      the selected modeling point. 
!      (to be registered using {\tt registergetnso} and later called
!      through {\tt getselctedobsnumber}
!   \end{description}
! 
!  The user-defined functions are included in the registry using a single 
!  index. For example, consider an instance where soil moisture
!  assimilation using TMI soil moisture data is conducted. The methods
!  should be defined in the registry as follows, if the index of the 
!  'assimilation set' (assimilating TMI to update soil moisture variables)
!  is defined to be 1
!
!  \begin{verbatim}
!    call registerdaobssetup(1,TMIsmobs_setup)  
!    call registerreaddaobs(1,read_TMIsm)
!    call registergetnso(1,getNSO_TMIsm)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as 
!   follows: 
!  
!  \begin{verbatim}
!    call daobssetup(1)       -  calls TMIsmobs_setup
!    call readobservations(1) -  calls read_TMIsm
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the 
!   following manner. 
!   \begin{verbatim}
!    call daobssetup(lis%daset, )       
!    call readobservations(lis%daset)
!   \end{verbatim}
!   where $lis\%daset$ is set through the configuration
!   utility, enabling the user to make selections at runtime. 
!
! !INTERFACE:
subroutine LIS_DAobs_plugin
!EOP
#if ( ( defined DA_DIRECT_INSERTION ) || \
      ( defined DA_ENKS )             || \
      ( defined DA_ENKF ) )

   use LIS_pluginIndices

#if ( defined DA_OBS_SYNTHETICSM )
   use syntheticsmobs_module,   only : syntheticsmobs_setup
#endif

#if ( defined DA_OBS_SYNTHETICWL )
   use syntheticwlobs_module,   only : syntheticwlobs_setup
#endif

#if ( defined DA_OBS_SYNTHETICSND )
   use syntheticsndobs_module,  only : syntheticsndobs_setup
#endif

#if ( defined DA_OBS_SYNTHETICSNOW )
   use syntheticsweobs_module,  only : syntheticsweobs_setup
#endif

#if ( defined DA_OBS_SYNTHETICSNOWTB )
   use syntheticSnowTbObs_Mod,  only : syntheticSnowTbobs_setup
#endif

#if 0 
   use syntheticlstobs_module,  only : syntheticlstobs_setup
   use multisynsmobs_Mod,       only : multisynsmobs_setup
   use ISCCP_Tskin_module,      only : ISCCP_Tskin_setup
   use MODISscaobs_module,      only : MODISscaobs_setup
#endif

#if ( defined DA_OBS_SNODEP )
   use SNODEPobs_Mod,           only : SNODEPobs_setup
#endif

#if ( defined DA_OBS_USAFSI )
   use USAFSIobs_Mod,           only : USAFSIobs_setup
#endif

#if 0 
   use NASA_AMSREsm_Mod,        only : NASA_AMSREsm_setup
#endif

#if ( defined DA_OBS_LPRM_AMSRESM )
   use LPRM_AMSREsm_Mod,        only : LPRM_AMSREsm_setup
#endif

#if ( defined DA_OBS_ESACCI_SM )
   use ESACCI_sm_Mod,           only : ESACCI_sm_setup
#endif

#if ( defined DA_OBS_SMOPSSM )
!   use SMOPSsm_Mod,             only : SMOPSsm_setup   
#endif

!MN
#if ( defined DA_OBS_SMOPS_ASCATSM )
   use SMOPS_ASCATsm_Mod,             only : SMOPS_ASCATsm_setup   
#endif

#if ( defined DA_OBS_SMOPS_SMOSSM )
   use SMOPS_SMOSsm_Mod,             only : SMOPS_SMOSsm_setup   
#endif

#if ( defined DA_OBS_SMOPS_AMSR2SM )
   use SMOPS_AMSR2sm_Mod,             only : SMOPS_AMSR2sm_setup   
#endif

#if ( defined DA_OBS_SMOPS_SMAPSM )
   use SMOPS_SMAPsm_Mod,             only : SMOPS_SMAPsm_setup   
#endif

#if 0
   use ANSASWEsnow_Mod,         only : ANSASWEsnow_setup
#endif

#if ( defined DA_OBS_ANSA_SCF )
   use ANSASCFsnow_Mod,         only : ANSASCFsnow_setup
#endif

#if ( defined DA_OBS_ANSA_SNWD )
   use ANSASNWDsnow_Mod,        only : ANSASNWDsnow_setup
#endif

#if 0  
   use AMSRE_SWE_Mod,           only : AMSRE_SWE_setup
!    use AMSRE_snow_Mod,          only : AMSRE_snow_setup !yliu
#endif

#if ( defined DA_OBS_PMW_SNOW )
   use PMW_snow_Mod,            only : PMW_snow_setup !yliu
#endif

#if ( defined DA_OBS_SMMR_SNWD )
   use SMMRSNWDsnow_Mod,        only : SMMRSNWDsnow_setup
#endif

#if ( defined DA_OBS_SSMI_SNWD )
   use SSMISNWDsnow_Mod,        only : SSMISNWDsnow_setup
#endif

#if ( defined DA_OBS_GCOMW_AMSR2L3SND )
   use GCOMW_AMSR2L3SND_Mod,    only : GCOMW_AMSR2L3SND_setup
#endif

#if ( defined DA_OBS_SMOS_NESDIS )
   use SMOSNESDISsm_Mod,        only : SMOSNESDISsm_setup
#endif

#if 0
    use WindSatsm_Mod,           only : WindSatsm_setup
    use WindSatCsm_Mod,          only : WindSatCsm_setup
    use simGRACEJPLobs_module,   only : simGRACEJPLobs_setup
    use SYN_LBAND_TB_Mod,        only : SYN_LBAND_TB_setup
    use ASCAT_TUWsm_Mod,         only : ASCAT_TUWsm_setup
    use GCOMW_AMSR2L3sm_Mod,     only : GCOMW_AMSR2L3sm_setup
    use SMOSL2sm_Mod,            only : SMOSL2sm_setup
#endif

#if ( defined DA_OBS_GRACE )
    use GRACEobs_module,         only : GRACEobs_setup
#endif

#if ( defined DA_OBS_PILDAS )
    use pildassmobs_module,      only : pildassmobs_setup
#endif

#if ( defined DA_OBS_NASA_SMAPSM )
    use NASASMAPsm_Mod,          only : NASASMAPsm_setup
#endif

#if ( defined DA_CDF_TRANSFER_NASA_SMAPSM )
    use cdfTransfer_NASASMAPsm_Mod,         only : cdfTransfer_NASASMAPsm_setup
#endif

!YK
#if ( defined DA_OBS_SMOS_NRT_NN )
    use SMOSNRTNNL2sm_Mod,       only : SMOSNRTNNL2sm_setup
#endif

!YK
#if ( defined DA_OBS_SMAP_E_OPL_SM )
    use SMAPEOPLsm_Mod,       only : SMAPEOPLsm_setup
#endif

#if ( defined DA_OBS_NASA_SMAPVOD )
    use NASASMAPvod_Mod,          only : NASASMAPvod_setup
#endif

#if ( defined DA_OBS_GLASS_LAI )
    use GLASSLAI_Mod,          only : GLASSlai_setup
#endif
#if ( defined DA_OBS_MCD15A2H_LAI )
    use MCD15A2HLAI_Mod,       only : MCD15A2Hlai_setup
#endif

!Y.Kwon
#if ( defined DA_OBS_VIIRS_GVF )
    use VIIRSgvf_Mod,       only : VIIRSgvf_setup
#endif
!Y.Kwon
#if ( defined DA_OBS_CDFS_GVF )
    use CDFSgvf_Mod,       only : CDFSgvf_setup
#endif

#if ( defined DA_OBS_NRT_SMAPSM )
    use SMAPNRTsm_Mod,           only : SMAPNRTsm_setup
#endif

#if ( defined DA_OBS_GLASS_Albedo )
    use GLASSAlbedo_Mod,       only : GLASSalbedo_setup
#endif

#if ( defined DA_OBS_MODISSPORT_LAI )
    use MODISsportLAI_Mod,          only : MODISsportLAI_setup
#endif

#if ( defined DA_OBS_ASO_SWE)
    use ASO_SWE_Mod,  only : ASO_SWE_setup
#endif

#if ( defined DA_OBS_THYSM)
    use THySM_Mod,    only : THySM_setup
#endif

#if ( defined DA_OBS_SNODAS )
   use SNODAS_Mod,    only : SNODAS_setup
#endif

#if ( defined DA_OBS_WUSUCLA )
   use WUS_UCLAsnowMod, only : WUS_UCLAsnow_setup
#endif

#if ( defined DA_OBS_HYDROWEBWL )
   use hydrowebWLobs_module,   only : hydrowebwlobs_setup
#endif
    
#if ( defined DA_OBS_SYNTHETICSM )
    external read_syntheticsmobs, write_syntheticsmobs
#endif

#if ( defined DA_OBS_SYNTHETICWL )
    external read_syntheticwlobs, write_syntheticwlobs
#endif

#if ( defined DA_OBS_THYSM )
    external read_THySM, write_THySM
#endif    

#if ( defined DA_OBS_SYNTHETICSND )
    external read_syntheticsndobs,write_syntheticsndobs
#endif

#if ( defined DA_OBS_SYNTHETICSNOW )
    external read_syntheticsweobs,write_syntheticsweobs
#endif

#if ( defined DA_OBS_SYNTHETICSNOWTB )
    external read_syntheticSnowTbObs,write_syntheticSnowTbObs
#endif


#if 0
   external read_syntheticlstobs
   external read_multisynsmobs
   external read_ISCCP_Tskin, write_ISCCP_Tskin
   external read_MODISscaobs, write_MODISsca
#endif

#if ( defined DA_OBS_SNODEP )
   external read_SNODEPobs, write_SNODEPobs
#endif

#if ( defined DA_OBS_USAFSI )
   external read_USAFSIobs, write_USAFSIobs
#endif

#if 0
   external read_NASA_AMSREsm, write_NASA_AMSREsmobs
#endif

#if ( defined DA_OBS_LPRM_AMSRESM )
   external read_LPRM_AMSREsm, write_LPRM_AMSREsmobs
#endif

#if ( defined DA_OBS_ESACCI_SM )
   external read_ESACCIsm, write_ESACCIsmobs
#endif


#if ( defined DA_OBS_SMOPSSM )
!   external read_SMOPSsm, write_SMOPSsmobs
#endif

!MN
#if ( defined DA_OBS_SMOPS_ASCATSM )
   external read_SMOPS_ASCATsm, write_SMOPS_ASCATsmobs
#endif

#if ( defined DA_OBS_SMOPS_SMOSSM )
   external read_SMOPS_SMOSsm, write_SMOPS_SMOSsmobs
#endif

#if ( defined DA_OBS_SMOPS_AMSR2SM )
   external read_SMOPS_AMSR2sm, write_SMOPS_AMSR2smobs
#endif

#if ( defined DA_OBS_SMOPS_SMAPSM )
   external read_SMOPS_SMAPsm, write_SMOPS_SMAPsmobs
#endif

#if ( defined DA_OBS_SMOS_NESDIS )
   external read_SMOSNESDISsm, write_SMOSNESDISsmobs
#endif

#if 0 
   external read_ANSASWEsnow, write_ANSASWEsnowobs
#endif

#if ( defined DA_OBS_ANSA_SCF )
   external read_ANSASCFsnow, write_ANSASCFsnowobs
#endif

#if ( defined DA_OBS_SMMR_SNWD )
   external read_SMMRSNWDsnow, write_SMMRSNWDsnowobs
#endif

#if ( defined DA_OBS_SSMI_SNWD )
   external read_SSMISNWDsnow, write_SSMISNWDsnowobs
#endif

#if ( defined DA_OBS_ANSA_SNWD )
   external read_ANSASNWDsnow, write_ANSASNWDsnowobs
#endif

#if ( defined DA_OBS_GCOMW_AMSR2L3SND )
   external read_GCOMW_AMSR2L3SND,  write_GCOMW_AMSR2L3sndobs
#endif

#if ( defined DA_OBS_HYDROWEBWL )
    external read_hydrowebWLobs, write_hydrowebWLobs
#endif

    
#if 0 
   external read_WindSatsm, write_WindSatsmobs
   external read_WindSatCsm, write_WindSatCsmobs
   external read_AMSRE_SWE, write_AMSRE_SWEobs
!    external read_AMSRE_snow, write_AMSRE_snowobs, getNSO_AMSRE_snow !yliu
#endif

#if ( defined DA_OBS_PMW_SNOW )
   external read_PMW_snow, write_PMW_snowobs
#endif

#if 0
    external read_simGRACEJPLobs, write_simGRACEJPLobs    
    external read_SYN_LBAND_TB, write_SYN_LBAND_TB
    external read_ASCAT_TUWsm, write_ASCAT_TUWsmobs
    external read_GCOMW_AMSR2L3sm,  write_GCOMW_AMSR2L3smobs
   external read_SMOSL2sm,write_SMOSL2smobs
#endif

#if ( defined DA_OBS_GRACE )
    external read_GRACEobs,  write_GRACEobs    
#endif

#if ( defined DA_OBS_PILDAS )
    external read_pildassmobs, write_pildassmobs
#endif

#if ( defined DA_OBS_NASA_SMAPSM )
    external read_NASASMAPsm, write_NASASMAPsmobs
#endif

#if ( defined DA_CDF_TRANSFER_NASA_SMAPSM )
    external read_cdfTransfer_NASASMAPsm, write_cdfTransfer_NASASMAPsmobs
#endif

!YK
#if ( defined DA_OBS_SMOS_NRT_NN )
    external read_SMOSNRTNNL2sm, write_SMOSNRTNNL2smobs
#endif

!YK
#if ( defined DA_OBS_SMAP_E_OPL_SM)
    external read_SMAPEOPLsm, write_SMAPEOPLsmobs
#endif

#if ( defined DA_OBS_NASA_SMAPVOD)
    external read_NASASMAPvod, write_NASASMAPvodobs
#endif

#if ( defined DA_OBS_GLASS_LAI)
    external read_GLASSlai, write_GLASSlai
#endif

#if ( defined DA_OBS_MCD15A2H_LAI)
    external read_MCD15A2Hlai, write_MCD15A2Hlai
#endif

#if ( defined DA_OBS_GLASS_Albedo)
    external read_GLASSalbedo, write_GLASSalbedo
#endif

#if ( defined DA_OBS_MODISSPORT_LAI )
    external read_MODISsportLAI, write_MODISsportLAI
#endif

!Y.Kwon
#if ( defined DA_OBS_VIIRS_GVF )
    external read_VIIRSgvf, write_VIIRSgvfobs
#endif
!Y.Kwon
#if ( defined DA_OBS_CDFS_GVF )
    external read_CDFSgvf, write_CDFSgvfobs
#endif

#if ( defined DA_OBS_NRT_SMAPSM )
    external read_SMAPNRTsm, write_SMAPNRTsmobs
#endif

#if ( defined DA_OBS_ASO_SWE)
    external read_ASO_SWE, write_ASO_SWEobs
#endif

#if ( defined DA_OBS_SNODAS )
   external read_SNODAS,  write_SNODAS
#endif

#if ( defined DA_OBS_WUSUCLA )    
   external read_WUS_UCLAsnow, write_WUS_UCLAsnow
#endif
   
    LIS_DAobsFuncEntry%head_daobsfunc_list => null()
    
#if ( defined DA_OBS_SYNTHETICSM )
!synthetic noah soil moisture    
   call registerdaobsclass(trim(LIS_synsmId),"LSM")
   call registerdaobssetup(trim(LIS_synsmId)//char(0),syntheticsmobs_setup)
   call registerreaddaobs(trim(LIS_synsmId)//char(0),read_syntheticsmobs)
   call registerwritedaobs(trim(LIS_synsmId)//char(0),write_syntheticsmobs)
#endif

#if ( defined DA_OBS_SYNTHETICWL )
!synthetic noah soil moisture    
   call registerdaobsclass(trim(LIS_synwlId),"Routing")
   call registerdaobssetup(trim(LIS_synwlId)//char(0),syntheticwlobs_setup)
   call registerreaddaobs(trim(LIS_synwlId)//char(0),read_syntheticwlobs)
   call registerwritedaobs(trim(LIS_synwlId)//char(0),write_syntheticwlobs)
#endif

#if ( defined DA_OBS_SYNTHETICSND )
   call registerdaobsclass(trim(LIS_synsndId),"LSM")
   call registerdaobssetup(trim(LIS_synsndId)//char(0),syntheticsndobs_setup)
   call registerreaddaobs(trim(LIS_synsndId)//char(0),read_syntheticsndobs)
   call registerwritedaobs(trim(LIS_synsndId)//char(0),write_syntheticsndobs)
#endif

#if ( defined DA_OBS_SYNTHETICSNOW )
   call registerdaobsclass(trim(LIS_synsweId),"LSM")
   call registerdaobssetup(trim(LIS_synsweId)//char(0),syntheticsweobs_setup)
   call registerreaddaobs(trim(LIS_synsweId)//char(0),read_syntheticsweobs)
   call registerwritedaobs(trim(LIS_synsweId)//char(0),write_syntheticsweobs)
#endif

#if ( defined DA_OBS_SYNTHETICSNOWTB )
   call registerdaobsclass(trim(LIS_synSnowTBId),"LSM")
   call registerdaobssetup(trim(LIS_synSnowTBId)//char(0),syntheticSnowTBobs_setup)
   call registerreaddaobs(trim(LIS_synSnowTBId)//char(0),read_syntheticSnowTBobs)
   call registerwritedaobs(trim(LIS_synSnowTBId)//char(0),write_syntheticSnowTBobs)
#endif

#if 0 
!synthetic lst
   call registerdaobsclass(trim(LIS_synlstId),"LSM")
   call registerdaobssetup(trim(LIS_synlstId)//char(0),syntheticlstobs_setup)
   call registerreaddaobs(trim(LIS_synlstId)//char(0),read_syntheticlstobs)

!multi layer synthetic sm obs 
   call registerdaobsclass(trim(LIS_multisynsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_multisynsmobsId)//char(0),multisynsmobs_setup)
   call registerreaddaobs(trim(LIS_multisynsmobsId)//char(0),read_multisynsmobs)


!ISSP Tskin 
   call registerdaobsclass(trim(LIS_isccpTskinId),"LSM")
   call registerdaobssetup(trim(LIS_isccpTskinId)//char(0),ISCCP_Tskin_setup)
   call registerreaddaobs(trim(LIS_isccpTskinId)//char(0),read_ISCCP_Tskin)
   call registerwritedaobs(trim(LIS_isccpTskinId)//char(0),write_ISCCP_Tskin)

!MODIS snow cover fraction
   call registerdaobsclass(trim(LIS_modisscfId),"LSM")
   call registerdaobssetup(trim(LIS_modisscfId)//char(0),MODISscaobs_setup)
   call registerreaddaobs(trim(LIS_modisscfId)//char(0),read_MODISscaobs)
   call registerwritedaobs(trim(LIS_modisscfId)//char(0),write_MODISsca)
#endif

#if ( defined DA_OBS_SNODEP )
!SNODEP obs 
   call registerdaobsclass(trim(LIS_snodepobsId),"LSM")
   call registerdaobssetup(trim(LIS_snodepobsId)//char(0),SNODEPobs_setup)
   call registerreaddaobs(trim(LIS_snodepobsId)//char(0),read_SNODEPobs)
   call registerwritedaobs(trim(LIS_snodepobsId)//char(0),write_SNODEPobs)
#endif


#if ( defined DA_OBS_USAFSI )
!USAFSI obs 
   call registerdaobsclass(trim(LIS_usafsiobsId),"LSM")
   call registerdaobssetup(trim(LIS_usafsiobsId)//char(0),USAFSIobs_setup)
   call registerreaddaobs(trim(LIS_usafsiobsId)//char(0),read_USAFSIobs)
   call registerwritedaobs(trim(LIS_usafsiobsId)//char(0),write_USAFSIobs)

#endif

#if 0
!NASA AMSRE obs 
   call registerdaobsclass(trim(LIS_NASA_AMSREsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_NASA_AMSREsmobsId)//char(0),NASA_AMSREsm_setup)
   call registerreaddaobs(trim(LIS_NASA_AMSREsmobsId)//char(0),read_NASA_AMSREsm)
   call registerwritedaobs(trim(LIS_NASA_AMSREsmobsId)//char(0),write_NASA_AMSREsmobs)
#endif

#if ( defined DA_OBS_LPRM_AMSRESM )
!LPRM AMSRE obs 
   call registerdaobsclass(trim(LIS_LPRM_AMSREsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_LPRM_AMSREsmobsId)//char(0), &
        LPRM_AMSREsm_setup)
   call registerreaddaobs(trim(LIS_LPRM_AMSREsmobsId)//char(0),  &
        read_LPRM_AMSREsm)
   call registerwritedaobs(trim(LIS_LPRM_AMSREsmobsId)//char(0), &
        write_LPRM_AMSREsmobs)
#endif

#if ( defined DA_OBS_ESACCI_SM )
!ESACCI sm obs
   call registerdaobsclass(trim(LIS_ESACCIsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_ESACCIsmobsId)//char(0),ESACCI_sm_setup)
   call registerreaddaobs(trim(LIS_ESACCIsmobsId)//char(0),read_ESACCIsm)
   call registerwritedaobs(trim(LIS_ESACCIsmobsId)//char(0),write_ESACCIsmobs)
#endif

#if 0
!ANSA SWE snow obs 
   call registerdaobsclass(trim(LIS_ANSASWEsnowobsId),"LSM")
   call registerdaobssetup(trim(LIS_ANSASWEsnowobsId)//char(0), &
        ANSASWEsnow_setup)
   call registerreaddaobs(trim(LIS_ANSASWEsnowobsId)//char(0),  &
        read_ANSASWEsnow)
   call registerwritedaobs(trim(LIS_ANSASWEsnowobsId)//char(0), &
        write_ANSASWEsnowobs)
#endif

#if ( defined DA_OBS_ANSA_SCF )
!ANSA SCF snow obs 
   call registerdaobsclass(trim(LIS_ANSASCFsnowobsId),"LSM")
   call registerdaobssetup(trim(LIS_ANSASCFsnowobsId)//char(0), &
        ANSASCFsnow_setup)
   call registerreaddaobs(trim(LIS_ANSASCFsnowobsId)//char(0),  &
        read_ANSASCFsnow)
   call registerwritedaobs(trim(LIS_ANSASCFsnowobsId)//char(0), &
        write_ANSASCFsnowobs)
#endif

#if ( defined DA_OBS_GCOMW_AMSR2L3SND )
!GCOMW AMSR2 L3 snow depth
   call registerdaobsclass(trim(LIS_GCOMW_AMSR2L3sndobsId),"LSM")
   call registerdaobssetup(trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0), &
        GCOMW_AMSR2L3snd_setup)
   call registerreaddaobs(trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),  &
        read_GCOMW_AMSR2L3snd)
   call registerwritedaobs(trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0), &
        write_GCOMW_AMSR2L3sndobs)
#endif

#if ( defined DA_OBS_ANSA_SNWD )
!ANSA SNWD snow obs 
   call registerdaobsclass(trim(LIS_ANSASNWDsnowobsId),"LSM")
   call registerdaobssetup(trim(LIS_ANSASNWDsnowobsId)//char(0), &
        ANSASNWDsnow_setup)
   call registerreaddaobs(trim(LIS_ANSASNWDsnowobsId)//char(0),  &
        read_ANSASNWDsnow)
   call registerwritedaobs(trim(LIS_ANSASNWDsnowobsId)//char(0), &
        write_ANSASNWDsnowobs)
#endif

#if ( defined DA_OBS_SMMR_SNWD )
!SMMR SNWD snow obs 
   call registerdaobsclass(trim(LIS_SMMRSNWDsnowobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMMRSNWDsnowobsId)//char(0), &
        SMMRSNWDsnow_setup)
   call registerreaddaobs(trim(LIS_SMMRSNWDsnowobsId)//char(0),  &
        read_SMMRSNWDsnow)
   call registerwritedaobs(trim(LIS_SMMRSNWDsnowobsId)//char(0), &
        write_SMMRSNWDsnowobs)
#endif

#if ( defined DA_OBS_SSMI_SNWD )
!SSMI SNWD snow obs 
   call registerdaobsclass(trim(LIS_SSMISNWDsnowobsId),"LSM")
   call registerdaobssetup(trim(LIS_SSMISNWDsnowobsId)//char(0), &
        SSMISNWDsnow_setup)
   call registerreaddaobs(trim(LIS_SSMISNWDsnowobsId)//char(0),  &
        read_SSMISNWDsnow)
   call registerwritedaobs(trim(LIS_SSMISNWDsnowobsId)//char(0), &
        write_SSMISNWDsnowobs)
#endif

#if 0

!AMSRE SWE
   call registerdaobsclass(trim(LIS_AMSREsweobsId),"LSM")
   call registerdaobssetup(trim(LIS_AMSREsweobsId)//char(0),AMSRE_SWE_setup)
   call registerreaddaobs(trim(LIS_AMSREsweobsId)//char(0),read_AMSRE_SWE)
   call registerwritedaobs(trim(LIS_AMSREsweobsId)//char(0),write_AMSRE_SWEobs)

!AMSRE Snow Depth or SWE, yliu
!    call registerdaobssetup(trim(LIS_AMSREsnowobsId)//char(0), &
!         AMSRE_snow_setup) !yliu
!    call registerreaddaobs(trim(LIS_AMSREsnowobsId)//char(0),  &
!         read_AMSRE_snow)   !yliu
!    call registerwritedaobs(trim(LIS_AMSREsnowobsId)//char(0), &
!         write_AMSRE_snowobs) !yliu

#endif

#if ( defined DA_OBS_PMW_SNOW )
!PMW Snow Depth or SWE, yliu
   call registerdaobsclass(trim(LIS_PMWsnowobsId),"LSM")
   call registerdaobssetup(trim(LIS_PMWsnowobsId)//char(0),PMW_snow_setup)
   call registerreaddaobs(trim(LIS_PMWsnowobsId)//char(0),read_PMW_snow)
   call registerwritedaobs(trim(LIS_PMWsnowobsId)//char(0),write_PMW_snowobs)
#endif

#if 0 
!WindSat soil moisture obs 
   call registerdaobsclass(trim(LIS_WindSatsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_WindSatsmobsId)//char(0),WindSatsm_setup)
   call registerreaddaobs(trim(LIS_WindSatsmobsId)//char(0),read_WindSatsm)
   call registerwritedaobs(trim(LIS_WindSatsmobsId)//char(0),write_WindSatsmobs)

!WindSat Cband (over Toulouse only) soil moisture obs 
   call registerdaobsclass(trim(LIS_WindSatCsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_WindSatCsmobsId)//char(0), &
        WindSatCsm_setup)
   call registerreaddaobs(trim(LIS_WindSatCsmobsId)//char(0),  &
        read_WindSatCsm)
   call registerwritedaobs(trim(LIS_WindSatCsmobsId)//char(0), &
        write_WindSatCsmobs)

!simulated GRACE
   call registerdaobsclass(trim(LIS_simGRACEJPLobsId),"LSM")
   call registerdaobssetup(trim(LIS_simGRACEJPLobsId)//char(0), &
        simGRACEJPLobs_setup)
   call registerreaddaobs(trim(LIS_simGRACEJPLobsId)//char(0),  &
        read_simGRACEJPLobs)
   call registerwritedaobs(trim(LIS_simGRACEJPLobsId)//char(0), &
        write_simGRACEJPLobs)

!Synthetic L-band Tb
   call registerdaobsclass(trim(LIS_synLbandTbobsId),"LSM")
   call registerdaobssetup(trim(LIS_synLbandTbobsId)//char(0), &
        SYN_LBAND_TB_setup)
   call registerreaddaobs(trim(LIS_synLbandTbobsId)//char(0),  &
        read_SYN_LBAND_TB)
   call registerwritedaobs(trim(LIS_synLbandTbobsId)//char(0), &
        write_SYN_LBAND_TB)
#endif

! MN : SMOPS soil moisture ! delete 
#if ( defined DA_OBS_SMOPSSM )
!   call registerdaobssetup(trim(LIS_SMOPSsmobsId)//char(0),SMOPSsm_setup)
!   call registerreaddaobs(trim(LIS_SMOPSsmobsId)//char(0),read_SMOPSsm)
!   call registerwritedaobs(trim(LIS_SMOPSsmobsId)//char(0),write_SMOPSsmobs)
#endif

! MN : SMOPS ASCAT soil moisture 
#if ( defined DA_OBS_SMOPS_ASCATSM )
   call registerdaobsclass(trim(LIS_SMOPS_ASCATsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMOPS_ASCATsmobsId)//char(0),SMOPS_ASCATsm_setup)
   call registerreaddaobs(trim(LIS_SMOPS_ASCATsmobsId)//char(0),read_SMOPS_ASCATsm)
   call registerwritedaobs(trim(LIS_SMOPS_ASCATsmobsId)//char(0),write_SMOPS_ASCATsmobs)
#endif

! MN : SMOPS SMOS soil moisture 
#if ( defined DA_OBS_SMOPS_SMOSSM )
   call registerdaobsclass(trim(LIS_SMOPS_SMOSsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMOPS_SMOSsmobsId)//char(0),SMOPS_SMOSsm_setup)
   call registerreaddaobs(trim(LIS_SMOPS_SMOSsmobsId)//char(0),read_SMOPS_SMOSsm)
   call registerwritedaobs(trim(LIS_SMOPS_SMOSsmobsId)//char(0),write_SMOPS_SMOSsmobs)
#endif

! MN : SMOPS AMSR2 soil moisture 
#if ( defined DA_OBS_SMOPS_AMSR2SM )
   call registerdaobsclass(trim(LIS_SMOPS_AMSR2smobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMOPS_AMSR2smobsId)//char(0),SMOPS_AMSR2sm_setup)
   call registerreaddaobs(trim(LIS_SMOPS_AMSR2smobsId)//char(0),read_SMOPS_AMSR2sm)
   call registerwritedaobs(trim(LIS_SMOPS_AMSR2smobsId)//char(0),write_SMOPS_AMSR2smobs)
#endif

! MN : SMOPS SMAP soil moisture 
#if ( defined DA_OBS_SMOPS_SMAPSM )
   call registerdaobsclass(trim(LIS_SMOPS_SMAPsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMOPS_SMAPsmobsId)//char(0),SMOPS_SMAPsm_setup)
   call registerreaddaobs(trim(LIS_SMOPS_SMAPsmobsId)//char(0),read_SMOPS_SMAPsm)
   call registerwritedaobs(trim(LIS_SMOPS_SMAPsmobsId)//char(0),write_SMOPS_SMAPsmobs)
#endif

#if ( defined DA_OBS_SMOS_NESDIS )
   call registerdaobsclass(trim(LIS_SMOSNESDISsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMOSNESDISsmobsId)//char(0), &
        SMOSNESDISsm_setup)
   call registerreaddaobs(trim(LIS_SMOSNESDISsmobsId)//char(0),  &
        read_SMOSNESDISsm)
   call registerwritedaobs(trim(LIS_SMOSNESDISsmobsId)//char(0), &
        write_SMOSNESDISsmobs)
#endif

#if 0
   call registerdaobssetup(trim(LIS_ASCAT_TUWsmobsId)//char(0), &
        ASCAT_TUWsm_setup)
   call registerreaddaobs(trim(LIS_ASCAT_TUWsmobsId)//char(0),  &
        read_ASCAT_TUWsm)
   call registerwritedaobs(trim(LIS_ASCAT_TUWsmobsId)//char(0), &
        write_ASCAT_TUWsmobs)

   call registerdaobssetup(trim(LIS_GCOMW_AMSR2L3smobsId)//char(0), &
        GCOMW_AMSR2L3sm_setup)
   call registerreaddaobs(trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),  &
        read_GCOMW_AMSR2L3sm)
   call registerwritedaobs(trim(LIS_GCOMW_AMSR2L3smobsId)//char(0), &
        write_GCOMW_AMSR2L3smobs)


   call registerdaobssetup(trim(LIS_SMOSL2smobsId)//char(0),SMOSL2sm_setup)
   call registerreaddaobs(trim(LIS_SMOSL2smobsId)//char(0),read_SMOSL2sm)
   call registerwritedaobs(trim(LIS_SMOSL2smobsId)//char(0),write_SMOSL2smobs)
#endif

#if ( defined DA_OBS_GRACE )
! GRACE TWS obs 
   call registerdaobsclass(trim(LIS_GRACEtwsobsId),"LSM")
    call registerdaobssetup(trim(LIS_GRACEtwsobsId)//char(0),GRACEobs_setup)
    call registerreaddaobs(trim(LIS_GRACEtwsobsId)//char(0),read_GRACEobs)
    call registerwritedaobs(trim(LIS_GRACEtwsobsId)//char(0),write_GRACEobs)
#endif

#if ( defined DA_OBS_PILDAS )
!pildas (synthetic soil moisture)    
   call registerdaobsclass(trim(LIS_pildassmobsId),"LSM")
   call registerdaobssetup(trim(LIS_pildassmobsId)//char(0),pildassmobs_setup)
   call registerreaddaobs(trim(LIS_pildassmobsId)//char(0),read_pildassmobs)
   call registerwritedaobs(trim(LIS_pildassmobsId)//char(0),write_pildassmobs)
#endif

#if ( defined DA_OBS_NASA_SMAPSM )
   call registerdaobsclass(trim(LIS_NASASMAPsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_NASASMAPsmobsId)//char(0),&
        NASASMAPsm_setup)
   call registerreaddaobs(trim(LIS_NASASMAPsmobsId)//char(0),&
        read_NASASMAPsm)
   call registerwritedaobs(trim(LIS_NASASMAPsmobsId)//char(0),&
        write_NASASMAPsmobs)
#endif

#if ( defined DA_CDF_TRANSFER_NASA_SMAPSM )
   call registerdaobsclass(trim(LIS_CDFTRANSFERNASASMAPsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_CDFTRANSFERNASASMAPsmobsId)//char(0),&
        cdfTransfer_NASASMAPsm_setup)
   call registerreaddaobs(trim(LIS_CDFTRANSFERNASASMAPsmobsId)//char(0),&
        read_cdfTransfer_NASASMAPsm)
   call registerwritedaobs(trim(LIS_CDFTRANSFERNASASMAPsmobsId)//char(0),&
        write_cdfTransfer_NASASMAPsmobs)
#endif

!YK
#if ( defined DA_OBS_SMOS_NRT_NN )
   call registerdaobsclass(trim(LIS_SMOSNRTNNL2smobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMOSNRTNNL2smobsId)//char(0),&
        SMOSNRTNNL2sm_setup)
   call registerreaddaobs(trim(LIS_SMOSNRTNNL2smobsId)//char(0),&
        read_SMOSNRTNNL2sm)
   call registerwritedaobs(trim(LIS_SMOSNRTNNL2smobsId)//char(0),&
        write_SMOSNRTNNL2smobs)
#endif

!YK
#if ( defined DA_OBS_SMAP_E_OPL_SM )
   call registerdaobsclass(trim(LIS_SMAPEOPLsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMAPEOPLsmobsId)//char(0),&
        SMAPEOPLsm_setup)
   call registerreaddaobs(trim(LIS_SMAPEOPLsmobsId)//char(0),&
        read_SMAPEOPLsm)
   call registerwritedaobs(trim(LIS_SMAPEOPLsmobsId)//char(0),&
        write_SMAPEOPLsmobs)
#endif

#if ( defined DA_OBS_NASA_SMAPVOD )
   call registerdaobsclass(trim(LIS_NASASMAPvodobsId),"LSM")
   call registerdaobssetup(trim(LIS_NASASMAPvodobsId)//char(0),&
        NASASMAPvod_setup)
   call registerreaddaobs(trim(LIS_NASASMAPvodobsId)//char(0),&
        read_NASASMAPvod)
   call registerwritedaobs(trim(LIS_NASASMAPvodobsId)//char(0),&
        write_NASASMAPvodobs)
#endif

#if ( defined DA_OBS_GLASS_LAI)
   call registerdaobsclass(trim(LIS_GLASSlaiobsId),"LSM")
   call registerdaobssetup(trim(LIS_GLASSlaiobsId)//char(0),&
        GLASSlai_setup)
   call registerreaddaobs(trim(LIS_GLASSlaiobsId)//char(0),&
        read_GLASSlai)
   call registerwritedaobs(trim(LIS_GLASSlaiobsId)//char(0),&
        write_GLASSlai)
#endif

#if ( defined DA_OBS_MCD15A2H_LAI)
   call registerdaobsclass(trim(LIS_MCD15A2HlaiobsId),"LSM")
   call registerdaobssetup(trim(LIS_MCD15A2HlaiobsId)//char(0),&
        MCD15A2Hlai_setup)
   call registerreaddaobs(trim(LIS_MCD15A2HlaiobsId)//char(0),&
        read_MCD15A2Hlai)
   call registerwritedaobs(trim(LIS_MCD15A2HlaiobsId)//char(0),&
        write_MCD15A2Hlai)
#endif

!Y.Kwon
#if ( defined DA_OBS_VIIRS_GVF )
   call registerdaobsclass(trim(LIS_VIIRSgvfobsId),"LSM")
   call registerdaobssetup(trim(LIS_VIIRSgvfobsId)//char(0),&
        VIIRSgvf_setup)
   call registerreaddaobs(trim(LIS_VIIRSgvfobsId)//char(0),&
        read_VIIRSgvf)
   call registerwritedaobs(trim(LIS_VIIRSgvfobsId)//char(0),&
        write_VIIRSgvfobs)
#endif
!Y.Kwon
#if ( defined DA_OBS_CDFS_GVF )
   call registerdaobsclass(trim(LIS_CDFSgvfobsId),"LSM")
   call registerdaobssetup(trim(LIS_CDFSgvfobsId)//char(0),&
        CDFSgvf_setup)
   call registerreaddaobs(trim(LIS_CDFSgvfobsId)//char(0),&
        read_CDFSgvf)
   call registerwritedaobs(trim(LIS_CDFSgvfobsId)//char(0),&
        write_CDFSgvfobs)
#endif

#if ( defined DA_OBS_NRT_SMAPSM )
   call registerdaobsclass(trim(LIS_SMAPNRTsmobsId),"LSM")
   call registerdaobssetup(trim(LIS_SMAPNRTsmobsId)//char(0),&
        SMAPNRTsm_setup)
   call registerreaddaobs(trim(LIS_SMAPNRTsmobsId)//char(0),&
        read_SMAPNRTsm)
   call registerwritedaobs(trim(LIS_SMAPNRTsmobsId)//char(0),&
        write_SMAPNRTsmobs)
#endif


#if ( defined DA_OBS_GLASS_Albedo)
   call registerdaobsclass(trim(LIS_GLASSalbedoobsId),"LSM")
   call registerdaobssetup(trim(LIS_GLASSalbedoobsId)//char(0),&
        GLASSalbedo_setup)
   call registerreaddaobs(trim(LIS_GLASSalbedoobsId)//char(0),&
        read_GLASSalbedo)
   call registerwritedaobs(trim(LIS_GLASSalbedoobsId)//char(0),&
        write_GLASSalbedo)
#endif

#if ( defined DA_OBS_MODISSPORT_LAI )
   call registerdaobsclass(trim(LIS_MODISsportLAIobsId),"LSM")
   call registerdaobssetup(trim(LIS_MODISsportLAIobsId)//char(0),&
        MODISsportLAI_setup)
   call registerreaddaobs(trim(LIS_MODISsportLAIobsId)//char(0),&
        read_MODISsportLAI)
   call registerwritedaobs(trim(LIS_MODISsportLAIobsId)//char(0),&
        write_MODISsportLAI)
#endif

#if ( defined DA_OBS_ASO_SWE )
   call registerdaobsclass(trim(LIS_ASOsweobsId),"LSM")
   call registerdaobssetup(trim(LIS_ASOsweobsId)//char(0),&
        ASO_SWE_setup)
   call registerreaddaobs(trim(LIS_ASOsweobsId)//char(0),&
        read_ASO_SWE)
   call registerwritedaobs(trim(LIS_ASOsweobsId)//char(0),&
        write_ASO_SWEobs)
#endif

#if ( defined DA_OBS_THYSM )
   call registerdaobsclass(trim(LIS_THySMid),"LSM")
   call registerdaobssetup(trim(LIS_THySMid)//char(0),&
        THySM_setup)
   call registerreaddaobs(trim(LIS_THySMid)//char(0),&
        read_THySM)
   call registerwritedaobs(trim(LIS_THySMid)//char(0),&
        write_THySM)
#endif

#if ( defined DA_OBS_SNODAS )
!MLW: SNODAS snow depth
   call registerdaobsclass(trim(LIS_SNODASobsId),"LSM")
   call registerdaobssetup(trim(LIS_SNODASobsId)//char(0), &
        SNODAS_setup)
   call registerreaddaobs(trim(LIS_SNODASobsId)//char(0),  &
        read_SNODAS)
   call registerwritedaobs(trim(LIS_SNODASobsId)//char(0), &
        write_SNODAS)
#endif

#if ( defined DA_OBS_HYDROWEBWL )
!synthetic noah soil moisture    
   call registerdaobsclass(trim(LIS_hydrowebwlId),"Routing")
   call registerdaobssetup(trim(LIS_hydrowebwlId)//char(0),hydrowebwlobs_setup)
   call registerreaddaobs(trim(LIS_hydrowebwlId)//char(0),read_hydrowebwlobs)
   call registerwritedaobs(trim(LIS_hydrowebwlId)//char(0),write_hydrowebwlobs)
#endif

#if ( defined DA_OBS_WUSUCLA )
!MLW: SNODAS snow depth
   call registerdaobsclass(trim(LIS_wusUCLAobsId),"LSM")
   call registerdaobssetup(trim(LIS_wusUCLAobsId)//char(0), &
        WUS_UCLAsnow_setup)
   call registerreaddaobs(trim(LIS_wusUCLAobsId)//char(0),  &
        read_WUS_UCLAsnow)
   call registerwritedaobs(trim(LIS_wusUCLAobsId)//char(0), &
        write_WUS_UCLAsnow)
#endif   
#endif
 end subroutine LIS_DAobs_plugin

 subroutine registerdaobsclass(funcName, funcClass)
   
   implicit none
   
   character(len=*) :: funcName
   character(len=*) :: funcClass
   
   type(DAobsFuncEntry), pointer :: dataEntry,current
   
   allocate(dataEntry)
   
   dataEntry%name = funcName
   dataEntry%class = funcClass
   
   if(.not.associated(LIS_DAobsFuncEntry%head_daobsfunc_list)) then 
      LIS_DAobsFuncEntry%head_daobsfunc_list => dataEntry
   else
      current => LIS_DAobsFuncEntry%head_daobsfunc_list
      do while(associated(current%next))
         current => current%next
      enddo
      current%next =>dataEntry
   endif
   dataEntry%next => null()
   
 end subroutine registerdaobsclass
 
!BOP
!
! !ROUTINE: LIS_isDAinstanceValid
! \label{LIS_isDAinstanceValid}
!
! !INTERFACE: 
 subroutine LIS_isDAinstanceValid(funcName, funcClass,flag)
!
! !DESCRIPTION: 
!  
!  This subroutine checks if the DA class instance is among
!  the supported list of DA class types. 
! 
!EOP
   implicit none
   
   character(len=*) :: funcName
   character(len=*) :: funcClass
   logical          :: flag
   
   type(DAobsFuncEntry), pointer :: dataEntry
   
   flag = .false. 
   dataEntry => LIS_DAobsFuncEntry%head_daobsfunc_list
   
   do while(associated(dataEntry))
      
      if(dataEntry%name.eq.funcName.and.&
           dataEntry%class.eq.funcClass) then 
         flag = .true. 
         exit
      endif
      dataEntry => dataEntry%next
   enddo
   
 end subroutine LIS_isDAinstanceValid
 
end module LIS_DAobs_pluginMod
