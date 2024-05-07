!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: LIS_initialize_registries
!  \label{LIS_initialize_registries}
!
! !REVISION HISTORY:
!  19 Aug 2010: Sujay Kumar; Initial specification
!  17 Jan 2011: David Mocko, added max/min greenness & slope type
!
! !INTERFACE:
subroutine LIS_initialize_registries()
! !USES:
  use ESMF 
  use LIS_param_pluginMod,           only : LIS_laisai_plugin,&
       LIS_alb_plugin, LIS_gfrac_plugin , LIS_roughness_plugin,&
       LIS_emissivity_plugin
  use LIS_runmode_pluginMod,         only : LIS_runmode_plugin
  use LIS_dataassim_pluginMod,       only : LIS_dataassim_plugin
  use LIS_biasEstimation_pluginMod,  only : LIS_biasEstimation_plugin
  use LIS_optUEAlgorithm_pluginMod,  only : LIS_optUEAlgorithm_plugin
!  use LIS_optUEType_pluginMod,      only : LIS_optUEType_plugin
  use LIS_lsmoptue_pluginMod,        only : LIS_lsmoptue_plugin
  use LIS_PEobs_pluginMod,           only : LIS_PEobs_plugin
  use LIS_ObjFunc_pluginMod,         only : LIS_ObjFunc_plugin
  use LIS_RTM_pluginMod,             only : LIS_RTM_plugin
  use LIS_rtmoptue_pluginMod,        only : LIS_rtmoptue_plugin
  use LIS_lsmrtm_pluginMod,          only : LIS_lsmrtm_plugin
  use LIS_landslidemodel_pluginMod,  only : LIS_landslidemodel_plugin
  use LIS_routing_pluginMod,         only : LIS_routing_plugin
  use LIS_lsmrouting_pluginMod,      only : LIS_lsmrouting_plugin
  use LIS_irrigationmodel_pluginMod, only : LIS_irrigationmodel_plugin    
  use LIS_lsmirrigation_pluginMod,   only : LIS_lsmirrigation_plugin 
  use LIS_runoffdata_pluginMod,      only : LIS_runoffdata_plugin
  use LIS_DAobs_pluginMod,           only : LIS_DAobs_plugin
  use LIS_lakemodel_pluginMod,       only : LIS_lakemodel_plugin
  use LIS_forecastAlg_pluginMod
  use LIS_lsm_pluginMod
  use LIS_sublsm_pluginMod
  use LIS_glaciermodel_pluginMod
  use LIS_glacierrouting_pluginMod
  use LIS_lsmcpl_pluginMod
  use LIS_lsmda_pluginMod
  use LIS_routingda_pluginMod
  use LIS_perturb_pluginMod

!
! !DESCRIPTION:
!
!  The code in this file initializes registries that set up the component plugin
!  definitions.
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_runmode\_plugin](\ref{LIS_runmode_plugin}) \newline
!    sets up function table registries for implemented runmodes
!   \item[LIS\_laisai\_plugin](\ref{LIS_laisai_plugin}) \newline
!    sets up function table registries for implemented LAI/SAI data sources
!   \item[LIS\_alb\_plugin](\ref{LIS_alb_plugin}) \newline
!    sets up function table registries for implemented
!    albedo data sources
!   \item[LIS\_gfrac\_plugin](\ref{LIS_gfrac_plugin}) \newline
!    sets up function table registries for implemented
!    greenness fraction data sources
!   \item[LIS\_optUEAlgorithm\_plugin](\ref{LIS_optUEAlgorithm_plugin}) \newline
!    sets up function table registries for implemented 
!    parameter estimation algorithms
!   %\item[LIS\_optUEType\_plugin](\ref{LIS_optUEType_plugin}) \newline
!   % sets up function table registries for implemented 
!   % parameter estimation types
!   \item[LIS\_lsmoptue\_plugin](\ref{LIS_lsmoptue_plugin}) \newline
!    sets up function table registries for implemented 
!    lsm plugins related to parameter estimation
!   \item[LIS\_PEobs\_plugin](\ref{LIS_PEobs_plugin}) \newline
!    sets up function table registries for implemented 
!    plugins for objective space handling related to parameter estimation
!  \item[LIS\_DAobs\_plugin](\ref{LIS_DAobs_plugin}) \newline
!    sets up function table registries for implemented 
!    observation sources to be assimilated.
!  \item[LIS\_lakemodel\_plugin] (\ref{LIS_lakemodel_plugin}) \newline
!    sets up function table registries for implemented land surface models
!  \item[LIS\_lsm\_plugin] (\ref{LIS_lsm_plugin}) \newline
!    sets up function table registries for implemented land surface models
!  \item[LIS\_lsmcpl\_plugin] (\ref{LIS_lsmcpl_plugin}) \newline
!    sets up function table registries for implemented land surface models
!    used in coupled modes
!  \item[LIS\_lsmda\_plugin] (\ref{LIS_lsmda_plugin}) \newline
!    sets up function table registries for land surface models used
!    in data assimilation
!  \item[LIS\_perturb\_plugin](\ref{LIS_perturb_plugin}) \newline
!    sets up the function table registries for implemented
!    perturbation algorithms
!  \end{description}
!EOP
  implicit none

  call LIS_RTM_plugin
  call LIS_lsm_plugin
  call LIS_sublsm_plugin
  call LIS_lsmcpl_plugin
  call LIS_lsmda_plugin
  call LIS_routingda_plugin
  call LIS_lsmrtm_plugin
  call LIS_landslidemodel_plugin
  call LIS_routing_plugin
  call LIS_lsmrouting_plugin

  call LIS_runmode_plugin
  call LIS_laisai_plugin
  call LIS_alb_plugin
  call LIS_gfrac_plugin
  call LIS_roughness_plugin
  call LIS_emissivity_plugin

  call LIS_dataassim_plugin    
  call LIS_biasestimation_plugin  

  call LIS_optUEAlgorithm_plugin
!  call LIS_optUEType_plugin
  call LIS_lsmoptue_plugin
  call LIS_rtmoptue_plugin
  call LIS_PEobs_plugin
  call LIS_ObjFunc_plugin

  call LIS_DAobs_plugin       

  call LIS_irrigationmodel_plugin    
  call LIS_lsmirrigation_plugin 

  call LIS_runoffdata_plugin
  call LIS_lakemodel_plugin
  call LIS_glaciermodel_plugin
  call LIS_glacierrouting_plugin

  call LIS_perturb_plugin
  call LIS_forecastAlg_plugin

end subroutine LIS_initialize_registries

