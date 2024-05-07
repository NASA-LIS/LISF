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
module LIS_param_pluginMod
!BOP
!
! !MODULE: LIS_param_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   defining routines to read various sources of parameters maps.
!   The user defined functions are incorporated into
!   the appropriate registry to be later invoked through generic calls.
!
!
! !REVISION HISTORY:
!  11 Dec 03    Sujay Kumar  Initial Specification
!  17 Jan 2011: David Mocko, added max/min greenness & slope type
!   4 Nov 2014: Jonathan Case, added support for daily NESDIS/VIIRS GVF for Noah
!
!EOP
  use LIS_pluginIndices
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_laisai_plugin
  PUBLIC :: LIS_gfrac_plugin
  PUBLIC :: LIS_alb_plugin
  PUBLIC :: LIS_roughness_plugin
  PUBLIC :: LIS_emissivity_plugin

contains

!BOP
! !ROUTINE: LIS_laisai_plugin
!  \label{LIS_laisai_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new LAI/SAI datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the LAI data]
!      Routines to retrieve the LAI data
!      (to be registered using {\tt registerreadlai} and later called
!       using the generic {\tt readlai} method)
!  \item[read the SAI data]
!      Routines to retrieve the SAI data
!      (to be registered using {\tt registerreadsai} and later called
!       using the generic {\tt readsai} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the LAI/SAI
!  datasets from the AVHRR source with LIS using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LIS uses and index of 1
!
!  \begin{verbatim}
!    call registerreadlai(1,1,read_ll_avhrrlai)
!    call registerreadsai(1,1,read_ll_avhrrsai)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readlai(1,1)  - calls read_ll_avhrrlai
!    call readsai(1,1)  - calls read_ll_avhrrsai
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call readlai(lis%domain,lislaisrc)
!    call readsai(lis%domain,lislaisrc)
!   \end{verbatim}
!   where $lis\%domain$ and $lis\%laisrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_laisai_plugin
!EOP

!  external setup_AVHRRlai, read_AVHRRlai,read_AVHRRsai
!  external setup_MODISlai, read_MODISlai,read_MODISsai
!  external setup_tiled_AVHRRlai, read_tiled_AVHRRlai, read_tiled_AVHRRsai
!  external setup_MMFlai, read_MMFlai, read_MMFsai

#if ( defined PARAM_MODIS_REAL_TIME )
  ! real-time MODIS data
  external setup_MODIS_RT_lai, read_MODIS_RT_lai,read_MODIS_RT_sai
#endif

#if ( defined PARAM_ALMIPII_LAI )
   external setup_ALMIPIIlai
   external read_ALMIPIIlai
#endif

!  call registerlaisetup(LIS_avhrrlaiId, setup_AVHRRlai)
!  call registerreadlai(LIS_avhrrlaiId,read_AVHRRlai)
!  call registerreadsai(LIS_avhrrlaiId,read_AVHRRsai)
!
!  call registerlaisetup(LIS_modislaiId, setup_MODISlai)
!  call registerreadlai(LIS_modislaiId,read_MODISlai)
!
!  call registerreadsai(LIS_modislaiId,read_MODISsai)
!  call registerlaisetup(LIS_tiledavhrrlaiId, setup_tiled_AVHRRlai)
!  call registerreadlai(LIS_tiledavhrrlaiId,read_tiled_AVHRRlai)
!  call registerreadsai(LIS_tiledavhrrlaiId,read_tiled_AVHRRsai)
!
!  call registerlaisetup(LIS_mmfdataId, setup_MMFlai)
!  call registerreadlai(LIS_mmfdataId,read_MMFlai)
!  call registerreadsai(LIS_mmfdataId,read_MMFsai)

#if ( defined PARAM_ALMIPII_LAI )
   call registerlaisetup(trim(LIS_ALMIPIIlaiId)//char(0),setup_ALMIPIIlai)
   call registerreadlai(trim(LIS_ALMIPIIlaiId)//char(0),read_ALMIPIIlai)
#endif

#if ( defined PARAM_MODIS_REAL_TIME )
   call registerlaisetup(trim(LIS_modis_RT_laiId)//char(0),setup_MODIS_RT_lai)
   call registerreadlai(trim(LIS_modis_RT_laiId)//char(0),read_MODIS_RT_lai)
   call registerreadsai(trim(LIS_modis_RT_laiId)//char(0),read_MODIS_RT_sai)
#endif

end subroutine LIS_laisai_plugin


!BOP
! !ROUTINE: LIS_gfrac_plugin
!  \label{LIS_gfrac_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new greenness datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the greenness data]
!      Routines to retrieve the greenness data
!      (to be registered using {\tt registerreadgfrac} and later called
!       using the generic {\tt readgfrac} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the greenness
!  datasets from the NCEP source with LIS using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LIS uses and index of 1
!
!  \begin{verbatim}
!    call registerreadgfrac(1,1,read_llgfrac)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readgfrac(1,1)  - calls read_llgfrac
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!     call readgfrac(lis%domain,lis%gfracsrc)
!   \end{verbatim}
!   where $lis\%domain$ and $lis\%gfracsrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_gfrac_plugin
!EOP

#if ( defined PARAM_NESDIS_WEEKLY )
   external setup_NESDISgfrac
   external read_NESDISgfrac
#endif

#if ( defined PARAM_SPORT )
   external setup_SPORTgfrac
   external read_SPORTgfrac
#endif

#if ( defined PARAM_VIIRS )
   external setup_VIIRSgfrac
   external read_VIIRSgfrac
#endif

#if ( defined PARAM_ALMIPII_GFRAC )
   external setup_ALMIPIIgfrac
   external read_ALMIPIIgfrac
#endif

#if ( defined PARAM_NESDIS_WEEKLY )
   call registergfracsetup(trim(LIS_NESDISgfracId)//char(0),setup_NESDISgfrac)
   call registerreadgfrac(trim(LIS_NESDISgfracId)//char(0),read_NESDISgfrac)
#endif

#if ( defined PARAM_SPORT )
   call registergfracsetup(trim(LIS_SPORTgfracId)//char(0),setup_SPORTgfrac)
   call registerreadgfrac(trim(LIS_SPORTgfracId)//char(0),read_SPORTgfrac)
#endif

#if ( defined PARAM_VIIRS )
   call registergfracsetup(trim(LIS_VIIRSgfracId)//char(0),setup_VIIRSgfrac)
   call registerreadgfrac(trim(LIS_VIIRSgfracId)//char(0),read_VIIRSgfrac)
#endif

#if ( defined PARAM_ALMIPII_GFRAC )
   call registergfracsetup(trim(LIS_ALMIPIIgfracId)//char(0),setup_ALMIPIIgfrac)
   call registerreadgfrac(trim(LIS_ALMIPIIgfracId)//char(0),read_ALMIPIIgfrac)
#endif
end subroutine LIS_gfrac_plugin


!BOP
! !ROUTINE: LIS_roughness_plugin
!  \label{LIS_roughness_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new greenness datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the greenness data]
!      Routines to retrieve the greenness data
!      (to be registered using {\tt registerreadroughness} and later called
!       using the generic {\tt readroughness} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the greenness
!  datasets from the NCEP source with LIS using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LIS uses and index of 1
!
!  \begin{verbatim}
!    call registerreadroughness(1,1,read_llroughness)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readroughness(1,1)  - calls read_llroughness
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!     call readroughness(lis%domain,lis%roughnesssrc)
!   \end{verbatim}
!   where $lis\%domain$ and $lis\%roughnesssrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_roughness_plugin
!EOP

#if ( defined PARAM_ALMIPII_ROUGHNESS )
   external setup_ALMIPIIroughness
   external read_ALMIPIIroughness
#endif

#if ( defined PARAM_ALMIPII_ROUGHNESS )
   call registerroughnesssetup(trim(LIS_ALMIPIIroughnessId)//char(0),&
                               setup_ALMIPIIroughness)
   call registerreadroughness(trim(LIS_ALMIPIIroughnessId)//char(0),&
                              read_ALMIPIIroughness)
#endif
end subroutine LIS_roughness_plugin

!BOP
! !ROUTINE: LIS_emissivity_plugin
!  \label{LIS_emissivity_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new greenness datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the greenness data]
!      Routines to retrieve the greenness data
!      (to be registered using {\tt registerreademissivity} and later called
!       using the generic {\tt reademissivity} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the greenness
!  datasets from the NCEP source with LIS using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LIS uses and index of 1
!
!  \begin{verbatim}
!    call registerreademissivity(1,1,read_llemissivity)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call reademissivity(1,1)  - calls read_llemissivity
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!     call reademissivity(lis%domain,lis%emissivitysrc)
!   \end{verbatim}
!   where $lis\%domain$ and $lis\%emissivitysrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_emissivity_plugin
!EOP

#if ( defined PARAM_ALMIPII_EMISSIVITY )
   external setup_ALMIPIIemiss
   external read_ALMIPIIemiss
#endif

#if ( defined PARAM_ALMIPII_EMISSIVITY )
   call registeremissivitysetup(trim(LIS_ALMIPIIemissId)//char(0), &
                                setup_ALMIPIIemiss)
   call registerreademissivity(trim(LIS_ALMIPIIemissId)//char(0), &
                               read_ALMIPIIemiss)
#endif
end subroutine LIS_emissivity_plugin

!BOP
! !ROUTINE: LIS_alb_plugin
!  \label{LIS_alb_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new albedo datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the max snow albedo data]
!      Routines to retrieve the max snow albedo data
!      (to be registered using {\tt registerreadmxsnoalb} and later called
!       using the generic {\tt readmxsnoalb} method)
!  \item[read the albedo climatology data]
!      Routines to retrieve the albedo climatology data
!      (to be registered using {\tt registerreadalbedo} and later called
!       using the generic {\tt readalbedo} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the albedo
!  datasets from the NCEP source with LIS using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LIS uses and index of 1
!
!  \begin{verbatim}
!    call registerreadmxsnoalb(1,1,read_llmxsnoalb)
!    call registerreadalbedo(1,1,read_llalbedo)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readmxsnoalb(1,1)  - calls read_llmxsnoalb
!    call readalbedo(1,1)    - calls read_llalbedo
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!     call readmxsnoalb(lis%domain,lis%albedosrc)
!     call readalbedo(lis%domain,lis%albedoarc)
!   \end{verbatim}
!   where $lis\%domain$ and $lis\%albedosrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_alb_plugin
!EOP

#if ( defined PARAM_ALMIPII_ALBEDO )
   external setup_ALMIPIIalbedo
   external read_ALMIPIIalbedo
#endif

#if ( defined PARAM_ALMIPII_ALBEDO )
   call registeralbedosetup(trim(LIS_ALMIPIIalbedoId)//char(0), &
                            setup_ALMIPIIalbedo)
   call registerreadalbedo(trim(LIS_ALMIPIIalbedoId)//char(0), &
                           read_ALMIPIIalbedo)
#endif

end subroutine LIS_alb_plugin

end module LIS_param_pluginMod
