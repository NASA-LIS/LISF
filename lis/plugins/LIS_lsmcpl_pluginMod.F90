!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_lsmcpl_pluginMod
!BOP
!
! !MODULE: LIS_lsmcpl_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions that
!   need to be defined to enable the use of a land surface
!   model in coupled land-atmosphere setups.
!
! !REVISION HISTORY:
!  09 Oct 07    Sujay Kumar  Initial Specification
!  14 Dec 15    Eric Kemp    Added Noah 2.7.1, Noah 3.2, Noah 3.6 WRF coupling.
!  22 Feb 19    Chandana Gangodagamage  Added NoahMP 3.6 for WRFHydro coupling 
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_lsmcpl_plugin
contains
!BOP
! !ROUTINE: LIS_lsmcpl_plugin
!  \label{LIS_lsmcpl_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a LSM
! in a coupled land-atmosphere setup. The interface mandates that the
! following routines be implemented and registered for each of the LSM
! that is used.
!
!  If used for a coupled run with an atmospheric component:
!  \begin{description}
!  \item[Specify an export state]
!      Routines to specify the variables to be sent to the atmos. component
!      (to be registered using {\tt registerlsmsetexport} and later called
!       using {\tt lsmsetexport})
!  \end{description}
!
!  The user-defined functions are included in the registry using two
!  indices. For example, consider the Noah LSM used in a mode coupled
!  to WRF, the setexport method is  incorporated in the registry with
!  an index of 1 corresponding to Noah and 3 corresponding to the runmode
!  coupled to WRF as follows:
!
!  \begin{verbatim}
!    call registerlsmsetexport(1,3,noah_setwrfexport)
!  \end{verbatim}
!
!   The function registered above is invoked using a generic call as
!   follows:
!
!  \begin{verbatim}
!    call lsmsetexport(1,3)      - calls noah_setwrfexport
!  \end{verbatim}
!
!   In the LIS code, the above call is typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call lsmsetexport(lis%lsm, lis%runmode)
!   \end{verbatim}
!   where $lis\%lsm$ and %lis\%runmode% are  set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_lsmcpl_plugin
!EOP

#if ( defined RM_WRF_COUPLING )
   use LIS_pluginIndices

#if ( defined SM_NOAH_2_7_1 )
   external noah271_wrf_f2t
   external noah271_setwrfexport
#endif

#if ( defined SM_NOAH_3_2 )
   external noah32_wrf_f2t
   external noah32_setwrfexport
#endif

#if ( defined SM_NOAH_3_3 )
   external noah33_wrf_f2t
   external noah33_setwrfexport
#endif

#if ( defined SM_NOAH_3_6 )
   external noah36_wrf_f2t
   external noah36_setwrfexport
#endif

#if ( defined SM_NOAHMP_3_6 )
   external noahMP36_setwrfexport
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   external noahMP401_setwrfexport
#endif

#if 0
!tight coupling interfaces: no ESMF
   external noah271_wrf_f2t
   external noah271_setwrfexport

   external clm2_wrf_f2t
   external clm2_setwrfexport

   external cable_wrf_f2t
   external cable_setwrfexport

   external tess_wrf_f2t
   external tess_setwrfexport

!loose coupling : using ESMF
   external noah271_wrf_esmff2t
   external noah271_setwrfesmfexport

   external noah271_gce_f2t
   external noah271_setgceexport

   external noah31_wrf_f2t
   external noah31_setwrfexport

   external noah31_gce_f2t
   external noah31_setgceexport

   external noah32_wrf_f2t
   external noah32_setwrfexport

!   external clm2_wrf_esmff2t
!   external clm2_setwrfesmfexport


! for coupled runs
!    call registerlsmf2t(trim(LIS_noah271Id)//"+"//&
!         trim(LIS_wrfcpl1Id)//char(0),noah271_wrf_esmff2t)
!    call registerlsmcplsetexport(trim(LIS_noah271Id)//"+"//&
!         trim(LIS_wrfcpl1Id)//char(0),noah271_setwrfesmfexport)

   call registerlsmf2t(trim(LIS_noah271Id)//"+"//&
        trim(LIS_wrfcplId)//char(0),noah271_wrf_f2t)
   call registerlsmcplsetexport(trim(LIS_noah271Id)//"+"//&
        trim(LIS_wrfcplId)//char(0),noah271_setwrfexport)

   call registerlsmf2t(trim(LIS_clm2id)//"+"//&
        trim(LIS_wrfcplId)//char(0),clm2_wrf_f2t)
   call registerlsmcplsetexport(trim(LIS_clm2id)//"+"//&
        trim(LIS_wrfcplId)//char(0),clm2_setwrfexport)

   call registerlsmf2t(trim(LIS_cableId)//"+"//&
        trim(LIS_wrfcpl2Id)//char(0),cable_wrf_f2t)
   call registerlsmcplsetexport(trim(LIS_cableId)//"+"//&
        trim(LIS_wrfcpl2Id)//char(0),cable_setwrfexport)

   call registerlsmf2t(trim(LIS_tessId)//"+"//&
        trim(LIS_wrfcplId)//char(0),tess_wrf_f2t)
   call registerlsmcplsetexport(trim(LIS_tessId)//"+"//&
        trim(LIS_wrfcplId)//char(0),tess_setwrfexport)

   call registerlsmf2t(trim(LIS_noah271Id)//"+"//&
        trim(LIS_gcecplId)//char(0),noah271_gce_f2t)
   call registerlsmcplsetexport(trim(LIS_noah271Id)//"+"//&
        trim(LIS_gcecplId)//char(0),noah271_setgceexport)

   call registerlsmf2t(trim(LIS_noah32Id)//"+"//&
        trim(LIS_wrfcplId)//char(0),noah32_wrf_f2t)
   call registerlsmcplsetexport(trim(LIS_noah32Id)//"+"//&
        trim(LIS_wrfcplId)//char(0),noah32_setwrfexport)

   call registerlsmf2t(trim(LIS_noah31Id)//"+"//&
        trim(LIS_gcecplId)//char(0),noah31_gce_f2t)
   call registerlsmcplsetexport(trim(LIS_noah31Id)//"+"//&
        trim(LIS_gcecplId)//char(0),noah31_setgceexport)
#endif

#if ( defined SM_NOAH_2_7_1 )
   call registerlsmf2t(trim(LIS_noah271Id)//"+"//&
                       trim(LIS_wrfcplId)//char(0), &
                       noah271_wrf_f2t)
   call registerlsmcplsetexport(trim(LIS_noah271Id)//"+"//&
                                trim(LIS_wrfcplId)//char(0), &
                                noah271_setwrfexport)
#endif

#if ( defined SM_NOAH_3_2 )
   call registerlsmf2t(trim(LIS_noah32Id)//"+"//&
                       trim(LIS_wrfcplId)//char(0), &
                       noah32_wrf_f2t)
   call registerlsmcplsetexport(trim(LIS_noah32Id)//"+"//&
                                trim(LIS_wrfcplId)//char(0), &
                                noah32_setwrfexport)
#endif

#if ( defined SM_NOAH_3_3 )
    call registerlsmf2t(trim(LIS_noah33Id)//"+"//&
                        trim(LIS_wrfcplId)//char(0), &
                        noah33_wrf_f2t)
    call registerlsmcplsetexport(trim(LIS_noah33Id)//"+"//&
                                 trim(LIS_wrfcplId)//char(0), &
                                 noah33_setwrfexport)
    call registerlsmcplsetexport(trim(LIS_noah33Id)//"+"//&
                                 trim(LIS_nuopccplId)//char(0), &
                                 noah33_setwrfexport)
    call registerlsmcplsetexport(trim(LIS_noah33Id)//"+"//&
                                 trim(LIS_retroId)//char(0), &
                                 noah33_setwrfexport)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmcplsetexport(trim(LIS_noahmp36Id)//"+"//&
                                 trim(LIS_wrfcplId)//char(0), &
                                 noahMP36_setwrfexport)
    call registerlsmcplsetexport(trim(LIS_noahmp36Id)//"+"//&
                                 trim(LIS_nuopccplId)//char(0), &
                                 noahMP36_setwrfexport)
    call registerlsmcplsetexport(trim(LIS_noahmp36Id)//"+"//&
                                 trim(LIS_retroId)//char(0), &
                                 noahMP36_setwrfexport)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmcplsetexport(trim(LIS_noahmp401Id)//"+"//&
                                 trim(LIS_wrfcplId)//char(0), &
                                 noahMP401_setwrfexport)
    call registerlsmcplsetexport(trim(LIS_noahmp401Id)//"+"//&
                                 trim(LIS_nuopccplId)//char(0), &
                                 noahMP401_setwrfexport)
    call registerlsmcplsetexport(trim(LIS_noahmp401Id)//"+"//&
                                 trim(LIS_retroId)//char(0), &
                                 noahMP401_setwrfexport)
#endif

#if ( defined SM_NOAH_3_6 )
    call registerlsmf2t(trim(LIS_noah36Id)//"+"//&
                        trim(LIS_wrfcplId)//char(0), &
                        noah36_wrf_f2t)
    call registerlsmcplsetexport(trim(LIS_noah36Id)//"+"//&
                                 trim(LIS_wrfcplId)//char(0), &
                                 noah36_setwrfexport)
#endif
#endif

end subroutine LIS_lsmcpl_plugin
end module LIS_lsmcpl_pluginMod
