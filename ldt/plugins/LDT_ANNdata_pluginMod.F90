!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_ANNdata_pluginMod
!BOP
!
! !MODULE: LDT_ANNdata_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various LDT-obss. 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  17 Feb 2004;   Sujay Kumar  Initial Specification
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LDT_ANNinputdata_plugin  
  PUBLIC :: LDT_ANNoutputdata_plugin  

contains
!BOP
! !ROUTINE: LDT_ANNinputdata_plugin
!  \label{LDT_ANNinputdata_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LDT-obs. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LDT-obs. 
!
!
! !INTERFACE:
  subroutine LDT_ANNinputdata_plugin

    use LDT_pluginIndices
!EOP
    use LISlsmSM_ANNdataMod,         only : LISlsmSM_ANNdataInit
    use MODISlst_ANNdataMod,         only : MODISlst_ANNdataInit
    use MOD10A1_ANNdataMod,          only : MOD10A1_ANNdataInit
    use GCOMWAMSR2TB_ANNdataMod,     only : GCOMWAMSR2TB_ANNdatainit

    external readLISlsmSMANNdata
    external readMODISlstANNdata
    external readMOD10A1ANNdata
    external readGCOMWAMSR2TbANNdata

    call registerANNinputsourcesetup(trim(LDT_ANNLISlsmSMId)//char(0), &
         LISlsmSM_ANNdataInit)
    call registerreadANNinputsource(trim(LDT_ANNLISlsmSMId)//char(0), &
         readLISlsmSMANNdata)

    call registerANNinputsourcesetup(trim(LDT_ANNMODISlstobsId)//char(0), &
         MODISlst_ANNdataInit)
    call registerreadANNinputsource(trim(LDT_ANNMODISlstobsId)//char(0), &
         readMODISlstANNdata)


    call registerANNinputsourcesetup(trim(LDT_ANNMOD10A1obsId)//char(0), &
         MOD10A1_ANNdataInit)
    call registerreadANNinputsource(trim(LDT_ANNMOD10A1obsId)//char(0), &
         readMOD10A1ANNdata)

    call registerANNinputsourcesetup(trim(LDT_ANNGCOMWAMSR2TbobsId)//char(0), &
         GCOMWAMSR2TB_ANNdatainit)
    call registerreadANNinputsource(trim(LDT_ANNGCOMWAMSR2TbobsId)//char(0), &
         readGCOMWAMSR2TbANNdata)

  end subroutine LDT_ANNinputdata_plugin

!BOP
! !ROUTINE: LDT_ANNoutputdata_plugin
!  \label{LDT_ANNoutputdata_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LDT-obs. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LDT-obs. 
!
!
! !INTERFACE:
  subroutine LDT_ANNoutputdata_plugin

! !USES:
    use LDT_pluginIndices
    use syntheticsm_ANNdataMod,      only : syntheticsm_ANNdatainit
    use LPRM_AMSREsm_ANNdataMod,     only : LPRM_AMSRE_ANNdataInit
    use GHCN_ANNdataMod,             only : GHCN_ANNdatainit
!EOP

    external readsyntheticsmANNdata
    external readLPRM_AMSRE_ANNdata
    external readGHCNANNdata

    call registerANNoutputsourcesetup(trim(LDT_ANNsynSMId)//char(0), &
         syntheticsm_ANNdatainit)
    call registerreadANNoutputsource(trim(LDT_ANNsynSMId)//char(0),&
         readsyntheticsmANNdata)

    call registerANNoutputsourcesetup(trim(LDT_ANNLPRMAMSREsmobsId)//char(0), &
         LPRM_AMSRE_ANNdataInit)
    call registerreadANNoutputsource(trim(LDT_ANNLPRMAMSREsmobsId)//char(0),&
         readLPRM_AMSRE_ANNdata)

    call registerANNoutputsourcesetup(trim(LDT_ANNGHCNsnwdobsId)//char(0), &
         GHCN_ANNdatainit)
    call registerreadANNoutputsource(trim(LDT_ANNGHCNsnwdobsId)//char(0), &
         readGHCNANNdata)

  end subroutine LDT_ANNoutputdata_plugin
end module LDT_ANNdata_pluginMod
