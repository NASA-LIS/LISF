!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module LDT_NatureRunData_pluginMod
!BOP
!
! !MODULE: LDT_NatureRunData_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining sources of nature runs to be used in OSSEs. 
! 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  08 Aug 2019;   Sujay Kumar  Initial Specification
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LDT_NatureRunData_plugin  

contains
!BOP
! !ROUTINE: LDT_NatureRunData_plugin
!  \label{LDT_NatureRunData_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new nature run source
!  The interface mandates that the following interfaces be implemented
!  and registered for each data source
!
!
! !INTERFACE:
  subroutine LDT_NatureRunData_plugin

    use LDT_pluginIndices
!EOP
    use LISoutNatureRun_Mod,         only : LISoutNatureRun_init

    external readLISoutNatureRun

    call registernaturerunsourcesetup(trim(LDT_LISoutNatureRunDataId)//char(0), &
         LISoutNatureRun_init)
    call registerreadNatureRunSource(trim(LDT_LISoutNatureRunDataId)//char(0), &
         readLISoutNatureRun)

  end subroutine LDT_NatureRunData_plugin

end module LDT_NatureRunData_pluginMod
