!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_NatureRunData_pluginMod
!BOP
!
! !MODULE: LDT_NatureRunData_pluginMod
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
  PUBLIC :: LDT_NatureRunData_plugin  

contains
!BOP
! !ROUTINE: LDT_NatureRunData_plugin
!  \label{LDT_NatureRunData_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LDT-obs. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LDT-obs. 
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
