!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_OSSEmaskData_pluginMod
!BOP
!
! !MODULE: LDT_OSSEmaskData_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   external masks to be used in the OSSEs 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  26 Sep 2019;   Sujay Kumar  Initial Specification
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LDT_OSSEmaskData_plugin  

contains
!BOP
! !ROUTINE: LDT_OSSEmaskData_plugin
!  \label{LDT_OSSEmaskData_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new data mask source for OSSEs
!  The interface mandates that the following interfaces be implemented
!  and registered for each source. 
!
!
! !INTERFACE:
  subroutine LDT_OSSEmaskData_plugin

    use LDT_pluginIndices
!EOP
    use LISoutOSSEmask_Mod,         only : LISoutOSSEmask_init

    external readLISoutOSSEmask

    call registerossemasksourcesetup(trim(LDT_LISoutOSSEmaskDataId)//char(0), &
         LISoutOSSEmask_init)
    call registerreadOssemaskSource(trim(LDT_LISoutOSSEmaskDataId)//char(0), &
         readLISoutOSSEmask)

  end subroutine LDT_OSSEmaskData_plugin

end module LDT_OSSEmaskData_pluginMod
