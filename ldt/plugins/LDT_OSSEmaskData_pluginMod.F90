!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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
    use AMSR2OSSEmask_Mod,          only : AMSR2OSSEmask_init
    use TSMMOSSEmask_Mod,           only : TSMMOSSEmask_init
    use MODISOSSEmask_Mod,          only : MODISOSSEmask_init
    use Sentinel1AOSSEmask_Mod,     only : Sentinel1AOSSEmask_init

    external readLISoutOSSEmask
    external readAMSR2OSSEmask
    external readTSMMOSSEmask
    external readMODISOSSEmask
    external readSentinel1AOSSEmask

    call registerossemasksourcesetup(trim(LDT_LISoutOSSEmaskDataId)//char(0), &
         LISoutOSSEmask_init)
    call registerreadOssemaskSource(trim(LDT_LISoutOSSEmaskDataId)//char(0), &
         readLISoutOSSEmask)

    call registerossemasksourcesetup(trim(LDT_AMSR2OSSEmaskDataId)//char(0), &
         AMSR2OSSEmask_init)
    call registerreadOssemaskSource(trim(LDT_AMSR2OSSEmaskDataId)//char(0), &
         readAMSR2OSSEmask)

    call registerossemasksourcesetup(trim(LDT_TSMMOSSEmaskDataId)//char(0), &
         TSMMOSSEmask_init)
    call registerreadOssemaskSource(trim(LDT_TSMMOSSEmaskDataId)//char(0), &
         readTSMMOSSEmask)

    call registerossemasksourcesetup(trim(LDT_MODISOSSEmaskDataId)//char(0), &
         MODISOSSEmask_init)
    call registerreadOssemaskSource(trim(LDT_MODISOSSEmaskDataId)//char(0), &
         readMODISOSSEmask)    

    call registerossemasksourcesetup(trim(LDT_Sentinel1AOSSEmaskDataId)//char(0), &
         Sentinel1AOSSEmask_init)
    call registerreadOssemaskSource(trim(LDT_Sentinel1AOSSEmaskDataId)//char(0), &
         readSentinel1AOSSEmask)
    
  end subroutine LDT_OSSEmaskData_plugin

end module LDT_OSSEmaskData_pluginMod
