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
! !MODULE: LVT_domain_pluginMod
!  \label(LVT_domain_pluginMod)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various LVT-domains. 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  17 Feb 2004;   Sujay Kumar  Initial Specification
! 
!EOP
module LVT_domain_pluginMod

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LVT_domain_plugin  
contains
!BOP
! 
! !ROUTINE: LVT_domain_plugin
!  \label{LVT_domain_plugin}
!
! !INTERFACE:
  subroutine LVT_domain_plugin
! 
! !USES:   
    use LVT_pluginIndices
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LVT-domain. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LVT-domain. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    external readinput_latlon
    external readinput_lambert
    external readinput_UTM

    call registerinput(trim(LVT_latlonId)//char(0),readinput_latlon) !regular lat/lon
    call registerinput(trim(LVT_lambertId)//char(0),readinput_lambert)
    call registerinput(trim(LVT_utmId)//char(0), readinput_UTM)

  end subroutine LVT_domain_plugin
end module LVT_domain_pluginMod
