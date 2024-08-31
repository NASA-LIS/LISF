!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_domain_pluginMod
!BOP
!
! !MODULE: LDT_domain_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various LDT-domains. 
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
  PUBLIC :: LDT_domain_plugin  
contains
!BOP
! !ROUTINE: LDT_domain_plugin
!  \label{LDT_domain_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LDT-domain. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LDT-domain. 
!
!
! !INTERFACE:
  subroutine LDT_domain_plugin
    use LDT_pluginIndices
!EOP
    external readinput_latlon
    external readinput_hrap
    external readinput_gaussian
    external readinput_lambert
    external readinput_polar
    external readinput_merc
    external readinput_easev2

    call registerinput(trim(LDT_latlonId)//char(0),readinput_latlon)   ! regular lat/lon
    call registerinput(trim(LDT_hrapId)//char(0),readinput_hrap)
    call registerinput(trim(LDT_gaussId)//char(0),readinput_gaussian)
    call registerinput(trim(LDT_lambertId)//char(0),readinput_lambert)
    call registerinput(trim(LDT_mercId)//char(0),readinput_merc)
    call registerinput(trim(LDT_polarId)//char(0),readinput_polar)
    call registerinput(trim(LDT_easev2Id)//char(0),readinput_easev2)

  end subroutine LDT_domain_plugin
end module LDT_domain_pluginMod
