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
! !MODULE: LVT_runmode_pluginMod
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
!   defining routines that initialize various LVT-runmodes. 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  17 Feb 2004;   Sujay Kumar  Initial Specification
! 
!EOP
module LVT_runmode_pluginMod

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LVT_runmode_plugin  
contains
!BOP
! !ROUTINE: LVT_runmode_plugin
!  \label{LVT_runmode_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LVT-runmode. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LVT-runmode. 
!
!
! !INTERFACE:
  subroutine LVT_runmode_plugin
    use LVT_pluginIndices

    implicit none
!EOP

    external LVT_init_DataComp
    external LVT_run_DataComp

    external LVT_init_DAstats
    external LVT_run_DAstats

!    external LVT_init_DAobs
!    external LVT_run_DAobs

    external LVT_init_optUE
    external LVT_run_optUE

    external LVT_init_Benchmarking
    external LVT_run_Benchmarking

    external LVT_init_557post
    external LVT_run_557post

    external LVT_init_USAFSIpost
    external LVT_run_USAFSIpost

    external LVT_init_LISpost
    external LVT_run_LISpost

    call registerlvtinit(trim(LVT_DataCompId)//char(0),LVT_init_DataComp) 
    call registerlvtrun(trim(LVT_DataCompId)//char(0),LVT_run_DataComp)

    call registerlvtinit(trim(LVT_dastatId)//char(0), LVT_init_DAstats)
    call registerlvtrun(trim(LVT_dastatId)//char(0), LVT_run_DAstats)

!    call registerlvtinit(trim(LVT_daobsId)//char(0), LVT_init_DAobs)
!    call registerlvtrun(trim(LVT_daobsId)//char(0), LVT_run_DAobs)

    call registerlvtinit(trim(LVT_optUEId)//char(0),LVT_init_optUE) 
    call registerlvtrun(trim(LVT_optUEId)//char(0),LVT_run_optUE)

    call registerlvtinit(trim(LVT_benchMarkId)//char(0),LVT_init_Benchmarking) 
    call registerlvtrun(trim(LVT_benchMarkId)//char(0),LVT_run_Benchmarking)

    call registerlvtinit(trim(LVT_557postId)//char(0),LVT_init_557post)
    call registerlvtrun(trim(LVT_557postId)//char(0),LVT_run_557post)

    call registerlvtinit(trim(LVT_USAFSIpostId)//char(0),LVT_init_USAFSIpost)
    call registerlvtrun(trim(LVT_USAFSIpostId)//char(0),LVT_run_USAFSIpost)

    call registerlvtinit(trim(LVT_LISpostId)//char(0),LVT_init_LISpost)
    call registerlvtrun(trim(LVT_LISpostId)//char(0),LVT_run_LISpost)

  end subroutine LVT_runmode_plugin
end module LVT_runmode_pluginMod
