!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_statDscaleMod
!BOP
!
! !MODULE: LDT_statDscaleMod
! 
! !DESCRIPTION:
! 
! !REVISION HISTORY: 
!  3 Feb 2016:  Sujay Kumar; initial specification
!

  use ESMF
  use LDT_coreMod
  use LDT_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LDT_statDscaleInit
  PUBLIC :: LDT_diagnoseForcStatDscale
  PUBLIC :: LDT_computeForcStatDscaleParams
  PUBLIC :: LDT_outputForcStatDscaleParams
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!----------------------------------------------------------------------------- 

!EOP

contains

!BOP
! !ROUTINE: LDT_statDscaleInit
! \label{LDT_statDscaleInit}
!
! !REVISION HISTORY:
!
! !INTERFACE:
  subroutine LDT_statDscaleInit()

! !USES:
    use LDT_statDscale_pluginMod

    implicit none
! !ARGUMENTS: 

! !DESCRIPTION:
!
!
!EOP

    integer          :: rc
! ___________________________________________

    write(LDT_logunit,*)" "
    write(LDT_logunit,*)" - - - - - - Statistical Downscaling - - - - - -"
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%statDscaleType,&
         label="Statistical downscaling method:",&
         rc=rc)
    call LDT_verify(rc,'Statistical downscaling method: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%statDscaleMode,&
         label="Statistical downscaling mode:",&
         rc=rc)
    call LDT_verify(rc,'Statistical downscaling mode: not defined')
    
    write(LDT_logunit,*) "-- Running Metforcing Downscaling Method: ",&
         trim(LDT_rc%statDscaleType)//" -- "
    
    call LDT_statdscale_plugin

    call initstatdscale(trim(LDT_rc%statDscaleType)//char(0))

  end subroutine LDT_statDscaleInit


  subroutine LDT_diagnoseForcStatDscale(n,pass)

    integer              :: n 
    integer              :: pass

    call diagnosestatdscale(trim(LDT_rc%statDscaleType)//char(0),&
         n, pass)
    
  end subroutine LDT_diagnoseForcStatDscale

  subroutine LDT_computeForcStatDscaleParams(pass)

    integer              :: pass

    call computestatdscale(trim(LDT_rc%statDscaleType)//char(0),&
         pass)
    
  end subroutine LDT_computeForcStatDscaleParams


  subroutine LDT_outputForcStatDscaleParams(pass)

    integer              :: pass

    call outputstatdscale(trim(LDT_rc%statDscaleType)//char(0),&
         pass)
    
  end subroutine LDT_outputForcStatDscaleParams


end module LDT_statDscaleMod

