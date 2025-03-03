!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_statDscale_pluginMod
!BOP
!
! !MODULE: LDT_statDscale_pluginModMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines to scale meteorological forcing datasets 
!   either temporally and/or spatially (up or down). 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  11 Dec 2003  Sujay Kumar   Initial Specification
!  11 Dec 2014  KR Arsenault  Added C-function tables for downscaling
! 
!EOP  

  use LDT_pluginIndices

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LDT_statdscale_plugin

 contains

!
!BOP
! !ROUTINE: LDT_statdscale_plugin
!  \label{LDT_statdscale_plugin}
!
! !DESCRIPTION:
!
  subroutine LDT_statdscale_plugin
!
!EOP
    use forcingClimoMod
    use BayesianMergingMod

! !USES:

    call registerinitstatdscale(trim(LDT_bayesianMergeId)//char(0),&
         BayesianMerging_init)
    call registerdiagnosestatdscale(trim(LDT_bayesianMergeId)//char(0),&
         BayesianMerging_diagnose)

    call registerinitstatdscale(trim(LDT_forcingClimoId)//char(0),&
         forcingClimo_init)
    call registerdiagnosestatdscale(trim(LDT_forcingClimoId)//char(0),&
         forcingClimo_diagnose)
    call registercomputestatdscale(trim(LDT_forcingClimoId)//char(0),&
         forcingClimo_compute)
    call registeroutputstatdscale(trim(LDT_forcingClimoId)//char(0),&
         forcingClimo_output)

  end subroutine LDT_statdscale_plugin

end module LDT_statDscale_pluginMod
