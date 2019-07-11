!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_runmode_pluginMod
!BOP
!
! !MODULE: LDT_runmode_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various LDT-runmodes. 
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
  PUBLIC :: LDT_runmode_plugin  
contains
!BOP
! !ROUTINE: LDT_runmode_plugin
!  \label{LDT_runmode_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LDT-runmode. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LDT-runmode. 
!
!
! !INTERFACE:
  subroutine LDT_runmode_plugin
    use LDT_pluginIndices

    implicit none
!EOP
    external LDT_init_LSMparamproc
    external LDT_run_LSMparamproc

    external LDT_init_DApreproc
    external LDT_run_DApreproc

    external LDT_init_EnsRstpreproc
    external LDT_run_EnsRstpreproc

    external LDT_init_Rstproc
    external LDT_run_Rstproc

    external LDT_init_MetforcProc
    external LDT_run_MetforcProc
    external LDT_final_MetforcProc

    external LDT_init_MetTimeDScale
    external LDT_run_MetTimeDScale
!    external LDT_final_MetTimeDScale

    external LDT_init_StatDscaleMetForc
    external LDT_run_StatDscaleMetForc
    
    external LDT_init_NUWRFpreproc
    external LDT_run_NUWRFpreproc

    external LDT_init_ANNproc
    external LDT_run_ANNproc

    external LDT_init_ldtsi
    external LDT_run_ldtsi

    external LDT_init_OPTUEparamproc
    external LDT_run_OPTUEparamproc

  ! Parameter Preprocessing:
    call registerldtinit(trim(LDT_LSMparamprocId)//char(0), &
         LDT_init_LSMparamproc)
    call registerldtrun(trim(LDT_LSMparamprocId)//char(0), &
         LDT_run_LSMparamproc)

  ! Data Assimilation CDF-Stats Preprocessing:
    call registerldtinit(trim(LDT_DApreprocId)//char(0), LDT_init_DApreproc)
    call registerldtrun(trim(LDT_DApreprocId)//char(0), LDT_run_DApreproc)

  ! Ensemble Restart Preprocessing:
    call registerldtinit(trim(LDT_EnsRstpreprocId)//char(0), &
         LDT_init_EnsRstpreproc)
    call registerldtrun(trim(LDT_EnsRstpreprocId)//char(0), &
         LDT_run_EnsRstpreproc)

  ! Restart processing:
    call registerldtinit(trim(LDT_rstProcId)//char(0), &
         LDT_init_Rstproc)
    call registerldtrun(trim(LDT_rstProcId)//char(0), &
         LDT_run_Rstproc)

  ! Meteorological Forcing Processing Only:
    call registerldtinit(trim(LDT_MetForcprocId)//char(0), &
         LDT_init_MetforcProc)
    call registerldtrun(trim(LDT_MetForcprocId)//char(0), &
         LDT_run_MetforcProc)
    call registerldtrun(trim(LDT_MetForcprocId)//char(0), &
         LDT_final_MetforcProc)

  ! Meteorological Forcing Temporal Downscaling Preprocessing:
    call registerldtinit(trim(LDT_MetTDscaleprocId)//char(0), &
         LDT_init_MetTimeDScale)
    call registerldtrun(trim(LDT_MetTDscaleprocId)//char(0), &
         LDT_run_MetTimeDScale)

    call registerldtinit(trim(LDT_StatDscaleMetforcprocId)//char(0), &
         LDT_init_StatDscaleMetForc)
    call registerldtrun(trim(LDT_StatDscaleMetforcprocId)//char(0), &
         LDT_run_StatDscaleMetForc)

  ! NUWRF Preprocessing:
    call registerldtinit(trim(LDT_NUWRFpreprocId)//char(0), &
         LDT_init_NUWRFpreproc)
    call registerldtrun(trim(LDT_NUWRFpreprocId)//char(0), &
         LDT_run_NUWRFpreproc)
  
  ! Artificial Neural Network Approach:
    call registerldtinit(trim(LDT_ANNprocId)//char(0), &
         LDT_init_ANNproc)
    call registerldtrun(trim(LDT_ANNprocId)//char(0), &
         LDT_run_ANNproc)

    ! LDTSI analysis
    call registerldtinit(trim(LDT_ldtsiId)//char(0), &
         LDT_init_ldtsi)
    call registerldtrun(trim(LDT_ldtsiId)//char(0), &
         LDT_run_ldtsi)

    ! OPTUE processing
    call registerldtinit(trim(LDT_OPTUEparamprocId)//char(0), &
         LDT_init_OPTUEparamproc)
    call registerldtrun(trim(LDT_OPTUEparamprocId)//char(0), &
         LDT_run_OPTUEparamproc)

  end subroutine LDT_runmode_plugin

end module LDT_runmode_pluginMod
