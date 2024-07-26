!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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

    external LDT_init_climoRstproc
    external LDT_run_climoRstproc

    external LDT_init_rstTransformProc
    external LDT_run_rstTransformProc

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

    external LDT_init_usafsi
    external LDT_run_usafsi

    external LDT_init_OPTUEparamproc
    external LDT_run_OPTUEparamproc

    external LDT_init_obsSim
    external LDT_run_obsSim

    external LDT_init_LISHydropreproc
    external LDT_run_LISHydropreproc

    external LDT_init_smap_e_opl     !Y.Kwon
    external LDT_run_smap_e_opl      !Y.Kwon

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

  ! climatological Restart processing:
    call registerldtinit(trim(LDT_climoRstProcId)//char(0), &
         LDT_init_climoRstproc)
    call registerldtrun(trim(LDT_climorstProcId)//char(0), &
         LDT_run_climoRstproc)

  ! Restart transformation processing:
    call registerldtinit(trim(LDT_rstTransformProcId)//char(0), &
         LDT_init_rstTransformproc)
    call registerldtrun(trim(LDT_rstTransformProcId)//char(0), &
         LDT_run_rstTransformproc)

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

    ! USAFSI analysis
    call registerldtinit(trim(LDT_usafsiId)//char(0), &
         LDT_init_usafsi)
    call registerldtrun(trim(LDT_usafsiId)//char(0), &
         LDT_run_usafsi)

    ! OPTUE processing
    call registerldtinit(trim(LDT_OPTUEparamprocId)//char(0), &
         LDT_init_OPTUEparamproc)
    call registerldtrun(trim(LDT_OPTUEparamprocId)//char(0), &
         LDT_run_OPTUEparamproc)

    ! obs simulator
    call registerldtinit(trim(LDT_obsSimprocId)//char(0), &
         LDT_init_obsSim)
    call registerldtrun(trim(LDT_obsSimprocId)//char(0), &
         LDT_run_obsSim)

  ! LISHydro Preprocessing for WRFHydro:
    call registerldtinit(trim(LDT_LISHydropreprocId)//char(0), &
         LDT_init_LISHydropreproc)
    call registerldtrun(trim(LDT_LISHydropreprocId)//char(0), &
         LDT_run_LISHydropreproc)

  ! OPL E SMAP soil moisture retrieval  (Y.Kwon)
    call registerldtinit(trim(LDT_SMAP_E_OPLId)//char(0), &
         LDT_init_smap_e_opl)
    call registerldtrun(trim(LDT_SMAP_E_OPLId)//char(0), &
         LDT_run_smap_e_opl)

  end subroutine LDT_runmode_plugin

end module LDT_runmode_pluginMod
