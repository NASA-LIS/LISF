!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_appMod
!BOP
!
! !MODULE: LIS_appMod
! 
! !DESCRIPTION:
!  
!  This module controls the operation of different application models 
!  within LIS. An application model is defined one that does not have 
!  any feedback to the land surface states. 
!
! !REVISION HISTORY: 
! 
!  
  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_appModel_init
  PUBLIC :: LIS_runAppModel
  PUBLIC :: LIS_outputAppModel
  PUBLIC :: LIS_appModel_finalize ! SY
  PUBLIC :: LIS_app_param_reset
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

!EOP  
contains
!BOP
! !ROUTINE: LIS_appModel_init
! \label{LIS_appModel_init}
! 
! !INTERFACE:
  subroutine LIS_appModel_init()
! !USES:
    use ESMF
    use LIS_coreMod,    only : LIS_rc, LIS_config
    use LIS_logMod,     only : LIS_verify
    use LIS_landslideMod, only : LIS_initLandSlideModel
    implicit none
! !DESCRIPTION:
!
! Reads the configuration options related to application models, 
! invokes the initialization of specific application models, 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[LIS\_initLandSlideModel](\ref{LIS_initLandSlideModel}) \newline
!    invokes the initialization of landslide models. 
!  \end{description}
!EOP
    integer :: rc

    TRACE_ENTER("appMod_init")
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nappmodel,&
         label="Number of application models:",rc=rc)
    call LIS_verify(rc,&
         'Number of application models: option not specified in the config file')
    if(LIS_rc%nappmodel.gt.0) then 
! call each application model types. 

       call LIS_initLandSlideModel
       
    endif
    TRACE_EXIT("appMod_init")
  end subroutine LIS_appModel_init



!BOP
! !ROUTINE: LIS_runAppModel
! \label{LIS_runAppModel}
! 
! !INTERFACE:
  subroutine LIS_runAppModel(n)
! !USES:
    use LIS_coreMod,   only : LIS_rc
    use LIS_landslideMod,     only : LIS_runLandSlideModel
#ifdef ESMF_TRACE
    use ESMF
#endif

! !ARGUMENTS: 
    integer, intent(in) :: n

! !DESCRIPTION:
!   This routine invokes the execution of specified application 
!   models. 
!   
!EOP
    TRACE_ENTER("appMod_run")
    if(LIS_rc%nappmodel.gt.0) then 
       call LIS_runLandslideModel(n)
    endif
    TRACE_EXIT("appMod_run")

  end subroutine LIS_runAppModel

!BOP
! !ROUTINE: LIS_outputAppModel(n)
! \label{LIS_outputAppModel}
! 
! !INTERFACE:
  subroutine LIS_outputAppModel(n)
! !USES:
    use LIS_coreMod,   only : LIS_rc
    use LIS_landslideMod,     only : LIS_outputLandSlideModel
#ifdef ESMF_TRACE
    use ESMF
#endif

! !ARGUMENTS: 
    integer, intent(in) :: n

! !DESCRIPTION:
!
!EOP
    TRACE_ENTER("appMod_out")
    if(LIS_rc%nappmodel.gt.0) then 
       call LIS_outputLandslideModel(n)
    endif
    TRACE_EXIT("appMod_out")

  end subroutine LIS_outputAppModel

! SY: Begin
!BOP
! !ROUTINE: LIS_appModel_finalize()
! \label{LIS_appModel_finalize}
! 
! !INTERFACE:
  subroutine LIS_appModel_finalize()
! !USES:
    use LIS_coreMod,   only : LIS_rc
    use LIS_landslideMod,     only : LIS_finalizeLandSlideModel
! !DESCRIPTION:
!
!EOP
    if(LIS_rc%nappmodel.gt.0) then 
       call LIS_finalizeLandSlideModel
    endif

  end subroutine LIS_appModel_finalize
! SY: End

  subroutine LIS_app_param_reset
    use LIS_coreMod, only : LIS_rc
    use LIS_landslideMod, only : LIS_landslide_param_reset
#ifdef ESMF_TRACE
    use ESMF
#endif

    TRACE_ENTER("appMod_reset")
    if(LIS_rc%nappmodel.gt.0) then 
       call LIS_landslide_param_reset
    endif
    TRACE_EXIT("appMod_reset")

  end subroutine LIS_app_param_reset

end module LIS_appMod
