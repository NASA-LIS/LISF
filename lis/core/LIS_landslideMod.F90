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
module LIS_landslideMod
!BOP
!
! !MODULE: LIS_landslideMod
! 
! !DESCRIPTION:
!  
!
! !REVISION HISTORY: 
! 
!  

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_initLandSlideModel
  PUBLIC :: LIS_runLandSlideModel
  PUBLIC :: LIS_outputLandSlideModel
  PUBLIC :: LIS_finalizeLandSlideModel ! SY
  PUBLIC :: LIS_landslide_param_reset
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

!EOP  
contains
!BOP
! !ROUTINE: LIS_initLandSlideModel
! \label{LIS_initLandSlideModel}
! 
! !INTERFACE:
  subroutine LIS_initLandSlideModel()
! !USES:
    use ESMF
    use LIS_coreMod,    only : LIS_rc, LIS_config
    use LIS_logMod,     only : LIS_verify

    implicit none
! !DESCRIPTION:
!EOP
    integer :: rc

    TRACE_ENTER("lslide_init")
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%landslidemodel,&
         label="land slide model:",rc=rc)
    call LIS_verify(rc,'land slide model: option not specified in the config file')

    !call initializelandslidemodel(LIS_rc%landslidemodel) ! SY
    call initializelandslidemodel(trim(LIS_rc%landslidemodel)//char(0)) ! SY
    TRACE_EXIT("lslide_init")

  end subroutine LIS_initLandSlideModel



!BOP
! !ROUTINE: LIS_runLandSlideModel
! \label{LIS_runLandSlideModel}
! 
! !INTERFACE:
  subroutine LIS_runLandSlideModel(n)
! !USES:
    use LIS_coreMod,   only : LIS_rc
#ifdef ESMF_TRACE
    use ESMF
#endif

! !ARGUMENTS: 
    integer, intent(in) :: n

! !DESCRIPTION:
!
!EOP
    TRACE_ENTER("lslide_run")
    !call runlandslidemodel(LIS_rc%landslidemodel, n) ! SY
    call runlandslidemodel(trim(LIS_rc%landslidemodel)//char(0), n) ! SY
    TRACE_EXIT("lslide_run")

  end subroutine LIS_runLandSlideModel

!BOP
! !ROUTINE: LIS_outputLandSlideModel
! \label{LIS_outputLandSlideModel}
! 
! !INTERFACE:
  subroutine LIS_outputLandSlideModel(n)
! !USES:
    use LIS_coreMod,   only : LIS_rc
#ifdef ESMF_TRACE
    use ESMF
#endif

! !ARGUMENTS: 
    integer, intent(in) :: n

! !DESCRIPTION:
!
!EOP
    TRACE_ENTER("lslide_out")
    !call outputlandslidemodel(LIS_rc%landslidemodel, n) ! SY
    call outputlandslidemodel(trim(LIS_rc%landslidemodel)//char(0), n) ! SY
    TRACE_EXIT("lslide_out")

  end subroutine LIS_outputLandSlideModel

! SY: Begin
!BOP
! !ROUTINE: LIS_finalizeLandSlideModel
! \label{LIS_finalizeLandSlideModel}
! 
! !INTERFACE:
  subroutine LIS_finalizeLandSlideModel()
! !USES:
    use LIS_coreMod,   only : LIS_rc
! !DESCRIPTION:
!
!EOP
    !call landslidemodelfinalize(LIS_rc%landslidemodel) ! SY
    call landslidemodelfinalize(trim(LIS_rc%landslidemodel)//char(0)) ! SY

  end subroutine LIS_finalizeLandSlideModel
! SY: End

  subroutine LIS_landslide_param_reset
    use LIS_coreMod, only : LIS_rc
#ifdef ESMF_TRACE
    use ESMF
#endif

    TRACE_ENTER("lslide_reset")
    !call resetlandslidemodel(LIS_rc%landslidemodel) ! SY
    call resetlandslidemodel(trim(LIS_rc%landslidemodel)//char(0)) ! SY
    TRACE_EXIT("lslide_reset")

  end subroutine LIS_landslide_param_reset

end module LIS_landslideMod
