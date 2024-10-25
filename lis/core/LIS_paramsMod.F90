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
module LIS_paramsMod
!BOP
!
! !MODULE: LIS_paramsMod
! 
! !DESCRIPTION:
!   The code in this file provides interfaces to manage different
!   sources of parameter datasets.
!   
!   \subsubsection{Overview}
!   This module contains interface plugins for the incorporating I/O 
!   associated with various sources of land surface parameter maps. 
!   These interfaces act as the entry points for handing issues such 
!   as format, ordering,and projection associated with a dataset. The
!   parameter reading routines are expected to perform all these 
!   transformations and convert the data to same grid, domain, 
!   and data ordering as the ones used in the running domain of LIS. 
!   The {\tt param\_module} provides interfaces to incorporate the 
!   following sources of parameter maps. 
!   \begin{description}
!   \item[soils]
!     sand, silt, clay fractions, soil color, soil texture
!   \item[albedo]
!     climatology and max albedo over deep snow
!   \item[greenness]
!     climatology
!   \item[Leaf Area Index]
!     climatology
!   \item[Stem Area Index]
!     climatology
!   \item[Topography]
!     static elevation, slope, aspect, curvature data
!   \item[Bottom Temperature]
!     static
!    \end{description}
!
! !REVISION HISTORY: 
!  26 Oct 2005    Sujay Kumar  Initial Specification
!  
  use LIS_LMLCMod
  use LIS_topoMod
  use LIS_vegDataMod
  use LIS_albedoMod
  use LIS_emissMod
  use LIS_soilsMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_param_init  !initializes structures, read static data
  public :: LIS_setDynParams !read time dependent data
  public :: LIS_param_finalize !cleanup allocated structures
  public :: LIS_param_reset
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

!EOP  
contains
!BOP
! !ROUTINE: LIS_param_init
! \label{LIS_param_init}
! 
! !INTERFACE:
  subroutine LIS_param_init()
! !USES:
    use LIS_coreMod,    only : LIS_rc
#ifdef ESMF_TRACE
    use ESMF
#endif
    
! !DESCRIPTION:
!  This interface provides the entry point for the initialization of 
!  data structures required for reading in parameter datasets. It 
!  also invokes routines to read static datasets such as soils. 
!
!  The calling sequence is: 
!  \begin{description}
!   \item[LIS\_albedo\_setup](\ref{LIS_albedo_setup}) \newline
!     call to setup albedo parameter options
!   \item[LIS\_greenness\_setup](\ref{LIS_greenness_setup}) \newline
!     call to setup greenness parameter options
!   \item[LIS\_lai\_setup](\ref{LIS_lai_setup}) \newline
!     call to setup LAI parameter options
!   \item[LIS\_sai\_setup](\ref{LIS_sai_setup}) \newline
!     call to setup SAI parameter options
!  \end{description}
!EOP
    TRACE_ENTER("param_init")
    call LIS_greenness_setup
    call LIS_roughness_setup
    call LIS_emiss_setup
    call LIS_albedo_setup
    call LIS_lai_setup
    call LIS_sai_setup
!    call LIS_prism_init()
    TRACE_EXIT("param_init")

  end subroutine LIS_param_init



!BOP
! !ROUTINE: LIS_setDynparams
! \label{LIS_setDynparams}
! 
! !INTERFACE:
  subroutine LIS_setDynparams(n)
! !USES:
    use LIS_coreMod,   only : LIS_rc
    use LIS_lsmMod,    only : LIS_setLSMDynparams
!    use prism_module, only : prism_dataops
#ifdef ESMF_TRACE
    use ESMF
#endif

! !ARGUMENTS: 
    integer, intent(in) :: n

! !DESCRIPTION:
!  This interface provides the entry point for routines to update
!  any time dependent land surface parameters. Current implementation
!  includes datasets such as albedo, greenness, LAI, and SAI. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!  The calling sequence is: 
!  \begin{description}
!   \item[LIS\_read\_greenness](\ref{LIS_read_greenness}) \newline
!     call to read the greenness data source
!   \item[LIS\_read\_albedo](\ref{LIS_read_albedo}) \newline
!     call to read the albedo data source
!   \item[LIS\_read\_lai](\ref{LIS_read_lai}) \newline
!     call to read the LAI data source
!   \item[LIS\_read\_sai](\ref{LIS_read_sai}) \newline
!     call to read the SAI data source
!  \end{description}
!EOP
    integer :: m 
    logical :: lsm_enabled

    TRACE_ENTER("param_dynsetup")
    lsm_enabled = .false. 
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          lsm_enabled = .true. 
       endif
    enddo

    call LIS_read_greenness(n)
    call LIS_read_emiss(n)
    call LIS_read_roughness(n)
    call LIS_read_albedo(n)
    call LIS_read_lai(n)
    call LIS_read_sai(n)
    if(lsm_enabled) then 
       call LIS_setLSMDynparams(n)
    endif
!    call prism_dataops(n)
    call diagnoseOutputparams(n)
    TRACE_EXIT("param_dynsetup")

  end subroutine LIS_setDynparams

!BOP
! !ROUTINE: diagnoseOutputparams
!  \label{diagnoseOutputparams}
! 
! !INTERFACE: 
  subroutine diagnoseOutputparams(n)
! !USES: 
    use LIS_coreMod,      only : LIS_rc
#ifdef ESMF_TRACE
    use ESMF
#endif

! !ARGUMENTS: 
    implicit none
    integer, intent(in) :: n 

! !DESCRIPTION: 
! This routine is used to enable the output of 
! land surface parameters to a file, once a month. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \end{description}
! 
!  The routines called are: 
!  \begin{description}
!   \item[LIS\_diagnoselandmask](\ref{LIS_diagnoselandmask}) \newline
!     maps the landmask data to the LIS history writer
!   \item[LIS\_diagnoseLandcover](\ref{LIS_diagnoseLandcover}) \newline
!     maps the landcover data to the LIS history writer
!   \item[LIS\_diagnoseTopography](\ref{LIS_diagnoseTopography}) \newline
!     maps the topography data to the LIS history writer
!   \item[LIS\_diagnosesoils](\ref{LIS_diagnosesoils}) \newline
!     maps the soils data to the LIS history writer
!   \item[LIS\_diagnosegfrac](\ref{LIS_diagnosegfrac}) \newline
!     maps the greenness data to the LIS history writer
!   \item[LIS\_diagnoseLAI](\ref{LIS_diagnoseLAI}) \newline
!     maps the LAI data to the LIS history writer
!   \item[LIS\_diagnoseSAI](\ref{LIS_diagnoseSAI}) \newline
!     maps the SAI data to the LIS history writer
!   \item[LIS\_diagnosealbedo](\ref{LIS_diagnosealbedo}) \newline
!     maps the albedo data to the LIS history writer
!  \end{description}
!EOP

    if(LIS_rc%wout.ne."none") then 
       call LIS_diagnoselandmask(n)
       call LIS_diagnoselandcover(n)
       call LIS_diagnosetopography(n)
       call LIS_diagnosesoils(n)
       call LIS_diagnosegfrac(n)
       call LIS_diagnoseemiss(n)
       call LIS_diagnoseroughness(n)
       call LIS_diagnoseLAI(n)
       call LIS_diagnoseSAI(n)
       call LIS_diagnosealbedo(n)
    endif

  end subroutine diagnoseOutputparams
    
    

!BOP
! !ROUTINE: LIS_param_finalize
! \label{LIS_param_finalize}
! 
! !INTERFACE:
  subroutine LIS_param_finalize()
! !USES:
    use LIS_coreMod, only  : LIS_rc
! !DESCRIPTION:
! This routine issues the invocation to deallocate and cleanup 
! any allocated data structures used in the parameter dataset
! implementations. 
!
!  The calling sequence is: 
!  \begin{description}
!   \item[LIS\_albedo\_finalize](\ref{LIS_albedo_finalize}) \newline
!    call to cleanup albedo related structures
!   \item[LIS\_greenness\_finalize](\ref{LIS_greenness_finalize}) \newline
!    call to cleanup greenness related structures
!   \item[LIS\_lai\_finalize](\ref{LIS_lai_finalize}) \newline
!    call to LAI related structures
!   \item[LIS\_sai\_finalize](\ref{LIS_sai_finalize}) \newline
!    call to SAI related structures
!  \end{description}
!EOP
    call LIS_albedo_finalize
    call LIS_greenness_finalize
    call LIS_emiss_finalize
    call LIS_roughness_finalize
    call LIS_lai_finalize
    call LIS_sai_finalize
  end subroutine LIS_param_finalize


!BOP
! !ROUTINE: LIS_param_reset
! \label{LIS_param_reset}
! 
! !INTERFACE:
  subroutine LIS_param_reset()
! !USES:
    use LIS_coreMod,   only : LIS_rc
#ifdef ESMF_TRACE
    use ESMF
#endif

! !ARGUMENTS: 

! !DESCRIPTION:
!  This interface provides the entry point for routines to update
!  any time dependent land surface parameters. Current implementation
!  includes datasets such as albedo, greenness, LAI, and SAI. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!  The calling sequence is: 
!  \begin{description}
!   \item[LIS\_greenness\_reset](\ref{LIS_greenness_reset}) \newline
!    resets the greenness data structures
!   \item[LIS\_albedo\_reset](\ref{LIS_albedo_reset}) \newline
!    resets the albedo data structures
!   \item[LIS\_lai\_reset](\ref{LIS_lai_reset}) \newline
!    resets the lai data structures
!   \item[LIS\_sai\_reset](\ref{LIS_sai_reset}) \newline
!    resets the sai data structures
!  \end{description}
!EOP
    TRACE_ENTER("param_reset")
    LIS_rc%rstflag = 1
    call LIS_greenness_reset()
    call LIS_emiss_reset
    call LIS_roughness_reset
    call LIS_albedo_reset()
    call LIS_lai_reset()
    call LIS_sai_reset()
    TRACE_EXIT("param_reset")

  end subroutine LIS_param_reset

end module LIS_paramsMod
