!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_openwatermodelMod
!BOP
!
! !MODULE: LIS_openwatermodelMod
! 
! !DESCRIPTION:
!  The code in this file provides interfaces to manage the operation
!  of different ocean models
! 
! !REVISION HISTORY: 
!  16 Jul 2012    Sujay Kumar  Initial Specification
! 
  use ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_openwatermodel_init       ! initialize openwater variables, memory
  public :: LIS_setupopenwatermodel       ! set land surface parameters
  public :: LIS_openwatermodel_readrestart    ! read the restart file
  public :: LIS_openwatermodel_f2t     ! transfer forcing to model tiles
  public :: LIS_openwatermodel_run       ! execute the land model 
  public :: LIS_openwatermodel_output     ! write model output
  public :: LIS_openwatermodel_setdynparams   ! set the time dependent parameters
  public :: LIS_openwatermodel_writerestart   ! write the restart file
  public :: LIS_openwatermodel_finalize   ! cleanup allocated structures
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------  
!EOP


contains

!BOP
! 
! !ROUTINE: LIS_openwatermodel_init
! \label{LIS_openwatermodel_init}
! 
! !INTERFACE:
  subroutine LIS_openwatermodel_init()
! !USES:
    use LIS_surfaceModelDataMod
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecTile
    use LIS_openwatermodel_pluginMod,    only : LIS_openwatermodel_plugin
    use LIS_logMod,       only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
!
! !DESCRIPTION:
! Interface for initializing the land model. The intialization includes
! the allocation of memory for OPENWATER specific variables and datastructures
! and specification of any runtime specific options. 
!
! The calling sequence is: 
! \begin{description}
!  \item[LIS\_openwatermodel\_plugin] (\ref{LIS_openwatermodel_plugin}) \newline
!    sets up function table registries for implemented land surface models
!  \item[openwaterinit] (\ref{openwaterinit}) \newline
!    invokes the generic method in the registry to initialize the 
!    land surface model
! \end{description}
!EOP
    integer       :: n 
    integer       :: status
    character*1   :: nestid(2)
    character*1   :: caseid(3)
    character*100 :: temp
    integer       :: ftn

    integer              :: i,j,k
    integer              :: rc

    TRACE_ENTER("owater_init")
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%openwatermodel,&
         label="Open water model:",rc=rc)
    call LIS_verify(rc,'Open water model: option not specified in the config file')
    
    call LIS_openwatermodel_plugin

    call openwaterinit(trim(LIS_rc%openwatermodel)//char(0), LIS_rc%openwater_index)
    
    do n=1,LIS_rc%nnest
       LIS_sfmodel_struc(n)%models_used = &
            trim(LIS_sfmodel_struc(n)%models_used)//&
            "+"//trim(LIS_rc%openwatermodel)
    enddo
    TRACE_EXIT("owater_init")

  end subroutine LIS_openwatermodel_init

!BOP
! !ROUTINE: LIS_setupopenwatermodel
! \label{LIS_setupopenwatermodel}
!
! !INTERFACE:
  subroutine  LIS_setupopenwatermodel
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
! The setup interfaces are used for the OPENWATER specification of
! OPENWATER bparameters. If a parameter falls outside the LIS-specified
! generic parameter list, the OPENWATER is expected to provide routines for the 
! handling of any input data, specific to those parameters.
! 
! The calling sequence is: 
! \begin{description}
!  \item[openwatersetup] (\ref{openwatersetup}) \newline
!    invokes the generic method in the registry to set up the 
!    land surface model
! \end{description}
!EOP
    TRACE_ENTER("owater_setup")
    call openwatersetup(trim(LIS_rc%openwatermodel)//char(0),LIS_rc%openwater_index)
    TRACE_EXIT("owater_setup")
  end subroutine LIS_setupopenwatermodel
  
!BOP
! !ROUTINE: LIS_openwatermodel_run
! \label{LIS_openwatermodel_run}
!
! !INTERFACE:            
  subroutine LIS_openwatermodel_run(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n 
!
! !DESCRIPTION:
! This interface provides the entry point to the OPENWATER routines that
! invokes the land surface model physics. 
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[openwaterrun] (\ref{openwaterrun}) \newline
!    invokes the generic method in the registry to run the
!    land surface model
! \end{description}
!EOP

    TRACE_ENTER("owater_run")
    call openwaterrun(trim(LIS_rc%openwatermodel)//char(0), n,LIS_rc%openwater_index)
    TRACE_EXIT("owater_run")

  end subroutine LIS_openwatermodel_run

!BOP
! !ROUTINE: LIS_openwatermodel_readrestart
! \label{LIS_openwatermodel_readrestart}
!
! !INTERFACE:
  subroutine LIS_openwatermodel_readrestart
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
!  This interface provides the entry point to read OPENWATER specific 
!  model restart files. 
!
! The calling sequence is: 
! \begin{description}
!  \item[openwaterrestart] (\ref{openwaterrestart}) \newline
!    invokes the generic method in the registry to read
!    restart files for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("owater_readrst")
    call openwaterrestart(trim(LIS_rc%openwatermodel)//char(0),LIS_rc%openwater_index)
    TRACE_EXIT("owater_readrst")
  end subroutine LIS_openwatermodel_readrestart
!BOP
! !ROUTINE: LIS_openwatermodel_output
! \label{LIS_openwatermodel_output}
!
! !INTERFACE:
  subroutine LIS_openwatermodel_output(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This interface provides the entry point to write OPENWATER specific 
!  model output files. 
!
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[openwateroutput] (\ref{openwateroutput}) \newline
!    invokes the generic method in the registry to write the
!    land surface model output
! \end{description}
!EOP
    TRACE_ENTER("owater_out")
    call openwateroutput(trim(LIS_rc%openwatermodel)//char(0),n,LIS_rc%openwater_index)
    TRACE_EXIT("owater_out")
  end subroutine LIS_openwatermodel_output

!BOP
! !ROUTINE: LIS_openwatermodel_setdynparams
! \label{LIS_openwatermodel_setdynparams}
!
! !INTERFACE:
  subroutine LIS_openwatermodel_setdynparams(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n

! 
! !DESCRIPTION:
! This interface provides an entry point to update any 
! time dependent land surface model parameters in an OPENWATER. 
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[openwaterdynsetup] (\ref{openwaterdynsetup}) \newline
!    invokes the generic method in the registry to set up 
!    time dependent parameters for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("owater_dynsetup")
    call openwaterdynsetup(trim(LIS_rc%openwatermodel)//char(0),n,LIS_rc%openwater_index)
    TRACE_EXIT("owater_dynsetup")
  end subroutine LIS_openwatermodel_setdynparams

!BOP
! !ROUTINE: LIS_openwatermodel_f2t
! \label{LIS_openwatermodel_f2t}
!
! !INTERFACE:
  subroutine LIS_openwatermodel_f2t(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n
! 
! !DESCRIPTION:
!  
!  This interface is used to transfer the forcing variables to the actual 
!  model tile space. Any forcing perturbations that need to be applied to 
!  the input forcing is applied at this stage as well. 
!
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[openwaterf2t] (\ref{openwaterf2t}) \newline
!    invokes the generic method in the registry to transfer the
!    forcing to the land surface model tiles
! \end{description}
!EOP

    TRACE_ENTER("owater_f2t")
    call openwaterf2t(trim(LIS_rc%openwatermodel)//"+"//trim(LIS_rc%runmode)//char(0),&
         n,LIS_rc%openwater_index)
    TRACE_EXIT("owater_f2t")

  end subroutine LIS_openwatermodel_f2t


!BOP
! !ROUTINE: LIS_openwatermodel_writerestart
! \label{LIS_openwatermodel_writerestart}
! 
! !INTERFACE:
  subroutine LIS_openwatermodel_writerestart(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n 
! !DESCRIPTION:
!  This interface provides the entry point to read the OPENWATER specific 
!  model restart files. 
!
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[openwaterwrst] (\ref{openwaterwrst}) \newline
!    invokes the generic method in the registry to write
!    restart files for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("owater_writerst")
    call openwaterwrst(trim(LIS_rc%openwatermodel)//char(0),n)
    TRACE_EXIT("owater_writerst")
  end subroutine LIS_openwatermodel_writerestart

!BOP
! !ROUTINE: LIS_openwatermodel_finalize
! \label{LIS_openwatermodel_finalize}
!
! !INTERFACE:
  subroutine LIS_openwatermodel_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures in the specific instance of a 
!  land surface model 
!
! The calling sequence is: 
! \begin{description}
!  \item[openwaterfinalize] (\ref{openwaterfinalize}) \newline
!    invokes the generic method in the registry to cleanup the 
!    OPENWATER related datastructures    
! \end{description}
!EOP
    call openwaterfinalize(trim(LIS_rc%openwatermodel)//char(0))
  end subroutine LIS_openwatermodel_finalize

end module LIS_openwatermodelMod
