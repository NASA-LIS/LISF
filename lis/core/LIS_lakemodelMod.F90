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
module LIS_lakemodelMod
!BOP
!
! !MODULE: LIS_lakemodelMod
! 
! !DESCRIPTION:
!  The code in this file provides interfaces to manage the operation
!  of different lake models
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
  public :: LIS_lakemodel_init       ! initialize lakemodel variables, memory
  public :: LIS_setuplakemodel       ! set land surface parameters
  public :: LIS_lakemodel_readrestart    ! read the restart file
  public :: LIS_lakemodel_f2t     ! transfer forcing to model tiles
  public :: LIS_lakemodel_run       ! execute the land model 
  public :: LIS_lakemodel_output     ! write model output
  public :: LIS_lakemodel_setdynparams   ! set the time dependent parameters
  public :: LIS_lakemodel_writerestart   ! write the restart file
  public :: LIS_lakemodel_finalize   ! cleanup allocated structures
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------  
!EOP


contains

!BOP
! 
! !ROUTINE: LIS_lakemodel_init
! \label{LIS_lakemodel_init}
! 
! !INTERFACE:
  subroutine LIS_lakemodel_init()
! !USES:
    use LIS_surfaceModelDataMod
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecTile
    use LIS_logMod,       only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
!
! !DESCRIPTION:
! Interface for initializing the land model. The intialization includes
! the allocation of memory for LAKEMODEL specific variables and datastructures
! and specification of any runtime specific options. 
!
! The calling sequence is: 
! \begin{description}
!  \item[lakemodelinit] (\ref{lakemodelinit}) \newline
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

    TRACE_ENTER("lake_init")
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%lakemodel,&
         label="Lake model:",rc=rc)
    call LIS_verify(rc,'Lake model: option not specified in the config file')
    
    call lakemodelinit(trim(LIS_rc%lakemodel)//char(0), LIS_rc%lake_index)

    do n=1,LIS_rc%nnest
       LIS_sfmodel_struc(n)%models_used = &
            trim(LIS_sfmodel_struc(n)%models_used)//&
            "+"//trim(LIS_rc%lakemodel)
    enddo
    TRACE_EXIT("lake_init")
  end subroutine LIS_lakemodel_init

!BOP
! !ROUTINE: LIS_setuplakemodel
! \label{LIS_setuplakemodel}
!
! !INTERFACE:
  subroutine  LIS_setuplakemodel
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
! The setup interfaces are used for the LAKEMODEL specification of
! LAKEMODEL bparameters. If a parameter falls outside the LIS-specified
! generic parameter list, the LAKEMODEL is expected to provide routines for the 
! handling of any input data, specific to those parameters.
! 
! The calling sequence is: 
! \begin{description}
!  \item[lakemodelsetup] (\ref{lakemodelsetup}) \newline
!    invokes the generic method in the registry to set up the 
!    land surface model
! \end{description}
!EOP
    TRACE_ENTER("lake_setup")
    call lakemodelsetup(trim(LIS_rc%lakemodel)//char(0),LIS_rc%lake_index)
    TRACE_EXIT("lake_setup")
  end subroutine LIS_setuplakemodel

!BOP
! !ROUTINE: LIS_lakemodel_run
! \label{LIS_lakemodel_run}
!
! !INTERFACE:            
  subroutine LIS_lakemodel_run(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n 
!
! !DESCRIPTION:
! This interface provides the entry point to the LAKEMODEL routines that
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
!  \item[lakemodelrun] (\ref{lakemodelrun}) \newline
!    invokes the generic method in the registry to run the
!    land surface model
! \end{description}
!EOP
    TRACE_ENTER("lake_run")
    call lakemodelrun(trim(LIS_rc%lakemodel)//char(0), n,LIS_rc%lake_index)
    TRACE_EXIT("lake_run")

  end subroutine LIS_lakemodel_run

!BOP
! !ROUTINE: LIS_lakemodel_readrestart
! \label{LIS_lakemodel_readrestart}
!
! !INTERFACE:
  subroutine LIS_lakemodel_readrestart
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
!  This interface provides the entry point to read LAKEMODEL specific 
!  model restart files. 
!
! The calling sequence is: 
! \begin{description}
!  \item[lakemodelrestart] (\ref{lakemodelrestart}) \newline
!    invokes the generic method in the registry to read
!    restart files for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("lake_readrst")
    call lakemodelrestart(trim(LIS_rc%lakemodel)//char(0),LIS_rc%lake_index)
    TRACE_EXIT("lake_readrst")
  end subroutine LIS_lakemodel_readrestart
!BOP
! !ROUTINE: LIS_lakemodel_output
! \label{LIS_lakemodel_output}
!
! !INTERFACE:
  subroutine LIS_lakemodel_output(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This interface provides the entry point to write LAKEMODEL specific 
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
!  \item[lakemodeloutput] (\ref{lakemodeloutput}) \newline
!    invokes the generic method in the registry to write the
!    land surface model output
! \end{description}
!EOP
    TRACE_ENTER("lake_out")
    call lakemodeloutput(trim(LIS_rc%lakemodel)//char(0),n,LIS_rc%lake_index)
    TRACE_EXIT("lake_out")
  end subroutine LIS_lakemodel_output

!BOP
! !ROUTINE: LIS_lakemodel_setdynparams
! \label{LIS_lakemodel_setdynparams}
!
! !INTERFACE:
  subroutine LIS_lakemodel_setdynparams(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n

! 
! !DESCRIPTION:
! This interface provides an entry point to update any 
! time dependent land surface model parameters in an LAKEMODEL. 
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[lakemodeldynsetup] (\ref{lakemodeldynsetup}) \newline
!    invokes the generic method in the registry to set up 
!    time dependent parameters for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("lake_dynsetup")
    call lakemodeldynsetup(trim(LIS_rc%lakemodel)//char(0),n,LIS_rc%lake_index)
    TRACE_EXIT("lake_dynsetup")
  end subroutine LIS_lakemodel_setdynparams

!BOP
! !ROUTINE: LIS_lakemodel_f2t
! \label{LIS_lakemodel_f2t}
!
! !INTERFACE:
  subroutine LIS_lakemodel_f2t(n)
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
!  \item[lakemodelf2t] (\ref{lakemodelf2t}) \newline
!    invokes the generic method in the registry to transfer the
!    forcing to the land surface model tiles
! \end{description}
!EOP
    TRACE_ENTER("lake_f2t")
    call lakemodelf2t(trim(LIS_rc%lakemodel)//"+"//trim(LIS_rc%runmode)//char(0),&
         n,LIS_rc%lake_index)
    TRACE_EXIT("lake_f2t")

  end subroutine LIS_lakemodel_f2t


!BOP
! !ROUTINE: LIS_lakemodel_writerestart
! \label{LIS_lakemodel_writerestart}
! 
! !INTERFACE:
  subroutine LIS_lakemodel_writerestart(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n 
! !DESCRIPTION:
!  This interface provides the entry point to read the LAKEMODEL specific 
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
!  \item[lakemodelwrst] (\ref{lakemodelwrst}) \newline
!    invokes the generic method in the registry to write
!    restart files for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("lake_writerst")
    call lakemodelwrst(trim(LIS_rc%lakemodel)//char(0),n)
    TRACE_EXIT("lake_writerst")
  end subroutine LIS_lakemodel_writerestart

!BOP
! !ROUTINE: LIS_lakemodel_finalize
! \label{LIS_lakemodel_finalize}
!
! !INTERFACE:
  subroutine LIS_lakemodel_finalize()
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
!  \item[lakemodelfinalize] (\ref{lakemodelfinalize}) \newline
!    invokes the generic method in the registry to cleanup the 
!    LAKEMODEL related datastructures    
! \end{description}
!EOP
    call lakemodelfinalize(trim(LIS_rc%lakemodel)//char(0))
  end subroutine LIS_lakemodel_finalize

end module LIS_lakemodelMod
