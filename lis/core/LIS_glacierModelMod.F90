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
module LIS_glaciermodelMod
!BOP
!
! !MODULE: LIS_glaciermodelMod
! 
! !DESCRIPTION:
!  The code in this file provides interfaces to manage the operation
!  of different glacier models
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
  public :: LIS_glaciermodel_init       ! initialize glaciermodel variables, memory
  public :: LIS_setupglaciermodel       ! set land surface parameters
  public :: LIS_glaciermodel_readrestart    ! read the restart file
  public :: LIS_glaciermodel_f2t     ! transfer forcing to model tiles
  public :: LIS_glaciermodel_run       ! execute the land model 
  public :: LIS_glaciermodel_output     ! write model output
  public :: LIS_glaciermodel_setdynparams   ! set the time dependent parameters
  public :: LIS_glaciermodel_writerestart   ! write the restart file
  public :: LIS_glaciermodel_finalize   ! cleanup allocated structures
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------  
!EOP


contains

!BOP
! 
! !ROUTINE: LIS_glaciermodel_init
! \label{LIS_glaciermodel_init}
! 
! !INTERFACE:
  subroutine LIS_glaciermodel_init()
! !USES:
    use LIS_surfaceModelDataMod
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecTile
    use LIS_logMod,       only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
!
! !DESCRIPTION:
! Interface for initializing the land model. The intialization includes
! the allocation of memory for GLACIERMODEL specific variables and datastructures
! and specification of any runtime specific options. 
!
! The calling sequence is: 
! \begin{description}
!  \item[glaciermodelinit] (\ref{glaciermodelinit}) \newline
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

    TRACE_ENTER("glacier_init")
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%glaciermodel,&
         label="Glacier model:",rc=rc)
    call LIS_verify(rc,'Glacier model: option not specified in the config file')
    
    call glaciermodelinit(trim(LIS_rc%glaciermodel)//char(0), LIS_rc%glacier_index)

    do n=1,LIS_rc%nnest
       LIS_sfmodel_struc(n)%models_used = &
            trim(LIS_sfmodel_struc(n)%models_used)//&
            "+"//trim(LIS_rc%glaciermodel)
    enddo
    TRACE_EXIT("glacier_init")
  end subroutine LIS_glaciermodel_init

!BOP
! !ROUTINE: LIS_setupglaciermodel
! \label{LIS_setupglaciermodel}
!
! !INTERFACE:
  subroutine  LIS_setupglaciermodel
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
! The setup interfaces are used for the GLACIERMODEL specification of
! GLACIERMODEL bparameters. If a parameter falls outside the LIS-specified
! generic parameter list, the GLACIERMODEL is expected to provide routines for the 
! handling of any input data, specific to those parameters.
! 
! The calling sequence is: 
! \begin{description}
!  \item[glaciermodelsetup] (\ref{glaciermodelsetup}) \newline
!    invokes the generic method in the registry to set up the 
!    land surface model
! \end{description}
!EOP
    TRACE_ENTER("glacier_setup")
    call glaciermodelsetup(trim(LIS_rc%glaciermodel)//char(0),LIS_rc%glacier_index)
    TRACE_EXIT("glacier_setup")
  end subroutine LIS_setupglaciermodel

!BOP
! !ROUTINE: LIS_glaciermodel_run
! \label{LIS_glaciermodel_run}
!
! !INTERFACE:            
  subroutine LIS_glaciermodel_run(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n 
!
! !DESCRIPTION:
! This interface provides the entry point to the GLACIERMODEL routines that
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
!  \item[glaciermodelrun] (\ref{glaciermodelrun}) \newline
!    invokes the generic method in the registry to run the
!    land surface model
! \end{description}
!EOP
    TRACE_ENTER("glacier_run")
    call glaciermodelrun(trim(LIS_rc%glaciermodel)//char(0), n,LIS_rc%glacier_index)
    TRACE_EXIT("glacier_run")

  end subroutine LIS_glaciermodel_run

!BOP
! !ROUTINE: LIS_glaciermodel_readrestart
! \label{LIS_glaciermodel_readrestart}
!
! !INTERFACE:
  subroutine LIS_glaciermodel_readrestart
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
!  This interface provides the entry point to read GLACIERMODEL specific 
!  model restart files. 
!
! The calling sequence is: 
! \begin{description}
!  \item[glaciermodelrestart] (\ref{glaciermodelrestart}) \newline
!    invokes the generic method in the registry to read
!    restart files for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("glacier_readrst")
    call glaciermodelrestart(trim(LIS_rc%glaciermodel)//char(0),LIS_rc%glacier_index)
    TRACE_EXIT("glacier_readrst")
  end subroutine LIS_glaciermodel_readrestart
!BOP
! !ROUTINE: LIS_glaciermodel_output
! \label{LIS_glaciermodel_output}
!
! !INTERFACE:
  subroutine LIS_glaciermodel_output(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This interface provides the entry point to write GLACIERMODEL specific 
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
!  \item[glaciermodeloutput] (\ref{glaciermodeloutput}) \newline
!    invokes the generic method in the registry to write the
!    land surface model output
! \end{description}
!EOP
    TRACE_ENTER("glacier_out")
    call glaciermodeloutput(trim(LIS_rc%glaciermodel)//char(0),n,LIS_rc%glacier_index)
    TRACE_EXIT("glacier_out")
  end subroutine LIS_glaciermodel_output

!BOP
! !ROUTINE: LIS_glaciermodel_setdynparams
! \label{LIS_glaciermodel_setdynparams}
!
! !INTERFACE:
  subroutine LIS_glaciermodel_setdynparams(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n

! 
! !DESCRIPTION:
! This interface provides an entry point to update any 
! time dependent land surface model parameters in an GLACIERMODEL. 
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[glaciermodeldynsetup] (\ref{glaciermodeldynsetup}) \newline
!    invokes the generic method in the registry to set up 
!    time dependent parameters for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("glacier_dynsetup")
    call glaciermodeldynsetup(trim(LIS_rc%glaciermodel)//char(0),n,LIS_rc%glacier_index)
    TRACE_EXIT("glacier_dynsetup")
  end subroutine LIS_glaciermodel_setdynparams

!BOP
! !ROUTINE: LIS_glaciermodel_f2t
! \label{LIS_glaciermodel_f2t}
!
! !INTERFACE:
  subroutine LIS_glaciermodel_f2t(n)
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
!  \item[glaciermodelf2t] (\ref{glaciermodelf2t}) \newline
!    invokes the generic method in the registry to transfer the
!    forcing to the land surface model tiles
! \end{description}
!EOP
    TRACE_ENTER("glacier_f2t")
    call glaciermodelf2t(trim(LIS_rc%glaciermodel)//"+"//trim(LIS_rc%runmode)//char(0),&
         n,LIS_rc%glacier_index)
    TRACE_EXIT("glacier_f2t")

  end subroutine LIS_glaciermodel_f2t


!BOP
! !ROUTINE: LIS_glaciermodel_writerestart
! \label{LIS_glaciermodel_writerestart}
! 
! !INTERFACE:
  subroutine LIS_glaciermodel_writerestart(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer, intent(in) :: n 
! !DESCRIPTION:
!  This interface provides the entry point to read the GLACIERMODEL specific 
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
!  \item[glaciermodelwrst] (\ref{glaciermodelwrst}) \newline
!    invokes the generic method in the registry to write
!    restart files for the land surface model 
! \end{description}
!EOP
    TRACE_ENTER("glacier_writerst")
    call glaciermodelwrst(trim(LIS_rc%glaciermodel)//char(0),n)
    TRACE_EXIT("glacier_writerst")
  end subroutine LIS_glaciermodel_writerestart

!BOP
! !ROUTINE: LIS_glaciermodel_finalize
! \label{LIS_glaciermodel_finalize}
!
! !INTERFACE:
  subroutine LIS_glaciermodel_finalize()
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
!  \item[glaciermodelfinalize] (\ref{glaciermodelfinalize}) \newline
!    invokes the generic method in the registry to cleanup the 
!    GLACIERMODEL related datastructures    
! \end{description}
!EOP
    call glaciermodelfinalize(trim(LIS_rc%glaciermodel)//char(0))
  end subroutine LIS_glaciermodel_finalize

end module LIS_glaciermodelMod
