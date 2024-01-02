!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module directInsertion_Mod
!BOP
!
! !MODULE: directInsertion_Mod
! \label{directInsertion_Mod}
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines that control
!   the incorporation of a data set using direct insertion in a 
!   land surface model. 
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
!  21Jun06    Sujay Kumar:   Updated implementation with ESMF 
!                            structures
! 
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_DAobservationsMod
  use LIS_surfaceModelMod
  use LIS_logMod

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: di_init  ! Initialization for Direct Insertion  
  public :: di_setup
  public :: di_increments ! Routine to compute DI increments
  public :: di_update     ! Routine to apply DI increments
  public :: di_diagnostics ! write DI related diagnostics 
  public :: di_final ! Finalization for Direct Insertion    
!EOP

contains
!BOP
! 
! !ROUTINE: di_init
! \label{di_init}
!  
! !INTERFACE: 
  subroutine di_init
! 
! !DESCRIPTION: 
!  This method performs the required initializations for the direct
!  insertion method. 
!
!EOP

  end subroutine di_init

  subroutine di_setup(k)
    
    integer, intent(in) :: k 

  end subroutine di_setup


!BOP
! 
! !ROUTINE: di_increments
! \label{di_increments}
!  
! !INTERFACE: 
  subroutine di_increments(n,k)
! 
! !DESCRIPTION: 
!  This method computes the analysis increments from direct insertion.
!  The direct insertion method simply retrieves 
!  the state vector, prompts the model to transform 
!  the OBS state to a model state vector. The transformed state is then 
!  used to overwrite the current model state. 
! 
!  The methods invoked are: 
!  \begin{description}
!   \item[LIS\_surfaceModel_DAgetstatevar]\ref{LIS_surfaceModel_DAgetstatevar}
!    obtain the specified prognostic variables. 
!   \item[LIS\_surfaceModel_DAobsTransform](\ref{LIS_surfacemodel_DAobsTransform}) \newline
!    transform the observation state into a model space
!   \item[LIS\_surfaceModel\_DAmapobsToModel](\ref{LIS_surfaceModel_DAmapObsToModel}) \newline
!    map the observation state to generate the analysis increments
!  \end{description} 
!
! !USES:     

!EOP
    integer, intent(in)      :: n 
    integer, intent(in)      :: k
    logical                  :: data_status
    integer                  :: status

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Update Status",&
         value=data_status,rc=status)
    call LIS_verify(status)

    call LIS_surfaceModel_DASetFreshIncrementsStatus(n,k,.false.)
    
    if(data_status) then 
       write(LIS_logunit,*) '[INFO] Using Direct Insertion for Assimilation',k
       call LIS_surfaceModel_DAGetStateVar(n,k)

       call LIS_surfaceModel_DAobsTransform(n,k)
       call LIS_surfaceModel_DAmapObsToModel(n,k)

       call LIS_surfaceModel_DASetFreshIncrementsStatus(n,k,.true.)
    endif

  end subroutine di_increments


!BOP
! !ROUTINE: di_update
! \label{di_update}
! 
! !INTERFACE:
  subroutine di_update(n,k)
! !USES: 
    use LIS_coreMod, only         : LIS_rc
! !ARGUMENTS: 
    integer, intent(in)    :: n
    integer, intent(in)    :: k

! !DESCRIPTION: 
!  This method applies the analysis increments computed by the 
!  direct insertion method and updates the model prognostic 
!  variables. 
! 
!  The methods invoked are: 
!  \begin{description}
!   \item[LIS_surfaceModel\_DAqcstate]\ref{LIS_surfaceModel_DAqcstate}
!    QC the updated state state
!   \item[LIS\_surfaceModel_DAsetstatevar](\ref{LIS_surfaceModel_DAsetstatevar}) \newline
!    assigns the specified state prognostic variablse
!   \end{description}
!EOP


    logical       :: fresh_incr
    integer                  :: status

    call LIS_surfaceModel_DAGetFreshIncrementsStatus(n,k,fresh_incr)

    if(fresh_incr) then 
!-----------------------------------------------------------------
! apply increments
!-----------------------------------------------------------------
       call LIS_surfaceModel_DAUpdateState(n,k)
!----------------------------------------------------------------------
!  Update the state variables
!----------------------------------------------------------------------       
       call LIS_surfaceModel_DAQCstate(n,k)
       
       call LIS_surfaceModel_DASetStateVar(n,k)

    endif

  end subroutine di_update
!BOP
! 
! !ROUTINE: di_diagnostics
! \label{di_diagnostics}
! 
! !INTERFACE:
  subroutine di_diagnostics(n,k)
! 
! !DESCRIPTION: 
!  This subroutine generates the direct insertion related diagnostics 
!
!EOP
    integer, intent(IN)    :: n 
    integer, intent(IN)    :: k
  end subroutine di_diagnostics

!BOP
! 
! !ROUTINE: di_final
! \label{di_final}
!  
! !INTERFACE: 
  subroutine di_final
! 
! !DESCRIPTION: 
!  This method performs the finalization for all direct insertion
!  related structures and subroutines. 
!
!EOP

  end subroutine di_final

end module directInsertion_Mod
