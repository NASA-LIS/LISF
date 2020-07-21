!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
!  the LSM state, prompts the land surface model to transform 
!  the OBS state to an LSM state. The transformed LSM state is then 
!  used to overwrite the current LSM state. 
! 
!  The methods invoked are: 
!  \begin{description}
!   \item[lsmdagetstatevar](\ref{lsmdagetstatevar}) \newline
!    obtain the specified LSM state prognostic variables
!   \item[Lsmdaobstransform](\ref{lsmdaobstransform}) \newline
!    transform the observation state into an LSM space
!   \item[lsmdamapobstolsm](\ref{lsmdamapobstolsm}) \newline
!    map the observation state to generate the analysis increments
!  \end{description} 
!
! !USES:     
    use ESMF
    use LIS_coreMod,only          : LIS_rc
    use LIS_lsmMod, only          : LIS_LSM_State, LIS_LSM_Incr_State
    use LIS_DAobservationsMod, only : LIS_OBS_State
    use LIS_logMod, only          : LIS_logunit, LIS_verify
!EOP
    integer, intent(in)      :: n 
    integer, intent(in)      :: k
    logical                  :: data_status
    integer                  :: status

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Update Status",&
         value=data_status,rc=status)
    call LIS_verify(status)

    call ESMF_AttributeSet(LIS_LSM_Incr_State(n,k),&
         "Fresh Increments Status", .false., rc=status)
    call LIS_verify(status)    
    
    if(data_status) then 
       write(LIS_logunit,*) 'Using Direct Insertion for Assimilation',n
       call lsmdagetstatevar(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k))
       call Lsmdaobstransform(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_OBS_State(n,k))
       call lsmdamapobstolsm(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, k, LIS_OBS_State(n,k),&
            LIS_LSM_Incr_State(n,k))
       
       call ESMF_AttributeSet(LIS_LSM_Incr_State(n,k), &
            "Fresh Increments Status",&
            value = .true., rc=status)
    endif

  end subroutine di_increments


!BOP
! !ROUTINE: di_update
! \label{di_update}
! 
! !INTERFACE:
  subroutine di_update(n,k)
! !USES: 
    use ESMF
    use LIS_coreMod, only         : LIS_rc
    use LIS_lsmMod, only          : LIS_LSM_State, LIS_LSM_Incr_State
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
!   \item[lsmdasetstatevar](\ref{lsmdasetstatevar}) \newline
!    set the specified LSM state prognostic variables
!   \end{description}
!EOP


    logical       :: fresh_incr
    integer                  :: status

    call ESMF_AttributeGet(LIS_LSM_Incr_State(n,k),"Fresh Increments Status",&
         value = fresh_incr, rc=status)

    if(fresh_incr) then 
!-----------------------------------------------------------------
! apply increments
!-----------------------------------------------------------------
       call lsmdaupdatestate(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k),&
            LIS_LSM_Incr_State(n,k))       
!----------------------------------------------------------------------
!  Update the LSM's state variables
!----------------------------------------------------------------------       
       call lsmdaqcstate(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0),n,LIS_LSM_State(n,k))
       call lsmdasetstatevar(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0),n, LIS_LSM_State(n,k))
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
