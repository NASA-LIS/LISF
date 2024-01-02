!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_optUEMod

!BOP
!
! !MODULE: LIS_optUEMod
!
! !DESCRIPTION:
!  The code in this file provides interfaces to manage optimization and
!  uncertainty modeling implementations in LIS 
! 
!  \subsubsection{Overview}
!  The module provides interfaces for incorporating optimization and uncertainty
!  modeling implementations into LIS. The abstract representation is to have the 
!  algorithm work between a certain objective function and a certain decision 
!  space using a specified opt/UE algorithm. The algorithm attempts to find
!  potential solutions by optimizating the objective 
!  function with respect to the decision space constraints. 
!  
!  An example of an optimization algorithm is
!  the suite of Genetic Algorithm-based heuristic optimization techniques
!  
! !REVISION HISTORY:
!
!  15 May 2007: Sujay Kumar; Initial implementation
!
!EOP
  use ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_optUE_init    ! initialization for perturbation routines
  public :: LIS_objectiveFunc_init
  public :: LIS_optUEAlg_init
  public :: LIS_isOptStopCriterionTrue ! method to check if the stopping criteria
                                      ! for the algorithm is met 
  public :: LIS_runOptUEAlg  ! run the optimization/uncertainty estimation algorithm
  public :: LIS_OptUEAlg_reset
  public :: LIS_optUEAlg_readrestart
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LIS_DecisionSpace   ! LIS object to store the decision space
  public :: LIS_FeasibleSpace   ! LIS object to store feasibility checks
  public :: LIS_ObjectiveFunc   ! LIS object to store the objective function
!EOP

  type(ESMF_State), save  :: LIS_DecisionSpace
  type(ESMF_State), save  :: LIS_FeasibleSpace
  type(ESMF_State), save  :: LIS_ObjectiveFunc
!EOP


contains
!BOP
! 
! !ROUTINE: LIS_optUE_init
! \label{LIS_optUE_init}
! 
! !INTERFACE:  
  subroutine LIS_optUE_init
! !USES: 
    use LIS_coreMod,  only : LIS_rc, LIS_config, LIS_vecTile
    use LIS_logMod,   only : LIS_verify, LIS_logunit
    
! 
! !DESCRIPTION: 
! 
!  This routine intializes the optimization/uncertainty estimation
!  mode, the objective function and decision space objects. 
!
! 
!EOP
    integer :: status
    type(ESMF_Field)            :: varField
    type(ESMF_ArraySpec)        :: arrspec1
    integer, pointer            :: mod_flag(:)
    integer                     :: n 

    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%optUEAlg,&
         label="Optimization/Uncertainty Estimation Algorithm:",rc=status)
    call LIS_verify(status,'Optimization/Uncertainty Estimation Algorithm: not defined')

!    LIS_rc%optUEType = 0 
!    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%optUEtype,&
!         label="Optimization/Uncertainty Estimation Type:",rc=status)
!    call LIS_verify(status,'Optimization/Uncertainty Estimation Type: not defined')

    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%optUEset,&
         label="Optimization/Uncertainty Estimation Set:",rc=status)
    call LIS_verify(status,'Optimization/Uncertainty Estimation Set: not defined')

    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%objfuncmethod,&
         label="Objective Function Method:",rc=status)
    call LIS_verify(status,'Objective Function Method: not defined')
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%wpeobs,&
         label="Write PE Observations:",rc=status)
    call LIS_verify(status,'Write PE Observations: not defined')

    LIS_decisionSpace = ESMF_StateCreate(name="Decision Space",rc=status)
    call LIS_verify(status)


    LIS_FeasibleSpace = ESMF_StateCreate(name="Feasible Space",rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(arrspec1, rank=1,typekind = ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)
    n = 1

    varField = ESMF_FieldCreate(arrayspec=arrspec1, grid=LIS_vecTile(n),&
         name = "Feasibility Flag",rc=status)
    call LIS_verify(status)
    call ESMF_StateAdd(LIS_feasibleSpace, (/varField/),rc=status)
    call LIS_verify(status)

!initialize the feasibility flag
    call ESMF_FieldGet(varField,localDE=0,farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)
    mod_flag = 0
!    write(LIS_logunit,*) 'LIS_optUE_init mod_flag:', mod_flag

    LIS_objectiveFunc = ESMF_StateCreate(name="Objective Function",rc=status)
    call LIS_verify(status)
    
  end subroutine LIS_optUE_init

!BOP
! 
! !ROUTINE: LIS_objectiveFunc_init
! \label{LIS_objectiveFunc_init}
! 
! !INTERFACE:  
  subroutine LIS_objectiveFunc_init
! !USES: 
    use LIS_coreMod,  only : LIS_rc
    
! 
! !DESCRIPTION: 
! 
!  This routine intializes the optimization/uncertainty estimation
!  algorithm, the objective  function 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[objectivefunctypeinit](\ref{objectivefunctypeinit}) \newline
!    invokes the initialization routine for the specified objection
!    function metric/method 
!  \end{description}
! 
!EOP
    call objectivefunctypeinit(trim(LIS_rc%objfuncmethod)//char(0))

  end subroutine LIS_objectiveFunc_init

!BOP
! 
! !ROUTINE: LIS_optUEAlg_init
! \label{LIS_optUEAlg_init}
! 
! !INTERFACE:  
  subroutine LIS_optUEAlg_init
! !USES: 
    use LIS_coreMod,  only : LIS_rc
    
! 
! !DESCRIPTION: 
! 
!  This routine intializes the optimization/uncertainty estimation
!  algorithm. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[optUEalginit](\ref{optuealginit}) \newline
!    invokes the init routine for the specified optimization
!    or uncertainty estimation algorithm.  
!  \end{description}
! 
!EOP
    call optuealginit(trim(LIS_rc%optUEAlg)//char(0))

  end subroutine LIS_optUEAlg_init

!BOP
! 
! !ROUTINE: LIS_isOptStopCriterionTrue
! \label{LIS_isOptStopCriterionTrue}
!
! !INTERFACE:  
  function LIS_isOptStopCriterionTrue() result(finish)
! !USES: 
    use LIS_coreMod,   only : LIS_rc
    use LIS_logMod,          only : LIS_logunit
! !ARGUMENTS: 
    logical :: finish
! 
! !DESCRIPTION: 
! Invokes the appropriate method from the registry to check 
! if the convergence criteria for the specified optimization 
! /uncertainty estimation algorithm is met. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[checkconvergence](\ref{checkconvergence}) \newline
!    invokes the method from the registry that specifies the 
!    stopping criteria of the algorithm
!  \end{description}
! 
!EOP        

!    write(LIS_logunit,*) 'begin checkconvergence in LIS_isOptStopCriterionTrue'
    call checkconvergence(trim(LIS_rc%optuealg)//char(0),finish)
!    write(LIS_logunit,*) 'end checkconvergence in LIS_isOptStopCriterionTrue:', finish

  end function LIS_isOptStopCriterionTrue

!BOP
! !ROUTINE: LIS_runoptUE
! \label{LIS_runoptUE}
! 
! !INTERFACE: 
  subroutine LIS_runOptUEAlg()
! !USES: 
    use LIS_coreMod,   only : LIS_rc
! 
! !DESCRIPTION: 
! 
!  invokes the run method of the selected optimization/uncertainty 
!  estimation algorithm from the registry 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[runoptue](\ref{runoptue}) \newline
!    invokes the specific optue algorithm related operations
!    to solve for potential solutions
!  \item[resetobjectivefunctype](\ref{resetobjectivefunctype}) \newline
!    invokes the call to reset the objective function objects based
!    on the specific method
!  \end{description}
!EOP
    call runoptue(trim(LIS_rc%optuealg)//char(0))
    
  end subroutine LIS_runOptUEAlg

!BOP
! !ROUTINE: LIS_optUEAlg_readrestart
! \label{LIS_optUEAlg_readrestart}
! 
! !INTERFACE: 
  subroutine LIS_optUEAlg_readrestart()
! !USES: 
    use LIS_coreMod, only : LIS_rc
! 
! !DESCRIPTION: 
!  invokes the read restart method for the selected optimization/uncertainty
!  estimation algorithm from the registry
! 
!  The methods invoked are: 
!  \begin{description}
!   \item[optuereadrestart](\ref{optuereadrestart})
!    invokes the call to read the restart file from the OPT/UE algorithm.
!  \end{description}
!EOP    
    implicit none 
    
    logical         :: rstflag

    call optuereadrestart(trim(LIS_rc%optuealg)//char(0), rstflag)

  end subroutine LIS_optUEAlg_readrestart


!BOP
! 
! !ROUTINE: LIS_optUEAlg_reset
! \label{LIS_optUEAlg_reset}
! 
! !INTERFACE:  
  subroutine LIS_optUEAlg_reset
! !USES: 
    use LIS_logMod
    
! 
! !DESCRIPTION: 
! 
! 
!EOP
    type(ESMF_Field)  :: feasField
    integer           :: status
    integer, pointer  :: modflag(:)

    call ESMF_StateGet(LIS_FeasibleSpace, "Feasibility Flag", feasField, rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField,localDE=0,farrayPtr=modflag,rc=status)
    call LIS_verify(status)

    modflag = 0
  end subroutine LIS_optUEAlg_reset
end module LIS_optUEMod
