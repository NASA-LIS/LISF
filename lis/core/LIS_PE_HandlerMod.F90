!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_PE_HandlerMod

!BOP
!
! !MODULE: LIS_PE_HandlerMod
!
! !DESCRIPTION:
!  The code in this file implements methods to handle the program flow related
!  to parameter estimation (PE). The module supports parameter estimation
!  of different model types (LSMs, RTMs, Routing, Landslide models) using
!  any of the optimization and uncertainty estimation algorithms in LIS.
!   
! !REVISION HISTORY:
!
!  21 Jun 2009: Sujay Kumar; Initial implementation
!  12 Jan 2012: Sujay Kumar; Implemented an updated version that includes
!                            generic support for different model types.
! 
!
! !NOTES: TBD - There has to be one additional step where we assemble the 
! individual patches from the surface models and send to the optimization 
! algorithm. The OPT/UE algorithm is assumed to work in the tile space. 
!
  use ESMF

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public ::  LIS_PE_init                !initializes the PE handler
  public ::  LIS_readPEobs              !read observations for PE
  public ::  LIS_updatePEObjectiveFunc  !updates the objective function
  public ::  LIS_computePEObjectiveFunc !computes objective function 
  public ::  LIS_setPEdecisionSpace     !sets the decision space
  public ::  LIS_resetPEObjectiveFunc   !resets objective function
  public ::  LIS_resetPEobs             !resets the PE observations
  public ::  LIS_PE_restart             !reads a PE restart file
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_PEOBS_State     !object to store observations for PE
  public :: LIS_PEOBSPred_State !object to store the model simulated values
                                !that correspond to the PE observations
!EOP
  type(ESMF_State), save :: LIS_PEOBS_State
  type(ESMF_State), save :: LIS_PEOBSPred_State

  integer :: nmodel_opts
  logical :: lsm_opt_enabled
  logical :: rtm_opt_enabled
  logical :: landslide_opt_enabled
  logical :: routing_opt_enabled

  integer :: nmodel_obspreds
  logical :: lsm_obspred_enabled
  logical :: rtm_obspred_enabled
  logical :: landslide_obspred_enabled
  logical :: routing_obspred_enabled
  
  real*8  :: stime        ! time to begin writing output
  integer :: syear        ! year to begin writing output
  integer :: smonth       ! month to begin writing output
  integer :: sday         ! day to begin writing output
  integer :: shour        ! hour to begin writing output
  integer :: smin         ! minutes to begin writing output
  integer :: ssec         ! seconds to begin writing output
  integer :: doy
  real    :: gmt

contains


!BOP
! 
! !ROUTINE: LIS_PE_init
! \label{LIS_PE_init}
! 
! !INTERFACE: 
  subroutine LIS_PE_init
! 
! !DESCRIPTION: 
!  This routine allocates the structures required for managing 
!  observational data, the corresponding ``obs predictor (obspred)''
!  from the models and the decision space (list of variables being
!  adjusted) for parameter estimation. 
! 
!  The methods invoked are: 
!  \begin{description}
!  \item[setupPEOBSSpace](\ref{setuppeobsspace}) \newline
!    invokes the setup method for the specified PE 
!    observational source
!  \item[setupPEOBSPredSpace](\ref{setuppeobspredspace}) \newline
!    invokes the setup method for the specified  
!    PE ``obspred'' - the simulated value of observation
!  \item[setupDecSpaceVars](\ref{setupDecSpaceVars}) \newline
!    invokes the methods for setting up the decision
!    space variables \newline
!  \end{description}
!
! !USES:     
    use LIS_coreMod, only : LIS_rc, LIS_config
    use LIS_logMod,  only : LIS_logunit, LIS_verify
    use LIS_timeMgrMod, only : LIS_date2time
!EOP
    integer       :: n 
    integer       :: status
    character*1   :: nestid(2)
    character*1   :: caseid(3)
    character*100 :: temp
    integer              :: max_index
    integer              :: i, k
    integer,    allocatable  :: insts(:)
    character*100        :: type1
    integer :: rc

    call ESMF_ConfigGetAttribute(LIS_config,nmodel_opts,&
         label="Number of model types subject to parameter estimation:",rc=status)
    call LIS_verify(status,&
         'Number of model types subject to parameter estimation: not defined')
    
    call ESMF_ConfigFindLabel(LIS_config,&
         label="Model types subject to parameter estimation:",rc=status)
    call LIS_verify(status,&
         'Model types subject to parameter estimation: not defined')    

    do i=1,nmodel_opts
       call ESMF_ConfigGetAttribute(LIS_config,type1, rc=status)
       if(trim(type1).eq."LSM") then 
          lsm_opt_enabled =.true. 
       elseif(trim(type1).eq."RTM") then 
          rtm_opt_enabled = .true.
       elseif(trim(type1).eq."Routing") then 
          routing_opt_enabled = .true.
       elseif(trim(type1).eq."Landslide") then 
          landslide_opt_enabled = .true. 
       endif
       
       call ESMF_ConfigGetAttribute(LIS_config,nmodel_obspreds,&
            label="Number of model types with observation predictors for parameter estimation:",rc=status)
       call LIS_verify(status,&
            'Number of model types with observation predictors for parameter estimation: not defined')
       
       call ESMF_ConfigFindLabel(LIS_config,&
            label="Model types with observation predictors for parameter estimation:",rc=status)
       call LIS_verify(status,&
            'Model types with observation predictors for parameter estimation: not defined')    
    enddo

    do i=1,nmodel_obspreds
       call ESMF_ConfigGetAttribute(LIS_config,type1, rc=status)
       if(trim(type1).eq."LSM") then 
          lsm_obspred_enabled =.true. 
       elseif(trim(type1).eq."RTM") then 
          rtm_obspred_enabled = .true.
       elseif(trim(type1).eq."Routing") then 
          routing_obspred_enabled = .true.
       elseif(trim(type1).eq."Landslide") then 
          landslide_obspred_enabled = .true. 
       endif
    enddo

! initmode = 1 means initialize with random values, 
! initmode = 0 means initialize with default values only. 
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%decSpaceInitMode,&
         label="Initialize decision space with default values:",rc=status)
    call LIS_verify(status, &
         "Initialize decision space with default values: not defined")

    call ESMF_ConfigGetAttribute(LIS_config,syear,&
         label="Calibration period start year:",default=LIS_rc%syr,rc=rc)
    call LIS_verify(rc,'Calibration period start year: not specified')
    call ESMF_ConfigGetAttribute(LIS_config,smonth,&
         label="Calibration period start month:",default=LIS_rc%smo,rc=rc)
    call LIS_verify(rc,'Calibration period start month: not specified')
    call ESMF_ConfigGetAttribute(LIS_config,sday,&
         label="Calibration period start day:",default=LIS_rc%sda,rc=rc)
    call LIS_verify(rc,'Calibration period start day: not specified')
    call ESMF_ConfigGetAttribute(LIS_config,shour,&
         label="Calibration period start hour:",default=LIS_rc%shr,rc=rc)
    call LIS_verify(rc,'Calibration period start hour: not specified')
    call ESMF_ConfigGetAttribute(LIS_config,smin,&
         label="Calibration period start minutes:",default=LIS_rc%smn,rc=rc)
    call LIS_verify(rc,'Calibration period start minutes: not specified')
    call ESMF_ConfigGetAttribute(LIS_config,ssec,&
         label="Calibration period start seconds:",default=LIS_rc%sss,rc=rc)
    call LIS_verify(rc,'Calibration period start seconds: not specified')
    
    call LIS_date2time(stime, doy, gmt, &
         syear,           &
         smonth,          &
         sday,            &
         shour,           &
         smin,            &
         ssec)
    
    write(unit=LIS_logunit,FMT=*) '[INFO] Calibration period start time:',  &
         stime,  &
         syear,  &
         smonth, &
         sday,   &
         shour,  &
         smin,   &
         ssec
    
    LIS_PEOBS_State = ESMF_StateCreate(name="PE OBS Space",rc=status)
    call LIS_verify(status)

    LIS_PEOBSPred_State = ESMF_StateCreate(name="PE OBS Pred Space",rc=status)
    call LIS_verify(status)

    call setupPEOBSSpace(trim(LIS_rc%optueset)//char(0),LIS_PEOBS_State)
    call setupPEOBSPredSpace()

    call setupDecSpaceVars()

  end subroutine LIS_PE_init

!BOP
! 
! !ROUTINE: setupPEOBSPredSpace
! \label{setuppeobspredspace}
! 
! !INTERFACE: 
  subroutine setupPEOBSPredSpace()
! !USES:
    use LIS_coreMod, only : LIS_rc    
! 
! !DESCRIPTION: 
!  This routine cycles through different model types to 
!  setup the OBSPred object. Depending on which model types
!  are chosen in the optimization, the setup method for 
!  those model types are invoked. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[lsmpesetupobspredspace](\ref{lsmpesetupobspredspace}) \newline
!    sets up the LSM's ObsPred object
!  \item[rtmpesetupobspredspace](\ref{rtmpesetupobspredspace}) \newline
!    sets up the RTM's ObsPred object
!  \end{description}
!EOP    
    implicit none

    if(lsm_obspred_enabled) then
       call lsmpesetupobspredspace(trim(LIS_rc%lsm)//"+"//trim(LIS_rc%optueset)//char(0), &
            LIS_PEOBSPred_State)
    endif
    if(rtm_obspred_enabled) then
       call rtmpesetupobspredspace(trim(LIS_rc%rtm)//"+"//trim(LIS_rc%optueset)//char(0), &
            LIS_PEOBSPred_State)
    endif
  end subroutine setupPEOBSPredSpace
      
!BOP
! 
! !ROUTINE: setupDecSpaceVars
! \label{setupDecSpaceVars}
! 
! !INTERFACE: 
  subroutine setupDecSpaceVars()
! !USES:     
    use LIS_coreMod, only : LIS_rc
    use LIS_optUEMod, only : LIS_decisionSpace,LIS_feasibleSpace
!
! !DESCRIPTION: 
!  This routine cycles through different model types to setup the 
!  decision space object. Depending on which model types are
!  chosen in the optimization, the methods to setup decision 
!  space variables for those model types are invoked. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[lsmpesetupdecisionspace](\ref{lsmpesetupdecisionspace}) \newline
!    sets up the LSM's decision space object
!  \end{description}
!EOP    
    if(lsm_opt_enabled) then 
       call lsmpesetupdecisionspace(trim(LIS_rc%lsm)//char(0), &
            LIS_DecisionSpace, LIS_feasibleSpace)
    endif
    if(rtm_opt_enabled) then 
       call rtmpesetupdecisionspace(trim(LIS_rc%rtm)//char(0), &
            LIS_DecisionSpace, LIS_feasibleSpace)
    endif
  end subroutine setupDecSpaceVars

!BOP
! 
! !ROUTINE: LIS_readPEobs
! \label{LIS_readPEobs}
! 
! !INTERFACE: 
  subroutine LIS_readPEobs(n)
! !USES: 
    use LIS_coreMod,   only : LIS_rc

! !ARGUMENTS: 
    integer, intent(in) :: n

! 
! !DESCRIPTION: 
!  This routines calls the method to read the specific observation
!  source. 
! 
!  The methods invoked are: 
!  \begin{description}
!  \item[getPEObs](\ref{getpeobs}) \newline
!   invokes the method to read the observations and incorporate into 
!   the observation state object
!  \item[writePEObs](\ref{writepeobs}) \newline
!   invokes the method to write the processed PE observations to disk
!  \end{description}
!EOP
    integer          :: k 
    if(LIS_rc%time .ge. stime) then
       call getPEObs(trim(LIS_rc%optueset)//char(0), LIS_PEOBS_State)
       
       if ( LIS_rc%wpeobs .ge. 1 ) then 
          call writePEObs(trim(LIS_rc%optueset)//char(0), LIS_PEOBS_State)
       endif
    end if
  end subroutine LIS_readPEobs

!BOP
! !ROUTINE: LIS_updatePEObjectiveFunc
! \label{LSMPE_LIS_updatePEObjectiveFunc}
! 
! !INTERFACE: 
  subroutine LIS_updatePEObjectiveFunc(n)
! !USES: 
    use LIS_coreMod,   only : LIS_rc
    use LIS_logMod,    only : LIS_verify
    use LIS_optUEMod,  only : LIS_objectiveFunc
! !ARGUMENTS: 
    integer, intent(in)   :: n 
! 
! !DESCRIPTION: 
!  Invokes the call to update the specified objective function
! 
!  The methods invoked are: 
!  \begin{description}
!  \item[getPEObsPred](\ref{getpeobspred}) \newline
!   invokes the method to retrieve the model simulated value of the 
!   PE observations
!  \item[updateObjectiveFunc](\ref{updateobjectivefunctype}) \newline
!   invokes the method to update the objective function based on 
!   the PE observations and the corresponding obspreds. 
!  \end{description}
!EOP
    logical               :: data_status
    integer               :: status

    if(LIS_rc%time .ge. stime) then
       call ESMF_AttributeGet(LIS_PEOBS_State,name="Data Update Status",&
            value=data_status,rc=status)
       call LIS_verify(status)
       
       if(data_status) then 
          call getPEObsPred()
          call updateObjectiveFunc(trim(LIS_rc%objfuncmethod)//char(0))
       endif
    endif
  end subroutine LIS_updatePEObjectiveFunc

!BOP
! 
! !ROUTINE: getPEOBSPred
! \label{getpeobspred}
!
! !INTERFACE: 
  subroutine getPEOBSPred()
! !USES: 
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION: 
!  This routine cycles through different model types to generate the 
!  obspred object, depending on which model types are chosen 
!  in the optimization. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[lsmpegetobspred](\ref{lsmpegetobspred}) \newline
!   invokes the method to retrieve obsPred object from the LSM. 
!  \item[rtmpegetobspred](\ref{rtmpegetobspred}) \newline
!   invokes the method to retrieve obsPred object from the RTM. 
!  \end{description}
!EOP    
       if(lsm_obspred_enabled) & 
            call lsmpegetobspred(trim(LIS_rc%lsm)//"+"//trim(LIS_rc%optueset)//char(0), &
            LIS_PEOBSPred_State)       
       if(rtm_obspred_enabled) & 
            call rtmpegetobspred(trim(LIS_rc%rtm)//"+"//trim(LIS_rc%optueset)//char(0), &
            LIS_PEOBSPred_State)       
    
  end subroutine getPEOBSPred
     
     !BOP
! 
! !ROUTINE: LIS_computePEObjectiveFunc
! \label{LIS_computePEObjectiveFunc}
!
! !INTERFACE: 
  subroutine LIS_computePEObjectiveFunc
! !USES:     
    use LIS_coreMod,    only : LIS_rc
!
! !DESCRIPTION: 
!   This routine calls the method to compute the objective function
!   prior to calling the optimization/uncertainty estimation algorithm
! 
!  The methods invoked are: 
!  \begin{description}
!  \item[computeObjectiveFuncType](\ref{computeobjectivefunctype}) \newline
!   invokes the objective function computation based on the 
!   specified method
!  \end{description}
!EOP
    call computeObjectiveFuncType(trim(LIS_rc%objfuncmethod)//char(0))

  end subroutine LIS_computePEObjectiveFunc

!BOP
! 
! !ROUTINE: LIS_setPEDecisionSpace
! \label{lis_setpedecisionspace}
! 
! !INTERFACE: 
  subroutine LIS_setPEDecisionSpace()
! !USES:     
    use LIS_coreMod,  only : LIS_rc

!
! !DESCRIPTION: 
!  This routine invokes the calls to retrieve the decision space variables
!  from the optimization/uncertainty estimation algorithm and assigns 
!  them to the model structures. 
!  
!  The methods invoked are: 
!  \begin{description}
!  \item[getOptUEAlgnparam](\ref{getoptuealgnparam}) \newline
!   invokes the method to return the number of parameters in the decision space
!   from the optimization/UE algorithm
!  \item[setDecisionSpace](\ref{setdecisionspace}) \newline
!   invokes the method to assign the reconciled decision space 
!   variables to the respective model types.
!  \end{description}
!EOP

    implicit none
    
    integer                 :: n 
    integer                 :: t, i

    n = 1

    !get the decision space values from the algorithm datastructures
    call getOptuealgDecSpace(trim(LIS_rc%optuealg)//char(0), n)   
    call setDecisionSpace()

  end subroutine LIS_setPEDecisionSpace

!BOP
!
! !ROUTINE: setDecisionSpace
! \label{setdecisionspace}
! 
! !INTERFACE:  
  subroutine setDecisionSpace()
! !USES: 
    use LIS_coreMod, only : LIS_rc
    use LIS_optUEMod, only : LIS_decisionSpace, LIS_feasibleSpace
    
!
! !DESCRIPTION: 
!  This routine cycles through different model types to set the 
!  decision space variables back to the respective models.
!
!  The methods invoked are: 
!  \begin{description}
!  \item[lsmpesetdecisionspace](\ref{lsmpesetdecisionspace}) \newline
!   invokes the method to return decision space object from the 
!   LSM. 
!  \end{description}
!EOP
    implicit none

    integer               :: status

    if(lsm_opt_enabled) then 
       call lsmpesetdecisionspace(trim(LIS_rc%lsm)//char(0),  &
            LIS_decisionSpace, LIS_feasibleSpace)
    endif
    if(rtm_opt_enabled) then 
       call rtmpesetdecisionspace(trim(LIS_rc%rtm)//char(0),  &
            LIS_decisionSpace, LIS_feasibleSpace)
    endif
      
  end subroutine setDecisionSpace

!BOP
! 
! !ROUTINE: LIS_resetPEobjectiveFunc
! \label{LIS_resetPEobjectiveFunc}
!
! !INTERFACE: 
  subroutine LIS_resetPEobjectiveFunc
! !USES:     
    use LIS_coreMod,     only : LIS_rc
!
! !DESCRIPTION:
!   This routine calls the method to reset the objective function
!   for the next iteration of the optimization/uncertainty estimation
!   algorithm. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[resetObjectiveFuncType](\ref{resetobjectivefunctype}) \newline
!   invokes the method to reset the objective function objects
!  \end{description}
!EOP 
    call resetObjectiveFuncType(trim(LIS_rc%objfuncmethod)//char(0))

  end subroutine LIS_resetPEobjectiveFunc

!BOP
! 
! !ROUTINE: LIS_resetPEobs
! \label{LIS_resetPEobs}
!
! !INTERFACE: 
  subroutine LIS_resetPEobs
! !USES:     
    use LIS_coreMod,     only : LIS_rc
!
! !DESCRIPTION:
!   This routine calls the method to reset the observation data structures
!   for the next iteration of the optimization/uncertainty estimation
!   algorithm. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[resetpeobsspace](\ref{resetpeobsspace}) \newline
!   invokes the method to reset the obs state objects
!  \end{description}
!EOP 
    call resetpeobsspace(trim(LIS_rc%optueset)//char(0),LIS_PEOBS_State)

  end subroutine LIS_resetPEobs

!BOP
! 
! !ROUTINE: LIS_PE_restart
! \label{LIS_PE_restart}
!
! !INTERFACE: 
  subroutine LIS_PE_restart
! !USES:     
    use LIS_optUEMod, only : LIS_optUEAlg_readrestart
!
! !DESCRIPTION: 
!  This routine invokes the method to read the restart file for
!  the optimization algorithm. If a restart file is read, then 
!  the call to update the model's decision space variables is 
!  also invoked. 
!
!  The methods invoked are: 
!  \begin{description}
!  \item[LIS\_optUEAlg\_readrestart](\ref{LIS_optUEAlg_readrestart}) \newline
!   invokes the method to read the restart file of the OPT/UE algorithm. 
!  \item[LIS\_setPEdecisionSpace](\ref{lis_setpedecisionspace})
!   calls the method to set the decision space variables in the 
!   model. 
!  \end{description}
!EOP

    implicit none
!    logical          :: rstflag

    call LIS_optUEAlg_readrestart

!    if(rstflag) then 
!       call LIS_setPEdecisionSpace 
!    endif

  end subroutine LIS_PE_restart

end module LIS_PE_HandlerMod
