!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module paramEstim_runMod
!BOP
!
! !MODULE: paramEstim_runMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   LIS initialization, execution, and finalization
!   for parameter estimation instances. 
!   
! !REVISION HISTORY: 
!  21Oct05    Sujay Kumar  Initial Specification
! 
!
  implicit none
  
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: lis_init_paramEstim  !init method for optimization-uncertainty modeling mode
  public :: lis_run_paramEstim   !run method for optimization-uncertainty modeling mode
  public :: lis_final_paramEstim !finalize method for optimization-uncertainty modeling mode

!EOP  
contains
!BOP
! !ROUTINE: lis_init_paramEstim
!
! 
! !INTERFACE:
  subroutine lis_init_paramEstim
! !USES: 
    use LIS_coreMod,         only : LIS_core_init
    use LIS_domainMod,       only : LIS_domain_init
    use LIS_surfaceModelMod, only : LIS_surfaceModel_init, &
         LIS_surfaceModel_setup, LIS_surfaceModel_readrestart
    use LIS_metforcingMod,  only : LIS_metforcing_init
    use LIS_perturbMod,     only : LIS_perturb_init, LIS_perturb_readrestart
    use LIS_DAobservationsMod, only : LIS_initDAObservations
    use LIS_dataAssimMod,   only : LIS_dataassim_init
    use LIS_paramsMod,      only : LIS_param_init
    use LIS_routingMod,     only : LIS_routing_init, LIS_routing_readrestart
    use LIS_RTMMod,          only : LIS_RTM_init
    use LIS_appMod,          only : LIS_appModel_init
    use LIS_optUEMod, only : LIS_optUE_init, &
         LIS_objectiveFunc_init,LIS_optUEAlg_init
    use LIS_PE_HandlerMod,    only : LIS_PE_init, LIS_PE_restart, &
         LIS_setPEDecisionSpace

! !DESCRIPTION:
!  This is the initialize method for LIS in a optUE running mode. 
!  The following calls are invoked from this method. 
! \begin{description} 
!  \item[LIS\_domain\_init] (\ref{LIS_domain_init}) \newline
!    initialize the LIS domains
!  \item[LIS\_param\_init] (\ref{LIS_param_init}) \newline
!    initialize parameters
!  \item[LIS\_lsm\_init] (\ref{LIS_lsm_init}) \newline
!    initialize the land surface model. 
!  \item[LIS\_metforcing\_init] (\ref{LIS_metforcing_init}) \newline
!    initialize the met forcing
!  \item[LIS\_setuplsm] (\ref{LIS_setuplsm}) \newline
!    complete the LSM setups
!  \item[LIS\_lsm\_readrestart] (\ref{LIS_lsm_readrestart}) \newline
!    read the restart files 
! \end{description} 
!
!EOP

    call LIS_domain_init
    call LIS_param_init
    call LIS_perturb_init
    call LIS_surfaceModel_init
    call LIS_metforcing_init
    call LIS_initDAObservations
    call LIS_dataassim_init
    call LIS_surfaceModel_setup
    call LIS_routing_init
    call LIS_routing_readrestart
    call LIS_RTM_init
    call LIS_appModel_init
    call LIS_optUE_init
    call LIS_PE_init            ! switch because P obj fxn needs number of dec var's ...
    call LIS_objectiveFunc_init ! ... to perform check that prior is placed over all dec var's
    call LIS_optUEAlg_init 
    call LIS_PE_restart
    call LIS_surfaceModel_readrestart  !moved below for opt of init conditions
    call LIS_setPEDecisionSpace
    call LIS_perturb_readrestart
    call LIS_core_init   

  end subroutine lis_init_paramEstim

!BOP
! !ROUTINE: lis_run_paramEstim
!
! !INTERFACE:
  subroutine lis_run_paramEstim
! !USES:
    use LIS_coreMod,  only : LIS_rc, LIS_endofrun, LIS_timetoRunNest, &
         LIS_ticktime, LIS_core_init
    use LIS_surfaceModelMod, only : LIS_surfaceModel_readrestart, LIS_surfaceModel_init, &
         LIS_surfaceModel_setup, LIS_surfaceModel_finalize, LIS_surfaceModel_reset
    use LIS_metforcingMod,  only :  LIS_metforcing_init, LIS_metforcing_reset, LIS_metforcing_finalize 
    use LIS_perturbMod,     only : LIS_perturb_init, LIS_perturb_readrestart
    use LIS_DAobservationsMod, only : LIS_initDAObservations
    use LIS_dataAssimMod,   only : LIS_dataassim_init
    use LIS_routingMod,     only : LIS_routing_init, LIS_routing_readrestart
    use LIS_RTMMod,          only : LIS_RTM_init
    use LIS_appMod,          only : LIS_appModel_init
    use LIS_paramsMod,       only : LIS_param_reset, LIS_param_init, LIS_param_finalize
    use LIS_timeMgrMod,      only : LIS_resetClock, LIS_timemgr_init
    use LIS_logMod,          only : LIS_logunit
    use LIS_optUEMod
    use LIS_PE_HandlerMod,   only : LIS_readPEobs, LIS_computePEobjectiveFunc,&
         LIS_resetPEObjectiveFunc, LIS_resetPEobs,LIS_updatePEObjectiveFunc,&
         LIS_setPEDecisionSpace
    use LIS_logMod,          only : LIS_logunit
!
! !DESCRIPTION:
! 
!  This is the run method for LIS in a optUE running mode. 
!  The following calls are invoked from this method. 
!
!EOP
    integer           :: n 

    do while (.NOT. LIS_isOptStopCriterionTrue())
       call LIS_resetPEObjectiveFunc

       do while(.not.LIS_endofrun())       
          call LIS_ticktime
          do n=1,LIS_rc%nnest
             if(LIS_timeToRunNest(n)) then
                call LIS_run_step(n)
                call LIS_readPEobs(n)
                call LIS_updatePEobjectivefunc(n)
             endif
          enddo
       enddo
       
       call LIS_computePEobjectiveFunc
       call LIS_resetClock(LIS_rc)
       call LIS_param_reset
       call LIS_metforcing_reset
       call LIS_surfaceModel_reset
       call LIS_runOptUEAlg()
       call LIS_OptUEAlg_reset()
       call LIS_setPEDecisionSpace()
       call LIS_resetPEobs
       call LIS_surfaceModel_readrestart
       call LIS_perturb_readrestart
       call LIS_core_init   
    enddo

  end subroutine lis_run_paramEstim

  subroutine lis_run_step(n)

    use LIS_coreMod,         only : LIS_rc, LIS_endofrun, LIS_timetoRunNest, &
         LIS_ticktime
    use LIS_surfaceModelMod, only : LIS_surfaceModel_f2t, LIS_surfaceModel_run,&
         LIS_surfaceModel_output, LIS_surfaceModel_writerestart, &
         LIS_surfaceModel_perturb_states
    use LIS_paramsMod,        only : LIS_setDynparams
    use LIS_metforcingMod,  only : LIS_get_met_forcing, LIS_perturb_forcing
    use LIS_perturbMod,      only : LIS_perturb_writerestart
    use LIS_DAobservationsMod, only : LIS_readDAobservations, &
         LIS_perturb_DAobservations
    use LIS_dataAssimMod,    only : LIS_dataassim_run, LIS_dataassim_output
    use LIS_routingMod,      only : LIS_routing_run, LIS_routing_writeoutput, &
         LIS_routing_writerestart
    use LIS_RTMMod,          only : LIS_RTM_run,LIS_RTM_output
    use LIS_logMod,          only : LIS_logunit

    integer, intent(in) :: n

     call LIS_setDynparams(n)
     call LIS_get_met_forcing(n)
     call LIS_perturb_forcing(n)
     call LIS_surfaceModel_f2t(n)  

     call LIS_surfaceModel_run(n)

     call LIS_routing_run(n)
     call LIS_routing_writeoutput(n)
     call LIS_routing_writerestart(n)

     call LIS_RTM_run(n)
     call LIS_RTM_output(n)

     call LIS_surfaceModel_perturb_states(n)
     call LIS_readDAobservations(n)
     call LIS_perturb_DAobservations(n)   
     call LIS_perturb_writerestart(n)
     call LIS_dataassim_run(n)
     call LIS_dataassim_output(n)

     call LIS_surfaceModel_output(n)  
     call LIS_surfaceModel_writerestart(n)

  end subroutine lis_run_step

!BOP
! !ROUTINE: lis_final_paramEstim
!
! !INTERFACE:
  subroutine lis_final_paramEstim
! !USES:
    use LIS_coreMod,       only : lis_finalize
    use LIS_surfaceModelMod,          only : LIS_surfaceModel_finalize
    use LIS_paramsMod,        only : LIS_param_finalize
    use LIS_metforcingMod,  only : LIS_metforcing_finalize

! !DESCRIPTION:
! 
!  This is the finalize method for LIS in a optUE running mode. 
!  The following calls are invoked from this method. 
! \begin{description} 
!  \item[LIS\_finalize] (\ref{LIS_finalize}) \newline
!    cleanup LIS generic structures
!  \item[LIS\_surfaceModel\_finalize] (\ref{LIS_surfaceModel_finalize}) \newline
!    cleanup land surface model specific structures
!  \item[LIS\_param\_finalize] (\ref{LIS_param_finalize}) \newline
!    cleanup parameter specific structures
!  \item[LIS\_metforcing\_finalize] (\ref{LIS_metforcing_finalize}) \newline
!    cleanup metforcing specific structures
!  \item[LIS\_dataassim\_finalize] (\ref{LIS_dataassim_finalize}) \newline
!    cleanup data assimilation specific structures
! \end{description} 
!
!EOP

    call lis_finalize()
    call LIS_surfaceModel_finalize()
    call LIS_param_finalize()
    call LIS_metforcing_finalize()
!    call LIS_dataassim_finalize()

  end subroutine lis_final_paramEstim
end module paramEstim_runMod
