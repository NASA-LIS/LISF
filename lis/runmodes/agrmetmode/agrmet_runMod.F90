!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module agrmet_runMod
!BOP
!
! !MODULE: agrmet_runMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   LIS initialization, execution, and finalization
!   for the forecast runmode used by AGRMET in LIS.
!   
! !REVISION HISTORY: 
!  21Oct05    Sujay Kumar  Initial Specification
!  08Nov07    Yudong Tian  Simplified precipitation processing logic. 
!			   No adjustment of cycle time here. It just marches along, 
!			   and metforcing controls when to do precip processing. 
!  25 Apr 2023 Eric Kemp  Add routing support.
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: lis_init_agrmet  !init method for AFWA forecast mode
  public :: lis_run_agrmet   !run method for AFWA forecast mode
  public :: lis_final_agrmet !finalize method for AFWA forecast mode

!EOP  
contains
!BOP
! !ROUTINE: lis_init_agrmet
!
! 
! !INTERFACE:
  subroutine lis_init_agrmet
! !USES: 
    use LIS_coreMod
    use LIS_domainMod,         only : LIS_domain_init
    use LIS_surfaceModelMod,   only : LIS_surfaceModel_init,  &
                                      LIS_surfaceModel_setup, &
                                      LIS_surfaceModel_readrestart
    use LIS_metforcingMod,     only : LIS_metforcing_init
    use LIS_DAobservationsMod, only : LIS_initDAObservations
    use LIS_perturbMod,        only : LIS_perturb_init, LIS_perturb_readrestart
    use LIS_dataAssimMod,      only : LIS_dataassim_init
    use LIS_paramsMod,         only : LIS_param_init
    use LIS_routingMod,        only : LIS_routing_init, &
         LIS_routing_readrestart ! EMK

! !DESCRIPTION:
!  This is the initialize method for LIS in the AGRMET running mode. 
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
!  \item[LIS\_initDAObservations] (\ref{LIS_initDAobservations}) \newline
!    initialize structures needed to read observations for 
!    data assimilation
!  \item[LIS\_setuplsm] (\ref{LIS_setuplsm}) \newline
!    complete the LSM setups
!  \item[LIS\_lsm\_readrestart] (\ref{LIS_lsm_readrestart}) \newline
!    read the restart files 
! \end{description} 
!
!EOP

!YDT,11/8/07:  when precip is not ready, agrmet_mod%pcp_ready = .false.  
! precip processing is triggered by: 
!  pcp_ready = .false. or the first time step past 0Z or 12Z 
!  once triggered, the processing will finish the whole 
!  (0Z, 12Z] or (12Z, 24Z] segment, depending on the current time. 

    LIS_rc%run_model = .true.      !always march the model

    call LIS_domain_init
    call LIS_param_init
    call LIS_perturb_init
    call LIS_surfaceModel_init
    call LIS_metforcing_init
    call LIS_initDAObservations
    call LIS_routing_init ! EMK
    call LIS_routing_readrestart ! EMK
    call LIS_dataassim_init
    call LIS_surfaceModel_setup
    call LIS_surfaceModel_readrestart
    call LIS_perturb_readrestart
    call LIS_core_init

  end subroutine lis_init_agrmet

!BOP
! !ROUTINE: lis_run_agrmet
!
! !INTERFACE:
  subroutine lis_run_agrmet
! !USES:
    use LIS_coreMod,           only : LIS_rc, LIS_endofrun, LIS_timeToRunNest, &
                                      LIS_ticktime
    use LIS_surfaceModelMod,   only : LIS_surfaceModel_f2t, LIS_surfaceModel_run,&
                                      LIS_surfaceModel_output, LIS_surfaceModel_writerestart, &
                                      LIS_surfaceModel_perturb_states
    use LIS_paramsMod,         only : LIS_setDynparams
    use LIS_metforcingMod,     only : LIS_get_met_forcing, LIS_perturb_forcing
    use LIS_perturbMod,        only : LIS_perturb_writerestart
    use LIS_DAobservationsMod, only : LIS_readDAobservations, &
                                      LIS_perturb_DAobservations
    use LIS_dataAssimMod,      only : LIS_dataassim_run, LIS_dataassim_output
    use LIS_logMod,            only : LIS_logunit
    use LIS_routingMod,        only : LIS_routing_run, &
         LIS_routing_writeoutput, LIS_routing_writerestart ! EMK
!
! !DESCRIPTION:
! 
!  This is the run method for LIS in the AGRMET running mode. 
!  The following calls are invoked from this method. 
! \begin{description} 
!  \item[LIS\_endofrun] (\ref{LIS_endofrun}) \newline
!    check to see if the end of simulation has reached
!  \item[LIS\_ticktime] (\ref{LIS_ticktime}) \newline
!    advance model clock
!  \item[LIS\_timeToRunNest] (\ref{LIS_timeToRunNest}) \newline
!    check to see if the current nest needs to be run. 
!  \item[LIS\_setDynparams] (\ref{LIS_setDynparams}) \newline
!    set the time dependent parameters
!  \item[LIS\_get\_met\_forcing] (\ref{LIS_get_met_forcing}) \newline
!    retrieve the met forcing
!  \item[LIS\_perturb\_forcing] (\ref{LIS_perturb_forcing}) \newline
!    perturbs the met forcing
!  \item[LIS\_surfaceModel\_f2t] (\ref{LIS_surfaceModel_f2t}) \newline
!    transfer forcing to model tiles
!  \item[LIS\_surfaceModel\_run] (\ref{LIS_surfaceModel_run}) \newline
!    run the land surface model
!  \item[LIS\_surfaceModel\_perturb\_states] (\ref{LIS_surfaceModel_perturb_states}) \newline
!    perturb the land surface model states
!  \item[LIS\_readDAobservations] (\ref{LIS_readDAobservations}) \newline
!    read observations to be used for data assimilation
!  \item[LIS\_perturb\_DAobservations] (\ref{LIS_perturb_DAobservations}) \newline
!    perturb observations to be used for data assimilation
!  \item[LIS\_dataassim\_run] (\ref{LIS_dataassim_run}) \newline
!    run the data assimilation algorithm
!  \item[LIS\_dataassim\_output] (\ref{LIS_dataassim_output}) \newline
!    write da output
!  \item[LIS\_surfaceModel\_output] (\ref{LIS_surfaceModel_output}) \newline
!    write surface model output
!  \item[LIS\_surfaceModel\_writerestart] (\ref{LIS_surfaceModel_writerestart}) \newline
!    write surface model restart files
! \end{description} 
!
!EOP
    integer :: n,ierr

    do while (.NOT. LIS_endofrun())

       do n=1,LIS_rc%nnest
          if(LIS_timeToRunNest(n)) then 
             call LIS_setDynparams(n)
             call LIS_get_met_forcing(n)
             call LIS_perturb_forcing(n)
             call LIS_surfaceModel_f2t(n)                  
             call LIS_surfaceModel_run(n)
             call LIS_surfaceModel_perturb_states(n)
             call LIS_readDAobservations(n)
             call LIS_perturb_DAobservations(n)
             call LIS_perturb_writerestart(n)
             call LIS_dataassim_run(n)
             call LIS_dataassim_output(n)
             call LIS_surfaceModel_output(n)
             call LIS_surfaceModel_writerestart(n)
             call LIS_routing_run(n) ! EMK
             call LIS_routing_writeoutput(n) ! EMK
             call LIS_routing_writerestart(n) ! EMK
          endif
       enddo
       call LIS_ticktime 
       flush(LIS_logunit)
    enddo
  end subroutine lis_run_agrmet

  !BOP
! !ROUTINE: lis_final_agrmet
!
! !INTERFACE:
  subroutine lis_final_agrmet
! !USES:
    use LIS_coreMod,         only : LIS_finalize
    use LIS_surfaceModelMod, only : LIS_surfaceModel_finalize
    use LIS_paramsMod,       only : LIS_param_finalize
    use LIS_metforcingMod,   only : LIS_metforcing_finalize
    use LIS_dataAssimMod,    only : LIS_dataassim_finalize
! !DESCRIPTION:
! 
!  This is the finalize method for LIS in the AGRMET running mode. 
!  The following calls are invoked from this method. 
! \begin{description} 
!  \item[LIS\_finalize] (\ref{LIS_finalize}) \newline
!    cleanup LIS generic structures
!  \item[LIS\_lsm\_finalize] (\ref{LIS_lsm_finalize}) \newline
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
    call LIS_dataassim_finalize()

  end subroutine lis_final_agrmet
end module agrmet_runMod
