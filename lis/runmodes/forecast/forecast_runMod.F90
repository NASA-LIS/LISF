!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module forecast_runMod
!BOP
!
! !MODULE: forecast_runMod
! 
! !DESCRIPTION: 
!   This is a typical running mode used in LIS to perform 
!   forecast simulations. The mode assumes that all the meteorological
!   forcing analyes required for the specified time period is archived 
!   for the simulation. 
! 
!   This module contains the definition of the functions used for
!   LIS initialization, execution, and finalization
!   for a forecast runmode in LIS.
!   
! !REVISION HISTORY: 
!  19 Feb 2016   Sujay Kumar  Initial Specification
! 
!
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_init_forecast  !init method for forecast mode
  public :: LIS_run_forecast   !run method for forecast mode
  public :: LIS_final_forecast !finalize method for forecast mode

!EOP  
contains
!BOP
! !ROUTINE: LIS_init_forecast
!
! 
! !INTERFACE:
  subroutine LIS_init_forecast
! !USES: 
    use LIS_coreMod
    use LIS_domainMod
    use LIS_surfaceModelMod
    use LIS_metforcingMod
    use LIS_perturbMod
    use LIS_paramsMod
    use LIS_routingMod
    use LIS_irrigationMod
    use LIS_appMod
    use LIS_forecastMod
    use LIS_tbotAdjustMod
!
! !DESCRIPTION:
!  This is the initialize method for LIS in a forecast running mode. 
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
    call LIS_createTmnUpdate
    call LIS_param_init
    call LIS_perturb_init
    call LIS_surfaceModel_init
    call LIS_forecast_init
    call LIS_metforcing_init
    call LIS_irrigation_init
    call LIS_surfaceModel_setup
    call LIS_surfaceModel_readrestart
    call LIS_perturb_readrestart
    call LIS_routing_init
    call LIS_routing_readrestart
    call LIS_appModel_init
    call LIS_core_init

  end subroutine lis_init_forecast

!BOP
! !ROUTINE: lis_run_forecast
!
! !INTERFACE:
  subroutine lis_run_forecast
! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_surfaceModelMod
    use LIS_paramsMod
    use LIS_metforcingMod
    use LIS_perturbMod
    use LIS_routingMod
    use LIS_irrigationMod
    use LIS_appMod
    use LIS_RTMMod
    use LIS_logMod
    use LIS_forecastMod
!
! !DESCRIPTION:
! 
!  This is the run method for LIS in a forecast running mode. 
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
!  \item[LIS\_irrigation\_run] (\ref{LIS_irrigation_run}) \newline
!    run the irrigation model
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
!  \item[LIS\_irrigation\_output] (\ref{LIS_irrigation_output}) \newline
!    write irrigation model output
!  \item[LIS\_routing\_run] (\ref{LIS_routing_run}) \newline
!    run the routing model
!  \item[LIS\_routing\_writeoutput] (\ref{LIS_routing_writeoutput}) \newline
!    write routing model output
!  \item[LIS\_routing\_writerestart] (\ref{LIS_routing_writerestart}) \newline
!    write routing model restart data
!  \item[LIS\_RTM\_run] (\ref{LIS_RTM_run}) \newline
!    run the radiative transfer model
!  \item[LIS\_RTM\_output] (\ref{LIS_RTM_output}) \newline
!    write radiative transfer model output
!  \item[LIS\_runAppModel] (\ref{LIS_runAppModel}) \newline
!    run the application model
!  \item[LIS\_outputAppModel] (\ref{LIS_outputAppModel}) \newline
!    write the application model output
! \end{description} 
!
!EOP
    integer    :: n
    integer    :: i 
    
!    do i=LIS_forecast_struc(1)%st_iterId,LIS_forecast_struc(1)%niterations
!       LIS_forecast_struc(:)%iterId = i
       do while (.NOT. LIS_endofrun())
          call LIS_ticktime 
          do n=1,LIS_rc%nnest
             if(LIS_timeToRunNest(n)) then 
                call LIS_setDynparams(n)
                call LIS_get_met_forcing(n)
                call LIS_perturb_forcing(n)
                call LIS_irrigation_run(n)
                call LIS_surfaceModel_f2t(n)  
                call LIS_surfaceModel_run(n)
                call LIS_surfaceModel_perturb_states(n)
                call LIS_perturb_writerestart(n)
                call LIS_surfaceModel_output(n)
                call LIS_surfaceModel_writerestart(n)
                call LIS_irrigation_output(n)
                call LIS_routing_run(n)
                call LIS_routing_writeoutput(n)
                call LIS_routing_writerestart(n)
                call LIS_runAppModel(n)             
                call LIS_outputAppModel(n)
             endif
          enddo
          flush(LIS_logunit)
       enddo

 !      call LIS_forecast_writerestart
 !      call LIS_resetClock(LIS_rc)
 !      call LIS_param_reset
 !      call LIS_metforcing_reset
 !      call LIS_surfaceModel_reset
 !      call LIS_surfaceModel_readrestart
 !      call LIS_perturb_readrestart
 !      call LIS_core_init   
 !
 !   enddo
  end subroutine lis_run_forecast

!BOP
! !ROUTINE: lis_final_forecast
!
! !INTERFACE:
  subroutine lis_final_forecast
! !USES:
    use LIS_coreMod
    use LIS_logMod
    use LIS_surfaceModelMod
    use LIS_paramsMod
    use LIS_metforcingMod
    use LIS_RTMMod
    use LIS_appMod

! !DESCRIPTION:
! 
!  This is the finalize method for LIS in a forecast running mode. 
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
!  \item[LIS\_RTM\_finalize] (\ref{LIS_RTM_finalize}) \newline
!    cleanup radiative transfer model specific structures
!  \item[LIS\_dataassim\_finalize] (\ref{LIS_dataassim_finalize}) \newline
!    cleanup data assimilation specific structures
! \end{description} 
!
!EOP

    call lis_finalize()
    call LIS_surfaceModel_finalize()
    call LIS_param_finalize()
    call LIS_metforcing_finalize()
    call LIS_RTM_finalize()
    call LIS_appModel_finalize()     

    write(LIS_logunit,*) "[INFO] LIS Run completed. "

  end subroutine lis_final_forecast
end module forecast_runMod
