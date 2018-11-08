!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module smootherDA_runMod
!BOP
!
! !MODULE: smootherDA_runMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   LIS initialization, execution, and finalization
!   for an optimization and uncertainty estimation runmode in LIS.
!   
! !REVISION HISTORY: 
!  21Oct05    Sujay Kumar  Initial Specification
! 
!
  use ESMF

  implicit none
  
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: lis_init_smootherDA  !init method for optimization-uncertainty modeling mode
  public :: lis_run_smootherDA   !run method for optimization-uncertainty modeling mode
  public :: lis_final_smootherDA !finalize method for optimization-uncertainty modeling mode
!EOP  
contains
!BOP
! !ROUTINE: lis_init_smootherDA
!
! 
! !INTERFACE:
  subroutine lis_init_smootherDA
! !USES: 
    use LIS_coreMod
    use LIS_domainMod
    use LIS_surfaceModelMod
    use LIS_metforcingMod
    use LIS_DAobservationsMod
    use LIS_perturbMod
    use LIS_dataAssimMod
    use LIS_paramsMod
    use LIS_irrigationMod
    use LIS_routingMod
    use LIS_logMod

! !DESCRIPTION:
!  This is the initialize method for LIS in a retrospective running mode. 
!  The following calls are invoked from this method. 
! \begin{description} 
!  \item[LIS\_domain\_init] (\ref{LIS_domain_init}) \newline
!    initialize the LIS domains
!  \item[LIS\_param\_init] (\ref{LIS_param_init}) \newline
!    initialize parameters
!  \item[LIS\_surfaceModel\_init] (\ref{LIS_surfaceModel_init}) \newline
!    initialize the land surface model
!  \item[LIS\_metforcing\_init] (\ref{LIS_metforcing_init}) \newline
!    initialize the met forcing
!  \item[LIS\_initDAObservations] (\ref{LIS_initDAobservations}) \newline
!    initialize structures needed to read observations for 
!    data assimilation
!  \item[LIS\_surfaceModel\_setup] (\ref{LIS_surfaceModel_setup}) \newline
!    complete the LSM setups
!  \item[LIS\_surfaceModel\_readrestart] (\ref{LIS_surfaceModel_readrestart}) \newline
!    read the surface model restart files 
!  \item[LIS\_perturb\_readrestart] (\ref{LIS_perturb_readrestart}) \newline
!    read the perturbation restart files 
!  \item[LIS\_routing\_init] (\ref{LIS_routing_init}) \newline
!    initialize the routing model
!  \item[LIS\_routing\_readrestart] (\ref{LIS_routing_readrestart}) \newline
!    read the routing model restart files 
!  \item[LIS\_core\_init] (\ref{LIS_core_init}) \newline
!    completes the setting of LIS' clocks and alarms
! \end{description} 
!
!EOP
    call LIS_domain_init
    call LIS_param_init
    call LIS_perturb_init
    call LIS_surfaceModel_init
    call LIS_metforcing_init
    call LIS_irrigation_init
    call LIS_initDAObservations
    call LIS_dataassim_init
    call LIS_surfaceModel_setup
    call LIS_surfaceModel_readrestart
    call LIS_perturb_readrestart
    call LIS_routing_init
    call LIS_routing_readrestart
    call LIS_core_init

  end subroutine lis_init_smootherDA

!BOP
! !ROUTINE: lis_run_smootherDA
!
! !INTERFACE:
  subroutine lis_run_smootherDA
! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_domainMod
    use LIS_surfaceModelMod
    use LIS_paramsMod
    use LIS_metforcingMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_dataAssimMod
    use LIS_routingMod
    use LIS_irrigationMod
    use LIS_logMod

!
! !DESCRIPTION:
! 
!  This is the run method for LIS in a smootherDA running mode. 
!
!  NOTES: 
!  If the endtime is set to exactly the beginning of a month, 
!  the GRACE analysis will not be applied to the last month
!  of the simulation. The recommended workaround is to specify
!  the end time to be not exactly at the end of the month. 
!
!EOP
    integer           :: n 

    LIS_rc%DAincrMode(:) = 0 
    
    call LIS_surfaceModel_readrestart
    do while (.NOT. LIS_endofrun())
       do while(.NOT.LIS_endofTimeWindow())
          call LIS_ticktime
          do n=1,LIS_rc%nnest
             if(LIS_timeToRunNest(n)) then
                call LIS_setDynparams(n)
                call LIS_get_met_forcing(n)
                call LIS_perturb_forcing(n)
                if(LIS_rc%DAincrMode(n).eq.1) then
                   call LIS_irrigation_run(n)
                endif
                call LIS_surfaceModel_f2t(n)  
                call LIS_surfaceModel_run(n)
                call LIS_surfaceModel_diagnoseVarsForDA(n)
                call LIS_surfaceModel_perturb_states(n)
                if(LIS_rc%DAincrMode(n).eq.1) then
                  call LIS_readDAobservations(n)
                  call LIS_perturb_DAobservations(n)   
                end if
                call LIS_perturb_writerestart(n)
                call LIS_dataassim_run(n)
                call LIS_dataassim_output(n)
                call LIS_surfaceModel_output(n)  
                call LIS_surfaceModel_writerestart(n)                
                if(LIS_rc%DAincrMode(n).eq.1) then
                   call LIS_irrigation_output(n)
                endif
                call LIS_routing_run(n)
                call LIS_routing_writeoutput(n)
                call LIS_routing_writerestart(n)
                call updateIncrementsFlag(n)
             endif
          enddo
          call LIS_flush(LIS_logunit)
       enddo

       if(LIS_rc%endtime.ne.1) then 
          call LIS_resetClockForTimeWindow(LIS_rc)
          call LIS_param_reset
          call LIS_metforcing_reset
          call LIS_surfaceModel_readrestart
          call LIS_routing_readrestart
          call LIS_perturb_readrestart 

          LIS_rc%iterationId(:) = LIS_rc%iterationId(:) + 1
       endif
    enddo

  end subroutine lis_run_smootherDA

  subroutine updateIncrementsFlag(n)
    use LIS_timeMgrMod
    use LIS_coreMod,  only: LIS_rc
    use LIS_logMod

    type(ESMF_time) :: currTime
    integer         :: status, n
    integer         :: yr,mo,da,hr,mn,ss

    call ESMF_ClockGet(LIS_clock, currTime = currTime, rc=status)
    if(currTime .gt. LIS_twMidTime) then 
       LIS_rc%DAincrMode(n) = 0
    else
       LIS_rc%DAincrMode(n) = 1
    endif
  end subroutine updateIncrementsFlag

 

!BOP
! !ROUTINE: lis_final_smootherDA
!
! !INTERFACE:
  subroutine lis_final_smootherDA
! !USES:
    use LIS_coreMod
    use LIS_surfaceModelMod
    use LIS_paramsMod
    use LIS_metforcingMod


! !DESCRIPTION:
! 
!  This is the finalize method for LIS in a smootherDA running mode. 
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
! \end{description} 
!
!EOP

    call lis_finalize()
    call LIS_surfaceModel_finalize
    call LIS_param_finalize()
    call LIS_metforcing_finalize()

  end subroutine lis_final_smootherDA
end module smootherDA_runMod
