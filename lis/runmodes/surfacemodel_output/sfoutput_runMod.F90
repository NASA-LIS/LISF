!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module sfoutput_runMod
!BOP
!
! !MODULE: sfoutput_runMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   LIS initialization, execution, and finalization
!   for the surface model runmode used in LIS.
!   
! !REVISION HISTORY: 
!  30May19    Daniel Rosen Initial Specification
!
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: lis_init_sfoutput  !init method for sfoutput run mode
  public :: lis_run_sfoutput   !run method for sfoutput run mode
  public :: lis_final_sfoutput !finalize method for sfoutput run mode

!EOP  
contains
!BOP
! !ROUTINE: lis_init_sfoutput
!
! 
! !INTERFACE:
  subroutine lis_init_sfoutput
! !USES: 
    use LIS_coreMod
    use LIS_domainMod,         only : LIS_domain_init
    use LIS_surfaceModelMod,   only : LIS_surfaceModel_init,  &
                                      LIS_surfaceModel_setup
    use LIS_metforcingMod,     only : LIS_metforcing_init
    use LIS_DAobservationsMod, only : LIS_initDAObservations
    use LIS_perturbMod,        only : LIS_perturb_init
    use LIS_dataAssimMod,      only : LIS_dataassim_init
    use LIS_paramsMod,         only : LIS_param_init

! !DESCRIPTION:
!  This is the initialize method for LIS in the sfoutput running mode. 
!  The following calls are invoked from this method. 
! \begin{description} 
!  \item[LIS\_domain\_init] (\ref{LIS_domain_init}) \newline
!    initialize the LIS domains
!  \item[LIS\_param\_init] (\ref{LIS_param_init}) \newline
!    initialize parameters
!  \item[LIS\_perturb\_init] (\ref{LIS_perturb_init}) \newline
!    initialize perturbations
!  \item[LIS\_surfaceModel\_init] (\ref{LIS_surfaceModel_init}) \newline
!    initialize the surface model
!  \item[LIS\_metforcing\_init] (\ref{LIS_metforcing_init}) \newline
!    initialize the met forcing
!  \item[LIS\_initDAObservations] (\ref{LIS_initDAobservations}) \newline
!    initialize structures needed to read observations for 
!    data assimilation
!  \item[LIS\_dataassim\_init] (\ref{LIS_dataassim_init}) \newline
!    initialize structures needed for data assimilation
!  \item[LIS\_surfaceModel\_setup] (\ref{LIS_surfaceModel_setup}) \newline
!    complete the surface model setup
!  \item[LIS\_core\_init] (\ref{LIS_core_init}) \newline
!    complete the LIS setup
! \end{description} 
!
!EOP

    LIS_rc%run_model = .true.      !always march the model

    call LIS_domain_init
    call LIS_param_init
    call LIS_perturb_init
    call LIS_surfaceModel_init
    call LIS_metforcing_init
    call LIS_initDAObservations
    call LIS_dataassim_init
    call LIS_surfaceModel_setup
    call LIS_core_init

  end subroutine lis_init_sfoutput

!BOP
! !ROUTINE: lis_run_sfoutput
!
! !INTERFACE:
  subroutine lis_run_sfoutput
! !USES:
    use LIS_coreMod,           only : LIS_rc, LIS_endofrun, LIS_timeToRunNest, &
                                      LIS_ticktime
    use LIS_surfaceModelMod,   only : LIS_surfaceModel_output
    use LIS_logMod,            only : LIS_flush, LIS_logunit
!
! !DESCRIPTION:
! 
!  This is the run method for LIS in the sfoutput running mode. 
!  The following calls are invoked from this method. 
! \begin{description} 
!  \item[LIS\_endofrun] (\ref{LIS_endofrun}) \newline
!    check to see if the end of simulation has reached
!  \item[LIS\_ticktime] (\ref{LIS_ticktime}) \newline
!    advance model clock
!  \item[LIS\_timeToRunNest] (\ref{LIS_timeToRunNest}) \newline
!    check to see if the current nest needs to be run. 
!  \item[LIS\_surfaceModel\_output] (\ref{LIS_surfaceModel_output}) \newline
!    write surface model output
! \end{description} 
!
!EOP
    integer :: n,ierr

    do while (.NOT. LIS_endofrun())

       do n=1,LIS_rc%nnest
          if(LIS_timeToRunNest(n)) then 
             call LIS_surfaceModel_output(n)
          endif
       enddo
       call LIS_ticktime 
       call LIS_flush(LIS_logunit)
    enddo
  end subroutine lis_run_sfoutput

  !BOP
! !ROUTINE: lis_final_sfoutput
!
! !INTERFACE:
  subroutine lis_final_sfoutput
! !USES:
    use LIS_coreMod,         only : LIS_finalize
    use LIS_surfaceModelMod, only : LIS_surfaceModel_finalize
    use LIS_paramsMod,       only : LIS_param_finalize
    use LIS_metforcingMod,   only : LIS_metforcing_finalize
    use LIS_dataAssimMod,    only : LIS_dataassim_finalize
! !DESCRIPTION:
! 
!  This is the finalize method for LIS in the sfoutput running mode. 
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

  end subroutine lis_final_sfoutput
end module sfoutput_runMod
