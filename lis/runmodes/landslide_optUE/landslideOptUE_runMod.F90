!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module landslideOptUE_runMod
!BOP
!
! !MODULE: landslideOptUE_runMod
! 
! !DESCRIPTION: 
!   This is a typical running mode used in LIS to perform 
!   landslide simulations. The mode assumes that all the meteorological
!   forcing analyes required for the specified time period is archived 
!   for the simulation. 
! 
!   This module contains the definition of the functions used for
!   LIS initialization, execution, and finalization
!   for a landslide runmode in LIS.
!   
! !REVISION HISTORY: 
!  21Oct05    Sujay Kumar  Initial Specification
! 
!
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_init_landslideOptUE  !init method for landslide mode
  public :: LIS_run_landslideOptUE   !run method for landslide mode
  public :: LIS_final_landslideOptUE !finalize method for landslide mode

!EOP  
contains
!BOP
! !ROUTINE: LIS_init_landslideOptUE
!
! 
! !INTERFACE:
  subroutine LIS_init_landslideOptUE
! !USES: 
    use LIS_domainMod,       only : LIS_domain_init
    use LIS_lsmMod,          only : LIS_lsm_init, LIS_setuplsm, &
         LIS_readrestart
    use LIS_metforcingMod,  only : LIS_metforcing_init
    use LIS_DAobservationsMod, only : LIS_initDAObservations
    use LIS_perturbMod,      only : LIS_perturb_init
    use LIS_dataAssimMod,    only : LIS_dataassim_init
    use LIS_paramsMod,        only : LIS_param_init
    use LIS_landslideMod,     only : LIS_initLandSlideModel
    use LIS_optUEMod, only : LIS_optUE_init, &
         LIS_optUE_setup

! !DESCRIPTION:
!  This is the initialize method for LIS in a landslide running mode. 
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
    call LIS_domain_init
    call LIS_param_init
    call LIS_perturb_init
    call LIS_lsm_init
    call LIS_metforcing_init
    call LIS_initDAObservations
    call LIS_dataassim_init
    call LIS_setuplsm
    call LIS_readrestart
    call LIS_initLandSlideModel
    call LIS_optUE_init
    call LIS_optUE_setup

  end subroutine lis_init_landslideOptUE

!BOP
! !ROUTINE: lis_run_landslideOptUE
!
! !INTERFACE:
  subroutine lis_run_landslideOptUE
! !USES:
    use LIS_coreMod,         only : LIS_rc
    use LIS_optUEMod, only : LIS_isOptStopCriterionTrue,&
         LIS_runOptUE
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
       call LIS_runOptUE()
    enddo

  end subroutine lis_run_landslideOptUE
!BOP
! !ROUTINE: lis_final_landslideOptUE
!
! !INTERFACE:
  subroutine lis_final_landslideOptUE
! !USES:
    use LIS_coreMod,         only : LIS_finalize
    use LIS_lsmMod,          only : LIS_lsm_finalize
    use LIS_paramsMod,        only : LIS_param_finalize
    use LIS_metforcingMod,  only : LIS_metforcing_finalize

! !DESCRIPTION:
! 
!  This is the finalize method for LIS in a landslide running mode. 
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
    call LIS_lsm_finalize()
    call LIS_param_finalize()
    call LIS_metforcing_finalize()
!    call LIS_dataassim_finalize()

  end subroutine lis_final_landslideOptUE
end module landslideOptUE_runMod
