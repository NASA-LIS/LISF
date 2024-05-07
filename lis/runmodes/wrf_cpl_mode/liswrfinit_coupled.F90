!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine liswrfinit_coupled(nx,ny,wrf_mpi_comm_compute_tasks) 
!BOP
! !ROUTINE: liswrfinit_coupled
! 
!  !REVISION HISTORY: 
!  14Nov02    Sujay Kumar  Initial Specification
!  21Oct05    Sujay Kumar  Modified to include the runmodes. Switched
!                          to a init,run,finalize mode
! !USES: 
  use LIS_coreMod
  use LIS_domainMod,       only : LIS_domain_init
  use LIS_surfaceModelMod, only : LIS_surfaceModel_init, &
       LIS_surfaceModel_setup, LIS_surfaceModel_readrestart
  use LIS_metforcingMod,  only : LIS_metforcing_init
  use LIS_DAobservationsMod, only : LIS_initDAObservations
  use LIS_perturbMod,      only : LIS_perturb_init, LIS_perturb_readrestart
  use LIS_dataAssimMod,    only : LIS_dataassim_init
  use LIS_paramsMod,        only : LIS_param_init
  use LISWRFGridCompMod,    only : LISWRF_alloc_states
  use LIS_logMod,           only : LIS_logunit

  use LIS_tbotAdjustMod,    only : LIS_createTmnUpdate
  use LIS_irrigationMod,    only : LIS_irrigation_init ! EMK

    implicit none
! !ARGUMENTS: 
    integer, intent(in)     :: nx, ny
    integer, intent(in)     :: wrf_mpi_comm_compute_tasks
!
! !DESCRIPTION: 
!  This routine defines the set of steps required from LIS during the 
!  initialization of the LIS-WRF system.  
!
! The arguments are:
! \begin{description}
! \item[nx]
!   number of processes in the x-direction.
! \item[ny]
!   number of processes in the y-direction.
! \item[comm]
!   MPI communicator to use.
! \end{description}
! 
!EOP  
    integer                 :: coupled
    

    if(.not.LIS_initialized) then
       call LIS_config_init(nx=nx,ny=ny,comm=wrf_mpi_comm_compute_tasks)

       call LIS_domain_init
       call LIS_createTmnUpdate
       call LIS_param_init
       call LIS_perturb_init
       call LIS_surfaceModel_init

       LIS_rc%met_nf(:) = 17 

       call LIS_metforcing_init(coupled)
       call LIS_irrigation_init ! EMK
       call LIS_initDAObservations
       call LIS_dataassim_init
       call LIS_surfaceModel_setup
       call LIS_surfaceModel_readrestart
       call LIS_perturb_readrestart
       call LISWRF_alloc_states
       call LIS_core_init

       flush(LIS_logunit)
       LIS_initialized = .true.

    endif
end subroutine liswrfinit_coupled
