!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine lis_gce_init(clock, vm)
  
  use ESMF
  use LIS_coreMod,         only : LIS_config_init!, lisGrid
  use LIS_domainMod,       only : LIS_domain_init
  use LIS_surfaceModelMod, only : LIS_surfaceModel_init, &
       LIS_surfaceModel_setup, LIS_surfaceModel_readrestart
  use LIS_metforcingMod,  only : LIS_metforcing_init
  use LIS_DAobservationsMod, only : LIS_initDAObservations
  use LIS_perturbMod,      only : LIS_perturb_init
  use LIS_dataAssimMod,    only : LIS_dataassim_init
  use LIS_paramsMod,        only : LIS_param_init
  use lisgceGridCompMod,   only : lisgce_alloc_states
  implicit none

  type(ESMF_VM)        :: vm
  type(ESMF_Clock)     :: clock
  integer              :: n
  integer              :: coupled
    
  coupled = 1

  call LIS_Config_init(vm=vm, clock=clock)
  call LIS_domain_init
  call LIS_param_init
  call LIS_perturb_init
  call LIS_surfaceModel_init
  call LIS_metforcing_init(coupled)
  call LIS_initDAObservations
  call LIS_dataassim_init
  call LIS_surfaceModel_setup
  call LIS_surfaceModel_readrestart
  call lisgce_alloc_states
  print*, 'Reached at the end of LIS_init'

  
end subroutine lis_gce_init
