!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: clsmf25_lsmMod.F90
!
! !DESCRIPTION:
!  Module for 1-D catchment land model driver variable initialization
!
! !REVISION HISTORY:
! 15 Dec 2005; Sujay Kumar, Initial Code
! 23 Nov 2012: David Mocko, Additional configs for Catchment Fortuna-2.5
!
! !INTERFACE:
module clsmf25_lsmMod
! !USES:        
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use clsmf25_types
  use clsmf25_drv_types
  use clsmf25_modis_alb_types

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: clsmf25_varalloc
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: clsmf25_struc

  type, public ::  cat_type_dec 
     integer :: catopen             ! Keeps track of opening files
     integer :: numout              ! Counts number of output times for CLSM
     real    :: rstInterval         ! CLSM restart interval (seconds)
     integer :: forc_count
     real    :: ts
     integer :: varid(39)           ! For netcdf output
     real    :: dzsfcrd             ! top soil layer depth from lis.config
     real    :: initSM              ! initial soil moisture
     real    :: initST              ! initial soil temperature
     real    :: RefHfixed           ! reference height [m]
     integer :: turbscheme          ! Catchment turbulence scheme option
     integer :: uselaiflag
     integer :: usemodisalbflag
     integer :: usegreennessflag
     logical :: modelstart

     character(len=LIS_CONST_PATH_LEN) :: rfile         ! restart file
     logical, allocatable               :: good_forcing_mask(:)
     type(cat_param_type), allocatable  :: cat_param(:)
     type(met_force_type), allocatable  :: met_force(:)
     type(cat_progn_type), allocatable  :: cat_progn(:)
     !ag(01Jan2021)
     type(cat_route_type), allocatable  :: cat_route(:)

     type(cat_diagn_type), allocatable  :: cat_diagn(:)
     type(cat_output_type), allocatable :: cat_output(:)
     type(modis_alb_type), allocatable  :: modis_alb_param(:,:)

     
  end type cat_type_dec
  
  type(cat_type_dec), allocatable :: clsmf25_struc(:)
!EOP  

  SAVE
contains
!BOP
! 
! !ROUTINE: clsmf25_varalloc
! 
! !DESCRIPTION:        
! Reads in runtime catchment parameters, allocates memory for variables
! 
! !INTERFACE:
  subroutine clsmf25_varalloc()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use LIS_timeMgrMod
    use LIS_logMod
!EOP
    implicit none
    integer :: n,i
    integer                 :: yr, mo, da, hr, mn, ss
    integer                 :: status
    character*3        :: fnest

    allocate(clsmf25_struc(LIS_rc%nnest))
    call clsmf25_readcrd()
    do n=1,LIS_rc%nnest
       allocate(clsmf25_struc(n)%cat_param(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(clsmf25_struc(n)%met_force(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(clsmf25_struc(n)%cat_progn(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(clsmf25_struc(n)%cat_diagn(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(clsmf25_struc(n)%cat_output(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(clsmf25_struc(n)%good_forcing_mask(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(clsmf25_struc(n)%modis_alb_param(LIS_rc%npatch(n,LIS_rc%lsm_index),12))
       !ag(01Jan2021)
       allocate(clsmf25_struc(n)%cat_route(LIS_rc%npatch(n,LIS_rc%lsm_index)))

       clsmf25_struc(n)%numout = 0 
       clsmf25_struc(n)%forc_count = 0 

       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          clsmf25_struc(n)%met_force(i)%tair = 0 
          clsmf25_struc(n)%met_force(i)%qair = 0 
          clsmf25_struc(n)%met_force(i)%swdown = 0 
          clsmf25_struc(n)%met_force(i)%lwdown = 0
          clsmf25_struc(n)%met_force(i)%wind = 0
          clsmf25_struc(n)%met_force(i)%psurf = 0 
          clsmf25_struc(n)%met_force(i)%rainf = 0 
          clsmf25_struc(n)%met_force(i)%snowf = 0 
          clsmf25_struc(n)%met_force(i)%rainf_c = 0 
          clsmf25_struc(n)%met_force(i)%refh = 0 
          clsmf25_struc(n)%met_force(i)%swnet = 0 
          clsmf25_struc(n)%met_force(i)%pardr = 0 
          clsmf25_struc(n)%met_force(i)%pardf = 0 
       enddo

       !ag(01Jan2021)
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          clsmf25_struc(n)%cat_route(i)%rivsto = 0
          clsmf25_struc(n)%cat_route(i)%fldsto = 0
          clsmf25_struc(n)%cat_route(i)%fldfrc = 0
       enddo

       clsmf25_struc(n)%cat_diagn%totlwc_prev = 0 
       clsmf25_struc(n)%modelstart = .true. 
!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, clsmf25_struc(n)%ts)

       write(fnest,'(i3.3)') n

       call LIS_registerAlarm("CLSM F2.5 model alarm "//trim(fnest),&
            clsmf25_struc(n)%ts,&
            clsmf25_struc(n)%ts)

       call LIS_registerAlarm("CLSM F2.5 restart alarm "//trim(fnest),&
            clsmf25_struc(n)%ts,&
            clsmf25_struc(n)%rstInterval)

       LIS_sfmodel_struc(n)%nsm_layers = 3
       LIS_sfmodel_struc(n)%nst_layers = 6

       allocate(LIS_sfmodel_struc(n)%lyrthk(3))
       LIS_sfmodel_struc(n)%lyrthk(1) = 2         !depth in cm
       LIS_sfmodel_struc(n)%lyrthk(2) = 100
       LIS_sfmodel_struc(n)%lyrthk(3) = 200

       allocate(LIS_sfmodel_struc(n)%lyrthk2(6))  !dzgt in cm
       LIS_sfmodel_struc(n)%lyrthk2(1) = 9.88
       LIS_sfmodel_struc(n)%lyrthk2(2) = 19.52
       LIS_sfmodel_struc(n)%lyrthk2(3) = 38.59
       LIS_sfmodel_struc(n)%lyrthk2(4) = 76.26
       LIS_sfmodel_struc(n)%lyrthk2(5) = 150.71
       LIS_sfmodel_struc(n)%lyrthk2(6) = 1000.0

       LIS_sfmodel_struc(n)%ts = clsmf25_struc(n)%ts

    enddo

  end subroutine clsmf25_varalloc

end module clsmf25_lsmMod



