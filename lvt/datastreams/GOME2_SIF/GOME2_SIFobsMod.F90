!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: GOME2_SIFobsMod
! \label(GOME2_SIFobsMod)
!
! !INTERFACE:
module GOME2_SIFobsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This plugin supports the processing of Solar Induced Fluorescence (SIF)
!  data from GOME-2 (Global Ozone Monitoring Experiment 2) instrument
!  flying on the METOP-A series of satellites (launched Oct 2006). 
! !FILES USED:
!
! !REVISION HISTORY: 
!  6 Jan 2016   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GOME2_SIFobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GOME2SIFobs !Object to hold GOME2SIF observation attributes
!EOP

  type, public :: gome2sifdec
     character*100           :: odir
     character*20            :: version
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: yr
     integer                 :: mo
     logical                 :: startFlag
     real                    :: datares
  end type gome2sifdec
     
  type(gome2sifdec), allocatable :: GOME2SIFObs(:)

contains
  
!BOP
! 
! !ROUTINE: GOME2_SIFobsInit
! \label{GOME2_SIFobsInit}
!
! !INTERFACE: 
  subroutine GOME2_SIFobsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the GOME-2 SIF and NDVI data, 
!   including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    if(.not.allocated(GOME2SIFobs)) then 
       allocate(GOME2SIFobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, GOME2SIFobs(i)%odir, &
         label='GOME2 SIF data directory:',rc=status)
    call LVT_verify(status, 'GOME2 SIF data directory: not defined')
    call ESMF_ConfigGetAttribute(LVT_Config, GOME2SIFobs(i)%version, &
         label='GOME2 SIF data version:',rc=status)
    if(status.ne.0) then 
       write(LVT_logunit,*) '[ERR] GOME2 SIF data version: not defined'
       write(LVT_logunit,*) "[ERR] supported options are 'v26', 'v27' or 'v28'"
       call LVT_endrun()
    endif

    call LVT_update_timestep(LVT_rc, 2592000)

    gome2sifobs(i)%gridDesc = 0
        
    gome2sifobs(i)%nc = 720
    gome2sifobs(i)%nr = 360

    !filling the items needed by the interpolation library
    gome2sifobs(i)%gridDesc(1) = 0  
    gome2sifobs(i)%gridDesc(2) = gome2sifobs(i)%nc
    gome2sifobs(i)%gridDesc(3) = gome2sifobs(i)%nr
    gome2sifobs(i)%gridDesc(4) = -89.75
    gome2sifobs(i)%gridDesc(5) = -179.75
    gome2sifobs(i)%gridDesc(7) = 89.75
    gome2sifobs(i)%gridDesc(8) = 179.75
    gome2sifobs(i)%gridDesc(6) = 128
    gome2sifobs(i)%gridDesc(9) = 0.5
    gome2sifobs(i)%gridDesc(10) = 0.5
    gome2sifobs(i)%gridDesc(20) = 64

    gome2sifobs(i)%datares  = 0.5

    if(LVT_isAtAfinerResolution(gome2sifobs(i)%datares)) then
       
       allocate(gome2sifobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gome2sifobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gome2sifobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(gome2sifobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            gome2sifobs(i)%rlat, &
            gome2sifobs(i)%rlon, &
            gome2sifobs(i)%n11)
    else
       allocate(gome2sifobs(i)%n11(gome2sifobs(i)%nc*gome2sifobs(i)%nr))
       call upscaleByAveraging_input(gome2sifobs(i)%gridDesc,&
            LVT_rc%gridDesc,gome2sifobs(i)%nc*gome2sifobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,gome2sifobs(i)%n11)
    endif

    gome2sifobs(i)%yr = -1
    gome2sifobs(i)%mo = LVT_rc%mo
    gome2sifobs(i)%startFlag = .false. 

  end subroutine GOME2_SIFobsinit


end module GOME2_SIFobsMod
