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
! !MODULE: GLDAS1obsMod
! \label(GLDAS1obsMod)
!
! !INTERFACE:
module GLDAS1obsMod
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLDAS1obsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GLDAS1obs !Object to hold GLDAS1 observation attributes
!EOP

  type, public :: gldas1dec
     character*100           :: odir
     character*20            :: model_name
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: yr
     integer                 :: mo
     logical                 :: startFlag
     real                    :: datares
  end type gldas1dec
     
  type(gldas1dec), allocatable :: GLDAS1Obs(:)

contains
  
!BOP
! 
! !ROUTINE: GLDAS1obsInit
! \label{GLDAS1obsInit}
!
! !INTERFACE: 
  subroutine GLDAS1obsinit(i)
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
!   for reading the GIMMS NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: k, status
    if(.not.allocated(GLDAS1obs)) then 
       allocate(GLDAS1obs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigFindLabel(LVT_Config, &
         label='GLDAS1 data directory:',rc=status)
    do k=1,i
       call ESMF_ConfigGetAttribute(LVT_Config, GLDAS1obs(k)%odir, &
            rc=status)
       call LVT_verify(status, 'GLDAS1 data directory: not defined')
    enddo
       
    call ESMF_ConfigFindLabel(LVT_Config, &
         label='GLDAS1 data model name:',rc=status)
    do k=1,i
       call ESMF_ConfigGetAttribute(LVT_Config, GLDAS1obs(k)%model_name, &
            rc=status)
       call LVT_verify(status, 'GLDAS1 data model name: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_Config, &
         label='GLDAS1 data spatial resolution (degree):',rc=status)
    do k=1,i
       call ESMF_ConfigGetAttribute(LVT_Config, GLDAS1obs(k)%datares, &
            rc=status)
       call LVT_verify(status, 'GLDAS1 data spatial resolution (degree): not defined')
    enddo

    call LVT_update_timestep(LVT_rc, 10800)

    gldas1obs(i)%gridDesc = 0

    if(gldas1obs(i)%datares .eq. 0.25) then 
        
       gldas1obs(i)%nc = 1440
       gldas1obs(i)%nr = 600
       
       gldas1obs(i)%gridDesc(1) = 0  
       gldas1obs(i)%gridDesc(2) = gldas1obs(i)%nc
       gldas1obs(i)%gridDesc(3) = gldas1obs(i)%nr
       gldas1obs(i)%gridDesc(4) = -59.875
       gldas1obs(i)%gridDesc(5) = -179.875
       gldas1obs(i)%gridDesc(7) = 89.875
       gldas1obs(i)%gridDesc(8) = 179.875
       gldas1obs(i)%gridDesc(6) = 128
       gldas1obs(i)%gridDesc(9) = 0.25
       gldas1obs(i)%gridDesc(10) = 0.25
       gldas1obs(i)%gridDesc(20) = 64

       
    elseif(gldas1obs(i)%datares.eq. 1.0) then 

       gldas1obs(i)%nc = 360
       gldas1obs(i)%nr = 150

       gldas1obs(i)%gridDesc(1) = 0  
       gldas1obs(i)%gridDesc(2) = gldas1obs(i)%nc
       gldas1obs(i)%gridDesc(3) = gldas1obs(i)%nr
       gldas1obs(i)%gridDesc(4) = -59.5
       gldas1obs(i)%gridDesc(5) = -179.5
       gldas1obs(i)%gridDesc(7) = 89.5
       gldas1obs(i)%gridDesc(8) = 179.5
       gldas1obs(i)%gridDesc(6) = 128
       gldas1obs(i)%gridDesc(9) = 1.0
       gldas1obs(i)%gridDesc(10) = 1.0
       gldas1obs(i)%gridDesc(20) = 64

    endif

    if(LVT_isAtAfinerResolution(gldas1obs(i)%datares)) then
          
       allocate(gldas1obs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gldas1obs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gldas1obs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(gldas1obs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            gldas1obs(i)%rlat, &
            gldas1obs(i)%rlon, &
            gldas1obs(i)%n11)
    else
       allocate(gldas1obs(i)%n11(gldas1obs(i)%nc*gldas1obs(i)%nr))
       call upscaleByAveraging_input(gldas1obs(i)%gridDesc,&
            LVT_rc%gridDesc,gldas1obs(i)%nc*gldas1obs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,gldas1obs(i)%n11)
    endif

    gldas1obs(i)%yr = -1
    gldas1obs(i)%mo = LVT_rc%mo
    gldas1obs(i)%startFlag = .false. 

  end subroutine GLDAS1obsinit


end module GLDAS1obsMod
