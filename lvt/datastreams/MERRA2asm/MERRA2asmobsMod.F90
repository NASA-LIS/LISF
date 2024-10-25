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
! !MODULE: MERRA2asmobsMod
! \label(MERRA2asmobsMod)
!
! !INTERFACE:
module MERRA2asmobsMod
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
!  3 June 2017   Eric Kemp  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MERRA2asmobsinit !Initializes structures for MERRA2 asm data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MERRA2asmobs !Object to hold MERRA2 asm observation attributes
!EOP

  type, public :: merra2asmdec
     character*100           :: odir
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: da
     real                    :: datares
     logical                 :: startFlag
     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: ts
     real, allocatable       :: t2m(:,:,:)
     real, allocatable       :: qv2m(:,:,:)
     real, allocatable       :: u10m(:,:,:)
     real, allocatable       :: v10m(:,:,:)
  end type merra2asmdec
     
  type(merra2asmdec), allocatable :: MERRA2asmObs(:)
  
contains
  
!BOP
! 
! !ROUTINE: MERRA2asmobsInit
! \label{MERRA2asmobsInit}
!
! !INTERFACE: 
  subroutine MERRA2asmobsinit(i)
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
!   for reading MERRA2 asm (single level diagnostic) data, including the
!   computation of spatial interpolation weights.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: k,status

    if(.not.allocated(MERRA2asmobs)) then 
       allocate(MERRA2asmobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigFindLabel(LVT_Config, &
         label='MERRA2 ASM data directory:',rc=status)
    do k=1,i
       call ESMF_ConfigGetAttribute(LVT_Config, MERRA2asmobs(k)%odir, &
            rc=status)
       call LVT_verify(status, 'MERRA2 ASM data directory: not defined')
    enddo

    merra2asmobs(i)%gridDesc = 0
        
    merra2asmobs(i)%nc = 576
    merra2asmobs(i)%nr = 361

    !filling the items needed by the interpolation library
    merra2asmobs(i)%gridDesc(1) = 0  
    merra2asmobs(i)%gridDesc(2) = merra2asmobs(i)%nc
    merra2asmobs(i)%gridDesc(3) = merra2asmobs(i)%nr
    merra2asmobs(i)%gridDesc(4) = -90.000
    merra2asmobs(i)%gridDesc(5) = -180.000
    merra2asmobs(i)%gridDesc(7) = 90.000
    merra2asmobs(i)%gridDesc(8) = 179.375
    merra2asmobs(i)%gridDesc(6) = 128
    merra2asmobs(i)%gridDesc(9) = 0.625
    merra2asmobs(i)%gridDesc(10) = 0.5
    merra2asmobs(i)%gridDesc(20) = 0

    merra2asmobs(i)%datares  = 0.625

    if(LVT_isAtAfinerResolution(merra2asmobs(i)%datares)) then
       
       allocate(merra2asmobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(merra2asmobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(merra2asmobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(merra2asmobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            merra2asmobs(i)%rlat, &
            merra2asmobs(i)%rlon, &
            merra2asmobs(i)%n11)
    else
       allocate(merra2asmobs(i)%n11(merra2asmobs(i)%nc*merra2asmobs(i)%nr))
       call upscaleByAveraging_input(merra2asmobs(i)%gridDesc,&
            LVT_rc%gridDesc,merra2asmobs(i)%nc*merra2asmobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,merra2asmobs(i)%n11)
    endif

    call ESMF_TimeIntervalSet(merra2asmobs(i)%ts, s=3600,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: MERRA2asmobsMod ')

    call LVT_update_timestep(LVT_rc, 3600)

    merra2asmobs(i)%da = -1
    merra2asmobs(i)%startFlag = .true.

    allocate(merra2asmobs(i)%t2m(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2asmobs(i)%qv2m(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2asmobs(i)%u10m(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2asmobs(i)%v10m(LVT_rc%lnc,LVT_rc%lnr,24))

    merra2asmobs(i)%t2m = LVT_rc%udef
    merra2asmobs(i)%qv2m = LVT_rc%udef
    merra2asmobs(i)%u10m = LVT_rc%udef
    merra2asmobs(i)%v10m = LVT_rc%udef
    
 end subroutine MERRA2asmobsinit
 
end module MERRA2asmobsMod
