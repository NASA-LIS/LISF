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
! !MODULE: MERRAlandobsMod
! \label(MERRAlandobsMod)
!
! !INTERFACE:
module MERRAlandobsMod
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
  PUBLIC :: MERRAlandobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MERRAlandobs !Object to hold MERRAland observation attributes
!EOP

  type, public :: merralanddec
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
     real, allocatable       :: qs(:,:,:)
     real, allocatable       :: qsb(:,:,:)
     real, allocatable       :: swnet(:,:,:)
     real, allocatable       :: qle(:,:,:)
     real, allocatable       :: qh(:,:,:)
     real, allocatable       :: frsno(:,:,:)
     real, allocatable       :: snod(:,:,:)
     real, allocatable       :: swe(:,:,:)
     real, allocatable       :: qg(:,:,:)
     real, allocatable       :: sfsm(:,:,:)
     real, allocatable       :: rzsm(:,:,:)
     real, allocatable       :: prcp(:,:,:)
     real, allocatable       :: tskin(:,:,:)
  end type merralanddec
     
  type(merralanddec), allocatable :: MERRAlandObs(:)

contains
  
!BOP
! 
! !ROUTINE: MERRAlandobsInit
! \label{MERRAlandobsInit}
!
! !INTERFACE: 
  subroutine MERRAlandobsinit(i)
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

    integer               :: status

    if(.not.allocated(MERRAlandobs)) then 
       allocate(MERRAlandobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, MERRAlandobs(i)%odir, &
         label='MERRA-Land data directory:',rc=status)
    call LVT_verify(status, 'MERRA-Land data directory: not defined')

    merralandobs(i)%gridDesc = 0
        
    merralandobs(i)%nc = 540
    merralandobs(i)%nr = 361

    !filling the items needed by the interpolation library
    merralandobs(i)%gridDesc(1) = 0  
    merralandobs(i)%gridDesc(2) = merralandobs(i)%nc
    merralandobs(i)%gridDesc(3) = merralandobs(i)%nr
    merralandobs(i)%gridDesc(4) = -90.000
    merralandobs(i)%gridDesc(5) = -180.000
    merralandobs(i)%gridDesc(6) = 128
    merralandobs(i)%gridDesc(7) = 90.000
    merralandobs(i)%gridDesc(8) = 179.33333
    merralandobs(i)%gridDesc(9) =  0.66666666667
    merralandobs(i)%gridDesc(10) = 0.5
    merralandobs(i)%gridDesc(20) = 0

    merralandobs(i)%datares  =  0.66666666667

    if(LVT_isAtAfinerResolution(merralandobs(i)%datares)) then
       
       allocate(merralandobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(merralandobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(merralandobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(merralandobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            merralandobs(i)%rlat, &
            merralandobs(i)%rlon, &
            merralandobs(i)%n11)
    else
       allocate(merralandobs(i)%n11(merralandobs(i)%nc*merralandobs(i)%nr))
       call upscaleByAveraging_input(merralandobs(i)%gridDesc,&
            LVT_rc%gridDesc,merralandobs(i)%nc*merralandobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,merralandobs(i)%n11)
    endif

    call ESMF_TimeIntervalSet(merralandobs(i)%ts, s=3600,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')

    call LVT_update_timestep(LVT_rc, 3600)

    merralandobs(i)%da = -1
    merralandobs(i)%startFlag = .true.

    allocate(merralandobs(i)%qs(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%qsb(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%swnet(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%qle(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%qh(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%frsno(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%snod(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%swe(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%qg(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%sfsm(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%rzsm(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%prcp(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merralandobs(i)%tskin(LVT_rc%lnc,LVT_rc%lnr,24))

    merralandobs(i)%qs = LVT_rc%udef
    merralandobs(i)%qsb = LVT_rc%udef
    merralandobs(i)%swnet = LVT_rc%udef
    merralandobs(i)%qle = LVT_rc%udef
    merralandobs(i)%qh = LVT_rc%udef
    merralandobs(i)%frsno = LVT_rc%udef
    merralandobs(i)%snod = LVT_rc%udef
    merralandobs(i)%swe = LVT_rc%udef
    merralandobs(i)%qg = LVT_rc%udef
    merralandobs(i)%sfsm = LVT_rc%udef
    merralandobs(i)%rzsm = LVT_rc%udef
    merralandobs(i)%prcp = LVT_rc%udef
    merralandobs(i)%tskin = LVT_rc%udef
    
  end subroutine MERRAlandobsinit


end module MERRAlandobsMod
