!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: SMAP_L3TBMod
! \label(SMAP_L3TBMod)
!
! !INTERFACE:
module SMAP_L3TBMod
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
!  This module handles the observation plugin for the 
!  retrieval products from the Soil Moisture Active Passive
!  (SMAP) mission
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  19 Nov 2018: Mahdi Navari , Sujay Kumar, Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_L3TBinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_L3TB
!EOP
  type, public :: smapobsdec

     character*100        :: odir
     character*20         :: data_designation

     integer              :: nc
     integer              :: nr
     integer              :: mo
     integer,allocatable  :: n112(:)
     real,allocatable     :: rlat2(:)
     real,allocatable     :: rlon2(:)
     real                 :: gridDesci(50)

     real,    allocatable    :: L3TB(:,:)
     real,    allocatable    :: Tbtime(:,:)
     integer*2, allocatable  :: Tbqc(:,:)     
     logical                 :: startflag     

  end type smapobsdec

  type(smapobsdec), allocatable:: SMAP_L3TB(:)

contains
  
!BOP
! 
! !ROUTINE: SMAP_L3TBInit
! \label{SMAP_L3TBInit}
!
! !INTERFACE: 
  subroutine SMAP_L3TBinit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading NASA AMSR-E soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer         ezlh_convert
    integer            :: npts
    integer                 :: ease_nc,ease_nr
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc

    if(.not.allocated(SMAP_L3TB)) then 
       allocate(SMAP_L3TB(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_L3TB(i)%odir, &
         label='SMAP soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'SMAP soil moisture observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_L3TB(i)%data_designation, &
         label='SMAP soil moisture data designation:', rc=status)
    call LVT_verify(status, &
         'SMAP soil moisture data designation: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    if(SMAP_L3TB(i)%data_designation.eq."SPL3SMAP") then 
       !SMAP L3 radar/radiometer daily 9km
       SMAP_L3TB(i)%gridDesci=0.0 

       ease_nr=1624
       ease_nc=3856       

    !filling the items needed by the interpolation library
       SMAP_L3TB(i)%gridDesci(1) = 9  !input is EASE grid
    !these  corner coordinates were calculated based on ezlh_convert
       SMAP_L3TB(i)%gridDesci(4) = -84.6564  !lat
       SMAP_L3TB(i)%gridDesci(5) = -179.953 !lon
       SMAP_L3TB(i)%gridDesci(7) = 84.6564  !lat
       SMAP_L3TB(i)%gridDesci(8) = 179.953  !lon
       SMAP_L3TB(i)%gridDesci(9) = 5 !M09 grid
       SMAP_L3TB(i)%gridDesci(20) = 64
       
       SMAP_L3TB(i)%gridDesci(2) = ease_nc  !nx
       SMAP_L3TB(i)%gridDesci(3) = ease_nr  !ny

       SMAP_L3TB(i)%nc=ease_nc
       SMAP_L3TB(i)%nr=ease_nr       
    elseif(SMAP_L3TB(i)%data_designation.eq."SPL3SMP") then 
       !SMAP L3 radiometer daily 36km 
    !filling the items needed by the interpolation library
       ease_nc=964       
       ease_nr=406
       SMAP_L3TB(i)%gridDesci = 0
       SMAP_L3TB(i)%gridDesci(1) = 9  !input is EASE grid
       SMAP_L3TB(i)%gridDesci(2) = ease_nc  !nx
       SMAP_L3TB(i)%gridDesci(3) = ease_nr  !ny
       SMAP_L3TB(i)%gridDesci(9) = 4 !M36 grid
       SMAP_L3TB(i)%gridDesci(20) = 64
       SMAP_L3TB(i)%gridDesci(10) = 0.36
       SMAP_L3TB(i)%gridDesci(11) = 1

       SMAP_L3TB(i)%nc=ease_nc
       SMAP_L3TB(i)%nr=ease_nr       
    elseif(SMAP_L3TB(i)%data_designation.eq."SPL3SMP_E") then 
       ease_nc=3856       
       ease_nr=1624
       SMAP_L3TB(i)%gridDesci = 0
       SMAP_L3TB(i)%gridDesci(1) = 9  
       SMAP_L3TB(i)%gridDesci(2) = ease_nc  !nx
       SMAP_L3TB(i)%gridDesci(3) = ease_nr  !ny
       SMAP_L3TB(i)%gridDesci(9) = 5 !M09 grid
       SMAP_L3TB(i)%gridDesci(20) = 64
       SMAP_L3TB(i)%gridDesci(10) = 0.09
       SMAP_L3TB(i)%gridDesci(11) = 1

       SMAP_L3TB(i)%nc=ease_nc
       SMAP_L3TB(i)%nr=ease_nr       
    endif

    npts= LVT_rc%lnc*LVT_rc%lnr
    SMAP_L3TB(i)%mo=npts
    
    
    allocate(SMAP_L3TB(i)%rlat2(npts))
    allocate(SMAP_L3TB(i)%rlon2(npts))
    allocate(SMAP_L3TB(i)%n112(npts))
    
    SMAP_L3TB(i)%rlat2=0.0
    SMAP_L3TB(i)%rlon2=0.0
    SMAP_L3TB(i)%n112=0.0
    call neighbor_interp_input(SMAP_L3TB(i)%gridDesci,LVT_rc%gridDesc,&
         npts,SMAP_L3TB(i)%rlat2,&
         SMAP_L3TB(i)%rlon2,SMAP_L3TB(i)%n112)
    

    SMAP_L3TB(i)%startflag = .true. 

    allocate(SMAP_L3TB(i)%L3TB(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAP_L3TB(i)%Tbtime(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAP_L3TB(i)%Tbqc(LVT_rc%lnc*LVT_rc%lnr,2))

!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------

  end subroutine SMAP_L3TBinit


end module SMAP_L3TBMod
