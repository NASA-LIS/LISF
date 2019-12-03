!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: SMAP_smobsMod
! \label(SMAP_smobsMod)
!
! !INTERFACE:
module SMAP_smobsMod
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
!  21 July 2016: Sujay Kumar, Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_smobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_smobs
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

     real,    allocatable    :: smobs(:,:)
     real,    allocatable    :: smtime(:,:)
     integer*2, allocatable  :: smqc(:,:)     
     logical                 :: startflag     

  end type smapobsdec

  type(smapobsdec), allocatable:: SMAP_smobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMAP_smobsInit
! \label{SMAP_smobsInit}
!
! !INTERFACE: 
  subroutine SMAP_smobsinit(i)
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

    if(.not.allocated(SMAP_smobs)) then 
       allocate(SMAP_smobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_smobs(i)%odir, &
         label='SMAP soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'SMAP soil moisture observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_smobs(i)%data_designation, &
         label='SMAP soil moisture data designation:', rc=status)
    call LVT_verify(status, &
         'SMAP soil moisture data designation: not defined')


    if(SMAP_smobs(i)%data_designation.eq."SPL3SMAP") then 
       !SMAP L3 radar/radiometer daily 9km
       SMAP_smobs(i)%gridDesci=0.0 

       ease_nr=1624
       ease_nc=3856       

    !filling the items needed by the interpolation library
       SMAP_smobs(i)%gridDesci(1) = 9  !input is EASE grid
    !these  corner coordinates were calculated based on ezlh_convert
       SMAP_smobs(i)%gridDesci(4) = -84.6564  !lat
       SMAP_smobs(i)%gridDesci(5) = -179.953 !lon
       SMAP_smobs(i)%gridDesci(7) = 84.6564  !lat
       SMAP_smobs(i)%gridDesci(8) = 179.953  !lon
       SMAP_smobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_smobs(i)%gridDesci(20) = 64
       
       SMAP_smobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_smobs(i)%gridDesci(3) = ease_nr  !ny

       SMAP_smobs(i)%nc=ease_nc
       SMAP_smobs(i)%nr=ease_nr    

       call LVT_update_timestep(LVT_rc, 86400)
   
    elseif(SMAP_smobs(i)%data_designation.eq."SPL3SMP") then 
       !SMAP L3 radiometer daily 36km 
    !filling the items needed by the interpolation library
       ease_nc=964       
       ease_nr=406
       SMAP_smobs(i)%gridDesci = 0
       SMAP_smobs(i)%gridDesci(1) = 9  !input is EASE grid
       SMAP_smobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_smobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_smobs(i)%gridDesci(9) = 4 !M36 grid
       SMAP_smobs(i)%gridDesci(20) = 64
       SMAP_smobs(i)%gridDesci(10) = 0.36
       SMAP_smobs(i)%gridDesci(11) = 1

       SMAP_smobs(i)%nc=ease_nc
       SMAP_smobs(i)%nr=ease_nr       
       
       call LVT_update_timestep(LVT_rc, 86400)
       
    elseif(SMAP_smobs(i)%data_designation.eq."SPL3SMP_E") then 
       ease_nc=3856       
       ease_nr=1624
       SMAP_smobs(i)%gridDesci = 0
       SMAP_smobs(i)%gridDesci(1) = 9  
       SMAP_smobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_smobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_smobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_smobs(i)%gridDesci(20) = 64
       SMAP_smobs(i)%gridDesci(10) = 0.09
       SMAP_smobs(i)%gridDesci(11) = 1

       SMAP_smobs(i)%nc=ease_nc
       SMAP_smobs(i)%nr=ease_nr   

       call LVT_update_timestep(LVT_rc, 86400)

    elseif(SMAP_smobs(i)%data_designation.eq."SPL2SMP") then 
       ease_nc=964       
       ease_nr=406
       SMAP_smobs(i)%gridDesci = 0
       SMAP_smobs(i)%gridDesci(1) = 9  !input is EASE grid
       SMAP_smobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_smobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_smobs(i)%gridDesci(9) = 4 !M36 grid
       SMAP_smobs(i)%gridDesci(20) = 64
       SMAP_smobs(i)%gridDesci(10) = 0.36
       SMAP_smobs(i)%gridDesci(11) = 1

       SMAP_smobs(i)%nc=ease_nc
       SMAP_smobs(i)%nr=ease_nr
       
       call LVT_update_timestep(LVT_rc, 3600)
    
    elseif(SMAP_smobs(i)%data_designation.eq."SPL2SMP_E") then 

       ease_nc=3856       
       ease_nr=1624
       SMAP_smobs(i)%gridDesci = 0
       SMAP_smobs(i)%gridDesci(1) = 9  
       SMAP_smobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_smobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_smobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_smobs(i)%gridDesci(20) = 64
       SMAP_smobs(i)%gridDesci(10) = 0.09
       SMAP_smobs(i)%gridDesci(11) = 1

       SMAP_smobs(i)%nc=ease_nc
       SMAP_smobs(i)%nr=ease_nr   
       
       call LVT_update_timestep(LVT_rc, 3600)

    endif

    npts= LVT_rc%lnc*LVT_rc%lnr
    SMAP_smobs(i)%mo=npts
    
    
    allocate(SMAP_smobs(i)%rlat2(npts))
    allocate(SMAP_smobs(i)%rlon2(npts))
    allocate(SMAP_smobs(i)%n112(npts))
    
    SMAP_smobs(i)%rlat2=0.0
    SMAP_smobs(i)%rlon2=0.0
    SMAP_smobs(i)%n112=0.0
    call neighbor_interp_input(SMAP_smobs(i)%gridDesci,LVT_rc%gridDesc,&
         npts,SMAP_smobs(i)%rlat2,&
         SMAP_smobs(i)%rlon2,SMAP_smobs(i)%n112)
    

    SMAP_smobs(i)%startflag = .true. 

    allocate(SMAP_smobs(i)%smobs(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAP_smobs(i)%smtime(LVT_rc%lnc,LVT_rc%lnr))
    allocate(SMAP_smobs(i)%smqc(LVT_rc%lnc*LVT_rc%lnr,2))

!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------

  end subroutine SMAP_smobsinit


end module SMAP_smobsMod
