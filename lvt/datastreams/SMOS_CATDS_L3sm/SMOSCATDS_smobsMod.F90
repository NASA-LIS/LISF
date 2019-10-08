!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: SMOSCATDS_smobsMod
! \label(SMOSCATDS_smobsMod)
!
! !INTERFACE:
module SMOSCATDS_smobsMod
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
!  retrieval products from the Soil Moisture Ocean Salinity Mission
!  (SMOS) mission. Specifically, this plugin handles the Level 3
!  CATDS product, which contains the filtered data processed
!  from the ESA L1B product. 
!  
!  website: http://www.catds.fr/Products/Available-products-from-CPDC
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  28 Jan 2017: Sujay Kumar, Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSCATDS_smobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSCATDS_smobs
!EOP
  type, public :: smoscatdsobsdec

     character*100        :: odir

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

  end type smoscatdsobsdec

  type(smoscatdsobsdec), allocatable:: SMOSCATDS_smobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOSCATDS_smobsInit
! \label{SMOSCATDS_smobsInit}
!
! !INTERFACE: 
  subroutine SMOSCATDS_smobsinit(i)
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

    if(.not.allocated(SMOSCATDS_smobs)) then 
       allocate(SMOSCATDS_smobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMOSCATDS_smobs(i)%odir, &
         label='SMOS (CATDS) soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'SMOS (CATDS) soil moisture observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    
    SMOSCATDS_smobs(i)%gridDesci(1) = 9  !input is EASE grid
!    ease_nr=586
!    ease_nc=1383
    ease_nr=584
    ease_nc=1388
    
    SMOSCATDS_smobs(i)%gridDesci=0.0 
    
    SMOSCATDS_smobs(i)%gridDesci(1) = 9  !input is EASE grid

    SMOSCATDS_smobs(i)%gridDesci(4) = -83.51714
    SMOSCATDS_smobs(i)%gridDesci(5) = -179.8703
    SMOSCATDS_smobs(i)%gridDesci(7) = 83.51714
    SMOSCATDS_smobs(i)%gridDesci(8) = 179.8703

    SMOSCATDS_smobs(i)%gridDesci(9) = 1
    SMOSCATDS_smobs(i)%gridDesci(10) = 0.25 ! MN: based on Kristi's email 1/16/2019    
    SMOSCATDS_smobs(i)%gridDesci(20) = 64
    
    SMOSCATDS_smobs(i)%gridDesci(2) = ease_nc  !nx
    SMOSCATDS_smobs(i)%gridDesci(3) = ease_nr  !ny
    
    SMOSCATDS_smobs(i)%nc=ease_nc
    SMOSCATDS_smobs(i)%nr=ease_nr       

    npts= LVT_rc%lnc*LVT_rc%lnr

!    npts =1440*720

    SMOSCATDS_smobs(i)%mo=npts
    
    
    allocate(SMOSCATDS_smobs(i)%rlat2(npts))
    allocate(SMOSCATDS_smobs(i)%rlon2(npts))
    allocate(SMOSCATDS_smobs(i)%n112(npts))
    
    SMOSCATDS_smobs(i)%rlat2=0.0
    SMOSCATDS_smobs(i)%rlon2=0.0
    SMOSCATDS_smobs(i)%n112=0.0

 !   LVT_rc%gridDesc(3) = 720
 !   LVT_rc%gridDesc(4) = -89.875

    call neighbor_interp_input(SMOSCATDS_smobs(i)%gridDesci,LVT_rc%gridDesc,&
         npts,SMOSCATDS_smobs(i)%rlat2,&
         SMOSCATDS_smobs(i)%rlon2,SMOSCATDS_smobs(i)%n112)
    

    SMOSCATDS_smobs(i)%startflag = .true. 

    allocate(SMOSCATDS_smobs(i)%smobs(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMOSCATDS_smobs(i)%smtime(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMOSCATDS_smobs(i)%smqc(LVT_rc%lnc*LVT_rc%lnr,2))

  end subroutine SMOSCATDS_smobsinit


end module SMOSCATDS_smobsMod
