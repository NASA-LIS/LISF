!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: SMAPEOPLSMobsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! SMAP_E_OPL soil moisture retrievals
!
! !REVISION HISTORY: 
!  26 Apr 2023: Mahdi Navari, Initial Specification
!
module SMAPEOPLSMobsMod
! !USES:
  use ESMF
  use map_utils

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAPEOPLSMobsinit  !Initializes structures for reading SMAP_E_OPL SM data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAPEOPLsmobs !Object to hold SMAPEOPLsm observation attributes
!EOP
  type, public :: smapeoplsmdec

     character*100        :: odir
     integer              :: nc, nr
     real                 :: gridDesci(50)
     real,    allocatable :: smobs(:,:)
     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)
     real,  allocatable   :: rlat(:)
     real,  allocatable   :: rlon(:)

  end type smapeoplsmdec

  type(smapeoplsmdec),allocatable :: SMAPEOPLsmobs(:)

contains

!BOP
! 
! !ROUTINE: SMAPEOPLSMobsinit
! \label{SMAPEOPLSMobsinit}
!
! !INTERFACE: 
  subroutine SMAPEOPLSMobsinit(i)
! 
! !USES:
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i

! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading the SMAP_E_OPL soil moisture data.
! 
!EOP

    integer                 :: status, rc
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: n

    allocate(SMAPEOPLsmobs(LVT_rc%nDataStreams))

    call ESMF_ConfigFindLabel(LVT_config, &
         'SMAP_E_OPL soil moisture data directory:', rc=status)
    call ESMF_ConfigGetAttribute(LVT_Config, SMAPEOPLsmobs(i)%odir, &
         rc=status)
    call LVT_verify(status, &
         'SMAP_E_OPL soil moisture data directory: not defined')
    call LVT_update_timestep(LVT_rc, 3600)

    allocate(SMAPEOPLsmobs(i)%smobs(LVT_rc%lnc,LVT_rc%lnr))

    SMAPEOPLsmobs(i)%smobs = -9999.0

    SMAPEOPLsmobs(i)%nc = 2560
    SMAPEOPLsmobs(i)%nr = 1920

    SMAPEOPLsmobs(i)%gridDesci(1) = 0
    SMAPEOPLsmobs(i)%gridDesci(2) = SMAPEOPLsmobs(i)%nc
    SMAPEOPLsmobs(i)%gridDesci(3) = SMAPEOPLsmobs(i)%nr
    SMAPEOPLsmobs(i)%gridDesci(4) = -89.9531250
    SMAPEOPLsmobs(i)%gridDesci(5) = -179.9296875
    SMAPEOPLsmobs(i)%gridDesci(6) = 128
    SMAPEOPLsmobs(i)%gridDesci(7) = 89.9531250
    SMAPEOPLsmobs(i)%gridDesci(8) = 179.9296875
    SMAPEOPLsmobs(i)%gridDesci(9) = 0.1406250  !dlon
    SMAPEOPLsmobs(i)%gridDesci(10) = 0.0937500 !dlat
    SMAPEOPLsmobs(i)%gridDesci(20) = 64

    if(LVT_isAtAfinerResolution(0.0937500)) then

       allocate(SMAPEOPLsmobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(SMAPEOPLsmobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

       call bilinear_interp_input ( &
            SMAPEOPLsmobs(i)%gridDesci,LVT_rc%gridDesc(:),&
            LVT_rc%lnc*LVT_rc%lnr,&
            SMAPEOPLsmobs(i)%rlat, &
            SMAPEOPLsmobs(i)%rlon, &
            SMAPEOPLsmobs(i)%n11,&
            SMAPEOPLsmobs(i)%n12,&
            SMAPEOPLsmobs(i)%n21,&
            SMAPEOPLsmobs(i)%n22,&
            SMAPEOPLsmobs(i)%w11,&
            SMAPEOPLsmobs(i)%w12,&
            SMAPEOPLsmobs(i)%w21,&
            SMAPEOPLsmobs(i)%w22)

    else

       allocate(SMAPEOPLsmobs(i)%n11(SMAPEOPLsmobs(i)%nc*&
            SMAPEOPLsmobs(i)%nr))

       call upscaleByAveraging_input (&
            SMAPEOPLsmobs(i)%gridDesci,&
            LVT_rc%gridDesc,&
            SMAPEOPLsmobs(i)%nc*&
            SMAPEOPLsmobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,&
            SMAPEOPLsmobs(i)%n11)

    endif

   end subroutine SMAPEOPLSMobsinit

end module SMAPEOPLSMobsMod
