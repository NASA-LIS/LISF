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
! !MODULE: GIMMSAVHRR_NDVIobsMod
! \label(GIMMSAVHRR_NDVIobsMod)
!
! !INTERFACE:
module GIMMSAVHRR_NDVIobsMod
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
  PUBLIC :: GIMMSAVHRR_NDVIobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GIMMSAVHRRNDVIobs !Object to hold GIMMSAVHRRNDVI observation attributes
!EOP

  type, public :: gimmsavhrrndvidec
     character*100           :: odir
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: yr
     integer                 :: mo
     logical                 :: startFlag
     real                    :: datares
  end type gimmsavhrrndvidec
     
  type(gimmsavhrrndvidec), save :: GIMMSAVHRRNDVIObs(2)

contains
  
!BOP
! 
! !ROUTINE: GIMMSAVHRR_NDVIobsInit
! \label{GIMMSAVHRR_NDVIobsInit}
!
! !INTERFACE: 
  subroutine GIMMSAVHRR_NDVIobsinit(i)
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
!   for reading the GIMMSAVHRR NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    call ESMF_ConfigGetAttribute(LVT_Config, GIMMSAVHRRNDVIobs(i)%odir, &
         label='GIMMS AVHRR NDVI data directory:',rc=status)
    call LVT_verify(status, 'GIMMS AVHRR NDVI data directory: not defined')

    call LVT_update_timestep(LVT_rc, 2592000)

    gimmsavhrrndviobs(i)%gridDesc = 0
        
    gimmsavhrrndviobs(i)%nc = 4320
    gimmsavhrrndviobs(i)%nr = 2160

    !filling the items needed by the interpolation library
    gimmsavhrrndviobs(i)%gridDesc(1) = 0  
    gimmsavhrrndviobs(i)%gridDesc(2) = gimmsavhrrndviobs(i)%nc
    gimmsavhrrndviobs(i)%gridDesc(3) = gimmsavhrrndviobs(i)%nr
    gimmsavhrrndviobs(i)%gridDesc(4) = -89.9583 ! dx/2,dy/2 = 1/24 = 0.04167
    gimmsavhrrndviobs(i)%gridDesc(5) = -179.9583
    gimmsavhrrndviobs(i)%gridDesc(7) = 89.9583
    gimmsavhrrndviobs(i)%gridDesc(8) = 179.9583
    gimmsavhrrndviobs(i)%gridDesc(6) = 128
    gimmsavhrrndviobs(i)%gridDesc(9) = 0.0833
    gimmsavhrrndviobs(i)%gridDesc(10) = 0.0833
    gimmsavhrrndviobs(i)%gridDesc(20) = 64

    gimmsavhrrndviobs(i)%datares  = 0.0833

    if(LVT_isAtAfinerResolution(gimmsavhrrndviobs(i)%datares)) then
       
       allocate(gimmsavhrrndviobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gimmsavhrrndviobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gimmsavhrrndviobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(gimmsavhrrndviobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            gimmsavhrrndviobs(i)%rlat, &
            gimmsavhrrndviobs(i)%rlon, &
            gimmsavhrrndviobs(i)%n11)
    else
       allocate(gimmsavhrrndviobs(i)%n11(gimmsavhrrndviobs(i)%nc*gimmsavhrrndviobs(i)%nr))
       call upscaleByAveraging_input(gimmsavhrrndviobs(i)%gridDesc,&
            LVT_rc%gridDesc,gimmsavhrrndviobs(i)%nc*gimmsavhrrndviobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,gimmsavhrrndviobs(i)%n11)
    endif

    gimmsavhrrndviobs(i)%yr = -1
    gimmsavhrrndviobs(i)%mo = LVT_rc%mo
    gimmsavhrrndviobs(i)%startFlag = .false. 

  end subroutine GIMMSAVHRR_NDVIobsinit


end module GIMMSAVHRR_NDVIobsMod
