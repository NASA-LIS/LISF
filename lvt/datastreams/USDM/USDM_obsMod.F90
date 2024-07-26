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
! !MODULE: USDM_obsMod
! \label(USDM_obsMod)
!
! !INTERFACE:
module USDM_obsMod
! 
! !USES:   
  use ESMF
  use map_utils

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
!  8 Mar 2017   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: USDM_obsinit !Initializes structures for reading USDM data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: USDMobs !Object to hold USDM observation attributes
!EOP

  type, public :: usdmdec
     character*100           :: odir
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
     real,    allocatable    :: drcategory(:)
     real                    :: gridDesc(50)
     logical                 :: startFlag
     type(proj_info)         :: map_proj
  end type usdmdec
     
  type(usdmdec), allocatable :: USDMObs(:)

contains
  
!BOP
! 
! !ROUTINE: USDM_obsInit
! \label{USDM_obsInit}
!
! !INTERFACE: 
  subroutine USDM_obsinit(i)
! 
! !USES: 
    use LVT_coreMod,   only : LVT_rc, LVT_Config
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
!   for reading the USDM data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status
    real                  :: cornerlat1, cornerlat2
    real                  :: cornerlon1, cornerlon2

    if(.not.allocated(USDMobs)) then 
       allocate(USDMobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, USDMObs(i)%odir, &
         label='USDM data directory: ',rc=status)
    call LVT_verify(status, 'USDM data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    allocate(USDMobs(i)%drcategory(LVT_rc%lnc*LVT_rc%lnr))

    usdmobs(i)%gridDesc = 0
    
    usdmobs(i)%startFlag = .true.

  end subroutine USDM_obsinit


end module USDM_obsMod
