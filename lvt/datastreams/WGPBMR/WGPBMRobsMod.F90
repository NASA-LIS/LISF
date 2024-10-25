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
! !MODULE: WGPBMRobsMod
! \label(WGPBMRobsMod)
!
! !INTERFACE:
module WGPBMRobsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: WGPBMRobsinit
  PUBLIC :: WGPBMRobs
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
! 
!EOP

  type, public :: wgpbmrobsdec
     character*100        :: odir
     integer              :: stnid
  end type wgpbmrobsdec

  type(wgpbmrobsdec), allocatable :: wgpbmrobs(:)

contains
  
!BOP
! 
! !ROUTINE: WGPBMRobsInit
! \label{WGPBMRobsInit}
!
! !INTERFACE: 
  subroutine WGPBMRobsinit(i)
! 
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_config
    use LVT_histDataMod
    use LVT_timeMgrMod, only : LVT_calendar
    use LVT_logMod,     only : LVT_verify, LVT_logunit, &
         LVT_getNextUnitNumber, LVT_releaseUnitNumber

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading walnut gulch PBMR data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: status, rc

    if(.not.allocated(wgpbmrobs)) then 
       allocate(wgpbmrobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, wgpbmrobs(i)%odir, &
         label='WG PBMR observation directory:', rc=status)
    call LVT_verify(status, 'WG PBMR observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, wgpbmrobs(i)%stnid, &
         label='WG PBMR site index:', rc=status)
    call LVT_verify(status, 'WG PBMR site index: not defined')


  end subroutine WGPBMRobsinit


end module WGPBMRobsMod
