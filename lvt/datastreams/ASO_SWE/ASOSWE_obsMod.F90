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
! !MODULE: ASOSWE_obsMod
! \label(ASOSWE_obsMod)
!
! !INTERFACE:
module ASOSWE_obsMod
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
!  This module handles the observation plugin for the NASA Airborne
!  Snow Observatory (ASO) snow water equivalent (SWE) data. 
! 
!   Website: https://aso.jpl.nasa.gov
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Aug 2018   Sujay Kumar  Initial Specification
! 
!EOP
! !USES: 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ASOSWE_obsinit !Initializes structures for reading ASOSWE data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ASOSWEobs !Object to hold ASOSWE observation attributes
!EOP
  type, public :: asosweobsdec
     character*100           :: odir
  end type asosweobsdec

  type(asosweobsdec), allocatable :: asosweobs(:)

contains
  
!BOP
! 
! !ROUTINE: ASOSWE_obsInit
! \label{ASOSWE_obsInit}
!
! !INTERFACE:
  subroutine ASOSWE_obsinit(i)
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
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading ASOSWE data. The ASOSWE data is provides in the WGS 1984 grid
!  with 30 arc second resolution. The domain extents are from (24.9504, -124.7337) 
!  to (52.8754, -66.9421). The data entries are 16-bit signed integers. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer            :: status, rc
    integer            :: ftn, k
    real*8             :: tdur
    integer            :: syr, smo, sda, shr, smn, sss
    integer            :: eyr, emo, eda, ehr, emn, ess
    integer            :: ts
    character*100      :: coordfile
    character*100      :: mdata
    real               :: xi1,xj1,xmesh,orient,alat1,alon1
    integer            :: t
    real               :: gridDesci(50)

    if(.not.allocated(asosweobs)) then 
       allocate(asosweobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, asosweobs(i)%odir, &
         label='ASO SWE observation directory:',rc=status)
    call LVT_verify(status, 'ASO SWE observation directory: not defined')


    ts = 86400
    call LVT_update_timestep(LVT_rc, 86400)

  end subroutine ASOSWE_obsinit


end module ASOSWE_obsMod
