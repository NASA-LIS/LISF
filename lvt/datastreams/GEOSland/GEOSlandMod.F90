!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: GEOSlandMod
! \label(GEOSlandMod)
!
! !INTERFACE:
module GEOSlandMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the use of a GEOS land outputs 
!  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  20 Oct 2021    Sujay Kumar  Initial Specification
! 
!EOP

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GEOSlandInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: geoslandoutput
!
!EOP
  
  type, public :: geoslandobsdec
     character*100  :: odir
  end type geoslandobsdec

  type(geoslandobsdec), allocatable :: geoslandoutput(:)

contains

!BOP
! 
! !ROUTINE: GEOSlandInit
! \label{GEOSlandInit}
!
! !INTERFACE: 
  subroutine GEOSlandInit(i)
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
    integer,     intent(IN) :: i   ! index of the observation type
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of
! GEOS land outputs    
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: t
    integer                 :: ts
    integer                 :: rc

    if(.not.allocated(geoslandoutput)) then 
       allocate(geoslandoutput(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config,geoslandoutput(i)%odir, &
         label="GEOS land output directory:",rc=rc)
    call LVT_verify(rc,'GEOS land output directory: not defined')

    ts = 86400 !daily
    call LVT_update_timestep(LVT_rc, ts)

  end subroutine GEOSlandInit
  
end module GEOSlandMod
