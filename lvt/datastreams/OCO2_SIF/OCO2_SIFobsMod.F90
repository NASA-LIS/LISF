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
! !MODULE: OCO2_SIFobsMod
! \label(OCO2_SIFobsMod)
!
! !INTERFACE:
module OCO2_SIFobsMod
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
!  This plugin supports the processing of Solar Induced Fluorescence (SIF)
!  data from OCO-2. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  22 Apr 2018   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: OCO2_SIFobsinit !Initializes structures for reading OCO2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: OCO2SIFobs !Object to hold OCO2SIF observation attributes
!EOP

  type, public :: oco2sifdec
     character*100           :: odir
     character*20            :: channel
     logical                 :: startFlag
  end type oco2sifdec
     
  type(oco2sifdec), allocatable :: OCO2SIFObs(:)

contains
  
!BOP
! 
! !ROUTINE: OCO2_SIFobsInit
! \label{OCO2_SIFobsInit}
!
! !INTERFACE: 
  subroutine OCO2_SIFobsinit(i)
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
!   for reading the OCO2 data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    if(.not.allocated(OCO2SIFobs)) then 
       allocate(OCO2SIFobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, OCO2SIFobs(i)%odir, &
         label='OCO2 SIF data directory:',rc=status)
    call LVT_verify(status, 'OCO2 SIF data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, OCO2SIFobs(i)%channel, &
         label='OCO2 SIF data channel frequency:',rc=status)
    if(status.ne.0) then 
       write(LVT_logunit,*) '[ERR] OCO2 SIF data channel frequency: not defined'
       write(LVT_logunit,*) '[ERR] The options are : '
       write(LVT_logunit,*) '[ERR] "757nm" or "771nm"'
       call LVT_endrun()
    endif

    call LVT_update_timestep(LVT_rc, 86400)

    oco2sifobs(i)%startFlag = .true. 

  end subroutine OCO2_SIFobsinit


end module OCO2_SIFobsMod
