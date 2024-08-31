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
! !MODULE: MODIS_LSTobsMod
! \label(MODIS_LSTobsMod)
!
! !INTERFACE:
module MODIS_LSTobsMod
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
  PUBLIC :: MODIS_LSTobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MODISLSTobs !Object to hold MODISLST observation attributes
!EOP

  type, public :: modislstdec
     character*100           :: odir
     integer                 :: nc, nr
     integer, allocatable        :: n11(:)
     real,    allocatable        :: lst(:)
     real                    :: gridDesc(50)
  end type modislstdec
     
  type(modislstdec), allocatable :: MODISLSTObs(:)

contains
  
!BOP
! 
! !ROUTINE: MODIS_LSTobsInit
! \label{MODIS_LSTobsInit}
!
! !INTERFACE: 
  subroutine MODIS_LSTobsinit(i)
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
!   for reading the MODIS LST data, including the computation of spatial 
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

    if(.not.allocated(MODISLSTobs)) then 
       allocate(MODISLSTobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, MODISLSTobs(i)%odir, &
         label='MODIS LST data directory:',rc=status)
    call LVT_verify(status, 'MODIS LST data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    allocate(Modislstobs(I)%lst(LVT_rc%lnc*LVT_rc%lnr))

    modislstobs(i)%gridDesc = 0
    
    cornerlat1 = max(-59.995,nint((LVT_rc%gridDesc(4)+59.995)/0.01)*0.01-59.995-2*0.01)
    cornerlon1 = max(-179.995,nint((LVt_rc%gridDesc(5)+179.995)/0.01)*0.01-179.995-2*0.01)
    cornerlat2 = min(89.995,nint((LVT_rc%gridDesc(7)+59.995)/0.01)*0.01-59.995+2*0.01)
    cornerlon2 = min(179.995,nint((LVT_rc%gridDesc(8)+179.995)/0.01)*0.01-179.995+2*0.01)
    
    modislstobs(i)%nr = nint((cornerlat2-cornerlat1)/0.01)+1
    modislstobs(i)%nc = nint((cornerlon2-cornerlon1)/0.01)+1

    allocate(Modislstobs(i)%n11(modislstobs(i)%nc*modislstobs(i)%nr))

    !filling the items needed by the interpolation library
    modislstobs(i)%gridDesc(1) = 0  !input is EASE grid
    modislstobs(i)%gridDesc(2) = modislstobs(i)%nc
    modislstobs(i)%gridDesc(3) = modislstobs(i)%nr
    modislstobs(i)%gridDesc(4) = cornerlat1
    modislstobs(i)%gridDesc(5) = cornerlon1
    modislstobs(i)%gridDesc(7) = cornerlat2
    modislstobs(i)%gridDesc(8) = cornerlon2
    modislstobs(i)%gridDesc(6) = 128
    modislstobs(i)%gridDesc(9) = 0.01
    modislstobs(i)%gridDesc(10) = 0.01
    modislstobs(i)%gridDesc(20) = 64

    call upscaleByAveraging_input(modislstobs(i)%gridDesc,&
         LVT_rc%gridDesc,modislstobs(i)%nc*modislstobs(i)%nr,&
         LVT_rc%lnc*LVT_rc%lnr,modislstobs(i)%n11)

  end subroutine MODIS_LSTobsinit


end module MODIS_LSTobsMod
