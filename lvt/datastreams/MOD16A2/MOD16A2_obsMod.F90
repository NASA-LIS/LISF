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
! !MODULE: MOD16A2_obsMod
! \label(MOD16A2_obsMod)
!
! !INTERFACE:
module MOD16A2_obsMod
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
!  18 May 2011   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MOD16A2_obsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MOD16A2obs !Object to hold MOD16A2 observation attributes
!EOP

  type, public :: mod16a2dec
     character*100           :: odir
     integer                 :: nc, nr
     integer, allocatable        :: n11(:)
     real,    allocatable        :: qle(:)
     integer                 :: yr
     integer                 :: mo
     real                    :: gridDesc(50)
     logical                 :: startFlag
  end type mod16a2dec
     
  type(mod16a2dec), allocatable :: MOD16A2Obs(:)

contains
  
!BOP
! 
! !ROUTINE: MOD16A2_obsInit
! \label{MOD16A2_obsInit}
!
! !INTERFACE: 
  subroutine MOD16A2_obsinit(i)
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
!   for reading the MOD16A2 data, including the computation of spatial 
!   interpolation weights. The MOD16A2 data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status
    real                  :: cornerlat1, cornerlat2
    real                  :: cornerlon1, cornerlon2

    if(.not.allocated(MOD16A2obs)) then 
       allocate(MOD16A2obs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, MOD16A2Obs(i)%odir, &
         label='MOD16A2 data directory: ',rc=status)
    call LVT_verify(status, 'MOD16A2 data directory: not defined')

    call LVT_update_timestep(LVT_rc, 2592000)

    allocate(MOD16A2obs(i)%qle(LVT_rc%lnc*LVT_rc%lnr))

    mod16a2obs(i)%gridDesc = 0
    
    cornerlat1 = max(-59.995,nint((LVT_rc%gridDesc(4)+59.995)/0.01)*0.01-59.995-2*0.01)
    cornerlon1 = max(-179.995,nint((LVt_rc%gridDesc(5)+179.995)/0.01)*0.01-179.995-2*0.01)
    cornerlat2 = min(89.995,nint((LVT_rc%gridDesc(7)+59.995)/0.01)*0.01-59.995+2*0.01)
    cornerlon2 = min(179.995,nint((LVT_rc%gridDesc(8)+179.995)/0.01)*0.01-179.995+2*0.01)
    
    mod16a2obs(i)%nr = nint((cornerlat2-cornerlat1)/0.01)+1
    mod16a2obs(i)%nc = nint((cornerlon2-cornerlon1)/0.01)+1

    allocate(MOD16A2Obs(i)%n11(mod16a2obs(i)%nc*mod16a2obs(i)%nr))

    !filling the items needed by the interpolation library
    mod16a2obs(i)%gridDesc(1) = 0  !input is EASE grid
    mod16a2obs(i)%gridDesc(2) = mod16a2obs(i)%nc
    mod16a2obs(i)%gridDesc(3) = mod16a2obs(i)%nr
    mod16a2obs(i)%gridDesc(4) = cornerlat1
    mod16a2obs(i)%gridDesc(5) = cornerlon1
    mod16a2obs(i)%gridDesc(7) = cornerlat2
    mod16a2obs(i)%gridDesc(8) = cornerlon2
    mod16a2obs(i)%gridDesc(6) = 128
    mod16a2obs(i)%gridDesc(9) = 0.01
    mod16a2obs(i)%gridDesc(10) = 0.01
    mod16a2obs(i)%gridDesc(20) = 64

    call upscaleByAveraging_input(mod16a2obs(i)%gridDesc,&
         LVT_rc%gridDesc,mod16a2obs(i)%nc*mod16a2obs(i)%nr,&
         LVT_rc%lnc*LVT_rc%lnr,mod16a2obs(i)%n11)

    MOD16A2obs(i)%mo = LVT_rc%mo
    MOD16A2obs(i)%yr = -1
    mod16a2obs(i)%startFlag = .true.

    if(LVT_rc%tavgInterval.lt.2592000) then 
       write(LVT_logunit,*) '[ERR] The time averaging interval must be greater than'
       write(LVT_logunit,*) '[ERR] equal to a month since the MOD16A2 data is monthly'
       call LVT_endrun()
    endif

  end subroutine MOD16A2_obsinit


end module MOD16A2_obsMod
