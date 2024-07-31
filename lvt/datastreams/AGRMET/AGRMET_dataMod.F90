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
! !MODULE: AGRMET_dataMod
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This subroutine provides the observation plugin for reading the 
!   operational AGRMET output from the Air Force Weather Agency (AFWA)
!   This plugin only handles the LIS-style outputs and not the old
!   AGRMET grib files at 1/2 deg.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 Dec 2010   Sujay Kumar  Initial Specification
!  02 Nov 2018   Eric Kemp    Added support for n1280e configuration
! 
!EOP
!
module AGRMET_dataMod
! !USES: 
  use ESMF

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: AGRMET_datainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: AGRMETdata
!EOP
  type, public :: agrmetdatadec
     character*255           :: odir
     real*8                  :: changetime1,changetime2
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
     integer, allocatable        :: n11(:)
     integer, allocatable        :: n12(:)
     integer, allocatable        :: n21(:)
     integer, allocatable        :: n22(:)     
     real,    allocatable        :: w11(:)
     real,    allocatable        :: w12(:)
     real,    allocatable        :: w21(:)
     real,    allocatable        :: w22(:)

     character*20           :: security_class
     character*20           :: distribution_class
     character*20           :: data_category
     character*20           :: area_of_data

     integer                 :: nc
     integer                 :: nr
     type(ESMF_TimeInterval) :: ts

     character*20           :: gridname
  end type agrmetdatadec

  type(agrmetdatadec),allocatable:: agrmetdata(:)

contains
  
!BOP
! 
! !ROUTINE: AGRMET_dataInit
! \label{AGRMET_dataInit}
!
! !INTERFACE: 
  subroutine AGRMET_datainit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_logMod
    use LVT_histDataMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This subroutine initializes and sets up the data structures required
!   for reading the AGRMET data, including the setup of spatial interpolation
!   weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer,   intent(IN) :: i 
!EOP
    integer              :: status
    real                 :: gridDesci(50)
    integer              :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real                 :: upgmt
    character*10         :: time
    integer              :: ts

    if(.not.allocated(agrmetdata)) then 
       allocate(agrmetdata(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, agrmetdata(i)%odir, &
         label='AGRMET data directory:', rc=status)
    call LVT_verify(status, 'AGRMET data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, agrmetdata(i)%security_class, &
         label='AGRMET data security class name:', rc=status)
    call LVT_verify(status, 'AGRMET data security class name: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, agrmetdata(i)%distribution_class, &
         label='AGRMET data distribution class name:', rc=status)
    call LVT_verify(status, 'AGRMET data distribution class name: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, agrmetdata(i)%data_category, &
         label='AGRMET data category name:', rc=status)
    call LVT_verify(status, 'AGRMET data category name: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, agrmetdata(i)%area_of_data, &
         label='AGRMET data area of data:', rc=status)
    call LVT_verify(status, 'AGRMET data area of data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, agrmetdata(i)%gridname, &
         label='AGRMET data gridname:', rc=status)
    call LVT_verify(status, 'AGRMET data gridname: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, time, &
         label='AGRMET data output interval:', rc=status)
    call LVT_verify(status, 'AGRMET data output interval: not defined')

    call LVT_parseTimeString(time,ts)

    call LVT_update_timestep(LVT_rc,ts)

    allocate(agrmetdata(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(agrmetdata(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(agrmetdata(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    ! Sanity check the grid type
    ! NOTE:  557WW refers to the 0.25 deg deterministic product as "GLOBAL"
    ! even though it excludes Antarctica
    if (trim(agrmetdata(i)%gridname) == "GLOBAL") then
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 1440
       gridDesci(3) = 600
       gridDesci(4) = -59.875
       gridDesci(5) = -179.875
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(6) = 128
       gridDesci(9) = 0.250
       gridDesci(10) = 0.250
       gridDesci(20) = 64

       agrmetdata(i)%nc = 1440
       agrmetdata(i)%nr =  600
       
    else if (trim(agrmetdata(i)%gridname) == "n1280e") then
       ! NOTE:  This is the in-house Bratseth reanalysis produced at rough
       ! 0.09 deg resolution.  The name is taken from the global configuration
       ! for GALWEM.
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 2560
       gridDesci(3) = 1920
       gridDesci(4) =  -89.9531250
       gridDesci(5) = -179.9296875
       gridDesci(7) =   89.9531250
       gridDesci(8) =  179.9296875
       gridDesci(6) = 128
       gridDesci(9) =  0.140625 ! dlon
       gridDesci(10) = 0.093750 ! dlat
       gridDesci(20) = 64

       agrmetdata(i)%nc = 2560
       agrmetdata(i)%nr = 1920
    else
       write(LVT_logunit,*)'[ERR] Invalid AGRMET grid type specified!'
       write(LVT_logunit,*)'Currently supports GLOBAL and n1280e'
       write(LVT_logunit,*)'Found ',trim(agrmetdata(i)%gridname)
       call LVT_endrun()
    end if

    write(LVT_logunit,*)'EMK: AGRMET calling bilinear_interp_input...'
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         agrmetdata(i)%rlat, agrmetdata(i)%rlon,&
         agrmetdata(i)%n11, agrmetdata(i)%n12, &
         agrmetdata(i)%n21, agrmetdata(i)%n22, & 
         agrmetdata(i)%w11, agrmetdata(i)%w12, &
         agrmetdata(i)%w21, agrmetdata(i)%w22)

    call ESMF_TimeIntervalSet(agrmetdata(i)%ts, s = 10800, &
         rc=status)
    call LVT_verify(status)

  end subroutine AGRMET_datainit


end module AGRMET_dataMod
