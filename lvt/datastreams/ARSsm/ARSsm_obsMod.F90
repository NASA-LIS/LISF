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
! !MODULE: ARSsm_obsMod
! \label(ARSsm_obsMod)
!
! !INTERFACE:
module ARSsm_obsMod
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
!  29 Mar 2012   Sujay Kumar  Initial Specification
! 
!EOP


  PUBLIC :: ARSsm_obsInit
  PUBLIC :: ARSsmobs
  
  type, public :: ismnobsdec

     character*100         :: odir
     integer               :: yr
     integer               :: n_stns
     character*50, allocatable :: stn_name(:)
     real,         allocatable :: stn_lat(:)
     real,         allocatable :: stn_lon(:)
     integer,      allocatable :: stn_col(:)
     integer,      allocatable :: stn_row(:)
     type(ESMF_Time)       :: startTime
     type(ESMF_TimeInterval) :: timestep
     real,         allocatable :: sm(:,:)
     real,         allocatable :: sm_std(:,:)
     real,         allocatable :: rootsm(:,:)
     real,         allocatable :: st(:,:)
  end type ismnobsdec

  type(ismnobsdec), allocatable :: ARSsmobs(:)

contains

!BOP
! 
! !ROUTINE: ARSsm_obsInit
! \label(ARSsm_obsInit)
!
! !INTERFACE:
  subroutine ARSsm_obsInit(i)
! 
! !USES:   
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod,     only : LVT_verify, LVT_logunit, &
         LVT_getNextUnitNumber, LVT_releaseUnitNumber
    use map_utils


    implicit none
!
! !INPUT PARAMETERS: 
    integer,     intent(IN) :: i   ! index of the observation type
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
    integer                 :: k 
    integer                 :: ftn
    integer                 :: status
    real                  :: col,row
    character*100           :: stnlist_file
    type(LVT_metadataEntry),    pointer :: sm
    
    if(.not.allocated(ARSsmobs)) then 
       allocate(ARSsmobs(LVT_rc%nDataStreams))
    endif

    write(LVT_logunit,*) '[INFO] Initializing ARSsm data reader....'
    call ESMF_ConfigGetAttribute(LVT_config, ARSsmobs(i)%odir, &
         label='ARS soil moisture observation directory:',rc=status)
    call LVT_verify(status, 'ARS soil moisture observation directory: not defined')
  
    call ESMF_ConfigGetAttribute(LVT_config, stnlist_file, &
         label='ARS soil moisture station list file:',rc=status)
    call LVT_verify(status, 'ARS soil moisture station list file: not defined')

    ftn = LVT_getNextUnitNumber()

    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) ARSsmobs(i)%n_stns
    read(ftn,*) 

    allocate(ARSsmobs(i)%stn_name(ARSsmobs(i)%n_stns))
    allocate(ARSsmobs(i)%stn_lat(ARSsmobs(i)%n_stns))
    allocate(ARSsmobs(i)%stn_lon(ARSsmobs(i)%n_stns))

    allocate(ARSsmobs(i)%stn_col(ARSsmobs(i)%n_stns))
    allocate(ARSsmobs(i)%stn_row(ARSsmobs(i)%n_stns))
    
    write(LVT_logunit,*) '[INFO] ARSsm station list .. '
    do k=1,ARSsmobs(i)%n_stns
       read(ftn,*) ARSsmobs(i)%stn_name(k), ARSsmobs(i)%stn_lat(k), &
            ARSsmobs(i)%stn_lon(k)
       write(LVT_logunit,*) '[INFO] ',trim(ARSsmobs(i)%stn_name(k)), &
            ARSsmobs(i)%stn_lat(k), ARSsmobs(i)%stn_lon(k)

       call latlon_to_ij(LVT_domain%lvtproj, ARSsmobs(i)%stn_lat(k), &
            ARSsmobs(i)%stn_lon(k), col, row)
       ARSsmobs(i)%stn_col(k) = nint(col)
       ARSsmobs(i)%stn_row(k) = nint(row)
       
    enddo
     
    call LVT_releaseUnitNumber(ftn)

    ARSsmobs(i)%yr = -1

    call ESMF_TimeIntervalSet(ARSsmobs(i)%timestep,s=1800,rc=status)
    call LVT_verify(status,"ESMF_TimeIntervalSet failed in ARSsm_obsInit")

    call LVT_update_timestep(LVT_rc, 1800)

    allocate(ARSsmobs(i)%sm(ARSsmobs(i)%n_stns,366*48))
    allocate(ARSsmobs(i)%sm_std(ARSsmobs(i)%n_stns,366*48))

    ARSsmobs(i)%sm = -9999.0
    ARSsmobs(i)%sm_std = -9999.0

! The ARS data includes reported standard deviation values. For LVT to log them, 
! the data structure must be appropriately initialized

    if(i.eq.1) then 
       sm => LVT_histData%ptr_into_ds1_list(&
            LVT_MOC_SOILMOIST(i))%dataEntryPtr
    else
       sm => LVT_histData%ptr_into_ds2_list(&
            LVT_MOC_SOILMOIST(i))%dataEntryPtr
    endif
! Turning off for testing: SVK
!    sm%stdev_flag  = .true. 

    write(LVT_logunit,*) '[INFO] Finished initializing ARSsm data reader....'

  end subroutine ARSsm_obsInit


end module ARSsm_obsMod
