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
! !MODULE: USGSSF_obsMod
! \label(USGSSF_obsMod)
!
! !INTERFACE:
module USGSSF_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: USGSSF_obsinit
  PUBLIC :: USGSSFobs

  type, public :: USGSSFobsdec
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
!  13 May 2011   Sujay Kumar  Initial Specification
! 
!EOP

     character*100           :: odir
     integer                 :: n_stns
     integer                 :: nts
     integer                 :: version
     character*100, allocatable  :: stn_name(:)
     real,          allocatable  :: stnlat(:)
     real,          allocatable  :: stnlon(:)
     integer                 :: yr
     real,          allocatable  :: q(:,:)
     type(ESMF_Time)         :: startTime
     type(ESMF_TimeInterval) :: timestep
  end type USGSSFobsdec

  type(USGSSFobsdec), allocatable :: USGSSFobs(:)

contains
  
!BOP
! 
! !ROUTINE: USGSSF_obsInit
! \label{USGSSF_obsInit}
!
! !INTERFACE: 
 subroutine USGSSF_obsinit(i)
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    character*100 :: stnlist_file
    integer       :: ftn, k, status

    if(.not.allocated(USGSSFobs)) then 
       allocate(USGSSFobs(LVT_rc%nDataStreams))
    endif
!------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------

    write(LVT_logunit,*) '[INFO] Initializing USGS streamflow data reader....'
    call ESMF_ConfigGetAttribute(LVT_config, USGSSFobs(i)%odir, &
         label='USGS streamflow observation directory:',rc=status)
    call LVT_verify(status, 'USGS streamflow observation directory: not defined')
    call ESMF_ConfigGetAttribute(LVT_config, USGSSFobs(i)%version, &
         label='USGS streamflow observation version:',rc=status)
    call LVT_verify(status, 'USGS streamflow observation version: not defined')
  
    call ESMF_ConfigGetAttribute(LVT_config, stnlist_file, &
         label='USGS streamflow station list file:',rc=status)
    call LVT_verify(status, 'USGS streamflow station list file: not defined')

    ftn = LVT_getNextUnitNumber()

    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) USGSSFobs(i)%n_stns
    read(ftn,*) 

    allocate(USGSSFobs(i)%stn_name(USGSSFobs(i)%n_stns))
    allocate(USGSSFobs(i)%stnlat(USGSSFobs(i)%n_stns))
    allocate(USGSSFobs(i)%stnlon(USGSSFobs(i)%n_stns))

!------------------------------------------------------------------------------
! For each station, this reads the site name, station name, and station
! position in lat, lon coordinates.
!------------------------------------------------------------------------------
    do k=1,USGSSFobs(i)%n_stns
       read(ftn,*) USGSSFobs(i)%stn_name(k), USGSSFobs(i)%stnlat(k), &
            USGSSFobs(i)%stnlon(k)
    end do
    call LVT_releaseUnitNumber(ftn)

    if(USGSSFobs(i)%version.eq.1) then 
       call ESMF_TimeIntervalSet(USGSSFobs(i)%timestep,s=86400,rc=status)
       call LVT_verify(status,"ESMF_TimeIntervalSet failed in USGSSF_obsInit")

       USGSSFobs(i)%nts = 366
       call LVT_update_timestep(LVT_rc, 86400)
    else
       call ESMF_TimeIntervalSet(USGSSFobs(i)%timestep,s=900,rc=status)
       call LVT_verify(status,"ESMF_TimeIntervalSet failed in USGSSF_obsInit")

       USGSSFobs(i)%nts = 35136 ! 366x24x4
       call LVT_update_timestep(LVT_rc, 900)
    endif

    allocate(USGSSFobs(i)%q(USGSSFobs(i)%n_stns,USGSSFobs(i)%nts))
    USGSSFobs(i)%q = -9999.0

    USGSSFobs(i)%yr = -1

  end subroutine USGSSF_obsinit


end module USGSSF_obsMod
