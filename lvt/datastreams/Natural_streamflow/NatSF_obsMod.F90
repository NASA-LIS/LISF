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
! !MODULE: NatSF_obsMod
! \label(NatSF_obsMod)
!
! !INTERFACE:
module NatSF_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: NatSF_obsinit
  PUBLIC :: NatSFobs

  type, public :: NatSFobsdec
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
     character*100, allocatable  :: stn_name(:)
     real,          allocatable  :: stnlat(:)
     real,          allocatable  :: stnlon(:)
     integer                 :: yr,mo
     logical                 :: startFlag
     real,          allocatable  :: q(:)
  end type NatSFobsdec

  type(NatSFobsdec), allocatable :: NatSFobs(:)

contains
  
!BOP
! 
! !ROUTINE: NatSF_obsInit
! \label{NatSF_obsInit}
!
! !INTERFACE: 
 subroutine NatSF_obsinit(i)
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

    if(.not.allocated(NatSFobs)) then 
       allocate(NatSFobs(LVT_rc%nDataStreams))
    endif

!------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------

    write(LVT_logunit,*) '[INFO] Initializing Naturalized streamflow data reader....'
    call ESMF_ConfigGetAttribute(LVT_config, NatSFobs(i)%odir, &
         label='Naturalized streamflow observation directory:',rc=status)
    call LVT_verify(status, 'Naturalized streamflow observation directory: not defined')
  
    call ESMF_ConfigGetAttribute(LVT_config, stnlist_file, &
         label='Naturalized streamflow station list file:',rc=status)
    call LVT_verify(status, 'Naturalized streamflow station list file: not defined')

    ftn = LVT_getNextUnitNumber()

    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) NatSFobs(i)%n_stns
    read(ftn,*) 

    allocate(NatSFobs(i)%stn_name(NatSFobs(i)%n_stns))
    allocate(NatSFobs(i)%stnlat(NatSFobs(i)%n_stns))
    allocate(NatSFobs(i)%stnlon(NatSFobs(i)%n_stns))

!------------------------------------------------------------------------------
! For each station, this reads the site name, station name, and station
! position in lat, lon coordinates.
!------------------------------------------------------------------------------
    do k=1,NatSFobs(i)%n_stns
       read(ftn,*) NatSFobs(i)%stn_name(k), NatSFobs(i)%stnlat(k), &
            NatSFobs(i)%stnlon(k)
!       print*, NatSFobs(i)%stn_name(k), NatSFobs(i)%stnlat(k), &
!            NatSFobs(i)%stnlon(k)
    end do
    call LVT_releaseUnitNumber(ftn)
    
    allocate(NatSFobs(i)%q(NatSFobs(i)%n_stns))
    NatSFobs(i)%q = -9999.0

    NatSFobs(i)%yr = -1
    NatSFobs(i)%mo = -1
    NatSFobs(i)%startFlag = .true.
    call LVT_update_timestep(LVT_rc, 2592000)

  end subroutine NatSF_obsinit
end module NatSF_obsMod
