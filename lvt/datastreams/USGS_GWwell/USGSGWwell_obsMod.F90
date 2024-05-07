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
! !MODULE: USGSGWwell_obsMod
! \label(USGSGWwell_obsMod)
!
! !INTERFACE:
module USGSGWwell_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: USGSGWwell_obsinit
  PUBLIC :: USGSGWwellobs
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for
!  the USGSGW depth-to-water well observations.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Aug 2009   Sujay Kumar  Initial Specification
!  10 Dec 2015   David Mocko  Added depth to water
!
!EOP

  type, public :: usgsgwwellobsdec
     character*100        :: odir
     integer              :: nstns
     integer              :: nstates
     real                 :: udef
     integer              :: nts
     character*30,allocatable :: stnid(:)
     real,        allocatable :: stnlat(:)
     real,        allocatable :: stnlon(:)
     real,        allocatable :: stnqa(:)
     real,        allocatable :: stnsy(:)
     real,        allocatable :: gw(:,:)
     real,        allocatable :: wtdepth(:,:)
     integer              :: yr
     type(ESMF_TimeInterval) :: timestep
     type(ESMF_Time)         :: startTime

  end type usgsgwwellobsdec

  type(usgsgwwellobsdec), save :: usgsgwwellobs(2)

contains
  
!BOP
! 
! !ROUTINE: USGSGWwell_obsInit
! \label{USGSGWwell_obsInit}
!
! !INTERFACE: 
  subroutine USGSGWwell_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in) :: i
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading USGSGWwell data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer            :: status, rc
    integer            :: ftn,k
    character*100      :: coordfile
    
    call ESMF_ConfigGetAttribute(LVT_config, usgsgwwellobs(i)%odir, &
         label='USGS ground water (well data) observation directory:',rc=status)
    call LVT_verify(status, 'USGS ground water (well data) observation directory: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, coordfile, &
         label='USGS ground water (well data) coord file:',rc=status)
    call LVT_verify(status, 'USGS ground water (well data) coord file: not defined')

    ftn=LVT_getNextUnitNumber()
    open(ftn,file=trim(coordfile),status='old')
    write(LVT_logunit,*) '[INFO] Reading USGSGWwell metadata file ',trim(coordfile)
    read(ftn,*)
    read(ftn,*) usgsgwwellobs(i)%nstns
    read(ftn,*) 

    allocate(usgsgwwellobs(i)%stnid(usgsgwwellobs(i)%nstns))
    allocate(usgsgwwellobs(i)%stnlat(usgsgwwellobs(i)%nstns))
    allocate(usgsgwwellobs(i)%stnlon(usgsgwwellobs(i)%nstns))
    allocate(usgsgwwellobs(i)%stnqa(usgsgwwellobs(i)%nstns))
    allocate(usgsgwwellobs(i)%stnsy(usgsgwwellobs(i)%nstns))

    do k=1,usgsgwwellobs(i)%nstns
       read(ftn,*) usgsgwwellobs(i)%stnid(k),usgsgwwellobs(i)%stnlat(k), &
            usgsgwwellobs(i)%stnlon(k),usgsgwwellobs(i)%stnqa(k), &
            usgsgwwellobs(i)%stnsy(k)
       write(LVT_logunit,*) '[INFO] ',&
            usgsgwwellobs(i)%stnid(k),usgsgwwellobs(i)%stnlat(k), &
            usgsgwwellobs(i)%stnlon(k),usgsgwwellobs(i)%stnqa(k), &
            usgsgwwellobs(i)%stnsy(k)
    enddo

    call LVT_releaseUnitNumber(ftn)

    call ESMF_TimeIntervalSet(usgsgwwellobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (gw well obs)')

    call LVT_update_timestep(LVT_rc, 86400)
    usgsgwwellobs(i)%nts = 366
    allocate(usgsgwwellobs(i)%gw(usgsgwwellobs(i)%nstns,usgsgwwellobs(i)%nts))
    allocate(usgsgwwellobs(i)%wtdepth(usgsgwwellobs(i)%nstns,usgsgwwellobs(i)%nts))
    usgsgwwellobs(i)%gw = LVT_rc%udef

  end subroutine USGSGWwell_obsinit


end module USGSGWwell_obsMod
