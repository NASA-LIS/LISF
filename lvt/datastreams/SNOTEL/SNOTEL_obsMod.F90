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
! !MODULE: SNOTEL_obsMod
! \label(SNOTEL_obsMod)
!
! !INTERFACE:
module SNOTEL_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: SNOTEL_obsinit
  PUBLIC :: SNOTELobs
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the SNOTEL SWE 
!  data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Aug 2009   Sujay Kumar  Initial Specification
! 
!EOP

  type, public :: snotelobsdec
     character*100        :: odir
     integer              :: nstns
     integer              :: nstates
     real                 :: udef
     integer              :: nts
     character*2, allocatable :: statename(:)
     character*37,allocatable :: stnname(:)
     character*6, allocatable :: stnid(:)
     real,        allocatable :: stnlat(:)
     real,        allocatable :: stnlon(:)
     type(ESMF_Clock)     :: clock
     type(ESMF_Time)      :: startTime, stopTime
     type(ESMF_TimeInterval) :: timestep
     real,  allocatable          :: swe(:,:)
     real,  allocatable          :: prcp(:,:)
     integer              :: yr
  end type snotelobsdec

  type(snotelobsdec), allocatable :: snotelobs(:)

contains
  
!BOP
! 
! !ROUTINE: SNOTEL_obsInit
! \label{SNOTEL_obsInit}
!
! !INTERFACE: 
  subroutine SNOTEL_obsinit(i)
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
!  for reading SNOTEL data. 
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
    integer            :: siteid
    character*100      :: coordfile
    character*100      :: mdata
    
    if(.not.allocated(snotelobs)) then 
       allocate(snotelobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, snotelobs(i)%odir, &
         label='SNOTEL observation directory:',rc=status)
    call LVT_verify(status, 'SNOTEL observation directory: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, coordfile, &
         label='SNOTEL coord file:',rc=status)
    call LVT_verify(status, 'SNOTEL coord file: not defined')

!    call ESMF_ConfigGetAttribute(LVT_config, mdata, &
!         label='SNOTEL metadata file:',rc=status)
!    call LVT_verify(status, 'SNOTEL metadata file: not defined')

    snotelobs(i)%nstns = 813

!    ftn=LVT_getNextUnitNumber()
!    open(ftn,file=trim(mdata),status='old')
!    write(LVT_logunit,*) '[INFO] Reading SNOTEL metadata file ',trim(mdata)
!    read(ftn,*)
!    read(ftn,*) snotelobs(i)%nstns, snotelobs(i)%udef, syr, smo, &
!         sda, shr, smn, eyr, emo, &
!         eda, ehr, emn, ts   
!
!    write(LVT_logunit,*) '[INFO] ',snotelobs(i)%nstns, snotelobs(i)%udef, syr, smo, &
!         sda, shr, smn, eyr, emo, &
!         eda, ehr, emn, ts    

!    call LVT_releaseUnitNumber(ftn)

    ts = 86400
    call ESMF_TimeIntervalSet(snotelobs(i)%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (snotelobs(i))')
    call LVT_update_timestep(LVT_rc, 86400)

    snotelobs(i)%nts = 366 !yearly data

    allocate(snotelobs(i)%statename(snotelobs(i)%nstns))
    allocate(snotelobs(i)%stnid(snotelobs(i)%nstns))
    allocate(snotelobs(i)%stnname(snotelobs(i)%nstns))
    allocate(snotelobs(i)%stnlat(snotelobs(i)%nstns))
    allocate(snotelobs(i)%stnlon(snotelobs(i)%nstns))

    ftn = LVT_getNextUnitNumber()
    open(ftn,file=trim(coordfile),form='formatted')
    
    write(LVT_logunit,*) '[INFO] opening coord file ',trim(coordfile)
    do k=1,snotelobs(i)%nstns
       read(ftn,222) snotelobs(i)%statename(k), snotelobs(i)%stnname(k), &
            snotelobs(i)%stnid(k),siteid, &
            snotelobs(i)%stnlat(k), snotelobs(i)%stnlon(k)
       write (LVT_logunit,*) '[INFO] ',&
            snotelobs(i)%stnid(k), snotelobs(i)%stnname(k), &
            siteid,snotelobs(i)%stnlat(k), snotelobs(i)%stnlon(k)
222 format(A2,A30,A6,I14,2F10.3)
    enddo

    call LVT_releaseUnitNumber(ftn)

    allocate(snotelobs(i)%swe(snotelobs(i)%nstns,snotelobs(i)%nts))
    allocate(snotelobs(i)%prcp(snotelobs(i)%nstns,snotelobs(i)%nts))
    snotelobs(i)%swe = LVT_rc%udef
    snotelobs(i)%prcp = LVT_rc%udef
    snotelobs(i)%yr = -1
  end subroutine SNOTEL_obsinit


end module SNOTEL_obsMod
