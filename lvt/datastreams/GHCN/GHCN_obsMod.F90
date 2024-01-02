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
! !MODULE: GHCN_obsMod
!  \label(GHCN_obsMod)
!
! !INTERFACE:
module GHCN_obsMod
! 
! !USES: 
  use ESMF

  implicit none
  PRIVATE 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the GHCN 
!  data (snow depth and precip). The data download link is: 
!  
!  ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Apr 2009   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GHCN_obsinit !Initializes structures for reading GHCN data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GHCNobs !Object to hold GHCN observation attributes
!EOP
  type, public :: ghcnobsdec
     character*100        :: odir
     integer              :: nstns
     real                 :: udef
     integer              :: nts
     character*100, allocatable :: stnid(:)
     real,        allocatable :: stnlat(:)
     real,        allocatable :: stnlon(:)
     real,        allocatable :: stnelev(:)
     type(ESMF_Clock)     :: clock
     type(ESMF_Time)      :: startTime
     type(ESMF_TimeInterval) :: timestep
     real,  allocatable          :: snod(:,:)
     real,  allocatable          :: prcp(:,:)
     logical              :: startflag
     integer              :: yr,syr
     integer              :: mo,smo
  end type ghcnobsdec

  type(ghcnobsdec), save :: ghcnobs(2)

contains
  
!BOP
! 
! !ROUTINE: GHCN_obsInit
! \label{GHCN_obsInit}
!
! !INTERFACE: 
  subroutine GHCN_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod
    use map_utils

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading GHCN data. 
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
    character*50       :: stnid
    real               :: stnlat, stnlon,stnelev
    character*100      :: stationfile
    integer            :: c,r
    real               :: col,row
    
    call ESMF_ConfigGetAttribute(LVT_config, ghcnobs(i)%odir, &
         label='GHCN observation directory:',rc=status)
    call LVT_verify(status, 'GHCN observation directory: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, stationfile, &
         label='GHCN station file:',rc=status)
    call LVT_verify(status, 'GHCN station file: not defined')

    ftn=LVT_getNextUnitNumber()
    open(ftn,file=trim(stationfile),status='old')
    write(LVT_logunit,*) '[INFO] Reading GHCN station file ',trim(stationfile)
! first read to figure out the number of stations included in the 
! LVT/LIS domain
    status = 0 
    do while(status.eq.0) 
       read(ftn,*,iostat=status) stnid, stnlat, stnlon,stnelev
       if(status.eq.0) then 
          call latlon_to_ij(LVT_domain%lvtproj, stnlat, stnlon, &
               col, row)
          c = nint(col)
          r = nint(row)
          if(c.ge.1.and.c.le.LVT_rc%lnc.and.&
               r.ge.1.and.r.le.LVT_rc%lnr) then 
             ghcnobs(i)%nstns =  ghcnobs(i)%nstns +1
          endif
       else
          exit
       endif
    enddo
    call LVT_releaseUnitNumber(ftn)

    allocate(ghcnobs(i)%stnid(ghcnobs(i)%nstns))
    allocate(ghcnobs(i)%stnlat(ghcnobs(i)%nstns))
    allocate(ghcnobs(i)%stnlon(ghcnobs(i)%nstns))
    allocate(ghcnobs(i)%stnelev(ghcnobs(i)%nstns))

    ftn=LVT_getNextUnitNumber()
    open(ftn,file=trim(stationfile),status='old')

    status = 0 
    k = 0 
    do while(status.eq.0) 
       read(ftn,*,iostat=status) stnid, stnlat, stnlon,stnelev
       if(status.eq.0) then 
          call latlon_to_ij(LVT_domain%lvtproj, stnlat, stnlon, &
               col, row)
          c = nint(col)
          r = nint(row)
          if(c.ge.1.and.c.le.LVT_rc%lnc.and.&
               r.ge.1.and.r.le.LVT_rc%lnr) then 
             k = k + 1
             ghcnobs(i)%stnid(k) = stnid
             ghcnobs(i)%stnlat(k) = stnlat
             ghcnobs(i)%stnlon(k) = stnlon
             ghcnobs(i)%stnelev(k) = stnelev
          endif
       else
          exit
       endif
    enddo
    call LVT_releaseUnitNumber(ftn)

222 format(A11,2F10.4,F7.1)
       
    ghcnobs(i)%nts = 366 

    call ESMF_TimeIntervalSet(ghcnobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (ghcnobs)')

    allocate(ghcnobs(i)%snod(ghcnobs(i)%nstns, ghcnobs(i)%nts))
    allocate(ghcnobs(i)%prcp(ghcnobs(i)%nstns, ghcnobs(i)%nts))
    ghcnobs(i)%snod = LVT_rc%udef
    ghcnobs(i)%prcp = LVT_rc%udef

    ghcnobs(i)%yr = -1
    ghcnobs(i)%mo = -1

    ghcnobs(i)%startflag = .true.
    call LVT_update_timestep(LVT_rc, 86400)
 
  end subroutine GHCN_obsinit

end module GHCN_obsMod
