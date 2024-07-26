!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: simGRACE_obsMod
! \label(simGRACE_obsMod)
!
! !INTERFACE:
module simGRACE_obsMod
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
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: simGRACE_obsinit !Initializes structures for reading simGRACE data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: simGRACEobs !Object to hold simGRACE observation attributes
!EOP

  type, public :: simgracedec
     character*100           :: odir
     character*50            :: config
     logical                 :: startFlag
     integer                 :: useRawData
     integer                 :: yr
     integer                 :: mo
     integer                 :: da

     integer                 :: nc
     integer                 :: nr
     real, allocatable       :: rlat(:)
     real, allocatable       :: rlon(:)
     integer, allocatable    :: n11(:)
     integer, allocatable    :: n12(:)
     integer, allocatable    :: n21(:)
     integer, allocatable    :: n22(:)     
     real,    allocatable    :: w11(:)
     real,    allocatable    :: w12(:)
     real,    allocatable    :: w21(:)
     real,    allocatable    :: w22(:)

     real, allocatable       :: rlat1(:)
     real, allocatable       :: rlon1(:)
     integer, allocatable    :: n111(:)
  end type simgracedec
     
  type(simgracedec), allocatable :: simGRACEObs(:)

contains
  
!BOP
! 
! !ROUTINE: simGRACE_obsInit
! \label{simGRACE_obsInit}
!
! !INTERFACE: 
  subroutine simGRACE_obsinit(i)
! 
! !USES: 
    use LVT_coreMod,   only : LVT_rc, LVT_Config
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the simGRACE data, including the computation of spatial 
!   interpolation weights. The simGRACE data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: i
    integer               :: status
    real                  :: gridDesci(50)
    character*100           :: domFile
    character*100           :: map_proj
    logical                 :: file_exists
    real                    :: stlat, stlon, dx, dy
    integer                 :: ftn,ncId, nrId,ios

    if(.not.allocated(simGRACEobs)) then 
       allocate(simGRACEobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, simGRACEObs(i)%odir, &
         label='simulated GRACE data directory: ',rc=status)
    call LVT_verify(status, 'simulated GRACE data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, simGRACEObs(i)%config, &
         label='simulated GRACE configuration: ',rc=status)
    call LVT_verify(status, 'simulated GRACE configuration: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,domFile, &
         label="simulated GRACE observation domain file:",rc=status)
    call LVT_verify(status,'simulated GRACE observation domain file: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, simGRACEObs(i)%useRawData, &
         label='simulated GRACE process raw (anomaly) data: ',rc=status)
    call LVT_verify(status, 'simulated GRACE process raw (anomaly) data: not defined')

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    inquire(file=trim(domFile), exist=file_exists)
    if(file_exists) then 
       ios = nf90_open(path=domFile,mode=NF90_NOWRITE,ncid=ftn)
       call LVT_verify(ios,'Error in nf90_open in readObsDomainInput')
       
       ios = nf90_inq_dimid(ftn,"east_west",ncId)
       call LVT_verify(ios,&
            'Error in nf90_inq_dimid in readObsDomainInput:east_west')
       
       ios = nf90_inq_dimid(ftn,"north_south",nrId)
       call LVT_verify(ios,&
            'Error in nf90_inq_dimid in readObsDomainInput:north_south')
       
       ios = nf90_inquire_dimension(ftn,ncId, len=simGRACEobs(i)%nc)
       call LVT_verify(ios,&
            'Error in nf90_inquire_dimension in readObsDomainInput:ncId')
       
       ios = nf90_inquire_dimension(ftn,nrId, len=simGRACEobs(i)%nr)
       call LVT_verify(ios,&
            'Error in nf90_inquire_dimension in readObsDomainInput:nrId')
       
       ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',map_proj)
       call LVT_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')
       
       ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
            stlat)
       call LVT_verify(ios, &
            'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
       
       ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
            stlon)
       call LVT_verify(ios, &
            'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')
       
       ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',dx)
       call LVT_verify(ios, 'Error in nf90_get_att: DX')
       
       ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',dy)
       call LVT_warning(ios, 'Error in nf90_get_att: DY')
       
       if(map_proj.eq."EQUIDISTANT CYLINDRICAL") then              

          gridDesci(1) = 0
          gridDesci(2) = simGRACEobs(i)%nc
          gridDesci(3) = simGRACEobs(i)%nr                    
          gridDesci(4) = stlat
          gridDesci(5) = stlon
          gridDesci(7) = stlat + (simGRACEobs(i)%nr-1)*dy
          gridDesci(8) = stlon + (simGRACEobs(i)%nc-1)*dx          
          gridDesci(9) = dx
          
          if(gridDesci(1).eq.0) then 
             gridDesci(10) = dy
             gridDesci(6) = 128
             gridDesci(11) = 64
             gridDesci(20) = 64
          endif
          if(gridDesci(7).lt.gridDesci(4)) then
             write(LVT_logunit,*) '[ERR] lat2 must be greater than lat1'
             write(LVT_logunit,*) '[ERR] ',gridDesci(7),&
                  gridDesci(4)
             write(LVT_logunit,*) '[ERR] Stopping run...'
             call LVT_endrun
          endif
          if(gridDesci(8).lt.gridDesci(5)) then
             write(LVT_logunit,*) '[ERR] lon2 must be greater than lon1'
             write(LVT_logunit,*) '[ERR] ',gridDesci(8),&
                  gridDesci(5)
             write(LVT_logunit,*) '[ERR] Stopping run...'
             call LVT_endrun
          endif
          
       else
          write(LVT_logunit,*) '[ERR] Map projection ',trim(map_proj) 
          write(LVT_logunit,*) '[ERR] not currently supported for LIS DAOBS plugin'
          call LVT_endrun()
       endif
                
    else
       write(LVT_logunit,*) '[ERR] LIS DAOBS domain file ',trim(domFile)
       write(LVT_logunit,*) '[ERR] does not exist..'
       call LVT_endrun()
    endif
#endif

    allocate(simGRACEobs(i)%rlat1(LVT_rc%lnc*LVT_rc%lnr))
    allocate(simGRACEobs(i)%rlon1(LVT_rc%lnc*LVT_rc%lnr))
    allocate(simGRACEobs(i)%n111(LVT_rc%lnc*LVT_rc%lnr))

    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         simGRACEobs(i)%rlat1, simGRACEobs(i)%rlon1,&
         simGRACEobs(i)%n111)

    if(simGRACEobs(i)%userawdata.gt.0) then 

       allocate(simGRACEObs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(simGRACEObs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       
       gridDesci    = 0
       gridDesci(1) = 0
       gridDesci(2) = 360
       gridDesci(3) = 180
       gridDesci(4) = -89.5
       gridDesci(5) = -179.5
       gridDesci(7) = 89.5
       gridDesci(8) = 179.5
       gridDesci(6) = 128
       gridDesci(9) = 1.0
       gridDesci(10) = 1.0
       gridDesci(20) = 64
       
       simGRACEObs(i)%nc = 360
       simGRACEObs(i)%nr = 180
       
       call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,            &
            LVT_rc%lnc*LVT_rc%lnr,                &
            simGRACEObs(i)%rlat, simGRACEObs(i)%rlon,     &
            simGRACEObs(i)%n11, simGRACEObs(i)%n12,       &
            simGRACEObs(i)%n21, simGRACEObs(i)%n22,       &
            simGRACEObs(i)%w11, simGRACEObs(i)%w12,       &
            simGRACEObs(i)%w21, simGRACEObs(i)%w22)
       
    endif

    simGRACEobs(i)%yr = -1
!    simGRACEobs(i)%mo = LVT_rc%mo
    simGRACEobs(i)%mo = -1
    simGRACEobs(i)%da = -1
    simGRACEobs(i)%startFlag = .true. 

!    if(LVT_rc%tavgInterval.lt.2592000) then 
!       write(LVT_logunit,*) 'The time averaging interval must be greater than'
!       write(LVT_logunit,*) 'equal to a month since the simGRACE data is monthly'
!       call LVT_endrun()
!    endif

    call LVT_update_timestep(LVT_rc, 86400)

  end subroutine simGRACE_obsinit


end module simGRACE_obsMod
