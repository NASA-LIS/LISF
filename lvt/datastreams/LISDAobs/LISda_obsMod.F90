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
! !MODULE: LISda_obsMod
! \label(LISda_obsMod)
!
! !INTERFACE:
module LISda_obsMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation output as 
!  "observations". 
!  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
!  July 18 2022   Madi Navari Added support for Lambert projection.
! 
!EOP

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISda_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lisdaobs
!
!EOP
  
  type, public :: daobsdec
     character*100               :: odir
     integer                     :: scal
     integer                     :: obstype
     integer                     :: obsinstance
     integer                     :: nc,nr
     real,    allocatable        :: rlat(:)
     real,    allocatable        :: rlon(:)
     integer, allocatable        :: n11(:)

  end type daobsdec

  type(daobsdec), allocatable  :: lisdaobs(:)

contains

!BOP
! 
! !ROUTINE: LISda_obsInit
! \label{LISda_obsInit}
!
! !INTERFACE: 
  subroutine LISda_obsInit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    implicit none
!
! !INPUT PARAMETERS: 
    integer,     intent(IN) :: i   ! index of the observation type
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model output (from a LIS simulation) as observations.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real                    :: run_dd(8)
    integer                 :: t
    integer                 :: ts
    integer                 :: ios
    integer                 :: ftn
    integer                 :: ncId, nrId
    type(ESMF_Config)       :: modelSpecConfig
    character*20            :: domain
    character*100           :: domFile
    character*100           :: map_proj
    character*10            :: time
    character*20            :: obstype
    real                    :: stlat, stlon, dx, dy
    real                    :: gridDesci(50)
    logical                 :: file_exists
    character*20            :: gridtype
    integer                 :: rc
    integer                 :: truelat1, truelat2 
    integer                 :: domain_resolution 
    integer                 :: standard_lon 

    if(.not.allocated(lisdaobs)) then 
       allocate(lisdaobs(LVT_rc%nDataStreams))
    endif

    gridDesci = 0.0

    call ESMF_ConfigGetAttribute(LVT_config,lisdaobs(i)%odir, &
         label="LIS DAOBS output directory:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS output directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,time,&
         label="LIS DAOBS output interval:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS output interval: not defined')

    call LVT_parseTimeString(time,ts)

    call LVT_update_timestep(LVT_rc, ts)

    call ESMF_ConfigGetAttribute(LVT_config,lisdaobs(i)%scal, &
         label="LIS DAOBS use scaled obs:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS use scaled obs: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,obstype,&
         label="LIS DAOBS observation type:",rc=rc)
    if(rc.ne.0) then 
       write(LVT_logunit,*) '[ERR] LIS DAOBS observation type: not defined'
       write(LVT_logunit,*) '[ERR] Supported options are: '
       write(LVT_logunit,*) "[ERR] 'soil moisture' 'snowdepth' 'SWE' 'LAI'"
       call LVT_endrun()
    endif

    call ESMF_ConfigGetAttribute(LVT_config,lisdaobs(i)%obsinstance, &
         label="LIS DAOBS instance index:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS instance index: not defined')

    if(obstype.eq."soil moisture") then 
       lisdaobs(i)%obstype = 1
    elseif(obstype.eq."snowdepth") then 
       lisdaobs(i)%obstype = 2
    elseif(obstype.eq."SWE") then 
       lisdaobs(i)%obstype = 3
    elseif(obstype.eq."LAI") then 
       lisdaobs(i)%obstype = 4
    endif

    call ESMF_ConfigGetAttribute(LVT_config,domFile, &
         label="LIS DAOBS domain file:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS domain file: not defined')

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    inquire(file=trim(domFile), exist=file_exists)
    if(file_exists) then 
       ios = nf90_open(path=domFile,mode=NF90_NOWRITE,ncid=ftn)
       call LVT_verify(ios,'Error in nf90_open in readObsDomainInput')

       ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',map_proj)
       call LVT_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')

       if(map_proj.eq."EQUIDISTANT CYLINDRICAL" .or. map_proj.eq."EASE V2") then
       
          ios = nf90_inq_dimid(ftn,"east_west",ncId)
          call LVT_verify(ios,&
               'Error in nf90_inq_dimid in readObsDomainInput:east_west')
       
          ios = nf90_inq_dimid(ftn,"north_south",nrId)
          call LVT_verify(ios,&
               'Error in nf90_inq_dimid in readObsDomainInput:north_south')
       
          ios = nf90_inquire_dimension(ftn,ncId, len=lisdaobs(i)%nc)
          call LVT_verify(ios,&
               'Error in nf90_inquire_dimension in readObsDomainInput:ncId')
       
          ios = nf90_inquire_dimension(ftn,nrId, len=lisdaobs(i)%nr)
          call LVT_verify(ios,&
               'Error in nf90_inquire_dimension in readObsDomainInput:nrId')
       
          !ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',map_proj)
          !call LVT_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')
       
          ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
               stlat)
          call LVT_verify(ios, &
               'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
       
          ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
               stlon)
          call LVT_verify(ios, &
               'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')
       
          ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',dx)
           call LVT_warning(ios, 'Error in nf90_get_att: DX')
          if(ios.ne.0) then 
             dx = 0 
          endif

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',dy)
          call LVT_warning(ios, 'Error in nf90_get_att: DY')
          if(ios.ne.0) then 
             dy = 0 
          endif

          ios = nf90_get_att(ftn,NF90_GLOBAL,'GRIDTYPE',gridtype)
          call LVT_warning(ios,'Error in nf90_get_att: GRIDTYPE')
       
          if(map_proj.eq."EQUIDISTANT CYLINDRICAL") then              

             gridDesci(1) = 0
             gridDesci(2) = lisdaobs(i)%nc
             gridDesci(3) = lisdaobs(i)%nr                    
             gridDesci(4) = stlat
             gridDesci(5) = stlon
             gridDesci(7) = stlat + (lisdaobs(i)%nr-1)*dy
             gridDesci(8) = stlon + (lisdaobs(i)%nc-1)*dx          
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
          elseif(map_proj.eq."EASE V2") then 

             gridDesci(1) = 9
             gridDesci(2) = lisdaobs(i)%nc
             gridDesci(3) = lisdaobs(i)%nr     
             gridDesci(4) = stlat
             gridDesci(5) = stlon
             gridDesci(6) = 128

             gridDesci(20) = 64

             if(gridtype.eq."M36") then 
                gridDesci(9) = 4
                gridDesci(10) = 0.36
                dx = 0.36
                dy = 0.36
             elseif(gridtype.eq."M09") then 
                gridDesci(9) = 5
                gridDesci(10) = 0.09
                dx = 0.09
                dy = 0.09
             endif
          endif

       elseif(map_proj.eq."LAMBERT CONFORMAL") then   
          ios = nf90_inq_dimid(ftn,"east_west",ncId)
          call LVT_verify(ios,&
               'Error in nf90_inq_dimid in readObsDomainInput:east_west')

          ios = nf90_inq_dimid(ftn,"north_south",nrId)
          call LVT_verify(ios,&
               'Error in nf90_inq_dimid in readObsDomainInput:north_south')

          ios = nf90_inquire_dimension(ftn,ncId, len=lisdaobs(i)%nc)
          call LVT_verify(ios,&
               'Error in nf90_inquire_dimension in readObsDomainInput:ncId')

          ios = nf90_inquire_dimension(ftn,nrId, len=lisdaobs(i)%nr)
          call LVT_verify(ios,&
               'Error in nf90_inquire_dimension in readObsDomainInput:nrId')

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
               stlat)
          call LVT_verify(ios, &
               'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
               stlon)
          call LVT_verify(ios, &
               'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',&
               truelat1)
          call LVT_verify(ios, &
               'Error in nf90_get_att: TRUELAT1')

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT2',&
               truelat2)
          call LVT_verify(ios, &
               'Error in nf90_get_att: TRUELAT1')

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',&
               standard_lon)
          call LVT_verify(ios, &
               'Error in nf90_get_att: STANDARD_LON')

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',&
               domain_resolution)
          call LVT_verify(ios, &
               'Error in nf90_get_att: DX')

          ios = nf90_get_att(ftn,NF90_GLOBAL,'GRIDTYPE',gridtype)
          call LVT_warning(ios,'Error in nf90_get_att: GRIDTYPE')

       gridDesci = 0
       gridDesci(1) = 3                   ! Lambert conic conformal grid
       gridDesci(2) = lisdaobs(i)%nc
       gridDesci(3) = lisdaobs(i)%nr
       gridDesci(4) = stlat               ! latitude of origin -- LL Lat
       gridDesci(5) = stlon               ! longitude of origin -- LL Lon
       gridDesci(6) = 8                   ! Set for Lambert in core/LDT_domainMod.F90 (line 2658)
       gridDesci(7) = truelat2            ! true lat2
       gridDesci(8) = domain_resolution   ! grid spacing in km
       gridDesci(9) = domain_resolution   ! grid spacing in km
       gridDesci(10) = truelat1           ! true lat1
       gridDesci(11) = standard_lon       ! standard long
       gridDesci(20) = 0.0


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

    allocate(lisdaobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(lisdaobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(lisdaobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))

    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         lisdaobs(i)%rlat, lisdaobs(i)%rlon,&
         lisdaobs(i)%n11)

  end subroutine LISda_obsInit
  
end module LISda_obsMod
