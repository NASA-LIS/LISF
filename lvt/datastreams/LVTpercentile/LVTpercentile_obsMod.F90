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
! !MODULE: LVTpercentile_obsMod
! \label(LVTpercentile_obsMod)
!
! !INTERFACE:
module LVTpercentile_obsMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the use of a LVT percentile output as 
!  "observations". 
!  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  08 Mar 2017   Sujay Kumar  Initial Specification
! 
!EOP

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LVTpercentile_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lvtpercobs
!
!EOP
  
  type, public :: lvtpercobsdec
     character*100               :: odir
     integer                     :: nc,nr
     character*100               :: var_name
     real,    allocatable        :: rlat(:)
     real,    allocatable        :: rlon(:)
     integer, allocatable        :: n11(:)

  end type lvtpercobsdec

  type(lvtpercobsdec), allocatable :: lvtpercobs(:)

contains

!BOP
! 
! !ROUTINE: LVTpercentile_obsInit
! \label{LVTpercentile_obsInit}
!
! !INTERFACE: 
  subroutine LVTpercentile_obsInit(i)
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
    integer                 :: k,rc

    if(.not.allocated(lvtpercobs)) then 
       allocate(lvtpercobs(LVT_rc%nDataStreams))
    endif

    gridDesci = 0.0

    call ESMF_ConfigFindLabel(LVT_config,  & 
         label='LVT percentile output directory:', rc=rc)
    do k=1,i
       call ESMF_ConfigGetAttribute(LVT_Config, lvtpercobs(i)%odir, &
            rc=rc)
       call LVT_verify(rc,'LVT percentile output directory: not defined')
    enddo

    ts = 86400
    call LVT_update_timestep(LVT_rc, ts)

    call ESMF_ConfigFindLabel(LVT_config,  & 
         label='LVT percentile variable name:', rc=rc)
    do k=1,i
       call ESMF_ConfigGetAttribute(LVT_Config, lvtpercobs(i)%var_name, &
            rc=rc)
       call LVT_verify(rc,'LVT percentile variable name: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config,  & 
         label='LVT percentile domain file:', rc=rc)
    do k=1,i
       call ESMF_ConfigGetAttribute(LVT_Config, domFile,&
            rc=rc)
       call LVT_verify(rc,'LVT percentile domain file: not defined')
    enddo

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
       
       ios = nf90_inquire_dimension(ftn,ncId, len=lvtpercobs(i)%nc)
       call LVT_verify(ios,&
            'Error in nf90_inquire_dimension in readObsDomainInput:ncId')
       
       ios = nf90_inquire_dimension(ftn,nrId, len=lvtpercobs(i)%nr)
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
          gridDesci(2) = lvtpercobs(i)%nc
          gridDesci(3) = lvtpercobs(i)%nr                    
          gridDesci(4) = stlat
          gridDesci(5) = stlon
          gridDesci(7) = stlat + (lvtpercobs(i)%nr-1)*dy
          gridDesci(8) = stlon + (lvtpercobs(i)%nc-1)*dx          
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

    allocate(lvtpercobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(lvtpercobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(lvtpercobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))

    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         lvtpercobs(i)%rlat, lvtpercobs(i)%rlon,&
         lvtpercobs(i)%n11)

  end subroutine LVTpercentile_obsInit
  
end module LVTpercentile_obsMod
