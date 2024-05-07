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
! !MODULE: ERA5obsMod
! \label(ERA5obsMod)
!
! !INTERFACE:
module ERA5obsMod
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
!  5 Dec 2020   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ERA5obsinit !Initializes structures for reading ERA5 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ERA5obs !Object to hold ERA5 observation attributes
!EOP

  type, public :: era5dec
     character*100           :: odir
     character*100           :: mapfile
     integer                 :: nc, nr,npts,ntimes
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: mo
     real                    :: datares
     logical                 :: startFlag
     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: ts
     real, allocatable       :: prcp(:,:,:)
     real, allocatable       :: tair(:,:,:)   
     integer, allocatable    :: G2P(:,:)
  end type era5dec
     
  type(era5dec), allocatable :: ERA5Obs(:)

contains
  
!BOP
! 
! !ROUTINE: ERA5obsInit
! \label{ERA5obsInit}
!
! !INTERFACE: 
  subroutine ERA5obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod    
    use LVT_logMod
    use LVT_timeMgrMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the ERA5 data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer               :: c,r,k,status
    integer               :: ftn
    integer               :: G2Pid

    integer                 :: iret     
    real, allocatable       :: var_inp_1d(:)
    logical*1, allocatable  :: input_bitmap(:)
    logical*1               :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)

    if(.not.allocated(ERA5obs)) then 
       allocate(ERA5obs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigFindLabel(LVT_Config, &
         label='ERA5 data directory:',rc=status)
    call ESMF_ConfigGetAttribute(LVT_Config, ERA5obs(i)%odir, &
         rc=status)
    call LVT_verify(status, 'ERA5 data directory: not defined')

    call ESMF_ConfigFindLabel(LVT_Config, &
         label='ERA5 forcing tile to grid mapping file:',rc=status)
    call ESMF_ConfigGetAttribute(LVT_Config, ERA5obs(i)%mapfile, &
         rc=status)
    call LVT_verify(status, 'ERA5 forcing tile to grid mapping file: not defined')

    era5obs(i)%gridDesc = 0
        
    era5obs(i)%nc = 1440
    era5obs(i)%nr =  720
    era5obs(i)%npts = 340819
    !filling the items needed by the interpolation library
    era5obs(i)%gridDesc(1) = 0  
    era5obs(i)%gridDesc(2) = era5obs(i)%nc
    era5obs(i)%gridDesc(3) = era5obs(i)%nr
    era5obs(i)%gridDesc(4) = -89.875
    era5obs(i)%gridDesc(5) = -179.875
    era5obs(i)%gridDesc(7) = 89.875
    era5obs(i)%gridDesc(8) = 179.875
    era5obs(i)%gridDesc(6) = 128
    era5obs(i)%gridDesc(9) = 0.25
    era5obs(i)%gridDesc(10) = 0.25
    era5obs(i)%gridDesc(20) = 0

    era5obs(i)%datares  = 0.25

    if(LVT_isAtAfinerResolution(era5obs(i)%datares)) then
       
       allocate(era5obs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(era5obs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(era5obs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(era5obs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            era5obs(i)%rlat, &
            era5obs(i)%rlon, &
            era5obs(i)%n11)
    else
       allocate(era5obs(i)%n11(era5obs(i)%nc*era5obs(i)%nr))
       call upscaleByAveraging_input(era5obs(i)%gridDesc,&
            LVT_rc%gridDesc,era5obs(i)%nc*era5obs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,era5obs(i)%n11)
    endif


    call ESMF_TimeIntervalSet(era5obs(i)%ts, s=3600,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')

    call LVT_update_timestep(LVT_rc, 3600)

    era5obs(i)%mo = -1
    era5obs(i)%startFlag = .true.

    era5obs(i)%ntimes = 745

    allocate(era5obs(i)%prcp(LVT_rc%lnc,LVT_rc%lnr,era5obs(i)%ntimes))
    allocate(era5obs(i)%tair(LVT_rc%lnc,LVT_rc%lnr,era5obs(i)%ntimes))

    era5obs(i)%prcp = LVT_rc%udef
    era5obs(i)%tair = LVT_rc%udef
    
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 

    allocate(era5obs(i)%G2P(era5obs(i)%nc,&
         era5obs(i)%nr))
    
    call LVT_verify(nf90_open(path=trim(era5obs(i)%mapfile), &
         mode=NF90_NOWRITE, &
         ncid=ftn), 'nf90_open failed for '//trim(era5obs(i)%mapfile))
    
    call LVT_verify(nf90_inq_varid(ftn,'G2P',G2PId), &
         'nf90_inq_varid failed for G2P in read_era5')
    
    call LVT_verify(nf90_get_var(ftn,G2PId, era5obs(i)%G2P),&
         'nf90_get_var failed for G2P in read_era5') 
    call LVT_verify(nf90_close(ftn))
#endif

  end subroutine ERA5obsinit


end module ERA5obsMod
