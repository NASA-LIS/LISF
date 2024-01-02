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
! !MODULE: GRACE_obsMod
! \label(GRACE_obsMod)
!
! !INTERFACE:
module GRACE_obsMod
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
!  18 May 2011: Sujay Kumar,  Initial Specification
!   8 Jun 2018: Kristi Arsenault, Updated to support newer GRACE datasets
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACE_obsinit !Initializes structures for reading GRACE data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACEobs !Object to hold GRACE observation attributes
!EOP

  type, public :: gracedec
     logical                 :: startFlag
     integer                 :: yr
     integer                 :: mo
     integer                 :: da
     integer                 :: nc
     integer                 :: nr
     real, allocatable       :: tvals(:)
     real, allocatable       :: time_bounds(:,:)
     real, allocatable       :: lwe_thickness(:,:,:)
     real, allocatable       :: scalefactor(:,:)  !BZ
     integer                 :: reftime
     integer                 :: tdims
     character*120           :: filename
     character*150           :: gracescalefile
     character*100           :: datasource
     character*20            :: gridtransformopt

     ! Bilinear interpolation
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
     ! Nearest-neighbor
     integer, allocatable    :: n111(:)

  end type gracedec
     
  type(gracedec), allocatable :: GRACEObs(:)

contains
  
!BOP
! 
! !ROUTINE: GRACE_obsInit
! \label{GRACE_obsInit}
!
! !INTERFACE: 
  subroutine GRACE_obsinit(i)
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
!   for reading the GRACE data, including the computation of spatial 
!   interpolation weights. The GRACE data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification 
!   8 Jun 2018: Kristi Arsenault, Updated to support newer GRACE datasets
! 
!EOP

    integer              :: i
    integer              :: status
    real                 :: gridDesci(50)
    character*100        :: domFile
    character*100        :: map_proj
    logical              :: file_exists

    if(.not.allocated(GRACEobs)) then 
       allocate(GRACEobs(LVT_rc%nDataStreams))
    endif

    ! LVT config options for GRACE Obs:
    call ESMF_ConfigGetAttribute(LVT_Config, GRACEObs(i)%filename, &
         label='GRACE data filename: ',rc=status)
    call LVT_verify(status, 'GRACE data filename: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, GRACEObs(i)%gracescalefile, &
         label='GRACE scale factor filename: ',rc=status)
    call LVT_verify(status, 'GRACE scale factor filename: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, GRACEObs(i)%datasource, &
         label="GRACE data source:",rc=status)
    if(status.ne.0) then
       write(LVT_logunit,*) '[ERR] GRACE data source: not defined'
       write(LVT_logunit,*) '[ERR] The options are ..'
       write(LVT_logunit,*) "[ERR] 'GRACE TWS Mascon 0.5 deg' or "
       write(LVT_logunit,*) "[ERR] 'GRACE TWS Original 1 deg'"
       call LVT_endrun()
    endif
    ! Option for modifying input data to final output grid:
    call ESMF_ConfigGetAttribute(LVT_config, GRACEObs(i)%gridtransformopt, &
         label="GRACE data grid transform option:",rc=status)
    call LVT_verify(status,'GRACE data grid transform option: not defined')
    if(status.ne.0) then
       write(LVT_logunit,*) '[ERR] GRACE data grid transform: not defined'
       write(LVT_logunit,*) '[ERR] The current options are ..'
       write(LVT_logunit,*) "[ERR] 'neighbor' or "
       write(LVT_logunit,*) "[ERR] 'bilinear'"
       call LVT_endrun()
    endif

    ! Set the GRACE data file domain: 

       ! Original 1.0 deg SH version:
    if( GRACEObs(i)%datasource.eq."GRACE TWS Original 1 deg" ) then
       GRACEObs(i)%nc = 360
       GRACEObs(i)%nr = 180
       gridDesci    = 0
       gridDesci(1) = 0
       gridDesci(2) = GRACEObs(i)%nc
       gridDesci(3) = GRACEObs(i)%nr
       gridDesci(4) = -89.5
       gridDesci(5) = -179.5
       gridDesci(7) = 89.5
       gridDesci(8) = 179.5
       gridDesci(6) = 128
       gridDesci(9) = 1.0
       gridDesci(10) = 1.0
       gridDesci(20) = 64

    elseif( GRACEObs(i)%datasource.eq."GRACE TWS Mascon 0.5 deg" ) then
       GRACEObs(i)%nc = 720
       GRACEObs(i)%nr = 360
       gridDesci    = 0
       gridDesci(1) = 0
       gridDesci(2) = GRACEObs(i)%nc
       gridDesci(3) = GRACEObs(i)%nr
       gridDesci(4) = -89.75
       gridDesci(5) = -179.75
       gridDesci(6) = 128
       gridDesci(7) = 89.75
       gridDesci(8) = 179.75
       gridDesci(9) = 0.5
       gridDesci(10) = 0.5
       gridDesci(20) = 64
    endif
    
! ---
    
    allocate(GRACEObs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(GRACEObs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    ! Bilinear interpolation option:
    if( trim(GRACEObs(i)%gridtransformopt) == "bilinear" ) then
      allocate(GRACEObs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
      allocate(GRACEObs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
      allocate(GRACEObs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
      allocate(GRACEObs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
      allocate(GRACEObs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
      allocate(GRACEObs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
      allocate(GRACEObs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
      allocate(GRACEObs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
      call bilinear_interp_input(gridDesci,LVT_rc%gridDesc, &
           LVT_rc%lnc*LVT_rc%lnr,                  &
           GRACEObs(i)%rlat, GRACEObs(i)%rlon,     &
           GRACEObs(i)%n11, GRACEObs(i)%n12,       &
           GRACEObs(i)%n21, GRACEObs(i)%n22,       &
           GRACEObs(i)%w11, GRACEObs(i)%w12,       &
           GRACEObs(i)%w21, GRACEObs(i)%w22)

    elseif( trim(GRACEObs(i)%gridtransformopt) == "neighbor" ) then
       ! Use nearest neighbor to have product on GRACE product grid:
       allocate( GRACEObs(i)%n111(LVT_rc%lnc*LVT_rc%lnr) )
       call neighbor_interp_input( gridDesci,&
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            GRACEObs(i)%rlat, &
            GRACEObs(i)%rlon, &
            GRACEObs(i)%n111 )
    endif

    GRACEobs(i)%yr = -1
!    GRACEobs(i)%mo = LVT_rc%mo
    GRACEobs(i)%mo = -1
    GRACEobs(i)%da = LVT_rc%da
    GRACEobs(i)%startFlag = .true. 

!    if(LVT_rc%tavgInterval.lt.2592000) then 
!       write(LVT_logunit,*) 'The time averaging interval must be greater than'
!       write(LVT_logunit,*) 'equal to a month since the GRACE data is monthly'
!       call LVT_endrun()
!    endif

    call LVT_get_julhr(2002,1,1,0,0,0,GRACEobs(i)%reftime)

    call LVT_update_timestep(LVT_rc, 86400)

  end subroutine GRACE_obsinit

end module GRACE_obsMod
