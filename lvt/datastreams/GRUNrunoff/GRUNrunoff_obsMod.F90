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
! !MODULE: GRUNrunoff_obsMod
! \label(GRUNrunoff_obsMod)
!
! !INTERFACE:
module GRUNrunoff_obsMod
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
!  This module adds the support for GRUN runoff data, covering a time
!  period of 1902 - 2004.
!  https://essd.copernicus.org/articles/11/1655/2019/
!  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  30 Jul 2022   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GRUNrunoffinit !Initializes structures for reading GRUN data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GRUNobs !Object to hold GRUN observation attributes
!EOP

  type, public :: grundec
     character*100           :: odir
     character*100           :: mapfile
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: mo
     real                    :: datares
     logical                 :: startFlag
     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: ts

  end type grundec
     
  type(grundec), allocatable :: GRUNObs(:)

contains
  
!BOP
! 
! !ROUTINE: GRUNrunoffInit
! \label{GRUNrunoffInit}
!
! !INTERFACE: 
  subroutine GRUNrunoffinit(i)
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
!   for reading the GRUN data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer               :: c,r,k,status
    integer               :: ftn

    integer                 :: iret     
    real, allocatable       :: var_inp_1d(:)
    logical*1, allocatable  :: input_bitmap(:)
    logical*1               :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)

    if(.not.allocated(GRUNobs)) then 
       allocate(GRUNobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigFindLabel(LVT_Config, &
         label='GRUN data directory:',rc=status)
    call ESMF_ConfigGetAttribute(LVT_Config, GRUNobs(i)%odir, &
         rc=status)
    call LVT_verify(status, 'GRUN data directory: not defined')

    grunobs(i)%gridDesc = 0
        
    grunobs(i)%nc = 720
    grunobs(i)%nr = 360

    !filling the items needed by the interpolation library
    grunobs(i)%gridDesc(1) = 0  
    grunobs(i)%gridDesc(2) = grunobs(i)%nc
    grunobs(i)%gridDesc(3) = grunobs(i)%nr
    grunobs(i)%gridDesc(4) = -89.75
    grunobs(i)%gridDesc(5) = -179.75
    grunobs(i)%gridDesc(7) = 89.75
    grunobs(i)%gridDesc(8) = 179.75
    grunobs(i)%gridDesc(6) = 128
    grunobs(i)%gridDesc(9) = 0.5
    grunobs(i)%gridDesc(10) = 0.5
    grunobs(i)%gridDesc(20) = 0

    grunobs(i)%datares  = 0.5

    if(LVT_isAtAfinerResolution(grunobs(i)%datares)) then
       
       allocate(grunobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(grunobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(grunobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(grunobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            grunobs(i)%rlat, &
            grunobs(i)%rlon, &
            grunobs(i)%n11)
    else
       allocate(grunobs(i)%n11(grunobs(i)%nc*grunobs(i)%nr))
       call upscaleByAveraging_input(grunobs(i)%gridDesc,&
            LVT_rc%gridDesc,grunobs(i)%nc*grunobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,grunobs(i)%n11)
    endif

    call ESMF_TimeIntervalSet(grunobs(i)%ts, s=86400,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: GRUNrunoff_obsMod ')

    call LVT_update_timestep(LVT_rc, 86400)

    grunobs(i)%mo = -1
    grunobs(i)%startFlag = .true.

    call ESMF_TimeSet(GRUNobs(i)%starttime, yy=1902, &
         mm = 1,&
         dd = 1,&
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status,'error in timeset: GRUNrunoff_obsMod')
 
  end subroutine GRUNrunoffinit
  
end module GRUNrunoff_obsMod
