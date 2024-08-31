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
! !MODULE: MERRA2obsMod
! \label(MERRA2obsMod)
!
! !INTERFACE:
module MERRA2obsMod
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
!  7 Mar 2015   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MERRA2obsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MERRA2obs !Object to hold MERRA2 observation attributes
!EOP

  type, public :: merra2dec
     character*100           :: odir
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: da
     real                    :: datares
     logical                 :: startFlag
     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: ts
     real, allocatable       :: qs(:,:,:)
     real, allocatable       :: qsb(:,:,:)
     real, allocatable       :: swnet(:,:,:)
     real, allocatable       :: qle(:,:,:)
     real, allocatable       :: qh(:,:,:)
     real, allocatable       :: frsno(:,:,:)
     real, allocatable       :: snod(:,:,:)
     real, allocatable       :: swe(:,:,:)
     real, allocatable       :: qg(:,:,:)
     real, allocatable       :: sfsm(:,:,:)
     real, allocatable       :: rzsm(:,:,:)
     real, allocatable       :: prcp(:,:,:)
     real, allocatable       :: tair(:,:,:)
     real, allocatable       :: tskin(:,:,:)
     
     integer                 :: usecorr   !Yeosang Yoon
     integer                 :: usescalef
     character*100           :: scaleffile
     real, allocatable       :: refxrange(:,:,:,:)
     real, allocatable       :: refcdf(:,:,:,:)
     real, allocatable       :: merraxrange(:,:,:,:)
     real, allocatable       :: merracdf(:,:,:,:)

  end type merra2dec
     
  type(merra2dec), allocatable :: MERRA2Obs(:)

contains
  
!BOP
! 
! !ROUTINE: MERRA2obsInit
! \label{MERRA2obsInit}
!
! !INTERFACE: 
  subroutine MERRA2obsinit(i)
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
!   for reading the GIMMS NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer,  parameter   :: nt = 12, nbins = 50
    integer               :: c,r,k,status
    integer               :: ftn
    integer               :: pcp1id,pcp2id,pcp3id,pcp4id

    integer                 :: iret     
    real, allocatable       :: var_inp_1d(:)
    logical*1, allocatable  :: input_bitmap(:)
    logical*1               :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)

    if(.not.allocated(MERRA2obs)) then 
       allocate(MERRA2obs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigFindLabel(LVT_Config, &
         label='MERRA2 data directory:',rc=status)
!    do k=1,i
    call ESMF_ConfigGetAttribute(LVT_Config, MERRA2obs(i)%odir, &
         rc=status)
    call LVT_verify(status, 'MERRA2 data directory: not defined')

    ! Yeosang Yoon (check for using bias corrected precip (=1) or uncorrected precip (=0))
    call ESMF_ConfigFindLabel(LVT_Config, &
         label='MERRA2 use corrected total precipitation:',rc=status)
    call ESMF_ConfigGetAttribute(LVT_Config, MERRA2obs(i)%usecorr, rc=status)
    call LVT_verify(status, 'MERRA2 use corrected total precipitation: not defined')

    ! Yeosang Yoon (check for using bias corrected precip (=1) or uncorrected precip (=0))
    call ESMF_ConfigFindLabel(LVT_Config, &
         label='MERRA2 apply precip scaling factors:',rc=status)
    call ESMF_ConfigGetAttribute(LVT_Config, MERRA2obs(i)%usescalef, rc=status)
    call LVT_verify(status, 'MERRA2 apply precip scaling factors: not defined')
   
    merra2obs(i)%gridDesc = 0
        
    merra2obs(i)%nc = 576
    merra2obs(i)%nr = 361

    !filling the items needed by the interpolation library
    merra2obs(i)%gridDesc(1) = 0  
    merra2obs(i)%gridDesc(2) = merra2obs(i)%nc
    merra2obs(i)%gridDesc(3) = merra2obs(i)%nr
    merra2obs(i)%gridDesc(4) = -90.000
    merra2obs(i)%gridDesc(5) = -180.000
    merra2obs(i)%gridDesc(7) = 90.000
    merra2obs(i)%gridDesc(8) = 179.375
    merra2obs(i)%gridDesc(6) = 128
    merra2obs(i)%gridDesc(9) = 0.625
    merra2obs(i)%gridDesc(10) = 0.5
    merra2obs(i)%gridDesc(20) = 0

    merra2obs(i)%datares  = 0.625

    if(MERRA2obs(i)%usescalef.eq.1) then 
       allocate(merra2obs(i)%refxrange(merra2obs(i)%nc,merra2obs(i)%nr,1,nbins))
       allocate(merra2obs(i)%refcdf(merra2obs(i)%nc,merra2obs(i)%nr,1,nbins))
       allocate(merra2obs(i)%merraxrange(merra2obs(i)%nc,merra2obs(i)%nr,1,nbins))
       allocate(merra2obs(i)%merracdf(merra2obs(i)%nc,merra2obs(i)%nr,1,nbins))

       call ESMF_ConfigFindLabel(LVT_Config, &
            label='MERRA2 precip scaling factor input file:',rc=status)
       call ESMF_ConfigGetAttribute(LVT_Config, merra2obs(i)%scaleffile, rc=status)
       call LVT_verify(status, 'MERRA2 precip scaling factor input file: not defined')

    endif

    if(LVT_isAtAfinerResolution(merra2obs(i)%datares)) then
       
       allocate(merra2obs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(merra2obs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(merra2obs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(merra2obs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            merra2obs(i)%rlat, &
            merra2obs(i)%rlon, &
            merra2obs(i)%n11)
    else
       allocate(merra2obs(i)%n11(merra2obs(i)%nc*merra2obs(i)%nr))
       call upscaleByAveraging_input(merra2obs(i)%gridDesc,&
            LVT_rc%gridDesc,merra2obs(i)%nc*merra2obs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,merra2obs(i)%n11)
    endif


    call ESMF_TimeIntervalSet(merra2obs(i)%ts, s=3600,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')

    call LVT_update_timestep(LVT_rc, 3600)

    merra2obs(i)%da = -1
    merra2obs(i)%startFlag = .true.

    allocate(merra2obs(i)%qs(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%qsb(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%swnet(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%qle(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%qh(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%frsno(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%snod(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%swe(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%qg(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%sfsm(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%rzsm(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%prcp(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%tair(LVT_rc%lnc,LVT_rc%lnr,24))
    allocate(merra2obs(i)%tskin(LVT_rc%lnc,LVT_rc%lnr,24))

    merra2obs(i)%qs = LVT_rc%udef
    merra2obs(i)%qsb = LVT_rc%udef
    merra2obs(i)%swnet = LVT_rc%udef
    merra2obs(i)%qle = LVT_rc%udef
    merra2obs(i)%qh = LVT_rc%udef
    merra2obs(i)%frsno = LVT_rc%udef
    merra2obs(i)%snod = LVT_rc%udef
    merra2obs(i)%swe = LVT_rc%udef
    merra2obs(i)%qg = LVT_rc%udef
    merra2obs(i)%sfsm = LVT_rc%udef
    merra2obs(i)%rzsm = LVT_rc%udef
    merra2obs(i)%prcp = LVT_rc%udef
    merra2obs(i)%tair = LVT_rc%udef
    merra2obs(i)%tskin = LVT_rc%udef
    
  end subroutine MERRA2obsinit


end module MERRA2obsMod
