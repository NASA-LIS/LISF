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
! !MODULE: UAsnow_obsMod
! 
! !DESCRIPTION: 
!  This module handles the observation plugin for the
!  University of Arizona (UA) SWE/Snow Depth data v01
!  for use within the LIS OPT/UE subsystem. 
!  The UA SNOW data is provided in the NAD 1983 grid
!  with ~4-km resolution. The domain extents are from
!  approximately (24N, -125W) to (50N, -66.5W).
!  The data entries are 16-bit signed integers.
!  
!  https://nsidc.org/data/nsidc-0719
!
!   
! !REVISION HISTORY: 
!  2 May 2020 Sujay Kumar;   Initial Specification 
! 
module UAsnow_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: UAsnow_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: UAsnow_obs_struc

  type, public ::  UAsnow_obs_data_dec

     integer                    :: nc
     integer                    :: nr
     type(ESMF_Time)            :: starttime
     type(ESMF_TimeInterval)    :: timestep

     integer, allocatable       :: n11(:)
     integer, allocatable       :: n12(:)
     integer, allocatable       :: n21(:)
     integer, allocatable       :: n22(:)

     real, allocatable       :: w11(:)
     real, allocatable       :: w12(:)
     real, allocatable       :: w21(:)
     real, allocatable       :: w22(:)

  end type UAsnow_obs_data_dec

  type(UAsnow_obs_data_dec), allocatable :: UAsnow_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: UAsnow_obs_setup
! \label{UAsnow_obs_setup}
! 
! !INTERFACE: 
  subroutine UAsnow_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use map_utils
    use LIS_timeMgrMod
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!  This routines completes the setup of the UA snow plugin
!  for OPTUE. This includes the definition of the data grid and the 
!  setup of the interpolation weights. 
!   
!   The arguments are: 
!   \begin{description}
!    \item[Obs\_State]   observation state object 
!   \end{description}
!EOP
    integer                   ::  n 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField1,obsField2
    character(len=LIS_CONST_PATH_LEN) ::  obsdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                 :: k
    integer                 :: ftn
    integer                 :: status
    real                   :: gridDesci(50)


    allocate(UAsnow_obs_struc(LIS_rc%nnest))

    write(LIS_logunit,*) '[INFO] Setting up UA snow data reader....'

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status, 'Error in ESMF_ArraySpecSet')
    
    call ESMF_ConfigFindLabel(LIS_config,"UA snow data directory:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,obsdir,&
         rc=status)
    call LIS_verify(status, 'UA snow data directory: not defined')
    
    do n=1,LIS_rc%nnest

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status, 'Error in ESMF_AttributeSet: Data Directory')

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status, 'Error in ESMF_AttributeSet: Data Update Status')
    enddo

    do n=1,LIS_rc%nnest
!----------------------------------------------------------------------
! Describes the native grid projection, spatial extent, and resolution
! of the UA snow data  (lat/lon, CONUS domain, at ~4km resolution)
!----------------------------------------------------------------------
       gridDesci = 0
       gridDesci(1) = 0
       gridDesci(2) = 1405
       gridDesci(3) = 621
       gridDesci(4) = 24.0833340
       gridDesci(5) = -125.0000
       gridDesci(6) = 128
       gridDesci(7) = 49.9166679
       gridDesci(8) = -66.5000
       gridDesci(9) = 0.04166662697178698
       gridDesci(10) = 0.04166662697178698
       gridDesci(20) = 64
       
       UAsnow_obs_struc(n)%nc = 1405
       UAsnow_obs_struc(n)%nr = 621
       

       allocate(UAsnow_obs_struc(n)%n11(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       allocate(UAsnow_obs_struc(n)%n12(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       allocate(UAsnow_obs_struc(n)%n21(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       allocate(UAsnow_obs_struc(n)%n22(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       allocate(UAsnow_obs_struc(n)%w11(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       allocate(UAsnow_obs_struc(n)%w12(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       allocate(UAsnow_obs_struc(n)%w21(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       allocate(UAsnow_obs_struc(n)%w22(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))

       call bilinear_interp_input (n, gridDesci,&
            UAsnow_obs_struc(n)%n11,&
            UAsnow_obs_struc(n)%n12,&
            UAsnow_obs_struc(n)%n21,&
            UAsnow_obs_struc(n)%n22,&
            UAsnow_obs_struc(n)%w11,&
            UAsnow_obs_struc(n)%w12,&
            UAsnow_obs_struc(n)%w21,&
            UAsnow_obs_struc(n)%w22)
      
       obsField1 = ESMF_FieldCreate(arrayspec=realarrspec, &
            grid=LIS_vecGrid(n), &
            name="UA_SWE", rc=status)
       call LIS_verify(status, 'Error in ESMF_FieldCreate: UA_SWE ')

       obsField2 = ESMF_FieldCreate(arrayspec=realarrspec, &
            grid=LIS_vecGrid(n), &
            name="UA_SNOD", rc=status)
       call LIS_verify(status, 'Error in ESMF_FieldCreate: UA_SNOD')

       call ESMF_StateAdd(Obs_State(n),(/obsField1/),rc=status)
       call LIS_verify(status, 'Error in ESMF_StateAdd: obsField1')

       call ESMF_StateAdd(Obs_State(n),(/obsField2/),rc=status)
       call LIS_verify(status, 'Error in ESMF_StateAdd: obsField2')

       call ESMF_TimeIntervalSet(UAsnow_obs_struc(n)%timestep, &
            s=86400,rc=status)
       call LIS_verify(status, 'Error in ESMF_TimeIntervalSet: UAsnow_obs_struc(n)%timestep')
    enddo

    write(LIS_logunit,*) '[INFO] created the States to hold the UA snow data'
    
  end subroutine UAsnow_obs_setup
  
end module UAsnow_obsMod
