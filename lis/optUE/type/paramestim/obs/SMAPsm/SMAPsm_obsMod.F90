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
! !MODULE: SMAPsm_obsMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  21 Sep 2018 Sujay Kumar;   Initial Specification 
! 
module SMAPsm_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SMAPsm_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAPsm_obs_struc

  type, public ::  SMAPsm_obs_data_dec

     integer                    :: nc
     integer                    :: nr
     real,     allocatable      :: smobs(:,:)
     real,     allocatable      :: smtime(:,:)

     integer, allocatable       :: n11(:)

     character*20               :: data_designation
     logical                    :: startMode

  end type SMAPsm_obs_data_dec

  type(SMAPsm_obs_data_dec), allocatable :: SMAPsm_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: SMAPsm_obs_setup
! \label{SMAPsm_obs_setup}
! 
! !INTERFACE: 
  subroutine SMAPsm_obs_setup(Obs_State)
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
!   The arguments are: 
!   \begin{description}
!    \item[Obs\_State]   observation state object 
!   \end{description}
!EOP
    integer                   ::  n 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField
    character(len=LIS_CONST_PATH_LEN) ::  obsdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                 :: k
    integer                 :: ftn
    integer                 :: status
    real                   :: gridDesci(50)


    allocate(SMAPsm_obs_struc(LIS_rc%nnest))

    write(LIS_logunit,*) '[INFO] Setting up SMAP sm data reader....'

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"SMAP soil moisture data directory:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,obsdir,&
         rc=status)
    call LIS_verify(status, 'SMAP soil moisture data directory: not defined')
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "SMAP soil moisture data designation:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            SMAPsm_obs_struc(n)%data_designation,&
            rc=status)
       call LIS_verify(status, 'SMAP soil moisture data designation: is missing')
    enddo
    do n=1,LIS_rc%nnest

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    do n=1,LIS_rc%nnest
       if(SMAPsm_obs_struc(n)%data_designation.eq."SPL3SMP") then 
          gridDesci = 0
          gridDesci(1) = 9
          gridDesci(2) = 964
          gridDesci(3) = 406
          gridDesci(9) = 4 !M36 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.36 
          gridDesci(11) = 1 !for the global switch

          SMAPsm_obs_struc(n)%nc = 964
          SMAPsm_obs_struc(n)%nr = 406

       elseif(SMAPsm_obs_struc(n)%data_designation.eq."SPL3SMP_E") then 
          gridDesci = 0
          gridDesci(1) = 9
          gridDesci(2) = 3856
          gridDesci(3) = 1624
          gridDesci(9) = 5 !M09 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.09 
          gridDesci(11) = 1 !for the global switch

          SMAPsm_obs_struc(n)%nc = 3856
          SMAPsm_obs_struc(n)%nr = 1624

       endif

       allocate(SMAPsm_obs_struc(n)%smobs(LIS_rc%lnc(n), &
            LIS_rc%lnr(n)))
       SMAPsm_obs_struc(n)%smobs = LIS_rc%udef
       allocate(SMAPsm_obs_struc(n)%smtime(LIS_rc%lnc(n), &
            LIS_rc%lnr(n)))


       allocate(SMAPsm_obs_struc(n)%n11(LIS_rc%lnc(n)*&
            LIS_rc%lnr(n)))
       
       call neighbor_interp_input(n, gridDesci,&
            SMAPsm_obs_struc(n)%n11)
       
       call LIS_registerAlarm("SMAP data read alarm",&
            86400.0, 86400.0)

       SMAPsm_obs_struc(n)%startMode = .true. 

       obsField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecGrid(n), &
            name="SMAP_sm", rc=status)
       call LIS_verify(status)
       
       call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
       call LIS_verify(status)

    enddo

    write(LIS_logunit,*) '[INFO] created the States to hold the SMAP soil moisture data'
    
  end subroutine SMAPsm_obs_setup
  
end module SMAPsm_obsMod
