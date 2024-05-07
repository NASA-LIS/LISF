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
! !MODULE: ISCCP_Tskin_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle land surface temperature retrievals from International 
!   Satellite Cloud Climatology Project (ISCCP; http://isccp.giss.nasa.gov,
!   Rossow and Schiffer 1991, 1992). The product used here includes
!   LST retrievals from Geostationary Operational Environmental 
!   Satellites (GOES) series, the European Meteosat series and the 
!   Japanese Geostationary Meteorological Satellite (GMS) series. 
!   The data aggregated to a global latitude-longitude grid at 1 degree 
!   resolution is used here for data assimilation 
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
! 
module ISCCP_Tskin_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ISCCP_Tskin_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: isccp_tskin_struc

  type, public ::  isccp_tskin_dec
     integer             :: mi
     integer, allocatable    :: n11(:)

!for scaling runs
     integer             :: scal
     character(len=LIS_CONST_PATH_LEN) :: modelmean
     character(len=LIS_CONST_PATH_LEN) :: modelstd
     character(len=LIS_CONST_PATH_LEN) :: obsmean
     character(len=LIS_CONST_PATH_LEN) :: obsstd
  end type isccp_tskin_dec

  type(isccp_tskin_dec), allocatable :: isccp_tskin_struc(:)

contains
!BOP
! 
! !ROUTINE: ISCCP_Tskin_setup
! \label{ISCCP_Tskin_setup}
! 
! !INTERFACE: 
  subroutine ISCCP_Tskin_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod
    use LIS_logMod, only : LIS_logunit, LIS_verify
    use LIS_perturbMod
    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data structures required for ISCCP Tskin
!   assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n 
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  synsmobsdir
    character*100          ::  temp
    real, pointer          :: obs_temp(:,:)
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40           ::  vname(1)
    real                   ::  varmin(1)
    real                   ::  varmax(1)
    type(pert_dec_type)    :: obs_pert
    real                   :: gridDesci(LIS_rc%nnest,50)

    allocate(isccp_tskin_struc(LIS_rc%nnest))
    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,synsmobsdir,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            synsmobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin scale data:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,isccp_tskin_struc(n)%scal,&
            rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin model mean data file:",&
         rc=status)
    call LIS_verify(status, 'ISCCP Tskin model mean data file: not defined')
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,isccp_tskin_struc(n)%modelmean,&
            rc=status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin model std data file:",&
         rc=status)
    call LIS_verify(status, 'ISCCP Tskin model std data file: not defined')
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,isccp_tskin_struc(n)%modelstd,&
            rc=status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin obs mean data file:",&
         rc=status)
    call LIS_verify(status, 'ISCCP Tskin obs mean data file: not defined')
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,isccp_tskin_struc(n)%obsmean,&
            rc=status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin obs std data file:",&
         rc=status)
    call LIS_verify(status, 'ISCCP Tskin obs std data file: not defined')
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,isccp_tskin_struc(n)%obsstd,&
            rc=status)
    enddo

    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
            -99.0, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
            LIS_rc%ngrid(n),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*)'read ISCCP Tskin data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_vecGrid(n),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)
       
!Perturbations State
       call LIS_readDAObsAttributes(k,vname,varmin,varmax)
              
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 

          allocate(obs_pert%vname(1))
          allocate(obs_pert%perttype(1))
          allocate(obs_pert%ssdev(1))
          allocate(obs_pert%stdmax(1))
          allocate(obs_pert%zeromean(1))
          allocate(obs_pert%tcorr(1))
          allocate(obs_pert%xcorr(1))
          allocate(obs_pert%ycorr(1))
          allocate(obs_pert%ccorr(1,1))

          call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
               obs_pert)

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_ensOnGrid(n),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               obs_pert%ssdev(1), rc=status)
          call LIS_verify(status)

          call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
               obs_pert%stdmax(1), rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
               obs_pert%zeromean(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
               obs_pert%tcorr(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
               obs_pert%xcorr(1),rc=status)
          
          call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
               obs_pert%ycorr(1),rc=status)

          call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
               obs_pert%ccorr(1,:),itemCount=1,rc=status)

          call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
          call LIS_verify(status)
       endif
          
!to obsarr

    enddo
!-------------------------------------------------------------
! setting up interpolation weights for neighbor search
!-------------------------------------------------------------
    gridDesci = 0 

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0 
       gridDesci(n,2) = 360
       gridDesci(n,3) = 180
       gridDesci(n,4) = -89.50
       gridDesci(n,5) = -179.50
       gridDesci(n,6) = 128
       gridDesci(n,7) = 89.50
       gridDesci(n,8) = 179.50
       gridDesci(n,9) = 1.000
       gridDesci(n,10) = 1.000
       gridDesci(n,20) = 64.0

       isccp_tskin_struc(n)%mi = 360*180
       allocate(isccp_tskin_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       call neighbor_interp_input(n,gridDesci(n,:),&
            isccp_tskin_struc(n)%n11)
    enddo

    write(LIS_logunit,*) 'Created the States to hold the observations data'
    
  end subroutine ISCCP_Tskin_setup
  
end module ISCCP_Tskin_module
