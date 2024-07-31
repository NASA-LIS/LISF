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
! !MODULE: multisynsmobs_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle synthetic soil moisture observations generated 
!   from a previous LIS-Noah simulation. In this example, 
!   synthetic observations from multiple soil layers 
!   (instead of just the surface soil moisture) are 
!   processed and demonstrates the use of observation 
!   interfaces for multiple observation types. 
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
! 
module multisynsmobs_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: multisynsmobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------


contains
!BOP
! 
! !ROUTINE: multisynsmobs_setup
! \label{multisynsmobs_setup}
! 
! !INTERFACE: 
  subroutine multisynsmobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid, LIS_ensOnGrid
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
!   creation of data strctures required for handling
!   multi-layer soil moisture observations
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n 
    integer                ::  m
    integer                ::  status
    type(ESMF_Field)       ::  obsField
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  synsmobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    type(pert_dec_type)    ::  obs_pert
    character*1            ::  vid(2)
    character*40           ::  vname(LIS_rc%nobtypes(k))
    real                   ::  varmin(LIS_rc%nobtypes(k))
    real                   ::  varmax(LIS_rc%nobtypes(k))
    real, pointer          ::  obs_temp(:,:)

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"Synthetic multi-sm data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synsmobsdir,&
            rc=status)
       call LIS_verify(status,'Synthetic multi-sm data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            synsmobsdir, rc=status)
       call LIS_verify(status)
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

    write(LIS_logunit,*)'read Synthetic soil moisture data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
!Perturbations State
       call LIS_readDAObsAttributes(k,vname,varmin,varmax)

       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
          allocate(obs_pert%vname(LIS_rc%nobtypes(k)))
          allocate(obs_pert%perttype(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ssdev(LIS_rc%nobtypes(k)))
          allocate(obs_pert%stdmax(LIS_rc%nobtypes(k)))
          allocate(obs_pert%zeromean(LIS_rc%nobtypes(k)))
          allocate(obs_pert%tcorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%xcorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ycorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ccorr(LIS_rc%nobtypes(k),LIS_rc%nobtypes(k)))

          call LIS_readPertAttributes(LIS_rc%nobtypes(k), &
               LIS_rc%obspertAttribFile(k),obs_pert)
       endif

       do m=1,LIS_rc%nobtypes(k)

          write(unit=temp,fmt='(i2.2)') m
          read(unit=temp,fmt='(2a1)') vid

          obsField = ESMF_FieldCreate(arrayspec=realarrspec,grid=LIS_vecGrid(n),&
               name="Observation"//vid(1)//vid(2), &
               rc=status)
          call LIS_verify(status, 'error in obsField create')
                     
          call ESMF_StateAdd(OBS_State(n),(/obsField/),rc=status)
          call LIS_verify(status, 'field add to obs_state')

          if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
             
             pertField = ESMF_FieldCreate(arrayspec=pertArrSpec,&
                  grid=LIS_ensOnGrid(n),name="Observation"&
                  //vid(1)//vid(2), rc=status)
             call LIS_verify(status, 'pert field create error')
             
! initializing the perturbations to be zero 
             call ESMF_FieldGet(pertField,localDE=0,farrayPtr=obs_temp,rc=status)
             call LIS_verify(status)
             obs_temp(:,:) = 0 

             call ESMF_AttributeSet(pertField,"Perturbation Type",obs_pert%perttype(m),&
                  rc=status)
             call LIS_verify(status, 'Perturbation Type attribute add to pertField')

             call ESMF_AttributeSet(pertField,"Standard Deviation",obs_pert%ssdev(m),&
                  rc=status)
             call LIS_verify(status, 'Standard Deviation attribute add to pertField')

             call ESMF_AttributeSet(pertField,"Std Normal Max",obs_pert%stdmax(m),&
                  rc=status)
             call LIS_verify(status, 'Standard Normal Max attribute add to pertField')
             
             call ESMF_AttributeSet(pertField,"Ensure Zero Mean",obs_pert%zeromean(m),&
               rc=status)
             call LIS_verify(status, 'Ensure Zero Mean attribute add to pertField')
             
             call ESMF_AttributeSet(pertField,"Temporal Correlation Scale",&
                  obs_pert%tcorr(m),rc=status)
             call LIS_verify(status, 'Temporal Correlation Scale attribute add to pertField')
          
             call ESMF_AttributeSet(pertField,"X Correlation Scale",&
                  obs_pert%xcorr(m),rc=status)
             
             call ESMF_AttributeSet(pertField,"Y Correlation Scale",&
                  obs_pert%ycorr(m),rc=status)
             
             call ESMF_AttributeSet(pertField,"Cross Correlation Strength",&
                  obs_pert%ccorr(m,:),itemCount=LIS_rc%nobtypes(k),rc=status)
             
             call ESMF_StateAdd(OBS_Pert_State(n),(/pertField/),rc=status)
             call LIS_verify(status, 'error in pertfield add to obs_pert_state')
          endif
       enddo
    enddo
   

    write(LIS_logunit,*) 'Created the States to hold the observations data'
    
  end subroutine multisynsmobs_setup
  
end module multisynsmobs_Mod
