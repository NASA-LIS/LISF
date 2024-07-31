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
! !MODULE: SMMRSNWDsnow_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle SMMR snow depth retrievals. 
!   
! !REVISION HISTORY: 
!  16 Oct 2012   Sujay Kumar;   Initial Specification
! 
module SMMRSNWDsnow_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SMMRSNWDsnow_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMMRSNWDsnow_struc

  type, public ::  SMMRSNWDsnow_dec
     logical             :: startMode
     real                :: ssdev
     integer             :: snowfield
     integer             :: mi
     integer             :: nc,nr

     integer, allocatable    :: n11(:)
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     real,    allocatable    :: snwd(:)
     real,    allocatable    :: snwdtime(:,:)
  end type SMMRSNWDsnow_dec

  type(SMMRSNWDsnow_dec), allocatable :: SMMRSNWDsnow_struc(:)

contains
!BOP
! 
! !ROUTINE: SMMRSNWDsnow_setup
! \label{SMMRSNWDsnow_setup}
! 
! !INTERFACE: 
  subroutine SMMRSNWDsnow_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_logMod
   
    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for SMMR snow depth data. 
!
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  smmrsnowobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real,         allocatable  ::  varmin(:)
    real,         allocatable  ::  varmax(:)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    real                   :: gridDesci(LIS_rc%nnest,50)
    real, allocatable          :: ssdev(:)
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2    
    real                   :: minlat, minlon, maxlon, maxlat, dx, dy

    allocate(SMMRSNWDsnow_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SMMR snow depth data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,smmrsnowobsdir,&
            rc=status)
       call LIS_verify(status,'SMMR snow depth data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            smmrsnowobsdir, rc=status)
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
            LIS_rc%obs_ngrid(k),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*)'[INFO] read SMMR snow depth data specifications'

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

       obsField(n) = ESMF_FieldCreate(grid=LIS_obsvecGrid(n,k),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)
       
!Perturbations State
       write(LIS_logunit,*) '[INFO] Opening attributes for observations ',&
            trim(LIS_rc%obsattribfile(k))
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
       read(ftn,*)
       read(ftn,*) LIS_rc%nobtypes(k)
       read(ftn,*)
    
       allocate(vname(LIS_rc%nobtypes(k)))
       allocate(varmax(LIS_rc%nobtypes(k)))
       allocate(varmin(LIS_rc%nobtypes(k)))
       
       do i=1,LIS_rc%nobtypes(k)
          read(ftn,fmt='(a40)') vname(i)
          read(ftn,*) varmin(i),varmax(i)
          write(LIS_logunit,*) '[INFO] ',vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)  
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       allocate(ssdev(LIS_rc%obs_ngrid(k)))

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

          ssdev = obs_pert%ssdev(1)
          SMMRSNWDsnow_struc(n)%ssdev =obs_pert%ssdev(1) 

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status)
          
          if(LIS_rc%obs_ngrid(k).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
             call LIS_verify(status)
          endif

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
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)
    enddo
!-------------------------------------------------------------
! set up the SMMR domain and interpolation weights. 
!-------------------------------------------------------------
    gridDesci = 0 
    do n=1,LIS_rc%nnest
       SMMRSNWDsnow_struc(n)%nc = 1440
       SMMRSNWDsnow_struc(n)%nr = 360 

       gridDesci(n,1) = 0 
       gridDesci(n,2) = SMMRSNWDsnow_struc(n)%nc
       gridDesci(n,3) = SMMRSNWDsnow_struc(n)%nr
       gridDesci(n,4) = 0.125
       gridDesci(n,5) = -179.875
       gridDesci(n,6) = 128
       gridDesci(n,7) = 89.875
       gridDesci(n,8) = 179.875
       gridDesci(n,9) = 0.25
       gridDesci(n,10) = 0.25
       gridDesci(n,20) = 64.0

       SMMRSNWDsnow_struc(n)%mi = &
            SMMRSNWDsnow_struc(n)%nc*SMMRSNWDsnow_struc(n)%nr

       allocate(SMMRSNWDsnow_struc(n)%n11(SMMRSNWDsnow_struc(n)%mi))
       allocate(SMMRSNWDsnow_struc(n)%rlat(SMMRSNWDsnow_struc(n)%mi))
       allocate(SMMRSNWDsnow_struc(n)%rlon(SMMRSNWDsnow_struc(n)%mi))

       allocate(SMMRSNWDsnow_struc(n)%snwd(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(SMMRSNWDsnow_struc(n)%snwdtime(&
            LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k)))
       SMMRSNWDsnow_struc(n)%snwd = LIS_rc%udef
       SMMRSNWDsnow_struc(n)%snwdtime = -1

       call neighbor_interp_input_withgrid(gridDesci(n,:), & 
            LIS_rc%obs_gridDesc(k,:),&
            SMMRSNWDsnow_struc(n)%nc*SMMRSNWDsnow_struc(n)%nr,&
            SMMRSNWDsnow_struc(n)%rlat, SMMRSNWDsnow_struc(n)%rlon,&
            SMMRSNWDsnow_struc(n)%n11)

       call LIS_registerAlarm("SMMR snow depth read alarm",&
            86400.0, 86400.0)

       SMMRSNWDsnow_struc(n)%startMode = .true. 
    enddo

    write(LIS_logunit,*) '[INFO] Created ESMF States to hold SMMR observations data'

  end subroutine SMMRSNWDsnow_setup
  
end module SMMRSNWDsnow_Mod

