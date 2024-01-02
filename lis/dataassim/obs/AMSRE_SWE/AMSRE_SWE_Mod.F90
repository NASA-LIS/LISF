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
! !MODULE: AMSRE_SWE_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle the Level 3 AMSR-E SWE retrievals. The documentation of the 
!   data can be found at : 
!   
!    http://nsidc.org/data/docs/daac/ae\_swe\_ease-grids.gd.html
!   
! !REVISION HISTORY: 
!  01Jul10   Sujay Kumar; Initial Specification
! 
module AMSRE_SWE_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: AMSRE_SWE_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: AMSRE_SWE_struc
!EOP
  type, public:: AMSRE_SWE_dec

     integer         :: mo
     integer,allocatable :: n112(:,:)
     real,allocatable    :: rlat2(:,:)
     real,allocatable    :: rlon2(:,:)
     real            :: gridDesci(50,2)

     integer             :: ihemi, nhemi
     real,    allocatable    :: sweobs(:)
     real,    allocatable    :: sweqc(:)

  end type AMSRE_SWE_dec
  
  type(AMSRE_SWE_dec),allocatable :: AMSRE_SWE_struc(:)
  
contains

!BOP
! 
! !ROUTINE: AMSRE_SWE_setup
! \label{AMSRE_SWE_setup}
! 
! !INTERFACE: 
  subroutine AMSRE_SWE_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_perturbMod
    use LIS_logmod
    use LIS_DAobservationsMod
    use map_utils

    implicit none 

! !ARGUMENTS: 
    integer                ::  k
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for the AMSR-E SWE
!   observation plugin for data assimilation. The routine also 
!   sets up interpolation weights (using neighbor search) for 
!   reprojecting the AMSR-E data from 2-hemisphere EASE grids
!   to the LIS grid. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer, parameter     ::  nbins = 1000 
    integer                ::  n
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  amsresweobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40           ::  vname(1)
    real                   ::  varmin(1)
    real                   ::  varmax(1)
    real                   ::  ssdev(1)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real                   ::  max_lat, min_lat, rlat, rlon
    integer                ::  c,r

    allocate(AMSRE_SWE_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"AMSR-E SWE data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,amsresweobsdir,&
            rc=status)
       call LIS_verify(status, 'AMSR-E SWE data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            trim(amsresweobsdir), rc=status)
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

    write(LIS_logunit,*)'read AMSR-E SWE data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. amsr-e 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_vecGrid(n),&
            name="Observation"//vid(1)//vid(2),rc=status)
       call LIS_verify(status)

!Perturbations State
       call LIS_readDAObsAttributes(k,vname,varmin,varmax)
       
       allocate(obsstd(LIS_rc%ngrid(n)))
       obsstd = 0 ! ssdev(1)
       call ESMF_AttributeSet(obsField(n),"Standard Deviation",&
            obsstd, itemCount=LIS_rc%ngrid(n),rc=status)
       call LIS_verify(status)
       deallocate(obsstd)
       
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
               grid=LIS_ensOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,&
               farrayPtr=obs_temp,rc=status)
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


    write(LIS_logunit,*) 'Created the States to hold the AMSR-E observations data'
!--------------------------------------------------------------------------------
! figure out the grid span and whether we need to process both southern 
! and northern hemispheres
!--------------------------------------------------------------------------------
    do n=1, LIS_rc%nnest
       max_lat = -10000
       min_lat = 10000 
       
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
             call ij_to_latlon(LIS_domain(n)%lisproj,float(c),float(r),&
                  rlat,rlon)
             
             if(rlat.gt.max_lat) max_lat = rlat
             if(rlat.lt.min_lat) min_lat = rlat
          enddo
       enddo

       if(max_lat.gt.0.and.min_lat.lt.0) then ! domain split in 2 hemispheres
          AMSRE_SWE_struc(n)%ihemi = 1
          AMSRE_SWE_struc(n)%nhemi = 2
       elseif(max_lat.ge.0.and.min_lat.ge.0) then !all northern hemisphere
          AMSRE_SWE_struc(n)%ihemi = 1
          AMSRE_SWE_struc(n)%nhemi = 1
       else !all southern hemisphere
          AMSRE_SWE_struc(n)%ihemi = 2
          AMSRE_SWE_struc(n)%nhemi = 2
       endif

    enddo

!--------------------------------------------------------------------------------
! set up interpolation weights for neighbor search
!--------------------------------------------------------------------------------
    call computeInterpWeights()

    do n=1,LIS_rc%nnest
       allocate(AMSRE_SWE_struc(n)%sweobs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(AMSRE_SWE_struc(n)%sweqc(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       AMSRE_SWE_struc(n)%sweobs = LIS_rc%udef
       AMSRE_SWE_struc(n)%sweqc = LIS_rc%udef

       call LIS_registerAlarm("AMSR-E SWE read alarm", 86400.0, 86400.0)

    enddo
    
  end subroutine AMSRE_SWE_setup

!BOP
! !ROUTINE: computeInterpWeights
! \label{AMSRE_SWE_computeInterpWeights}
! 
! !INTERFACE: 
  subroutine computeInterpWeights()
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config
    use LIS_logMod, only : LIS_logunit

    implicit none

! 
! !DESCRIPTION: 
!   This subroutine sets up the interpolation weights to transform the 
!   AMSRE data in EASE grid to the grid  used in the LIS simulation. 
!   The code employs a neighbor-based interpolation
! 
!EOP
    integer,parameter  :: ease_nr=721
    integer,parameter  :: ease_nc=721
    integer         :: npts,n
    integer         :: hemi

    !read the domain specs from LIS_rc% struct

    do n =1, LIS_rc%nnest

       npts= LIS_rc%lnc(n)*LIS_rc%lnr(n)
       AMSRE_SWE_struc(n)%mo=npts
       allocate(AMSRE_SWE_struc(n)%rlat2(npts,2))
       allocate(AMSRE_SWE_struc(n)%rlon2(npts,2))
       allocate(AMSRE_SWE_struc(n)%n112(npts,2))

       do hemi=AMSRE_SWE_struc(n)%ihemi, AMSRE_SWE_struc(n)%nhemi 
          !initialize the entire array
          AMSRE_SWE_struc(n)%gridDesci(:,hemi) =0.0 
          
      !filling the items needed by the interpolation library
          AMSRE_SWE_struc(n)%gridDesci(1,hemi) = 9  !input is EASE grid
          !these  corner coordinates were calculated based on ezlh_convert
          AMSRE_SWE_struc(n)%gridDesci(4,hemi) = -90.0  !lat
          AMSRE_SWE_struc(n)%gridDesci(5,hemi) = -179.6096 !lon
          AMSRE_SWE_struc(n)%gridDesci(7,hemi) = 83.33788  !lat
          AMSRE_SWE_struc(n)%gridDesci(8,hemi) = 180.1301  !lon


          AMSRE_SWE_struc(n)%gridDesci(2,hemi) = ease_nc  !nx
          AMSRE_SWE_struc(n)%gridDesci(3,hemi) = ease_nr  !ny
          if(hemi.eq.1) then 
             AMSRE_SWE_struc(n)%gridDesci(9,hemi)  = 2  !Northern hemi
          else
             AMSRE_SWE_struc(n)%gridDesci(9,hemi)  = 3  !Southern hemi
          endif
          
          AMSRE_SWE_struc(n)%rlat2(:,hemi)=0.0
          AMSRE_SWE_struc(n)%rlon2(:,hemi)=0.0
          AMSRE_SWE_struc(n)%n112(:,hemi)=0.0
          call neighbor_interp_input_withgrid(AMSRE_SWE_struc(n)%gridDesci(:,hemi),&
               LIS_rc%gridDesc(n,:),&
               npts,AMSRE_SWE_struc(n)%rlat2(:,hemi),&
               AMSRE_SWE_struc(n)%rlon2(:,hemi),AMSRE_SWE_struc(n)%n112(:,hemi))
       enddo
       
    end do
      
  end subroutine computeInterpWeights

end module AMSRE_SWE_Mod
