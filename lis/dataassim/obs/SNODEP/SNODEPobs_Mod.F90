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
! !MODULE: SNODEPobs_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle SNODEP (SWE) observations 
!   
! !REVISION HISTORY: 
!  21Aug2008: Sujay Kumar; Initial Specification
! 
module SNODEPobs_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SNODEPobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!----------------------------------------------------------------------------- 
  PUBLIC :: SNODEP_obs_obj

!EOP
  type, public :: SNODEPobsdec
     real, allocatable      :: rlat1_nh(:)
     real, allocatable      :: rlat1_sh(:)
     real, allocatable      :: rlon1_nh(:)
     real, allocatable      :: rlon1_sh(:)
     integer, allocatable   :: n111_nh(:)
     integer, allocatable   :: n111_sh(:)
     integer, allocatable   :: n121_nh(:)
     integer, allocatable   :: n121_sh(:)
     integer, allocatable   :: n211_nh(:)
     integer, allocatable   :: n211_sh(:)
     integer, allocatable   :: n221_nh(:)
     integer, allocatable   :: n221_sh(:)
     real, allocatable      :: w111_nh(:)
     real, allocatable      :: w121_nh(:)
     real, allocatable      :: w111_sh(:)
     real, allocatable      :: w121_sh(:)
     real, allocatable      :: w211_nh(:)
     real, allocatable      :: w221_nh(:)
     real, allocatable      :: w211_sh(:)
     real, allocatable      :: w221_sh(:)     
     integer                :: gridspan
     integer                :: shemi
     integer                :: nhemi
     integer                :: mo1
     integer                :: mo2
     integer                :: hemi_nc(2)
     integer                :: hemi_nr(2)
     integer                :: mi
     integer                :: pmax
     integer                :: mesh
     character(len=5)       :: conv
    real                   ::  gridDesco(2,50)

  end type SNODEPobsdec

  type(SNODEPobsdec), allocatable :: SNODEP_obs_obj(:)
contains

!BOP
! 
! !ROUTINE: SNODEPobs_setup
! \label{SNODEPobs_setup}
! 
! !INTERFACE: 
  subroutine SNODEPobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_DAobservationsMod
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
!   creation of data strctures required for SNODEP assimilation
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
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  SNODEPobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40,allocatable   ::  vname(:)
    real  , allocatable        ::  varmin(:)
    real  , allocatable        ::  varmax(:)
    real                   ::  gridDesci(50)
    integer                ::  i,ihemi
    real, pointer          ::  obs_temp(:,:)
    real,  allocatable         ::  ssdev(:)
    type(pert_dec_type)    ::  obs_pert
    integer                ::  ftn
    real, parameter :: xmeshl1 = 47.625
    real, parameter :: xmeshl2 = 23.812
    real, parameter :: xpnmcaf1 = 257
    real, parameter :: ypnmcaf1 = 257
    real, parameter :: xpnmcaf2 = 513
    real, parameter :: ypnmcaf2 = 513
    real :: xmesh, orient,xi1,xj1
    real :: alat1,alon1
    
    allocate(SNODEP_obs_obj(LIS_rc%nnest))
    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SNODEP data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SNODEPobsdir,&
            rc=status)
       if(status /= ESMF_SUCCESS)then
         write(LIS_logunit,*)"[ERR] SNODEP data directory is missing"
       end if
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            SNODEPobsdir, rc=status)
       call LIS_verify(status)

    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SNODEP mesh resolution:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SNODEP_obs_obj(n)%mesh,&
            rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SNODEP naming convention:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SNODEP_obs_obj(n)%conv,&
            rc=status)
       call LIS_verify(status, "SNODEP naming convention: not defined")
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

    write(LIS_logunit,*)'[INFO] read SNODEP data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. modis 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_obsVecGrid(n,k),&
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
          
          ! Set obs err to be uniform (will be rescaled later for each grid point). 
          ssdev = obs_pert%ssdev(1)
          
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
       
       if(SNODEP_obs_obj(n)%mesh.eq.8) then 
          SNODEP_obs_obj(n)%pmax = 512
          SNODEP_obs_obj(n)%mi   = 512*512
       elseif(SNODEP_obs_obj(n)%mesh.eq.16) then 
          SNODEP_obs_obj(n)%pmax = 1024
          SNODEP_obs_obj(n)%mi   = 1024*1024
       endif
       if(LIS_rc%obs_gridDesc(k,1) .eq.0) then !latlon domain
          if(LIS_rc%obs_gridDesc(k,4).ge.0.and.LIS_rc%obs_gridDesc(k,7).ge.0) then 
             SNODEP_obs_obj(n)%gridspan = 1 ! NH only
             SNODEP_obs_obj(n)%shemi = 1
             SNODEP_obs_obj(n)%nhemi = 1
          elseif(LIS_rc%obs_gridDesc(k,4).le.0.and.LIS_rc%obs_gridDesc(k,7).le.0) then 
             SNODEP_obs_obj(n)%gridspan = 2 ! SH only 
             SNODEP_obs_obj(n)%shemi = 2
             SNODEP_obs_obj(n)%nhemi = 2
          else
             SNODEP_obs_obj(n)%gridspan = 3 ! NH and SH
             SNODEP_obs_obj(n)%shemi = 1
             SNODEP_obs_obj(n)%nhemi = 2
          endif
          SNODEP_obs_obj(n)%hemi_nc = nint((LIS_rc%obs_gridDesc(k,8)-&
               LIS_rc%obs_gridDesc(k,5))&
               /LIS_rc%obs_gridDesc(k,9))+1
          if(SNODEP_obs_obj(n)%gridspan.eq.1) then 
             SNODEP_obs_obj(n)%hemi_nr(1) = nint((LIS_rc%obs_gridDesc(k,7)-&
                  LIS_rc%obs_gridDesc(k,4))/LIS_rc%obs_gridDesc(k,10))+1
             SNODEP_obs_obj(n)%hemi_nr(2) = 0 
             SNODEP_obs_obj(n)%mo1 = SNODEP_obs_obj(n)%hemi_nc(1)*&
                  SNODEP_obs_obj(n)%hemi_nr(1)
             SNODEP_obs_obj(n)%mo2 = 0 
          elseif(SNODEP_obs_obj(n)%gridspan.eq.2) then 
             SNODEP_obs_obj(n)%hemi_nr(1) = 0 
             SNODEP_obs_obj(n)%hemi_nr(2) = nint((LIS_rc%obs_gridDesc(k,7)-&
                  LIS_rc%obs_gridDesc(k,4))/LIS_rc%obs_gridDesc(k,10)+1)
             SNODEP_obs_obj(n)%mo1 = 0 
             SNODEP_obs_obj(n)%mo2 = SNODEP_obs_obj(n)%hemi_nc(2)*&
                  SNODEP_obs_obj(n)%hemi_nr(2)
          else
             SNODEP_obs_obj(n)%hemi_nr(1) = nint((LIS_rc%obs_gridDesc(k,7)-&
                  LIS_rc%obs_gridDesc(k,10)/2)/LIS_rc%obs_gridDesc(k,10)+1)
             SNODEP_obs_obj(n)%hemi_nr(2) = nint((-LIS_rc%obs_gridDesc(k,10)/2-&
                  LIS_rc%obs_gridDesc(k,4))/LIS_rc%obs_gridDesc(k,10)+1) 
             SNODEP_obs_obj(n)%mo1 = SNODEP_obs_obj(n)%hemi_nc(1)*&
                  SNODEP_obs_obj(n)%hemi_nr(1)
             SNODEP_obs_obj(n)%mo2 = SNODEP_obs_obj(n)%hemi_nc(2)*&
                  SNODEP_obs_obj(n)%hemi_nr(2)
          endif
          
          SNODEP_obs_obj(n)%gridDesco = 0 
          gridDesci = 0 
          do ihemi = SNODEP_obs_obj(n)%shemi,SNODEP_obs_obj(n)%nhemi
             SNODEP_obs_obj(n)%gridDesco(ihemi,1) = 0 
             SNODEP_obs_obj(n)%gridDesco(ihemi,2) = SNODEP_obs_obj(n)%hemi_nc(ihemi)
             SNODEP_obs_obj(n)%gridDesco(ihemi,3) = SNODEP_obs_obj(n)%hemi_nr(ihemi)
             SNODEP_obs_obj(n)%gridDesco(ihemi,5) = LIS_rc%obs_gridDesc(k,5)
             SNODEP_obs_obj(n)%gridDesco(ihemi,8) = LIS_rc%obs_gridDesc(k,8)
             SNODEP_obs_obj(n)%gridDesco(ihemi,6) = LIS_rc%obs_gridDesc(k,6)
             SNODEP_obs_obj(n)%gridDesco(ihemi,9) = LIS_rc%obs_gridDesc(k,9)
             SNODEP_obs_obj(n)%gridDesco(ihemi,10) = LIS_rc%obs_gridDesc(k,10)
             SNODEP_obs_obj(n)%gridDesco(ihemi,20) = 255
             if(SNODEP_obs_obj(n)%gridspan.eq.1.or.&
                  SNODEP_obs_obj(n)%gridspan.eq.2) then 
                SNODEP_obs_obj(n)%gridDesco(ihemi,4) = LIS_rc%obs_gridDesc(k,4)
                SNODEP_obs_obj(n)%gridDesco(ihemi,7) = LIS_rc%obs_gridDesc(k,7)
             else
                if(ihemi.eq.1) then 
                   SNODEP_obs_obj(n)%gridDesco(ihemi,4) = LIS_rc%obs_gridDesc(k,9)/2
                   SNODEP_obs_obj(n)%gridDesco(ihemi,7) = LIS_rc%obs_gridDesc(k,7)
                else
                   SNODEP_obs_obj(n)%gridDesco(ihemi,4) = LIS_rc%obs_gridDesc(k,4)
                   SNODEP_obs_obj(n)%gridDesco(ihemi,7) = -LIS_rc%obs_gridDesc(k,9)/2
                endif
             endif
             
             if(SNODEP_obs_obj(n)%mesh.eq.8) then 
                if(ihemi.eq.1) then 
                   xmesh = xmeshl1
                   orient = 100.0
                else
                   xmesh = -xmeshl1
                   orient = 280.0
                endif
                xj1 = float(1)-ypnmcaf1
                xi1 = float(1)-xpnmcaf1
                
                call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
             elseif(SNODEP_obs_obj(n)%mesh.eq.16) then 
                if(ihemi.eq.1) then 
                   xmesh = xmeshl2
                   orient = 100.0
                else
                   xmesh = -xmeshl2
                   orient = 280.0
                endif
                xj1 = float(1)-ypnmcaf2
                xi1 = float(1)-xpnmcaf2
                
                call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
             endif

             gridDesci = 0 
             gridDesci(1) = 5
             gridDesci(2) = SNODEP_obs_obj(n)%pmax
             gridDesci(3) = SNODEP_obs_obj(n)%pmax
             gridDesci(4) = alat1
             gridDesci(5) = alon1
             gridDesci(6) = 8
             gridDesci(7) = orient
             gridDesci(8) = xmesh
             gridDesci(9) = xmesh
             gridDesci(10) = 0.0       
             if(ihemi .eq.2) then 
                gridDesci(20) = 128
                gridDesci(11) = 128
             endif
             gridDesci(13) = 1  !global grid
             gridDesci(20) = 0 
             
             if(SNODEP_obs_obj(n)%gridspan.eq.1) then 
                allocate(SNODEP_obs_obj(n)%rlat1_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%rlon1_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n111_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n121_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n211_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n221_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w111_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w121_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w211_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w221_nh(SNODEP_obs_obj(n)%mo1))
#if 0
                call bilinear_interp_input_withgrid(gridDesci&
                     ,SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo1,SNODEP_obs_obj(n)%rlat1_nh,&
                     SNODEP_obs_obj(n)%rlon1_nh,SNODEP_obs_obj(n)%n111_nh,&
                     SNODEP_obs_obj(n)%n121_nh,SNODEP_obs_obj(n)%n211_nh,&
                     SNODEP_obs_obj(n)%n221_nh,SNODEP_obs_obj(n)%w111_nh,&
                     SNODEP_obs_obj(n)%w121_nh,SNODEP_obs_obj(n)%w211_nh,&
                     SNODEP_obs_obj(n)%w221_nh)
#endif

                call neighbor_interp_input_withgrid(gridDesci&
                     ,SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo1,SNODEP_obs_obj(n)%rlat1_nh,&
                     SNODEP_obs_obj(n)%rlon1_nh,SNODEP_obs_obj(n)%n111_nh)


             elseif(SNODEP_obs_obj(n)%gridspan.eq.2) then 
                allocate(SNODEP_obs_obj(n)%rlat1_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%rlon1_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n111_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n121_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n211_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n221_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w111_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w121_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w211_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w221_sh(SNODEP_obs_obj(n)%mo2))
#if 0 
                call bilinear_interp_input_withgrid(gridDesci,&
                     SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo2,SNODEP_obs_obj(n)%rlat1_sh,&
                     SNODEP_obs_obj(n)%rlon1_sh,&
                     SNODEP_obs_obj(n)%n111_sh,SNODEP_obs_obj(n)%n121_sh,&
                     SNODEP_obs_obj(n)%n211_sh,SNODEP_obs_obj(n)%n221_sh,&
                     SNODEP_obs_obj(n)%w111_sh,SNODEP_obs_obj(n)%w121_sh,&
                     SNODEP_obs_obj(n)%w211_sh,SNODEP_obs_obj(n)%w221_sh)
#endif
                call neighbor_interp_input_withgrid(gridDesci,&
                     SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo2,SNODEP_obs_obj(n)%rlat1_sh,&
                     SNODEP_obs_obj(n)%rlon1_sh,&
                     SNODEP_obs_obj(n)%n111_sh)
             else
                if(ihemi.eq.1) then 
                   allocate(SNODEP_obs_obj(n)%rlat1_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%rlon1_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n111_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n121_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n211_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n221_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w111_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w121_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w211_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w221_nh(SNODEP_obs_obj(n)%mo1))

#if 0 
                   call bilinear_interp_input_withgrid(gridDesci,&
                        SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                        SNODEP_obs_obj(n)%mo1,SNODEP_obs_obj(n)%rlat1_nh,&
                        SNODEP_obs_obj(n)%rlon1_nh,SNODEP_obs_obj(n)%n111_nh,&
                        SNODEP_obs_obj(n)%n121_nh,SNODEP_obs_obj(n)%n211_nh,&
                        SNODEP_obs_obj(n)%n221_nh,SNODEP_obs_obj(n)%w111_nh,&
                        SNODEP_obs_obj(n)%w121_nh,SNODEP_obs_obj(n)%w211_nh,&
                        SNODEP_obs_obj(n)%w221_nh)
#endif
                   call neighbor_interp_input_withgrid(gridDesci,&
                        SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                        SNODEP_obs_obj(n)%mo1,SNODEP_obs_obj(n)%rlat1_nh,&
                        SNODEP_obs_obj(n)%rlon1_nh,SNODEP_obs_obj(n)%n111_nh)
                elseif(ihemi.eq.2) then 
                   allocate(SNODEP_obs_obj(n)%rlat1_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%rlon1_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n111_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n121_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n211_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n221_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w111_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w121_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w211_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w221_sh(SNODEP_obs_obj(n)%mo2))

                   call neighbor_interp_input_withgrid(gridDesci,&
                        SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                        SNODEP_obs_obj(n)%mo2,SNODEP_obs_obj(n)%rlat1_sh,&
                        SNODEP_obs_obj(n)%rlon1_sh,SNODEP_obs_obj(n)%n111_sh)
#if 0 
                   call bilinear_interp_input_withgrid(gridDesci,&
                        SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                        SNODEP_obs_obj(n)%mo2,SNODEP_obs_obj(n)%rlat1_sh,&
                        SNODEP_obs_obj(n)%rlon1_sh,SNODEP_obs_obj(n)%n111_sh,&
                        SNODEP_obs_obj(n)%n121_sh,SNODEP_obs_obj(n)%n211_sh,&
                        SNODEP_obs_obj(n)%n221_sh,SNODEP_obs_obj(n)%w111_sh,&
                        SNODEP_obs_obj(n)%w121_sh,SNODEP_obs_obj(n)%w211_sh,&
                        SNODEP_obs_obj(n)%w221_sh)
#endif
                endif
             endif
          enddo
       elseif(LIS_rc%obs_gridDesc(k,1).ne.0) then 
          if(LIS_rc%obs_gridDesc(k,4).ge.0.and.LIS_rc%obs_gridDesc(k,10).ge.0) then 
             SNODEP_obs_obj(n)%gridspan = 1 ! NH only
             SNODEP_obs_obj(n)%shemi = 1
             SNODEP_obs_obj(n)%nhemi = 1
          elseif(LIS_rc%obs_gridDesc(k,4).le.0.and.LIS_rc%obs_gridDesc(k,10).ge.0) then 
             SNODEP_obs_obj(n)%gridspan = 2 ! SH only 
             SNODEP_obs_obj(n)%shemi = 2
             SNODEP_obs_obj(n)%nhemi = 2
          elseif(LIS_rc%obs_gridDesc(k,10).eq.-100.and.LIS_rc%obs_gridDesc(k,20).eq.-100) then 
!-----------------------------------------------------------------------------
! Global grid in polar stereographic projection. No interpolation 
! will be done.
!-----------------------------------------------------------------------------
             SNODEP_obs_obj(n)%gridspan = 3
             SNODEP_obs_obj(n)%shemi = 1
             SNODEP_obs_obj(n)%nhemi = 2
          else
             write(LIS_logunit,*) '[ERR] Currently spanning across hemispheres is'
             write(LIS_logunit,*) '[ERR] not supported for SNODEP data'
             call LIS_endrun()
          endif
          SNODEP_obs_obj(n)%hemi_nc = LIS_rc%obs_gridDesc(k,2)
          if(SNODEP_obs_obj(n)%gridspan.eq.1) then 
             SNODEP_obs_obj(n)%hemi_nr(1) = LIS_rc%obs_gridDesc(k,3)
             SNODEP_obs_obj(n)%hemi_nr(2) = 0 
             SNODEP_obs_obj(n)%mo1 = SNODEP_obs_obj(n)%hemi_nc(1)*&
                  SNODEP_obs_obj(n)%hemi_nr(1)
             SNODEP_obs_obj(n)%mo2 = 0 
          elseif(SNODEP_obs_obj(n)%gridspan.eq.2) then 
             SNODEP_obs_obj(n)%hemi_nr(1) = 0 
             SNODEP_obs_obj(n)%hemi_nr(2) = LIS_rc%obs_gridDesc(k,3)
             SNODEP_obs_obj(n)%mo1 = 0 
             SNODEP_obs_obj(n)%mo2 = SNODEP_obs_obj(n)%hemi_nc(2)*&
                  SNODEP_obs_obj(n)%hemi_nr(2)
          else
             SNODEP_obs_obj(n)%hemi_nr(1) = LIS_rc%obs_gridDesc(k,3)
             SNODEP_obs_obj(n)%hemi_nr(2) = LIS_rc%obs_gridDesc(k,3)
             SNODEP_obs_obj(n)%mo1 = SNODEP_obs_obj(n)%hemi_nc(1)*&
                  SNODEP_obs_obj(n)%hemi_nr(1)
             SNODEP_obs_obj(n)%mo2 = SNODEP_obs_obj(n)%hemi_nc(2)*&
                  SNODEP_obs_obj(n)%hemi_nr(2)             
          endif
          SNODEP_obs_obj(n)%gridDesco = 0 
          gridDesci = 0 
          do ihemi = SNODEP_obs_obj(n)%shemi,SNODEP_obs_obj(n)%nhemi
             SNODEP_obs_obj(n)%gridDesco(ihemi,1) = LIS_rc%obs_gridDesc(k,1)
             SNODEP_obs_obj(n)%gridDesco(ihemi,2) = SNODEP_obs_obj(n)%hemi_nc(ihemi)
             SNODEP_obs_obj(n)%gridDesco(ihemi,3) = SNODEP_obs_obj(n)%hemi_nr(ihemi)
             SNODEP_obs_obj(n)%gridDesco(ihemi,5) = LIS_rc%obs_gridDesc(k,5)
             SNODEP_obs_obj(n)%gridDesco(ihemi,8) = LIS_rc%obs_gridDesc(k,8)
             SNODEP_obs_obj(n)%gridDesco(ihemi,6) = LIS_rc%obs_gridDesc(k,6)
             SNODEP_obs_obj(n)%gridDesco(ihemi,9) = LIS_rc%obs_gridDesc(k,9)
             SNODEP_obs_obj(n)%gridDesco(ihemi,10) = LIS_rc%obs_gridDesc(k,10)
             SNODEP_obs_obj(n)%gridDesco(ihemi,11) = LIS_rc%obs_gridDesc(k,11)
             SNODEP_obs_obj(n)%gridDesco(ihemi,20) = 255
             if(SNODEP_obs_obj(n)%gridspan.eq.1.or.&
                  SNODEP_obs_obj(n)%gridspan.eq.2) then 
                SNODEP_obs_obj(n)%gridDesco(ihemi,4) = LIS_rc%obs_gridDesc(k,4)
                SNODEP_obs_obj(n)%gridDesco(ihemi,7) = LIS_rc%obs_gridDesc(k,7)
             endif
             if(SNODEP_obs_obj(n)%mesh.eq.8) then 
                if(ihemi.eq.1) then 
                   xmesh = xmeshl1
                   orient = 100.0
                else
                   xmesh = -xmeshl1
                   orient = 280.0
                endif
                xj1 = float(1)-ypnmcaf1
                xi1 = float(1)-xpnmcaf1
                
                call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
             elseif(SNODEP_obs_obj(n)%mesh.eq.16) then 
                if(ihemi.eq.1) then 
                   xmesh = xmeshl2
                   orient = 100.0
                else
                   xmesh = -xmeshl2
                   orient = 280.0
                endif
                xj1 = float(1)-ypnmcaf2
                xi1 = float(1)-xpnmcaf2
                
                call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
             endif

             gridDesci = 0 
             gridDesci(1) = 5
             gridDesci(2) = SNODEP_obs_obj(n)%pmax
             gridDesci(3) = SNODEP_obs_obj(n)%pmax
             gridDesci(4) = alat1
             gridDesci(5) = alon1
             gridDesci(6) = 8
             gridDesci(7) = orient
             gridDesci(8) = xmesh
             gridDesci(9) = xmesh
             gridDesci(10) = 0.0
             if(ihemi .eq.2) then 
                gridDesci(20) = 128
                gridDesci(11) = 128
             endif
             
             gridDesci(20) = 0 
             gridDesci(13) = 1  !global grid
             
             if(SNODEP_obs_obj(n)%gridspan.eq.1) then 
                allocate(SNODEP_obs_obj(n)%rlat1_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%rlon1_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n111_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n121_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n211_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%n221_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w111_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w121_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w211_nh(SNODEP_obs_obj(n)%mo1))
                allocate(SNODEP_obs_obj(n)%w221_nh(SNODEP_obs_obj(n)%mo1))
                call neighbor_interp_input_withgrid(gridDesci,&
                     SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo1,SNODEP_obs_obj(n)%rlat1_nh,&
                     SNODEP_obs_obj(n)%rlon1_nh,SNODEP_obs_obj(n)%n111_nh)
#if 0 
                call bilinear_interp_input_withgrid(gridDesci,&
                     SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo1,SNODEP_obs_obj(n)%rlat1_nh,&
                     SNODEP_obs_obj(n)%rlon1_nh,SNODEP_obs_obj(n)%n111_nh,&
                     SNODEP_obs_obj(n)%n121_nh,SNODEP_obs_obj(n)%n211_nh,&
                     SNODEP_obs_obj(n)%n221_nh,SNODEP_obs_obj(n)%w111_nh,&
                     SNODEP_obs_obj(n)%w121_nh,SNODEP_obs_obj(n)%w211_nh,&
                     SNODEP_obs_obj(n)%w221_nh)
#endif
             elseif(SNODEP_obs_obj(n)%gridspan.eq.2) then 
                
                allocate(SNODEP_obs_obj(n)%rlat1_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%rlon1_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n111_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n121_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n211_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%n221_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w111_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w121_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w211_sh(SNODEP_obs_obj(n)%mo2))
                allocate(SNODEP_obs_obj(n)%w221_sh(SNODEP_obs_obj(n)%mo2))

                call neighbor_interp_input_withgrid(gridDesci,&
                     SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo2,SNODEP_obs_obj(n)%rlat1_sh,&
                     SNODEP_obs_obj(n)%rlon1_sh,&
                     SNODEP_obs_obj(n)%n111_sh)
#if 0
                call bilinear_interp_input_withgrid(gridDesci,&
                     SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                     SNODEP_obs_obj(n)%mo2,SNODEP_obs_obj(n)%rlat1_sh,&
                     SNODEP_obs_obj(n)%rlon1_sh,&
                     SNODEP_obs_obj(n)%n111_sh,SNODEP_obs_obj(n)%n121_sh,&
                     SNODEP_obs_obj(n)%n211_sh,SNODEP_obs_obj(n)%n221_sh,&
                     SNODEP_obs_obj(n)%w111_sh,SNODEP_obs_obj(n)%w121_sh,&
                     SNODEP_obs_obj(n)%w211_sh,SNODEP_obs_obj(n)%w221_sh)
#endif
             elseif(SNODEP_obs_obj(n)%gridspan.eq.3) then 
!-------------------------------------------------------------------------
!   No interpolation is being done. So no weights will be computed. 
!-------------------------------------------------------------------------
                if(ihemi.eq.1) then 
                   allocate(SNODEP_obs_obj(n)%rlat1_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%rlon1_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n111_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n121_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n211_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%n221_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w111_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w121_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w211_nh(SNODEP_obs_obj(n)%mo1))
                   allocate(SNODEP_obs_obj(n)%w221_nh(SNODEP_obs_obj(n)%mo1))
                elseif(ihemi.eq.2) then 
                   allocate(SNODEP_obs_obj(n)%rlat1_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%rlon1_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n111_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n121_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n211_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%n221_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w111_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w121_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w211_sh(SNODEP_obs_obj(n)%mo2))
                   allocate(SNODEP_obs_obj(n)%w221_sh(SNODEP_obs_obj(n)%mo2))
                endif
             endif
          enddo
       endif

       call LIS_registerAlarm("SNODEP read alarm", 3600.0, 3600.0)

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

    enddo
    write(LIS_logunit,*) '[INFO] Created the States to hold the observations data'
    
  end subroutine SNODEPobs_setup

  subroutine SNODEPobs_finalize()

  implicit none

  end subroutine SNODEPobs_finalize
  
end module SNODEPobs_Mod
