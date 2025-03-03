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
! !MODULE: LPRM_AMSREsm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle Land Parameter Retrieval Model (LPRM) AMSR-E soil moisture
!   retrievals
! 
!   Ref: Owe et al. 2008; Multi-sensor historical climatology of satellite-
!   derived global land surface moisture. 
!   
! !REVISION HISTORY: 
!  17 Jun 10    Sujay Kumar; Updated for use with LPRM retrieval version 5. 
!  20 Sep 12    Sujay Kumar; Updated to the NETCDF version from GES-DISC. 
! 
module LPRM_AMSREsm_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LPRM_AMSREsm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LPRM_AMSREsm_struc
!EOP
  type, public:: LPRM_AMSREsm_dec
     
     logical                :: startMode
     integer                :: useSsdevScal
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: smobs(:,:)
     real,     allocatable      :: smtime(:,:)
     
     real                   :: ssdev_inp
     integer                :: lprmnc, lprmnr
     type(proj_info)        :: lprmproj
     integer                :: rawdata
     integer, allocatable       :: n11(:)
     real,    allocatable       :: rlat(:)
     real,    allocatable       :: rlon(:)

     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: obs_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: obs_cdf(:,:,:)
     real,    allocatable       :: model_mu(:,:)
     real,    allocatable       :: obs_mu(:,:)
     real,    allocatable       :: model_sigma(:,:)
     real,    allocatable       :: obs_sigma(:,:)

     integer                :: nbins
     integer                :: ntimes

  end type LPRM_AMSREsm_dec
  
  type(LPRM_AMSREsm_dec),allocatable :: LPRM_AMSREsm_struc(:)
  
contains

!BOP
! 
! !ROUTINE: LPRM_AMSREsm_setup
! \label{LPRM_AMSREsm_setup}
! 
! !INTERFACE: 
  subroutine LPRM_AMSREsm_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_dataAssimMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_logmod

    implicit none 

! !ARGUMENTS: 
    integer                ::  k
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for handling LPRM 
!   AMSR-E soil moisture data. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
!based on Liu et al. JHM 2011                   
!    real, parameter        ::  minssdev = 0.01
    real, parameter        ::  minssdev = 0.001
    real, parameter        ::  maxssdev = 0.11
    integer                ::  n,i,t,kk,c,r
    real, allocatable          ::  obserr(:,:)
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  amsresmobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real                   :: gridDesci(50)
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
    character(len=LIS_CONST_PATH_LEN) :: obscdffile(LIS_rc%nnest)
    real, allocatable          :: ssdev(:)
    integer                :: jj
    integer                :: ngrid

    allocate(LPRM_AMSREsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"AMSR-E(LPRM) soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,amsresmobsdir,&
            rc=status)
       call LIS_verify(status, 'AMSR-E(LPRM) soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            amsresmobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR-E(LPRM) soil moisture use raw data:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_struc(n)%rawdata,&
            rc=status)
       call LIS_verify(status, 'AMSR-E(LPRM) soil moisture use raw data: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR-E(LPRM) use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_struc(n)%useSsdevScal, &
            rc=status)
       call LIS_verify(status, "AMSR-E(LPRM) use scaled standard deviation model: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR-E(LPRM) model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'AMSR-E(LPRM) model CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR-E(LPRM) observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
          call LIS_verify(status, 'AMSR-E(LPRM) observation CDF file: not defined')
       endif
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "AMSR-E(LPRM) soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_struc(n)%nbins, rc=status)
          call LIS_verify(status, "AMSR-E(LPRM) soil moisture number of bins in the CDF: not defined")
       endif
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

    write(LIS_logunit,*)'[INFO] read AMSR-E(LPRM) soil moisture data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. amsr-e 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%obs_ngrid(k). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_obsvecGrid(n,k),&
            name="Observation"//vid(1)//vid(2),rc=status)
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
          LPRM_AMSREsm_struc(n)%ssdev_inp = obs_pert%ssdev(1)

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsensOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
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

       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)

    enddo
    write(LIS_logunit,*) &
         '[INFO] Created the States to hold the AMSR-E(LPRM) observations data'
    do n=1,LIS_rc%nnest
       LPRM_AMSREsm_struc(n)%nc = 1440
       LPRM_AMSREsm_struc(n)%nr = 720
       allocate(LPRM_AMSREsm_struc(n)%smobs(&
            LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(LPRM_AMSREsm_struc(n)%smtime(&
            LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))

    enddo
    
    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 

          call LIS_getCDFattributes(k,modelcdffile(n),&
               LPRM_AMSREsm_struc(n)%ntimes,ngrid)

          allocate(ssdev(LIS_rc%obs_ngrid(k)))
          ssdev = obs_pert%ssdev(1)

          allocate(LPRM_AMSREsm_struc(n)%model_xrange(&
               LIS_rc%obs_ngrid(k),  LPRM_AMSREsm_struc(n)%ntimes,&
               LPRM_AMSREsm_struc(n)%nbins))
          allocate(LPRM_AMSREsm_struc(n)%obs_xrange(&
               LIS_rc%obs_ngrid(k),  LPRM_AMSREsm_struc(n)%ntimes, &
               LPRM_AMSREsm_struc(n)%nbins))
          allocate(LPRM_AMSREsm_struc(n)%model_cdf(&
               LIS_rc%obs_ngrid(k),  LPRM_AMSREsm_struc(n)%ntimes, &
               LPRM_AMSREsm_struc(n)%nbins))
          allocate(LPRM_AMSREsm_struc(n)%obs_cdf(&
               LIS_rc%obs_ngrid(k),  LPRM_AMSREsm_struc(n)%ntimes, &
               LPRM_AMSREsm_struc(n)%nbins))
          allocate(LPRM_AMSREsm_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
               LPRM_AMSREsm_struc(n)%ntimes))
          allocate(LPRM_AMSREsm_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
               LPRM_AMSREsm_struc(n)%ntimes))
          allocate(LPRM_AMSREsm_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
               LPRM_AMSREsm_struc(n)%ntimes))
          allocate(LPRM_AMSREsm_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
               LPRM_AMSREsm_struc(n)%ntimes))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
          call LIS_readMeanSigmaData(n,k,&
               LPRM_AMSREsm_struc(n)%ntimes,& 
               LIS_rc%obs_ngrid(k), &
               modelcdffile(n), &
               "SoilMoist",&
               LPRM_AMSREsm_struc(n)%model_mu,&
               LPRM_AMSREsm_struc(n)%model_sigma)

          call LIS_readMeanSigmaData(n,k,&
               LPRM_AMSREsm_struc(n)%ntimes,& 
               LIS_rc%obs_ngrid(k), &
               obscdffile(n), &
               "SoilMoist",&
               LPRM_AMSREsm_struc(n)%obs_mu,&
               LPRM_AMSREsm_struc(n)%obs_sigma)

          call LIS_readCDFdata(n,k,&
               LPRM_AMSREsm_struc(n)%nbins,&
               LPRM_AMSREsm_struc(n)%ntimes,& 
               LIS_rc%obs_ngrid(k), &
               modelcdffile(n), &
               "SoilMoist",&
               LPRM_AMSREsm_struc(n)%model_xrange,&
               LPRM_AMSREsm_struc(n)%model_cdf)

          call LIS_readCDFdata(n,k,&
               LPRM_AMSREsm_struc(n)%nbins,&
               LPRM_AMSREsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               obscdffile(n), &
               "SoilMoist",&
               LPRM_AMSREsm_struc(n)%obs_xrange,&
               LPRM_AMSREsm_struc(n)%obs_cdf)

          if(LPRM_AMSREsm_struc(n)%useSsdevScal.eq.1) then 
             if(LPRM_AMSREsm_struc(n)%ntimes.eq.1) then 
                jj = 1
             else
                jj = LIS_rc%mo
             endif
             do t=1,LIS_rc%obs_ngrid(k)
                if(LPRM_AMSREsm_struc(n)%obs_sigma(t,jj).gt.0) then 
                   ssdev(t) = ssdev(t)*LPRM_AMSREsm_struc(n)%model_sigma(t,jj)/&
                        LPRM_AMSREsm_struc(n)%obs_sigma(t,jj)

                   if(ssdev(t).gt.maxssdev) ssdev(t) = maxssdev
                   if(ssdev(t).lt.minssdev) then 
                      ssdev(t) = minssdev
                   endif
                endif
             enddo
          endif

#if 0           
          allocate(obserr(LIS_rc%obs_gnc(n),LIS_rc%obs_gnr(n)))
          obserr = -9999.0

!          lobserr(:,:) = obserr(&
!               LIS_ews_halo_ind(n,LIS_localPet+1):&         
!               LIS_ewe_halo_ind(n,LIS_localPet+1), &
!               LIS_nss_halo_ind(n,LIS_localPet+1): &
!               LIS_nse_halo_ind(n,LIS_localPet+1))

          do r=1,LIS_rc%obs_lnr(k)
             do c=1,LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                   obserr(c,r)  =  ssdev(LIS_obs_domain(n,k)%gindex(c,r)) 
                   
                endif
             enddo
          enddo

!          do r=1,LIS_rc%obs_lnr(k)
!             do c=1,LIS_rc%obs_lnc(k)
!                if(LIS_domain(n)%gindex(c,r).ne.-1) then 
!                   lobserr(c,r) = ssdev(LIS_domain(n)%gindex(c,r)) 
!                   
!                endif
!             enddo
!          enddo
          print*, 'LPRM ',LIS_rc%obs_gnc(n),LIS_rc%obs_gnr(n)
          open(100,file='test.bin',form='unformatted')
          write(100) obserr
          close(100)
          stop
          deallocate(obserr)
#endif
         
          if(LIS_rc%obs_ngrid(k).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
             call LIS_verify(status)
          endif

          deallocate(ssdev)
       endif
    enddo

    do n=1,LIS_rc%nnest
       LPRM_AMSREsm_struc(n)%lprmnc = 1440
       LPRM_AMSREsm_struc(n)%lprmnr = 720

       call map_set(PROJ_LATLON, -89.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            LPRM_AMSREsm_struc(n)%lprmnc,LPRM_AMSREsm_struc(n)%lprmnr,&
            LPRM_AMSREsm_struc(n)%lprmproj)
       
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 1440
       gridDesci(3) = 720
       gridDesci(4) = -89.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
       
       allocate(LPRM_AMSREsm_struc(n)%n11(LPRM_AMSREsm_struc(n)%lprmnc*&
            LPRM_AMSREsm_struc(n)%lprmnr))
       allocate(LPRM_AMSREsm_struc(n)%rlat(LPRM_AMSREsm_struc(n)%lprmnc*&
            LPRM_AMSREsm_struc(n)%lprmnr))
       allocate(LPRM_AMSREsm_struc(n)%rlon(LPRM_AMSREsm_struc(n)%lprmnc*&
            LPRM_AMSREsm_struc(n)%lprmnr))
       
       call neighbor_interp_input_withgrid(gridDesci, &
            LIS_rc%obs_gridDesc(k,:),&
            LPRM_AMSREsm_struc(n)%lprmnc*&
            LPRM_AMSREsm_struc(n)%lprmnr,&
            LPRM_AMSREsm_struc(n)%rlat, &
            LPRM_AMSREsm_struc(n)%rlon, &
            LPRM_AMSREsm_struc(n)%n11)
       
       call LIS_registerAlarm("AMSR-E(LPRM) read alarm",&
            86400.0, 86400.0)
       LPRM_AMSREsm_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine LPRM_AMSREsm_setup
end module LPRM_AMSREsm_Mod
