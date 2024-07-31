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
! !MODULE: pildassmobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle Pildas soil moisture observations (such as the 
!   one generated from a previous LIS-Noah simulation)
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
!  9Sep2016  Mahdi Navari;  Modified for pildas

module pildassmobs_module
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN


  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: pildassmobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: pildas_sm_struc

  type, public :: Pildas_sm_dec
! MN
    logical                     :: startMode
    real                        :: version
     integer                    :: useSsdevScal
     integer                    :: nc
     integer                    :: nr

     real,     allocatable      :: smobs(:,:)
     real,     allocatable      :: smtime(:,:)

     real                       :: ssdev_inp
     integer                    :: pil_nc, pil_nr
     type(proj_info)            :: ecvproj
     integer, allocatable       :: n11(:)
     real,    allocatable       :: rlat(:)
     real,    allocatable       :: rlon(:)
! MN end
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
     integer                :: ngrid

  end type Pildas_sm_dec
  type(Pildas_sm_dec), allocatable :: pildas_sm_struc(:)

contains

!BOP
! 
! !ROUTINE: pildassmobs_setup
! \label{pildassmobs_setup}
! 
! !INTERFACE: 
  subroutine pildassmobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
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
!   creation of data strctures required for soil moisture
!    assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    real, parameter        ::  minssdev = 0.035
    integer                ::  n,i,t,kk,jj
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  pildassmobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
   integer                 ::  ngrid
    real, pointer          ::  obs_temp(:,:)
    real, allocatable          :: xrange(:), cdf(:)
    real                   :: gridDesci(50)
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
    character(len=LIS_CONST_PATH_LEN) :: obscdffile(LIS_rc%nnest)
    real, allocatable          ::  ssdev(:)

    allocate(pildas_sm_struc(LIS_rc%nnest))
    ! MN: Creates a description of the data – the typekind, the rank, and the dimensionality
    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)  !MN: I think we don't need this 
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    ! MN: Get the path from LIS_config file
    call ESMF_ConfigFindLabel(LIS_config,"PILDAS soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,pildassmobsdir,&
            rc=status)
       call LIS_verify(status)
     ! MN: save the path in the state var.
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            pildassmobsdir, rc=status)
       call LIS_verify(status)
    enddo

!    call ESMF_ConfigFindLabel(LIS_config,"PILDAS soil moisture data version:",&
!         rc=status)
!    do n=1,LIS_rc%nnest
!       call ESMF_ConfigGetAttribute(LIS_config,pildas_sm_struc(n)%version,&
!            rc=status)
!       call LIS_verify(status, 'PILDAS soil moisture data version: is missing')

!    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PILDAS use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,pildas_sm_struc(n)%useSsdevScal, &
            rc=status)
       call LIS_verify(status, "PILDAS use scaled standard deviation model: not defined")
    enddo


    ! MN: Get model CDF filename from LIS_config
    call ESMF_ConfigFindLabel(LIS_config,"PILDAS soil moisture model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'PILDAS soil moisture model CDF file: not defined')
       endif
    enddo
    ! MN: Get pildas CDF filename from LIS_config
    call ESMF_ConfigFindLabel(LIS_config,"PILDAS soil moisture observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
          call LIS_verify(status, 'PILDAS soil moisture observation CDF file: not defined')
       endif
    enddo
     ! MN: Get number of bins in the CDF filename from LIS_config    
    call ESMF_ConfigFindLabel(LIS_config, "PILDAS soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,pildas_sm_struc(n)%nbins, rc=status)
          call LIS_verify(status, "PILDAS soil moisture number of bins in the CDF: not defined")
       endif
    enddo
    ! MN: Initialize the state variable (OBS_State) attribute 
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

    write(LIS_logunit,*)'[INFO] read PILDAS soil moisture data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%obs_ngrid(k). 
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
          ! MN assign value to obs_pert from LIS
          call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
               obs_pert)

! Set obs err to be uniform (will be rescaled later for each grid point). 
          ssdev = obs_pert%ssdev(1)
          pildas_sm_struc(n)%ssdev_inp = obs_pert%ssdev(1)

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 
          ! MN: Set the ESMF attributes
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
    write(LIS_logunit,*) '[INFO] Created the States to hold the observations data'

    do n=1,LIS_rc%nnest
       pildas_sm_struc(n)%pil_nc = 64    !????????  64 or 1440 
       pildas_sm_struc(n)%pil_nr = 44    !????????  44 0r 720 
       allocate(pildas_sm_struc(n)%smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(pildas_sm_struc(n)%smtime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       pildas_sm_struc(n)%smtime = -1

    enddo

    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          
          call LIS_getCDFattributes(k, modelcdffile(n),&
               pildas_sm_struc(n)%ntimes, pildas_sm_struc(n)%ngrid)

          allocate(ssdev(LIS_rc%obs_ngrid(k)))
          ssdev = obs_pert%ssdev(1)

          allocate(pildas_sm_struc(n)%model_mu(LIS_rc%obs_ngrid(k), &
               pildas_sm_struc(n)%ntimes))
          allocate(pildas_sm_struc(n)%model_sigma(&
               LIS_rc%obs_ngrid(k),&
               pildas_sm_struc(n)%ntimes))
          allocate(pildas_sm_struc(n)%obs_mu(&
               LIS_rc%obs_ngrid(k),pildas_sm_struc(n)%ntimes))
          allocate(pildas_sm_struc(n)%obs_sigma(&
               LIS_rc%obs_ngrid(k),pildas_sm_struc(n)%ntimes))
          allocate(pildas_sm_struc(n)%model_xrange(&
               LIS_rc%obs_ngrid(k), pildas_sm_struc(n)%ntimes, &
               pildas_sm_struc(n)%nbins))
          allocate(pildas_sm_struc(n)%obs_xrange(&
               LIS_rc%obs_ngrid(k), pildas_sm_struc(n)%ntimes, &
               pildas_sm_struc(n)%nbins))
          allocate(pildas_sm_struc(n)%model_cdf(&
               LIS_rc%obs_ngrid(k), pildas_sm_struc(n)%ntimes, &
               pildas_sm_struc(n)%nbins))
          allocate(pildas_sm_struc(n)%obs_cdf(&
               LIS_rc%obs_ngrid(k), pildas_sm_struc(n)%ntimes, &
               pildas_sm_struc(n)%nbins))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
          call LIS_readMeanSigmaData(n,k, &
               pildas_sm_struc(n)%ntimes, &
               LIS_rc%obs_ngrid(k), &
               modelcdffile(n), &
               "SoilMoist",&
               pildas_sm_struc(n)%model_mu,&
               pildas_sm_struc(n)%model_sigma)

          call LIS_readMeanSigmaData(n,k,&
               pildas_sm_struc(n)%ntimes, &
               LIS_rc%obs_ngrid(k), &
               obscdffile(n), &
               "SoilMoist",&
               pildas_sm_struc(n)%obs_mu,&
               pildas_sm_struc(n)%obs_sigma)

          call LIS_readCDFdata(n,k,&
               pildas_sm_struc(n)%nbins, &
               pildas_sm_struc(n)%ntimes, &
               LIS_rc%obs_ngrid(k), &
               modelcdffile(n), &
               "SoilMoist",&
               pildas_sm_struc(n)%model_xrange,&
               pildas_sm_struc(n)%model_cdf)

          call LIS_readCDFdata(n,k,&
               pildas_sm_struc(n)%nbins,&
               pildas_sm_struc(n)%ntimes, &
               LIS_rc%obs_ngrid(k), &
               obscdffile(n), &
               "SoilMoist",&
               pildas_sm_struc(n)%obs_xrange,&
               pildas_sm_struc(n)%obs_cdf)

          if(pildas_sm_struc(n)%useSsdevScal.eq.1) then 
             if(pildas_sm_struc(n)%ntimes.eq.1) then 
                jj = 1
             else
                jj = LIS_rc%mo
             endif
             do t=1,LIS_rc%obs_ngrid(k)
                if(pildas_sm_struc(n)%obs_sigma(t,jj).ne.LIS_rc%udef) then 
                   ssdev(t) = ssdev(t)*pildas_sm_struc(n)%model_sigma(t,jj)/&
                        pildas_sm_struc(n)%obs_sigma(t,jj)
                   if(ssdev(t).lt.minssdev) then 
                      ssdev(t) = minssdev
                   endif
                endif
             enddo
          endif

          if(LIS_rc%obs_ngrid(k).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
             call LIS_verify(status)
          endif

          deallocate(ssdev)

       endif
    enddo

   do n=1,LIS_rc%nnest

       pildas_sm_struc(n)%pil_nc = 64   
       pildas_sm_struc(n)%pil_nr = 44  

       call map_set(PROJ_LATLON, -89.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            pildas_sm_struc(n)%pil_nc,pildas_sm_struc(n)%pil_nr,&
            pildas_sm_struc(n)%ecvproj)
       
       gridDesci = 0 ! no change
       gridDesci(1) = 0 
       gridDesci(2) = 64
       gridDesci(3) = 44
       gridDesci(4) = 33.5625
       gridDesci(5) = -102.9375
       gridDesci(6) = 128 ! no change
       gridDesci(7) = 38.9375
       gridDesci(8) = -95.0625
       gridDesci(9) = 0.125
       gridDesci(10) = 0.125
       gridDesci(20) = 64 ! no change
       
       allocate(pildas_sm_struc(n)%n11(pildas_sm_struc(n)%pil_nc*&
            pildas_sm_struc(n)%pil_nr))
       allocate(pildas_sm_struc(n)%rlat(pildas_sm_struc(n)%pil_nc*&
            pildas_sm_struc(n)%pil_nr))
       allocate(pildas_sm_struc(n)%rlon(pildas_sm_struc(n)%pil_nc*&
            pildas_sm_struc(n)%pil_nr))
       ! compute interpolation weight 
       call neighbor_interp_input_withgrid(&
            gridDesci,&
            LIS_rc%obs_gridDesc(k,:),&
            pildas_sm_struc(n)%pil_nc*pildas_sm_struc(n)%pil_nr,&
            pildas_sm_struc(n)%rlat, &
            pildas_sm_struc(n)%rlon, &
            pildas_sm_struc(n)%n11)       

       !call LIS_registerAlarm("ESACCI read alarm",&
       !     86400.0, 86400.0)
!       pildas_sm_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)       

    enddo
  end subroutine pildassmobs_setup  
end module pildassmobs_module
