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
! !MODULE: syntheticsmobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle synthetic soil moisture observations (such as the 
!   one generated from a previous LIS-Noah simulation)
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
! 
module syntheticsmobs_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: syntheticsmobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: synthetic_sm_struc

  type, public :: synthetic_sm_dec
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

  end type synthetic_sm_dec
  type(synthetic_sm_dec), allocatable :: synthetic_sm_struc(:)
contains
!BOP
! 
! !ROUTINE: syntheticsmobs_setup
! \label{syntheticsmobs_setup}
! 
! !INTERFACE: 
  subroutine syntheticsmobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_dataAssimMod
    use LIS_perturbMod
    use LIS_DAobservationsMod

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
    integer                ::  n 
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  synsmobsdir
    character*100          ::  temp
    integer                ::  ftn,i,t
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real, allocatable          ::  ssdev(:)
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
    character(len=LIS_CONST_PATH_LEN) :: obscdffile(LIS_rc%nnest)

    allocate(synthetic_sm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"Synthetic soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synsmobsdir,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            synsmobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Synthetic soil moisture model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'Synthetic soil moisture model CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Synthetic soil moisture observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
          call LIS_verify(status, 'Synthetic soil moisture observation CDF file: not defined')
       endif
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "Synthetic soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       if(LIS_rc%dascaloption(k).ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config,synthetic_sm_struc(n)%nbins, rc=status)
          call LIS_verify(status, "Synthetic soil moisture number of bins in the CDF: not defined")
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

    write(LIS_logunit,*)'[INFO] read Synthetic soil moisture data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. 
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
!to obsarr

    enddo
    write(LIS_logunit,*) '[INFO] Created the States to hold the observations data'
    

    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       ssdev = obs_pert%ssdev(1)
       
       if(LIS_rc%dascaloption(k).ne."none") then 
          
          call LIS_getCDFattributes(k, modelcdffile(n),&
               synthetic_sm_struc(n)%ntimes, synthetic_sm_struc(n)%ngrid)

          synthetic_sm_struc(n)%ngrid = LIS_rc%obs_ngrid(k)

          allocate(synthetic_sm_struc(n)%model_mu(synthetic_sm_struc(n)%ngrid, &
               synthetic_sm_struc(n)%ntimes))
          allocate(synthetic_sm_struc(n)%model_sigma(&
               synthetic_sm_struc(n)%ngrid,&
               synthetic_sm_struc(n)%ntimes))
          allocate(synthetic_sm_struc(n)%obs_mu(&
               synthetic_sm_struc(n)%ngrid,synthetic_sm_struc(n)%ntimes))
          allocate(synthetic_sm_struc(n)%obs_sigma(&
               synthetic_sm_struc(n)%ngrid,synthetic_sm_struc(n)%ntimes))

          allocate(synthetic_sm_struc(n)%model_xrange(&
               synthetic_sm_struc(n)%ngrid, synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%nbins))
          allocate(synthetic_sm_struc(n)%obs_xrange(&
               synthetic_sm_struc(n)%ngrid, synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%nbins))
          allocate(synthetic_sm_struc(n)%model_cdf(&
               synthetic_sm_struc(n)%ngrid, synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%nbins))
          allocate(synthetic_sm_struc(n)%obs_cdf(&
               synthetic_sm_struc(n)%ngrid, synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%nbins))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
          call LIS_readMeanSigmaData(n,k, &
               synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               synthetic_sm_struc(n)%model_mu,&
               synthetic_sm_struc(n)%model_sigma)

          call LIS_readMeanSigmaData(n,k,&
               synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%ngrid, &
               obscdffile(n), &
               "SoilMoist",&
               synthetic_sm_struc(n)%obs_mu,&
               synthetic_sm_struc(n)%obs_sigma)

          call LIS_readCDFdata(n,k,&
               synthetic_sm_struc(n)%nbins, &
               synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               synthetic_sm_struc(n)%model_xrange,&
               synthetic_sm_struc(n)%model_cdf)

          call LIS_readCDFdata(n,k,&
               synthetic_sm_struc(n)%nbins,&
               synthetic_sm_struc(n)%ntimes, &
               synthetic_sm_struc(n)%ngrid, &
               obscdffile(n), &
               "SoilMoist",&
               synthetic_sm_struc(n)%obs_xrange,&
               synthetic_sm_struc(n)%obs_cdf)
!          do t=1,synthetic_sm_struc(n)%ngrid
!             if(synthetic_sm_struc(n)%obs_sigma(t).ne.LIS_rc%udef) then 
!                ssdev(t) = ssdev(t)*synthetic_sm_struc(n)%model_sigma(t)/&
!                     synthetic_sm_struc(n)%obs_sigma(t)
!             endif
!          enddo
       else
          synthetic_sm_struc(n)%ngrid = LIS_rc%obs_ngrid(k)

       endif

       if(synthetic_sm_struc(n)%ngrid.gt.0) then 
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=synthetic_sm_struc(n)%ngrid,rc=status)
          call LIS_verify(status, 'Error in AttributeSet: Standard Deviation')
       endif
       
       deallocate(ssdev)

    enddo
    
    do n=1,LIS_rc%nnest
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)
    enddo

  end subroutine syntheticsmobs_setup
  
end module syntheticsmobs_module
