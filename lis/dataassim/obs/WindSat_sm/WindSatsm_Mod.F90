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
! !MODULE: WindSatsm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle WindSat soil moisture retrievals. This plugin processes
!   the retrievals based on 10, 18.7 and 37GHz channels. For more
!   details, please see: 
! 
!   Li et al, ``WindSat global soil moisture retrieval and validation'',
!   IEEE Transactions on Geoscience and Remote Sensing, 2009
! 
! !REVISION HISTORY: 
!  22 Dec 09    Sujay Kumar;   Initial Specification
! 
module WindSatsm_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: WindSatsm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: WindSatsm_struc

  type, public ::  WindSatsm_dec

     logical             :: startMode
     integer             :: mi
     integer, allocatable    :: n11(:)
     
     real,    allocatable    :: model_xrange(:,:)
     real,    allocatable    :: obs_xrange(:,:)
     real,    allocatable    :: model_cdf(:,:)
     real,    allocatable    :: obs_cdf(:,:)
     real,    allocatable    :: model_mu(:)
     real,    allocatable    :: obs_mu(:)
     real,    allocatable    :: model_sigma(:)
     real,    allocatable    :: obs_sigma(:)
     real,    allocatable    :: smobs(:)
     real,    allocatable    :: tmobs(:)

     integer             :: nbins
     integer             :: scal
  end type WindSatsm_dec

  type(WindSatsm_dec), allocatable :: WindSatsm_struc(:)

contains
!BOP
! 
! !ROUTINE: WindSatsm_setup
! \label{WindSatsm_setup}
! 
! !INTERFACE: 
  subroutine WindSatsm_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_dataAssimMod
    use LIS_perturbMod
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
!   creation of data strctures required for WindSat soil
!   moisture data. 
!
!   Please note that currently only neighbor interpolation is supported 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n,i,t
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  synsmobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    real                   :: gridDesci(LIS_rc%nnest,50)
    real, allocatable          :: xrange(:), cdf(:)
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
    character(len=LIS_CONST_PATH_LEN) :: obscdffile(LIS_rc%nnest)
    real,  allocatable         :: ssdev(:)

    allocate(WindSatsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"WindSat soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,synsmobsdir,&
            rc=status)
       call LIS_verify(status,'WindSat soil moisture data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            synsmobsdir, rc=status)
       call LIS_verify(status)
    enddo
    call ESMF_ConfigFindLabel(LIS_config,"WindSat scale observations:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,WindSatsm_struc(n)%scal, &
            rc=status)
       call LIS_verify(status, "WindSat scale observations: not defined")
    enddo
   
    call ESMF_ConfigFindLabel(LIS_config,"WindSat model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(WindSatsm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'WindSat model CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WindSat observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(WindSatsm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
          call LIS_verify(status, 'WindSat observation CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WindSat number of bins in the CDF:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(WindSatsm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,WindSatsm_struc(n)%nbins,rc=status)
          call LIS_verify(status, 'WindSat number of bins in the CDF: not defined')
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
            LIS_rc%ngrid(n),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*)'read WindSat soil moisure data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       allocate(WindSatsm_struc(n)%smobs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(WindSatsm_struc(n)%tmobs(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       WindSatsm_struc(n)%smobs = LIS_rc%udef
       WindSatsm_struc(n)%tmobs = LIS_rc%udef

       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_vecGrid(n),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)
       
!Perturbations State
       write(LIS_logunit,*) 'Opening attributes for observations ',&
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
          write(LIS_logunit,*) vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)  
              
       allocate(ssdev(LIS_rc%ngrid(n)))

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

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_ensOnGrid(n),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1),&
               rc=status)
          call LIS_verify(status)

          if(LIS_rc%ngrid(n).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
             call LIS_verify(status)
          endif

          call ESMF_AttributeSet(pertField(n),"Std Normal Max", obs_pert%stdmax(1),&
               rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean", obs_pert%zeromean(1),&
               rc=status)
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

!-------------------------------------------------------------
! setting up interpolation weights for neighbor search
!-------------------------------------------------------------
    gridDesci = 0 

    do n=1,LIS_rc%nnest

       gridDesci(n,1) = 9
       gridDesci(n,2) = 1383
       gridDesci(n,3) = 586
       gridDesci(n,4) = -90.0
       gridDesci(n,5) = -179.6096
       gridDesci(n,7) = 83.33788
       gridDesci(n,8) = 180.1301
       gridDesci(n,9)  = 1  !Ml

       WindSatsm_struc(n)%mi = 1383*586
       allocate(WindSatsm_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       call neighbor_interp_input(n,gridDesci(n,:),&
            WindSatsm_struc(n)%n11)
    enddo

!--------------------------------------------------------------------------------
! Read the model and observation CDF data
!--------------------------------------------------------------------------------
    do n=1,LIS_rc%nnest
       if(WindSatsm_struc(n)%scal.eq.1) then 

          allocate(ssdev(LIS_rc%ngrid(n)))
          ssdev =obs_pert%ssdev(1)
          
          allocate(WindSatsm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), WindSatsm_struc(n)%nbins))
          allocate(WindSatsm_struc(n)%obs_xrange(&
               LIS_rc%ngrid(n), WindSatsm_struc(n)%nbins))
          allocate(WindSatsm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), WindSatsm_struc(n)%nbins))
          allocate(WindSatsm_struc(n)%obs_cdf(&
               LIS_rc%ngrid(n), WindSatsm_struc(n)%nbins))
          allocate(WindSatsm_struc(n)%model_mu(LIS_rc%ngrid(n)))
          allocate(WindSatsm_struc(n)%model_sigma(LIS_rc%ngrid(n)))
          allocate(WindSatsm_struc(n)%obs_mu(LIS_rc%ngrid(n)))
          allocate(WindSatsm_struc(n)%obs_sigma(LIS_rc%ngrid(n)))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
          call LIS_readMeanSigmaData(n,&
               modelcdffile(n), &
               "SoilMoist",&
               WindSatsm_struc(n)%model_mu,&
               WindSatsm_struc(n)%model_sigma)

          call LIS_readMeanSigmaData(n,&
               obscdffile(n), &
               "SoilMoist",&
               WindSatsm_struc(n)%obs_mu,&
               WindSatsm_struc(n)%obs_sigma)

          call LIS_readCDFdata(n,&
               WindSatsm_struc(n)%nbins,&
               modelcdffile(n), &
               "SoilMoist",&
               WindSatsm_struc(n)%model_xrange,&
               WindSatsm_struc(n)%model_cdf)

          call LIS_readCDFdata(n,&
               WindSatsm_struc(n)%nbins,&
               obscdffile(n), &
               "SoilMoist",&
               WindSatsm_struc(n)%obs_xrange,&
               WindSatsm_struc(n)%obs_cdf)

          do t=1,LIS_rc%ngrid(n)
             if(WindSatsm_struc(n)%obs_sigma(t).gt.0) then 
!                ssdev(t) = ssdev(t)*WindSatsm_struc(n)%model_sigma(t)/&
!                     WindSatsm_struc(n)%obs_sigma(t)
!                print*, ssdev(t), WindSatsm_struc(n)%model_sigma(t), &
!                     WindSatsm_struc(n)%obs_sigma(t)
!                if(ssdev(t).gt.0.20) ssdev(t) = 0.20 !setting as the max obs err
             endif
          enddo

          if(LIS_rc%ngrid(n).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
             call LIS_verify(status)
          endif

          deallocate(ssdev)
       endif
    enddo

    do n=1,LIS_rc%nnest       

       WindSatsm_struc(n)%startMode = .true.

       call LIS_registerAlarm("WindSat sm read alarm", &
            86400.0, 86400.0)
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)
       

    enddo

    write(LIS_logunit,*) 'Created ESMF States to hold the WindSat observations data'

  end subroutine WindSatsm_setup
  
end module WindSatsm_Mod

