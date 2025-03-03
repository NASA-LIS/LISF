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
! !MODULE: WindSatCsm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle WindSat soil moisture retrievals from the C-band. 
!   
! !REVISION HISTORY: 
!  22 Dec 09    Sujay Kumar;   Initial Specification
! 
module WindSatCsm_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: WindSatCsm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: WindSatCsm_struc

  type, public ::  WindSatCsm_dec
     integer             :: mi
     
     real,    allocatable    :: model_xrange(:,:)
     real,    allocatable    :: obs_xrange(:,:)
     real,    allocatable    :: model_cdf(:,:)
     real,    allocatable    :: obs_cdf(:,:)

     real,    allocatable    :: smobs(:)
     integer,    allocatable    :: tmobs(:)

     integer             :: nbins
     integer             :: scal
  end type WindSatCsm_dec

  type(WindSatCsm_dec), allocatable :: WindSatCsm_struc(:)

contains
!BOP
! 
! !ROUTINE: WindSatCsm_setup
! \label{WindSatCsm_setup}
! 
! !INTERFACE: 
  subroutine WindSatCsm_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid, LIS_ensOnGrid, &
         LIS_masterproc    
    use LIS_timeMgrMod, only : LIS_clock, LIS_calendar,LIS_get_julss
    use LIS_historyMod, only : LIS_readvar_gridded
    use LIS_perturbMod
    use LIS_logMod, only : LIS_logunit, LIS_verify, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber
   
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
    integer                ::  n,i
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  synsmobsdir
    character*100          ::  temp
    integer                ::  yr,mo,da,hr,mn,ss
    real                   ::  smvalue
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40           ::  vname(1)
    real                   ::  varmin(1)
    real                   ::  varmax(1)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    real                   :: gridDesci(LIS_rc%nnest,50)
    real, allocatable          :: xrange(:), cdf(:)
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
    character(len=LIS_CONST_PATH_LEN) :: obscdffile(LIS_rc%nnest)
    integer                :: count,dummy, ios, obstime

    allocate(WindSatCsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"WindSat C-band soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,synsmobsdir,&
            rc=status)
       call LIS_verify(status,'WindSat C-band soil moisture data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            synsmobsdir, rc=status)
       call LIS_verify(status)
    enddo
    call ESMF_ConfigFindLabel(LIS_config,"WindSat C-band scale observations:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,WindSatCsm_struc(n)%scal, &
            rc=status)
       call LIS_verify(status, "WindSat C-band scale observations: not defined")
    enddo
   
    call ESMF_ConfigFindLabel(LIS_config,"WindSat C-band model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(WindSatCsm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'WindSat C-band model CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WindSat C-band observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(WindSatCsm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
          call LIS_verify(status, 'WindSat C-band observation CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WindSat C-band number of bins in the CDF:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(WindSatCsm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,WindSatCsm_struc(n)%nbins,rc=status)
          call LIS_verify(status, 'WindSat C-band number of bins in the CDF: not defined')
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
       
       allocate(WindSatCsm_struc(n)%smobs(182))
       allocate(WindSatCsm_struc(n)%tmobs(182))

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
    ios = 0 
    count = 1
    do n=1,LIS_rc%nnest

       ftn = LIS_getNextUnitNumber()
       write(LIS_logunit,*) 'Reading Windsat C-band data ', synsmobsdir
       open(ftn,file=synsmobsdir,form='formatted')
       do while(ios.eq.0) 
          read(ftn,*,iostat=ios) yr, mo, da, hr, mn, ss
          read(ftn,*,iostat=ios) smvalue
          call LIS_get_julss(yr, mo, da, hr, mn, ss, obstime)
          WindSatCsm_struc(n)%smobs(count) = smvalue
          WindSatCsm_struc(n)%tmobs(count) = obstime
          count = count + 1
       enddo

       call LIS_releaseUnitNumber(ftn)
    enddo

!--------------------------------------------------------------------------------
! Read the model and observation CDF data
!--------------------------------------------------------------------------------
    do n=1,LIS_rc%nnest
       if(WindSatCsm_struc(n)%scal.eq.1) then 

          allocate(WindSatCsm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), WindSatCsm_struc(n)%nbins))
          allocate(WindSatCsm_struc(n)%obs_xrange(&
               LIS_rc%ngrid(n), WindSatCsm_struc(n)%nbins))
          allocate(WindSatCsm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), WindSatCsm_struc(n)%nbins))
          allocate(WindSatCsm_struc(n)%obs_cdf(&
               LIS_rc%ngrid(n), WindSatCsm_struc(n)%nbins))

          allocate(xrange(LIS_rc%ngrid(n)))
          allocate(cdf(LIS_rc%ngrid(n)))


          ftn = LIS_getNextUnitNumber()
          write(LIS_logunit,*) 'Reading model CDF file ', modelcdffile(n)
          open(ftn,file=modelcdffile(n),form='unformatted')
          do i=1,WindSatCsm_struc(n)%nbins
             call LIS_readvar_gridded(ftn,n,xrange,dummy)
             call LIS_readvar_gridded(ftn,n,cdf,dummy)
             
             WindSatCsm_struc(n)%model_xrange(:,i) = xrange(:)
             WindSatCsm_struc(n)%model_cdf(:,i) = cdf(:)
          enddo
          
          call LIS_releaseUnitNumber(ftn)
          
          ftn = LIS_getNextUnitNumber()
          write(LIS_logunit,*) 'Reading obs CDF file ', obscdffile(n)
          open(ftn,file=obscdffile(n),form='unformatted')
          
          do i=1,WindSatCsm_struc(n)%nbins
             call LIS_readvar_gridded(ftn,n,xrange,dummy)
             call LIS_readvar_gridded(ftn,n,cdf,dummy)
             
             WindSatCsm_struc(n)%obs_xrange(:,i) = xrange(:)
             WindSatCsm_struc(n)%obs_cdf(:,i) = cdf(:)
          enddo
          
          call LIS_releaseUnitNumber(ftn)
       endif
    enddo

    write(LIS_logunit,*) 'Created ESMF States to hold the WindSat observations data'

  end subroutine WindSatCsm_setup
  
end module WindSatCsm_Mod

