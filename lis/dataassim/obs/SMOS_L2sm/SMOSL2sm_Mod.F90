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
! !MODULE: SMOSL2sm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle SMOS L2 soil moisture retrievals
! 
!   
! !REVISION HISTORY: 
!  16 Dec 14    Sujay Kumar; Initial specification
! 
module SMOSL2sm_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SMOSL2sm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: SMOSL2sm_struc
!EOP
  type, public:: SMOSL2sm_dec
     character(len=LIS_CONST_PATH_LEN) :: odir
     logical                :: startMode
     integer                :: useSsdevScal
     integer                :: nc
     integer                :: nr
     real,            allocatable      :: smobs(:,:)
     type(ESMF_TIme), allocatable      :: smtime(:,:)

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
     integer                :: scal

  end type SMOSL2sm_dec
  
  type(SMOSL2sm_dec),allocatable :: SMOSL2sm_struc(:)
  
contains

!BOP
! 
! !ROUTINE: SMOSL2sm_setup
! \label{SMOSL2sm_setup}
! 
! !INTERFACE: 
  subroutine SMOSL2sm_setup(k, OBS_State, OBS_Pert_State)
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
    real, parameter        ::  minssdev = 0.01
    real, parameter        ::  maxssdev = 0.11
    integer                ::  n,i,t,kk,c,r,jj
    real, allocatable          ::  obserr(:,:)
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
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
    integer                    :: ngrid

    allocate(SMOSL2sm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SMOS L2 soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOSL2sm_struc(n)%odir,&
            rc=status)
       call LIS_verify(status, 'SMOS L2 soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            SMOSL2sm_struc(n)%odir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOS L2 scale observations:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOSL2sm_struc(n)%scal, &
            rc=status)
       call LIS_verify(status, "SMOS L2 scale observations: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOS L2 use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOSL2sm_struc(n)%useSsdevScal, &
            rc=status)
       call LIS_verify(status, "SMOS L2 use scaled standard deviation model: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOS L2 model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(SMOSL2sm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'SMOS L2 model CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOS L2 observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(SMOSL2sm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
          call LIS_verify(status, 'SMOS L2 observation CDF file: not defined')
       endif
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "SMOS L2 soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       if(SMOSL2sm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,SMOSL2sm_struc(n)%nbins, rc=status)
          call LIS_verify(status, "SMOS L2 soil moisture number of bins in the CDF: not defined")
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

    write(LIS_logunit,*)'read SMOS L2 soil moisture data specifications'       

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

! Set obs err to be uniform (will be rescaled later for each grid point). 
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
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status)

          if(LIS_rc%ngrid(n).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
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
         'Created the States to hold the SMOS L2 observations data'
    
    do n=1,LIS_rc%nnest
       if(SMOSL2sm_struc(n)%scal.eq.1) then 
          
          allocate(ssdev(LIS_rc%ngrid(n)))
          ssdev = obs_pert%ssdev(1)

          call LIS_getCDFattributes(modelcdffile(n),&
               SMOSL2sm_struc(n)%ntimes, ngrid)
          
          allocate(SMOSL2sm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), SMOSL2sm_struc(n)%ntimes, &
               SMOSL2sm_struc(n)%nbins))
          allocate(SMOSL2sm_struc(n)%obs_xrange(&
               LIS_rc%ngrid(n), SMOSL2sm_struc(n)%ntimes, &
               SMOSL2sm_struc(n)%nbins))
          allocate(SMOSL2sm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), SMOSL2sm_struc(n)%ntimes, &
               SMOSL2sm_struc(n)%nbins))
          allocate(SMOSL2sm_struc(n)%obs_cdf(&
               LIS_rc%ngrid(n), SMOSL2sm_struc(n)%ntimes, &
               SMOSL2sm_struc(n)%nbins))

          allocate(SMOSL2sm_struc(n)%model_mu(LIS_rc%ngrid(n),&
               SMOSL2sm_struc(n)%ntimes))
          allocate(SMOSL2sm_struc(n)%model_sigma(LIS_rc%ngrid(n),&
               SMOSL2sm_struc(n)%ntimes))
          allocate(SMOSL2sm_struc(n)%obs_mu(LIS_rc%ngrid(n),&
               SMOSL2sm_struc(n)%ntimes))
          allocate(SMOSL2sm_struc(n)%obs_sigma(LIS_rc%ngrid(n),&
               SMOSL2sm_struc(n)%ntimes))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
          call LIS_readMeanSigmaData(n,&
               SMOSL2sm_struc(n)%ntimes, &
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               SMOSL2sm_struc(n)%model_mu,&
               SMOSL2sm_struc(n)%model_sigma)

          call LIS_readMeanSigmaData(n,&
               SMOSL2sm_struc(n)%ntimes, &
               ngrid, &
               obscdffile(n), &
               "SoilMoist",&
               SMOSL2sm_struc(n)%obs_mu,&
               SMOSL2sm_struc(n)%obs_sigma)

          call LIS_readCDFdata(n,&
               SMOSL2sm_struc(n)%nbins,&
               SMOSL2sm_struc(n)%ntimes, &
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               SMOSL2sm_struc(n)%model_xrange,&
               SMOSL2sm_struc(n)%model_cdf)

          call LIS_readCDFdata(n,&
               SMOSL2sm_struc(n)%nbins,&
               SMOSL2sm_struc(n)%ntimes, &
               ngrid, &
               obscdffile(n), &
               "SoilMoist",&
               SMOSL2sm_struc(n)%obs_xrange,&
               SMOSL2sm_struc(n)%obs_cdf)

          if(SMOSL2sm_struc(n)%useSsdevScal.eq.1) then

             if(SMOSL2sm_struc(n)%ntimes.eq.1) then 
                jj = 1
             else
                jj = LIS_rc%mo
             endif

             do t=1,LIS_rc%ngrid(n)
                if(SMOSL2sm_struc(n)%obs_sigma(t,jj).gt.0) then 
                   ssdev(t) = ssdev(t)*SMOSL2sm_struc(n)%model_sigma(t,jj)/&
                        SMOSL2sm_struc(n)%obs_sigma(t,jj)

                   if(ssdev(t).gt.maxssdev) ssdev(t) = maxssdev
                   if(ssdev(t).lt.minssdev) then 
                      ssdev(t) = minssdev
                   endif
                endif
             enddo
          endif
         
          if(LIS_rc%ngrid(n).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
             call LIS_verify(status)
          endif

          deallocate(ssdev)
       endif
    enddo

    do n=1,LIS_rc%nnest
       allocate(SMOSL2sm_struc(n)%smobs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(SMOSL2sm_struc(n)%smtime(LIS_rc%lnc(n),LIS_rc%lnr(n)))

       call LIS_registerAlarm("SMOS L2 read alarm",&
            86400.0, 86400.0)
       SMOSL2sm_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine SMOSL2sm_setup
end module SMOSL2sm_Mod
