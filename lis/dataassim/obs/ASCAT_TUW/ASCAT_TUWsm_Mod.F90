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
! !MODULE: ASCAT_TUWsm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle 
! 
! !REVISION HISTORY: 
!  8 May 2013    Sujay Kumar; initial specification
! 
module ASCAT_TUWsm_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ASCAT_TUWsm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ASCAT_TUWsm_struc
!EOP
  type, public:: ASCAT_TUWsm_dec 
     
     logical                :: startMode
     integer                :: useSsdevScal
     integer                :: nc
     integer                :: nr
     type(proj_info)        :: ascattuwproj

     integer, allocatable       :: n11(:)
     real,    allocatable       :: model_xrange(:,:)
     real,    allocatable       :: obs_xrange(:,:)
     real,    allocatable       :: model_cdf(:,:)
     real,    allocatable       :: obs_cdf(:,:)
     real,    allocatable       :: model_mu(:)
     real,    allocatable       :: obs_mu(:)
     real,    allocatable       :: model_sigma(:)
     real,    allocatable       :: obs_sigma(:)

     integer                :: nbins

  end type ASCAT_TUWsm_dec
  
  type(ASCAT_TUWsm_dec),allocatable :: ASCAT_TUWsm_struc(:)
  
contains

!BOP
! 
! !ROUTINE: ASCAT_TUWsm_setup
! \label{ASCAT_TUWsm_setup}
! 
! !INTERFACE: 
  subroutine ASCAT_TUWsm_setup(k, OBS_State, OBS_Pert_State)
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
!   creation of data strctures required for handling RT SMOPS 
!   soil moisture data. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    real, parameter        ::  minssdev = 0.04
    integer                ::  n,i,t,kk
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  ascattuwsmobsdir
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

!    integer                :: c,r
!    real, allocatable          :: ssdev_grid(:,:)


    allocate(ASCAT_TUWsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"ASCAT (TUW) soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ascattuwsmobsdir,&
            rc=status)
       call LIS_verify(status, 'ASCAT (TUW) soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            ascattuwsmobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ASCAT (TUW) use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ASCAT_TUWsm_struc(n)%useSsdevScal,rc=status)
       call LIS_verify(status, 'ASCAT (TUW) use scaled standard deviation model: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ASCAT (TUW) model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
       call LIS_verify(status, 'ASCAT (TUW) model CDF file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ASCAT (TUW) observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
       call LIS_verify(status, 'ASCAT (TUW) observation CDF file: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "ASCAT (TUW) soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ASCAT_TUWsm_struc(n)%nbins, rc=status)
       call LIS_verify(status, "ASCAT (TUW) soil moisture number of bins in the CDF: not defined")
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

    write(LIS_logunit,*)'read ASCAT (TUW) soil moisture data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. ASCAT (TUW)
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
         'Created the States to hold the ASCAT (TUW) observations data'
    do n=1,LIS_rc%nnest
       ASCAT_TUWsm_struc(n)%nc = 1440
       ASCAT_TUWsm_struc(n)%nr = 600

       call map_set(PROJ_LATLON, -59.875, -179.875, &
            0.0, 0.25, 0.25, 0.0, & 
            ASCAT_TUWsm_struc(n)%nc, ASCAT_TUWsm_struc(n)%nr, &
            ASCAT_TUWsm_struc(n)%ascattuwproj)
    enddo
    
    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%ngrid(n)))
       ssdev = obs_pert%ssdev(1)
       
       allocate(ASCAT_TUWsm_struc(n)%model_xrange(&
            LIS_rc%ngrid(n), ASCAT_TUWsm_struc(n)%nbins))
       allocate(ASCAT_TUWsm_struc(n)%obs_xrange(&
            LIS_rc%ngrid(n), ASCAT_TUWsm_struc(n)%nbins))
       allocate(ASCAT_TUWsm_struc(n)%model_cdf(&
            LIS_rc%ngrid(n), ASCAT_TUWsm_struc(n)%nbins))
       allocate(ASCAT_TUWsm_struc(n)%obs_cdf(&
            LIS_rc%ngrid(n), ASCAT_TUWsm_struc(n)%nbins))
       allocate(ASCAT_TUWsm_struc(n)%model_mu(LIS_rc%ngrid(n)))
       allocate(ASCAT_TUWsm_struc(n)%model_sigma(LIS_rc%ngrid(n)))
       allocate(ASCAT_TUWsm_struc(n)%obs_mu(LIS_rc%ngrid(n)))
       allocate(ASCAT_TUWsm_struc(n)%obs_sigma(LIS_rc%ngrid(n)))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
       call LIS_readMeanSigmaData(n,&
            modelcdffile(n), &
            "SoilMoist",&
            ASCAT_TUWsm_struc(n)%model_mu,&
            ASCAT_TUWsm_struc(n)%model_sigma)
       
       call LIS_readMeanSigmaData(n,&
            obscdffile(n), &
            "SoilMoist",&
            ASCAT_TUWsm_struc(n)%obs_mu,&
            ASCAT_TUWsm_struc(n)%obs_sigma)
       
       call LIS_readCDFdata(n,&
            ASCAT_TUWsm_struc(n)%nbins,&
            modelcdffile(n), &
            "SoilMoist",&
            ASCAT_TUWsm_struc(n)%model_xrange,&
            ASCAT_TUWsm_struc(n)%model_cdf)
       
       call LIS_readCDFdata(n,&
            ASCAT_TUWsm_struc(n)%nbins,&
            obscdffile(n), &
            "SoilMoist",&
            ASCAT_TUWsm_struc(n)%obs_xrange,&
            ASCAT_TUWsm_struc(n)%obs_cdf)
       
!       allocate(ssdev_grid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!       ssdev_grid = LIS_rc%udef

       if(ASCAT_TUWsm_struc(n)%useSsdevScal.eq.1) then 
          do t=1,LIS_rc%ngrid(n)
             if(ASCAT_TUWsm_struc(n)%obs_sigma(t).gt.0) then 
                !             print*, t, ssdev(t), ASCAT_TUWsm_struc(n)%model_sigma(t), & 
!                  ASCAT_TUWsm_struc(n)%obs_sigma(t)
                ssdev(t) = ssdev(t)*ASCAT_TUWsm_struc(n)%model_sigma(t)/&
                     ASCAT_TUWsm_struc(n)%obs_sigma(t)
!             c = LIS_domain(n)%grid(t)%col
!             r = LIS_domain(n)%grid(t)%row
!             ssdev_grid(c,r) = ssdev(t)             
                if(ssdev(t).gt.minssdev) then 
                   ssdev(t) = minssdev
                endif
             endif
          enddo
       endif
!       open(100,file='ssdev.bin',form='unformatted')
!       write(100) ssdev_grid
!       close(100)
!       stop
      
       if(LIS_rc%ngrid(n).gt.0) then 
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
          call LIS_verify(status)
       endif

       deallocate(ssdev)
    enddo
    
    do n=1,LIS_rc%nnest
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 1440
       gridDesci(3) = 600
       gridDesci(4) = -59.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
       
       allocate(ASCAT_TUWsm_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       
       call neighbor_interp_input(n, gridDesci,&
            ASCAT_TUWsm_struc(n)%n11)
       
       call LIS_registerAlarm("ASCAT TUW read alarm",&
            3600.0,3600.0)
       ASCAT_TUWsm_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine ASCAT_TUWsm_setup
end module ASCAT_TUWsm_Mod
