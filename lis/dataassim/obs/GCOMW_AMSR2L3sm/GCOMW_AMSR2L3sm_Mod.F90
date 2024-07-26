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
! !MODULE: GCOMW_AMSR2L3sm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle AMSR2 retrievals
! 
! !REVISION HISTORY: 
!  12 Jan 15    Sujay Kumar; Initial version
! 
module GCOMW_AMSR2L3sm_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: GCOMW_AMSR2L3sm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: GCOMW_AMSR2L3sm_struc
!EOP
  type, public:: GCOMW_AMSR2L3sm_dec
     
     logical                :: startMode
     integer                :: useSsdevScal
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: smobs(:,:)
     real,     allocatable      :: smtime(:,:)

     real                   :: datares
     real                   :: ssdev_inp
     integer                :: amsr2nc, amsr2nr
     type(proj_info)        :: amsr2proj
     integer, allocatable       :: n11(:)
     integer, allocatable       :: n12(:)
     integer, allocatable       :: n21(:)
     integer, allocatable       :: n22(:)
     real,    allocatable       :: w11(:)
     real,    allocatable       :: w12(:)
     real,    allocatable       :: w21(:)
     real,    allocatable       :: w22(:)

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

  end type GCOMW_AMSR2L3sm_dec
  
  type(GCOMW_AMSR2L3sm_dec),allocatable :: GCOMW_AMSR2L3sm_struc(:)
  
contains

!BOP
! 
! !ROUTINE: GCOMW_AMSR2L3sm_setup
! \label{GCOMW_AMSR2L3sm_setup}
! 
! !INTERFACE: 
  subroutine GCOMW_AMSR2L3sm_setup(k, OBS_State, OBS_Pert_State)
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
!   creation of data strctures required for handling AMSR2 
!   AMSR-E soil moisture data. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP

    real, parameter        ::  minssdev = 0.01
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
    integer                :: jj,ngrid

    allocate(GCOMW_AMSR2L3sm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,amsresmobsdir,&
            rc=status)
       call LIS_verify(status, 'AMSR2(GCOMW) soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            amsresmobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) scale observations:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GCOMW_AMSR2L3sm_struc(n)%scal, &
            rc=status)
       call LIS_verify(status, "AMSR2(GCOMW) scale observations: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GCOMW_AMSR2L3sm_struc(n)%useSsdevScal, &
            rc=status)
       call LIS_verify(status, "AMSR2(GCOMW) use scaled standard deviation model: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(GCOMW_AMSR2L3sm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'AMSR2(GCOMW) model CDF file: not defined')
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if(GCOMW_AMSR2L3sm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
          call LIS_verify(status, 'AMSR2(GCOMW) observation CDF file: not defined')
       endif
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "AMSR2(GCOMW) soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       if(GCOMW_AMSR2L3sm_struc(n)%scal.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,GCOMW_AMSR2L3sm_struc(n)%nbins, rc=status)
          call LIS_verify(status, "AMSR2(GCOMW) soil moisture number of bins in the CDF: not defined")
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

    write(LIS_logunit,*)'read AMSR2(GCOMW) soil moisture data specifications'       

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
          GCOMW_AMSR2L3sm_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
         'Created the States to hold the AMSR2(GCOMW) observations data'
    do n=1,LIS_rc%nnest
       GCOMW_AMSR2L3sm_struc(n)%nc = 3600
       GCOMW_AMSR2L3sm_struc(n)%nr = 1800
       allocate(GCOMW_AMSR2L3sm_struc(n)%smobs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(GCOMW_AMSR2L3sm_struc(n)%smtime(LIS_rc%lnc(n),LIS_rc%lnr(n)))

    enddo
    
    do n=1,LIS_rc%nnest
       if(GCOMW_AMSR2L3sm_struc(n)%scal.eq.1) then 

          call LIS_getCDFattributes(modelcdffile(n),&
               GCOMW_AMSR2L3sm_struc(n)%ntimes,ngrid)

          allocate(ssdev(LIS_rc%ngrid(n)))
          ssdev = obs_pert%ssdev(1)

          allocate(GCOMW_AMSR2L3sm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n),  GCOMW_AMSR2L3sm_struc(n)%ntimes,&
               GCOMW_AMSR2L3sm_struc(n)%nbins))
          allocate(GCOMW_AMSR2L3sm_struc(n)%obs_xrange(&
               LIS_rc%ngrid(n),  GCOMW_AMSR2L3sm_struc(n)%ntimes, &
               GCOMW_AMSR2L3sm_struc(n)%nbins))
          allocate(GCOMW_AMSR2L3sm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n),  GCOMW_AMSR2L3sm_struc(n)%ntimes, &
               GCOMW_AMSR2L3sm_struc(n)%nbins))
          allocate(GCOMW_AMSR2L3sm_struc(n)%obs_cdf(&
               LIS_rc%ngrid(n),  GCOMW_AMSR2L3sm_struc(n)%ntimes, &
               GCOMW_AMSR2L3sm_struc(n)%nbins))
          allocate(GCOMW_AMSR2L3sm_struc(n)%model_mu(LIS_rc%ngrid(n),&
               GCOMW_AMSR2L3sm_struc(n)%ntimes))
          allocate(GCOMW_AMSR2L3sm_struc(n)%model_sigma(LIS_rc%ngrid(n),&
               GCOMW_AMSR2L3sm_struc(n)%ntimes))
          allocate(GCOMW_AMSR2L3sm_struc(n)%obs_mu(LIS_rc%ngrid(n),&
               GCOMW_AMSR2L3sm_struc(n)%ntimes))
          allocate(GCOMW_AMSR2L3sm_struc(n)%obs_sigma(LIS_rc%ngrid(n),&
               GCOMW_AMSR2L3sm_struc(n)%ntimes))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
          call LIS_readMeanSigmaData(n,&
               GCOMW_AMSR2L3sm_struc(n)%ntimes,& 
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               GCOMW_AMSR2L3sm_struc(n)%model_mu,&
               GCOMW_AMSR2L3sm_struc(n)%model_sigma)

          call LIS_readMeanSigmaData(n,&
               GCOMW_AMSR2L3sm_struc(n)%ntimes,& 
               ngrid, &
               obscdffile(n), &
               "SoilMoist",&
               GCOMW_AMSR2L3sm_struc(n)%obs_mu,&
               GCOMW_AMSR2L3sm_struc(n)%obs_sigma)

          call LIS_readCDFdata(n,&
               GCOMW_AMSR2L3sm_struc(n)%nbins,&
               GCOMW_AMSR2L3sm_struc(n)%ntimes,& 
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               GCOMW_AMSR2L3sm_struc(n)%model_xrange,&
               GCOMW_AMSR2L3sm_struc(n)%model_cdf)

          call LIS_readCDFdata(n,&
               GCOMW_AMSR2L3sm_struc(n)%nbins,&
               GCOMW_AMSR2L3sm_struc(n)%ntimes,&
               ngrid, &
               obscdffile(n), &
               "SoilMoist",&
               GCOMW_AMSR2L3sm_struc(n)%obs_xrange,&
               GCOMW_AMSR2L3sm_struc(n)%obs_cdf)

          if(GCOMW_AMSR2L3sm_struc(n)%useSsdevScal.eq.1) then 
             if(GCOMW_AMSR2L3sm_struc(n)%ntimes.eq.1) then 
                jj = 1
             else
                jj = LIS_rc%mo
             endif
             do t=1,LIS_rc%ngrid(n)
                if(GCOMW_AMSR2L3sm_struc(n)%obs_sigma(t,jj).gt.0) then 
                   ssdev(t) = ssdev(t)*GCOMW_AMSR2L3sm_struc(n)%model_sigma(t,jj)/&
                        GCOMW_AMSR2L3sm_struc(n)%obs_sigma(t,jj)

                   if(ssdev(t).gt.maxssdev) ssdev(t) = maxssdev
                   if(ssdev(t).lt.minssdev) then 
                      ssdev(t) = minssdev
                   endif
                endif
             enddo
          endif

#if 0           
          allocate(obserr(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          obserr = -9999.0

!          lobserr(:,:) = obserr(&
!               LIS_ews_halo_ind(n,LIS_localPet+1):&         
!               LIS_ewe_halo_ind(n,LIS_localPet+1), &
!               LIS_nss_halo_ind(n,LIS_localPet+1): &
!               LIS_nse_halo_ind(n,LIS_localPet+1))

          do r=1,LIS_rc%lnr(n)
             do c=1,LIS_rc%lnc(n)
                if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                   obserr(c,r)  =  ssdev(LIS_domain(n)%gindex(c,r)) 
                   
                endif
             enddo
          enddo

!          do r=1,LIS_rc%lnr(n)
!             do c=1,LIS_rc%lnc(n)
!                if(LIS_domain(n)%gindex(c,r).ne.-1) then 
!                   lobserr(c,r) = ssdev(LIS_domain(n)%gindex(c,r)) 
!                   
!                endif
!             enddo
!          enddo
          open(100,file='test.bin',form='unformatted')
          write(100) obserr
          close(100)
          stop
          deallocate(obserr)
#endif
         
          if(LIS_rc%ngrid(n).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
             call LIS_verify(status)
          endif

          deallocate(ssdev)
       endif
    enddo

    do n=1,LIS_rc%nnest
       GCOMW_AMSR2L3sm_struc(n)%amsr2nc = 3600
       GCOMW_AMSR2L3sm_struc(n)%amsr2nr = 1800

       call map_set(PROJ_LATLON, -89.95,-179.95,&
            0.0, 0.10,0.10, 0.0,&
            GCOMW_AMSR2L3sm_struc(n)%amsr2nc,GCOMW_AMSR2L3sm_struc(n)%amsr2nr,&
            GCOMW_AMSR2L3sm_struc(n)%amsr2proj)
       
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 3600
       gridDesci(3) = 1800
       gridDesci(4) = -89.95
       gridDesci(5) = -179.95
       gridDesci(6) = 128
       gridDesci(7) = 89.95
       gridDesci(8) = 179.95
       gridDesci(9) = 0.10
       gridDesci(10) = 0.10
       gridDesci(20) = 64
       
       GCOMW_AMSR2L3sm_struc(n)%datares = 0.10
       
       if(LIS_isatAfinerResolution(n,GCOMW_AMSR2L3sm_struc(n)%datares)) then 
          allocate(GCOMW_AMSR2L3sm_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3sm_struc(n)%n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3sm_struc(n)%n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3sm_struc(n)%n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          allocate(GCOMW_AMSR2L3sm_struc(n)%w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3sm_struc(n)%w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3sm_struc(n)%w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3sm_struc(n)%w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call bilinear_interp_input(n, gridDesci, &
               GCOMW_AMSR2L3sm_struc(n)%n11, &
               GCOMW_AMSR2L3sm_struc(n)%n12, GCOMW_AMSR2L3sm_struc(n)%n21, &
               GCOMW_AMSR2L3sm_struc(n)%n22, GCOMW_AMSR2L3sm_struc(n)%w11, &
               GCOMW_AMSR2L3sm_struc(n)%w12, GCOMW_AMSR2L3sm_struc(n)%w21, &
               GCOMW_AMSR2L3sm_struc(n)%w22)
       
       else
          allocate(GCOMW_AMSR2L3sm_struc(n)%n11(&
               GCOMW_AMSR2L3sm_struc(n)%amsr2nc*GCOMW_AMSR2L3sm_struc(n)%amsr2nr))

          call upscaleByAveraging_input(gridDesci,&
               LIS_rc%gridDesc(n,:),&
               GCOMW_AMSR2L3sm_struc(n)%amsr2nc*GCOMW_AMSR2L3sm_struc(n)%amsr2nr,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n),&
               GCOMW_AMSR2L3sm_struc(n)%n11)
       endif
       call LIS_registerAlarm("AMSR2(GCOMW) read alarm",&
            86400.0, 86400.0)
       GCOMW_AMSR2L3sm_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine GCOMW_AMSR2L3sm_setup

end module GCOMW_AMSR2L3sm_Mod
