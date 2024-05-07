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
! !MODULE: SMOPS_AMSR2sm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle NOAA/NESDIS-based AMSR2 Soil Moisture (SM) Files (aka, SMOPS)
!
!  There are three versions of the SMOPS datasets.  According to the
!  use by the 557th Weather Wing:
!
!                         version_1.3 <  2016-10-31T12:00:00
!  2016-10-31T12:00:00 <= version_2.0 <  2017-08-24T12:00:00
!                         version_3.0 >= 2017-08-24T12:00:00
! 
! !REVISION HISTORY: 
!  8 May 2013    Sujay Kumar; initial specification
!  28 Sep 2017: Mahdi Navari; Updated to read AMSR2 from SMOPS V3
! 
module SMOPS_AMSR2sm_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SMOPS_AMSR2sm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: SMOPS_AMSR2sm_struc
!EOP
  type, public:: SMOPS_AMSR2sm_dec
     
     integer                :: useRealtime
     integer                :: useAMSR2
     integer                :: useSsdevScal
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: smobs(:,:)
     real,     allocatable      :: smtime(:,:)

     real                       :: ssdev_inp
     integer, allocatable       :: n11(:)
     real, allocatable          :: rlat(:)
     real, allocatable          :: rlon(:)

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
     real*8                 :: version2_time, version3_time
     character(len=5)       :: conv
  end type SMOPS_AMSR2sm_dec
  
  type(SMOPS_AMSR2sm_dec),allocatable :: SMOPS_AMSR2sm_struc(:)
  
contains

!BOP
! 
! !ROUTINE: SMOPS_AMSR2sm_setup
! \label{SMOPS_AMSR2sm_setup}
! 
! !INTERFACE: 
  subroutine SMOPS_AMSR2sm_setup(k, OBS_State, OBS_Pert_State)
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
!   creation of data strctures required for handling SMOPS 
!   soil moisture data. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    real, parameter        ::  minssdev =0.02
    integer                ::  n,i,t,kk,jj
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  rtsmopssmobsdir
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

    real, allocatable          ::  obserr(:,:)
    real, allocatable          ::  lobserr(:,:)
    integer                :: c,r
    real, allocatable          :: ssdev_grid(:,:)
    integer                :: ngrid
    integer                 :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real                    :: upgmt

    allocate(SMOPS_AMSR2sm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SMOPS soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,rtsmopssmobsdir,&
            rc=status)
       call LIS_verify(status, 'SMOPS soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            rtsmopssmobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOPS use realtime data:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOPS_AMSR2sm_struc(n)%useRealtime,&
            rc=status)
       call LIS_verify(status, 'SMOPS use realtime data: is missing')

    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOPS soil moisture use AMSR2 data:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOPS_AMSR2sm_struc(n)%useAMSR2,&
            rc=status)
       call LIS_verify(status, 'SMOPS soil moisture use AMSR2 data: is missing')
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOPS soil moisture use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOPS_AMSR2sm_struc(n)%useSsdevScal,&
            rc=status)
       call LIS_verify(status, 'SMOPS soil moisture use scaled standard deviation model: is missing')
       
    enddo
    call ESMF_ConfigFindLabel(LIS_config,"SMOPS model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
       call LIS_verify(status, 'SMOPS model CDF file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMOPS observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
       call LIS_verify(status, 'SMOPS observation CDF file: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "SMOPS soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOPS_AMSR2sm_struc(n)%nbins, rc=status)
       call LIS_verify(status, "SMOPS soil moisture number of bins in the CDF: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config, "SMOPS naming convention:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMOPS_AMSR2sm_struc(n)%conv, rc=status)
       call LIS_verify(status, "SMOPS naming convention: not defined")
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

    write(LIS_logunit,*)'[INFO] read SMOPS soil moisture data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. SMOPS
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
          SMOPS_AMSR2sm_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
         '[INFO] Created the States to hold the SMOPS observations data'
    do n=1,LIS_rc%nnest
       SMOPS_AMSR2sm_struc(n)%nc = 1440
       SMOPS_AMSR2sm_struc(n)%nr = 720
       allocate(SMOPS_AMSR2sm_struc(n)%smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(SMOPS_AMSR2sm_struc(n)%smtime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))

    enddo
    
    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       ssdev = obs_pert%ssdev(1)
       
       call LIS_getCDFattributes(k,modelcdffile(n),&
            SMOPS_AMSR2sm_struc(n)%ntimes,ngrid)
       
       allocate(SMOPS_AMSR2sm_struc(n)%model_xrange(&
            LIS_rc%obs_ngrid(k), SMOPS_AMSR2sm_struc(n)%ntimes, &
            SMOPS_AMSR2sm_struc(n)%nbins))
       allocate(SMOPS_AMSR2sm_struc(n)%obs_xrange(&
            LIS_rc%obs_ngrid(k), SMOPS_AMSR2sm_struc(n)%ntimes, &
            SMOPS_AMSR2sm_struc(n)%nbins))
       allocate(SMOPS_AMSR2sm_struc(n)%model_cdf(&
            LIS_rc%obs_ngrid(k), SMOPS_AMSR2sm_struc(n)%ntimes, &
            SMOPS_AMSR2sm_struc(n)%nbins))
       allocate(SMOPS_AMSR2sm_struc(n)%obs_cdf(&
            LIS_rc%obs_ngrid(k), SMOPS_AMSR2sm_struc(n)%ntimes, &
            SMOPS_AMSR2sm_struc(n)%nbins))
       allocate(SMOPS_AMSR2sm_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
            SMOPS_AMSR2sm_struc(n)%ntimes))
       allocate(SMOPS_AMSR2sm_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
            SMOPS_AMSR2sm_struc(n)%ntimes))
       allocate(SMOPS_AMSR2sm_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
            SMOPS_AMSR2sm_struc(n)%ntimes))
       allocate(SMOPS_AMSR2sm_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
            SMOPS_AMSR2sm_struc(n)%ntimes))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
       call LIS_readMeanSigmaData(n,k,&
            SMOPS_AMSR2sm_struc(n)%ntimes,&
            LIS_rc%obs_ngrid(k), &
            modelcdffile(n), &
            "SoilMoist",&
            SMOPS_AMSR2sm_struc(n)%model_mu,&
            SMOPS_AMSR2sm_struc(n)%model_sigma)
       
       call LIS_readMeanSigmaData(n,k,&
            SMOPS_AMSR2sm_struc(n)%ntimes,&
            LIS_rc%obs_ngrid(k), &
            obscdffile(n), &
            "SoilMoist",&
            SMOPS_AMSR2sm_struc(n)%obs_mu,&
            SMOPS_AMSR2sm_struc(n)%obs_sigma)
       
       call LIS_readCDFdata(n,k,&
            SMOPS_AMSR2sm_struc(n)%nbins,&
            SMOPS_AMSR2sm_struc(n)%ntimes,&
            LIS_rc%obs_ngrid(k), &
            modelcdffile(n), &
            "SoilMoist",&
            SMOPS_AMSR2sm_struc(n)%model_xrange,&
            SMOPS_AMSR2sm_struc(n)%model_cdf)
       
       call LIS_readCDFdata(n,k,&
            SMOPS_AMSR2sm_struc(n)%nbins,&
            SMOPS_AMSR2sm_struc(n)%ntimes,&
            LIS_rc%obs_ngrid(k), &
            obscdffile(n), &
            "SoilMoist",&
            SMOPS_AMSR2sm_struc(n)%obs_xrange,&
            SMOPS_AMSR2sm_struc(n)%obs_cdf)
       

       if(SMOPS_AMSR2sm_struc(n)%useSsdevScal.eq.1) then 
          if(SMOPS_AMSR2sm_struc(n)%ntimes.eq.1) then 
             jj = 1
          else
             jj = LIS_rc%mo
          endif
          do t=1,LIS_rc%obs_ngrid(k)
             if(SMOPS_AMSR2sm_struc(n)%obs_sigma(t,jj).gt.0) then 
                ssdev(t) = ssdev(t)*SMOPS_AMSR2sm_struc(n)%model_sigma(t,jj)/&
                     SMOPS_AMSR2sm_struc(n)%obs_sigma(t,jj)
!                c = LIS_domain(n)%grid(t)%col
!                r = LIS_domain(n)%grid(t)%row
!                ssdev_grid(c,r) = ssdev(t) 
                if(ssdev(t).lt.minssdev) then 
                   ssdev(t) = minssdev
                endif
                print*, t, ssdev(t)

             endif
          enddo
          
!          open(100,file='ssdev.bin',form='unformatted')
!          write(100) ssdev_grid
!          close(100)
!          stop
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
          open(100,file='obs_err.bin',form='unformatted')
          write(100) obserr
          close(100)
          stop
          deallocate(obserr)
#endif
#if 0           
          allocate(obserr(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(lobserr(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
          obserr = -9999.0
          lobserr = -9999.0
          open(100,file='ssdev.bin',form='unformatted',access='direct',&
               recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
          read(100,rec=1) obserr
          close(100)
!          stop
          lobserr(:,:) = obserr(&
               LIS_ews_halo_ind(n,LIS_localPet+1):&         
               LIS_ewe_halo_ind(n,LIS_localPet+1), &
               LIS_nss_halo_ind(n,LIS_localPet+1): &
               LIS_nse_halo_ind(n,LIS_localPet+1))

          do r=1,LIS_rc%obs_lnr(k)
             do c=1,LIS_rc%obs_lnc(k)
                if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                   if(lobserr(c,r).gt.0.001) then 
                      ssdev(LIS_domain(n)%gindex(c,r))  = &
                           ssdev(LIS_domain(n)%gindex(c,r)) *2
                   endif                    
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
!          open(100,file='obs_err.bin',form='unformatted')
!          write(100) obserr
!          close(100)
!          stop
          deallocate(obserr)
          deallocate(lobserr)

#endif

       if(LIS_rc%obs_ngrid(k).gt.0) then 
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
          call LIS_verify(status)
       endif

       deallocate(ssdev)
       
    enddo
    
    do n=1,LIS_rc%nnest
       yr1 = 2016; mo1 = 10; da1 = 31; hr1 = 12; mn1 = 0; ss1 = 0
       call LIS_date2time(SMOPS_AMSR2sm_struc(n)%version2_time,updoy,upgmt, &
                          yr1,mo1,da1,hr1,mn1,ss1)

       yr1 = 2017; mo1 = 8; da1 = 24; hr1 = 12; mn1 = 0; ss1 = 0
       call LIS_date2time(SMOPS_AMSR2sm_struc(n)%version3_time,updoy,upgmt, &
                          yr1,mo1,da1,hr1,mn1,ss1)

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
       
       allocate(SMOPS_AMSR2sm_struc(n)%n11(SMOPS_AMSR2sm_struc(n)%nc*SMOPS_AMSR2sm_struc(n)%nr))
       allocate(SMOPS_AMSR2sm_struc(n)%rlat(&
            SMOPS_AMSR2sm_struc(n)%nc*SMOPS_AMSR2sm_struc(n)%nr))
       allocate(SMOPS_AMSR2sm_struc(n)%rlon(&
            SMOPS_AMSR2sm_struc(n)%nc*SMOPS_AMSR2sm_struc(n)%nr))
       
       call neighbor_interp_input_withgrid(gridDesci,&
            LIS_rc%obs_gridDesc(k,:),&
            SMOPS_AMSR2sm_struc(n)%nc*&
            SMOPS_AMSR2sm_struc(n)%nr,&
            SMOPS_AMSR2sm_struc(n)%rlat, &
            SMOPS_AMSR2sm_struc(n)%rlon, &
            SMOPS_AMSR2sm_struc(n)%n11)
       
       if ( SMOPS_AMSR2sm_struc(n)%useRealtime == 1 ) then
          call LIS_registerAlarm("SMOPS read alarm",&
               10800.0, 10800.0)
       else
          call LIS_registerAlarm("SMOPS read alarm",&
               86400.0, 86400.0)
       endif
       SMOPS_AMSR2sm_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine SMOPS_AMSR2sm_setup
end module SMOPS_AMSR2sm_Mod
