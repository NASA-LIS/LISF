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
! !MODULE: SMAPNRTsm_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle 
! 
! !REVISION HISTORY: 
!  22 Aug 2016    Sujay Kumar; initial specification
!  1  Apr 2019  Yonghwan Kwon: Upated for reading monthy CDF for the current month
! 
module SMAPNRTsm_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SMAPNRTsm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: SMAPNRTsm_struc
!EOP
  type, public:: SMAPNRTsm_dec
     
     integer                :: useSsdevScal
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: smobs(:,:)
     real*8,   allocatable      :: smtime(:,:)

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
     logical                :: cdf_read_mon  !(for reading monthly CDF when
                                             !LIS_rc%da > 1 but the first model time step,
                                             !e.g., 4/29 13:00:00)
     integer                :: cdf_read_opt  ! 0: read all months at one time
                                             ! 1: read only the current month
     character(len=LIS_CONST_PATH_LEN) :: modelcdffile
     character(len=LIS_CONST_PATH_LEN) :: obscdffile
  end type SMAPNRTsm_dec
  
  type(SMAPNRTsm_dec),allocatable :: SMAPNRTsm_struc(:)
  
contains

!BOP
! 
! !ROUTINE: SMAPNRTsm_setup
! \label{SMAPNRTsm_setup}
! 
! !INTERFACE: 
  subroutine SMAPNRTsm_setup(k, OBS_State, OBS_Pert_State)
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
!   creation of data strctures required for handling NASASMAP 
!   soil moisture data. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    real, parameter        ::  minssdev =0.04 
    integer                ::  n,i,t,kk,jj
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  obsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real                   :: gridDesci(50)
    real, allocatable          :: ssdev(:)

    real, allocatable          ::  obserr(:,:)
    real, allocatable          ::  lobserr(:,:)
    integer                :: c,r
    real, allocatable          :: ssdev_grid(:,:)
    integer                :: ngrid

    allocate(SMAPNRTsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NRT) soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsdir,&
            rc=status)
       call LIS_verify(status, 'SMAP(NRT) soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status)
    enddo


    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NRT) soil moisture use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMAPNRTsm_struc(n)%useSsdevScal,&
            rc=status)
       call LIS_verify(status, 'SMAP(NRT) soil moisture use scaled standard deviation model: is missing')
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NRT) model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMAPNRTsm_struc(n)%modelcdffile,rc=status)   
       call LIS_verify(status, 'SMAP(NRT) model CDF file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NRT) observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMAPNRTsm_struc(n)%obscdffile,rc=status)     
       call LIS_verify(status, 'SMAP(NRT) observation CDF file: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "SMAP(NRT) soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SMAPNRTsm_struc(n)%nbins, rc=status)
       call LIS_verify(status, "SMAP(NRT) soil moisture number of bins in the CDF: not defined")
    enddo

   do n=1, LIS_rc%nnest
      SMAPNRTsm_struc(n)%cdf_read_mon = .false.

      call ESMF_ConfigFindLabel(LIS_config, "SMAP(NRT) CDF read option:", rc=status)    ! 0: read CDF for all months/year
                                                                                         ! 1: read CDF for current month
      call ESMF_ConfigGetAttribute(LIS_config, SMAPNRTsm_struc(n)%cdf_read_opt, rc=status)
      call LIS_verify(status, "SMAP(NRT) CDF read option: not defined")
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

    write(LIS_logunit,*)&
         '[INFO] read SMAP(NRT) soil moisture data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. NASASMAP
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
          SMAPNRTsm_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
         '[INFO] Created the States to hold the SMAP NASA observations data'
    do n=1,LIS_rc%nnest
       SMAPNRTsm_struc(n)%nc = 964
       SMAPNRTsm_struc(n)%nr = 406
       
       allocate(SMAPNRTsm_struc(n)%smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(SMAPNRTsm_struc(n)%smtime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       
       SMAPNRTsm_struc(n)%smtime = -1
      
    enddo
    
    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       ssdev = obs_pert%ssdev(1)
       
       if(LIS_rc%dascaloption(k).eq."CDF matching") then 
          call LIS_getCDFattributes(k,SMAPNRTsm_struc(n)%modelcdffile,&
               SMAPNRTsm_struc(n)%ntimes,ngrid)          
         
          if (SMAPNRTsm_struc(n)%cdf_read_opt.eq.0) then   
             allocate(SMAPNRTsm_struc(n)%model_xrange(&
                  LIS_rc%obs_ngrid(k), SMAPNRTsm_struc(n)%ntimes, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%obs_xrange(&
                  LIS_rc%obs_ngrid(k), SMAPNRTsm_struc(n)%ntimes, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%model_cdf(&
                  LIS_rc%obs_ngrid(k), SMAPNRTsm_struc(n)%ntimes, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%obs_cdf(&
                  LIS_rc%obs_ngrid(k), SMAPNRTsm_struc(n)%ntimes, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
                  SMAPNRTsm_struc(n)%ntimes))
             allocate(SMAPNRTsm_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
                  SMAPNRTsm_struc(n)%ntimes))
             allocate(SMAPNRTsm_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
                  SMAPNRTsm_struc(n)%ntimes))
             allocate(SMAPNRTsm_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
                  SMAPNRTsm_struc(n)%ntimes))
          else
             allocate(SMAPNRTsm_struc(n)%model_xrange(&
                  LIS_rc%obs_ngrid(k), 1, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%obs_xrange(&
                  LIS_rc%obs_ngrid(k), 1, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%model_cdf(&
                  LIS_rc%obs_ngrid(k), 1, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%obs_cdf(&
                  LIS_rc%obs_ngrid(k), 1, &
                  SMAPNRTsm_struc(n)%nbins))
             allocate(SMAPNRTsm_struc(n)%model_mu(LIS_rc%obs_ngrid(k),1))
             allocate(SMAPNRTsm_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),1))
             allocate(SMAPNRTsm_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),1))
             allocate(SMAPNRTsm_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),1))
          endif

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
         if (SMAPNRTsm_struc(n)%cdf_read_opt.eq.0) then  
          call LIS_readMeanSigmaData(n,k,&
               SMAPNRTsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               SMAPNRTsm_struc(n)%modelcdffile, &
               "SoilMoist",&
               SMAPNRTsm_struc(n)%model_mu,&
               SMAPNRTsm_struc(n)%model_sigma)  
          
          call LIS_readMeanSigmaData(n,k,&
               SMAPNRTsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               SMAPNRTsm_struc(n)%obscdffile, &
               "SoilMoist",&
               SMAPNRTsm_struc(n)%obs_mu,&
               SMAPNRTsm_struc(n)%obs_sigma)      
          
          call LIS_readCDFdata(n,k,&
               SMAPNRTsm_struc(n)%nbins,&
               SMAPNRTsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               SMAPNRTsm_struc(n)%modelcdffile, &
               "SoilMoist",&
               SMAPNRTsm_struc(n)%model_xrange,&
               SMAPNRTsm_struc(n)%model_cdf)      
          
          call LIS_readCDFdata(n,k,&
               SMAPNRTsm_struc(n)%nbins,&
               SMAPNRTsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               SMAPNRTsm_struc(n)%obscdffile, &
               "SoilMoist",&
               SMAPNRTsm_struc(n)%obs_xrange,&
               SMAPNRTsm_struc(n)%obs_cdf)    
          
          if(SMAPNRTsm_struc(n)%useSsdevScal.eq.1) then 
             if(SMAPNRTsm_struc(n)%ntimes.eq.1) then 
                jj = 1
             else
                jj = LIS_rc%mo
             endif
             do t=1,LIS_rc%obs_ngrid(k)
                if(SMAPNRTsm_struc(n)%obs_sigma(t,jj).gt.0) then 
                   ssdev(t) = ssdev(t)*SMAPNRTsm_struc(n)%model_sigma(t,jj)/&
                        SMAPNRTsm_struc(n)%obs_sigma(t,jj)
                   !                c = LIS_domain(n)%grid(t)%col
                   !                r = LIS_domain(n)%grid(t)%row
                   !                ssdev_grid(c,r) = ssdev(t) 
                   if(ssdev(t).lt.minssdev) then 
                      ssdev(t) = minssdev
                   endif
                endif
             enddo
             
!          open(100,file='ssdev.bin',form='unformatted')
!          write(100) ssdev_grid
!          close(100)
!          stop
          endif
         endif   
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
       endif
       if(LIS_rc%obs_ngrid(k).gt.0) then 
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
          call LIS_verify(status)
       endif

       deallocate(ssdev)
       
    enddo
    
    do n=1,LIS_rc%nnest
       gridDesci = 0
       gridDesci(1) = 9
       gridDesci(2) = 964
       gridDesci(3) = 406
       gridDesci(9) = 4 !M36 grid
       gridDesci(20) = 64
       gridDesci(10) = 0.36 
       gridDesci(11) = 1 !for the global switch

       allocate(SMAPNRTsm_struc(n)%n11(LIS_rc%obs_lnc(k)*&
            LIS_rc%obs_lnr(k)))
       allocate(SMAPNRTsm_struc(n)%rlat(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(SMAPNRTsm_struc(n)%rlon(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       
       call neighbor_interp_input_withgrid(gridDesci,&
            LIS_rc%obs_gridDesc(k,:),&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
            SMAPNRTsm_struc(n)%rlat, &
            SMAPNRTsm_struc(n)%rlon, &
            SMAPNRTsm_struc(n)%n11)
       
       call LIS_registerAlarm("SMAP NRT read alarm",&
            3600.0, 3600.0)

       SMAPNRTsm_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine SMAPNRTsm_setup
end module SMAPNRTsm_Mod
