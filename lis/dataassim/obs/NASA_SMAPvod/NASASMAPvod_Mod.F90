!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: NASASMAPvod_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle the processing of SMAP vegetation optical depth 
!   retrievals
! 
! !REVISION HISTORY: 
!  28 Mar 2019    Sujay Kumar; initial specification
! 
module NASASMAPvod_Mod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: NASASMAPvod_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: NASASMAPvod_struc
!EOP
  type, public:: NASASMAPvod_dec
     
     integer                :: useSsdevScal
     integer                :: qcFlag
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: vodobs(:,:)
     real,     allocatable      :: vodtime(:,:)

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
     character*20           :: data_designation

     integer                :: nbins
     integer                :: ntimes

  end type NASASMAPvod_dec
  
  type(NASASMAPvod_dec),allocatable :: NASASMAPvod_struc(:)
  
contains

!BOP
! 
! !ROUTINE: NASASMAPvod_setup
! \label{NASASMAPvod_setup}
! 
! !INTERFACE: 
  subroutine NASASMAPvod_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use ESMF
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
!   creation of data strctures required for handling NASA SMAP 
!   vegetation optical depth data. The code expects
!   inputs of a reference CDF (of Leaf Area Index) and the
!   an observation CDF (of VOD).
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    real, parameter        ::  minssdev =0.001
    integer                ::  n,i,t,kk,jj
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character*100          ::  smapvoddir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real                   :: gridDesci(50)
    character*100          :: modelcdffile(LIS_rc%nnest)
    character*100          :: obscdffile(LIS_rc%nnest)
    real, allocatable          :: ssdev(:)

    real, allocatable          ::  obserr(:,:)
    real, allocatable          ::  lobserr(:,:)
    integer                :: c,r
    real, allocatable          :: ssdev_grid(:,:)
    integer                :: ngrid

    allocate(NASASMAPvod_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,&
         "SMAP(NASA) vegetation optical depth data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,smapvoddir,&
            rc=status)
       call LIS_verify(status, &
            'SMAP(NASA) vegetation optical depth data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            smapvoddir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "SMAP(NASA) vegetation optical depth data designation:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            NASASMAPvod_struc(n)%data_designation,&
            rc=status)
       call LIS_verify(status, &
            'SMAP(NASA) vegetation optical depth data designation: is missing')

    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "SMAP(NASA) vegetation optical depth use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,NASASMAPvod_struc(n)%useSsdevScal,&
            rc=status)
       call LIS_verify(status, 'SMAP(NASA) vegetation optical depth use scaled standard deviation model: is missing')
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "SMAP(NASA) vegetation optical depth apply SMAP QC flags:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,NASASMAPvod_struc(n)%qcFlag,&
            rc=status)
       call LIS_verify(status, &
            'SMAP(NASA) vegetation optical depth apply SMAP QC flags: is missing')
       
    enddo
    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NASA) reference LAI CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
       call LIS_verify(status, 'SMAP(NASA) reference LAI CDF file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NASA) observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
       call LIS_verify(status, 'SMAP(NASA) observation CDF file: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, &
         "SMAP(NASA) vegetation optical depth number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,NASASMAPvod_struc(n)%nbins, rc=status)
       call LIS_verify(status, &
            "SMAP(NASA) vegetation optical depth number of bins in the CDF: not defined")
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
         '[INFO] read SMAP(NASA) vegetation optical depth data specifications'       

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
          NASASMAPvod_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
       if(NASASMAPvod_struc(n)%data_designation.eq."SPL3SMP") then 
          NASASMAPvod_struc(n)%nc = 964
          NASASMAPvod_struc(n)%nr = 406
       elseif(NASASMAPvod_struc(n)%data_designation.eq."SPL3SMP_E") then 
          NASASMAPvod_struc(n)%nc = 3856
          NASASMAPvod_struc(n)%nr = 1624
       endif

       allocate(NASASMAPvod_struc(n)%vodobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(NASASMAPvod_struc(n)%vodtime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       
       NASASMAPvod_struc(n)%vodtime = -1
      
    enddo
    
    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       ssdev = obs_pert%ssdev(1)
       
       if(LIS_rc%dascaloption(k).eq."CDF matching") then 

          call LIS_getCDFattributes(k,modelcdffile(n),&
               NASASMAPvod_struc(n)%ntimes,ngrid)
          
          allocate(NASASMAPvod_struc(n)%model_xrange(&
               LIS_rc%obs_ngrid(k), NASASMAPvod_struc(n)%ntimes, &
               NASASMAPvod_struc(n)%nbins))
          allocate(NASASMAPvod_struc(n)%obs_xrange(&
               LIS_rc%obs_ngrid(k), NASASMAPvod_struc(n)%ntimes, &
               NASASMAPvod_struc(n)%nbins))
          allocate(NASASMAPvod_struc(n)%model_cdf(&
               LIS_rc%obs_ngrid(k), NASASMAPvod_struc(n)%ntimes, &
               NASASMAPvod_struc(n)%nbins))
          allocate(NASASMAPvod_struc(n)%obs_cdf(&
               LIS_rc%obs_ngrid(k), NASASMAPvod_struc(n)%ntimes, &
               NASASMAPvod_struc(n)%nbins))
          allocate(NASASMAPvod_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
               NASASMAPvod_struc(n)%ntimes))
          allocate(NASASMAPvod_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
               NASASMAPvod_struc(n)%ntimes))
          allocate(NASASMAPvod_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
               NASASMAPvod_struc(n)%ntimes))
          allocate(NASASMAPvod_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
               NASASMAPvod_struc(n)%ntimes))

!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
          call LIS_readMeanSigmaData(n,k,&
               NASASMAPvod_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               modelcdffile(n), &
               "LAI",&
               NASASMAPvod_struc(n)%model_mu,&
               NASASMAPvod_struc(n)%model_sigma)
          
          call LIS_readMeanSigmaData(n,k,&
               NASASMAPvod_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               obscdffile(n), &
               "VOD",&
               NASASMAPvod_struc(n)%obs_mu,&
               NASASMAPvod_struc(n)%obs_sigma)
          
          call LIS_readCDFdata(n,k,&
               NASASMAPvod_struc(n)%nbins,&
               NASASMAPvod_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               modelcdffile(n), &
               "LAI",&
               NASASMAPvod_struc(n)%model_xrange,&
               NASASMAPvod_struc(n)%model_cdf)
          
          call LIS_readCDFdata(n,k,&
               NASASMAPvod_struc(n)%nbins,&
               NASASMAPvod_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               obscdffile(n), &
               "VOD",&
               NASASMAPvod_struc(n)%obs_xrange,&
               NASASMAPvod_struc(n)%obs_cdf)
          
          
          if(NASASMAPvod_struc(n)%useSsdevScal.eq.1) then 
             if(NASASMAPvod_struc(n)%ntimes.eq.1) then 
                jj = 1
             else
                jj = LIS_rc%mo
             endif
             do t=1,LIS_rc%obs_ngrid(k)
                if(NASASMAPvod_struc(n)%obs_sigma(t,jj).gt.0) then 

                   print*, ssdev(t), NASASMAPvod_struc(n)%model_sigma(t,jj),&
                        NASASMAPvod_struc(n)%obs_sigma(t,jj)


                   ssdev(t) = ssdev(t)*NASASMAPvod_struc(n)%model_sigma(t,jj)/&
                        NASASMAPvod_struc(n)%obs_sigma(t,jj)
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

#if 0           
          allocate(obserr(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
          obserr = -9999.0

          do r=1,LIS_rc%obs_lnr(k)
             do c=1,LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                   obserr(c,r)  =  ssdev(LIS_obs_domain(n,k)%gindex(c,r)) 
                   
                endif
             enddo
          enddo
          print*, 'SMAP ',LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)
          open(100,file='smap_obs_err.bin',form='unformatted')
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
       endif
       if(LIS_rc%obs_ngrid(k).gt.0) then 
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
          call LIS_verify(status)
       endif

       deallocate(ssdev)
       
    enddo
    
    do n=1,LIS_rc%nnest
       if(NASASMAPvod_struc(n)%data_designation.eq."SPL3SMP") then 
          gridDesci = 0
          gridDesci(1) = 9
          gridDesci(2) = 964
          gridDesci(3) = 406
          gridDesci(9) = 4 !M36 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.36 
          gridDesci(11) = 1 !for the global switch
       elseif(NASASMAPvod_struc(n)%data_designation.eq."SPL3SMP_E") then 
          gridDesci = 0
          gridDesci(1) = 9
          gridDesci(2) = 3856
          gridDesci(3) = 1624
          gridDesci(9) = 5 !M09 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.09 
          gridDesci(11) = 1 !for the global switch
       endif

       allocate(NASASMAPvod_struc(n)%n11(LIS_rc%obs_lnc(k)*&
            LIS_rc%obs_lnr(k)))
       allocate(NASASMAPvod_struc(n)%rlat(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(NASASMAPvod_struc(n)%rlon(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       
       call neighbor_interp_input_withgrid(gridDesci,&
            LIS_rc%obs_gridDesc(k,:),&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
            NASASMAPvod_struc(n)%rlat, &
            NASASMAPvod_struc(n)%rlon, &
            NASASMAPvod_struc(n)%n11)
       
       call LIS_registerAlarm("NASASMAP read alarm",&
            86400.0, 86400.0)

       NASASMAPvod_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine NASASMAPvod_setup
end module NASASMAPvod_Mod
