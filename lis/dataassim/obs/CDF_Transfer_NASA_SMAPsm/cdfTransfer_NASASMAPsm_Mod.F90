!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE:  cdfTransfer_NASASMAPsm_Mod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines to
!   handle
!
! !REVISION HISTORY:
!  2 Mar 2022    Mahdi Navari; initial specification
!
module cdfTransfer_NASASMAPsm_Mod
! !USES:
  use ESMF
  use map_utils

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: cdfTransfer_NASASMAPsm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: cdfT_SMAPsm_struc
!EOP
  type, public:: cdfT_SMAPsm_dec

     integer                :: useSsdevScal
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: smobs(:,:)
     real*8,     allocatable    :: smtime(:,:)

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
     character*20               :: data_designation
     character*3                :: release_number
     integer                    :: nbins
     integer                    :: ntimes

     logical                    :: cdf_read_mon  !(for reading monthly CDF when
                                                 !LIS_rc%da > 1 but the first model time step,
                                                 !e.g., 4/29 13:00:00)
     integer                    :: cdf_read_opt  ! 0: read all months at one time
                                                 ! 1: read only the current month
     character*100              :: modelcdffile
     character*100              :: obscdffile
     integer                    :: n_strat_bins
     integer                    :: useCDFtransfer
     character*100              :: ref_p_climo_file
     character*100              :: target_p_climo_file
     real,    allocatable       :: ref_p_climo_maxval(:)
     real,    allocatable       :: target_p_climo(:,:,:)

  end type cdfT_SMAPsm_dec

  type(cdfT_SMAPsm_dec),allocatable :: cdfT_SMAPsm_struc(:)

contains

!BOP
!
! !ROUTINE: cdfTransfer_NASASMAPsm_setup
! \label{cdfTransfer_NASASMAPsm_setup}
!
! !INTERFACE:
  subroutine cdfTransfer_NASASMAPsm_setup(k, OBS_State, OBS_Pert_State)
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
    real, parameter        ::  minssdev =0.001
    integer                ::  n,i,t,kk,jj
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character*100          ::  rtsmopssmobsdir
    character*100          ::  temp
    real,  allocatable     ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real                   :: gridDesci(50)
    real, allocatable      :: ssdev(:)

    real, allocatable      ::  obserr(:,:)
    real, allocatable      ::  lobserr(:,:)
    integer                :: c,r
    real, allocatable      :: ssdev_grid(:,:)
    integer                :: ngrid
    real, allocatable      :: target_precip_climo(:,:,:)

    allocate(cdfT_SMAPsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NASA) soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,rtsmopssmobsdir,&
            rc=status)
       call LIS_verify(status, 'SMAP(NASA) soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            rtsmopssmobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NASA) soil moisture data designation:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            cdfT_SMAPsm_struc(n)%data_designation,&
            rc=status)
       call LIS_verify(status, 'SMAP(NASA) soil moisture data designation: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NASA) soil moisture Composite Release ID:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            cdfT_SMAPsm_struc(n)%release_number,&
            rc=status)
       call LIS_verify(status, 'SMAP(NASA) soil moisture Composite Release ID: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SMAP(NASA) soil moisture use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,cdfT_SMAPsm_struc(n)%useSsdevScal,&
            rc=status)
       call LIS_verify(status, 'SMAP(NASA) soil moisture use scaled standard deviation model: is missing')

    enddo

    call ESMF_ConfigFindLabel(LIS_config, "SMAP(NASA) soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,cdfT_SMAPsm_struc(n)%nbins, rc=status)
       call LIS_verify(status, "SMAP(NASA) soil moisture number of bins in the CDF: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Use CDF transfer for soil moisture data assimilation:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            cdfT_SMAPsm_struc(n)%useCDFtransfer,&
            default=0,rc=status)
       call LIS_verify(status, 'Use CDF transfer for soil moisture data assimilation: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Reference domain model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            cdfT_SMAPsm_struc(n)%modelcdffile,&
            rc=status)
       call LIS_verify(status, 'Reference domain model CDF file: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Reference domain obs CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            cdfT_SMAPsm_struc(n)%obscdffile,&
            rc=status)
       call LIS_verify(status, 'Reference domain obs CDF file: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Reference domain precipitation climatology data source:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            cdfT_SMAPsm_struc(n)%ref_p_climo_file,&
            rc=status)
       call LIS_verify(status, 'Reference domain precipitation climatology data source: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Target domain precipitation climatology data source:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            cdfT_SMAPsm_struc(n)%target_p_climo_file,&
            rc=status)
       call LIS_verify(status, 'Target domain precipitation climatology data source: is missing')
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
         '[INFO] read reference and target domains precipitation climatology data'
!--------------------------------------------------------------------------------
! Data must be in the local domain
! The stratified CDF is geolocation independent
! Subset the precipition climotology to local space
!-------------------------------------------------------------------------------
    do n=1,LIS_rc%nnest
       allocate(cdfT_SMAPsm_struc(n)%ref_p_climo_maxval(12)) !cdfT_SMAPsm_struc(n)%ntimes))
       allocate(cdfT_SMAPsm_struc(n)%target_p_climo(LIS_rc%lnc(n), LIS_rc%lnr(n),12)) !cdfT_SMAPsm_struc(n)%ntimes))
       ! Note: we do not need to know the dimenesions of the reference domain precip climatology
       !       for allocation so we just extract max value for each month

       call read_Precip_climo_maxval (n,cdfT_SMAPsm_struc(n)%ref_p_climo_file,&
                               cdfT_SMAPsm_struc(n)%ref_p_climo_maxval)

       call read_Precip_climo (n,k,&
                               cdfT_SMAPsm_struc(n)%target_p_climo_file,&
                               target_precip_climo)

      ! subset from the global 2-d space to the local 2-d space
        do i=1,12
           cdfT_SMAPsm_struc(n)%target_p_climo(:,:,i) = target_precip_climo(&
                    LIS_ews_halo_ind(n,LIS_localPet+1):&
                    LIS_ewe_halo_ind(n,LIS_localPet+1), &
                    LIS_nss_halo_ind(n,LIS_localPet+1): &
                    LIS_nse_halo_ind(n,LIS_localPet+1),i)
        enddo
    enddo

    write(LIS_logunit,*)&
         '[INFO] read SMAP(NASA) soil moisture data specifications'

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
          cdfT_SMAPsm_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
       if(cdfT_SMAPsm_struc(n)%data_designation.eq."SPL3SMP") then
          cdfT_SMAPsm_struc(n)%nc = 964
          cdfT_SMAPsm_struc(n)%nr = 406
       elseif(cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP") then
          cdfT_SMAPsm_struc(n)%nc = 964
          cdfT_SMAPsm_struc(n)%nr = 406
       elseif(cdfT_SMAPsm_struc(n)%data_designation.eq."SPL3SMP_E") then
          cdfT_SMAPsm_struc(n)%nc = 3856
          cdfT_SMAPsm_struc(n)%nr = 1624
       elseif(cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP_E") then
          cdfT_SMAPsm_struc(n)%nc = 3856
          cdfT_SMAPsm_struc(n)%nr = 1624
       endif

       allocate(cdfT_SMAPsm_struc(n)%smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(cdfT_SMAPsm_struc(n)%smtime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))

       cdfT_SMAPsm_struc(n)%smtime = -1

    enddo

    write(LIS_logunit,*)&
         '[INFO] read cdf transfer data specifications'

    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       ssdev = obs_pert%ssdev(1)

       if(cdfT_SMAPsm_struc(n)%useCDFtransfer.gt.0) then ! LIS_rc%dascaloption(k).eq."CDF matching" .and. &

          call LIS_getCDFtransferattributes(k,cdfT_SMAPsm_struc(n)%modelcdffile,&
                   cdfT_SMAPsm_struc(n)%n_strat_bins , cdfT_SMAPsm_struc(n)%ntimes)

          allocate(cdfT_SMAPsm_struc(n)%model_xrange(&
               cdfT_SMAPsm_struc(n)%n_strat_bins, cdfT_SMAPsm_struc(n)%ntimes, &
               cdfT_SMAPsm_struc(n)%nbins))
          allocate(cdfT_SMAPsm_struc(n)%obs_xrange(&
               cdfT_SMAPsm_struc(n)%n_strat_bins, cdfT_SMAPsm_struc(n)%ntimes, &
               cdfT_SMAPsm_struc(n)%nbins))
          allocate(cdfT_SMAPsm_struc(n)%model_cdf(&
               cdfT_SMAPsm_struc(n)%n_strat_bins, cdfT_SMAPsm_struc(n)%ntimes, &
               cdfT_SMAPsm_struc(n)%nbins))
          allocate(cdfT_SMAPsm_struc(n)%obs_cdf(&
               cdfT_SMAPsm_struc(n)%n_strat_bins, cdfT_SMAPsm_struc(n)%ntimes, &
               cdfT_SMAPsm_struc(n)%nbins))
          write(LIS_logunit,*)&
               '[INFO] Successfully read cdf transfer data specifications'
       endif

       if(cdfT_SMAPsm_struc(n)%useSsdevScal.eq.1) then
          write(LIS_logunit,*) '[ERR] "use scaled standard deviation model" does not work '
          write(LIS_logunit,*) '[ERR] for CDF transfer method'
          call LIS_endrun()
       endif

       if(LIS_rc%obs_ngrid(k).gt.0) then
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
          call LIS_verify(status)
       endif

       deallocate(ssdev)
    enddo


    do n=1,LIS_rc%nnest
       if(cdfT_SMAPsm_struc(n)%data_designation.eq."SPL3SMP".or.&
            cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP") then
          gridDesci = 0
          gridDesci(1) = 9
          gridDesci(2) = 964
          gridDesci(3) = 406
          gridDesci(9) = 4 !M36 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.36
          gridDesci(11) = 1 !for the global switch
       elseif(cdfT_SMAPsm_struc(n)%data_designation.eq."SPL3SMP_E".or.&
            cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP_E") then
          gridDesci = 0
          gridDesci(1) = 9
          gridDesci(2) = 3856
          gridDesci(3) = 1624
          gridDesci(9) = 5 !M09 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.09
          gridDesci(11) = 1 !for the global switch
       endif

       allocate(cdfT_SMAPsm_struc(n)%n11(LIS_rc%obs_lnc(k)*&
            LIS_rc%obs_lnr(k)))
       allocate(cdfT_SMAPsm_struc(n)%rlat(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(cdfT_SMAPsm_struc(n)%rlon(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

       call neighbor_interp_input_withgrid(gridDesci,&
            LIS_rc%obs_gridDesc(k,:),&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
            cdfT_SMAPsm_struc(n)%rlat, &
            cdfT_SMAPsm_struc(n)%rlon, &
            cdfT_SMAPsm_struc(n)%n11)

       if(cdfT_SMAPsm_struc(n)%data_designation.eq."SPL3SMP".or.&
            cdfT_SMAPsm_struc(n)%data_designation.eq."SPL3SMP_E") then
          call LIS_registerAlarm("NASASMAP read alarm",&
               86400.0, 86400.0)
       else
          call LIS_registerAlarm("NASASMAP read alarm",&
               3600.0, 3600.0)
       endif

       cdfT_SMAPsm_struc(n)%startMode = .true.

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine cdfTransfer_NASASMAPsm_setup

end module cdfTransfer_NASASMAPsm_Mod
