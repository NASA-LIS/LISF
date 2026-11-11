!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: WSFsm_Mod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines to
!   handle WSF soil moisture retrievals
!
! !REVISION HISTORY:
!  2 Oct 2025 Ehsan Jalilvand; initial specification
!  Oct 2026   Ehsan Jalilvand; grid of the WSF files read from lis.config
!             ("WSF soil moisture grid ..." entries). All six entries
!             are required; there is no default grid.

module WSFsm_Mod
! !USES:
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: WSFsm_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: WSFsm_struc ! make it accessible by other module
!EOP
  type, public:: WSFsm_dec

     integer                :: useSsdevScal
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: smobs(:,:)
     real*8,     allocatable    :: smtime(:,:)

     real                       :: ssdev_inp
     real, allocatable          :: rlat(:)
     real, allocatable          :: rlon(:)
     real                       :: gridDesci(50)
     real                       :: lat_ll, lon_ll  ! lower-left cell centre
     real                       :: dlat, dlon      ! grid spacing (deg)
     logical                    :: use_bilinear
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

     integer                    :: nbins
     integer                    :: ntimes
     logical                    :: cdf_read_mon  !(for reading monthly CDF when
                                                 !LIS_rc%da > 1 but the first model time step,
                                                 !e.g., 4/29 13:00:00)
     integer                    :: cdf_read_opt  ! 0: read all months at one time
                                                 ! 1: read only the current month
     character(len=LIS_CONST_PATH_LEN) :: modelcdffile
     character(len=LIS_CONST_PATH_LEN) :: obscdffile

  end type WSFsm_dec

  type(WSFsm_dec),allocatable :: WSFsm_struc(:)

contains

!BOP
!
! !ROUTINE: WSFsm_setup
! \label{WSFsm_setup}
!
! !INTERFACE:
  subroutine WSFsm_setup(k, OBS_State, OBS_Pert_State)
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
!   creation of data strctures required for handling WSF
!   soil moisture data.
!
!   The arguments are:
!   \begin{description}
!    \item[OBS\_State]   observation state
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    real, parameter                   ::  minssdev =0.001
    integer                           ::  n,i,t,jj
    integer                           ::  ftn
    integer                           ::  status
    type(ESMF_Field)                  ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)              ::  intarrspec, realarrspec
    type(ESMF_Field)                  ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)              ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  smobsdir
    character*100                     ::  temp
    character*1                       ::  vid(2)
    character*40, allocatable         ::  vname(:)
    real        , allocatable         ::  varmin(:)
    real        , allocatable         ::  varmax(:)
    type(pert_dec_type)               ::  obs_pert
    real, pointer                     ::  obs_temp(:,:)
    real, allocatable                 ::  ssdev(:)
    integer                           ::  ngrid
    real, allocatable                 ::  grid_dd(:,:)
    character(len=48), parameter      ::  grid_labels(6) = [character(len=48) :: &
         'WSF soil moisture grid lower left lat:', &
         'WSF soil moisture grid lower left lon:', &
         'WSF soil moisture grid upper right lat:', &
         'WSF soil moisture grid upper right lon:', &
         'WSF soil moisture grid resolution (dx):', &
         'WSF soil moisture grid resolution (dy):' ]

    external :: bilinear_interp_input_withgrid
    external :: upscaleByAveraging_input

    allocate(WSFsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"WSF soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,smobsdir,&
            rc=status)
       call LIS_verify(status, 'WSF soil moisture data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            smobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WSF soil moisture use scaled standard deviation model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,WSFsm_struc(n)%useSsdevScal,&
            rc=status)
       call LIS_verify(status, 'WSF soil moisture use scaled standard deviation model: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WSF model CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,WSFsm_struc(n)%modelcdffile,rc=status)
       call LIS_verify(status, 'WSF model CDF file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WSF observation CDF file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,WSFsm_struc(n)%obscdffile,rc=status)
       call LIS_verify(status, 'WSF observation CDF file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config, "WSF soil moisture number of bins in the CDF:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,WSFsm_struc(n)%nbins, rc=status)
       call LIS_verify(status, "WSF soil moisture number of bins in the CDF: not defined")
    enddo

   do n=1, LIS_rc%nnest
      WSFsm_struc(n)%cdf_read_mon = .false.

      call ESMF_ConfigFindLabel(LIS_config, "WSF CDF read option:", rc=status)    ! 0: read CDF for all months/year
                                                                                         ! 1: read CDF for current month
      call ESMF_ConfigGetAttribute(LIS_config, WSFsm_struc(n)%cdf_read_opt, rc=status)
      call LIS_verify(status, "WSF CDF read option: not defined")
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
         '[INFO] read WSF soil moisture data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. WSF
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
          WSFsm_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
         '[INFO] Created the States to hold the WSF observations data'

    ! Grid of the WSF soil moisture files (cell centres). All six
    ! entries are required; there is no default grid.
    allocate(grid_dd(LIS_rc%nnest,6))
    do i = 1, 6
       call ESMF_ConfigFindLabel(LIS_config, trim(grid_labels(i)), rc=status)
       if (status /= 0) then
          write(LIS_logunit,*) '[ERR] ', trim(grid_labels(i)), &
               ' not defined in lis.config'
          write(LIS_logunit,*) '[ERR] Give the grid of the WSF soil ', &
               'moisture files with all six entries:'
          write(LIS_logunit,*) '[ERR]   ', trim(grid_labels(1)), ' ... ', &
               trim(grid_labels(6))
       endif
       call LIS_verify(status, trim(grid_labels(i))//' not defined')
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, grid_dd(n,i), rc=status)
          call LIS_verify(status, trim(grid_labels(i))//' not defined')
       enddo
    enddo
    do n=1,LIS_rc%nnest
       if (grid_dd(n,5) <= 0.0 .or. grid_dd(n,6) <= 0.0 .or. &
            grid_dd(n,3) <= grid_dd(n,1) .or. &
            grid_dd(n,4) <= grid_dd(n,2)) then
          write(LIS_logunit,*) '[ERR] Invalid WSF soil moisture grid in ', &
               'lis.config (lower left, upper right, dx, dy):'
          write(LIS_logunit,*) '[ERR] ', grid_dd(n,:)
          call LIS_endrun()
       endif
    enddo

    do n=1,LIS_rc%nnest
       WSFsm_struc(n)%lat_ll = grid_dd(n,1)
       WSFsm_struc(n)%lon_ll = grid_dd(n,2)
       WSFsm_struc(n)%dlon   = grid_dd(n,5)
       WSFsm_struc(n)%dlat   = grid_dd(n,6)
       WSFsm_struc(n)%nc = nint((grid_dd(n,4) - grid_dd(n,2))/grid_dd(n,5)) + 1
       WSFsm_struc(n)%nr = nint((grid_dd(n,3) - grid_dd(n,1))/grid_dd(n,6)) + 1

       WSFsm_struc(n)%gridDesci = 0
       WSFsm_struc(n)%gridDesci(1) = 0
       WSFsm_struc(n)%gridDesci(2) = WSFsm_struc(n)%nc
       WSFsm_struc(n)%gridDesci(3) = WSFsm_struc(n)%nr
       WSFsm_struc(n)%gridDesci(4) = grid_dd(n,1)
       WSFsm_struc(n)%gridDesci(5) = grid_dd(n,2)
       WSFsm_struc(n)%gridDesci(6) = 128
       WSFsm_struc(n)%gridDesci(7) = grid_dd(n,3)
       WSFsm_struc(n)%gridDesci(8) = grid_dd(n,4)
       WSFsm_struc(n)%gridDesci(9) = grid_dd(n,5)   !dlon
       WSFsm_struc(n)%gridDesci(10) = grid_dd(n,6)  !dlat
       WSFsm_struc(n)%gridDesci(20) = 64

       ! Bilinear if the obs domain is finer than the WSF grid, else
       ! averaging (identity when the two grids are the same). Same rule
       ! as the LDT reader used for the CDF.
       WSFsm_struc(n)%use_bilinear = &
            LIS_rc%obs_gridDesc(k,10) < WSFsm_struc(n)%dlat

       write(LIS_logunit,'(a,i6,a,i6,a,2f12.6,a,2f10.6)') &
            ' [INFO] WSF SM grid: ', WSFsm_struc(n)%nc, ' x ', WSFsm_struc(n)%nr, &
            ', lower left ', grid_dd(n,1), grid_dd(n,2), &
            ', dx dy ', grid_dd(n,5), grid_dd(n,6)

       allocate(WSFsm_struc(n)%smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(WSFsm_struc(n)%smtime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))

       WSFsm_struc(n)%smtime = -1

    enddo
    deallocate(grid_dd)

    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       ssdev = obs_pert%ssdev(1)

       if(LIS_rc%dascaloption(k).eq."CDF matching" .or. &
               LIS_rc%dascaloption(k).eq."Linear scaling" .or. &
               LIS_rc%dascaloption(k).eq."Anomaly scaling") then

          call LIS_getCDFattributes(k,WSFsm_struc(n)%modelcdffile,&
               WSFsm_struc(n)%ntimes,ngrid)

          if (WSFsm_struc(n)%cdf_read_opt.eq.0) then
             allocate(WSFsm_struc(n)%model_xrange(&
                  LIS_rc%obs_ngrid(k), WSFsm_struc(n)%ntimes, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%obs_xrange(&
                  LIS_rc%obs_ngrid(k), WSFsm_struc(n)%ntimes, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%model_cdf(&
                  LIS_rc%obs_ngrid(k), WSFsm_struc(n)%ntimes, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%obs_cdf(&
                  LIS_rc%obs_ngrid(k), WSFsm_struc(n)%ntimes, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
                  WSFsm_struc(n)%ntimes))
             allocate(WSFsm_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
                  WSFsm_struc(n)%ntimes))
             allocate(WSFsm_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
                  WSFsm_struc(n)%ntimes))
             allocate(WSFsm_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
                  WSFsm_struc(n)%ntimes))
          else
             allocate(WSFsm_struc(n)%model_xrange(&
                  LIS_rc%obs_ngrid(k), 1, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%obs_xrange(&
                  LIS_rc%obs_ngrid(k), 1, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%model_cdf(&
                  LIS_rc%obs_ngrid(k), 1, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%obs_cdf(&
                  LIS_rc%obs_ngrid(k), 1, &
                  WSFsm_struc(n)%nbins))
             allocate(WSFsm_struc(n)%model_mu(LIS_rc%obs_ngrid(k),1))
             allocate(WSFsm_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),1))
             allocate(WSFsm_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),1))
             allocate(WSFsm_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),1))
          endif
!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
         if (WSFsm_struc(n)%cdf_read_opt.eq.0) then
          call LIS_readMeanSigmaData(n,k,&
               WSFsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               WSFsm_struc(n)%modelcdffile, &
               "SoilMoist",&
               WSFsm_struc(n)%model_mu,&
               WSFsm_struc(n)%model_sigma)

          call LIS_readMeanSigmaData(n,k,&
               WSFsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               WSFsm_struc(n)%obscdffile, &
               "SoilMoist",&
               WSFsm_struc(n)%obs_mu,&
               WSFsm_struc(n)%obs_sigma)

          call LIS_readCDFdata(n,k,&
               WSFsm_struc(n)%nbins,&
               WSFsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               WSFsm_struc(n)%modelcdffile, &
               "SoilMoist",&
               WSFsm_struc(n)%model_xrange,&
               WSFsm_struc(n)%model_cdf)

          call LIS_readCDFdata(n,k,&
               WSFsm_struc(n)%nbins,&
               WSFsm_struc(n)%ntimes,&
               LIS_rc%obs_ngrid(k), &
               WSFsm_struc(n)%obscdffile, &
               "SoilMoist",&
               WSFsm_struc(n)%obs_xrange,&
               WSFsm_struc(n)%obs_cdf)

          if(WSFsm_struc(n)%useSsdevScal.eq.1) then
             if(WSFsm_struc(n)%ntimes.eq.1) then
                jj = 1
             else
                jj = LIS_rc%mo
             endif
             do t=1,LIS_rc%obs_ngrid(k)
                if(WSFsm_struc(n)%obs_sigma(t,jj).gt.0) then

                   print*, ssdev(t), WSFsm_struc(n)%model_sigma(t,jj),&
                        WSFsm_struc(n)%obs_sigma(t,jj)


                   ssdev(t) = ssdev(t)*WSFsm_struc(n)%model_sigma(t,jj)/&
                        WSFsm_struc(n)%obs_sigma(t,jj)
                   !                c = LIS_domain(n)%grid(t)%col
                   !                r = LIS_domain(n)%grid(t)%row
                   !                ssdev_grid(c,r) = ssdev(t)
                   if(ssdev(t).lt.minssdev) then
                      ssdev(t) = minssdev
                   endif
                endif
             enddo
          endif
         endif
       endif

       if(LIS_rc%obs_ngrid(k).gt.0) then
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
          call LIS_verify(status)
       endif

       deallocate(ssdev)

    enddo

    do n=1,LIS_rc%nnest

       if(WSFsm_struc(n)%use_bilinear) then

          allocate(WSFsm_struc(n)%rlat(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%rlon(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%n11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%n12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%n21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%n22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%w11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%w12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%w21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WSFsm_struc(n)%w22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

          call bilinear_interp_input_withgrid(WSFsm_struc(n)%gridDesci(:), &
               LIS_rc%obs_gridDesc(k,:),&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
               WSFsm_struc(n)%rlat, WSFsm_struc(n)%rlon,&
               WSFsm_struc(n)%n11, WSFsm_struc(n)%n12, &
               WSFsm_struc(n)%n21, WSFsm_struc(n)%n22, &
               WSFsm_struc(n)%w11, WSFsm_struc(n)%w12, &
               WSFsm_struc(n)%w21, WSFsm_struc(n)%w22)

       else

          allocate(WSFsm_struc(n)%n11(&
               WSFsm_struc(n)%nc*WSFsm_struc(n)%nr))

          call upscaleByAveraging_input(WSFsm_struc(n)%gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),&
               WSFsm_struc(n)%nc*WSFsm_struc(n)%nr, &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), WSFsm_struc(n)%n11)
       endif

       call LIS_registerAlarm("WSF read alarm",&
            3600.0, 3600.0)

       WSFsm_struc(n)%startMode = .true.

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo
  end subroutine WSFsm_setup
end module WSFsm_Mod
