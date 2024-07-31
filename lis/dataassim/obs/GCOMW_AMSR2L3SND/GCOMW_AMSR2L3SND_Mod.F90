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
! !MODULE: GCOMW_AMSR2L3SND_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle AMSR2 retrievals
! 
! !REVISION HISTORY: 
!  12 Jan 15    Sujay Kumar; Initial version
! 
module GCOMW_AMSR2L3SND_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: GCOMW_AMSR2L3SND_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: GCOMW_AMSR2L3SND_struc
!EOP
  type, public:: GCOMW_AMSR2L3SND_dec
     
     integer                :: bc_version
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: sndobs(:,:)
     real,     allocatable      :: sndtime(:,:)

     real                   :: datares
     real                   :: ssdev_inp
     integer                :: amsr2nc, amsr2nr
     type(proj_info)        :: amsr2proj
     integer, allocatable       :: n11(:)
     integer, allocatable       :: n12(:)
     integer, allocatable       :: n21(:)
     integer, allocatable       :: n22(:)
     real, allocatable          :: w11(:)
     real, allocatable          :: w12(:)
     real, allocatable          :: w21(:)
     real, allocatable          :: w22(:)
     real,    allocatable       :: rlat(:)
     real,    allocatable       :: rlon(:)

     integer             :: useIMS
     character(len=LIS_CONST_PATH_LEN) :: IMSdir
     integer             :: useMODIS
     character(len=LIS_CONST_PATH_LEN) :: MODISdir

     integer             :: ims_mi
     real                :: ims_gridDesci(50)
     real, allocatable    :: ims_rlat(:) ! EMK...Corrected type
     real, allocatable    :: ims_rlon(:) ! EMK...Corrected type
     integer, allocatable    :: ims_n11(:)
     integer, allocatable    :: ims_n12(:)
     integer, allocatable    :: ims_n21(:)
     integer, allocatable    :: ims_n22(:)
     real, allocatable    :: ims_w11(:) ! EMK...Corrected type
     real, allocatable    :: ims_w12(:) ! EMK...Corrected type
     real, allocatable    :: ims_w21(:) ! EMK...Corrected type
     real, allocatable    :: ims_w22(:) ! EMK...Corrected type
     
     real, allocatable       :: input_mask(:,:)
     integer             :: usr_input_mask
     integer             :: mod_mi
     real                :: mod_gridDesci(50)
     integer, allocatable    :: mod_n11(:)
     integer, allocatable    :: mod_n12(:)
     integer, allocatable    :: mod_n21(:)
     integer, allocatable    :: mod_n22(:)
     real, allocatable    :: mod_w11(:)
     real, allocatable    :: mod_w12(:)
     real, allocatable    :: mod_w21(:)
     real, allocatable    :: mod_w22(:)
     real, allocatable    :: mod_rlat(:)
     real, allocatable    :: mod_rlon(:)

     real, allocatable   :: IMSdata_obs(:)
     real, allocatable   :: MODISdata_obs(:)

  end type GCOMW_AMSR2L3SND_dec
  
  type(GCOMW_AMSR2L3SND_dec),allocatable :: GCOMW_AMSR2L3SND_struc(:)
  
contains

!BOP
! 
! !ROUTINE: GCOMW_AMSR2L3SND_setup
! \label{GCOMW_AMSR2L3SND_setup}
! 
! !INTERFACE: 
  subroutine GCOMW_AMSR2L3SND_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_dataAssimMod
    use LIS_perturbMod
    use LIS_logmod
    use LIS_DAobservationsMod

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

    integer                ::  n,i,t,kk,c,r
    real, allocatable          ::  obserr(:,:)
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  amsresndobsdir
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
    integer                :: jj
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2  
    real, allocatable      :: obsmask(:,:)
    character(len=LIS_CONST_PATH_LEN) :: input_mask_file
    integer                ::  modis_nc,modis_nr
    integer                ::  ims_nc, ims_nr


    allocate(GCOMW_AMSR2L3SND_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) snow depth data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,amsresndobsdir,&
            rc=status)
       call LIS_verify(status, 'AMSR2(GCOMW) snow depth data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            amsresndobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "AMSR2(GCOMW) snow depth use IMS data for snow detection:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,GCOMW_AMSR2L3SND_struc(n)%useIMS,&
            rc=status)
       call LIS_verify(status,&
            'AMSR2(GCOMW) snow depth use IMS data for snow detection: not defined')
       if(GCOMW_AMSR2L3SND_struc(n)%useIMS.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,&
               GCOMW_AMSR2L3SND_struc(n)%IMSdir,&
               label="AMSR2(GCOMW) snow depth IMS data directory:",rc=status)
          call LIS_verify(status,&
               'AMSR2(GCOMW) snow depth IMS data directory: option not specified')
       endif
    enddo
    

    call ESMF_ConfigFindLabel(LIS_config,&
         "AMSR2(GCOMW) snow depth use MODIS (MOD10C1) data for snow detection:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,GCOMW_AMSR2L3SND_struc(n)%useMODIS,&
            rc=status)
       call LIS_verify(status,&
            'AMSR2(GCOMW) snow depth use MODIS (MOD10C1) data for snow detection: not defined')
       if(GCOMW_AMSR2L3SND_struc(n)%useMODIS.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,&
               GCOMW_AMSR2L3SND_struc(n)%MODISdir,&
               label="AMSR2(GCOMW) snow depth MOD10C1 data directory:",rc=status)
          call LIS_verify(status,&
               'AMSR2(GCOMW) snow depth MOD10C1 data directory: option not specified')
       endif

    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) snow depth use bias corrected version:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GCOMW_AMSR2L3SND_struc(n)%bc_version,&
            rc=status)
       call LIS_verify(status, 'AMSR2(GCOMW) snow depth use bias corrected version: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"AMSR2(GCOMW) snow depth use input mask:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GCOMW_AMSR2L3SND_struc(n)%usr_input_mask,&
            default = 0, rc=status)
       call LIS_verify(status, 'AMSR2(GCOMW) snow depth use input mask: is missing')
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

    write(LIS_logunit,*)'[INFO] read AMSR2(GCOMW) snow depth data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. amsr-e 
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
          GCOMW_AMSR2L3SND_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
         '[INFO] Created the States to hold the AMSR2(GCOMW) observations data'
    do n=1,LIS_rc%nnest
       GCOMW_AMSR2L3SND_struc(n)%nc = 3600
       GCOMW_AMSR2L3SND_struc(n)%nr = 1800
       allocate(GCOMW_AMSR2L3SND_struc(n)%sndobs(&
            LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%sndtime(&
            LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))

       GCOMW_AMSR2L3SND_struc(n)%sndobs = LIS_rc%udef
       GCOMW_AMSR2L3SND_struc(n)%sndtime = -1
    enddo
    
    do n=1,LIS_rc%nnest
       GCOMW_AMSR2L3SND_struc(n)%amsr2nc = 3600
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr = 1800

       call map_set(PROJ_LATLON, -89.95,-179.95,&
            0.0, 0.10,0.10, 0.0,&
            GCOMW_AMSR2L3SND_struc(n)%amsr2nc,&
            GCOMW_AMSR2L3SND_struc(n)%amsr2nr,&
            GCOMW_AMSR2L3SND_struc(n)%amsr2proj)
       
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
       
       GCOMW_AMSR2L3SND_struc(n)%datares = 0.10
       
       allocate(GCOMW_AMSR2L3SND_struc(n)%n11(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%n12(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%n21(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%n22(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%w11(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%w12(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%w21(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%w22(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%rlat(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%rlon(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          
       call bilinear_interp_input_withgrid(gridDesci(:), & 
            LIS_rc%obs_gridDesc(k,:),&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
            GCOMW_AMSR2L3SND_struc(n)%rlat, GCOMW_AMSR2L3SND_struc(n)%rlon,&
            GCOMW_AMSR2L3SND_struc(n)%n11,&
            GCOMW_AMSR2L3SND_struc(n)%n12,&
            GCOMW_AMSR2L3SND_struc(n)%n21,&
            GCOMW_AMSR2L3SND_struc(n)%n22,&
            GCOMW_AMSR2L3SND_struc(n)%w11,&
            GCOMW_AMSR2L3SND_struc(n)%w12,&
            GCOMW_AMSR2L3SND_struc(n)%w21,&
            GCOMW_AMSR2L3SND_struc(n)%w22)

       call LIS_registerAlarm("AMSR2(GCOMW) read alarm",&
            86400.0, 86400.0)
       GCOMW_AMSR2L3SND_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo

    do n=1,LIS_rc%nnest
       GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci = 0
       if(GCOMW_AMSR2L3SND_struc(n)%useIMS.eq.1) then 
          
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(1) = 0 
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(2) = 1500
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(3) = 375
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(4) = 0.12
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(5) = -179.88
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(6) = 128
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(7) = 89.88
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(8) = 179.88
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(9) = 0.24
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(10) = 0.24
          GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(20) = 64.0
          
          ims_nc = 1500
          ims_nr = 375
          GCOMW_AMSR2L3SND_struc(n)%ims_mi = ims_nc*ims_nr
!-----------------------------------------------------------------------------
!   Set up interpolation weights to transform the IMS data to the AMSR2 grid
!-----------------------------------------------------------------------------

          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_rlat(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_rlon(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_n11(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_n12(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_n21(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_n22(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_w11(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_w12(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_w21(ims_nc*ims_nr))
          allocate(GCOMW_AMSR2L3SND_struc(n)%ims_w22(ims_nc*ims_nr))          

          call bilinear_interp_input_withgrid(&
               GCOMW_AMSR2L3SND_struc(n)%ims_gridDesci(:), & 
               LIS_rc%obs_gridDesc(k,:),&
               ims_nc*ims_nr,&
               GCOMW_AMSR2L3SND_struc(n)%ims_rlat, &
               GCOMW_AMSR2L3SND_struc(n)%ims_rlon, &
               GCOMW_AMSR2L3SND_struc(n)%ims_n11,&
               GCOMW_AMSR2L3SND_struc(n)%ims_n12,&
               GCOMW_AMSR2L3SND_struc(n)%ims_n21,&
               GCOMW_AMSR2L3SND_struc(n)%ims_n22,&
               GCOMW_AMSR2L3SND_struc(n)%ims_w11,&
               GCOMW_AMSR2L3SND_struc(n)%ims_w12,&
               GCOMW_AMSR2L3SND_struc(n)%ims_w21,&
               GCOMW_AMSR2L3SND_struc(n)%ims_w22)

       endif
    enddo

    do n=1,LIS_rc%nnest
       GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci = 0

       if(GCOMW_AMSR2L3SND_struc(n)%useMODIS.eq.1) then 
          
          cornerlat1 = max(-59.975, &
               nint((LIS_rc%gridDesc(n,4)+59.975)/0.05)*0.05-59.975-2*0.05)
          cornerlon1 = max(-179.975, &
               nint((LIS_rc%gridDesc(n,5)+179.975)/0.05)*0.05-179.975-2*0.05)
          cornerlat2 = min(89.975, &
               nint((LIS_rc%gridDesc(n,7)+59.975)/0.05)*0.05-59.975+2*0.05)
          cornerlon2 = min(179.975, &
               nint((LIS_rc%gridDesc(n,8)+179.975)/0.05)*0.05-179.975+2*0.05)
          
          modis_nr = nint((cornerlat2-cornerlat1)/0.05)+1
          modis_nc = nint((cornerlon2-cornerlon1)/0.05)+1
          
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(1) = 0 
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(2) = modis_nc
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(3) = modis_nr
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(4) = cornerlat1
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(5) = cornerlon1
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(6) = 128
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(7) = cornerlat2
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(8) = cornerlon2
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(9) = 0.05
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(10) = 0.05
          GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(20) = 64.0
          
          GCOMW_AMSR2L3SND_struc(n)%mod_mi = modis_nc*modis_nr
!-----------------------------------------------------------------------------
!   Use upscaling since SSMI data is coarser than MODIS
!-----------------------------------------------------------------------------
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_rlat(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_rlon(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_n11(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_n12(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_n21(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_n22(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_w11(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_w12(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_w21(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%mod_w22(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

          call bilinear_interp_input_withgrid(&
               GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
               GCOMW_AMSR2L3SND_struc(n)%mod_rlat, &
               GCOMW_AMSR2L3SND_struc(n)%mod_rlon,&
               GCOMW_AMSR2L3SND_struc(n)%mod_n11,&
               GCOMW_AMSR2L3SND_struc(n)%mod_n12,&
               GCOMW_AMSR2L3SND_struc(n)%mod_n21,&
               GCOMW_AMSR2L3SND_struc(n)%mod_n22,&
               GCOMW_AMSR2L3SND_struc(n)%mod_w11,&
               GCOMW_AMSR2L3SND_struc(n)%mod_w12,&
               GCOMW_AMSR2L3SND_struc(n)%mod_w21,&
               GCOMW_AMSR2L3SND_struc(n)%mod_w22)

!          call upscaleByAveraging_input(&
!               GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(:),&
!               LIS_rc%obs_gridDesc(k,:),modis_nc*modis_nr, &
!               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
!               GCOMW_AMSR2L3SND_struc(n)%mod_n11)
       endif

       if(GCOMW_AMSR2L3SND_struc(n)%usr_input_mask.eq.1) then 

          call ESMF_ConfigGetAttribute(LIS_config,input_mask_file,&
               label="AMSR2(GCOMW) snow depth input mask file:",&
               rc=status)
          allocate(obsmask(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
          allocate(GCOMW_AMSR2L3SND_struc(n)%input_mask(&
               LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
          obsmask = -9999.0

          write(LIS_logunit,*) '[INFO] Reading AMSR2 mask ',trim(input_mask_file)
          ftn = LIS_getNextUnitNumber()
          open(ftn,file=input_mask_file,form='unformatted')
          read(ftn) obsmask
          call LIS_releaseUnitNumber(ftn)

          GCOMW_AMSR2L3SND_struc(n)%input_mask(:,:) = obsmask(&
               LIS_ews_obs_halo_ind(k,LIS_localPet+1):&         
               LIS_ewe_obs_halo_ind(k,LIS_localPet+1), &
               LIS_nss_obs_halo_ind(k,LIS_localPet+1): &
               LIS_nse_obs_halo_ind(k,LIS_localPet+1))
          deallocate(obsmask)
       endif

       allocate(GCOMW_AMSR2L3SND_struc(n)%IMSdata_obs(LIS_rc%obs_ngrid(k)))
       allocate(GCOMW_AMSR2L3SND_struc(n)%MODISdata_obs(LIS_rc%obs_ngrid(k)))
    enddo

  end subroutine GCOMW_AMSR2L3SND_setup

end module GCOMW_AMSR2L3SND_Mod
