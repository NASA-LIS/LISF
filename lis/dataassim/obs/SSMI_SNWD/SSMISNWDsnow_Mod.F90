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
! !MODULE: SSMISNWDsnow_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle SSMI snow depth retrievals. 
!   
! !REVISION HISTORY: 
!  16 Oct 2012   Sujay Kumar;   Initial Specification
! 
module SSMISNWDsnow_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SSMISNWDsnow_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SSMISNWDsnow_struc

  type, public ::  SSMISNWDsnow_dec
     logical             :: startMode
     real                :: ssdev
     integer             :: snowfield
     integer             :: mi
     integer             :: nc,nr

     integer, allocatable    :: n11(:)
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     real,    allocatable    :: snwd(:)
     real,    allocatable   :: snwdtime(:,:)

     integer             :: useIMS
     character(len=LIS_CONST_PATH_LEN) :: IMSdir
     integer             :: useMODIS
     character(len=LIS_CONST_PATH_LEN) :: MODISdir

     integer             :: ims_mi
     real                :: ims_gridDesci(50)
     integer, allocatable    :: ims_n11(:)

     integer             :: mod_mi
     real                :: mod_gridDesci(50)
     integer, allocatable    :: mod_n11(:)
     
     real, allocatable   :: IMSdata_obs(:)
     real, allocatable   :: MODISdata_obs(:)

  end type SSMISNWDsnow_dec

  type(SSMISNWDsnow_dec), allocatable :: SSMISNWDsnow_struc(:)

contains
!BOP
! 
! !ROUTINE: SSMISNWDsnow_setup
! \label{SSMISNWDsnow_setup}
! 
! !INTERFACE: 
  subroutine SSMISNWDsnow_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_perturbMod
    use LIS_logMod
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
!   creation of data strctures required for processing 
!   SSMI snow depth data. 
!
!   The arguments are: 
!   \begin{description}
!    \item[k]    index of the assimilation instance
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  ssmisnowobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real,         allocatable  ::  varmin(:)
    real,         allocatable  ::  varmax(:)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    real                   :: gridDesci(LIS_rc%nnest,50)
    real, allocatable          :: ssdev(:)
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2    
    real                   :: minlat, minlon, maxlon, maxlat, dx, dy
    integer                ::  modis_nc,modis_nr
    integer                ::  ims_nc, ims_nr

    allocate(SSMISNWDsnow_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SSMI snow depth data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ssmisnowobsdir,&
            rc=status)
       call LIS_verify(status,'SSMI snow depth data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            ssmisnowobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "SSMI snow depth use IMS data for snow detection:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,SSMISNWDsnow_struc(n)%useIMS,&
            rc=status)
       call LIS_verify(status,&
            'SSMI snow depth use IMS data for snow detection: not defined')
       if(SSMISNWDsnow_struc(n)%useIMS.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,&
               SSMISNWDsnow_struc(n)%IMSdir,&
               label="SSMI snow depth IMS data directory:",rc=status)
          call LIS_verify(status,&
               'SSMI snow depth IMS data directory: option not specified')
       endif
    enddo
    

    call ESMF_ConfigFindLabel(LIS_config,&
         "SSMI snow depth use MODIS (MOD10C1) data for snow detection:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,SSMISNWDsnow_struc(n)%useMODIS,&
            rc=status)
       call LIS_verify(status,&
            'SSMI snow depth use MODIS (MOD10C1) data for snow detection: not defined')
       if(SSMISNWDsnow_struc(n)%useMODIS.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,&
               SSMISNWDsnow_struc(n)%MODISdir,&
               label="SSMI snow depth MOD10C1 data directory:",rc=status)
          call LIS_verify(status,&
               'SSMI snow depth MOD10C1 data directory: option not specified')
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
            LIS_rc%obs_ngrid(k),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*)'[INFO] read SSMI snow depth data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. 
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_obsvecGrid(n,k),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
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
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

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

          ssdev = obs_pert%ssdev(1)
          SSMISNWDsnow_struc(n)%ssdev =obs_pert%ssdev(1) 

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

          call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
          call LIS_verify(status)
       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)
    enddo
!-------------------------------------------------------------
! set up the SSMI domain and interpolation weights. 
!-------------------------------------------------------------
    gridDesci = 0 
    do n=1,LIS_rc%nnest
       SSMISNWDsnow_struc(n)%nc = 1440
       SSMISNWDsnow_struc(n)%nr = 360 

       gridDesci(n,1) = 0 
       gridDesci(n,2) = SSMISNWDsnow_struc(n)%nc
       gridDesci(n,3) = SSMISNWDsnow_struc(n)%nr
       gridDesci(n,4) = 0.125
       gridDesci(n,5) = -179.875
       gridDesci(n,6) = 128
       gridDesci(n,7) = 89.875
       gridDesci(n,8) = 179.875
       gridDesci(n,9) = 0.25
       gridDesci(n,10) = 0.25
       gridDesci(n,20) = 64.0

       SSMISNWDsnow_struc(n)%mi = &
            SSMISNWDsnow_struc(n)%nc*SSMISNWDsnow_struc(n)%nr

       allocate(SSMISNWDsnow_struc(n)%n11(SSMISNWDsnow_struc(n)%mi))
       allocate(SSMISNWDsnow_struc(n)%rlat(SSMISNWDsnow_struc(n)%mi))
       allocate(SSMISNWDsnow_struc(n)%rlon(SSMISNWDsnow_struc(n)%mi))

       allocate(SSMISNWDsnow_struc(n)%snwd(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(SSMISNWDsnow_struc(n)%snwdtime(&
            LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k)))

       SSMISNWDsnow_struc(n)%snwdtime = -1
       SSMISNWDsnow_struc(n)%snwd = LIS_rc%udef

       call neighbor_interp_input_withgrid(gridDesci(n,:), & 
            LIS_rc%obs_gridDesc(k,:),&
            SSMISNWDsnow_struc(n)%nc*SSMISNWDsnow_struc(n)%nr,&
            SSMISNWDsnow_struc(n)%rlat,SSMISNWDsnow_struc(n)%rlon,&
            SSMISNWDsnow_struc(n)%n11)

       call LIS_registerAlarm("SSMI snow depth read alarm",&
            86400.0, 86400.0)

       SSMISNWDsnow_struc(n)%startMode = .true. 
       
    enddo
    
    do n=1,LIS_rc%nnest
       SSMISNWDsnow_struc(n)%ims_gridDesci = 0
       if(SSMISNWDsnow_struc(n)%useIMS.eq.1) then 
          
          SSMISNWDsnow_struc(n)%ims_gridDesci(1) = 0 
          SSMISNWDsnow_struc(n)%ims_gridDesci(2) = 1500
          SSMISNWDsnow_struc(n)%ims_gridDesci(3) = 375
          SSMISNWDsnow_struc(n)%ims_gridDesci(4) = 0.12
          SSMISNWDsnow_struc(n)%ims_gridDesci(5) = -179.88
          SSMISNWDsnow_struc(n)%ims_gridDesci(6) = 128
          SSMISNWDsnow_struc(n)%ims_gridDesci(7) = 89.88
          SSMISNWDsnow_struc(n)%ims_gridDesci(8) = 179.88
          SSMISNWDsnow_struc(n)%ims_gridDesci(9) = 0.24
          SSMISNWDsnow_struc(n)%ims_gridDesci(10) = 0.24
          SSMISNWDsnow_struc(n)%ims_gridDesci(20) = 64.0
          
          ims_nc = 1500
          ims_nr = 375
          SSMISNWDsnow_struc(n)%ims_mi = ims_nc*ims_nr
!-----------------------------------------------------------------------------
!   Use upscaling since SSMI data is coarser than IMS
!-----------------------------------------------------------------------------
          allocate(SSMISNWDsnow_struc(n)%ims_n11(ims_nc*ims_nr))
          
          call upscaleByAveraging_input(&
               SSMISNWDsnow_struc(n)%ims_gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),ims_nc*ims_nr, &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), SSMISNWDsnow_struc(n)%ims_n11)

       endif
    enddo

    do n=1,LIS_rc%nnest
       SSMISNWDsnow_struc(n)%mod_gridDesci = 0

       if(SSMISNWDsnow_struc(n)%useMODIS.eq.1) then 
          
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
          
          SSMISNWDsnow_struc(n)%mod_gridDesci(1) = 0 
          SSMISNWDsnow_struc(n)%mod_gridDesci(2) = modis_nc
          SSMISNWDsnow_struc(n)%mod_gridDesci(3) = modis_nr
          SSMISNWDsnow_struc(n)%mod_gridDesci(4) = cornerlat1
          SSMISNWDsnow_struc(n)%mod_gridDesci(5) = cornerlon1
          SSMISNWDsnow_struc(n)%mod_gridDesci(6) = 128
          SSMISNWDsnow_struc(n)%mod_gridDesci(7) = cornerlat2
          SSMISNWDsnow_struc(n)%mod_gridDesci(8) = cornerlon2
          SSMISNWDsnow_struc(n)%mod_gridDesci(9) = 0.05
          SSMISNWDsnow_struc(n)%mod_gridDesci(10) = 0.05
          SSMISNWDsnow_struc(n)%mod_gridDesci(20) = 64.0
          
          SSMISNWDsnow_struc(n)%mod_mi = modis_nc*modis_nr
!-----------------------------------------------------------------------------
!   Use upscaling since SSMI data is coarser than MODIS
!-----------------------------------------------------------------------------
          allocate(SSMISNWDsnow_struc(n)%mod_n11(modis_nc*modis_nr))
          
          call upscaleByAveraging_input(&
               SSMISNWDsnow_struc(n)%mod_gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),modis_nc*modis_nr, &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), SSMISNWDsnow_struc(n)%mod_n11)
       endif

       allocate(SSMISNWDsnow_struc(n)%IMSdata_obs(LIS_rc%obs_ngrid(k)))
       allocate(SSMISNWDsnow_struc(n)%MODISdata_obs(LIS_rc%obs_ngrid(k)))
       
    enddo
    

    write(LIS_logunit,*) '[INFO] Created ESMF States to hold SSMI observations data'

  end subroutine SSMISNWDsnow_setup
  
end module SSMISNWDsnow_Mod

