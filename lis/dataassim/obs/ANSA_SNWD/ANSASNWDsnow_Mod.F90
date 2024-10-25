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
! !MODULE: ANSASNWDsnow_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle ANSA (AFWA NASA snow algorithm) snow depth retrievals. 
!   The snow depth data is primarily based on AMSR-E, with observation
!   gaps filled in using previous overpasses. For more details, please see 
!   the ANSA reference: 
!
!   Foster et al. "A blended global snow product using visible, 
!   passive microwave and scatterometer satellite data", International
!   Journal of Remote Sensing, 2010. 
!   
! !REVISION HISTORY: 
!  1 Jun 09   Sujay Kumar;   Initial Specification
! 
module ANSASNWDsnow_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN 
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ANSASNWDsnow_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ANSASNWDsnow_struc

  type, public ::  ANSASNWDsnow_dec
     logical             :: startMode
     real                :: ssdev
     integer             :: snowfield
     integer             :: mi
     real                :: cornerlat1,cornerlat2
     real                :: cornerlon1,cornerlon2
     real                :: minlat,minlon
     real                :: maxlat,maxlon
     integer             :: nc,nr
     integer             :: offset1, offset2
     real                :: gridDesc(6)
     integer             :: useIMS
     character(len=LIS_CONST_PATH_LEN) :: IMSdir, MODISdir
     integer             :: useMODIS
     integer, allocatable    :: n11(:)
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     real,    allocatable    :: snwd(:)
     real,    allocatable    :: snwdtime(:,:)

     real, allocatable    :: ims_rlat(:) ! EMK changed type to real
     real, allocatable    :: ims_rlon(:) ! EMK changed type to real
     integer, allocatable    :: ims_n11(:)
     integer, allocatable    :: ims_n12(:)
     integer, allocatable    :: ims_n21(:)
     integer, allocatable    :: ims_n22(:)
     real,    allocatable    :: ims_w11(:)
     real,    allocatable    :: ims_w12(:)
     real,    allocatable    :: ims_w21(:)
     real,    allocatable    :: ims_w22(:)

     integer             :: mod_mi
     integer             :: ims_mi
     real                :: mod_gridDesci(50)

     integer, allocatable    :: mod_n11(:)
     real, allocatable    :: mod_rlat(:)
     real, allocatable    :: mod_rlon(:)

     real,    allocatable    :: snwd_flag(:)
     real, allocatable   :: IMSdata_obs(:)
     real, allocatable   :: MODISdata_obs(:)
  end type ANSASNWDsnow_dec

  type(ANSASNWDsnow_dec), allocatable :: ANSASNWDsnow_struc(:)

contains
!BOP
! 
! !ROUTINE: ANSASNWDsnow_setup
! \label{ANSASNWDsnow_setup}
! 
! !INTERFACE: 
  subroutine ANSASNWDsnow_setup(k, OBS_State, OBS_Pert_State)
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
!   creation of data strctures required for ANSA snow depth data. 
!
!   The arguments are: 
!   \begin{description}
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
    character(len=LIS_CONST_PATH_LEN) ::  ansasnowobsdir
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
    real                   :: dx, dy
    integer                ::  modis_nc,modis_nr
    real                   ::  cornerlat1, cornerlon1
    real                   ::  cornerlat2, cornerlon2

    integer                 :: c,r
    real,       allocatable :: obserr(:,:), lobserr(:,:)

    allocate(ANSASNWDsnow_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ansasnowobsdir,&
            rc=status)
       call LIS_verify(status,'ANSA snow depth data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            ansasnowobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "ANSA snow depth use IMS data for snow detection:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%useIMS,&
            rc=status)
       call LIS_verify(status,&
            'ANSA snow depth use IMS data for snow detection: not defined')
       if(ANSASNWDsnow_struc(n)%useIMS.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,&
               ANSASNWDsnow_struc(n)%IMSdir,&
               label="ANSA snow depth IMS data directory:",rc=status)
          call LIS_verify(status,&
               'ANSA snow depth IMS data directory: option not specified')
       endif
    enddo
    

    call ESMF_ConfigFindLabel(LIS_config,&
         "ANSA snow depth use MODIS (MOD10C1) data for snow detection:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%useMODIS,&
            rc=status)
       call LIS_verify(status,&
            'ANSA snow depth use MODIS (MOD10C1) data for snow detection: not defined')
       if(ANSASNWDsnow_struc(n)%useMODIS.eq.1) then 
          call ESMF_ConfigGetAttribute(LIS_config,&
               ANSASNWDsnow_struc(n)%MODISdir,&
               label="ANSA snow depth MOD10C1 data directory:",rc=status)
          call LIS_verify(status,&
               'ANSA snow depth MOD10C1 data directory: option not specified')
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

    write(LIS_logunit,*)'[INFO] read ANSA snow depth data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%obs_ngrid(k). 
!   
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
          ANSASNWDsnow_struc(n)%ssdev =obs_pert%ssdev(1) 

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsensongrid(n,k),name="Observation"//vid(1)//vid(2),&
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

#if 0 
    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       allocate(obserr(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(lobserr(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       obserr = -9999.0
       lobserr = -9999.0
       
       open(100,file='PMSNWD_obserr.bin',form='unformatted')
       read(100) obserr
       close(100)
       
       lobserr(:,:) = obserr(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1))
    
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
             if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                ssdev(LIS_domain(n)%gindex(c,r))  = lobserr(c,r)                
             endif
          enddo
       enddo

       deallocate(obserr)
       deallocate(lobserr)
    enddo
#endif
!-------------------------------------------------------------
! set up the ANSA domain and interpolation weights. 
!-------------------------------------------------------------
    gridDesci = 0 
    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth lower left lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%gridDesc(1),&
            rc=status)
       call LIS_verify(status,'ANSA snow depth lower left lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth lower left lon:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%gridDesc(2),&
            rc=status)
       call LIS_verify(status,'ANSA snow depth lower left lon: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth upper right lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%gridDesc(3),&
            rc=status)
       call LIS_verify(status,'ANSA snow depth upper right lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth upper right lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%gridDesc(3),&
            rc=status)
       call LIS_verify(status,'ANSA snow depth upper right lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth upper right lon:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%gridDesc(4),&
            rc=status)
       call LIS_verify(status,'ANSA snow depth upper right lon: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth resolution (dx):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%gridDesc(5),&
            rc=status)
       call LIS_verify(status,'ANSA snow depth resolution (dx): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA snow depth resolution (dy):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASNWDsnow_struc(n)%gridDesc(6),&
            rc=status)
       call LIS_verify(status,'ANSA snow depth resolution (dy): not defined')
    enddo

    do n=1,LIS_rc%nnest

       ANSASNWDsnow_struc(n)%minlat = ANSASNWDsnow_struc(n)%gridDesc(1)
       ANSASNWDsnow_struc(n)%minlon = ANSASNWDsnow_struc(n)%gridDesc(2)
       ANSASNWDsnow_struc(n)%maxlat = ANSASNWDsnow_struc(n)%gridDesc(3)
       ANSASNWDsnow_struc(n)%maxlon = ANSASNWDsnow_struc(n)%gridDesc(4)
       dx = ANSASNWDsnow_struc(n)%gridDesc(5)
       dy = ANSASNWDsnow_struc(n)%gridDesc(6)
    enddo
    
    do n=1, LIS_rc%nnest
!sets the local domain corner points with additional buffer 
       ANSASNWDsnow_struc(n)%cornerlat1 = &
            max(ANSASNWDsnow_struc(n)%minlat, &
            nint((LIS_domain(n)%minlat-&
            ANSASNWDsnow_struc(n)%minlat)/dx)*dx+&
            ANSASNWDsnow_struc(n)%minlat-2*dx)
       ANSASNWDsnow_struc(n)%cornerlon1 = &
            max(ANSASNWDsnow_struc(n)%minlon, &
            nint((LIS_domain(n)%minlon-&
            ANSASNWDsnow_struc(n)%minlon)/dy)*dy+&
            ANSASNWDsnow_struc(n)%minlon-2*dy)
       ANSASNWDsnow_struc(n)%cornerlat2 = &
            min(ANSASNWDsnow_struc(n)%maxlat, &
            nint((LIS_domain(n)%maxlat-&
            ANSASNWDsnow_struc(n)%minlat)/dx)*dx+&
            ANSASNWDsnow_struc(n)%minlat+2*dx)
       ANSASNWDsnow_struc(n)%cornerlon2 = &
            min(ANSASNWDsnow_struc(n)%maxlon, &
            nint((LIS_domain(n)%maxlon-&
            ANSASNWDsnow_struc(n)%minlon)/dy)*dy+&
            ANSASNWDsnow_struc(n)%minlon+2*dy)


       ANSASNWDsnow_struc(n)%offset1 = &
            nint((ANSASNWDsnow_struc(n)%cornerlon1-&
            ANSASNWDsnow_struc(n)%minlon)/dy)
       ANSASNWDsnow_struc(n)%offset2 = &
            nint((ANSASNWDsnow_struc(n)%cornerlat1-&
            ANSASNWDsnow_struc(n)%minlat)/dx)

       ANSASNWDsnow_struc(n)%nr = &
            nint((ANSASNWDsnow_struc(n)%cornerlat2-&
            ANSASNWDsnow_struc(n)%cornerlat1)/dx)+1
       ANSASNWDsnow_struc(n)%nc = &
            nint((ANSASNWDsnow_struc(n)%cornerlon2-&
            ANSASNWDsnow_struc(n)%cornerlon1)/dy)+1

       gridDesci(n,1) = 0 
       gridDesci(n,2) = ANSASNWDsnow_struc(n)%nc
       gridDesci(n,3) = ANSASNWDsnow_struc(n)%nr
       gridDesci(n,4) = ANSASNWDsnow_struc(n)%cornerlat1
       gridDesci(n,5) = ANSASNWDsnow_struc(n)%cornerlon1
       gridDesci(n,6) = 128
       gridDesci(n,7) = ANSASNWDsnow_struc(n)%cornerlat2
       gridDesci(n,8) = ANSASNWDsnow_struc(n)%cornerlon2
       gridDesci(n,9) = ANSASNWDsnow_struc(n)%gridDesc(5)
       gridDesci(n,10) = ANSASNWDsnow_struc(n)%gridDesc(6)
       gridDesci(n,20) = 64.0

       ANSASNWDsnow_struc(n)%mi = &
            ANSASNWDsnow_struc(n)%nc*ANSASNWDsnow_struc(n)%nr

       allocate(ANSASNWDsnow_struc(n)%n11(ANSASNWDsnow_struc(n)%mi))

       allocate(ANSASNWDsnow_struc(n)%rlat(ANSASNWDsnow_struc(n)%mi))
       allocate(ANSASNWDsnow_struc(n)%rlon(ANSASNWDsnow_struc(n)%mi))

       allocate(ANSASNWDsnow_struc(n)%snwd(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(ANSASNWDsnow_struc(n)%snwdtime(&
            LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k)))
       
       ANSASNWDsnow_struc(n)%snwdtime = -1
       ANSASNWDsnow_struc(n)%snwd = LIS_rc%udef

       call neighbor_interp_input_withgrid(gridDesci(n,:),& 
            LIS_rc%obs_gridDesc(k,:),&
            ANSASNWDsnow_struc(n)%mi, &
            ANSASNWDsnow_struc(n)%rlat,ANSASNWDsnow_struc(n)%rlon,&
            ANSASNWDsnow_struc(n)%n11)

       gridDesci = 0
       if(ANSASNWDsnow_struc(n)%useIMS.eq.1) then 
          
          gridDesci(n,1) = 0 
          gridDesci(n,2) = 1500
          gridDesci(n,3) = 375
          gridDesci(n,4) = 0.12
          gridDesci(n,5) = -179.88
          gridDesci(n,6) = 128
          gridDesci(n,7) = 89.88
          gridDesci(n,8) = 179.88
          gridDesci(n,9) = 0.24
          gridDesci(n,10) = 0.24
          gridDesci(n,20) = 64.0

          ANSASNWDsnow_struc(n)%ims_mi = gridDesci(n,2)*gridDesci(n,3)
          
          allocate(ANSASNWDsnow_struc(n)%ims_rlat(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_rlon(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_n11(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_n12(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_n21(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_n22(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          
          allocate(ANSASNWDsnow_struc(n)%ims_w11(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_w12(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_w21(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(ANSASNWDsnow_struc(n)%ims_w22(&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          

          call bilinear_interp_input_withgrid(gridDesci(n,:),&
               LIS_rc%obs_gridDesc(k,:),&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
               ANSASNWDsnow_struc(n)%ims_rlat, &
               ANSASNWDsnow_struc(n)%ims_rlon, &
               ANSASNWDsnow_struc(n)%ims_n11,ANSASNWDsnow_struc(n)%ims_n12,&
               ANSASNWDsnow_struc(n)%ims_n21,ANSASNWDsnow_struc(n)%ims_n22,&
               ANSASNWDsnow_struc(n)%ims_w11,ANSASNWDsnow_struc(n)%ims_w12,&
               ANSASNWDsnow_struc(n)%ims_w21,ANSASNWDsnow_struc(n)%ims_w22)
       endif

       gridDesci = 0
       if(ANSASNWDsnow_struc(n)%useMODIS.eq.1) then 
          
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
          
          ANSASNWDsnow_struc(n)%mod_gridDesci(1) = 0 
          ANSASNWDsnow_struc(n)%mod_gridDesci(2) = modis_nc
          ANSASNWDsnow_struc(n)%mod_gridDesci(3) = modis_nr
          ANSASNWDsnow_struc(n)%mod_gridDesci(4) = cornerlat1
          ANSASNWDsnow_struc(n)%mod_gridDesci(5) = cornerlon1
          ANSASNWDsnow_struc(n)%mod_gridDesci(6) = 128
          ANSASNWDsnow_struc(n)%mod_gridDesci(7) = cornerlat2
          ANSASNWDsnow_struc(n)%mod_gridDesci(8) = cornerlon2
          ANSASNWDsnow_struc(n)%mod_gridDesci(9) = 0.05
          ANSASNWDsnow_struc(n)%mod_gridDesci(10) = 0.05
          ANSASNWDsnow_struc(n)%mod_gridDesci(20) = 64.0
          
          ANSASNWDsnow_struc(n)%mod_mi = modis_nc*modis_nr
!-----------------------------------------------------------------------------
!   Use neighbor interpolation since the ANSA grid and the obs grid
!   are at 5km
!-----------------------------------------------------------------------------
          
          allocate(ANSASNWDsnow_struc(n)%mod_n11(ANSASNWDsnow_struc(n)%mod_mi))
          allocate(ANSASNWDsnow_struc(n)%mod_rlat(ANSASNWDsnow_struc(n)%mod_mi))
          allocate(ANSASNWDsnow_struc(n)%mod_rlon(ANSASNWDsnow_struc(n)%mod_mi))
          
          call neighbor_interp_input_withgrid(&
               ANSASNWDsnow_struc(n)%mod_gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),&
               ANSASNWDsnow_struc(n)%mod_mi, &
               ANSASNWDsnow_struc(n)%mod_rlat,&
               ANSASNWDsnow_struc(n)%mod_rlon,&
               ANSASNWDsnow_struc(n)%mod_n11)
          
       endif
       call LIS_registerAlarm("ANSA snow depth read alarm",&
            86400.0, 86400.0)

       ANSASNWDsnow_struc(n)%startMode = .true. 

       allocate(ANSASNWDsnow_struc(n)%IMSdata_obs(LIS_rc%obs_ngrid(k)))
       allocate(ANSASNWDsnow_struc(n)%MODISdata_obs(LIS_rc%obs_ngrid(k)))
    enddo

    write(LIS_logunit,*) '[INFO] Created ESMF States to hold ANSA observations data'

  end subroutine ANSASNWDsnow_setup
  
end module ANSASNWDsnow_Mod

