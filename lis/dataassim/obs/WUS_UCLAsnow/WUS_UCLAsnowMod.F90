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
! !MODULE: WUS_UCLAsnowMod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle the Margulis Western US Snow Reanalysis dataset.
!   Available online at: https://nsidc.org/data/WUS_UCLA_SR/versions/1
! 
! !REVISION HISTORY: 
!  08 Jun 2022: Sujay Kumar; Initial version
! 
module WUS_UCLAsnowMod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: WUS_UCLAsnow_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: WUS_UCLAsnow_struc
!EOP
  type, public:: WUS_UCLAsnow_dec
     
     logical                :: startMode
     logical                :: obsMap
     integer                :: nc     
     integer                :: nr
     integer                :: mi
     integer                :: c_off, r_off
     real                   :: datares
     real                   :: ssdev_inp
     type(proj_info)        :: proj
     real                   :: gridDesci(50)     
     real,    allocatable :: rlat(:)
     real,    allocatable :: rlon(:)
     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)
  end type WUS_UCLAsnow_dec
  
  type(WUS_UCLAsnow_dec),allocatable :: WUS_UCLAsnow_struc(:)
  
contains

!BOP
! 
! !ROUTINE: WUS_UCLAsnow_setup
! \label{WUS_UCLAsnow_setup}
! 
! !INTERFACE: 
  subroutine WUS_UCLAsnow_setup(k, OBS_State, OBS_Pert_State)
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
!   creation of data strctures required for handling  
!   WUS_UCLAsnow data. 
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
    character(len=LIS_CONST_PATH_LEN) ::  snodasobsdir
    character*100          ::  temp
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real, allocatable          :: ssdev(:)
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2    

    allocate(WUS_UCLAsnow_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"WUS UCLA snow snow depth data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,snodasobsdir,&
            rc=status)
       call LIS_verify(status, 'WUS UCLA snow snow depth data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            snodasobsdir, rc=status)
       call LIS_verify(status)
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

    write(LIS_logunit,*)'[INFO] read WUS UCLA snow snow depth data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. WUS_UCLAsnow 
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

          WUS_UCLAsnow_struc(n)%obsMap = .true. 
! Set obs err to be uniform (will be rescaled later for each grid point). 
          ssdev = obs_pert%ssdev(1)
          WUS_UCLAsnow_struc(n)%ssdev_inp = obs_pert%ssdev(1)

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
         '[INFO] Created the States to hold the WUS_UCLAsnow observations data'
    
    do n=1,LIS_rc%nnest

       if(LIS_rc%lis_obs_map_proj(k).eq."latlon") then
          cornerlat1 = max(31.002222, &
               nint((LIS_rc%obs_gridDesc(k,4)-31.002222)/&
               0.0044444)*0.0044444+31.002222-50*0.0044444)
          cornerlon1 = max(-124.99777777778, &
               nint((LIS_rc%obs_gridDesc(k,5)+&
               124.997777778)/0.0044444)*0.0044444-124.99777777778-50*0.0044444)
          cornerlat2 = min(48.997777778, &
               nint((LIS_rc%obs_gridDesc(k,7)-31.002222)/&
               0.0044444)*0.0044444+31.002222+50*0.0044444)
          cornerlon2 = min(-102.002222222222, &
               nint((LIS_rc%obs_gridDesc(k,8)+124.997777778)/&
               0.0044444)*0.0044444-124.99777777778+50*0.0044444)
          if(cornerlat1.gt.48.997777778.or.&
               cornerlat1.lt.31.002222.or.&
               cornerlat2.gt.48.997777778.or.&
               cornerlat2.lt.31.002222.or.&               
               cornerlon1.gt.-102.002222222222.or.&
               cornerlon1.lt.-124.99777777778.or.&
               cornerlon2.gt.-102.002222222222.or.&
               cornerlon2.lt.-124.99777777778) then
             WUS_UCLAsnow_struc(n)%obsMap = .false.             
          endif
          
               
       elseif(LIS_rc%lis_obs_map_proj(k).eq."lambert") then
          cornerlat1 = max(31.002222, &
               nint((LIS_rc%minLat(n)-31.002222)/&
               0.0044444)*0.0044444+31.002222-50*0.0044444)
          cornerlon1 = max(-124.99777777778, &
               nint((LIS_rc%minLon(n)+124.997777778)/&
               0.0044444)*0.0044444-124.99777777778-50*0.0044444)
          cornerlat2 = min(48.997777778, &
               nint((LIS_rc%maxLat(n)-31.002222)/&
               0.0044444)*0.0044444+31.002222+50*0.0044444)
          cornerlon2 = min(-102.002222222222, &
               nint((LIS_rc%maxLon(n)+124.997777778)/&
               0.0044444)*0.0044444-124.99777777778+50*0.0044444)
       endif
       
       WUS_UCLAsnow_struc(n)%c_off = nint((cornerlon1 + 124.997778)/0.0044444)+1
       WUS_UCLAsnow_struc(n)%r_off = nint((cornerlat1 - 31.002222)/0.0044444)+1

              
       WUS_UCLAsnow_struc(n)%nc = nint((cornerlon2-cornerlon1)/0.0044444)+1
       WUS_UCLAsnow_struc(n)%nr = nint((cornerlat2-cornerlat1)/0.0044444)+1
       
       WUS_UCLAsnow_struc(n)%gridDesci(1) = 0 
       WUS_UCLAsnow_struc(n)%gridDesci(2) = WUS_UCLAsnow_struc(n)%nc
       WUS_UCLAsnow_struc(n)%gridDesci(3) = WUS_UCLAsnow_struc(n)%nr 
       WUS_UCLAsnow_struc(n)%gridDesci(4) = cornerlat1
       WUS_UCLAsnow_struc(n)%gridDesci(5) = cornerlon1
       WUS_UCLAsnow_struc(n)%gridDesci(6) = 128
       WUS_UCLAsnow_struc(n)%gridDesci(7) = cornerlat2
       WUS_UCLAsnow_struc(n)%gridDesci(8) = cornerlon2
       WUS_UCLAsnow_struc(n)%gridDesci(9) = 0.0044444
       WUS_UCLAsnow_struc(n)%gridDesci(10) = 0.0044444
       WUS_UCLAsnow_struc(n)%gridDesci(20) = 64

       WUS_UCLAsnow_struc(n)%mi = WUS_UCLAsnow_struc(n)%nc*WUS_UCLAsnow_struc(n)%nr

!-----------------------------------------------------------------------------
!   Use interpolation if LIS is running finer than 500 m. 
!-----------------------------------------------------------------------------
       if(LIS_rc%obs_gridDesc(k,10).le.0.0044444) then 

          allocate(WUS_UCLAsnow_struc(n)%rlat(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%rlon(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%n11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%n12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%n21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%n22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%w11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%w12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%w21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(WUS_UCLAsnow_struc(n)%w22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          
          call bilinear_interp_input_withgrid(WUS_UCLAsnow_struc(n)%gridDesci(:), &
               LIS_rc%obs_gridDesc(k,:),&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
               WUS_UCLAsnow_struc(n)%rlat, WUS_UCLAsnow_struc(n)%rlon,&
               WUS_UCLAsnow_struc(n)%n11, WUS_UCLAsnow_struc(n)%n12, &
               WUS_UCLAsnow_struc(n)%n21, WUS_UCLAsnow_struc(n)%n22, &
               WUS_UCLAsnow_struc(n)%w11, WUS_UCLAsnow_struc(n)%w12, &
               WUS_UCLAsnow_struc(n)%w21, WUS_UCLAsnow_struc(n)%w22)
       else
          
          allocate(WUS_UCLAsnow_struc(n)%n11(&
               WUS_UCLAsnow_struc(n)%nc*WUS_UCLAsnow_struc(n)%nr))

          call upscaleByAveraging_input(WUS_UCLAsnow_struc(n)%gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),&
               WUS_UCLAsnow_struc(n)%nc*WUS_UCLAsnow_struc(n)%nr, &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), WUS_UCLAsnow_struc(n)%n11)
       endif
                     
       WUS_UCLAsnow_struc(n)%datares = 0.0044444

       call LIS_registerAlarm("WUS_UCLAsnow read alarm",&
            86400.0, 86400.0)
       WUS_UCLAsnow_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)

    enddo

   
  end subroutine WUS_UCLAsnow_setup

end module WUS_UCLAsnowMod
