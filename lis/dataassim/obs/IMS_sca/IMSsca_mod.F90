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
! !MODULE: IMSsca_Mod
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
module IMSsca_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: IMSsca_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: IMSsca_struc

  type, public ::  IMSsca_dec
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

     integer, allocatable    :: n11(:)
     integer, allocatable    :: n12(:)
     integer, allocatable    :: n21(:)
     integer, allocatable    :: n22(:)
     real,    allocatable    :: w11(:)
     real,    allocatable    :: w12(:)
     real,    allocatable    :: w21(:)
     real,    allocatable    :: w22(:)
     real,    allocatable    :: snwd(:)

     integer, allocatable    :: ims_n11(:)
     integer, allocatable    :: ims_n12(:)
     integer, allocatable    :: ims_n21(:)
     integer, allocatable    :: ims_n22(:)
     real,    allocatable    :: ims_w11(:)
     real,    allocatable    :: ims_w12(:)
     real,    allocatable    :: ims_w21(:)
     real,    allocatable    :: ims_w22(:)

  end type IMSsca_dec

  type(IMSsca_dec), allocatable :: IMSsca_struc(:)

contains
!BOP
! 
! !ROUTINE: IMSsca_setup
! \label{IMSsca_setup}
! 
! !INTERFACE: 
  subroutine IMSsca_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_perturbMod
    use LIS_logMod
   
    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for IMS data. 
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
    character(len=LIS_CONST_PATH_LEN) ::  imsobsdir
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


    allocate(IMSsca_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"IMS data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,imsobsdir,&
            rc=status)
       call LIS_verify(status,'IMS data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            imsobsdir, rc=status)
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
            LIS_rc%ngrid(n),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*)'read IMS data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_vecGrid(n),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
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
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

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

          ssdev = obs_pert%ssdev(1)
          IMSsca_struc(n)%ssdev =obs_pert%ssdev(1) 

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

          call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
          call LIS_verify(status)
       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)
    enddo
!-------------------------------------------------------------
! set up the ANSA domain and interpolation weights. 
!-------------------------------------------------------------
    gridDesci = 0 
    
    do n=1,LIS_rc%nnest
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
       
       allocate(IMSsca_struc(n)%ims_n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IMSsca_struc(n)%ims_n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IMSsca_struc(n)%ims_n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IMSsca_struc(n)%ims_n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       
       allocate(IMSsca_struc(n)%ims_w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IMSsca_struc(n)%ims_w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IMSsca_struc(n)%ims_w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IMSsca_struc(n)%ims_w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       
       call bilinear_interp_input(n,gridDesci(n,:), & 
            IMSsca_struc(n)%ims_n11,IMSsca_struc(n)%ims_n12,&
            IMSsca_struc(n)%ims_n21,IMSsca_struc(n)%ims_n22,&
            IMSsca_struc(n)%ims_w11,IMSsca_struc(n)%ims_w12,&
            IMSsca_struc(n)%ims_w21,IMSsca_struc(n)%ims_w22)
       
       call LIS_registerAlarm("IMS read alarm",&
            86400.0, 86400.0)
       
       IMSsca_struc(n)%startMode = .true. 
    enddo
    
    write(LIS_logunit,*) 'Created ESMF States to hold ANSA observations data'

  end subroutine IMSsca_setup
  
end module IMSsca_Mod

