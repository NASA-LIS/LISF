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
! !MODULE: MODISscaobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle modis snow cover area retrievals. This implementation 
!   handles the 5km standard modis product and the 25km gap filled product
!   
! !REVISION HISTORY: 
!  31 Mar 08    Jiarui Dong; Initial Specification
!  10 Sep 08    Sujay Kumar; Adopted in LIS 6.0
!  17 Apr 09    Sujay Kumar; Added support for the gap-filled product
! 
module MODISscaobs_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MODISscaobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!----------------------------------------------------------------------------- 
  PUBLIC :: MODISsca_obs_obj

!EOP
  type, public :: MODISscaobsdec
     integer          :: mi
     real             :: gridDesci(50)    

     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

     integer, allocatable :: stc(:,:)
     integer, allocatable :: str(:,:)
     integer, allocatable :: enc(:,:)
     integer, allocatable :: enr(:,:)
     
     integer          :: use_fill
     real             :: cloud_thres  !choose everything below the threshold
     real             :: cloud_pers_thres !choose everything below the threshold
  end type MODISscaobsdec

  type(MODISscaobsdec), allocatable :: MODISsca_obs_obj(:)
contains

!BOP
! 
! !ROUTINE: MODISscaobs_setup
! \label{MODISscaobs_setup}
! 
! !INTERFACE: 
  subroutine MODISscaobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod,  only : LIS_logunit, LIS_verify
    use LIS_perturbMod

    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for snow cover area 
!   observations to be used in data assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[k]  index of the data assimilation instance
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n 
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  MODISscaobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40           ::  vname(1)
    real                   ::  varmin(1)
    real                   ::  varmax(1)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    integer                ::  modis_nc,modis_nr
    real                   ::  cornerlat1, cornerlon1
    real                   ::  cornerlat2, cornerlon2

    allocate(MODISsca_obs_obj(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"MODIS SCF data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       MODISsca_obs_obj(n)%gridDesci = 0 
       call ESMF_ConfigGetAttribute(LIS_config,MODISscaobsdir,&
            rc=status)
       if(status /= ESMF_SUCCESS)then
         write(LIS_logunit,*)"MODIS SCF data directory is missing"
       end if
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            MODISscaobsdir, rc=status)
       call LIS_verify(status)

    enddo

   call ESMF_ConfigFindLabel(LIS_config,"MODIS SCF use gap filled product:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MODISsca_obs_obj(n)%use_fill,&
            rc=status)
       call LIS_verify(status,'MODIS SCF use gap filled product: is not specified')
    enddo

   call ESMF_ConfigFindLabel(LIS_config,"MODIS SCF cloud threshold:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MODISsca_obs_obj(n)%cloud_thres,&
            rc=status)
       call LIS_verify(status,'MODIS SCF cloud threshold: is not specified')
    enddo

   call ESMF_ConfigFindLabel(LIS_config,"MODIS SCF cloud persistence threshold:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MODISsca_obs_obj(n)%cloud_pers_thres,&
            rc=status)
       call LIS_verify(status,'MODIS SCF cloud persistence threshold: is not specified')
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

    write(LIS_logunit,*)'read modis snow cover fraction data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. modis 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,grid=LIS_vecGrid(n),&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)

!Perturbations State
       call LIS_readDAObsAttributes(k,vname,varmin,varmax)
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)


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
          
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
               obs_pert%ssdev(1), rc=status)
          call LIS_verify(status)

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
          
       cornerlat1 = max(-59.975, nint((LIS_rc%gridDesc(n,4)+59.975)/0.05)*0.05-59.975-2*0.05)
       cornerlon1 = max(-179.975, nint((LIS_rc%gridDesc(n,5)+179.975)/0.05)*0.05-179.975-2*0.05)
       cornerlat2 = min(89.975, nint((LIS_rc%gridDesc(n,7)+59.975)/0.05)*0.05-59.975+2*0.05)
       cornerlon2 = min(179.975, nint((LIS_rc%gridDesc(n,8)+179.975)/0.05)*0.05-179.975+2*0.05)


       modis_nr = nint((cornerlat2-cornerlat1)/0.05)+1
       modis_nc = nint((cornerlon2-cornerlon1)/0.05)+1
          
       MODISsca_obs_obj(n)%gridDesci(1) = 0 
       MODISsca_obs_obj(n)%gridDesci(2) = modis_nc
       MODISsca_obs_obj(n)%gridDesci(3) = modis_nr
       MODISsca_obs_obj(n)%gridDesci(4) = cornerlat1
       MODISsca_obs_obj(n)%gridDesci(5) = cornerlon1
       MODISsca_obs_obj(n)%gridDesci(6) = 128
       MODISsca_obs_obj(n)%gridDesci(7) = cornerlat2
       MODISsca_obs_obj(n)%gridDesci(8) = cornerlon2
       MODISsca_obs_obj(n)%gridDesci(9) = 0.05
       MODISsca_obs_obj(n)%gridDesci(10) = 0.05
       MODISsca_obs_obj(n)%gridDesci(20) = 64

       MODISsca_obs_obj(n)%mi = modis_nc*modis_nr

!-----------------------------------------------------------------------------
!   Use interpolation if LIS is running finer than 5km. 
!-----------------------------------------------------------------------------
       if(LIS_rc%gridDesc(n,10).le.0.05) then 

          allocate(MODISsca_obs_obj(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MODISsca_obs_obj(n)%n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MODISsca_obs_obj(n)%n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MODISsca_obs_obj(n)%n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MODISsca_obs_obj(n)%w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MODISsca_obs_obj(n)%w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MODISsca_obs_obj(n)%w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MODISsca_obs_obj(n)%w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call bilinear_interp_input(n, MODISsca_obs_obj(n)%gridDesci(:), &
               MODISsca_obs_obj(n)%n11, MODISsca_obs_obj(n)%n12, &
               MODISsca_obs_obj(n)%n21, MODISsca_obs_obj(n)%n22, &
               MODISsca_obs_obj(n)%w11, MODISsca_obs_obj(n)%w12, &
               MODISsca_obs_obj(n)%w21, MODISsca_obs_obj(n)%w22)
       else
          
          allocate(MODISsca_obs_obj(n)%n11(modis_nc*modis_nr))

          call upscaleByAveraging_input(MODISsca_obs_obj(n)%gridDesci(:),&
               LIS_rc%gridDesc(n,:),modis_nc*modis_nr, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n), MODISsca_obs_obj(n)%n11)
       endif
    enddo


    write(LIS_logunit,*) 'Created the States to hold the observations data'
    
  end subroutine MODISscaobs_setup

end module MODISscaobs_module
