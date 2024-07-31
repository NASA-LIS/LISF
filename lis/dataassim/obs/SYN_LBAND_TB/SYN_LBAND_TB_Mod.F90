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
! !MODULE: SYN_LBAND_TB_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle synthetic (simulated) L-band brightness temperature observations
!   This plugin is based on a simulated data over CONUS at 1deg resolution
!   The domain size is assumed to be 58x29. 
!   
! !REVISION HISTORY: 
!  13 Sept 2012: Sujay Kumar, initial specification
! 
module SYN_LBAND_TB_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SYN_LBAND_TB_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: SYN_LBAND_TB_struc
!EOP
  type, public:: SYN_LBAND_TB_dec

!     integer              :: nc
!     integer              :: nr
     real,     allocatable    :: TbH(:,:)
     real,     allocatable    :: TbV(:,:)
  end type SYN_LBAND_TB_dec
  
  type(SYN_LBAND_TB_dec),allocatable :: SYN_LBAND_TB_struc(:)
  
contains

!BOP
! 
! !ROUTINE: SYN_LBAND_TB_setup
! \label{SYN_LBAND_TB_setup}
! 
! !INTERFACE: 
  subroutine SYN_LBAND_TB_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid, LIS_ensOnGrid
    use LIS_timeMgrMod, only : LIS_clock, LIS_calendar, LIS_registerAlarm
    use LIS_perturbMod
    use LIS_logmod, only : LIS_logunit, LIS_verify, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber

    implicit none 

! !ARGUMENTS: 
    integer                ::  k
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for synthetic L-band 
!   brightness temperature. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n, i,kk
    integer                ::  m
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  lbandtbobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40,allocatable   ::  vname(:)
    real        ,allocatable   ::  varmin(:)
    real        ,allocatable   ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real, allocatable          :: ssdev(:)

    allocate(SYN_LBAND_TB_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"Synthetic L-band Tb data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,lbandtbobsdir,&
            rc=status)
       call LIS_verify(status, 'Synthetic L-band Tb data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            lbandtbobsdir, rc=status)
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

    write(LIS_logunit,*)'read Synthetic L-band Tb data specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. Observations are in the grid space. 
!   There are two types of observations (for the H-pol and V-pol)
!   being assimilated.
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
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

       do m=1,LIS_rc%nobtypes(k)
          write(unit=temp,fmt='(i2.2)') m
          read(unit=temp,fmt='(2a1)') vid

          obsField = ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=LIS_vecGrid(n),&
               name="Observation"//vid(1)//vid(2),rc=status)
          call LIS_verify(status)
          
          call ESMF_StateAdd(OBS_State(n),(/obsField/),rc=status)
          call LIS_verify(status)
       enddo

       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
          allocate(obs_pert%vname(LIS_rc%nobtypes(k)))
          allocate(obs_pert%perttype(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ssdev(LIS_rc%nobtypes(k)))
          allocate(obs_pert%stdmax(LIS_rc%nobtypes(k)))
          allocate(obs_pert%zeromean(LIS_rc%nobtypes(k)))
          allocate(obs_pert%tcorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%xcorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ycorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ccorr(LIS_rc%nobtypes(k),LIS_rc%nobtypes(k)))

          call LIS_readPertAttributes(LIS_rc%nobtypes(k),&
               LIS_rc%obspertAttribfile(k),&
               obs_pert)

          allocate(ssdev(LIS_rc%ngrid(n)))
          do m=1,LIS_rc%nobtypes(k)
             write(unit=temp,fmt='(i2.2)') m
             read(unit=temp,fmt='(2a1)') vid

             pertField = ESMF_FieldCreate(arrayspec=pertArrSpec,&
                  grid=LIS_ensOnGrid(n),name="Observation"//vid(1)//vid(2),&
                  rc=status)
             call LIS_verify(status)
          
! initializing the perturbations to be zero 
             call ESMF_FieldGet(pertField,localDE=0,farrayPtr=obs_temp,rc=status)
             call LIS_verify(status)
             obs_temp(:,:) = 0 
             
             call ESMF_AttributeSet(pertField,"Perturbation Type",&
                  obs_pert%perttype(m), rc=status)
             call LIS_verify(status)
             

             ssdev = obs_pert%ssdev(m)

             if(LIS_rc%ngrid(n).gt.0) then 
                call ESMF_AttributeSet(pertField,"Standard Deviation",&
                     ssdev, itemCount=LIS_rc%ngrid(n),rc=status)
                call LIS_verify(status)
             endif

             call ESMF_AttributeSet(pertField,"Std Normal Max",&
                  obs_pert%stdmax(m), rc=status)
             call LIS_verify(status)
             
             call ESMF_AttributeSet(pertField,"Ensure Zero Mean",&
                  obs_pert%zeromean(m),rc=status)
             call LIS_verify(status)
             
             call ESMF_AttributeSet(pertField,"Temporal Correlation Scale",&
                  obs_pert%tcorr(m),rc=status)
             call LIS_verify(status)
             
             call ESMF_AttributeSet(pertField,"X Correlation Scale",&
                  obs_pert%xcorr(m),rc=status)
             
             call ESMF_AttributeSet(pertField,"Y Correlation Scale",&
                  obs_pert%ycorr(m),rc=status)
             
             call ESMF_AttributeSet(pertField,"Cross Correlation Strength",&
                  obs_pert%ccorr(m,:),itemCount=LIS_rc%nobtypes(k),rc=status)
             
             call ESMF_StateAdd(OBS_Pert_State(n),(/pertField/),rc=status)
             call LIS_verify(status)
          enddo
          deallocate(ssdev)
       endif

       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
!to obsarr

    enddo
    write(LIS_logunit,*) 'Created the States to hold the Synthetic L-band Tb observations data'
    do n=1,LIS_rc%nnest
!       SYN_LBAND_TB_struc(n)%nc = LIS_rc%gnc(n)
!       SYN_LBAND_TB_struc(n)%nr = LIS_rc%gnr(n)
       allocate(SYN_LBAND_TB_struc(n)%TbH(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(SYN_LBAND_TB_struc(n)%TbV(LIS_rc%gnc(n),LIS_rc%gnr(n)))

       call LIS_registerAlarm("SYN LBAND TB read alarm",&
            86400.0, 86400.0)
    enddo
  end subroutine SYN_LBAND_TB_setup
end module SYN_LBAND_TB_Mod
