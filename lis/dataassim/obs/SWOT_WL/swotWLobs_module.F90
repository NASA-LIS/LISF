!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
#include "LIS_NetCDF_inc.h"
!BOP
!
! !MODULE: swotWLobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle SWOT water level observations
!   
! !REVISION HISTORY: 
!  15 Apr 2024    Yeosang Yoon;   Initial Specification
!  25 Nov 2025: Yeosang Yoon; Clean up the code
! 
module swotWLobs_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: swotWLobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: swot_wl_struc

  type, public :: swot_wl_dec
     type(ESMF_TimeInterval) :: ts
     integer                 :: nlisgrids
     integer                 :: ntimes
     logical                 :: readflag
     real,   allocatable     :: time(:,:)
     real,   allocatable     :: WLobs(:,:)
     integer,allocatable     :: lisid(:,:)

     real,  allocatable      :: model_mu(:,:)
     real,  allocatable      :: obs_mu(:,:)
     real,  allocatable      :: model_sigma(:,:)
     real,  allocatable      :: obs_sigma(:,:)

     integer                 :: nt
     real                    :: daInterval
     integer                 :: time_file
  end type swot_wl_dec
  type(swot_wl_dec), allocatable :: swot_wl_struc(:)
contains
!BOP
! 
! !ROUTINE: swotWLobs_setup
! \label{swotWLobs_setup}
! 
! !INTERFACE: 
  subroutine swotWLobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_dataAssimMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_timeMgrMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for water level 
!   assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                           ::  n 
    integer                           ::  status
    type(ESMF_Field)                  ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)              ::  intarrspec, realarrspec
    type(ESMF_Field)                  ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)              ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  swotwlobsdir
    character(len=LIS_CONST_PATH_LEN) ::  wlmap
    character(len=LIS_CONST_PATH_LEN) ::  temp
    integer                           ::  ftn,i,t
    integer                           ::  varid
    integer                           ::  ngrid
    real,  allocatable                ::  obsstd(:)
    character*1                       ::  vid(2)
    character*40, allocatable         ::  vname(:)
    real, allocatable                 ::  varmin(:)
    real, allocatable                 ::  varmax(:)
    type(pert_dec_type)               ::  obs_pert
    real, pointer                     ::  obs_temp(:,:)
    real, allocatable                 ::  ssdev(:)
    character*10                      ::  time

    allocate(swot_wl_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SWOT water level data directory:", rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,swotwlobsdir,rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",swotwlobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"SWOT water level map:", rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wlmap,rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Data assimilation frequency:", rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,"Data assimilation frequency: not defined")

       call LIS_parseTimeString(time,swot_wl_struc(n)%daInterval)
    enddo

    do n=1,LIS_rc%nnest

       call ESMF_AttributeSet(OBS_State(n),"Data Update Status", .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Update Time", -99.0, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status", .false., rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
            LIS_rc%obs_ngrid(k),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*) '[INFO] Reading SWOT water level map ',trim(wlmap)
    do n=1,LIS_rc%nnest

       call ESMF_TimeIntervalSet(swot_wl_struc(n)%ts,d=1,rc=status)
       call LIS_verify(status, 'Error in ESMF_TimeIntervalSet in swotwlobs')

       allocate(swot_wl_struc(n)%lisid(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       swot_wl_struc(n)%lisid = -9999

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
       call LIS_verify(nf90_open(path=trim(wlmap),mode=NF90_NOWRITE,ncid=ftn),&
            '[ERR] Error opening file '//trim(wlmap))
       call LIS_verify(nf90_inq_varid(ftn,'lis_id',varid), '[ERR] Error with nf90_inq_varid: lis_id')
       call LIS_verify(nf90_get_var(ftn,varid, swot_wl_struc(n)%lisid, &
            start=(/LIS_ews_obs_halo_ind(n,LIS_localPet+1), LIS_nss_obs_halo_ind(n,LIS_localPet+1)/),&
            count = (/LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)/)),'[ERR] Error with nf90_get_var: lis_id')
       call LIS_verify(nf90_close(ftn))
#endif
    enddo
!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. 
!----------------------------------------------------------------------------
    write(LIS_logunit,*) '[INFO] Reading SWOT water level data specifications'

    do n=1,LIS_rc%nnest

       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid
       
       obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_obsVecGrid(n,k),&
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
               grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
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

       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
!to obsarr

    enddo
    write(LIS_logunit,*) '[INFO] Created the States to hold the observations data'
    

    do n=1,LIS_rc%nnest
       allocate(ssdev(LIS_rc%obs_ngrid(k)))
       ssdev = obs_pert%ssdev(1)
       
       if(LIS_rc%obs_ngrid(k).gt.0) then
          call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
          call LIS_verify(status, 'Error in AttributeSet: Standard Deviation')
       endif
 
       deallocate(ssdev)
    enddo
    
    do n=1,LIS_rc%nnest

       call LIS_registerAlarm("swotWL read alarm", 86400.0, 86400.0) !daily input

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status)
    enddo

  end subroutine swotWLobs_setup
  
end module swotWLobs_module
