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
! !MODULE: ANSASWEsnow_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle ANSA (AFWA NASA snow algorithm) SWE retrievals. 
!   The SWE data is primarily based on AMSR-E, with observation
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
module ANSASWEsnow_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ANSASWEsnow_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ANSASWEsnow_struc

  type, public ::  ANSASWEsnow_dec
     integer             :: snowfield
     integer             :: mi
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
     real,    allocatable    :: swe(:)
  end type ANSASWEsnow_dec

  type(ANSASWEsnow_dec), allocatable :: ANSASWEsnow_struc(:)

contains
!BOP
! 
! !ROUTINE: ANSASWEsnow_setup
! \label{ANSASWEsnow_setup}
! 
! !INTERFACE: 
  subroutine ANSASWEsnow_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_config, LIS_vecGrid, LIS_ensOnGrid, &
         LIS_masterproc    
    use LIS_timeMgrMod, only : LIS_clock, LIS_calendar, LIS_registerAlarm
    use LIS_historyMod, only : LIS_readvar_gridded
    use LIS_perturbMod
    use LIS_logMod, only : LIS_logunit, LIS_verify, LIS_getNextUnitNumber, &
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
!   creation of data strctures required for ANSA snow data. 
!
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer, parameter     ::  nbins = 1000 
    integer                ::  n
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  ansasnowobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    real                   ::  varmin(1)
    real                   ::  varmax(1)
    character*40           ::  vname(1)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          :: obs_temp(:,:)
    real                   :: gridDesci(LIS_rc%nnest,50)
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2
    real                   :: minlat, minlon, maxlon, maxlat, dx, dy

    allocate(ANSASWEsnow_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ansasnowobsdir,&
            rc=status)
       call LIS_verify(status,'ANSA snow data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            ansasnowobsdir, rc=status)
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

    write(LIS_logunit,*)'read ANSA SWE data specifications'

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
          
!to obsarr

    enddo
!-------------------------------------------------------------
! setting up interpolation we%ights for neighbor search
!-------------------------------------------------------------
    gridDesci = 0 
    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE lower left lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASWEsnow_struc(n)%gridDesc(1),&
            rc=status)
       call LIS_verify(status,'ANSA SWE lower left lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE lower left lon:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASWEsnow_struc(n)%gridDesc(2),&
            rc=status)
       call LIS_verify(status,'ANSA SWE lower left lon: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE upper right lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASWEsnow_struc(n)%gridDesc(3),&
            rc=status)
       call LIS_verify(status,'ANSA SWE upper right lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE upper right lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASWEsnow_struc(n)%gridDesc(3),&
            rc=status)
       call LIS_verify(status,'ANSA SWE upper right lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE upper right lon:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASWEsnow_struc(n)%gridDesc(4),&
            rc=status)
       call LIS_verify(status,'ANSA SWE upper right lon: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE resolution (dx):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASWEsnow_struc(n)%gridDesc(5),&
            rc=status)
       call LIS_verify(status,'ANSA SWE resolution (dx): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SWE resolution (dy):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASWEsnow_struc(n)%gridDesc(6),&
            rc=status)
       call LIS_verify(status,'ANSA SWE resolution (dy): not defined')
    enddo

    do n=1,LIS_rc%nnest

       minlat = ANSASWEsnow_struc(n)%gridDesc(1)
       minlon = ANSASWEsnow_struc(n)%gridDesc(2)
       maxlat = ANSASWEsnow_struc(n)%gridDesc(3)
       maxlon = ANSASWEsnow_struc(n)%gridDesc(4)
       dx = ANSASWEsnow_struc(n)%gridDesc(5)
       dy = ANSASWEsnow_struc(n)%gridDesc(6)

    enddo
    
    do n=1, LIS_rc%nnest
!sets the local domain corner points with additional buffer 
       cornerlat1 = max(minlat, nint((LIS_domain(n)%minlat-minlat)/dx)*dx+minlat-2*dx)
       cornerlon1 = max(minlon, nint((LIS_domain(n)%minlon-minlon)/dy)*dy+minlon-2*dy)
       cornerlat2 = min(maxlat, nint((LIS_domain(n)%maxlat-minlat)/dx)*dx+minlat+2*dx)
       cornerlon2 = min(maxlon, nint((LIS_domain(n)%maxlon-minlon)/dy)*dy+minlon+2*dy)

       ANSASWEsnow_struc(n)%offset1 = nint((cornerlon1-minlon)/dy)
       ANSASWEsnow_struc(n)%offset2 = nint((cornerlat1-minlat)/dx)

       ANSASWEsnow_struc(n)%nr = nint((cornerlat2-cornerlat1)/dx)+1
       ANSASWEsnow_struc(n)%nc = nint((cornerlon2-cornerlon1)/dy)+1

       gridDesci(n,1) = 0 
       gridDesci(n,2) = ANSASWEsnow_struc(n)%nc
       gridDesci(n,3) = ANSASWEsnow_struc(n)%nr
       gridDesci(n,4) = cornerlat1
       gridDesci(n,5) = cornerlon1
       gridDesci(n,6) = 128
       gridDesci(n,7) = cornerlat2
       gridDesci(n,8) = cornerlon2
       gridDesci(n,9) = ANSASWEsnow_struc(n)%gridDesc(5)
       gridDesci(n,10) = ANSASWEsnow_struc(n)%gridDesc(6)
       gridDesci(n,20) = 64.0

       ANSASWEsnow_struc(n)%mi = ANSASWEsnow_struc(n)%nc*ANSASWEsnow_struc(n)%nr
       allocate(ANSASWEsnow_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASWEsnow_struc(n)%n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASWEsnow_struc(n)%n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASWEsnow_struc(n)%n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(ANSASWEsnow_struc(n)%w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASWEsnow_struc(n)%w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASWEsnow_struc(n)%w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASWEsnow_struc(n)%w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(ANSASWEsnow_struc(n)%swe(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       ANSASWEsnow_struc(n)%swe = LIS_rc%udef

       call bilinear_interp_input(n, gridDesci(n,:),& 
            ANSASWEsnow_struc(n)%n11,ANSASWEsnow_struc(n)%n12,&
            ANSASWEsnow_struc(n)%n21,ANSASWEsnow_struc(n)%n22,&
            ANSASWEsnow_struc(n)%w11,ANSASWEsnow_struc(n)%w12,&
            ANSASWEsnow_struc(n)%w21,ANSASWEsnow_struc(n)%w22)
            
       call LIS_registerAlarm("ANSA SWE read alarm", 86400.0,86400.0)
    enddo

    write(LIS_logunit,*) 'Created ESMF States to hold ANSA observations data'

  end subroutine ANSASWEsnow_setup
  
end module ANSASWEsnow_Mod

