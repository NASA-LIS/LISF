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
! !MODULE: wgPBMRsmobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle the Walnut Gulch (WG) PBMR soil moisture data
!   
! !REVISION HISTORY: 
!  09 Jul 09    Sujay Kumar;   Initial Specification
! 
module wgPBMRsmobs_module
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: wgPBMRsmdata_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: wgPBMRsm_struc

  type, public ::  wgPBMRsm_data_dec
     integer             :: site
  end type wgPBMRsm_data_dec

  type(wgPBMRsm_data_dec), allocatable :: wgPBMRsm_struc(:)

contains
!BOP
! 
! !ROUTINE: wgPBMRsmdata_setup
! \label{wgPBMRsmdata_setup}
! 
! !INTERFACE: 
  subroutine wgPBMRsmdata_setup(Obs_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for reading the PBMR
!   soil moisture data for the Walnut Gulch domain. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[Obj\_Space]   observation/Objective space object 
!   \end{description}
!EOP
    integer                   ::  n 
    integer                   ::  status
    integer                   ::  i 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField
    character(len=LIS_CONST_PATH_LEN) ::  smobsdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    real                      ::  gridDesci(LIS_rc%nnest,50)

    allocate(wgPBMRsm_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"WG PBMR soil moisture data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,smobsdir,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            smobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WG PBMR observations attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"WG PBMR site index:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wgPBMRsm_struc(n)%site,&
            rc=status)
       call LIS_verify(status)
   
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(obsAttribFile(n)),status='old')
       read(ftn,*)
       
       read(ftn,fmt='(a100)') vname
       obsField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecGrid(n), &
            name=trim(vname), rc=status)
       call LIS_verify(status)
       
       call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
       call LIS_verify(status)

       call LIS_releaseUnitNumber(ftn)
    enddo


    write(LIS_logunit,*) 'Created the States to hold the WG PBMR observations'
    
  end subroutine wgPBMRsmdata_setup
  
end module WgPBMRsmobs_module
