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
! !MODULE: ISCCP_Tskinobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle synthetic ISCCP Tskin data
!   
! !REVISION HISTORY: 
!  02Feb08    Sujay Kumar;   Initial Specification
! 
module ISCCP_Tskinobs_module
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ISCCP_Tskinobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: isccp_tskin_struc

  type, public ::  isccp_ts_dec
     integer             :: mi
     integer, allocatable    :: n11(:)
  end type isccp_ts_dec

  type(isccp_ts_dec), allocatable :: isccp_tskin_struc(:)

contains
!BOP
! 
! !ROUTINE: ISCCP_Tskinobs_setup
! \label{ISCCP_Tskinobs_setup}
! 
! !INTERFACE: 
  subroutine ISCCP_Tskinobs_setup(Obj_Space)
! !USES: 
    use LIS_coreMod
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obj_Space(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for soil moisture
!    assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                   ::  n 
    integer                   ::  status
    integer                   ::  nobstypes
    integer                   ::  i 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field), allocatable ::  obsField(:)
    character(len=LIS_CONST_PATH_LEN) ::  synsmobsdir
    character*100,  allocatable   ::  vname(:)
    character*100             ::  objspaceAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    real                      ::  gridDesci(LIS_rc%nnest,50)

    allocate(isccp_tskin_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synsmobsdir,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obj_Space(n),"Data Directory",&
            synsmobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obj_Space(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin objective space attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,objspaceAttribFile(n),rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISCCP Tskin number of observation types:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,nobstypes,&
            rc=status)
       call LIS_verify(status)
   
       allocate(obsField(nobstypes))
       allocate(vname(nobstypes))

       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(objspaceAttribFile(n)),status='old')
       read(ftn,*)
       
       do i=1,nobstypes
          read(ftn,fmt='(a100)') vname(i)
          obsField(i) = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecGrid(n), &
               name=trim(vname(i)), rc=status)
          call LIS_verify(status)

          call ESMF_StateAdd(Obj_Space(n),(/obsField(i)/),rc=status)
          call LIS_verify(status)
       enddo
       call LIS_releaseUnitNumber(ftn)
       deallocate(vname)
    enddo

!-------------------------------------------------------------
! setting up interpolation weights for neighbor search
!-------------------------------------------------------------
    gridDesci = 0 

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0 
       gridDesci(n,2) = 360
       gridDesci(n,3) = 180
       gridDesci(n,4) = -89.50
       gridDesci(n,5) = -179.50
       gridDesci(n,6) = 128
       gridDesci(n,7) = 89.50
       gridDesci(n,8) = 179.50
       gridDesci(n,9) = 1.000
       gridDesci(n,10) = 1.000
       gridDesci(n,20) = 64.0

       isccp_tskin_struc(n)%mi = 360*180
       allocate(isccp_tskin_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       call neighbor_interp_input(gridDesci(n,:),&
            LIS_rc%lnc(n), LIS_rc%lnr(n),&
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            isccp_tskin_struc(n)%n11)
    enddo

    write(LIS_logunit,*) 'Created the States to hold the observations data'
    
  end subroutine ISCCP_Tskinobs_setup
  
end module ISCCP_Tskinobs_module
