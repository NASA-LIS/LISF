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
! !MODULE: FLUXNETdata_module
! 
! !DESCRIPTION: 
!  
!   
! !REVISION HISTORY: 
!  29 Mar 11    Sujay Kumar;   Initial Specification
! 
module FLUXNETdata_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: FLUXNETdata_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: FLUXNETdata_struc

  type, public ::  FLUXNETdata_data_dec

     character(len=LIS_CONST_PATH_LEN) :: odir
     integer                 :: nc,nr
     integer, allocatable        :: n11(:)
     real,    allocatable        :: qle(:,:)
     real,    allocatable        :: qh(:,:)
     integer                 :: yr,mo

  end type FLUXNETdata_data_dec

  type(FLUXNETdata_data_dec), allocatable :: FLUXNETdata_struc(:)

contains
!BOP
! 
! !ROUTINE: FLUXNETdata_setup
! \label{FLUXNETdata_setup}
! 
! !INTERFACE: 
  subroutine FLUXNETdata_setup(Obj_space)
! !USES: 
    use LIS_coreMod
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obj_space(LIS_rc%nnest)
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
    integer                   ::  i
    integer                   ::  ftn
    character*100, allocatable    ::  vname(:)
    integer                   ::  status
    integer                   ::  nobstypes
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field),allocatable  ::  obsField(:)
    real                      ::  gridDesci(50)
    character(len=LIS_CONST_PATH_LEN) ::  objspaceAttribFile(LIS_rc%nnest)

    allocate(FLUXNETdata_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"FLUXNET data directory:",&
         rc=status)

    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,FLUXNETdata_struc(n)%odir,&
            rc=status)
       call LIS_verify(status, "FLUXNET data directory: not defined")
       
       call ESMF_AttributeSet(Obj_Space(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "FLUXNET objective space attributes file:",rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,objspaceAttribFile(n),rc=status)
       call LIS_verify(status,"FLUXNET objective space attributes file: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"FLUXNET number of observation types:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,nobstypes,&
            rc=status)
       call LIS_verify(status, "FLUXNET number of observation types: not defined")
   
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

       FLUXNETdata_struc(n)%nc = 720
       FLUXNETdata_struc(n)%nr = 291

       FLUXNETdata_struc(n)%yr = -1
       FLUXNETdata_struc(n)%mo = LIS_rc%mo

       gridDesci = 0 
       gridDesci(1) = 0  
       gridDesci(2) = 720
       gridDesci(3) = 291
       gridDesci(4) = -55.25
       gridDesci(5) = -179.75
       gridDesci(7) = 89.75
       gridDesci(8) = 179.75
       gridDesci(6) = 128
       gridDesci(9) = 0.50
       gridDesci(10) = 0.50
       gridDesci(20) = 64

       allocate(FLUXNETdata_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(FLUXNETdata_struc(n)%qle(LIS_rc%lnc(n)*LIS_rc%lnr(n),12))
       allocate(FLUXNETdata_struc(n)%qh(LIS_rc%lnc(n)*LIS_rc%lnr(n),12))

       call neighbor_interp_input(n,gridDesci,&
            FLUXNETdata_struc(n)%n11)       
    enddo
    write(LIS_logunit,*) 'Created the States to hold FLUXNET observations'
    
  end subroutine FLUXNETdata_setup
  
end module FLUXNETdata_module
