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
! !MODULE: MaconLSDataMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  09 Jul 09    Sujay Kumar;   Initial Specification
! 
module MaconLSDataMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: maconlsobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: maconlsobs_struc

  type, public :: maconls_data_dec
     integer             :: site
  end type maconls_data_dec

  type(maconls_data_dec), allocatable :: maconlsobs_struc(:)

contains
!BOP
! 
! !ROUTINE: maconlsobs_setup
! \label{maconlsobs_setup}
! 
! !INTERFACE: 
  subroutine maconlsobs_setup(Obs_State)
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
    character(len=LIS_CONST_PATH_LEN) ::  landslideobsdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    real                      ::  gridDesci(LIS_rc%nnest,50)

    allocate(maconlsobs_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"Macon County Landslide Obs data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,landslideobsdir,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            landslideobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Macon County Landslide observations attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status)
    enddo

    do n=1,LIS_rc%nnest
   
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(obsAttribFile(n)),status='old')
       read(ftn,*)
       
       read(ftn,fmt='(a100)') vname
       obsField = ESMF_FieldCreate(arrayspec=realarrspec, &
            grid=LIS_vecGrid(n), &
            name=trim(vname), rc=status)
       call LIS_verify(status)
       
       call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
       call LIS_verify(status)

       call LIS_releaseUnitNumber(ftn)
    enddo


    write(LIS_logunit,*) 'Created the States to hold the landslide observations'
    
  end subroutine maconlsobs_setup
  
end module MaconLSDataMod
