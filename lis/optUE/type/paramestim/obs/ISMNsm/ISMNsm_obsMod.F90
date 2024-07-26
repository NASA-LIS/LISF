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
! !MODULE: ISMNsm_obsMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  21 Sep 2018 Sujay Kumar;   Initial Specification 
! 
module ISMNsm_obsMod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ISMNsm_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ISMNsm_obs_struc

  type, public :: ismnstn
     integer                    :: vlevels
     real                       :: lat,lon
     character(len=LIS_CONST_PATH_LEN), allocatable :: fname(:)
     real,         allocatable  :: sm(:,:)
     real,         allocatable  :: sfsm(:)
     real,         allocatable  :: rzsm(:)

  end type ismnstn

  type, public ::  ISMNsm_obs_data_dec

     character(len=LIS_CONST_PATH_LEN) :: odir
     integer                     :: yr
     integer                     :: n_stns
     integer                     :: nts 
     type(ESMF_Time)             :: startTime
     type(ESMF_TimeInterval)     :: timestep
     real                        :: lis_sf_d
     real                        :: lis_rz_d
     type(ismnstn), allocatable  :: stn(:)

  end type ISMNsm_obs_data_dec

  type(ISMNsm_obs_data_dec), allocatable :: ISMNsm_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: ISMNsm_obs_setup
! \label{ISMNsm_obs_setup}
! 
! !INTERFACE: 
  subroutine ISMNsm_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, &
         LIS_vecGrid, LIS_domain
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use map_utils, only : latlon_to_ij 
    use LIS_timeMgrMod, only : LIS_update_timestep

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   The arguments are: 
!   \begin{description}
!    \item[Obs\_State]   observation state object 
!   \end{description}
!EOP
    integer                   ::  n 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField
    character(len=LIS_CONST_PATH_LEN) ::  obsdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                 :: k
    integer                 :: ftn
    integer                 :: status


    allocate(ISMNsm_obs_struc(LIS_rc%nnest))

    write(LIS_logunit,*) '[INFO] Setting up ISMN sm data reader....'

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"ISMN soil moisture data directory:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,obsdir,&
         rc=status)
    call LIS_verify(status, 'ISMN soil moisture data directory: not defined')

    do n=1,LIS_rc%nnest

       ISMNsm_obs_struc(n)%odir = obsdir  

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    do n=1,LIS_rc%nnest
       ISMNsm_obs_struc(n)%yr = -1
       
       call ESMF_TimeIntervalSet(ISMNsm_obs_struc(n)%timestep,s=3600,rc=status)
       call LIS_verify(status,"ESMF_TimeIntervalSet failed in ISMN_obsInit")

       call LIS_update_timestep(LIS_rc, n, 3600.0)

       ISMNsm_obs_struc(n)%nts =8784 !24*366
    enddo 
    
    call ESMF_ConfigFindLabel(LIS_config,"ISMN soil moisture target surface soil layer thickness:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ISMNsm_obs_struc(n)%lis_sf_d,rc=status)
       call LIS_verify(status, 'ISMN soil moisture target surface soil layer thickness: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISMN soil moisture target root zone soil layer thickness:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ISMNsm_obs_struc(n)%lis_rz_d,rc=status)
       call LIS_verify(status, 'ISMN soil moisture target root zone soil layer thickness: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ISMN soil moisture data attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'ISMN soil moisture data attributes file: not defined')
   
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
    enddo !  do n=1,LIS_rc%nnest

    write(LIS_logunit,*) '[INFO] created the States to hold the ISMN soil moisture data'
    
  end subroutine ISMNsm_obs_setup
  
end module ISMNsm_obsMod
