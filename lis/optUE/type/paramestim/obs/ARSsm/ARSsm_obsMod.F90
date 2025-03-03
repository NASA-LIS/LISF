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
! !MODULE: ARSsm_obsMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  26 Jan 2018    Soni Yatheendradas;   Initial Specification based on LVT ARSsm datastream
! 
module ARSsm_obsMod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ARSsm_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ARSsm_obs_struc

  type, public ::  ARSsm_obs_data_dec

     character(len=LIS_CONST_PATH_LEN) :: odir
     integer               :: yr
     integer               :: n_stns
     character*50, allocatable :: stn_name(:)
     real,         allocatable :: stn_lat(:)
     real,         allocatable :: stn_lon(:)
     integer,      allocatable :: stn_col(:)
     integer,      allocatable :: stn_row(:)
     type(ESMF_Time)       :: startTime
     type(ESMF_TimeInterval), allocatable :: timestep(:)
     real,         allocatable :: sm(:,:)
     real,         allocatable :: sm_std(:,:)

  end type ARSsm_obs_data_dec

  type(ARSsm_obs_data_dec), allocatable :: ARSsm_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: ARSsm_obs_setup
! \label{ARSsm_obs_setup}
! 
! !INTERFACE: 
  subroutine ARSsm_obs_setup(Obs_State)
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
    real                  :: col,row
    character(len=LIS_CONST_PATH_LEN) :: stnlist_file

    allocate(ARSsm_obs_struc(LIS_rc%nnest))

    write(LIS_logunit,*) '[INFO] Setting up ARSsm data reader....'

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"ARS soil moisture data directory:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,obsdir,&
         rc=status)
    call LIS_verify(status, 'ARS soil moisture data directory: not defined')
    do n=1,LIS_rc%nnest

       ARSsm_obs_struc(n)%odir = obsdir  

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config, "ARS soil moisture station list file:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,stnlist_file,&
         rc=status)
    call LIS_verify(status, 'ARS soil moisture station list file: not defined')

    do n=1,LIS_rc%nnest

       ftn = LIS_getNextUnitNumber()
       write(LIS_logunit,*) 'Reading ARS soil moisture station list file ',&
            trim(stnlist_file)
       open(ftn, file=trim(stnlist_file), form='formatted')
       read(ftn,*)
       read(ftn,*) ARSsm_obs_struc(n)%n_stns
       read(ftn,*) 
       
       allocate(ARSsm_obs_struc(n)%stn_name(ARSsm_obs_struc(n)%n_stns))
       allocate(ARSsm_obs_struc(n)%stn_lat(ARSsm_obs_struc(n)%n_stns))
       allocate(ARSsm_obs_struc(n)%stn_lon(ARSsm_obs_struc(n)%n_stns))

       allocate(ARSsm_obs_struc(n)%stn_col(ARSsm_obs_struc(n)%n_stns))
       allocate(ARSsm_obs_struc(n)%stn_row(ARSsm_obs_struc(n)%n_stns))

       allocate(ARSsm_obs_struc(n)%timestep(ARSsm_obs_struc(n)%n_stns))

       write(LIS_logunit,*) '[INFO] ARSsm station list .. '
       do k=1,ARSsm_obs_struc(n)%n_stns
          read(ftn,*) ARSsm_obs_struc(n)%stn_name(k), ARSsm_obs_struc(n)%stn_lat(k), &
               ARSsm_obs_struc(n)%stn_lon(k)
          write(LIS_logunit,*) '[INFO] ',trim(ARSsm_obs_struc(n)%stn_name(k)), &
               ARSsm_obs_struc(n)%stn_lat(k), ARSsm_obs_struc(n)%stn_lon(k)
          call latlon_to_ij(LIS_domain(n)%lisproj, ARSsm_obs_struc(n)%stn_lat(k), &
               ARSsm_obs_struc(n)%stn_lon(k), col, row)
          ARSsm_obs_struc(n)%stn_col(k) = nint(col)
          ARSsm_obs_struc(n)%stn_row(k) = nint(row)
!          write(LIS_logunit,*) '[INFO] stn col and row are:', ARSsm_obs_struc(n)%stn_col(k),',',&
!                               ARSsm_obs_struc(n)%stn_row(k)
        

          if ( (ARSsm_obs_struc(n)%stn_name(k) .eq. 'Walnut_Gulch') .or. &
               (ARSsm_obs_struc(n)%stn_name(k) .eq. 'Little_River') .or. &
               (ARSsm_obs_struc(n)%stn_name(k) .eq. 'Little_Washita') .or. &
               (ARSsm_obs_struc(n)%stn_name(k) .eq. 'Fort_Cobb') ) then
              call ESMF_TimeIntervalSet(ARSsm_obs_struc(n)%timestep(k),s=1800,rc=status) ! SY: because different watersheds have different time interval files
              call LIS_verify(status,"ESMF_TimeIntervalSet failed in ARSsm_obs_setup")
              call LIS_update_timestep(LIS_rc, n, 1800.0)
          elseif ( (ARSsm_obs_struc(n)%stn_name(k) .eq. 'Reynolds_Creek') .or. &
               (ARSsm_obs_struc(n)%stn_name(k) .eq. 'South_Fork') .or. &
               (ARSsm_obs_struc(n)%stn_name(k) .eq. 'StJoseph') ) then
              call ESMF_TimeIntervalSet(ARSsm_obs_struc(n)%timestep(k),s=3600,rc=status) ! SY: because different watersheds have different time interval files
              call LIS_verify(status,"ESMF_TimeIntervalSet failed in ARSsm_obs_setup")
              call LIS_update_timestep(LIS_rc, n, 3600.0)
          endif 
       enddo !  do k=1,ARSsm_obs_struc(n)%n_stns

       call LIS_releaseUnitNumber(ftn)

       ARSsm_obs_struc(n)%yr = -1
   
       allocate(ARSsm_obs_struc(n)%sm(ARSsm_obs_struc(n)%n_stns, 366*48)) ! SY: 366*48 30-min timesteps for leap yr is maximum possible 30-min timesteps for any year during LIS simulation; Both 30-min and 1-hr files present, so taking the smallest timestep of 30-min 
       allocate(ARSsm_obs_struc(n)%sm_std(ARSsm_obs_struc(n)%n_stns, 366*48)) ! SY: 366*48 30-min timesteps for leap yr is maximum possible 30-min timesteps for any  year during LIS simulation; Both 30-min and 1-hr files present, so taking the smallest timestep of 30-min


       ARSsm_obs_struc(n)%sm = LIS_rc%udef
       ARSsm_obs_struc(n)%sm_std = LIS_rc%udef

    enddo !  do n=1,LIS_rc%nnest
    
    call ESMF_ConfigFindLabel(LIS_config,"ARS soil moisture data attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'ARS soil moisture data attributes file: not defined')
   
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

    write(LIS_logunit,*) 'Created the States to hold the ARS soil moisture data'
    
  end subroutine ARSsm_obs_setup
  
end module ARSsm_obsMod
