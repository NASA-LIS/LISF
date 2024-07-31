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
! !MODULE: GlobalLSDataMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  09 Jul 09    Sujay Kumar;   Initial Specification
! 
module GlobalLSDataMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: globallsobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: globallsobs_struc

  type, public :: landslide_data_dec
     integer             :: size
     real, allocatable       :: lat(:)
     real, allocatable       :: lon(:)
     logical, allocatable    :: flag(:)
     type(ESMF_Time), allocatable       :: time(:)
  end type landslide_data_dec

  type(landslide_data_dec), allocatable :: globallsobs_struc(:)

contains
!BOP
! 
! !ROUTINE: globallsobs_setup
! \label{globallsobs_setup}
! 
! !INTERFACE: 
  subroutine globallsobs_setup(Obs_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid
    use LIS_timeMgrMod, only : LIS_calendar
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
    integer                   ::  size, ios, id, yr, mo, da, hr, mn
    real                      ::  lat, lon
    type(ESMF_Time)           ::  time1
    real                      ::  gridDesci(LIS_rc%nnest,50)

    allocate(globallsobs_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"Global Landslide Obs data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,landslideobsdir,&
            rc=status)
       call LIS_verify(status, 'Err: Global Landslide Obs data directory: not defined')

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            landslideobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Global Landslide observations attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'Err: Global Landslide observations attributes file: not defined')
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

!first read to figure out the size of entries
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(landslideobsdir),status='old')

       ios = 0 
       size = 0 
       do while(ios.eq.0)
          read(ftn,*,iostat=ios) id, yr, mo, da, hr, mn, lat, lon
          size = size + 1
       enddo

       call LIS_releaseUnitNumber(ftn)
       
       globallsobs_struc(n)%size = size
       allocate(globallsobs_struc(n)%time(size))
       allocate(globallsobs_struc(n)%lat(size))
       allocate(globallsobs_struc(n)%lon(size))
       allocate(globallsobs_struc(n)%flag(size))

       globallsobs_struc(n)%flag = .false. 

!read to store the values
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(landslideobsdir),status='old')

       ios = 0 
       size = 0 
       do while(ios.eq.0)
          read(ftn,*,iostat=ios) id, yr, mo, da, hr, mn, lat, lon
          size = size + 1
          if(hr.le.0) then 
             hr = 0 
             mn = 0 
          endif
          
          call ESMF_TimeSet(time1, yy=yr, mm=mo, dd=da,h=hr,&
               m=mn,calendar=LIS_calendar,rc=status)

          globallsobs_struc(n)%time(size) = time1 
          globallsobs_struc(n)%lat(size) = lat
          globallsobs_struc(n)%lon(size) = lon
       enddo

       call LIS_releaseUnitNumber(ftn)
    enddo

    write(LIS_logunit,*) 'Created the States to hold the landslide observations'
    
  end subroutine globallsobs_setup
  
end module GlobalLSDataMod
