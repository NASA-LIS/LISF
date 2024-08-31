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
! !MODULE: GLERL_dataMod
! \label(GLERL_dataMod)
!
! !INTERFACE:
module GLERL_dataMod
! 
! !USES: 
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 May 2011   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLERL_obsinit !Initializes structures for reading GLERL data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GLERLobs !Object to hold GLERL observation attributes
!EOP

  type, public :: glerldec
     character*100               :: odir, locfile
     integer                     :: nlocs
     character*50,  allocatable  :: lake_locname(:)
     real,          allocatable  :: lake_lat(:)
     real,          allocatable  :: lake_lon(:)
     real,    allocatable        :: qle(:,:)
     real,    allocatable        :: watert(:,:)
     integer                     :: yr
     integer                     :: mo
  end type glerldec
     
  type(glerldec), allocatable :: GLERLobs(:)

contains
  
!BOP
! 
! !ROUTINE: GLERL_obsInit
! \label{GLERL_obsInit}
!
! !INTERFACE: 
  subroutine GLERL_obsinit(i)
! 
! !USES: 
    use LVT_coreMod,   only : LVT_rc, LVT_Config
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)   :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the GLERL data, including the computation of spatial 
!   interpolation weights. The GLERL data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer               :: ftn
    integer               :: status
    integer               :: k

    if(.not.allocated(GLERLobs)) then 
       allocate(GLERLobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, GLERLobs(i)%odir, &
         label='GLERL hydro data directory: ',rc=status)
    call LVT_verify(status, 'GLERL hydro data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, GLERLobs(i)%locfile, &
         label="GLERL hydro lake locations file: ",rc=status)
    call LVT_verify(status, "GLERL hydro lake locations file: not defined")

    call LVT_update_timestep(LVT_rc, 2592000)

    ftn = LVT_getNextUnitNumber()
    open(ftn, file=GLERLobs(i)%locfile,form='formatted')
    read(ftn,*)
    read(ftn,*) GLERLobs(i)%nlocs
    read(ftn,*)

    allocate(GLERLobs(i)%lake_locname(GLERLobs(i)%nlocs))
    allocate(GLERLobs(i)%lake_lat(GLERLobs(i)%nlocs))
    allocate(GLERLobs(i)%lake_lon(GLERLobs(i)%nlocs))
    
    write(LVT_logunit,*) '[INFO] Using GLERL hydro data ...'
    do k=1,GLERLobs(i)%nlocs
       read(ftn,*) GLERLobs(i)%lake_locname(k)
       read(ftn,*) GLERLobs(i)%lake_lat(k), GLERLobs(i)%lake_lon(k)
       write(LVT_logunit,*) '[INFO] ',trim(GLERLobs(i)%lake_locname(k)), &
            GLERLobs(i)%lake_lat(k), GLERLobs(i)%lake_lon(k)
    enddo
    
    call LVT_releaseUnitNumber(ftn)

    allocate(GLERLobs(i)%qle(GLERLobs(i)%nlocs,12))
    allocate(GLERLobs(i)%watert(GLERLobs(i)%nlocs,12))

    GLERLobs(i)%yr = -1
    GLERLobs(i)%mo = LVT_rc%mo

    if(LVT_rc%tavgInterval.lt.2592000) then 
       write(LVT_logunit,*) '[ERR] The time averaging interval must be greater than'
       write(LVT_logunit,*) '[ERR] equal to a month since the GLERL data is monthly'
       call LVT_endrun()
    endif
  end subroutine GLERL_obsinit


end module GLERL_dataMod
