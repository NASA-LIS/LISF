!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module PALSmetdata_forcingMod
!BOP
! !MODULE: PALSmetdata_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of station data from 
!  PALS (Protocol for the Analysis of Land Surface Models; 
!  http://www.pals.unsw.edu.au) used as forcing
!  within LIS. No spatial interpolation of the data is performed
!  as it is expected to be point data. The LIS domain must be configured
!  to exactly match the station location. 
!
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_PALSmetdata      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: PALSmetdata_struc
!EOP

  type, public ::  PALSmetdata_type_dec 

     real                       :: ts
     character(len=LIS_CONST_PATH_LEN) :: PALSmetdatadir     
     character*80               :: stn_name
     real*8                     :: fcsttime1,fcsttime2
     integer                    :: findtime1, findtime2
     integer                    :: yr,mo
     integer                    :: tindex1,tindex2
     integer                    :: syr, smo, sda, shr, smn, sss
     real,       allocatable        :: tair1(:)
     real,       allocatable        :: qair1(:)
     real,       allocatable        :: swdown1(:)
     real,       allocatable        :: lwdown1(:)
     real,       allocatable        :: wind1(:)
     real,       allocatable        :: psurf1(:)
     real,       allocatable        :: rainf1(:)

     real,       allocatable        :: tair2(:)
     real,       allocatable        :: qair2(:)
     real,       allocatable        :: swdown2(:)
     real,       allocatable        :: lwdown2(:)
     real,       allocatable        :: wind2(:)
     real,       allocatable        :: psurf2(:)
     real,       allocatable        :: rainf2(:)
     logical                    :: startFlag
     type(ESMF_Time)            :: reftime
     type(ESMF_TimeInterval)    :: timestep

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type PALSmetdata_type_dec

  type(PALSmetdata_type_dec), allocatable :: PALSmetdata_struc(:)

!EOP
contains
  
!BOP
!
! !ROUTINE: init_PALSmetdata
! \label{init_PALSmetdata}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 
! !INTERFACE:
  subroutine init_PALSmetdata(findex)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native domain of the input forcing for PALS 
!  data (which is assumed to be at the point scale)
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_PALSmetdata](\ref{readcrd_PALSmetdata}) \newline
!     reads the runtime options specified for PALS station data
!  \end{description}
!EOP
    
    integer :: n
    integer :: status

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the PALS met forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    LIS_rc%met_nf(findex) = 10
    
    allocate(PALSmetdata_struc(LIS_rc%nnest))
    call readcrd_PALSmetdata()
    
    do n=1, LIS_rc%nnest
       PALSmetdata_struc(n)%ts = 1800
       call LIS_update_timestep(LIS_rc, n, PALSmetdata_struc(n)%ts)
    enddo
    
    do n=1,LIS_rc%nnest

       allocate(PALSmetdata_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(PALSmetdata_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       PALSmetdata_struc(n)%metdata1 = 0
       PALSmetdata_struc(n)%metdata2 = 0

       PALSmetdata_struc(n)%findtime1 = 0 
       PALSmetdata_struc(n)%findtime2 = 0 
       
       call ESMF_TimeSet(PALSmetdata_struc(n)%reftime, &
            yy=PALSmetdata_struc(n)%syr, &
            mm=PALSmetdata_struc(n)%smo, &
            dd=PALSmetdata_struc(n)%sda, &
            h =PALSmetdata_struc(n)%shr, &
            m =PALSmetdata_struc(n)%smn, &
            s =PALSmetdata_struc(n)%sss,&
            calendar=LIS_calendar, rc=status)
       call LIS_verify(status, 'error in ESMF_TimeSet in get_PALSmetdata')

       call ESMF_TimeIntervalSet(PALSmetdata_struc(n)%timestep,&
            s=1800,rc=status)
       call LIS_verify(status, 'Error in timeintervalset in PALSmetdata')

       allocate(PALSmetdata_struc(n)%tair1(17568)) !366*24*2
       allocate(PALSmetdata_struc(n)%qair1(17568)) 
       allocate(PALSmetdata_struc(n)%swdown1(17568))
       allocate(PALSmetdata_struc(n)%lwdown1(17568))
       allocate(PALSmetdata_struc(n)%wind1(17568)) 
       allocate(PALSmetdata_struc(n)%psurf1(17568))
       allocate(PALSmetdata_struc(n)%rainf1(17568)) 

       allocate(PALSmetdata_struc(n)%tair2(17568)) 
       allocate(PALSmetdata_struc(n)%qair2(17568)) 
       allocate(PALSmetdata_struc(n)%swdown2(17568))
       allocate(PALSmetdata_struc(n)%lwdown2(17568))
       allocate(PALSmetdata_struc(n)%wind2(17568)) 
       allocate(PALSmetdata_struc(n)%psurf2(17568))
       allocate(PALSmetdata_struc(n)%rainf2(17568)) 

       PALSmetdata_struc(n)%tair1   = LIS_rc%udef
       PALSmetdata_struc(n)%qair1   = LIS_rc%udef
       PALSmetdata_struc(n)%swdown1 = LIS_rc%udef
       PALSmetdata_struc(n)%lwdown1 = LIS_rc%udef
       PALSmetdata_struc(n)%wind1   = LIS_rc%udef
       PALSmetdata_struc(n)%psurf1  = LIS_rc%udef
       PALSmetdata_struc(n)%rainf1  = LIS_rc%udef

       PALSmetdata_struc(n)%tair2   = LIS_rc%udef
       PALSmetdata_struc(n)%qair2   = LIS_rc%udef
       PALSmetdata_struc(n)%swdown2 = LIS_rc%udef
       PALSmetdata_struc(n)%lwdown2 = LIS_rc%udef
       PALSmetdata_struc(n)%wind2   = LIS_rc%udef
       PALSmetdata_struc(n)%psurf2  = LIS_rc%udef
       PALSmetdata_struc(n)%rainf2  = LIS_rc%udef
    enddo
  end subroutine init_PALSmetdata

end module PALSmetdata_forcingMod
