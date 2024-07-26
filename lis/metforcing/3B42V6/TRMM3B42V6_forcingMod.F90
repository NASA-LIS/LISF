!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module TRMM3B42V6_forcingMod
!BOP
! !MODULE: TRMM3B42V6_forcingMod
!
! !REVISION HISTORY:
! 21 Jun 2013: Soni Yatheendradas; changes from earlier code to avoid
!              alternate file skip, jump to previous day TRMM and
!              absence of rain rate weighting
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  NASA TRMM 3B42 (Version 6) merged analysis of precipitation (TRMM 3B42V6). 
!  TRMM 3B42 merges satellite passive microwave estimates including TMI, 
!  SSM/I, AMSR and AMSU to produce a quasi-global 0.25 degree 3-hourly
!  precipitation analysis. 
!
!  The implementation in LIS has the derived data type {\tt TRMM3B42V6\_struc}
!  that includes the variables to specify the runtime options, and
!  the weights and neighbor information for spatial interpolation
!
!  They are desribed below:
! \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[TRMM3B42V6dir]
!    Directory containing the input data
!  \item[TRMM3B42V6time]
!    The nearest instance of the incoming
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES:
    USE LIS_constantsMod, only: LIS_CONST_PATH_LEN

    implicit none
    PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_TRMM3B42V6     !defines the native resolution of
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: TRMM3B42V6_struc

!EOP

  type, public :: TRMM3B42V6_type_dec
     real                     :: ts 
     integer                  :: ncold
     integer                  :: nrold
     character(len=LIS_CONST_PATH_LEN) :: TRMM3B42V6dir
     real*8                   :: TRMM3B42V6time_TStepStart 
     integer                  :: TRMM3B42V6yr_TStepStart ! SY
     integer                  :: TRMM3B42V6mo_TStepStart ! SY
     integer                  :: TRMM3B42V6da_TStepStart ! SY
     integer                  :: TRMM3B42V6hr_TStepStart ! SY
     integer                  :: TRMM3B42V6yr_TStepStart_Previous ! SY
     integer                  :: TRMM3B42V6mo_TStepStart_Previous ! SY
     integer                  :: TRMM3B42V6da_TStepStart_Previous ! SY
     integer                  :: TRMM3B42V6hr_TStepStart_Previous ! SY
     real*8                   :: TRMM3B42V6time_TStepEnd ! SY
     integer                  :: TRMM3B42V6yr_TStepEnd ! SY
     integer                  :: TRMM3B42V6mo_TStepEnd ! SY
     integer                  :: TRMM3B42V6da_TStepEnd ! SY
     integer                  :: TRMM3B42V6hr_TStepEnd ! SY
     integer                  :: TRMM3B42V6yr_TStepEnd_Previous ! SY
     integer                  :: TRMM3B42V6mo_TStepEnd_Previous ! SY
     integer                  :: TRMM3B42V6da_TStepEnd_Previous ! SY
     integer                  :: TRMM3B42V6hr_TStepEnd_Previous ! SY
     real*8                   :: LIS_timeAtTStepStart 
     real*8                   :: LIS_timeAtTStepEnd 
     integer                  :: mi
     integer, allocatable     :: n112(:,:)
     integer, allocatable     :: n122(:,:)
     integer, allocatable     :: n212(:,:)
     integer, allocatable     :: n222(:,:)
     real,    allocatable     :: w112(:,:),w122(:,:)
     real,    allocatable     :: w212(:,:),w222(:,:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type TRMM3B42V6_type_dec

  type(TRMM3B42V6_type_dec), allocatable :: TRMM3B42V6_struc(:)
contains

!BOP
!
! !ROUTINE: init_TRMM3B42V6
! \label{init_TRMM3B42V6}
!
! !REVISION HISTORY:
! 11Dec2003: Sujay Kumar; Initial Specification 
! 25Aug2006: Yudong Tian; Implementation for TRMM 3B42 V6
!
! !INTERFACE:
  subroutine init_TRMM3B42V6(findex)
! !USES:
    use LIS_coreMod, only: LIS_rc, LIS_config, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep, &
                               LIS_tick, LIS_parseTimeString, &
                               LIS_registerAlarm ! SY
    use LIS_logMod,    only : LIS_logunit, LIS_endrun, &
                              LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use LIS_FORC_AttributesMod
    use ESMF

    implicit none
    
    integer,  intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for TRMM 3B42V6
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_TRMM3B42V6](\ref{readcrd_TRMM3B42V6}) \newline
!     reads the runtime options specified for TRMM 3B42V6 data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    real    :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: n
    integer :: ts1,ts2 
    integer :: LIS_syr,LIS_smo,LIS_sda,LIS_shr,LIS_smn,LIS_sss 

    real*8  :: DataAvailabilityStartTime 
    real*8  :: LIS_StartTime 

    integer :: rc 
    character (len=10):: time 

    allocate(TRMM3B42V6_struc(LIS_rc%nnest))

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the TRMM-3B42 version 6, forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode. May be'
       write(LIS_logunit,*) '[ERR]  added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    ! SY: Begin LIS time check against TRMM data availability time
    yr1 = 1998
    mo1 = 01
    da1 = 01
    hr1 = 0
    mn1 = 0; ss1 = 0
    ts1 = (-1) * 90 * 60 ! For 90 min behind 1st available TRMM data
    call LIS_tick(DataAvailabilityStartTime,updoy,upgmt,yr1,mo1,&
         da1,hr1,mn1,ss1,real(ts1) )
    LIS_syr = LIS_rc%syr
    LIS_smo = LIS_rc%smo
    LIS_sda = LIS_rc%sda
    LIS_shr = LIS_rc%shr
    LIS_smn = LIS_rc%smn
    LIS_sss = LIS_rc%sss
    ts2 = 0
    call LIS_tick(LIS_StartTime,updoy,upgmt,LIS_syr,LIS_smo,&
         LIS_sda,LIS_shr,LIS_smn,LIS_sss,real(ts2) )
    if (DataAvailabilityStartTime .gt. LIS_StartTime) then
       write(LIS_logunit,*) 'LIS start time is earlier than TRMM data availability time!'
       write(LIS_logunit,*) 'Program stopping ... '
       call LIS_endrun()
    endif
    ! SY: End LIS time check against TRMM data availability time

    call readcrd_TRMM3B42V6()

    ! SY: Start for obtaining intended TRMM 3B42V6 read time step
    do n=1, LIS_rc%nnest
       TRMM3B42V6_struc(n)%ts = 3*60*60 
    ! SY: Because a LIS 3-hr timestep is typically valid from Z to Z, 
    !     while TRMM 3-hr data validity is from 0.5 Z to 0.5 Z 
       call LIS_update_timestep(LIS_rc, n, TRMM3B42V6_struc(n)%ts)

       call LIS_registerAlarm("TRMM 3B42V6 alarm",&
            TRMM3B42V6_struc(n)%ts,&
            TRMM3B42V6_struc(n)%ts, alarm_offset=90*60)
    enddo

    LIS_rc%met_nf(findex) = 2 !number of met variables in TRMM forcing

    do n=1,LIS_rc%nnest

       allocate(TRMM3B42V6_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(TRMM3B42V6_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       TRMM3B42V6_struc(n)%metdata1 = 0
       TRMM3B42V6_struc(n)%metdata2 = 0

      gridDesci = 0
      gridDesci(1) = 0
      gridDesci(2) = 1440
      gridDesci(3) = 400
      gridDesci(4) = -49.875
      gridDesci(5) = -179.875
      gridDesci(6) = 128
      gridDesci(7) = 49.875
      gridDesci(8) = 179.875
      gridDesci(9) = 0.25
      gridDesci(10) = 0.25
      gridDesci(20) = 64

      TRMM3B42V6_struc(n)%ncold = 1440
      TRMM3B42V6_struc(n)%nrold = 400

      ! SY: Begin initialization of TRMM3B42V6_struc time components to 0
      TRMM3B42V6_struc(n)%TRMM3B42V6time_TStepStart = 0.0
      TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart_Previous = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart_Previous = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart_Previous = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart_Previous = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6time_TStepEnd = 0.0
      TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd_Previous = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd_Previous = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd_Previous = 0
      TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd_Previous = 0
      TRMM3B42V6_struc(n)%LIS_timeAtTStepStart = 0.0
      TRMM3B42V6_struc(n)%LIS_timeAtTStepEnd = 0.0
      ! SY: End initialization of TRMM3B42V6_struc time components to 0

       allocate(TRMM3B42V6_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(TRMM3B42V6_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(TRMM3B42V6_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(TRMM3B42V6_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(TRMM3B42V6_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(TRMM3B42V6_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(TRMM3B42V6_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(TRMM3B42V6_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

      if(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
         call conserv_interp_input(n,gridDesci,&
              TRMM3B42V6_struc(n)%n112,TRMM3B42V6_struc(n)%n122,&
              TRMM3B42V6_struc(n)%n212,TRMM3B42V6_struc(n)%n222,&
              TRMM3B42V6_struc(n)%w112,TRMM3B42V6_struc(n)%w122,&
              TRMM3B42V6_struc(n)%w212,TRMM3B42V6_struc(n)%w222)
      elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then ! SY
         ! SY: begin for nearest neighbor
         call neighbor_interp_input(n,gridDesci,&
              TRMM3B42V6_struc(n)%n112)
         ! SY: end for nearest neighbor
      else
         write(LIS_logunit,*) 'This interpolation not defined for TRMM data'
         write(LIS_logunit,*) 'Program stopping ... '
         call LIS_endrun()
      endif
      
    enddo

  end subroutine init_TRMM3B42V6

end module TRMM3B42V6_forcingMod
                                                                        
