!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module TRMM3B42RTV7_forcingMod
!BOP
! !MODULE: TRMM3B42RTV7_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation product derived from 
!  the TRMM real-time multi-satellite precipitation analysis (3B42RT)
!  from NASA GSFC (Huffman et al.2003)
!   
!  Huffman, G.J.,R.F. Adler, E.F.Stocker, D.T.Bolvin, and E.J.Nelkin
!  (2003): Analysis of TRMM 3-hourly multi-satellite precipitation
!  estimates computed in real and post-real time. Combined preprints
!  CD-ROM, 83rd AMS Annual Meeting, Paper P4.11 in 12th Conference 
!  on Sat. Meteor. and Oceanog. 9-13 Feb 2003, Long Beach CA. 6pp. 
! 
!  The implementation in LIS has the derived data type {\tt TRMM3B42RT\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
! 
!  They are desribed below: 
! \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[directory]
!    Directory containing the input data
!  \item[griduptime1]
!    The time to switch the input resolution or file naming convention
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Array containing the weights of the input grid
!    for each grid point in LIS, for neighbor search.
!  \end{description}
!
! !REVISION HISTORY:
!  Yudong Tian,  10/26/2010: Updated to support reading raw .gz files 
!  KR Arsenault, 01/06/2015: Updated to support latest V7 files

! !USES:
  use ESMF
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_TRMM3B42RTV7      !defines the native resolution of
                                 !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: TRMM3B42RTV7_struc

!EOP

  type, public ::  TRMM3B42RTV7_type_dec
     real                     :: ts
     integer                  :: nc
     integer                  :: nr
     character(len=LIS_CONST_PATH_LEN) :: directory  
     real*8                   :: time_TStepStart ! SY
     integer                  :: yr_TStepStart ! SY
     integer                  :: mo_TStepStart ! SY
     integer                  :: da_TStepStart ! SY
     integer                  :: hr_TStepStart ! SY
     integer                  :: yr_TStepStart_Previous ! SY
     integer                  :: mo_TStepStart_Previous ! SY
     integer                  :: da_TStepStart_Previous ! SY
     integer                  :: hr_TStepStart_Previous ! SY
     real*8                   :: time_TStepEnd ! SY
     integer                  :: yr_TStepEnd ! SY
     integer                  :: mo_TStepEnd ! SY
     integer                  :: da_TStepEnd ! SY
     integer                  :: hr_TStepEnd ! SY
     integer                  :: yr_TStepEnd_Previous ! SY
     integer                  :: mo_TStepEnd_Previous ! SY
     integer                  :: da_TStepEnd_Previous ! SY
     integer                  :: hr_TStepEnd_Previous ! SY
     real*8                   :: LIS_timeAtTStepStart ! SY
     real*8                   :: LIS_timeAtTStepEnd   ! SY
     real*8                   :: griduptime1
     logical                  :: gridchange1
     integer                  :: mi
   ! Budget-bilinear:
     integer, allocatable     :: n112(:,:)
     integer, allocatable     :: n122(:,:)
     integer, allocatable     :: n212(:,:)
     integer, allocatable     :: n222(:,:)
     real,    allocatable     :: w112(:,:),w122(:,:)
     real,    allocatable     :: w212(:,:),w222(:,:)
   ! Neighbor:
     integer, allocatable     :: n113(:)

     integer           :: nIter, st_iterid,en_iterid  ! Forecast mode

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type TRMM3B42RTV7_type_dec
  
  type(TRMM3B42RTV7_type_dec), allocatable :: TRMM3B42RTV7_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_TRMM3B42RTV7
! \label{init_TRMM3B42RTV7}
! 
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_TRMM3B42RTV7(findex)
! !USES: 
    use LIS_coreMod, only: LIS_rc, LIS_config, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep, &
                               LIS_tick, LIS_parseTimeString, &
                               LIS_registerAlarm ! SY
    use LIS_logMod,  only : LIS_logunit, LIS_endrun, &
                            LIS_getNextUnitNumber, &
                            LIS_releaseUnitNumber, LIS_verify 
    use LIS_FORC_AttributesMod
    use LIS_forecastMod

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for TRMM 3B42RTV7
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_TRMM3B42RTV7](\ref{readcrd_TRMM3B42RTV7}) \newline
!     reads the runtime options specified for TRMM 3B42RTV7 data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
    real    :: gridDesci(50)
    integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: n 
    integer :: ts1,ts2 
    integer :: LIS_syr,LIS_smo,LIS_sda,LIS_shr,LIS_smn,LIS_sss ! SY

    real*8  :: DataAvailabilityStartTime ! SY
    real*8  :: LIS_StartTime  ! SY
    integer :: rc             
    character (len=10):: time ! SY

    allocate(TRMM3B42RTV7_struc(LIS_rc%nnest))
! ___________

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    ! SY: Begin LIS time check against TRMM data availability time
    yr1 = 2000
    mo1 = 03
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

    ! SY: End LIS time check against TRMM data availability time
    if (DataAvailabilityStartTime .gt. LIS_StartTime) then
       write(LIS_logunit,*) "LIS start time is earlier than TRMM data availability time!"
       write(LIS_logunit,*) "Program stopping ... "
       call LIS_endrun
    endif

    call readcrd_TRMM3B42RTV7()

    ! SY: Start for obtaining intended TRMM 3B42RT read time step
    do n=1, LIS_rc%nnest

       TRMM3B42RTV7_struc(n)%ts = 3*60*60
       call LIS_update_timestep(LIS_rc, n, TRMM3B42RTV7_struc(n)%ts)

       call LIS_registerAlarm("TRMM 3B42RT alarm",&
            TRMM3B42RTV7_struc(n)%ts,&
            TRMM3B42RTV7_struc(n)%ts, alarm_offset=90*60)
    enddo

    LIS_rc%met_nf(findex) = 2 !number of met variables in TRMM forcing

    TRMM3B42RTV7_struc(:)%nc = 1440
    TRMM3B42RTV7_struc(:)%nr = 480

    do n=1,LIS_rc%nnest

       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then

         if(mod(LIS_rc%nensem(n),&
             LIS_forecast_struc(1)%niterations).ne.0) then
            write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
            write(LIS_logunit,*) '[ERR] of the number of iterations '
            write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
            write(LIS_logunit,*) '[ERR] niter  = ',LIS_forecast_struc(1)%niterations
            call LIS_endrun()
         endif

         TRMM3B42RTV7_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
         TRMM3B42RTV7_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
         TRMM3B42RTV7_struc(n)%nIter = LIS_forecast_struc(1)%niterations

         allocate(TRMM3B42RTV7_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
         allocate(TRMM3B42RTV7_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

       else  ! Regular retrospective or non-forecast mode:

          TRMM3B42RTV7_struc(n)%st_iterid = 1
          TRMM3B42RTV7_struc(n)%en_iterId = 1
          TRMM3B42RTV7_struc(n)%nIter = 1

          allocate(TRMM3B42RTV7_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(TRMM3B42RTV7_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
       endif

       TRMM3B42RTV7_struc(n)%metdata1 = 0
       TRMM3B42RTV7_struc(n)%metdata2 = 0

       gridDesci = 0
       gridDesci(1) = 0
       gridDesci(2) = TRMM3B42RTV7_struc(n)%nc
       gridDesci(3) = TRMM3B42RTV7_struc(n)%nr
!       gridDesci(4) = -49.875
       gridDesci(4) = -59.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128
!       gridDesci(7) = 49.875
       gridDesci(7) = 59.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
       
       ! SY: Begin initialization of TRMM3B42RTV7_struc time components to 0
       TRMM3B42RTV7_struc(n)%time_TStepStart = 0.0
       TRMM3B42RTV7_struc(n)%yr_TStepStart = 0
       TRMM3B42RTV7_struc(n)%mo_TStepStart = 0
       TRMM3B42RTV7_struc(n)%da_TStepStart = 0
       TRMM3B42RTV7_struc(n)%hr_TStepStart = 0
       TRMM3B42RTV7_struc(n)%yr_TStepStart_Previous = 0
       TRMM3B42RTV7_struc(n)%mo_TStepStart_Previous = 0
       TRMM3B42RTV7_struc(n)%da_TStepStart_Previous = 0
       TRMM3B42RTV7_struc(n)%hr_TStepStart_Previous = 0
       TRMM3B42RTV7_struc(n)%time_TStepEnd = 0.0
       TRMM3B42RTV7_struc(n)%yr_TStepEnd = 0
       TRMM3B42RTV7_struc(n)%mo_TStepEnd = 0
       TRMM3B42RTV7_struc(n)%da_TStepEnd = 0
       TRMM3B42RTV7_struc(n)%hr_TStepEnd = 0
       TRMM3B42RTV7_struc(n)%yr_TStepEnd_Previous = 0
       TRMM3B42RTV7_struc(n)%mo_TStepEnd_Previous = 0
       TRMM3B42RTV7_struc(n)%da_TStepEnd_Previous = 0
       TRMM3B42RTV7_struc(n)%hr_TStepEnd_Previous = 0
       TRMM3B42RTV7_struc(n)%LIS_timeAtTStepStart = 0.0
       TRMM3B42RTV7_struc(n)%LIS_timeAtTStepEnd = 0.0
       ! SY: End initialization of TRMM3B42RTV7_struc time components to 0
       
    !- Select and define interp input arrays:
       select case( LIS_rc%met_interp(findex) )

        case( "budget-bilinear" )

          allocate(TRMM3B42RTV7_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               TRMM3B42RTV7_struc(n)%n112,TRMM3B42RTV7_struc(n)%n122,&
               TRMM3B42RTV7_struc(n)%n212,TRMM3B42RTV7_struc(n)%n222,&
               TRMM3B42RTV7_struc(n)%w112,TRMM3B42RTV7_struc(n)%w122,&
               TRMM3B42RTV7_struc(n)%w212,TRMM3B42RTV7_struc(n)%w222)

        case( "neighbor" )   
          allocate(TRMM3B42RTV7_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n,gridDesci,&
               TRMM3B42RTV7_struc(n)%n113)

        case default
          write(LIS_logunit,*) "This interpolation option not defined yet for 3B42RT data"
          write(LIS_logunit,*) "Program stopping ... "
          call LIS_endrun

       end select

       yr1 = 2012     ! Grid update time 3B42RT V7 file name changes from 
       mo1 = 11       ! 3B42RT.2007101603.7R2.bin to 3B42RT.2008101603.7.bin
       da1 = 7 
       hr1 = 3
       mn1 = 0; ss1 = 0
       call LIS_date2time(TRMM3B42RTV7_struc(n)%griduptime1,&
            updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       TRMM3B42RTV7_struc(n)%gridchange1 = .true.

    enddo

  end subroutine init_TRMM3B42RTV7

end module TRMM3B42RTV7_forcingMod
