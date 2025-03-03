!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module TRMM3B42V7_forcingMod
!BOP
! !MODULE: TRMM3B42V7_forcingMod
!
! !REVISION HISTORY:
! 21 Jun 2013: Soni Yatheendradas; changes from earlier code to avoid
!              alternate file skip, jump to previous day TRMM and 
!              absence of rain rate weighting
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  NASA TRMM 3B42 (Version 7) merged analysis of precipitation (TRMM 3B42V7).
!  TRMM 3B42V7 algorithm consists of two separate steps to produce a
!  quasi-global 0.25 degree 3-hourly precipitation analysis. The first step
!  uses the TRMM VIRS and TMI orbit data (TRMM products 1B01 and 2A12)
!  and the monthly TMI/TRMM Combined Instrument (TCI) calibration parameters
!  (from TRMM product 3B31) to produce monthly IR calibration parameters. The
!  second step uses these derived monthly IR calibration parameters to adjust
!  the merged-IR precipitation data, which consists of GMS, GOES-E, GOES-W,
!  Meteosat-7, Meteosat-5, and NOAA-12 data.
!
!  The implementation in LDT has the derived data type {\tt TRMM3B42V7\_struc}
!  that includes the variables to specify the runtime options, and
!  the weights and neighbor information for spatial interpolation
!
!  They are desribed below:
! \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[TRMM3B42V7dir]
!    Directory containing the input data
!  \item[TRMM3B42V7time]
!    The nearest instance of the incoming
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \end{description}
!
! !USES:
    use LDT_constantsMod, only : LDT_CONST_PATH_LEN
    implicit none
    PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_TRMM3B42V7     !defines the native resolution of
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: TRMM3B42V7_struc

!EOP

  type, public :: TRMM3B42V7_type_dec
     real                     :: ts 
     integer                  :: nc
     integer                  :: nr
     character(len=LDT_CONST_PATH_LEN) :: TRMM3B42V7dir
     real*8                   :: TRMM3B42V7time_TStepStart ! SY
     integer                  :: TRMM3B42V7yr_TStepStart ! SY
     integer                  :: TRMM3B42V7mo_TStepStart ! SY
     integer                  :: TRMM3B42V7da_TStepStart ! SY
     integer                  :: TRMM3B42V7hr_TStepStart ! SY
     integer                  :: TRMM3B42V7yr_TStepStart_Previous ! SY
     integer                  :: TRMM3B42V7mo_TStepStart_Previous ! SY
     integer                  :: TRMM3B42V7da_TStepStart_Previous ! SY
     integer                  :: TRMM3B42V7hr_TStepStart_Previous ! SY
     real*8                   :: TRMM3B42V7time_TStepEnd ! SY
     integer                  :: TRMM3B42V7yr_TStepEnd ! SY
     integer                  :: TRMM3B42V7mo_TStepEnd ! SY
     integer                  :: TRMM3B42V7da_TStepEnd ! SY
     integer                  :: TRMM3B42V7hr_TStepEnd ! SY
     integer                  :: TRMM3B42V7yr_TStepEnd_Previous ! SY
     integer                  :: TRMM3B42V7mo_TStepEnd_Previous ! SY
     integer                  :: TRMM3B42V7da_TStepEnd_Previous ! SY
     integer                  :: TRMM3B42V7hr_TStepEnd_Previous ! SY
     real*8                   :: LDT_timeAtTStepStart ! SY
     real*8                   :: LDT_timeAtTStepEnd ! SY
     integer                  :: mi
     integer, allocatable     :: n112(:,:)
     integer, allocatable     :: n122(:,:)
     integer, allocatable     :: n212(:,:)
     integer, allocatable     :: n222(:,:)
     real,    allocatable     :: w112(:,:),w122(:,:)
     real,    allocatable     :: w212(:,:),w222(:,:)
  end type TRMM3B42V7_type_dec

  type(TRMM3B42V7_type_dec), allocatable :: TRMM3B42V7_struc(:)

contains

!BOP
!
! !ROUTINE: init_TRMM3B42V7
! \label{init_TRMM3B42V7}
!
! !REVISION HISTORY:
! 11Dec2003: Sujay Kumar; Initial Specification 
! 25Aug2006: Yudong Tian; Implementation for TRMM 3B42 V7
!
! !INTERFACE:
  subroutine init_TRMM3B42V7(findex)
! !USES:
    use LDT_coreMod, only: LDT_rc, LDT_config, LDT_domain
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep, &
                               LDT_tick, LDT_parseTimeString, &
                               LDT_registerAlarm ! SY
    use LDT_logMod,    only : LDT_logunit, LDT_endrun, &
                              LDT_getNextUnitNumber, LDT_releaseUnitNumber, LDT_verify 
    use ESMF

    implicit none
    
    integer,  intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for TRMM 3B42V7
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_TRMM3B42V7](\ref{readcrd_TRMM3B42V7}) \newline
!     reads the runtime options specified for TRMM 3B42V7 data
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: n
    integer :: ts1,ts2 ! SY
    integer :: LDT_syr,LDT_smo,LDT_sda,LDT_shr,LDT_smn,LDT_sss ! SY

    real*8  :: DataAvailabilityStartTime ! SY
    real*8  :: LDT_StartTime ! SY
    integer :: rc 
    character (len=10) :: time ! SY

    real    :: gridDesci(20)
! ________________________________________________

    allocate(TRMM3B42V7_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*) "MSG: Initializing TRMM 3B42V6 forcing grid ... "

    ! SY: Begin LDT time check against TRMM data availability time
    yr1 = 1998
    mo1 = 01
    da1 = 01
    hr1 = 0
    mn1 = 0; ss1 = 0
    ts1 = (-1) * 90 * 60 ! For 90 min behind 1st available TRMM data
    call LDT_tick(DataAvailabilityStartTime,updoy,upgmt,yr1,mo1,&
         da1,hr1,mn1,ss1,real(ts1) )
    LDT_syr = LDT_rc%syr
    LDT_smo = LDT_rc%smo
    LDT_sda = LDT_rc%sda
    LDT_shr = LDT_rc%shr
    LDT_smn = LDT_rc%smn
    LDT_sss = LDT_rc%sss
    ts2 = 0
    call LDT_tick(LDT_StartTime,updoy,upgmt,LDT_syr,LDT_smo,&
         LDT_sda,LDT_shr,LDT_smn,LDT_sss,real(ts2))
    if (DataAvailabilityStartTime .gt. LDT_StartTime) then
       write(LDT_logunit,*) 'LDT start time is earlier than TRMM data availability time!'
       write(LDT_logunit,*) 'Program stopping ... '
       call LDT_endrun()
    endif
    ! SY: End LDT time check against TRMM data availability time

    call readcrd_TRMM3B42V7()

    LDT_rc%met_nf(findex) = 2   ! number of met variables in TRMM forcing
    LDT_rc%met_ts(findex) = 3*3600
    LDT_rc%met_zterp(findex) = .false.
    LDT_rc%met_proj(findex)  = "latlon"

    TRMM3B42V7_struc%nc = 1440
    TRMM3B42V7_struc%nr = 400
    LDT_rc%met_nc(findex) = TRMM3B42V7_struc(1)%nc
    LDT_rc%met_nr(findex) = TRMM3B42V7_struc(1)%nr

 !- AGRMET Grid description:
    gridDesci(1)  = 0
    gridDesci(2)  = TRMM3B42V7_struc(1)%nc
    gridDesci(3)  = TRMM3B42V7_struc(1)%nr
    gridDesci(4)  = -49.875
    gridDesci(5)  = -179.875
    gridDesci(6)  = 128
    gridDesci(7)  = 49.875
    gridDesci(8)  = 179.875
    gridDesci(9)  = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64

    LDT_rc%met_gridDesc(findex,:)  = gridDesci(:)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    ! SY: Start for obtaining intended TRMM 3B42V7 read time step
    do n=1, LDT_rc%nnest

       TRMM3B42V7_struc(n)%ts = 3*60*60 
     ! SY: Because a LDT 3-hr timestep is typically valid from Z to Z, 
     !     while TRMM 3-hr data validity is from 0.5 Z to 0.5 Z 
       call LDT_update_timestep(LDT_rc, n, TRMM3B42V7_struc(n)%ts)

       call LDT_registerAlarm("TRMM 3B42V7 alarm",&
            TRMM3B42V7_struc(n)%ts,&
            TRMM3B42V7_struc(n)%ts, alarm_offset=90*60)

     ! Set local - 1 hour timestep (to replicate model timestep):
       call LDT_update_timestep(LDT_rc, n, 3600.)

    enddo

    do n=1,LDT_rc%nnest

      ! SY: Begin initialization of TRMM3B42V7_struc time components to 0
      TRMM3B42V7_struc(n)%TRMM3B42V7time_TStepStart = 0.0
      TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart_Previous = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart_Previous = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart_Previous = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart_Previous = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7time_TStepEnd = 0.0
      TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd_Previous = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd_Previous = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd_Previous = 0
      TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd_Previous = 0
      TRMM3B42V7_struc(n)%LDT_timeAtTStepStart = 0.0
      TRMM3B42V7_struc(n)%LDT_timeAtTStepEnd = 0.0
      ! SY: End initialization of TRMM3B42V7_struc time components to 0

       allocate(TRMM3B42V7_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(TRMM3B42V7_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(TRMM3B42V7_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(TRMM3B42V7_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(TRMM3B42V7_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(TRMM3B42V7_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(TRMM3B42V7_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(TRMM3B42V7_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

      select case( LDT_rc%met_gridtransform(findex) )

       case( "budget-bilinear" )
         call conserv_interp_input(n,LDT_rc%met_gridDesc(findex,:),&
              TRMM3B42V7_struc(n)%n112,TRMM3B42V7_struc(n)%n122,&
              TRMM3B42V7_struc(n)%n212,TRMM3B42V7_struc(n)%n222,&
              TRMM3B42V7_struc(n)%w112,TRMM3B42V7_struc(n)%w122,&
              TRMM3B42V7_struc(n)%w212,TRMM3B42V7_struc(n)%w222)

       case( "neighbor" )  ! SY
         call neighbor_interp_input(n,LDT_rc%met_gridDesc(findex,:),&
              TRMM3B42V7_struc(n)%n112)

       case default
         write(LDT_logunit,*) 'This interpolation not defined for TRMM data'
         write(LDT_logunit,*) 'Program stopping ... '
         call LDT_endrun()

       end select

    enddo
  end subroutine init_TRMM3B42V7

end module TRMM3B42V7_forcingMod
                                                                        
