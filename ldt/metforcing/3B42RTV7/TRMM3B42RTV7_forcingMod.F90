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
!  The implementation in LDT has the derived data type {\tt TRMM3B42RT\_struc}
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
!    for each grid point in LDT, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Array containing the weights of the input grid
!    for each grid point in LDT, for neighbor search.
!  \end{description}
!
! !REVISION HISTORY:
!  Yudong Tian,  10/26/2010: Updated to support reading raw .gz files 
!  KR Arsenault, 01/06/2015: Updated to support latest V7 files

! !USES:
  use ESMF
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

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
     character(len=LDT_CONST_PATH_LEN) :: directory  
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
     real*8                   :: LDT_timeAtTStepStart ! SY
     real*8                   :: LDT_timeAtTStepEnd   ! SY
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
    use LDT_coreMod, only: LDT_rc, LDT_config, LDT_domain
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep, &
                               LDT_tick, LDT_parseTimeString, &
                               LDT_registerAlarm ! SY
    use LDT_logMod,  only : LDT_logunit, LDT_endrun, &
                            LDT_getNextUnitNumber, &
                            LDT_releaseUnitNumber, LDT_verify 

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for TRMM 3B42RTV7
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_TRMM3B42RTV7](\ref{readcrd_TRMM3B42RTV7}) \newline
!     reads the runtime options specified for TRMM 3B42RTV7 data
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
    integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: n 
    integer :: ts1,ts2 ! SY
    integer :: LDT_syr,LDT_smo,LDT_sda,LDT_shr,LDT_smn,LDT_sss ! SY

    real*8  :: DataAvailabilityStartTime ! SY
    real*8  :: LDT_StartTime  ! SY
    integer :: rc             ! SY
    character (len=10):: time ! SY

    real :: gridDesci(20)

! _________________________________________________________________

    allocate(TRMM3B42RTV7_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing TRMM 3B42-RT V7 forcing grid ... "

    ! SY: Begin LDT time check against TRMM data availability time
    yr1 = 2000
    mo1 = 03
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
         LDT_sda,LDT_shr,LDT_smn,LDT_sss,real(ts2) )

    ! SY: End LDT time check against TRMM data availability time
    if (DataAvailabilityStartTime .gt. LDT_StartTime) then
       write(LDT_logunit,*) "LDT start time is earlier than TRMM data availability time!"
       write(LDT_logunit,*) "Program stopping ... "
       call LDT_endrun
    endif

    call readcrd_TRMM3B42RTV7()

    LDT_rc%met_nf(findex) = 2        ! number of met variables in TRMM forcing
    LDT_rc%met_ts(findex) = 3*3600  
    LDT_rc%met_zterp(findex) = .false.
    LDT_rc%met_proj(findex)  = "latlon"

    TRMM3B42RTV7_struc(:)%nc = 1440
    TRMM3B42RTV7_struc(:)%nr = 480
    LDT_rc%met_nc(findex) = TRMM3B42RTV7_struc(1)%nc
    LDT_rc%met_nr(findex) = TRMM3B42RTV7_struc(1)%nr

    gridDesci(:) = 0
    gridDesci(1) = 0
    gridDesci(2) = TRMM3B42RTV7_struc(1)%nc
    gridDesci(3) = TRMM3B42RTV7_struc(1)%nr
!    gridDesci(4) = -49.875
    gridDesci(4) = -59.875
    gridDesci(5) = -179.875
    gridDesci(6) = 128
!    gridDesci(7) = 49.875
    gridDesci(7) = 59.875
    gridDesci(8) = 179.875
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64

    LDT_rc%met_gridDesc(findex,:) = gridDesci(:)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return
       
    do  n=1,LDT_rc%nnest

       TRMM3B42RTV7_struc(n)%mi = TRMM3B42RTV7_struc(n)%nc*TRMM3B42RTV7_struc(n)%nr

       ! SY: Start for obtaining intended TRMM 3B42RT read time step
       TRMM3B42RTV7_struc(n)%ts = 3*60*60
       call LDT_update_timestep(LDT_rc, n, TRMM3B42RTV7_struc(n)%ts)

       call LDT_registerAlarm("TRMM 3B42RT alarm",&
            TRMM3B42RTV7_struc(n)%ts,&
            TRMM3B42RTV7_struc(n)%ts, alarm_offset=90*60)

     ! Set local - 1 hour timestep (to replicate model timestep):
       call LDT_update_timestep(LDT_rc, n, 3600.)

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
       TRMM3B42RTV7_struc(n)%LDT_timeAtTStepStart = 0.0
       TRMM3B42RTV7_struc(n)%LDT_timeAtTStepEnd = 0.0
       ! SY: End initialization of TRMM3B42RTV7_struc time components to 0
       
    !- Select and define interp input arrays:
       select case( LDT_rc%met_gridtransform(findex) )

        case( "budget-bilinear" )

          allocate(TRMM3B42RTV7_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(TRMM3B42RTV7_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(:),&
               TRMM3B42RTV7_struc(n)%n112,TRMM3B42RTV7_struc(n)%n122,&
               TRMM3B42RTV7_struc(n)%n212,TRMM3B42RTV7_struc(n)%n222,&
               TRMM3B42RTV7_struc(n)%w112,TRMM3B42RTV7_struc(n)%w122,&
               TRMM3B42RTV7_struc(n)%w212,TRMM3B42RTV7_struc(n)%w222)

        case( "neighbor" )   

          allocate(TRMM3B42RTV7_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n,gridDesci(:),&
               TRMM3B42RTV7_struc(n)%n113)

        case default
          write(LDT_logunit,*) "This interpolation option not defined yet for 3B42RT data"
          write(LDT_logunit,*) "Program stopping ... "
          call LDT_endrun

       end select

       yr1 = 2012     ! Grid update time 3B42RT V7 file name changes from 
       mo1 = 11       ! 3B42RT.2007101603.7R2.bin to 3B42RT.2008101603.7.bin
       da1 = 7 
       hr1 = 3
       mn1 = 0; ss1 = 0
       call LDT_date2time(TRMM3B42RTV7_struc(n)%griduptime1,&
            updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       TRMM3B42RTV7_struc(n)%gridchange1 = .true.

    enddo

  end subroutine init_TRMM3B42RTV7

end module TRMM3B42RTV7_forcingMod
