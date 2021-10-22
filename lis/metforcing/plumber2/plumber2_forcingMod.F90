!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module plumber2_forcingMod
!BOP
! !MODULE: plumber2_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various PLUMBER2 stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt plumber2\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 15 Sep 2021: Mark Beauharnois, Derived from Bondville and GSWP2 readers
! 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_plumber2 !defines the native resolution of 
                                       !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: plumber2_struc
!EOP

  type, public :: plumber2_type_dec

     real                     :: ts
     character*256            :: plumber2file
     real                     :: undef
     real*8                   :: starttime,plumber2time1,plumber2time2
     integer                  :: sYear, sMon, sDay, sHour, sMin, sSec
     integer                  :: eYear, eMon, eDay, eHour, eMin, eSec
     integer                  :: findtime1,findtime2
     integer                  :: utcoffset_sec
     integer                  :: nstns
     character*6,allocatable  :: stnid(:)
     character*3              :: veg_type
     real, allocatable        :: stnlat(:),stnlon(:)
     real, allocatable        :: stnwt(:,:)
     real, allocatable        :: metdata1(:,:) 
     real, allocatable        :: metdata2(:,:) 

     logical                  :: startRead

  end type plumber2_type_dec
  
  type(plumber2_type_dec), allocatable :: plumber2_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_plumber2
! \label{init_plumber2}
! 
! !INTERFACE:
  subroutine init_plumber2(findex)
! !USES:
    use LIS_coreMod,only     : LIS_rc
    use LIS_logMod, only     : LIS_logunit,LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    
    implicit none
    integer, intent(in) :: findex
    
! !DESCRIPTION:
!  This routines reads the runtime configurations for using the
!  PLUMBER2 station data. Using the metadata provided for the
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used.
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_plumber2](\ref{readcrd_plumber2}) \newline
!     reads the runtime options specified for PLUMBER2 station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP

    real    :: gmt
    integer :: i,n,doy

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the PLUMBER2 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
    
    allocate(plumber2_struc(LIS_rc%nnest))
    call readcrd_plumber2()

    do n=1, LIS_rc%nnest
       !! MCB next line commented for PLUMBER2 because 'ts' is now
       !! set in the 'readcrd_plumber2' routine for automation
       !! plumber2_struc(n)%ts = 1800
       call LIS_update_timestep(LIS_rc, n, plumber2_struc(n)%ts)
    enddo

    ! MCB NOTE: For PLUMBER2 the following met forcing variables are
    !           required:
    !  1. prec rain (total precip)
    !  2. psurf
    !  3. qair
    !  4. tair
    !  5. swdown
    !  6. lwdown
    !  7. wind u
    !  8. wind v
    !  9. lai    (currently not used)
    !  
    LIS_rc%met_nf(findex) = 8 !number of met variables from PLUMBER2

    do n = 1,LIS_rc%nnest
       allocate(plumber2_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(plumber2_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       plumber2_struc(n)%metdata1 = 0
       plumber2_struc(n)%metdata2 = 0

       plumber2_struc(n)%undef = -9999.0
       plumber2_struc(n)%nstns = 1

       allocate(plumber2_struc(n)%stnid(plumber2_struc(n)%nstns))
       allocate(plumber2_struc(n)%stnlat(plumber2_struc(n)%nstns))
       allocate(plumber2_struc(n)%stnlon(plumber2_struc(n)%nstns))

!!!!!! MCB NOTE: in the Bondville_forcingMod.F90 file
!!!!!! the call to LIS_date2time is outside of any NEST (n) loop, BUG?
       call LIS_date2time(plumber2_struc(n)%starttime,doy,gmt,     &
            plumber2_struc(n)%sYear, &
            plumber2_struc(n)%sMon,  &
            plumber2_struc(n)%sDay,  &
            plumber2_struc(n)%sHour, &
            plumber2_struc(n)%sMin,  &
            plumber2_struc(n)%sSec)

       allocate(plumber2_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            plumber2_struc(n)%nstns))

!      MCB Note: No spatial interpolation for PLUMBER2 test case(s)
!      call compute_stnwts(plumber2_struc(n)%nstns,LIS_rc%gridDesc,&
!           plumber2_struc(n)%stnlat,plumber2_struc(n)%stnlon,&
!           LIS_rc%lnc(n)*LIS_rc%lnr(n),plumber2_struc(n)%stnwt)
    enddo

  end subroutine init_plumber2

end module plumber2_forcingMod

