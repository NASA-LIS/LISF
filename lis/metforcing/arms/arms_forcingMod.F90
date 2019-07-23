!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module arms_forcingMod
!BOP
! !MODULE: arms_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various Walnut Gulch stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt arms\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 08 Jan 2004: Mike Tischler and Matt Garcia : Initial Specification
! 28 Jun 2009: Sujay Kumar: Reengineered for LIS version 6
! 
  use ESMF

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_ARMS      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: arms_struc
!EOP

  type, public :: arms_type_dec
     real          :: ts 
     character*40  :: armsfile
     character*40  :: armspcpfile
     character*40  :: mdatafile
     real          :: undef
     real*8        :: starttime,armstime1,armstime2
     logical       :: startRead
     integer       :: nstns
     real, allocatable :: stnwt(:,:)
     real, allocatable :: stnlat(:),stnlon(:)
     type(ESMF_Time) :: startDate
     type(ESMF_TimeInterval) :: timestep

!forcing data to be stored
     real,  allocatable :: precip(:,:)
     real,  allocatable :: t_a(:)
     real,  allocatable :: q_a(:)
     real,  allocatable :: swd(:)
     real,  allocatable :: lwd(:)
     real,  allocatable :: uwind(:)
     real,  allocatable :: psurf(:)

     integer        :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type arms_type_dec

  type(arms_type_dec), allocatable :: arms_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_ARMS
! \label{init_ARMS}
! 
! !REVISION HISTORY: 
! 28 Jun 2009: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_ARMS(findex)
! !USES:
    use UTM_utils,  only : UTM2geo
    use LIS_coreMod,only : LIS_rc
    use LIS_timeMgrMod, only : LIS_date2time, LIS_calendar, &
         LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun, LIS_verify, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber

    implicit none
    integer, intent(in) :: findex

! !DESCRIPTION: 
!  This routines reads the runtime configurations for using the
!  ARMS station data. Using the metadata provided for the 
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_arms](\ref{readcrd_arms}) \newline
!     reads the runtime options specified for ARMS station data
!   \item[utm2geo](\ref{utm2geo}) \newline
!     converts utm coordinates (easting, northing) to latitude and longitude values.
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP
    integer, parameter :: ndata = 577
    integer      :: n 
    integer      :: j
    integer      :: status
    integer      :: ftn
    real         :: stn_north, stn_east, dummy
    integer      :: stnid
!    integer      :: col, row
!    real         :: stnlocs(660,333)

!    stnlocs = -9999.0

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the ARMS forcing reader '
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(arms_struc(LIS_rc%nnest))
    call readcrd_arms()

    do n=1, LIS_rc%nnest
       arms_struc(n)%ts = 60*60 
       call LIS_update_timestep(LIS_rc, n, arms_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9 

    do n=1,LIS_rc%nnest

       allocate(arms_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(arms_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       arms_struc(n)%metdata1 = 0
       arms_struc(n)%metdata2 = 0

       arms_struc(n)%findtime1 = 0 
       arms_struc(n)%findtime2 = 0
       
       call ESMF_TimeSet(arms_struc(n)%startDate, yy=1990,mm=7,dd=23,&
         h=0,calendar=LIS_calendar,rc=status)
       call LIS_verify(status,'error in ESMF_TimeSet in init_ARMS')
    
       call ESMF_TimeIntervalSet(arms_struc(n)%timeStep,h=1,rc=status)
       call LIS_verify(status,'error in ESMF_TimeIntervalSet in init_ARMS')
    enddo
    
!----------------------------------------------------------------------------------
! Read the station metadata (northing, easting) information 
!----------------------------------------------------------------------------------
    do n=1, LIS_rc%nnest
       allocate(arms_struc(n)%stnlat(arms_struc(n)%nstns))
       allocate(arms_struc(n)%stnlon(arms_struc(n)%nstns))

       ftn = LIS_getNextUnitNumber()
       
       open(ftn,file=trim(arms_struc(n)%mdatafile),form='formatted')
       do j=1,arms_struc(n)%nstns
          read(ftn,*) stnid, stn_east, stn_north, dummy
          
!          col = nint((stn_east-LIS_rc%gridDesc(n,5))/40.0)+1
!          row = nint((stn_north-LIS_rc%gridDesc(n,4))/40.0)+1

!          stnlocs(col,row) = 1.0
          
          call utm2geo(nint(LIS_rc%gridDesc(n,10)), stn_north, stn_east, &
               arms_struc(n)%stnlat(j), arms_struc(n)%stnlon(j))
          write(LIS_logunit,*) 'stn : ',stnid, arms_struc(n)%stnlat(j), &
               arms_struc(n)%stnlon(j)
       enddo
       call LIS_releaseUnitNumber(ftn)

!       open(100,file='precip.bin',form='unformatted')
!       write(100) stnlocs
!       close(100)       
!       stop
!----------------------------------------------------------------------------------
! Compute the station weights
!----------------------------------------------------------------------------------
       allocate(arms_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),arms_struc(n)%nstns))
       call compute_stnwts(arms_struc(n)%nstns, LIS_rc%gridDesc, & 
            arms_struc(n)%stnlat, arms_struc(n)%stnlon, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), arms_struc(n)%stnwt)
!----------------------------------------------------------------------------------
! Allocate variables for future storage
!----------------------------------------------------------------------------------
       allocate(arms_struc(n)%precip(ndata, arms_struc(n)%nstns))
       allocate(arms_struc(n)%t_a(ndata))
       allocate(arms_struc(n)%q_a(ndata))
       allocate(arms_struc(n)%swd(ndata))
       allocate(arms_struc(n)%lwd(ndata))
       allocate(arms_struc(n)%uwind(ndata))
       allocate(arms_struc(n)%psurf(ndata))
    enddo
  end subroutine init_ARMS

end module arms_forcingMod
