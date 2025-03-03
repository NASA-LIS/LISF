!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nam242_forcingMod
!BOP
! !MODULE: nam242_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Data
!  Assimilation System (NAM) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. NAM forcing variables are produced
!  on a quadratic gaussian grid. LIS uses the 00, 03, 06 and as 
!  needed, the 09 forecasts. The forecasts are produced at 6 hr intervals. 
!
!  The implementation in LIS has the derived data type {\tt nam242\_struc} that
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the NAM data
!  \item[ts]
!    Frequency in seconds of the forcing data
!  \item[namtime1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[namtime2]
!    The nearest, next 3 hour instance of the incoming 
!    data (as a real time).
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for 
!    temporal interpolation.
!  \item[namdir]
!    Directory containing the input data
!  \item[elevfile]
!    File with the elevation definition for the input data. 
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_nam242    !defines the native resolution of 
                           !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nam242_struc
!EOP

  type, public        :: nam242_type_dec
     integer          :: ncold, nrold
     integer          :: nmif
     character(len=LIS_CONST_PATH_LEN) :: namdir   !NAM Forcing Directory
     character(len=LIS_CONST_PATH_LEN) :: elevfile
     real             :: ts
     real*8           :: namtime1,namtime2
     integer          :: findtime1,findtime2
     integer          :: mi
     integer, allocatable :: n111(:)
     integer, allocatable :: n121(:)
     integer, allocatable :: n211(:)
     integer, allocatable :: n221(:)
     real, allocatable    :: w111(:),w121(:)
     real, allocatable    :: w211(:),w221(:)
     integer, allocatable :: n112(:,:)
     integer, allocatable :: n122(:,:)
     integer, allocatable :: n212(:,:)
     integer, allocatable :: n222(:,:)
     real, allocatable    :: w112(:,:),w122(:,:)
     real, allocatable    :: w212(:,:),w222(:,:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type nam242_type_dec
  
  type(nam242_type_dec), allocatable :: nam242_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_nam242
!  \label{init_nam242}
!
! !REVISION HISTORY: 
!     Sep 2012: NOHRSC/NOAA: Initial specification
! 
! !INTERFACE:
  subroutine init_nam242(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex

! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NAM
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}). Based on the NAM data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the NAM resolution switches are also 
!  defined in this routine. 
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_nam242](\ref{readcrd_nam242}) \newline
!     reads the runtime options specified for NAM data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP

    integer :: n
    real    :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1,tdoy
    real    :: upgmt, tgmt

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the NAM-242 grid forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(nam242_struc(LIS_rc%nnest))

    call readcrd_nam242()

    do n=1, LIS_rc%nnest
       nam242_struc(n)%ts = 3*60*60 
       call LIS_update_timestep(LIS_rc, n, nam242_struc(n)%ts)
    enddo

    nam242_struc(:)%nmif  = 9
    LIS_rc%met_nf(findex) = 9
    
    do n=1,LIS_rc%nnest

       allocate(nam242_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(nam242_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       nam242_struc(n)%metdata1 = 0
       nam242_struc(n)%metdata2 = 0

       gridDesci = 0
       gridDesci(1) = 5
       gridDesci(2) = 553
       gridDesci(3) = 425
       gridDesci(4) = 30
       gridDesci(5) = -173
       gridDesci(6) = 0       ! not used
       gridDesci(7) = -135    ! Dagang question
       gridDesci(8) = 11.25
       gridDesci(9) = 11.25
       gridDesci(10) = 60     ! Dagang question
       gridDesci(11) = -135   ! Dagang question
       gridDesci(20) = 0      ! Dagang question

       nam242_struc(n)%ncold = 553
       nam242_struc(n)%nrold = 425
       nam242_struc(n)%mi    = nam242_struc(n)%ncold*nam242_struc(n)%nrold
       
!Setting up weights for Interpolation
       if(LIS_rc%met_interp(findex).eq."bilinear") then 
          allocate(nam242_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               nam242_struc(n)%n111,nam242_struc(n)%n121,&
               nam242_struc(n)%n211,nam242_struc(n)%n221,&
               nam242_struc(n)%w111,nam242_struc(n)%w121,&
               nam242_struc(n)%w211,nam242_struc(n)%w221)

       elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 
          allocate(nam242_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nam242_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               nam242_struc(n)%n111,nam242_struc(n)%n121,&
               nam242_struc(n)%n211,nam242_struc(n)%n221,&
               nam242_struc(n)%w111,nam242_struc(n)%w121,&
               nam242_struc(n)%w211,nam242_struc(n)%w221)

          allocate(nam242_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nam242_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nam242_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nam242_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nam242_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nam242_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nam242_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nam242_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               nam242_struc(n)%n112,nam242_struc(n)%n122,&
               nam242_struc(n)%n212,nam242_struc(n)%n222,&
               nam242_struc(n)%w112,nam242_struc(n)%w122,&
               nam242_struc(n)%w212,nam242_struc(n)%w222)
       endif

       if ( LIS_rc%met_ecor(findex).eq."lapse-rate" ) then 
          call read_nam242_elev(n, findex, 1)
       endif

    enddo
  end subroutine init_nam242
end module nam242_forcingMod

