!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ecmwf_forcingMod
!BOP
! !MODULE: ecmwf_forcingMod
! 
! !DESCRIPTION: 
!  The operational, global analysis products from the European Center for 
!  Medium-Range Weather Forecasts (ECMWF) is available on a $T_L$511 triangular
!  truncation, linear reduced gaussian grid for four synoptic hours:
!  00, 06, 12, and 18 UTC, which are used in LIS.  
!  The radiation products from the data assimilation system at the ECMWF is based
!  on the scheme simulated after Morcrette's work~\cite{morcrette}.
!  The estimation of shortwave and longwave radiation 
!  includes schemes to include effects of 
!  adsorption by water vapor and other gases and aerosol scattering based
!  on different parameterizations~\cite{rothman,morcrette2}.
! 
!  This module contains variables and data structures that are used
!  for the implementation of the ECMWF forcing data. 
!  The data is global 0.25 degree dataset in latlon
!  projection, and at 3 hourly intervals. The derived
!  data type {\tt ecmwf\_struc}
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[ecmwftime1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[ecmwftime2]
!    The nearest, next 3 hour instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch Geopotential field
!  \item[griduptime2]
!    The time to switch Geopotential field
!  \item[ecmwfdir]
!    Directory containing the input data
!  \item[elevfile]
!    File with the elevation definition for the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n12,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_ecmwf      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ecmwf_struc

!EOP

  type, public     :: ecmwf_type_dec
     real          :: ts
     integer       :: ncold, nrold    ! AWIPS 212 dimensions
     character(len=LIS_CONST_PATH_LEN) :: ecmwfdir ! ecmwf Forcing Directory
     !character*45 :: elevfileifs23r4  ! elevation for 2001-2002
     !character*45 :: elevfileifs25r1  ! elevation for 2003-2006
     !character*45 :: elevfileifs30r1  ! elevation for 2006 - jun 2008
     !character*45 :: elevfileifs33r1  ! elevation for jun 2008 - mar 2009
     !character*45 :: elevfileifs35r2  ! elevation for mar 2009 - sep 2009
     !character*45 :: elevfileifs35r3  ! elevation for sep 2009 - jan 2010
     !character*45 :: elevfileifs36r1  ! elevation for jan 2010 - may 2011
     !character*45 :: elevfileifs37r2  ! elevation for may 2011 onward
     real*8        :: ecmwftime1,ecmwftime2
     real*8        :: griduptime1,griduptime2,griduptime3,griduptime4, &
                      griduptime5,griduptime6,griduptime7
     logical       :: gridchange1,gridchange2,gridchange3,gridchange4, &
                      gridchange5,gridchange6,gridchange7
     integer       :: mi
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)
     
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
     integer                :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type ecmwf_type_dec

  type(ecmwf_type_dec), allocatable :: ecmwf_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_ecmwf
! \label{init_ecmwf}
! 
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_ecmwf(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
    
    integer, intent(in)      :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for ECMWF
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_ecmwf](\ref{readcrd_ecmwf}) \newline
!     reads the runtime options specified for ECMWF data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real    :: gridDesci(LIS_rc%nnest,50)
    integer :: n 
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the ECMWF forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(ecmwf_struc(LIS_rc%nnest))
    call readcrd_ecmwf()

    do n=1, LIS_rc%nnest
       ecmwf_struc(n)%ts = 3*3600 
       call LIS_update_timestep(LIS_rc, n, ecmwf_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9
    ecmwf_struc(:)%ncold  = 1440
    ecmwf_struc(:)%nrold  = 601

    gridDesci = 0

    do n=1,LIS_rc%nnest

       allocate(ecmwf_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(ecmwf_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       ecmwf_struc(n)%metdata1 = 0
       ecmwf_struc(n)%metdata2 = 0

       ecmwf_struc(n)%findtime1 = 0 
       ecmwf_struc(n)%findtime2 = 0 

       gridDesci(n,1) = 0
       gridDesci(n,2) = ecmwf_struc(n)%ncold
       gridDesci(n,3) = ecmwf_struc(n)%nrold
       gridDesci(n,4) = 90.000
       gridDesci(n,5) = -180.000
       gridDesci(n,6) = 128
       gridDesci(n,7) = -60.000
       gridDesci(n,8) = 179.750
       gridDesci(n,9) = 0.250
       gridDesci(n,10) = 0.250
       gridDesci(n,20) = 64

       ecmwf_struc(n)%mi = ecmwf_struc(n)%ncold*ecmwf_struc(n)%nrold

       yr1 = 2003     !grid update time to IFS25R1 after data gap
       mo1 = 01
       da1 = 08
       hr1 = 00
       mn1 = 0; ss1 = 0
       call LIS_date2time( ecmwf_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
       
       yr1 = 2006     !grid update time to IFS30R1
       mo1 = 02
       da1 = 01
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LIS_date2time(ecmwf_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2008     !grid update time to IFS33R1
       mo1 = 06
       da1 = 03
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LIS_date2time(ecmwf_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2009     !grid update time to IFS35R2
       mo1 = 03
       da1 = 10
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LIS_date2time(ecmwf_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2009     !grid update time to IFS35R3
       mo1 = 09
       da1 = 08
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LIS_date2time(ecmwf_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2010     !grid update time to IFS36R1
       mo1 = 01
       da1 = 26
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LIS_date2time(ecmwf_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2011     !grid update time to IFS37R2
       mo1 = 05
       da1 = 18
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LIS_date2time(ecmwf_struc(n)%griduptime7,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       ecmwf_struc(n)%gridchange1 = .true.
       ecmwf_struc(n)%gridchange2 = .true.
       ecmwf_struc(n)%gridchange3 = .true.
       ecmwf_struc(n)%gridchange4 = .true.
       ecmwf_struc(n)%gridchange5 = .true.
       ecmwf_struc(n)%gridchange6 = .true.
       ecmwf_struc(n)%gridchange7 = .true.

!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(ecmwf_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               ecmwf_struc(n)%n111,ecmwf_struc(n)%n121,&
               ecmwf_struc(n)%n211,ecmwf_struc(n)%n221,&
               ecmwf_struc(n)%w111,ecmwf_struc(n)%w121,&
               ecmwf_struc(n)%w211,ecmwf_struc(n)%w221)


       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(ecmwf_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               ecmwf_struc(n)%n111,ecmwf_struc(n)%n121,&
               ecmwf_struc(n)%n211,ecmwf_struc(n)%n221,&
               ecmwf_struc(n)%w111,ecmwf_struc(n)%w121,&
               ecmwf_struc(n)%w211,ecmwf_struc(n)%w221)

          allocate(ecmwf_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               ecmwf_struc(n)%n112,ecmwf_struc(n)%n122,&
               ecmwf_struc(n)%n212,ecmwf_struc(n)%n222,&
               ecmwf_struc(n)%w112,ecmwf_struc(n)%w122,&
               ecmwf_struc(n)%w212,ecmwf_struc(n)%w222)
       endif
       if(trim(LIS_rc%met_ecor(findex)).ne."none") &
            call read_ecmwf_elev(n,findex,0)
    enddo

  end subroutine init_ecmwf
end module ecmwf_forcingMod
