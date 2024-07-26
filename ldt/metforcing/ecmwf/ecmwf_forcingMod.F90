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
!  00, 06, 12, and 18 UTC, which are used in LDT.  
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
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n12,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !USES:
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
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

  type, public ::  ecmwf_type_dec 
     real         :: ts
     integer      :: ncold, nrold   !AWIPS 212 dimensions
     character(len=LDT_CONST_PATH_LEN) :: ecmwfdir   !ecmwf Forcing Directory
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs23r4  !elevation for 2001-2002
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs25r1  !elevation for 2003-2006
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs30r1  !elevation for 2006 - jun 2008
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs33r1  !elevation for jun 2008 - mar 2009
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs35r2  !elevation for mar 2009 - sep 2009
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs35r3  !elevation for sep 2009 - jan 2010
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs36r1  !elevation for jan 2010 - may 2011
     character(len=LDT_CONST_PATH_LEN) :: elevfileifs37r2  !elevation for may 2011 onward
     real*8       :: ecmwftime1,ecmwftime2
     real*8       :: griduptime1,griduptime2,griduptime3,griduptime4, &
                     griduptime5,griduptime6,griduptime7
     logical      :: gridchange1,gridchange2,gridchange3,gridchange4, &
                     gridchange5,gridchange6,gridchange7
     integer      :: mi
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
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_logMod,     only : LDT_logunit,LDT_endrun
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep

    implicit none
    
    integer, intent(in)      :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for ECMWF
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
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
    integer :: n 
    integer :: i
    integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    real    :: gridDesci(20)

    allocate(ecmwf_struc(LDT_rc%nnest))

    write(LDT_logunit,*)"MSG: Initializing ECMWF forcing grid ... "

 !- Read in config file entries:
    call readcrd_ecmwf(findex)

    do n=1, LDT_rc%nnest
       ecmwf_struc(n)%ts = 3*3600 
       call LDT_update_timestep(LDT_rc, n, ecmwf_struc(n)%ts)

     ! Set local - 1 hour timestep (to replicate model timestep):
       call LDT_update_timestep(LDT_rc, n, 3600.)
    enddo

  ! Metforcing and parameter grid info:
    ecmwf_struc(:)%ncold  = 1440
    ecmwf_struc(:)%nrold  = 601

    gridDesci = 0

    gridDesci(1) = 0
    gridDesci(2) = ecmwf_struc(1)%ncold
    gridDesci(3) = ecmwf_struc(1)%nrold
    gridDesci(4) = 90.000
    gridDesci(5) = -180.000
    gridDesci(6) = 128
    gridDesci(7) = -60.000
    gridDesci(8) = 179.750
    gridDesci(9) = 0.250
    gridDesci(10) = 0.250
    gridDesci(20) = 64

    LDT_rc%met_gridDesc(findex,:) = gridDesci(:) 

  ! Set number of cols and rows for all ECMWF elevation files:
    do i = 0, 7   ! Number of file revision changes
       LDT_rc%met_nc(findex+i) = ecmwf_struc(1)%ncold
       LDT_rc%met_nr(findex+i) = ecmwf_struc(1)%nrold
       if( i /= 0 ) then
         LDT_rc%met_gridDesc(findex+i,:) = LDT_rc%met_gridDesc(findex,:)
       endif
    end do

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 9   ! number of met variables in ECMWF forcing
    LDT_rc%met_ts(findex) = 3*3600
    LDT_rc%met_zterp(findex) = .true.

    do n=1,LDT_rc%nnest

       ecmwf_struc(n)%mi = ecmwf_struc(n)%ncold*ecmwf_struc(n)%nrold

       ecmwf_struc(n)%findtime1 = 0 
       ecmwf_struc(n)%findtime2 = 0 

       yr1 = 2003     !grid update time to IFS25R1 after data gap
       mo1 = 01
       da1 = 08
       hr1 = 00
       mn1 = 0; ss1 = 0
       call LDT_date2time( ecmwf_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
       
       yr1 = 2006     !grid update time to IFS30R1
       mo1 = 02
       da1 = 01
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LDT_date2time(ecmwf_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2008     !grid update time to IFS33R1
       mo1 = 06
       da1 = 03
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LDT_date2time(ecmwf_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2009     !grid update time to IFS35R2
       mo1 = 03
       da1 = 10
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LDT_date2time(ecmwf_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2009     !grid update time to IFS35R3
       mo1 = 09
       da1 = 08
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LDT_date2time(ecmwf_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2010     !grid update time to IFS36R1
       mo1 = 01
       da1 = 26
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LDT_date2time(ecmwf_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2011     !grid update time to IFS37R2
       mo1 = 05
       da1 = 18
       hr1 = 06
       mn1 = 0; ss1 = 0
       call LDT_date2time(ecmwf_struc(n)%griduptime7,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       ecmwf_struc(n)%gridchange1 = .true.
       ecmwf_struc(n)%gridchange2 = .true.
       ecmwf_struc(n)%gridchange3 = .true.
       ecmwf_struc(n)%gridchange4 = .true.
       ecmwf_struc(n)%gridchange5 = .true.
       ecmwf_struc(n)%gridchange6 = .true.
       ecmwf_struc(n)%gridchange7 = .true.

     ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )

        case( "bilinear" )

          allocate(ecmwf_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(:),&
               ecmwf_struc(n)%n111,ecmwf_struc(n)%n121,&
               ecmwf_struc(n)%n211,ecmwf_struc(n)%n221,&
               ecmwf_struc(n)%w111,ecmwf_struc(n)%w121,&
               ecmwf_struc(n)%w211,ecmwf_struc(n)%w221)

        case( "budget-bilinear" )

          allocate(ecmwf_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwf_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(:),&
               ecmwf_struc(n)%n111,ecmwf_struc(n)%n121,&
               ecmwf_struc(n)%n211,ecmwf_struc(n)%n221,&
               ecmwf_struc(n)%w111,ecmwf_struc(n)%w121,&
               ecmwf_struc(n)%w211,ecmwf_struc(n)%w221)

          allocate(ecmwf_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwf_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(:),&
               ecmwf_struc(n)%n112,ecmwf_struc(n)%n122,&
               ecmwf_struc(n)%n212,ecmwf_struc(n)%n222,&
               ecmwf_struc(n)%w112,ecmwf_struc(n)%w122,&
               ecmwf_struc(n)%w212,ecmwf_struc(n)%w222)

       case default
          write(LDT_logunit,*) "ERR MSG:  Meteorological grid transform option, ",&
                trim(LDT_rc%met_gridtransform(findex)), ", not known."
          write(LDT_logunit,*) " This interpolation scheme currently not supported for ECMWF reader."
          write(LDT_logunit,*) "Program stopping ..."
          call LDT_endrun()
       end select

       if( LDT_rc%met_ecor(findex) .ne. "none" ) then
          call read_ecmwfelev_ldtproc(n, findex, 0)
       endif

    enddo   ! End nest loop

  end subroutine init_ecmwf

end module ecmwf_forcingMod
