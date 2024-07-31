!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdas_forcingMod
!BOP
! !MODULE: gdas_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Data
!  Assimilation System (GDAS) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. GDAS forcing variables are produced
!  on a quadratic gaussian grid. LIS uses the 00, 03, 06 and as
!  needed, the 09 forecasts. The forecasts are produced at 6 hr intervals.
!  The resolution of GDAS forcing varies as follows:
!
!   upto 2000/1/24          :   T126 (384x190)  grid \newline
!   2001/01/24 - 2002/10/29 :   T170 (512x256)  grid \newline
!   2002/10/29 - 2005/05/31 :   T254 (768x384)  grid \newline
!   2005/05/31 - 2010/07/27 :   T382 (1152x576) grid \newline
!   2010/07/28 - 2015/01/14 :   T574 (1760x880) grid \newline
!   2015/01/14 - onwards    :  T1534 (3072x1536) grid
!
!  On 2019/06/12 12Z, GDAS removed precipitation fields from the f00 data
!  files. The data fields in these files are now all instantaneous values.
!  When the reader is using data files after this time, a new subroutine will
!  be used that excludes precipitation as well as reads in instantaneous radiation
!  data. For data files prior to the switch, the reader will use the old subroutine.
!
!  The implementation in LIS has the derived data type {\tt gdas\_struc} that
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the GDAS data
!  \item[gdastime1]
!    The nearest, previous 3 hour instance of the incoming
!    data (as a real time).
!  \item[gdastime2]
!    The nearest, next 3 hour instance of the incoming
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch GDAS resolution to T126
!  \item[griduptime2]
!    The time to switch GDAS resolution to T170
!  \item[griduptime3]
!    The time to switch GDAS resolution to T254
!  \item[griduptime4]
!    The time to switch GDAS resolution to T382
!  \item[griduptime5]
!    The time to switch GDAS resolution to T574
!  \item[griduptime6]
!    The time to switch GDAS resolution to T1534
!  \item[datastructime1]
!    The time to switch to new data structure for f00 files
!    that removed precipitation fields.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \item[gdasdir]
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
  public :: init_GDAS      !defines the native resolution of
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdas_struc
!EOP

  type, public :: gdas_type_dec
     real          :: ts
     integer       :: ncold, nrold   !AWIPS 212 dimensions
     integer       :: nmif
     character(len=LIS_CONST_PATH_LEN) :: gdasdir   !GDAS Forcing Directory
     character*50  :: met_interp

     real*8        :: gdastime1, gdastime2
     real*8        :: griduptime1, griduptime2, griduptime3
     real*8        :: griduptime4, griduptime5, griduptime6
     real*8        :: datastructime1
     logical       :: gridchange1, gridchange2, gridchange3
     logical       :: gridchange4, gridchange5, gridchange6
     logical       :: dstrucchange1
     integer       :: mi
     integer       :: findtime1, findtime2

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

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 


  end type gdas_type_dec

  type(gdas_type_dec), allocatable :: gdas_struc(:)
contains

!BOP
!
! !ROUTINE: init_GDAS
!  \label{init_GDAS}
!
! !REVISION HISTORY:
! 11Dec2003: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine init_GDAS(findex)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep

    implicit none
! !ARGUMENTS:
    integer,  intent(in)    :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for GDAS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}). Based on the GDAS data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the GDAS resolution switches are also
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
!   \item[readcrd\_gdas](\ref{readcrd_gdas}) \newline
!     reads the runtime options specified for GDAS data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!  \end{description}
!EOP
    real :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n

    allocate(gdas_struc(LIS_rc%nnest))

    do n=1, LIS_rc%nnest
       gdas_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, gdas_struc(n)%ts)
    enddo

    call readcrd_gdas()

    gdas_struc(:)%nmif    = 9
    LIS_rc%met_nf(findex) = 9 !number of met variables in GDAS forcing

    gdas_struc(:)%ncold = 192
    gdas_struc(:)%nrold = 94

    do n=1,LIS_rc%nnest

       allocate(gdas_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gdas_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gdas_struc(n)%metdata1 = 0
       gdas_struc(n)%metdata2 = 0

       gdas_struc(n)%findtime1 = 0
       gdas_struc(n)%findtime2 = 0

       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
       gridDesci = 0
       gridDesci(1) = 4
       gridDesci(2) = 192
       gridDesci(3) = 94
       gridDesci(4) = 88.542
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) =  -88.542
       gridDesci(8) = -1.875
       gridDesci(9) = 1.875
       gridDesci(10) = 47
       gridDesci(20) = 0
       gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold

       ! This grid is good for some time in the 1990's.
       ! Look up the exact dates.
       yr1 = 1991
       mo1 = 01
       da1 = 01
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time( gdas_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2000
       mo1 = 01
       da1 = 24
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time( gdas_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2002     !grid update time ~ 0.469
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2005     !grid update time ~ 0.313
       mo1 = 05
       da1 = 31
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2010     !grid update time ~ 0.205
       mo1 = 07
       da1 = 28
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2015     !grid update time ~ 0.117
       mo1 = 01
       da1 = 14
       hr1 = 6
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       ! Set time for f00 data structure change
       yr1 = 2019
       mo1 = 06
       da1 = 12
       hr1 = 9 !09Z is when the reader reads in the 12Zf00 file
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%datastructime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
 
       gdas_struc(n)%gridchange1 = .true.
       gdas_struc(n)%gridchange2 = .true.
       gdas_struc(n)%gridchange3 = .true.
       gdas_struc(n)%gridchange4 = .true.
       gdas_struc(n)%gridchange5 = .true.
       gdas_struc(n)%gridchange6 = .true.

       gdas_struc(n)%dstrucchange1 = .true.

       ! Setting up weights for Interpolation
       call gdas_reset_interp_input(n, findex, gridDesci)
    enddo

  end subroutine init_GDAS

end module gdas_forcingMod

