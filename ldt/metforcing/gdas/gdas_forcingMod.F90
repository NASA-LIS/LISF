!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
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
!  on a quadratic gaussian grid. LDT uses the 00, 03, 06 and as 
!  needed, the 09 forecasts. The forecasts are produced at 6 hr intervals. 
!  The resolution of GDAS forcing varies as follows:
!   
!   upto 2000/1/24          :   T126 (384x190)  grid \newline
!   2001/01/24 - 2002/10/29 :   T170 (512x256)  grid \newline
!   2002/10/29 - 2005/05/31 :   T254 (768x384)  grid \newline
!   2005/05/31 - 2010/07/27 :   T382 (1152x576) grid \newline
!   2010/07/28 - 2015/01/14 :   T574 (1760x880) grid
!   2015/01/14 - onwards    :  T1534 (1760x880) grid
!
!  On 2019/06/12 12Z, GDAS removed precipitation fields from the f00 data
!  files. The data fields in these files are now all instantaneous values.
!  When the reader is using data files after this time, the subroutine will
!  skip the precipitation fields and read in instantaneous radiation data 
!  from the f00 files.
!
!  The implementation in LDT has the derived data type {\tt gdas\_struc} that
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
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
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \end{description}
!
! !USES: 

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
     integer       :: nc, nr   ! GDAS dimensions
     integer       :: nmif
     character*100 :: gdasdir        ! GDAS Forcing Directory
     character*100 :: elevfile

     real*8        :: gdastime1, gdastime2
     real*8        :: griduptime1, griduptime2, griduptime3
     real*8        :: griduptime4, griduptime5, griduptime6
     real*8        :: datastructime1
     logical       :: gridchange1, gridchange2, gridchange3
     logical       :: gridchange4, gridchange5, gridchange6
     logical       :: dstrucchange1
     integer       :: findtime1, findtime2
     integer       :: mi

    !Suffixes 1 are for bilinear
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)
     
    !Suffixes 2 are for conservative
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)

     logical       :: reset_flag
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
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_logMod,     only : LDT_logunit, LDT_endrun
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
    
    implicit none
! !ARGUMENTS: 
    integer,  intent(in)    :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GDAS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
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
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    integer :: n 
    integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    real    :: gridDesci(20)

    allocate(gdas_struc(LDT_rc%nnest))

    write(LDT_logunit,*) "Reading the GDAS forcing data"

!-  Read in config GDAS inputs:
    call readcrd_gdas(findex)
    
    do n=1, LDT_rc%nnest
       gdas_struc(n)%ts = 21600
       call LDT_update_timestep(LDT_rc, n, gdas_struc(n)%ts)
    enddo

    gdas_struc%reset_flag = .false.
    gdas_struc(:)%nmif    = 9 

  ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex) = "gaussian"

 != T126:
    write(LDT_logunit,*)"MSG: Initializing GDAS grid-1 for: 1991-2000 grid (T126)"
    gdas_struc(:)%nc = 384
    gdas_struc(:)%nr = 190
    gridDesci = 0
    gridDesci(1) = 4
    gridDesci(2) = gdas_struc(1)%nc 
    gridDesci(3) = gdas_struc(1)%nr
    gridDesci(4) = 89.277
    gridDesci(5) = 0
    gridDesci(6) = 128
    gridDesci(7) = -89.277
    gridDesci(8) = -0.9375  !-0.938   ! -0.9375 from NOAA Grib table
    gridDesci(9) = 0.9375   !0.938    ! 0.9375 from NOAA Grib table
    gridDesci(10) = 95
    gridDesci(20) = 0

    LDT_rc%met_nc(findex) = gdas_struc(1)%nc
    LDT_rc%met_nr(findex) = gdas_struc(1)%nr
    LDT_rc%met_gridDesc(findex,1:20) = gridDesci(1:20)
   
 != T170:
    write(LDT_logunit,*)"MSG: Initializing GDAS grid-2 for: 2000-2002 grid (T170)"
    LDT_rc%met_nc(findex+1) = 512
    LDT_rc%met_nr(findex+1) = 256
    LDT_rc%met_proj(findex+1)  = "gaussian"
    LDT_rc%met_gridDesc(findex+1,1)  = 4
    LDT_rc%met_gridDesc(findex+1,2)  = 512
    LDT_rc%met_gridDesc(findex+1,3)  = 256
    LDT_rc%met_gridDesc(findex+1,4)  = 89.463
    LDT_rc%met_gridDesc(findex+1,5)  = 0
    LDT_rc%met_gridDesc(findex+1,6)  = 128
    LDT_rc%met_gridDesc(findex+1,7)  = -89.463
    LDT_rc%met_gridDesc(findex+1,8)  = -0.7025 !-0.703  ! -0.7025 from NOAA Grib table
    LDT_rc%met_gridDesc(findex+1,9)  = 0.7025 !0.703   ! 0.7025 from NOAA Grib table
    LDT_rc%met_gridDesc(findex+1,10) = 128
    LDT_rc%met_gridDesc(findex+1,20) = 0

 != T254:
    write(LDT_logunit,*)"MSG: Initializing GDAS grid-3 for: 2002-2005 grid (T254)"
    LDT_rc%met_nc(findex+2) = 768
    LDT_rc%met_nr(findex+2) = 384
    LDT_rc%met_proj(findex+2)  = "gaussian"
    LDT_rc%met_gridDesc(findex+2,1)  = 4
    LDT_rc%met_gridDesc(findex+2,2)  = 768
    LDT_rc%met_gridDesc(findex+2,3)  = 384
    LDT_rc%met_gridDesc(findex+2,4)  = 89.642
    LDT_rc%met_gridDesc(findex+2,5)  = 0
    LDT_rc%met_gridDesc(findex+2,6)  = 128
    LDT_rc%met_gridDesc(findex+2,7)  = -89.642
    LDT_rc%met_gridDesc(findex+2,8)  = -0.46875  !-0.469
    LDT_rc%met_gridDesc(findex+2,9)  = 0.46875   !0.469
    LDT_rc%met_gridDesc(findex+2,10) = 192
    LDT_rc%met_gridDesc(findex+2,20) = 0

 != T382:
    write(LDT_logunit,*)"MSG: Initializing GDAS grid-4 for: 2005-2010 grid (T382)"
    LDT_rc%met_nc(findex+3) = 1152
    LDT_rc%met_nr(findex+3) = 576
    LDT_rc%met_proj(findex+3)  = "gaussian"
    LDT_rc%met_gridDesc(findex+3,1)  = 4
    LDT_rc%met_gridDesc(findex+3,2)  = 1152
    LDT_rc%met_gridDesc(findex+3,3)  = 576
    LDT_rc%met_gridDesc(findex+3,4)  = 89.761
    LDT_rc%met_gridDesc(findex+3,5)  = 0
    LDT_rc%met_gridDesc(findex+3,6)  = 128
    LDT_rc%met_gridDesc(findex+3,7)  = -89.761
    LDT_rc%met_gridDesc(findex+3,8)  = -0.3125  !-0.313
    LDT_rc%met_gridDesc(findex+3,9)  = 0.3125   !0.313
    LDT_rc%met_gridDesc(findex+3,10) = 288
    LDT_rc%met_gridDesc(findex+3,20) = 0

 != T574:
    write(LDT_logunit,*)"MSG: Initializing GDAS grid-5 for: 2010-2015 grid (T574)"
    LDT_rc%met_nc(findex+4) = 1760
    LDT_rc%met_nr(findex+4) = 880
    LDT_rc%met_proj(findex+4)  = "gaussian"
    LDT_rc%met_gridDesc(findex+4,1)  = 4
    LDT_rc%met_gridDesc(findex+4,2)  = 1760
    LDT_rc%met_gridDesc(findex+4,3)  = 880
    LDT_rc%met_gridDesc(findex+4,4)  = 89.844
    LDT_rc%met_gridDesc(findex+4,5)  = 0
    LDT_rc%met_gridDesc(findex+4,6)  = 128
    LDT_rc%met_gridDesc(findex+4,7)  = -89.844
    LDT_rc%met_gridDesc(findex+4,8)  = -0.204545454545455  !-0.205
    LDT_rc%met_gridDesc(findex+4,9)  = 0.204545454545455   !0.205
    LDT_rc%met_gridDesc(findex+4,10) = 440
    LDT_rc%met_gridDesc(findex+4,20) = 0

 != T1534:
    write(LDT_logunit,*)"MSG: Initializing GDAS grid-6 for: 2015-present grid (T1534)"
    LDT_rc%met_nc(findex+5) = 3072
    LDT_rc%met_nr(findex+5) = 1536
    LDT_rc%met_proj(findex+5)  = "gaussian"
    LDT_rc%met_gridDesc(findex+5,1)  = 4
    LDT_rc%met_gridDesc(findex+5,2)  = 3072
    LDT_rc%met_gridDesc(findex+5,3)  = 1536
    LDT_rc%met_gridDesc(findex+5,4)  = 89.91000
    LDT_rc%met_gridDesc(findex+5,5)  = 0.
    LDT_rc%met_gridDesc(findex+5,6)  = 128
    LDT_rc%met_gridDesc(findex+5,7)  = -89.91000
    LDT_rc%met_gridDesc(findex+5,8)  = -0.1171875
    LDT_rc%met_gridDesc(findex+5,9)  = 0.1171875
    LDT_rc%met_gridDesc(findex+5,10) = 768.0
    LDT_rc%met_gridDesc(findex+5,20) = 0

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 9    ! number of met variables in GDAS forcing
    LDT_rc%met_ts(findex) = 21600
    LDT_rc%met_zterp(findex) = .true.

    do n=1,LDT_rc%nnest

       gdas_struc(n)%findtime1 = 0 
       gdas_struc(n)%findtime2 = 0 

       gdas_struc(n)%mi = gdas_struc(n)%nc*gdas_struc(n)%nr
       
     ! T126
       yr1 = 2000
       mo1 = 01
       da1 = 24
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LDT_date2time( gdas_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
       
     ! T170
       yr1 = 2002     !grid update time ~ 0.469
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LDT_date2time(gdas_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2005     !grid update time ~ 0.313
       mo1 = 05
       da1 = 31
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LDT_date2time(gdas_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2010     !grid update time ~ 0.205
       mo1 = 07
       da1 = 28
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LDT_date2time(gdas_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2015     !grid update time ~ 0.117
       mo1 = 01
       da1 = 14
       hr1 = 6
       mn1 = 0; ss1 = 0
       call LDT_date2time(gdas_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       ! Set time for f00 data structure change
       yr1 = 2019
       mo1 = 06
       da1 = 12
       hr1 = 9 !09Z is when the reader reads in the 12Zf00 file
       mn1 = 0; ss1 = 0
       call LDT_date2time(gdas_struc(n)%datastructime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
       
       gdas_struc(n)%gridchange1 = .true.
       gdas_struc(n)%gridchange2 = .true.
       gdas_struc(n)%gridchange3 = .true.
       gdas_struc(n)%gridchange4 = .true.
       gdas_struc(n)%gridchange5 = .true.
       gdas_struc(n)%gridchange6 = .true.

       gdas_struc(n)%dstrucchange1 = .true.
     ! Setting up weights for Interpolation
      
       select case( LDT_rc%met_gridtransform(findex) ) 

         case( "bilinear" )
          allocate(gdas_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               gdas_struc(n)%n111,gdas_struc(n)%n121,&
               gdas_struc(n)%n211,gdas_struc(n)%n221,&
               gdas_struc(n)%w111,gdas_struc(n)%w121,&
               gdas_struc(n)%w211,gdas_struc(n)%w221)

         case( "budget-bilinear" )

          allocate(gdas_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               gdas_struc(n)%n111,gdas_struc(n)%n121,&
               gdas_struc(n)%n211,gdas_struc(n)%n221,&
               gdas_struc(n)%w111,gdas_struc(n)%w121,&
               gdas_struc(n)%w211,gdas_struc(n)%w221)

          allocate(gdas_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gdas_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gdas_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gdas_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gdas_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gdas_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gdas_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gdas_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               gdas_struc(n)%n112,gdas_struc(n)%n122,&
               gdas_struc(n)%n212,gdas_struc(n)%n222,&
               gdas_struc(n)%w112,gdas_struc(n)%w122,&
               gdas_struc(n)%w212,gdas_struc(n)%w222)
       case default
          write(LDT_logunit,*) 'The specified spatial interpolation option not'
          write(LDT_logunit,*) 'supported for GDAS..'
          write(LDT_logunit,*) 'program stopping..'
          call LDT_endrun()
      end select
    enddo

  end subroutine init_GDAS

end module gdas_forcingMod

