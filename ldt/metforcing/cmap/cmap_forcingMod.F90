!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module cmap_forcingMod
!BOP
! !MODULE: cmap_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  Climate Prediction Center (CPC)'s Merged Analysis of Precipitation
!  (CMAP). CMAP merges gauge measurements and satellite estimates
!  including GPI, OPI, SSM/I to produce a global 2.5 degree pentad
!  precipitation analysis. The 6-hourly product obtained by disaggregating 
!  CMAP using the GDAS forcing is used in this implementation. 
!
!   upto 2000/1/24          :   T126 (384x190)  grid
!   2001/1/24  - 2002/10/29 :   T170 (512x256)  grid
!   2002/10/29 - 2005/5/31  :   T254 (768x384)  grid
!   2005/5/31  - 2012/9/30  :   T382 (1152x576) grid ~~ CMAP only
!   2012/10/01 onwards      :   T574 (1760x880) grid ~~ CMAP only
!  Original GDAS grid change:
!   2010/7/28  onwards      :   T574 (1760x880) grid ~~ Original GDAS
! 
!  The implementation in LDT has the derived data type {\tt cmap\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation.
! 
!  They are desribed below: 
! \begin{description}
!  \item[cmapdir]
!    Directory containing the input data
!  \item[cmaptime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution to T170
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
  public :: init_CMAP      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: cmap_struc

!EOP
 
  type, public ::  cmap_type_dec
     real                 :: ts
     integer              :: nc
     integer              :: nr
     integer              :: mi
     character(len=LDT_CONST_PATH_LEN) :: cmapdir  
     real*8               :: cmaptime
     real*8               :: griduptime1
     real*8               :: griduptime2
     real*8               :: griduptime3
     real*8               :: griduptime4
     real*8               :: griduptime5
     logical              :: gridchange1
     logical              :: gridchange2
     logical              :: gridchange3
     logical              :: gridchange4
     logical              :: gridchange5
     integer, allocatable :: n112(:,:)
     integer, allocatable :: n122(:,:)
     integer, allocatable :: n212(:,:)
     integer, allocatable :: n222(:,:)
     real,    allocatable :: w112(:,:),w122(:,:)
     real,    allocatable :: w212(:,:),w222(:,:)
  end type cmap_type_dec

  type(cmap_type_dec), allocatable :: cmap_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_CMAP
! \label{init_CMAP}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 10Oct2014: KR Arsenault: Updated in LDT
! 
! !INTERFACE:
  subroutine init_CMAP(findex)

! !USES: 
    use ESMF
    use LDT_coreMod, only : LDT_rc
    use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep

    implicit none
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for CMAP
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcmapcrd](\ref{readcmapcrd}) \newline
!     reads the runtime options specified for CMAP data
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    integer :: n 
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: status
    real    :: gridDesci(20)
!_________________________________________________

    allocate(cmap_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing CMAP forcing grid ... "

!- Read in LDT Config file entries:
    call readcrd_cmap()

    LDT_rc%met_nf(findex) = 2
    LDT_rc%met_ts(findex) = 21600
    LDT_rc%met_validhr(findex) = 0   ! CMAP/GDAS valid at 0,6,12,18Z
    LDT_rc%met_proj(findex) = "gaussian"

!- Initialize first available grid domain/res for CMAP/GDAS:
!   2001/1/24 - 2002/10/29 : T170 (512x256)  grid
    cmap_struc%nc = 512
    cmap_struc%nr = 256
    LDT_rc%met_nc(findex) = cmap_struc(1)%nc
    LDT_rc%met_nr(findex) = cmap_struc(1)%nr

 !- CMAP/GDAS Grid description:
    gridDesci(1)  = 4 
    gridDesci(2)  = cmap_struc(1)%nc
    gridDesci(3)  = cmap_struc(1)%nr
    gridDesci(4)  = 89.463
    gridDesci(5)  = 0
    gridDesci(6)  = 128
    gridDesci(7)  = -89.463
    gridDesci(8)  = -0.703
    gridDesci(9)  = 0.703
    gridDesci(10) = 128
    gridDesci(20) = 64

    LDT_rc%met_gridDesc(findex,:) = gridDesci(:)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

 !- Initialize and set timestep:
    do n=1, LDT_rc%nnest
       cmap_struc(n)%ts = 21600
       call LDT_update_timestep(LDT_rc, n, cmap_struc(n)%ts)
    enddo

 !- Set interp arrays for reinterpolation later:
    do n = 1, LDT_rc%nnest

       select case( LDT_rc%met_gridtransform(findex) ) 

        case( "budget-bilinear" ) 
          allocate(cmap_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(cmap_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(cmap_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(cmap_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(cmap_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(cmap_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(cmap_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(cmap_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n, gridDesci(:),&
               cmap_struc(n)%n112, cmap_struc(n)%n122,&
               cmap_struc(n)%n212, cmap_struc(n)%n222,&
               cmap_struc(n)%w112, cmap_struc(n)%w122,&
               cmap_struc(n)%w212, cmap_struc(n)%w222)

        case default
          write(LDT_logunit,*) "CMAP precip implementation only supports "
          write(LDT_logunit,*) "conservative implementation .... "
          write(LDT_logunit,*) "Program stopping ..."
!           call LDT_endrun()
          stop
       end select

!------------------------------------------------------------------------ 
! CMAP product start time to use Jon G.'s data. 1979010100 for 00Z-06Z.
!------------------------------------------------------------------------

     ! This grid is good for some time in the 1990's.
     ! Look up the exact dates.
       yr1 = 1991     !grid update time
       mo1 = 01
       da1 = 01
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LDT_date2time(cmap_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,&
               da1,hr1,mn1,ss1 )

!   2001/1/24 - 2002/10/29 : T170 (512x256)  grid
       yr1 = 2000     !grid update time
       mo1 = 01
       da1 = 24
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LDT_date2time(cmap_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! ==  2002/10/29 - 2005/5/31 : T254 (768x384)  grid
       yr1 = 2002     !grid update time
       mo1 = 10
       da1 = 29
       hr1 = 12      ! -- ORIGINAL HOUR SET IN CODE
!       hr1 = 18     ! When actual grid change occurs
       mn1 = 0; ss1 = 0
       call LDT_date2time(cmap_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )
       
! == 2005/5/31 - 2012/9/30  : T382 (1152x576) grid ~~ CMAP only
       yr1 = 2005     !grid update time
       mo1 = 05
       da1 = 31
       hr1 = 0      ! -- ORIGINAL HOUR SET IN CODE
!       hr1 = 6      ! When actual grid change occurs
       mn1 = 0; ss1 = 0
       call LDT_date2time(cmap_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! == 2012/10/01 onwards  : T574 (1760x880) grid ~~ CMAP only
       yr1 = 2012     !grid update time
       mo1 = 09
       da1 = 30
       hr1 = 18      ! -- ORIGINAL HOUR SET IN CODE
!       hr1 = 0      ! When actual grid change occurs
       mn1 = 0; ss1 = 0
       call LDT_date2time(cmap_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )


       cmap_struc(n)%gridchange1 = .true.
       cmap_struc(n)%gridchange2 = .true.
       cmap_struc(n)%gridchange3 = .true.
       cmap_struc(n)%gridchange4 = .true.
       cmap_struc(n)%gridchange5 = .true.

    enddo

  end subroutine init_CMAP

end module cmap_forcingMod
