!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module WRF_AKdom_forcingMod
!BOP
! !MODULE: WRF_AKdom_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data extracted from WRF output
!  files. 
!
!  The implementation in LDT has the derived data type {\tt WRFAK\_struc} that
!  includes the variables that specify the runtime options.
!  They are described below: 
!  \begin{description}
!  \item[nest\_id]
!    Value of the WRF nest
!  \item[WRFAKdir]
!    Directory containing the input data
!  \item[ts]
!    Frequency in seconds of the forcing data
!  \item[WRFouttime1]
!    The nearest, previous 1 hour instance of the incoming 
!    data (as a real time). 
!  \item[WRFouttime2]
!    The nearest, next 1 hour instance of the incoming 
!    data (as a real time).
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for 
!    temporal interpolation.
!  \end{description}

! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_WRF_AKdom      !defines the native resolution of 
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: WRFAK_struc
!EOP

  type, public    :: WRFAK_type_dec

     integer      :: nest_id
     real         :: ts
     integer      :: nc, nr
     integer      :: mi
     character(len=LDT_CONST_PATH_LEN) :: WRFAKdir
     real*8       :: WRFouttime1,WRFouttime2
     integer      :: findtime1,findtime2
     integer      :: nIter, st_iterid, en_iterid

     real         :: gridDesci(20)

     integer, allocatable  :: n111(:)
     integer, allocatable  :: n121(:)
     integer, allocatable  :: n211(:)
     integer, allocatable  :: n221(:)
     real, allocatable     :: w111(:),w121(:)
     real, allocatable     :: w211(:),w221(:)
     integer, allocatable  :: n112(:,:)
     integer, allocatable  :: n122(:,:)
     integer, allocatable  :: n212(:,:)
     integer, allocatable  :: n222(:,:)
     real, allocatable     :: w112(:,:),w122(:,:)
     real, allocatable     :: w212(:,:),w222(:,:)

     ! Neighbor
     integer, allocatable  :: n113(:)

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

     ! Elevation file:
     character(len=LDT_CONST_PATH_LEN) :: file_wrfelev

  end type WRFAK_type_dec

  type(WRFAK_type_dec), allocatable :: WRFAK_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_WRF_AKdom
! \label{init_WRF_AKdom}
!
! !REVISION HISTORY: 
! 21 Jun 2021; K.R. Arsenault, Updated for different WRF output files
! 
! !INTERFACE:
  subroutine init_WRF_AKdom(findex)
! !USES: 
  use ESMF
  use LDT_coreMod,    only : LDT_rc, LDT_domain
  use LDT_timeMgrMod, only : LDT_update_timestep, LDT_calendar
  use LDT_logMod,     only : LDT_logunit, LDT_endrun, LDT_verify
  use map_utils  

   implicit none
! !ARGUMENTS:  
   integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for LDT output
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readconfig\_WRF_AKdom](\ref{readconfig_WRF_AKdom}) \newline
!    reads the runtime options specified for WRF output data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    
   integer :: n
   integer :: status
   integer :: doy
   real    :: gmt
   type(ESMF_Time)  :: DatastartTime
   type(ESMF_Time)  :: LDTstartTime

   type(proj_info)   :: proj_temp        
   real              :: lat_str,lon_str  
   integer           :: c,r              
   real, allocatable :: xlat(:,:),xlon(:,:) 
    
    allocate(WRFAK_struc(LDT_rc%nnest))
    write(LDT_logunit,*)"[INFO] Initializing WRF-Alaska forcing grid ... "

    LDT_rc%met_nf(findex) = 8

    ! Read inputs from LDT config file:
    call readconfig_WRF_AKdom(findex)

    ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex) = "polar"
    LDT_rc%met_nc(findex) = 659
    LDT_rc%met_nr(findex) = 609

    ! Number of x- and y-dir points:
    WRFAK_struc(:)%nc = 659
    WRFAK_struc(:)%nr = 609

    WRFAK_struc%mi = WRFAK_struc%nc * WRFAK_struc%nr

    do n=1,LDT_rc%nnest

       ! Polar Stereographic projection
       ! See below link for more information on polar stereographic grids:
       !  https://apps.ecmwf.int/codes/grib/format/grib2/templates/3/20
       ! CEN_LAT = 64.f 
       ! CEN_LON = -150.f 
       ! TRUELAT1 = 64.f 
       ! TRUELAT2 = 64.f 
       ! MOAD_CEN_LAT = 64.f 
       ! STAND_LON = -150.f 
       ! POLE_LAT = 90.f 
       ! POLE_LON = 0.f 

       WRFAK_struc(n)%gridDesci = 0
       WRFAK_struc(n)%gridDesci(1) = 5          ! Polar-stereographic
       WRFAK_struc(n)%gridDesci(2) = WRFAK_struc(n)%nc
       WRFAK_struc(n)%gridDesci(3) = WRFAK_struc(n)%nr
       WRFAK_struc(n)%gridDesci(4) = 51.5418    ! latitude of origin -- LL Lat
       WRFAK_struc(n)%gridDesci(5) = -168.175   ! longitude of origin -- LL Lon
       WRFAK_struc(n)%gridDesci(6) = 0          ! ** Matches NAM242 AK grid entry
!       WRFAK_struc(n)%gridDesci(6) = 8          ! Set for polar-stereo in core/LIS_domainMod.F90 (line 2658)
       WRFAK_struc(n)%gridDesci(7) = -150.0     ! Orientation (LoV)
       WRFAK_struc(n)%gridDesci(8) = 4.         ! grid spacing in km
       WRFAK_struc(n)%gridDesci(9) = 4.         ! grid spacing in km
       WRFAK_struc(n)%gridDesci(10) = 64.0      ! true lat1
       WRFAK_struc(n)%gridDesci(11) = -150.0    ! standard long
       WRFAK_struc(n)%gridDesci(13) = 0         ! Regional
       WRFAK_struc(n)%gridDesci(20) = 0         ! ** Matches NAM242 AK grid entry

       LDT_rc%met_gridDesc(findex,1:20) = WRFAK_struc(n)%gridDesci(1:20)

    enddo

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    ! Start date of WRF data availability: Sept 1, 2002
    call ESMF_TimeSet(DatastartTime, yy=2002, &
           mm = 9, dd = 1, h=0, m = 0, s=0, calendar=LDT_calendar, &
           rc=status)
    call LDT_verify(status, &
            'WRF_AKdom_forcingMod: ESMF_TimeSet DatastartTime')
    call ESMF_TimeSet(LDTstartTime, yy=LDT_rc%syr, &
             mm = LDT_rc%smo, dd = LDT_rc%sda, h=LDT_rc%shr, &
             m = LDT_rc%smn, s=LDT_rc%sss, calendar=LDT_calendar, &
             rc=status)
    call LDT_verify(status,'WRF_AKdom_forcingMod: ESMF_TimeSet LDTstartTime')

    ! Check if LDT start date/time begins before Data start date/time :
    if( LDTstartTime < DatastartTime ) then
      write(LDT_logunit,*) &
            '[ERR] LDT start date+time is before WRF-AK data start date+time ...'
      write(LDT_logunit,*) ' WRF-AK data is available from Oct. 1, 2000, and onwards.'
      write(LDT_logunit,*) ' Please update your lis.config file ...'
      call LDT_endrun()
    endif


    do n=1,LDT_rc%nnest

       ! WRF data timestep: 1 hour  == 3600 sec (60 min * 60 sec)
       WRFAK_struc(n)%ts = 60*60
       call LDT_update_timestep(LDT_rc, n, WRFAK_struc(n)%ts)

       allocate(WRFAK_struc(n)%metdata1(1,&
               LDT_rc%met_nf(findex),&
               LDT_rc%ngrid(n)))
       allocate(WRFAK_struc(n)%metdata2(1,&
               LDT_rc%met_nf(findex),&
               LDT_rc%ngrid(n)))

       WRFAK_struc(n)%st_iterid = 1
       WRFAK_struc(n)%en_iterId = 1
       WRFAK_struc(n)%nIter = 1

       WRFAK_struc(n)%metdata1 = 0
       WRFAK_struc(n)%metdata2 = 0

     ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )

         case( "bilinear" )
          allocate(WRFAK_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,WRFAK_struc(n)%gridDesci,&
               WRFAK_struc(n)%n111,WRFAK_struc(n)%n121,&
               WRFAK_struc(n)%n211,WRFAK_struc(n)%n221,&
               WRFAK_struc(n)%w111,WRFAK_struc(n)%w121,&
               WRFAK_struc(n)%w211,WRFAK_struc(n)%w221)

         case( "budget-bilinear" )
          allocate(WRFAK_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,WRFAK_struc(n)%gridDesci,&
               WRFAK_struc(n)%n111,WRFAK_struc(n)%n121,&
               WRFAK_struc(n)%n211,WRFAK_struc(n)%n221,&
               WRFAK_struc(n)%w111,WRFAK_struc(n)%w121,&
               WRFAK_struc(n)%w211,WRFAK_struc(n)%w221)

          allocate(WRFAK_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,WRFAK_struc(n)%gridDesci,&
               WRFAK_struc(n)%n112,WRFAK_struc(n)%n122,&
               WRFAK_struc(n)%n212,WRFAK_struc(n)%n222,&
               WRFAK_struc(n)%w112,WRFAK_struc(n)%w122,&
               WRFAK_struc(n)%w212,WRFAK_struc(n)%w222)

         case( "neighbor" )

          allocate(WRFAK_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call neighbor_interp_input(n,WRFAK_struc(n)%gridDesci,&
                        WRFAK_struc(n)%n113)

       case default
         write(LDT_logunit,*) "[ERR] User-input issue with WRF AK domain forcing ..."
         write(LDT_logunit,*) " -- Currently only supported interpolation options include:"
         write(LDT_logunit,*) "  - bilinear, budget-bilinear, or neighbor - "
         call LDT_endrun
       end select

    enddo  ! End nest loop

  end subroutine init_WRF_AKdom

end module WRF_AKdom_forcingMod
