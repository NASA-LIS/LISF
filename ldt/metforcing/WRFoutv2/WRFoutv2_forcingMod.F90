!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module WRFoutv2_forcingMod
!BOP
! !MODULE: WRFoutv2_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data extracted from WRF output
!  files.  Here WRF output files are consider input data.
!
!  The implementation in LDT has the derived data type {\tt WRFoutv2\_struc} that
!  includes the variables that specify the runtime options.
!  They are described below: 
!  \begin{description}
!  \item[nest\_id]
!    Value of the WRF nest
!  \item[WRFoutv2dir]
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
  public :: init_WRFoutv2      !defines the native resolution of 
                               !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: WRFoutv2_struc
!EOP

  type, public    :: WRFoutv2_type_dec

     integer      :: nest_id
     real         :: ts
     integer      :: nc, nr
     integer      :: mi
     character(len=LDT_CONST_PATH_LEN) :: WRFoutv2dir
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

  end type WRFoutv2_type_dec

  type(WRFoutv2_type_dec), allocatable :: WRFoutv2_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_WRFoutv2
! \label{init_WRFoutv2}
!
! !REVISION HISTORY: 
! 20 Nov 2020; K.R. Arsenault, Updated for different WRF output files
! 
! !INTERFACE:
  subroutine init_WRFoutv2(findex)
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
!  Note that no interpolation is performed for this forcing. 
!  The data is expected to be in the same map projection
!  and resolution as that of the current LDT run.
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readconfig\_WRFoutv2](\ref{readconfig_WRFoutv2}) \newline
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
    
    allocate(WRFoutv2_struc(LDT_rc%nnest))
    write(LDT_logunit,*)"[INFO] Initializing WRFoutv2 forcing grid ... "

    ! Read inputs from LDT config file:
    call readconfig_WRFoutv2(findex)

  ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex) = "lambert"
    LDT_rc%met_nc(findex) = 1359
    LDT_rc%met_nr(findex) = 1015

    ! Number of x- and y-dir points:
    WRFoutv2_struc%nc = 1359
    WRFoutv2_struc%nr = 1015

    WRFoutv2_struc%mi = WRFoutv2_struc%nc*WRFoutv2_struc%nr

    ! Lambert CC grid ...
    WRFoutv2_struc%gridDesci(1) = 3          ! Lambert conic conformal grid
    WRFoutv2_struc%gridDesci(2) = WRFoutv2_struc%nc
    WRFoutv2_struc%gridDesci(3) = WRFoutv2_struc%nr
    WRFoutv2_struc%gridDesci(4) = 18.1363    ! latitude of origin -- LL Lat
    WRFoutv2_struc%gridDesci(5) = -122.8839  ! longitude of origin -- LL Lon
    WRFoutv2_struc%gridDesci(6) = 8          ! Set for Lambert in core/LDT_domainMod.F90 (line 2658)
    WRFoutv2_struc%gridDesci(7) = 50.0       ! true lat2
    WRFoutv2_struc%gridDesci(8) = 4.         ! grid spacing in km
    WRFoutv2_struc%gridDesci(9) = 4.         ! grid spacing in km
    WRFoutv2_struc%gridDesci(10) = 28.0      ! true lat1
    WRFoutv2_struc%gridDesci(11) = -98.0     ! standard long
    WRFoutv2_struc%gridDesci(20) = 0.0

    LDT_rc%met_gridDesc(findex,1:20) = WRFoutv2_struc(1)%gridDesci(1:20)


 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 8

    ! Start date of WRF data availability: October 1, 2000
    call ESMF_TimeSet(DatastartTime, yy=2000, &
           mm = 10, dd = 1, h=0, m = 0, s=0, calendar=LDT_calendar, &
           rc=status)
    call LDT_verify(status, &
            'WRFoutv2_forcingMod: ESMF_TimeSet DatastartTime')
    call ESMF_TimeSet(LDTstartTime, yy=LDT_rc%syr, &
             mm = LDT_rc%smo, dd = LDT_rc%sda, h=LDT_rc%shr, &
             m = LDT_rc%smn, s=LDT_rc%sss, calendar=LDT_calendar, &
             rc=status)
    call LDT_verify(status,'WRFoutv2_forcingMod: ESMF_TimeSet LDTstartTime')

    ! Check if LDT start date/time begins before Data start date/time :
    if( LDTstartTime < DatastartTime ) then
      write(LDT_logunit,*) &
            '[ERR] LDT start date+time is before WRF data start date+time ...'
      write(LDT_logunit,*) ' WRF data is available from Oct. 1, 2000, and onwards.'
      write(LDT_logunit,*) ' Please update your ldt.config file ...'
      call LDT_endrun()
    endif

    do n=1,LDT_rc%nnest

       ! WRF output timestep: 1 hour  == 3600 sec (60 min * 60 sec)
       WRFoutv2_struc(n)%ts = 60*60
       call LDT_update_timestep(LDT_rc, n, WRFoutv2_struc(n)%ts)

       ! Regular retrospective or non-forecast mode:
       allocate(WRFoutv2_struc(n)%metdata1(1,&
            LDT_rc%met_nf(findex),&
            LDT_rc%ngrid(n)))
       allocate(WRFoutv2_struc(n)%metdata2(1,&
            LDT_rc%met_nf(findex),&
            LDT_rc%ngrid(n)))

       WRFoutv2_struc(n)%st_iterid = 1
       WRFoutv2_struc(n)%en_iterId = 1
       WRFoutv2_struc(n)%nIter = 1

       WRFoutv2_struc(n)%metdata1 = 0
       WRFoutv2_struc(n)%metdata2 = 0

       ! Lambert CC grid ...
       WRFoutv2_struc(n)%gridDesci(1) = 3          ! Lambert conic conformal grid
       WRFoutv2_struc(n)%gridDesci(2) = WRFoutv2_struc(n)%nc
       WRFoutv2_struc(n)%gridDesci(3) = WRFoutv2_struc(n)%nr
       WRFoutv2_struc(n)%gridDesci(4) = 18.1363    ! latitude of origin -- LL Lat
       WRFoutv2_struc(n)%gridDesci(5) = -122.8839  ! longitude of origin -- LL Lon
       WRFoutv2_struc(n)%gridDesci(6) = 8          ! Set for Lambert in core/LDT_domainMod.F90 (line 2658)
       WRFoutv2_struc(n)%gridDesci(7) = 50.0       ! true lat2
       WRFoutv2_struc(n)%gridDesci(8) = 4.         ! grid spacing in km
       WRFoutv2_struc(n)%gridDesci(9) = 4.         ! grid spacing in km
       WRFoutv2_struc(n)%gridDesci(10) = 28.0      ! true lat1
       WRFoutv2_struc(n)%gridDesci(11) = -98.0     ! standard long
       WRFoutv2_struc(n)%gridDesci(20) = 0.0

       ! CEN_LAT = 39.70001f 
       ! CEN_LON = -98.f 
       ! TRUELAT1 = 28.f
       ! TRUELAT2 = 50.f 
       ! MOAD_CEN_LAT = 39.70001f 
       ! STAND_LON = -98.f 

     ! CHECK FULL DOMAIN LAT/LONG:
     call map_set( PROJ_LC, &
              WRFoutv2_struc(n)%gridDesci(4), WRFoutv2_struc(n)%gridDesci(5),&
              4000.0, WRFoutv2_struc(n)%gridDesci(11),&
              WRFoutv2_struc(n)%gridDesci(10), WRFoutv2_struc(n)%gridDesci(7),&
              WRFoutv2_struc(n)%nc, WRFoutv2_struc(n)%nr, proj_temp)

       allocate( xlat(WRFoutv2_struc(1)%nc,WRFoutv2_struc(1)%nr) ) ! KRA
       allocate( xlon(WRFoutv2_struc(1)%nc,WRFoutv2_struc(1)%nr) ) ! KRA

       do r=1,WRFoutv2_struc(1)%nr
          do c=1,WRFoutv2_struc(1)%nc
             call ij_to_latlon(proj_temp,&
                        real(c), real(r),&
                        xlat(c,r), xlon(c,r))
          enddo
       enddo


     ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )

         case( "bilinear" )
          allocate(WRFoutv2_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,WRFoutv2_struc(n)%gridDesci,&
               WRFoutv2_struc(n)%n111,WRFoutv2_struc(n)%n121,&
               WRFoutv2_struc(n)%n211,WRFoutv2_struc(n)%n221,&
               WRFoutv2_struc(n)%w111,WRFoutv2_struc(n)%w121,&
               WRFoutv2_struc(n)%w211,WRFoutv2_struc(n)%w221)

         case( "budget-bilinear" )
          allocate(WRFoutv2_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WRFoutv2_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,WRFoutv2_struc(n)%gridDesci,&
               WRFoutv2_struc(n)%n111,WRFoutv2_struc(n)%n121,&
               WRFoutv2_struc(n)%n211,WRFoutv2_struc(n)%n221,&
               WRFoutv2_struc(n)%w111,WRFoutv2_struc(n)%w121,&
               WRFoutv2_struc(n)%w211,WRFoutv2_struc(n)%w221)

          allocate(WRFoutv2_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFoutv2_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFoutv2_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFoutv2_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFoutv2_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFoutv2_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFoutv2_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(WRFoutv2_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,WRFoutv2_struc(n)%gridDesci,&
               WRFoutv2_struc(n)%n112,WRFoutv2_struc(n)%n122,&
               WRFoutv2_struc(n)%n212,WRFoutv2_struc(n)%n222,&
               WRFoutv2_struc(n)%w112,WRFoutv2_struc(n)%w122,&
               WRFoutv2_struc(n)%w212,WRFoutv2_struc(n)%w222)

         case( "neighbor" )

          allocate(WRFoutv2_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call neighbor_interp_input(n,WRFoutv2_struc(n)%gridDesci,&
                        WRFoutv2_struc(n)%n113)

       case default
         write(LDT_logunit,*) "[ERR] User-input issue with WRFoutv2 forcing ..."
         write(LDT_logunit,*) " -- Currently only supported interpolation options include:"
         write(LDT_logunit,*) "  - bilinear, budget-bilinear, or neighbor - "
         call LDT_endrun
       end select

    enddo

  end subroutine init_WRFoutv2

end module WRFoutv2_forcingMod
