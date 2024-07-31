!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
 module gefs_forcingMod
   use ESMF
!
!BOP
! !MODULE: gefs_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of GEFS ensemble forecast data used as forcing
!  within LIS. 
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GEFS        ! defines the native input resolution 

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gefs_struc
!EOP

  type, public ::  gefs_type_dec

    integer            :: max_ens_members
    character(len=LIS_CONST_PATH_LEN) :: gefs_dir
    character*20       :: gefs_fcsttype
    character*20       :: gefs_runmode
    character*20       :: gefs_proj
    character*20       :: gefs_preslevel
    real               :: gefs_res
    
    real               :: ts
    real*8             :: fcsttime1,fcsttime2

    type(ESMF_Time)    :: initTime
    integer            :: init_yr, init_mo, init_da, init_hr
    integer            :: gribrec, fcst_hour

    real               :: gridDesc(50)
    integer            :: nc, nr
    integer            :: mi
   
    ! Bilinear weights
    integer, allocatable  :: n111(:)
    integer, allocatable  :: n121(:)
    integer, allocatable  :: n211(:)
    integer, allocatable  :: n221(:)
    real, allocatable     :: w111(:),w121(:)
    real, allocatable     :: w211(:),w221(:)

    ! Budget-bilinear weights ("Conserve")
    integer, allocatable  :: n112(:,:)
    integer, allocatable  :: n122(:,:)
    integer, allocatable  :: n212(:,:)
    integer, allocatable  :: n222(:,:)
    real, allocatable     :: w112(:,:),w122(:,:)
    real, allocatable     :: w212(:,:),w222(:,:)

    real, allocatable  :: metdata1(:,:,:)
    real, allocatable  :: metdata2(:,:,:)

    integer            :: findtime1, findtime2

  end type gefs_type_dec

  type(gefs_type_dec), allocatable :: gefs_struc(:)
  logical, public :: gefs_initialized = .false.

!EOP
contains

!BOP
!
! !ROUTINE: init_GEFS
! \label{init_GEFS}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
! 28 Jan 2021: Sujay Kumar; Updated for GEFS operational data
! 
! !INTERFACE:
  subroutine init_GEFS(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GEFS
!  forecast data.
!  The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}). 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gefs](\ref{readcrd_gefs}) \newline
!     reads the runtime options specified for GEFS forecast data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!  \end{description}
!EOP

    integer :: n
    integer :: rc
    real :: lonEnd
    write(LIS_logunit,*) "[INFO] Initializing the GEFS forecast inputs "

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GEFS forecast forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in ESP-forecast mode.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
 
    ! 8 - key met fields:
    LIS_rc%met_nf(findex) = 8

    if(.not.gefs_initialized) then
       gefs_initialized = .true.

       allocate(gefs_struc(LIS_rc%nnest))

       ! Read LIS-GEFS config file entries:
       call readcrd_gefs()

       ! Allocate and initialize GEFS metforcing data structures:
       do n=1, LIS_rc%nnest
          LIS_rc%met_nensem(findex) = gefs_struc(n)%max_ens_members

          allocate(gefs_struc(n)%metdata1(LIS_rc%met_nf(findex),&
               gefs_struc(n)%max_ens_members,&
               LIS_rc%ngrid(n)))
          allocate(gefs_struc(n)%metdata2(LIS_rc%met_nf(findex),&
               gefs_struc(n)%max_ens_members,&
               LIS_rc%ngrid(n)))

          gefs_struc(n)%metdata1 = 0
          gefs_struc(n)%metdata2 = 0

          gefs_struc(n)%findtime1 = 0
          gefs_struc(n)%findtime2 = 0
       enddo

       ! Set time and gridDesc arrays for different GEFS products:
       do n=1,LIS_rc%nnest

          gefs_struc(n)%gridDesc  = 0
          ! Reforecast product:
          if( gefs_struc(n)%gefs_fcsttype .eq. "Reforecast2" ) then

            ! 6-hourly timestep:
            gefs_struc(n)%ts = 21600   ! 6 * 3600sec
            call LIS_update_timestep(LIS_rc, n, gefs_struc(n)%ts)

            ! Check if starting hour of LIS run matches 00Z:
            ! MAY CHANGE THIS LATER TO ACCOUNT FOR STARTING FORECAST TIME
            ! AT 23:45Z or 23:30Z THE DAY BEFORE TO ACCOUNT FOR THE 
            ! FULL FIRST 3-HR or 6-HR OUTPUT WINDOW ...
            if( LIS_rc%shr .ne. 0 ) then
               write(LIS_logunit,*) "[ERR] GEFS 'Reforecast' forecast type begins"
               write(LIS_logunit,*) "[ERR] at 00Z for a forecast window, so the "
               write(LIS_logunit,*) "[ERR] 'Starting hour:' should be set to 0 in"
               write(LIS_logunit,*) "[ERR]  your lis.config file.."
               call LIS_endrun()
            endif
         
            ! Initialize the forecast initial date-time and grib record:
            gefs_struc(n)%init_yr = LIS_rc%syr
            gefs_struc(n)%init_mo = LIS_rc%smo
            gefs_struc(n)%init_da = LIS_rc%sda
            gefs_struc(n)%init_hr = LIS_rc%shr
            
            gefs_struc(n)%gribrec = 1   ! Record is different per field
            gefs_struc(n)%fcst_hour = 0

            if( gefs_struc(n)%gefs_proj .eq. "latlon" ) then
              if( gefs_struc(n)%gefs_res .eq. 0.25 ) then
                 gefs_struc(n)%nc = 1440
                 gefs_struc(n)%nr = 721
                 lonEnd = 179.75
              else
                 gefs_struc(n)%nc = 360
                 gefs_struc(n)%nr = 181
                 lonEnd = 178.75
              endif

              gefs_struc(n)%gridDesc(1) = 0
              gefs_struc(n)%gridDesc(2) = gefs_struc(n)%nc
              gefs_struc(n)%gridDesc(3) = gefs_struc(n)%nr
              gefs_struc(n)%gridDesc(4) = 90.0    ! original GEFS
              gefs_struc(n)%gridDesc(5) = -180.0
!              gefs_struc(n)%gridDesc(5) = 0.0     ! original GEFS
              gefs_struc(n)%gridDesc(6) = 128
              gefs_struc(n)%gridDesc(7) = -90.0   ! original GEFS
              gefs_struc(n)%gridDesc(8) = lonEnd
!              gefs_struc(n)%gridDesc(8) = 359.00  ! original GEFS
              gefs_struc(n)%gridDesc(9) =  gefs_struc(n)%gefs_res
              gefs_struc(n)%gridDesc(10) = gefs_struc(n)%gefs_res
              gefs_struc(n)%gridDesc(20) = 64  ! -180 to 180?
!              gefs_struc(n)%gridDesc(20) = 0   ! for 0 to 360?

            elseif( gefs_struc(n)%gefs_proj .eq. "gaussian" ) then
              write(LIS_logunit,*) "[ERR] GEFS Reforecast2 -- 'gaussian' projection"
              write(LIS_logunit,*) "[ERR]  not currently supported ... "
              call LIS_endrun()
            else
              write(LIS_logunit,*) "[ERR] GEFS Reforecast2 projection selected"
              write(LIS_logunit,*) "[ERR]  is not available. Only 'latlon' is "
              write(LIS_logunit,*) "[ERR]  currently supported. "
              call LIS_endrun()
            endif

          elseif( gefs_struc(n)%gefs_fcsttype .eq. "Operational" ) then

             ! Initialize the forecast initial date-time and grib record:
             gefs_struc(n)%init_yr = LIS_rc%syr
             gefs_struc(n)%init_mo = LIS_rc%smo
             gefs_struc(n)%init_da = LIS_rc%sda
             gefs_struc(n)%init_hr = LIS_rc%shr
             
             call ESMF_TimeSet(gefs_struc(n)%initTime, &
                  yy = gefs_struc(n)%init_yr, &
                  mm = gefs_struc(n)%init_mo, &
                  dd = gefs_struc(n)%init_da, &
                  h  = gefs_struc(n)%init_hr, &
                  m  = 0,&
                  s  = 0,&
                  rc=rc) 
             call LIS_verify(rc,'ESMF_TimeSet failed in init_gefs')
             
             if( gefs_struc(n)%gefs_proj .eq. "latlon" ) then
                if( gefs_struc(n)%gefs_res .eq. 0.25 ) then
                   gefs_struc(n)%nc = 1440
                   gefs_struc(n)%nr = 721
                   lonEnd = 179.75
                else
                   gefs_struc(n)%nc = 720
                   gefs_struc(n)%nr = 361
                   lonEnd = 179.5
                endif
                
                gefs_struc(n)%gridDesc(1) = 0
                gefs_struc(n)%gridDesc(2) = gefs_struc(n)%nc
                gefs_struc(n)%gridDesc(3) = gefs_struc(n)%nr
                gefs_struc(n)%gridDesc(4) = 90.0    ! original GEFS
                gefs_struc(n)%gridDesc(5) = -180.0
!              gefs_struc(n)%gridDesc(5) = 0.0     ! original GEFS
                gefs_struc(n)%gridDesc(6) = 128
                gefs_struc(n)%gridDesc(7) = -90.0   ! original GEFS
                gefs_struc(n)%gridDesc(8) = lonEnd
!              gefs_struc(n)%gridDesc(8) = 359.00  ! original GEFS
                gefs_struc(n)%gridDesc(9) =  gefs_struc(n)%gefs_res
                gefs_struc(n)%gridDesc(10) = gefs_struc(n)%gefs_res
                gefs_struc(n)%gridDesc(20) = 64  ! -180 to 180?
!              gefs_struc(n)%gridDesc(20) = 0   ! for 0 to 360?

             elseif( gefs_struc(n)%gefs_proj .eq. "gaussian" ) then
                write(LIS_logunit,*) "[ERR] GEFS Reforecast2 -- 'gaussian' projection"
                write(LIS_logunit,*) "[ERR]  not currently supported ... "
                call LIS_endrun()
             else
                write(LIS_logunit,*) "[ERR] GEFS Reforecast2 projection selected"
                write(LIS_logunit,*) "[ERR]  is not available. Only 'latlon' is "
                write(LIS_logunit,*) "[ERR]  currently supported. "
                call LIS_endrun()
             endif
          endif

          write(LIS_logunit,*) "[INFO] GEFS number of cols: ", gefs_struc(n)%nc
          write(LIS_logunit,*) "[INFO] GEFS number of rows: ", gefs_struc(n)%nr

          gefs_struc(n)%mi = gefs_struc(n)%nc*gefs_struc(n)%nr

          !Setting up weights for Interpolation
          if( LIS_rc%met_interp(findex).eq."bilinear" ) then

             ! Initialize bilinear weights:
             allocate(gefs_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

             call bilinear_interp_input(n,gefs_struc(n)%gridDesc(:),&
                  gefs_struc(n)%n111,gefs_struc(n)%n121,&
                  gefs_struc(n)%n211,&
                  gefs_struc(n)%n221,gefs_struc(n)%w111,&
                  gefs_struc(n)%w121,&
                  gefs_struc(n)%w211,gefs_struc(n)%w221)

          elseif( LIS_rc%met_interp(findex).eq."budget-bilinear" ) then

             ! Initialize bilinear weights:
             allocate(gefs_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gefs_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

             call bilinear_interp_input(n,gefs_struc(n)%gridDesc(:),&
                  gefs_struc(n)%n111,gefs_struc(n)%n121,&
                  gefs_struc(n)%n211,&
                  gefs_struc(n)%n221,gefs_struc(n)%w111,&
                  gefs_struc(n)%w121,&
                  gefs_struc(n)%w211,gefs_struc(n)%w221)

             ! Initialize budget-bilinear weights:
             allocate(gefs_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gefs_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gefs_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gefs_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gefs_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gefs_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gefs_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gefs_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

             call conserv_interp_input(n, gefs_struc(n)%gridDesc(:), &
                  gefs_struc(n)%n112,gefs_struc(n)%n122, &
                  gefs_struc(n)%n212,gefs_struc(n)%n222, &
                  gefs_struc(n)%w112,gefs_struc(n)%w122, &
                  gefs_struc(n)%w212,gefs_struc(n)%w222)

          else
             write(LIS_logunit,*) &
                  '[ERR] Interpolation option not specified for GEFS ...'
             write(LIS_logunit,*) '[ERR] Program stopping ...'
             call LIS_endrun()
          endif
       enddo
    end if

  end subroutine init_GEFS

end module gefs_forcingMod

