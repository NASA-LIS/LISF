!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module chirps2_forcingMod
!BOP
! !MODULE: chirps2_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation product derived from 
! 
!  The implementation in LIS has the derived data type {\tt chirps2\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[directory]
!    Directory containing the input data
!  \item[griduptime1]
!    The time to switch the input resolution or file naming convention
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Array containing the weights of the input grid
!    for each grid point in LIS, for neighbor search.
!  \end{description}
!
! !REVISION HISTORY:
!  Yudong Tian,  10/26/2010: Updated to support reading raw .gz files 
!  KR Arsenault, 07/06/2015: Updated to support latest V7 files

! !USES:
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_chirps2      !defines the native resolution of
                                 !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: chirps2_struc

!EOP

  type, public ::  chirps2_type_dec
     real               :: xres, yres
     real               :: ts
     integer            :: nc
     integer            :: nr
     integer            :: mi
     character(len=LIS_CONST_PATH_LEN) :: directory  
     real*8             :: chirpstime1, chirpstime2
     logical            :: reset_flag
     integer            :: start_nc, start_nr
     real               :: gridDesc(50) ! EMK Corrected array length

   ! "Interp" or spatial transform arrays:
   ! Budget-bilinear:
     integer, allocatable  :: n112(:,:)
     integer, allocatable  :: n122(:,:)
     integer, allocatable  :: n212(:,:)
     integer, allocatable  :: n222(:,:)
     real,    allocatable  :: w112(:,:),w122(:,:)
     real,    allocatable  :: w212(:,:),w222(:,:)
   ! Neighbor:
     integer, allocatable  :: n113(:)
   ! Average:
     integer, allocatable  :: n111(:)
     
     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type chirps2_type_dec
  
  type(chirps2_type_dec), allocatable :: chirps2_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_chirps2
! \label{init_chirps2}
! 
!
! !REVISION HISTORY: 
! 11Dec2015: KR Arsenault; Initial Specification
! 
! !INTERFACE:
  subroutine init_chirps2(findex)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_domain, &
                            LIS_isatAfinerResolution
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep, &
                               LIS_tick, LIS_parseTimeString, &
                               LIS_registerAlarm 
    use LIS_logMod,  only : LIS_logunit, LIS_endrun, &
                            LIS_getNextUnitNumber, &
                            LIS_releaseUnitNumber, LIS_verify 
    use LIS_gridmappingMod
    use LIS_FORC_AttributesMod
    use LIS_forecastMod

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for CHIRPS 2.0
!  data. The grid description arrays are based on the netcdf files
!  (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readconfig\_chirps2](\ref{readconfig_chirps2}) \newline
!     reads the runtime options specified for CHIRPS 2.0 data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[neighbor\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
    integer :: n 
    integer :: rc             
    integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: ts1,ts2 

    integer :: LIS_syr,LIS_smo,LIS_sda,LIS_shr,LIS_smn,LIS_sss 
    real*8  :: FirstDateTime
    real*8  :: LIS_StartTime  
    character (len=10):: time 

    integer :: subpnc, subpnr, glpnc, glpnr
    real     :: sub_gridDesci(50) ! EMK Corrected array length
    real    :: fulldomain_gridDesci(50) ! EMK Corrected array length
    integer, allocatable  :: lat_line(:,:), lon_line(:,:)

! _________________________________________________________________

    allocate(chirps2_struc(LIS_rc%nnest))
    chirps2_struc%reset_flag = .false.

    write(LIS_logunit,fmt=*)"[INFO] Initializing CHIRPS v2.0 forcing parameters ... "

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

  ! Read CHIRPS 2.0 config file entries:
    call readconfig_chirps2()

    LIS_rc%met_nf(findex)   = 2          ! number of met variables in CHIRPS forcing
    LIS_rc%met_proj(findex) = "latlon"

 !- Select and define interp input arrays:
    if( chirps2_struc(1)%xres == 0.05 ) then
      glpnc = 7200
      glpnr = 2000
    elseif( chirps2_struc(1)%xres == 0.25 ) then
      glpnc = 1440
      glpnr = 400
    else
       write(LIS_logunit,*) "[ERR] Other resolutions are not available with CHIRPS 2.0" 
       write(LIS_logunit,*) "     precipitation reader, at this time.  Only"
       write(LIS_logunit,*) "     0.05 or 0.25 deg are available ..."
       write(LIS_logunit,*) "Program stopping ... "
       call LIS_endrun
    endif

  ! Grid description array:
    fulldomain_gridDesci(:) = 0
    fulldomain_gridDesci(1) = 0
    fulldomain_gridDesci(2) = glpnc
    fulldomain_gridDesci(3) = glpnr
    fulldomain_gridDesci(4) = -50.0 +  ( chirps2_struc(1)%yres/2 )
    fulldomain_gridDesci(5) = -180.0 + ( chirps2_struc(1)%xres/2 )
    fulldomain_gridDesci(6) = 128
    fulldomain_gridDesci(7) = 50.0 -  ( chirps2_struc(1)%yres/2 )
    fulldomain_gridDesci(8) = 180.0 + ( chirps2_struc(1)%xres/2 )
    fulldomain_gridDesci(9)  = chirps2_struc(1)%xres
    fulldomain_gridDesci(10) = chirps2_struc(1)%yres
    fulldomain_gridDesci(20) = 64

!-- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
    sub_gridDesci = 0.
    subpnc = 0; subpnr = 0
    chirps2_struc(:)%start_nc = 0
    chirps2_struc(:)%start_nr = 0

    do n=1,LIS_rc%nnest

       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then

         if(mod(LIS_rc%nensem(n),&
              LIS_forecast_struc(1)%niterations).ne.0) then
            write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
            write(LIS_logunit,*) '[ERR] of the number of iterations '
            write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
            write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
            call LIS_endrun()
         endif

         chirps2_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
         chirps2_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
         chirps2_struc(n)%nIter = LIS_forecast_struc(1)%niterations

         allocate(chirps2_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(chirps2_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))

       ! Regular retrospective or non-forecast mode:
       else

         chirps2_struc(n)%st_iterid = 1
         chirps2_struc(n)%en_iterId = 1
         chirps2_struc(n)%nIter = 1

         allocate(chirps2_struc(n)%metdata1(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(chirps2_struc(n)%metdata2(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))

       endif
       chirps2_struc(n)%metdata1 = 0
       chirps2_struc(n)%metdata2 = 0

#if 0
     ! Do not subset for speeding up runtime (read in full CHIRPS domain):
     ! if( optimize_runtime == 0 ) then
         chirps2_struc(n)%start_nc = 1
         chirps2_struc(n)%start_nr = 1
         subpnc = glpnc
         subpnr = glpnr
         sub_gridDesci(:) = fulldomain_gridDesci(:)
     ! elseif( optimize_runtime == 1 ) then
#endif 
     ! Subset option
       call LIS_RunDomainPts( n, LIS_rc%met_proj(findex), fulldomain_gridDesci(:), &
            glpnc, glpnr, subpnc, subpnr, sub_gridDesci, lat_line, lon_line )

       chirps2_struc(n)%start_nc = lon_line(1,1)
       chirps2_struc(n)%start_nr = lat_line(1,1)
   end do

! - Subsetted domains:
    chirps2_struc(:)%nc = subpnc
    chirps2_struc(:)%nr = subpnr

  ! First date available for the CHIRPS v2.0 dataset:
    yr1 = 1981
    mo1 = 01
    da1 = 01
    hr1 = 0
    mn1 = 0; ss1 = 0
    ts1 = 3600 * 24  ! Need to figure out when CHIRPS precip is valid for ...
    call LIS_tick(FirstDateTime,updoy,upgmt,yr1,mo1,&
         da1,hr1,mn1,ss1,real(ts1) )

  ! Starting date/time for user-specified inputs:
    LIS_syr = LIS_rc%syr
    LIS_smo = LIS_rc%smo
    LIS_sda = LIS_rc%sda
    LIS_shr = LIS_rc%shr
    LIS_smn = LIS_rc%smn
    LIS_sss = LIS_rc%sss
    ts2 = 0
    call LIS_tick(LIS_StartTime,updoy,upgmt,LIS_syr,LIS_smo,&
         LIS_sda,LIS_shr,LIS_smn,LIS_sss,real(ts2) )

  ! Starting date/time for user-specified inputs:
    if( FirstDateTime .gt. LIS_StartTime ) then
      write(LIS_logunit,*) "[WARN] LIS start time is earlier than CHIRPS 2.0 data availability time ..."
      write(LIS_logunit,*) " ... Relying on baseforcing precipitation to fill ..."
      if( LIS_rc%nmetforc == 1 ) then
        write(LIS_logunit,*) "[ERR] No other underlying forcings and no available CHIRPS precipitation files."
        write(LIS_logunit,*) " Program stopping. "
        call LIS_endrun
      endif
    endif

    do n=1,LIS_rc%nnest

       chirps2_struc(n)%gridDesc(:) = sub_gridDesci(:)
       chirps2_struc(n)%mi = chirps2_struc(n)%nc*chirps2_struc(n)%nr

       chirps2_struc(n)%ts = 24 * 3600.   ! seconds per day
       call LIS_update_timestep(LIS_rc, n, chirps2_struc(n)%ts)

       call LIS_registerAlarm("CHIRPS 2.0 alarm",&
            chirps2_struc(n)%ts,&
            chirps2_struc(n)%ts )

     ! Set local - 1 hour timestep (to replicate model timestep):
       call LIS_update_timestep(LIS_rc, n, 3600.) 

     ! Check if interp/upscale option is correct for resolution selected:
       ! LIS resolution < CHIRPS resolution
       if( LIS_isatAfinerResolution(n,chirps2_struc(n)%xres) ) then
         write(LIS_logunit,*) " The LIS Run domain is at a finer resolution than "
         write(LIS_logunit,*) "  CHIRPS 2.0 selected resolution, ",chirps2_struc(n)%xres

         if( trim(LIS_rc%met_interp(findex)) .eq. "average" ) then
            write(LIS_logunit,*) "[ERR] Given the selected LIS and CHIRPS resolutions, "
            write(LIS_logunit,*) "  user must select a downscaling option.  Current"
            write(LIS_logunit,*) "  options are: 'neighbor','budget-bilinear'.  Please select"
            write(LIS_logunit,*) "  one and run LIS again.  Program stopping ..."
            call LIS_endrun
         endif

       ! LIS resolution >= CHIRPS resolution
       else
         write(LIS_logunit,*) "[WARN] The LIS Run domain is at a coarser (or =) resolution than "
         write(LIS_logunit,*) "  CHIRPS 2.0 selected resolution, ",chirps2_struc(n)%xres

         LIS_rc%met_interp(findex) = LIS_rc%met_upscale(findex) 

         if( trim(LIS_rc%met_interp(findex)) .ne. "average" ) then
            write(LIS_logunit,*) "[ERR] Given the selected LIS and CHIRPS resolutions, "
            write(LIS_logunit,*) "  user must select an upscaling option.  Current"
            write(LIS_logunit,*) "  option available is:  'average'.  Please select"
            write(LIS_logunit,*) "  and run LIS again.  Program stopping ..."
            call LIS_endrun
         endif
       endif

    !- Select and define interp input arrays:
       select case( LIS_rc%met_interp(findex) )

        case( "budget-bilinear" )

          allocate(chirps2_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(chirps2_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(chirps2_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(chirps2_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,sub_gridDesci(:),&
               chirps2_struc(n)%n112,chirps2_struc(n)%n122,&
               chirps2_struc(n)%n212,chirps2_struc(n)%n222,&
               chirps2_struc(n)%w112,chirps2_struc(n)%w122,&
               chirps2_struc(n)%w212,chirps2_struc(n)%w222)

        case( "neighbor" )   

          allocate(chirps2_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n,sub_gridDesci(:),&
               chirps2_struc(n)%n113)

       ! Upscale option(s):
        case( "average" )

          allocate( chirps2_struc(n)%n111(chirps2_struc(n)%mi) )
          call upscaleByAveraging_input(sub_gridDesci,&
               LIS_rc%gridDesc(n,:),&
               chirps2_struc(n)%mi,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n),&
               chirps2_struc(n)%n111)

        case default
          write(LIS_logunit,*) "[ERR] This interpolation option not defined yet for CHIRPS data"
          write(LIS_logunit,*) "Program stopping ... "
          call LIS_endrun

       end select

    enddo

  end subroutine init_chirps2

end module chirps2_forcingMod
