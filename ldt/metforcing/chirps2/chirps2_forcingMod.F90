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
!  The implementation in LDT has the derived data type {\tt chirps2\_struc}
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
!    for each grid point in LDT, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Array containing the weights of the input grid
!    for each grid point in LDT, for neighbor search.
!  \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbor information for upscaling data
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!     computes the neighbor for neighbor interpolation
!  \end{description}
!
! !REVISION HISTORY:
!  KR Arsenault, 07/06/2015: Updated to support latest v2 files

! !USES:
  use ESMF
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_chirps2      ! defines the native resolution of
                              ! the input data
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
     character(len=LDT_CONST_PATH_LEN) :: directory  
     real*8             :: chirpstime1, chirpstime2
     logical            :: reset_flag 

     integer            :: start_nc, start_nr

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
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_chirps2(findex)
! !USES: 
    use LDT_coreMod, only : LDT_rc, LDT_config, LDT_domain, &
            LDT_isLDTatAfinerResolution
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep, &
                               LDT_tick, LDT_parseTimeString, &
                               LDT_registerAlarm 
    use LDT_logMod,  only : LDT_logunit, LDT_endrun, &
                            LDT_getNextUnitNumber, &
                            LDT_releaseUnitNumber, LDT_verify 
    use LDT_gridmappingMod

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
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
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

    integer :: LDT_syr,LDT_smo,LDT_sda,LDT_shr,LDT_smn,LDT_sss 
    real*8  :: FirstDateTime
    real*8  :: LDT_StartTime  
    character (len=10):: time 

    integer :: subpnc, subpnr, glpnc, glpnr
    real    :: sub_gridDesci(20)
    real    :: fulldomain_gridDesci(20)
    integer, allocatable  :: lat_line(:,:), lon_line(:,:)

! _________________________________________________________________

    allocate(chirps2_struc(LDT_rc%nnest))
    chirps2_struc%reset_flag = .false.

    write(LDT_logunit,fmt=*) "[INFO] Initializing CHIRPS v2.0 forcing parameters ... "

    ! Read CHIRPS 2.0 config file entries:
    call readconfig_chirps2()

    LDT_rc%met_nf(findex) = 2          ! number of met variables in CHIRPS forcing
    LDT_rc%met_ts(findex) = 24*3600    ! seconds per day
    LDT_rc%met_zterp(findex) = .false.
    LDT_rc%met_proj(findex)  = "latlon"

    ! Select and define interp input arrays:
    fulldomain_gridDesci(:) = 0
    if( chirps2_struc(1)%xres == 0.05 ) then
      glpnc = 7200  ! nc
      glpnr = 2000  ! nr
    elseif( chirps2_struc(1)%xres == 0.25 ) then
      glpnc = 1440  ! nc
      glpnr =  400  ! nr
    else
       write(LDT_logunit,*) "[ERR] Other resolutions are not available with CHIRPS 2.0" 
       write(LDT_logunit,*) "     precipitation reader, at this time.  Only"
       write(LDT_logunit,*) "     0.05 or 0.25 deg are available ..."
       write(LDT_logunit,*) "Program stopping ... "
       call LDT_endrun
    endif

    ! The Full Domain - Grid description array:
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

    do n=1,LDT_rc%nnest

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
       call LDT_RunDomainPts( n, LDT_rc%met_proj(findex), fulldomain_gridDesci(:), &
            glpnc, glpnr, subpnc, subpnr, sub_gridDesci, lat_line, lon_line )

       chirps2_struc(n)%start_nc = lon_line(1,1)
       chirps2_struc(n)%start_nr = lat_line(1,1)

!-- Code below is want to include more surrounding points for interp/averaging --
#if 0
       print *, " Before: "
       print *, subpnc, subpnr, lon_line(1,1), lat_line(1,1)
       print *, sub_gridDesci(1:10)
       if( lon_line(1,1) > 4 ) then
          chirps2_struc(n)%start_nc = lon_line(1,1)-2
          subpnc = subpnc + 4
          sub_gridDesci(2) = subpnc 
          sub_gridDesci(5) = sub_gridDesci(5) - 2*(chirps2_struc(n)%xres)
          sub_gridDesci(8) = sub_gridDesci(8) + 2*(chirps2_struc(n)%xres)
       endif
       if( lat_line(1,1) > 4 ) then
          chirps2_struc(n)%start_nr = lat_line(1,1)-2
          subpnr = subpnr + 4
          sub_gridDesci(3) = subpnr
          sub_gridDesci(4) = sub_gridDesci(4) - 2*(chirps2_struc(n)%yres)
          sub_gridDesci(7) = sub_gridDesci(7) + 2*(chirps2_struc(n)%yres)
       endif
       print *, " After: "
       print *, subpnc, subpnr, chirps2_struc(n)%start_nc, chirps2_struc(n)%start_nr
       print *, sub_gridDesci(1:10)
#endif
    end do

    ! Subsetted domains:
    chirps2_struc(:)%nc = subpnc
    chirps2_struc(:)%nr = subpnr
    LDT_rc%met_nc(findex) = subpnc
    LDT_rc%met_nr(findex) = subpnr
    LDT_rc%met_gridDesc(findex,:) = sub_gridDesci(:)

    ! Check for First date available for the CHIRPS v2.0 dataset:
    yr1 = 1981
    mo1 = 01
    da1 = 01
    hr1 = 0
    mn1 = 0; ss1 = 0
    ts1 = 3600 * 24  ! Need to figure out when CHIRPS precip is valid for ...
    call LDT_tick(FirstDateTime,updoy,upgmt,yr1,mo1,&
         da1,hr1,mn1,ss1,real(ts1) )

    ! Starting date/time for user-specified inputs:
    LDT_syr = LDT_rc%syr
    LDT_smo = LDT_rc%smo
    LDT_sda = LDT_rc%sda
    LDT_shr = LDT_rc%shr
    LDT_smn = LDT_rc%smn
    LDT_sss = LDT_rc%sss
    ts2 = 0
    call LDT_tick(LDT_StartTime,updoy,upgmt,LDT_syr,LDT_smo,&
         LDT_sda,LDT_shr,LDT_smn,LDT_sss,real(ts2) )

    ! Starting date/time for user-specified inputs:
    if( FirstDateTime .gt. LDT_StartTime ) then
       write(LDT_logunit,*) "[WARN] LDT start time is earlier than CHIRPS 2.0 data availability time ..."
       write(LDT_logunit,*) " ... Relying on baseforcing precipitation to fill ..."
      if( LDT_rc%nmetforc == 1 ) then
        write(LDT_logunit,*) "[ERR] No other underlying forcings and no available CHIRPS precipitation files."
        write(LDT_logunit,*) " Program stopping. "
        call LDT_endrun
      endif
    endif

    ! If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return
       
    ! Set up the main inputs:
    do  n=1,LDT_rc%nnest

       ! Time step settings:
       chirps2_struc(n)%ts = 24 * 3600.
       call LDT_update_timestep(LDT_rc, n, chirps2_struc(n)%ts)

       call LDT_registerAlarm("CHIRPS 2.0 alarm",&
            chirps2_struc(n)%ts,&
            chirps2_struc(n)%ts )

       ! Set local - 1 hour timestep (to replicate model timestep):
       call LDT_update_timestep(LDT_rc, n, 24*3600.) 

       ! Check and set interp/upscaling options:
       chirps2_struc(n)%mi = chirps2_struc(n)%nc*chirps2_struc(n)%nr

       ! Check if interp/upscale option is correct for resolution selected:
       ! LIS resolution < CHIRPS resolution
       if( LDT_isLDTatAfinerResolution(n,chirps2_struc(n)%xres) ) then
         write(LDT_logunit,*) "[WARN] The LIS Run domain is at a finer resolution than "
         write(LDT_logunit,*) " CHIRPS 2.0 selected resolution, ",chirps2_struc(n)%xres

         if( trim(LDT_rc%met_gridtransform(findex)) .eq. "average" ) then
            write(LDT_logunit,*) "[ERR] Given the selected LIS and CHIRPS resolutions, "
            write(LDT_logunit,*) " user must select a downscaling option.  Current"
            write(LDT_logunit,*) " options are: 'neighbor','budget-bilinear'.  Please select"
            write(LDT_logunit,*) " one and run LDT again.  Program stopping ..."
            call LDT_endrun
         endif
 
       ! LIS resolution >= CHIRPS resolution
       else
         write(LDT_logunit,*) "[WARN] The LIS Run domain is at a coarser (or =) resolution than "
         write(LDT_logunit,*) " CHIRPS 2.0 selected resolution, ",chirps2_struc(n)%xres

         if( LDT_rc%met_gridtransform(findex) .ne. "average" .and. &
             LDT_rc%met_gridtransform(findex) .ne. "neighbor" ) then
            write(LDT_logunit,*) "[ERR] Given the selected LIS and CHIRPS resolutions, "
            write(LDT_logunit,*) " user must select an upscaling option.  Current"
            write(LDT_logunit,*) " options available: 'average' | 'neighbor'."
            write(LDT_logunit,*) " Please select and run LDT again.  Program stopping ..."
            call LDT_endrun
         endif
       endif

       !- Select and define interp input arrays:
       select case( LDT_rc%met_gridtransform(findex) )

        ! Downscale (interp) option:
        case( "budget-bilinear" )

          allocate(chirps2_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(chirps2_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(chirps2_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(chirps2_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(chirps2_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,sub_gridDesci(:),&
               chirps2_struc(n)%n112,chirps2_struc(n)%n122,&
               chirps2_struc(n)%n212,chirps2_struc(n)%n222,&
               chirps2_struc(n)%w112,chirps2_struc(n)%w122,&
               chirps2_struc(n)%w212,chirps2_struc(n)%w222)

        ! Resample option:
        case( "neighbor" )   

          allocate(chirps2_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n,sub_gridDesci(:),&
               chirps2_struc(n)%n113)

        ! Upscale option(s):
        case( "average" )   

          allocate( chirps2_struc(n)%n111(chirps2_struc(n)%mi) )
          call upscaleByAveraging_input(sub_gridDesci,&
               LDT_rc%gridDesc(n,:),&
               chirps2_struc(n)%mi,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               chirps2_struc(n)%n111)

        case default
          write(LDT_logunit,*) "[ERR] This interpolation option not defined yet for CHIRPS data"
          write(LDT_logunit,*) " Program stopping ... "
          call LDT_endrun

       end select

    enddo

  end subroutine init_chirps2

end module chirps2_forcingMod
