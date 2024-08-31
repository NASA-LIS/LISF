!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"
module LVT_DataStreamsMod
!BOP
!
! !MODULE: LVT_DataStreamMod
! \label(LVT_DataStreamMod)
!
! !INTERFACE:
!
! !USES:
  use LVT_histDataMod
  use LVT_coreMod
  use LVT_logMod
  use LVT_LISoutputHandlerMod
  use LVT_navgemMod
  use LVT_timeMgrMod
  use map_utils
  use grib_api
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  The code in this file contains the basic datastructures and
!  control routines for handling the operations associated
!  with the datastreams. The invocation to read, perform
!  temporal averaging and resetting of the datastreams are
!  performed from this module. The calculations of derived
!  variables are also performed in this module.
!
! !FILES USED:
!
! !REVISION HISTORY:
!  02 Oct 2008    Sujay Kumar  Initial Specification
!
!EOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_DataStreamsInit
  public :: LVT_readDataStreams
  public :: LVT_writeDataStreams
  public :: LVT_tavgDataStreams
  public :: LVT_resetDataStreams
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

!EOP

contains

!BOP
!
! !ROUTINE: LVT_DataStreamsInit
! \label{LVT_DataStreamsInit}
!
! !INTERFACE:
  subroutine LVT_DataStreamsInit
!
! !USES:
    use LVT_datastream_pluginMod, only : LVT_datastream_plugin
    implicit none
!
! !DESCRIPTION:
!
!  This subroutine invokes the call to initialize each datastream.
!
!   The routines invoked are:
!   \begin{description}
!    \item[LVT\_datastream\_plugin] (\ref{LVT_datastream_plugin}) \newline
!      routine to register all the supported datastream plugin implementations
!    \item[LVT\_LISoutputInit] (\ref{LVT_LISoutputInit}) \newline
!      routine to initialize the handling of LIS output data (if LIS output
!      data is one of the datastreams)
!   \end{description}
!EOP
    integer                          :: kk
    type(LVT_metadataEntry), pointer :: ds1, ds2
    real                             :: gridDesci(50)

    call LVT_datastream_plugin

    call observationsetup(trim(LVT_rc%obssource(1))//char(0),1)
    call observationsetup(trim(LVT_rc%obssource(2))//char(0),2)

    if(LVT_rc%nDatastreams.gt.2) then
       call observationsetup(trim(LVT_rc%obssource(3))//char(0),3)
    endif

    call LVT_LISoutputInit

! checking for duplicate entries in a given datastream
! Note that this check is not enabled for three datastrems.
! The responsibility of ensuring non-duplicate entries is
! on the user.

    LVT_rc%ds1_dup = .false.
    ds1 => LVT_histData%head_ds1_list
    do while(associated(ds1))
       ds2 => ds1%next
       do while(associated(ds2))
          if(ds2%index.ne.ds1%index.and.&
               ds1%short_name.eq.ds2%short_name) then
             LVT_rc%ds1_dup = .true.
          endif
          ds2 => ds2%next
       enddo
       ds1 => ds1%next
    enddo

    LVT_rc%ds2_dup = .false.
    ds1 => LVT_histData%head_ds2_list
    do while(associated(ds1))
       ds2 => ds1%next
       do while(associated(ds2))
          if(ds2%index.ne.ds1%index.and.&
               ds1%short_name.eq.ds2%short_name) then
             LVT_rc%ds2_dup = .true.
          endif
          ds2 => ds2%next
       enddo
       ds1 => ds1%next
    enddo

!-------------------------------------------------------------------
! for 557 post, the HYCOM data is processed to include the water
! temperature fields
!-------------------------------------------------------------------
    if (LVT_rc%runmode .eq. "557 post") then


       ! EMK FIXME...Replace HYCOM with NAVGEM
       ! EMK 20220519...Reinstate HYCOM SST processing.
       if (LVT_rc%processHYCOM .eq. 1) then

          LVT_rc%HYCOM_proc_start = .true.

          ! First, handle water_temp
          LVT_rc%HYCOM_nc = 4500
          LVT_rc%HYCOM_nr = 2001

          gridDesci = 0
          gridDesci(1) = 0
          gridDesci(2) = LVT_rc%HYCOM_nc
          gridDesci(3) = LVT_rc%HYCOM_nr
          gridDesci(4) = -80.0
          gridDesci(5) = -180.0
          gridDesci(7) = 80.0
          gridDesci(8) = 180.0
          gridDesci(6) = 128
          gridDesci(9) = 0.08
          gridDesci(10) = 0.08
          gridDesci(20) = 64

          allocate(LVT_rc%HYCOM_n11(LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr))

          call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc,&
               LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%HYCOM_n11)

          LVT_histData%watertemp%short_name = "water_temp"
          LVT_histData%watertemp%long_name = "water_temp"
          LVT_histData%watertemp%standard_name = "water temperature"
          LVT_histData%watertemp%units = "K"
          LVT_histData%watertemp%nunits = 1
          LVT_histData%watertemp%format = 'F'
          LVT_histData%watertemp%vlevels = 1
          LVT_histData%watertemp%timeAvgOpt = 0
          LVT_histData%watertemp%startNlevs = 1
          LVT_histData%watertemp%endNlevs = 1
          allocate(LVT_histData%watertemp%value(LVT_rc%ngrid, &
               1, LVT_histData%watertemp%vlevels))
          allocate(LVT_histData%watertemp%unittypes(1))
          LVT_histData%watertemp%unittypes(1) = "K"

          ! Now handle Arctic sea ice fraction (aice)
          LVT_rc%HYCOM_aice_arc_nc = 4500
          LVT_rc%HYCOM_aice_arc_nr = 1251

          ! See LIS_PRIV_rcMod.F90 for documentation of gridDesc
          gridDesci = 0
          gridDesci(1) = 0 ! Lat/lon projection
          gridDesci(2) = LVT_rc%HYCOM_aice_arc_nc ! Number of columns
          gridDesci(3) = LVT_rc%HYCOM_aice_arc_nr ! Number of rows
          gridDesci(4) = 40.     ! Lower-left latitude (deg N)
          gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
          gridDesci(6) = 128     ! Not used
          gridDesci(7) = 90.0             ! Upper-right latitude (deg N)
          gridDesci(8) = 179.920043945312 ! Upper-right longitude (deg E)
          gridDesci(9) = 0.080017089844005795  ! delta-lon (deg)
          gridDesci(10) = 0.040000915527301117 ! delta-lat (deg)
          gridDesci(20) = 64  ! East-west ordering

          allocate(LVT_rc%HYCOM_aice_arc_n11( &
               LVT_rc%HYCOM_aice_arc_nc*LVT_rc%HYCOM_aice_arc_nr))

          call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc, &
               LVT_rc%HYCOM_aice_arc_nc*LVT_rc%HYCOM_aice_arc_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%HYCOM_aice_arc_n11)

          ! Now handle Antarctic sea ice fraction (aice)
          LVT_rc%HYCOM_aice_ant_nc = 4500
          LVT_rc%HYCOM_aice_ant_nr = 775

          ! See LIS_PRIV_rcMod.F90 for documentation of gridDesc
          gridDesci = 0
          gridDesci(1) = 0 ! Lat/lon projection
          gridDesci(2) = LVT_rc%HYCOM_aice_ant_nc ! Number of columns
          gridDesci(3) = LVT_rc%HYCOM_aice_ant_nr ! Number of rows
          gridDesci(4) = -80.4800033569336     ! Lower-left latitude (deg N)
          gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
          gridDesci(6) = 128     ! Not used
          gridDesci(7) = -49.5200004577637 ! Upper-right latitude (deg N)
          gridDesci(8) = 179.920043945312 ! Upper-right longitude (deg E)
          gridDesci(9) = 0.080017089844005795  ! delta-lon (deg)
          gridDesci(10) = 0.040000915527400593 ! delta-lat (deg)
          gridDesci(20) = 64  ! East-west ordering

          allocate(LVT_rc%HYCOM_aice_ant_n11( &
               LVT_rc%HYCOM_aice_ant_nc*LVT_rc%HYCOM_aice_ant_nr))

          call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc, &
               LVT_rc%HYCOM_aice_ant_nc*LVT_rc%HYCOM_aice_ant_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%HYCOM_aice_ant_n11)

          LVT_histData%aice%short_name = "aice"
          LVT_histData%aice%long_name = "aice"
          LVT_histData%aice%standard_name = "sea_ice_area_fraction"
          LVT_histData%aice%units = ""
          LVT_histData%aice%nunits = 1
          LVT_histData%aice%format = 'F'
          LVT_histData%aice%vlevels = 1
          LVT_histData%aice%timeAvgOpt = 0
          LVT_histData%aice%startNlevs = 1
          LVT_histData%aice%endNlevs = 1
          allocate(LVT_histData%aice%value(LVT_rc%ngrid, &
               1, LVT_histData%aice%vlevels))
          allocate(LVT_histData%aice%unittypes(1))
          LVT_histData%aice%unittypes(1) = ""

          ! Now handle Arctic sea ice thickness (hi)
          LVT_rc%HYCOM_hi_arc_nc = 4500
          LVT_rc%HYCOM_hi_arc_nr = 1251

          ! See LIS_PRIV_rcMod.F90 for documentation of gridDesc
          gridDesci = 0
          gridDesci(1) = 0 ! Lat/lon projection
          gridDesci(2) = LVT_rc%HYCOM_hi_arc_nc ! Number of columns
          gridDesci(3) = LVT_rc%HYCOM_hi_arc_nr ! Number of rows
          gridDesci(4) = 40.     ! Lower-left latitude (deg N)
          gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
          gridDesci(6) = 128     ! Not used
          gridDesci(7) = 90.0             ! Upper-right latitude (deg N)
          gridDesci(8) = 179.920043945312 ! Upper-right longitude (deg E)
          gridDesci(9) = 0.080017089844005795  ! delta-lon (deg)
          gridDesci(10) = 0.040000915527301117 ! delta-lat (deg)
          gridDesci(20) = 64  ! East-west ordering

          allocate(LVT_rc%HYCOM_hi_arc_n11( &
               LVT_rc%HYCOM_hi_arc_nc*LVT_rc%HYCOM_hi_arc_nr))

          call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc, &
               LVT_rc%HYCOM_hi_arc_nc*LVT_rc%HYCOM_hi_arc_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%HYCOM_hi_arc_n11)

          ! Now handle Antarctic sea ice thickness (hi)
          LVT_rc%HYCOM_hi_ant_nc = 4500
          LVT_rc%HYCOM_hi_ant_nr = 775

          ! See LIS_PRIV_rcMod.F90 for documentation of gridDesc
          gridDesci = 0
          gridDesci(1) = 0 ! Lat/lon projection
          gridDesci(2) = LVT_rc%HYCOM_hi_ant_nc ! Number of columns
          gridDesci(3) = LVT_rc%HYCOM_hi_ant_nr ! Number of rows
          gridDesci(4) = -80.4800033569336     ! Lower-left latitude (deg N)
          gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
          gridDesci(6) = 128     ! Not used
          gridDesci(7) = -49.5200004577637 ! Upper-right latitude (deg N)
          gridDesci(8) = 179.920043945312 ! Upper-right longitude (deg E)
          gridDesci(9) = 0.080017089844005795  ! delta-lon (deg)
          gridDesci(10) = 0.040000915527400593 ! delta-lat (deg)
          gridDesci(20) = 64  ! East-west ordering

          allocate(LVT_rc%HYCOM_hi_ant_n11( &
               LVT_rc%HYCOM_hi_ant_nc*LVT_rc%HYCOM_hi_ant_nr))

          call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc, &
               LVT_rc%HYCOM_hi_ant_nc*LVT_rc%HYCOM_hi_ant_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%HYCOM_hi_ant_n11)

          LVT_histData%hi%short_name = "hi"
          LVT_histData%hi%long_name = "hi"
          LVT_histData%hi%standard_name = "sea_ice_thickness"
          LVT_histData%hi%units = "m"
          LVT_histData%hi%nunits = 1
          LVT_histData%hi%format = 'F'
          LVT_histData%hi%vlevels = 1
          LVT_histData%hi%timeAvgOpt = 0
          LVT_histData%hi%startNlevs = 1
          LVT_histData%hi%endNlevs = 1
          allocate(LVT_histData%hi%value(LVT_rc%ngrid, &
               1, LVT_histData%hi%vlevels))
          allocate(LVT_histData%hi%unittypes(1))
          LVT_histData%hi%unittypes(1) = ""

       endif
    endif
  end subroutine LVT_DataStreamsInit

!BOP
!
! !ROUTINE: LVT_readDataStreams
! \label{LVT_readDataStreams}
!
! !INTERFACE:
  subroutine LVT_readDataStreams
!
! !USES:
    implicit none
!
!
! !DESCRIPTION:
!  This subroutine invokes the routines that read the datastreams
!
!EOP

    call readObservationSource(trim(LVT_rc%obssource(1))//char(0),1)
    call readObservationSource(trim(LVT_rc%obssource(2))//char(0),2)

    if(LVT_rc%nDatastreams.gt.2) then
       call readObservationSource(trim(LVT_rc%obssource(3))//char(0),3)
    endif
  end subroutine LVT_readDataStreams


!BOP
!
! !ROUTINE: LVT_writeDataStreams
! \label{LVT_writeDataStreams}
!
! !INTERFACE:
  subroutine LVT_writeDataStreams
!
! !USES:
    use LVT_logMod
    use LVT_coreMod, only: LVT_LIS_rc ! EMK
    use LVT_557post_ps41_snowMod ! EMK

    implicit none
!
!
! !DESCRIPTION:
!  This subroutine invokes the routines that writes the datastream values
!  to an external file. Currently this feature is only supported
!  in the '557 post' runmode, for the processing of LIS outputs to the
!  grib format. The datastream1 output must be set to 'LIS output'.
!
!EOP


    integer, parameter                   :: nsoillayers = 4
    character*200                        :: fname_mean,fname_ssdev
    character(len=8)                     :: cdate2
    character(len=4)                     :: cdate3
    integer                              :: ftn_mean,ftn_ssdev
    integer                              :: time_unit
    integer                              :: time_past
    integer                              :: time_curr
    integer                              :: iret
    integer                              :: timeRange
    integer                              :: yr,mo,da,hr,mn,ss
    character*10                         :: stepType
    type(LVT_metadataEntry), pointer     :: dataEntry
    type(LVT_LISmetadataEntry), pointer  :: lisdataEntry
    real, dimension(nsoillayers)         :: toplev, botlev, depscale
    real, dimension(1)                   :: toplev0, botlev0
    real                                 :: lyrthk(nsoillayers)
    integer                              :: i,k,t,m
    real*8                               :: time
    real                                 :: gmt
    integer                              :: doy
    integer                              :: pdTemplate
    integer                              :: c,r,gid
    real                                 :: gtmp1_1d_mem(LVT_rc%lnc*LVT_rc%lnr)
    real                                 :: gtmp1_1d(LVT_rc%lnc*LVT_rc%lnr)
    real                                 :: gtmp1_ss(LVT_rc%lnc*LVT_rc%lnr)
    integer                              :: ngtmp1_1d(LVT_rc%lnc*LVT_rc%lnr)
    real                                 :: variance
    real                                 :: lat(LVT_rc%lnc,LVT_rc%lnr)
    real                                 :: lon(LVT_rc%lnc,LVT_rc%lnr)
    character*20                         :: output_fmt
    integer                              :: shuffle, deflate, deflate_level
    integer                              :: xlatID,xlonID,xtimeID
    integer                              :: xlat_ss_ID,xlon_ss_ID,xtime_ss_ID
    integer                              :: dimID(4),tdimID
    character*8                          :: xtime_begin_date
    character*6                          :: xtime_begin_time
    character*50                         :: xtime_units
    character*50                         :: xtime_timeInc
    character(len=8)                     :: date
    character(len=10)                    :: stime
    character(len=5)                     :: zone
    integer, dimension(8)                :: values
    type(LVT_metadataEntry)              :: xlat,xlon
    character*10                         :: stepType_mean
    character*10                         :: stepType_ssdev

    ! EMK...For Welford algorithm
    integer :: count
    real :: mean, m2, stddev, new_value

    ! EMK...Special processing of some JULES PS41 multi-layer snow ensembles.
    integer :: count_jules_ps41_ens_snow_vars
    logical :: jules_ps41_ens_snow
    logical :: is_ps41_snow_var
    character(len=100) :: short_name
    type(LVT_lismetadataEntry), target :: GrndSnow
    type(LVT_lismetadataEntry), target :: LayerSnowDensity
    type(LVT_lismetadataEntry), target :: SWE
    type(LVT_lismetadataEntry), target :: SnowDepth
    type(LVT_lismetadataEntry), target :: SnowDensity
    type(LVT_lismetadataEntry), target :: SnowGrain
    type(LVT_lismetadataEntry), target :: SurftSnow

    character*20 :: model_name ! EMK

    ! EMK...This is only used when LVT is run in "557 post" mode.
    if (trim(LVT_rc%runmode) .ne. "557 post") return

!    output_fmt = "grib2"

!    lyrthk(1) = 0.1*100.0
!    lyrthk(2) = 0.3*100.0
!    lyrthk(3) = 0.6*100.0
!    lyrthk(4) = 1.0*100.0
    ! EMK...Use soil thicknesses read in from file.
    if (LVT_LIS_rc(1)%nsmlayers .ne. nsoillayers) then
       write(LVT_logunit,*) '[ERR] Internal error, bad value of soil layers!'
       write(LVT_logunit,*) 'Program failed in LVT_writeDataStreams'
       call LVT_endrun()
    end if
    if (LVT_LIS_rc(1)%nstlayers .ne. nsoillayers) then
       write(LVT_logunit,*) '[ERR] Internal error, bad value of soil layers!'
       write(LVT_logunit,*) 'Program failed in LVT_writeDataStreams'
       call LVT_endrun()
    end if

    ! EMK...Soil layers are in centimeters.  Technically GRIB2 requires
    ! meters, but in practice we keep as centimeters and just modify
    ! the scale factor by 100.
    if (LVT_rc%lvt_out_format .eq. "grib2") then
       depscale = 2
    else
       depscale = 0
    end if
    do k = 1, nsoillayers
       lyrthk(k) = LVT_LIS_rc(1)%smthick(k)
    end do

    if (LVT_557post_alarm_is_on()) then
       ! EMK...We need lat/lon for all grid points, not just for land.
       ! So we will recalculate here.
       ! FIXME...Add support for other projections, not just lat/lon.
       lat = LVT_rc%udef
       lon = LVT_rc%udef
       do r = 1, LVT_rc%gnr
          do c = 1, LVT_rc%gnc
             lat(c,r) = LVT_rc%gridDesc(4) + (r-1)*LVT_rc%gridDesc(10)
             lon(c,r) = LVT_rc%gridDesc(5) + (c-1)*LVT_rc%gridDesc(9)
          enddo
       enddo

       ! EMK...Special handling of JULES PS41 multi-layer snow physics
       ! when ensembles are processed.
       ! FIXME...Add LVT flag specifying PS41?
       dataEntry => LVT_histData%head_ds1_list
       jules_ps41_ens_snow = .false.
       if (trim(LVT_LIS_rc(1)%anlys_data_class) .eq. "LSM" .and. &
            trim(LVT_LIS_rc(1)%model_name) .eq. "JULES.5.0" .and. &
            LVT_rc%nensem .gt. 1) then

          ! It's possible JULES PS41 snow variables will be processed.
          ! Take some preliminary steps.
          call LVT_init_jules_PS41_ens_snow()

          count_jules_ps41_ens_snow_vars = 0
          do while(associated(dataEntry))

             if (dataEntry%timeAvgOpt .ne. 0) then
                dataEntry => dataEntry%next
                cycle
             end if

             short_name = trim(dataEntry%short_name)//"_inst"

             call LVT_prep_jules_ps41_ens_snow_var(short_name, &
                  dataEntry%vlevels, dataEntry%value, is_ps41_snow_var)

             if (is_ps41_snow_var) then
                count_jules_ps41_ens_snow_vars = &
                     count_jules_ps41_ens_snow_vars + 1
             end if

             if (count_jules_ps41_ens_snow_vars .eq. 6) exit
             dataEntry => dataEntry%next
          end do

          if (count_jules_ps41_ens_snow_vars .eq. 0) then
             jules_ps41_ens_snow = .false.
          else if (count_jules_ps41_ens_snow_vars .ne. 6) then
             write(LVT_logunit,*) &
                  '[ERR] Cannot process JULES PS41 multi-layer snow'
             write(LVT_logunit,*) &
                  '[ERR] Missing some snow variables for ensemble processing'
             flush(LVT_logunit)
             call LVT_endrun()
          else
             jules_ps41_ens_snow = .true.
          end if

          ! Go to head of list for later processing
          dataEntry => LVT_histData%head_ds1_list
       end if

       if (jules_ps41_ens_snow) then
          ! We now have all the PS41 snow variables needed for ensemble
          ! processing.  We invoke the relayer algorithm.  The output variables
          ! will be fetched further down.
          call LVT_proc_jules_ps41_ens_snow()

          ! Set metadata for new derived snow variables not copied from
          ! LIS file.
          call LVT_set_SWE_metadata(SWE)
          call LVT_set_SnowDensity_metadata(SnowDensity)
          call LVT_set_LayerSnowDensity_metadata(LayerSnowDensity)
          call LVT_set_SnowGrain_metadata(SnowGrain)
          call LVT_set_SnowDepth_metadata(SnowDepth)
          call LVT_set_GrndSnow_metadata(GrndSnow)
          call LVT_set_SurftSnow_metadata(SurftSnow)
       end if
       ! EMK END JULES PS41 Snow

       if (LVT_rc%lvt_out_format .eq. "grib1") then

          write(unit=cdate2, fmt='(i4.4,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(unit=cdate3, fmt='(i2.2,i2.2)') &
               LVT_rc%hr, LVT_rc%mn

          ! EMK...Include LSM in GP section
          if (trim(LVT_LIS_rc(1)%model_name) == "NOAH.3.9") then
             model_name = "LIS-NOAH"
          else if (trim(LVT_LIS_rc(1)%model_name) == "NOAHMP.4.0.1") then
             model_name = "LIS-NOAHMP"
          else if (trim(LVT_LIS_rc(1)%model_name) == "JULES.5.0") then
             model_name = "LIS-JULES"
          else
             write(LVT_logunit,*)'[ERR] Unknown LSM selected'
             write(LVT_logunit,*)&
                  '[ERR] Must be NOAH.3.9, NOAHMP.4.0.1, or JULES.5.0'
             write(LVT_logunit,*) &
                  "[ERR] Update 'LIS output model name:' in lis.config" // &
                  " and try again!"
             call LVT_endrun()
          end if

          ! EMK...Different file name convention for 24-hr data
          if (LVT_rc%tavgInterval == 86400) then

             fname_mean = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.LIS24' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.GR1'

             fname_ssdev = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.LIS24-SSDEV' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.GR1'
          else

             fname_mean = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.LIS' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.GR1'

             fname_ssdev = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.SSDEV' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.GR1'

          end if

          ! Setup of GRIB-1 and GRIB-2 Metadata Section

          ! toplev is the depth of the top of each soil layer
          ! botlev is the depth of the bottom of each soil layer
          toplev(1) = 0.0
          botlev(1) = lyrthk(1)

          ! determine bounding levels for each soil moisture layer
          do i = 2, nsoillayers
             toplev(i) = toplev(i-1) + lyrthk(i-1)
             botlev(i) = botlev(i-1) + lyrthk(i)
          enddo
          !hardcoded to zero for now
          !depscale = 0

          ! Set values for non layered fields (Fluxes, Sfc Fields, etc.)
          toplev0 = 0
          botlev0 = 0

          yr = LVT_rc%yr
          mo = LVT_rc%mo
          da = LVT_rc%da
          hr = LVT_rc%hr
          mn = LVT_rc%mn
          ss = LVT_rc%ss

          call LVT_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,-1*LVT_rc%statswriteint)

          if(LVT_rc%statswriteint .GT. 0) then
             time_unit = 254     ! seconds
             time_curr = 0
             time_past = LVT_rc%statswriteint
          endif
          if(LVT_rc%statswriteint .GE. 60) then
             time_unit = 0      ! minutes
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 60)
          endif
          if(LVT_rc%statswriteint .GE. 3600) then
             time_unit = 1    ! hours
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 3600)
          endif
          if(LVT_rc%statswriteint .GE. 86400) then
             time_unit = 2   ! days
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 86400)
          endif

          !time_past: from LVT_grib1_finalize
          !time_P1 (Negative Time Unit for avg, or 0 for analysis)
          !According to the in-line comments, time_past must be negative or 0.
          !Here we are setting it to a positive value.  This produces bad output.
          !Setting it to a negative value also produces bad output.
          !So I am resetting it to zero.  This produces output that matches
          !the binary output.
          !    time_past=0

          call grib_open_file(ftn_mean, fname_mean, 'w', iret)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call grib_open_file(ftn_ssdev, fname_ssdev, 'w', iret)
             call LVT_verify(iret, &
                  'failed to open grib file '//trim(fname_ssdev))
          end if

       elseif (LVT_rc%lvt_out_format .eq. "grib2") then
          write(unit=cdate2, fmt='(i4.4,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(unit=cdate3, fmt='(i2.2,i2.2)') &
               LVT_rc%hr, LVT_rc%mn

          ! EMK...Include LSM in GP section
          if (trim(LVT_LIS_rc(1)%model_name) == "NOAH.3.9") then
             model_name = "LIS-NOAH"
          else if (trim(LVT_LIS_rc(1)%model_name) == "NOAHMP.4.0.1") then
             model_name = "LIS-NOAHMP"
          else if (trim(LVT_LIS_rc(1)%model_name) == "JULES.5.0") then
             model_name = "LIS-JULES"
          else
             write(LVT_logunit,*)'[ERR] Unknown LSM selected'
             write(LVT_logunit,*)&
                  '[ERR] Must be NOAH.3.9, NOAHMP.4.0.1, or JULES.5.0'
             write(LVT_logunit,*) &
                  "[ERR] Update 'LIS output model name:' in lis.config" // &
                  " and try again!"
             call LVT_endrun()
          end if

          ! EMK...Different file name convention for 24-hr data
          if (LVT_rc%tavgInterval == 86400) then
             fname_mean = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.LIS24' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.GR2'

             if (LVT_rc%nensem > 1) then
                fname_ssdev = trim(LVT_rc%statsodir) &
                     //'/PS.557WW' &
                     //'_SC.'//trim(LVT_rc%security_class) &
                     //'_DI.'//trim(LVT_rc%data_category) &
                     //'_GP.'//trim(model_name) &
                     //'_GR.C0P09DEG' &
                     //'_AR.'//trim(LVT_rc%area_of_data) &
                     //'_PA.LIS24-SSDEV' &
                     //'_DD.'//trim(cdate2) &
                     //'_DT.'//trim(cdate3) &
                     //'_DF.GR2'
             end if
          else
             ! EMK...Assume 3-hr
             fname_mean = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.LIS' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.GR2'

             if (LVT_rc%nensem > 1) then
                fname_ssdev = trim(LVT_rc%statsodir) &
                     //'/PS.557WW' &
                     //'_SC.'//trim(LVT_rc%security_class) &
                     //'_DI.'//trim(LVT_rc%data_category) &
                     //'_GP.'//trim(model_name) &
                     //'_GR.C0P09DEG' &
                     //'_AR.'//trim(LVT_rc%area_of_data) &
                     //'_PA.SSDEV' &
                     //'_DD.'//trim(cdate2) &
                     //'_DT.'//trim(cdate3) &
                     //'_DF.GR2'
             end if
          end if
          ! Setup of GRIB-1 and GRIB-2 Metadata Section

          ! toplev is the depth of the top of each soil layer
          ! botlev is the depth of the bottom of each soil layer
          toplev(1) = 0.0
          botlev(1) = lyrthk(1)

          ! determine bounding levels for each soil moisture layer
          do i = 2, nsoillayers
             toplev(i) = toplev(i-1) + lyrthk(i-1)
             botlev(i) = botlev(i-1) + lyrthk(i)
          enddo
          !hardcoded to zero for now
          !depscale = 0

          ! Set values for non layered fields (Fluxes, Sfc Fields, etc.)
          toplev0 = 0
          botlev0 = 0

          yr = LVT_rc%yr
          mo = LVT_rc%mo
          da = LVT_rc%da
          hr = LVT_rc%hr
          mn = LVT_rc%mn
          ss = LVT_rc%ss

          call LVT_tick(time, doy, gmt, yr, mo, da, hr, mn, ss, &
               -1*LVT_rc%statswriteint)

          if (LVT_rc%statswriteint .GT. 0) then
             time_unit = 254     ! seconds
             time_curr = 0
             time_past = LVT_rc%statswriteint
          endif
          if (LVT_rc%statswriteint .GE. 60) then
             time_unit = 0      ! minutes
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 60)
          endif
          if (LVT_rc%statswriteint .GE. 3600) then
             time_unit = 1    ! hours
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 3600)
          endif
          if (LVT_rc%statswriteint .GE. 86400) then
             time_unit = 2   ! days
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 86400)
          endif

          !time_past: from LVT_grib1_finalize
          !time_P1 (Negative Time Unit for avg, or 0 for analysis)
          !According to the in-line comments, time_past must be negative or 0.
          !Here we are setting it to a positive value.  This produces bad
          !output. Setting it to a negative value also produces bad output.
          !So I am resetting it to zero.  This produces output that matches
          !the binary output.
          !    time_past=0

          call grib_open_file(ftn_mean,fname_mean, 'w', iret)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call grib_open_file(ftn_ssdev, fname_ssdev, 'w', iret)
             call LVT_verify(iret, &
                  'failed to open grib file '//trim(fname_ssdev))
          end if

       elseif (LVT_rc%lvt_out_format .eq. "netcdf") then

          call date_and_time(date, stime, zone, values)

          write(unit=cdate2, fmt='(i4.4,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(unit=cdate3, fmt='(i2.2,i2.2)') &
               LVT_rc%hr, LVT_rc%mn

          ! EMK...Include LSM in GP section
          if (trim(LVT_LIS_rc(1)%model_name) == "NOAH.3.9") then
             model_name = "LIS-NOAH"
          else if (trim(LVT_LIS_rc(1)%model_name) == "NOAHMP.4.0.1") then
             model_name = "LIS-NOAHMP"
          else if (trim(LVT_LIS_rc(1)%model_name) == "JULES.5.0") then
             model_name = "LIS-JULES"
          else
             write(LVT_logunit,*)'[ERR] Unknown LSM selected'
             write(LVT_logunit,*)&
                  '[ERR] Must be NOAH.3.9, NOAHMP.4.0.1, or JULES.5.0'
             write(LVT_logunit,*) &
                  "[ERR] Update 'LIS output model name:' in lis.config" // &
                  " and try again!"
             call LVT_endrun()
          end if

          ! EMK...Different file name convention for 24-hr data
          if (LVT_rc%tavgInterval == 86400) then
             fname_mean = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.LIS24' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.nc'

             if (LVT_rc%nensem > 1) then
                fname_ssdev = trim(LVT_rc%statsodir) &
                     //'/PS.557WW' &
                     //'_SC.'//trim(LVT_rc%security_class) &
                     //'_DI.'//trim(LVT_rc%data_category) &
                     //'_GP.'//trim(model_name) &
                     //'_GR.C0P09DEG' &
                     //'_AR.'//trim(LVT_rc%area_of_data) &
                     //'_PA.LIS24-SSDEV' &
                     //'_DD.'//trim(cdate2) &
                     //'_DT.'//trim(cdate3) &
                     //'_DF.nc'
             end if
          else

             fname_mean = trim(LVT_rc%statsodir) &
                  //'/PS.557WW' &
                  //'_SC.'//trim(LVT_rc%security_class) &
                  //'_DI.'//trim(LVT_rc%data_category) &
                  //'_GP.'//trim(model_name) &
                  //'_GR.C0P09DEG' &
                  //'_AR.'//trim(LVT_rc%area_of_data) &
                  //'_PA.LIS' &
                  //'_DD.'//trim(cdate2) &
                  //'_DT.'//trim(cdate3) &
                  //'_DF.nc'

             if (LVT_rc%nensem > 1) then
                fname_ssdev = trim(LVT_rc%statsodir) &
                     //'/PS.557WW' &
                     //'_SC.'//trim(LVT_rc%security_class) &
                     //'_DI.'//trim(LVT_rc%data_category) &
                     //'_GP.'//trim(model_name) &
                     //'_GR.C0P09DEG' &
                     //'_AR.'//trim(LVT_rc%area_of_data) &
                     //'_PA.SSDEV' &
                     //'_DD.'//trim(cdate2) &
                     //'_DT.'//trim(cdate3) &
                     //'_DF.nc'
             end if
          end if
          ! Setup of GRIB-1 and GRIB-2 Metadata Section

          ! toplev is the depth of the top of each soil layer
          ! botlev is the depth of the bottom of each soil layer
          toplev(1) = 0.0
          botlev(1) = lyrthk(1)

          ! determine bounding levels for each soil moisture layer
          do i = 2, nsoillayers
             toplev(i) = toplev(i-1) + lyrthk(i-1)
             botlev(i) = botlev(i-1) + lyrthk(i)
          enddo
          !hardcoded to zero for now
          !depscale = 0

          ! Set values for non layered fields (Fluxes, Sfc Fields, etc.)
          toplev0 = 0
          botlev0 = 0

          yr = LVT_rc%yr
          mo = LVT_rc%mo
          da = LVT_rc%da
          hr = LVT_rc%hr
          mn = LVT_rc%mn
          ss = LVT_rc%ss

          call LVT_tick(time, doy, gmt, yr, mo, da, hr, mn, ss, &
               -1*LVT_rc%statswriteint)

          if (LVT_rc%statswriteint .GT. 0) then
             time_unit = 254     ! seconds
             time_curr = 0
             time_past = LVT_rc%statswriteint
          endif
          if (LVT_rc%statswriteint .GE. 60) then
             time_unit = 0      ! minutes
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 60)
          endif
          if (LVT_rc%statswriteint .GE. 3600) then
             time_unit = 1    ! hours
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 3600)
          endif
          if (LVT_rc%statswriteint .GE. 86400) then
             time_unit = 2   ! days
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 86400)
          endif

          !time_past: from LVT_grib1_finalize
          !time_P1 (Negative Time Unit for avg, or 0 for analysis)
          !According to the in-line comments, time_past must be negative or 0.
          !Here we are setting it to a positive value.  This produces bad
          !output. Setting it to a negative value also produces bad output.
          !So I am resetting it to zero.  This produces output that matches
          !the binary output.
          !    time_past=0

          shuffle = NETCDF_shuffle
          deflate = NETCDF_deflate
          deflate_level =NETCDF_deflate_level

          xlat%short_name = "latitude"
          xlat%long_name = "latitude"
          xlat%standard_name = "latitude"
          xlat%units = "degree_north"
          xlat%nunits = 1
          xlat%format = 'F'
          xlat%vlevels = 1
          xlat%timeAvgOpt = 0
          xlat%startNlevs = 1
          xlat%endNlevs = 1
          allocate(xlat%value(LVT_rc%ngrid, &
               1, xlat%vlevels))
          allocate(xlat%unittypes(1))
          xlat%unittypes(1) = "degree_north"

          xlon%short_name = "longitude"
          xlon%long_name = "longitude"
          xlon%standard_name = "longitude"
          xlon%units = "degree_east"
          xlon%nunits = 1
          xlon%format = 'F'
          xlon%vlevels = 1
          xlon%timeAvgOpt = 0
          xlon%startNlevs = 1
          xlon%endNlevs = 1
          allocate(xlon%value(LVT_rc%ngrid, &
               1, xlon%vlevels))
          allocate(xlon%unittypes(1))
          xlon%unittypes(1) = "degree_east"

#if (defined USE_NETCDF4)
          iret = nf90_create(path=trim(fname_mean), cmode=nf90_hdf5, &
               ncid=ftn_mean)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             iret = nf90_create(path=trim(fname_ssdev), cmode=nf90_hdf5, &
                  ncid=ftn_ssdev)
             call LVT_verify(iret, &
                  'failed to open grib file '//trim(fname_ssdev))
          end if
#endif
#if (defined USE_NETCDF3)
          iret = nf90_create(path=trim(fname_mean), cmode=nf90_clobber, &
               ncid=ftn_mean)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             iret = nf90_create(path=trim(fname_ssdev), cmode=nf90_clobber, &
                  ncid=ftn_ssdev)
             call LVT_verify(iret, &
                  'failed to open grib file '//trim(fname_ssdev))
          end if
#endif
          !Headers
          call LVT_verify(nf90_def_dim(ftn_mean, 'east_west', &
               LVT_rc%gnc, dimID(1)))
          call LVT_verify(nf90_def_dim(ftn_mean, 'north_south', &
               LVT_rc%gnr, dimID(2)))

          call LVT_verify(nf90_def_dim(ftn_mean, 'time', 1, tdimID))
          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
               "missing_value", &
               LVT_rc%udef))

          call LVT_verify(nf90_def_var(ftn_mean, &
               trim(xlat%short_name), &
               nf90_float, &
               dimids = dimID(1:2), varID=xlatID))
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_deflate(ftn_mean, &
               xlatID, &
               shuffle, deflate, deflate_level))
#endif
          call LVT_verify(nf90_def_var(ftn_mean, &
               trim(xlon%short_name), &
               nf90_float, &
               dimids = dimID(1:2), varID=xlonID))
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_deflate(ftn_mean, &
               xlonID, &
               shuffle, deflate, deflate_level))
#endif
          call LVT_verify(nf90_put_att(ftn_mean, xlatID, &
               "units", trim(xlat%units)))
          call LVT_verify(nf90_put_att(ftn_mean, xlatID, &
               "standard_name", trim(xlat%standard_name)))
          call LVT_verify(nf90_put_att(ftn_mean, xlatID, &
               "long_name", trim(xlat%long_name)))
          call LVT_verify(nf90_put_att(ftn_mean, xlatID, &
               "scale_factor", 1.0))
          call LVT_verify(nf90_put_att(ftn_mean, xlatID, &
               "add_offset", 0.0))
          call LVT_verify(nf90_put_att(ftn_mean, xlatID, &
               "missing_value", LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn_mean, xlatID, &
               "_FillValue", LVT_rc%udef))

          call LVT_verify(nf90_put_att(ftn_mean, xlonID, &
               "units", trim(xlon%units)))
          call LVT_verify(nf90_put_att(ftn_mean, xlonID, &
               "standard_name", trim(xlon%standard_name)))
          call LVT_verify(nf90_put_att(ftn_mean, xlonID, &
               "long_name", trim(xlon%long_name)))
          call LVT_verify(nf90_put_att(ftn_mean, xlonID, &
               "scale_factor", 1.0))
          call LVT_verify(nf90_put_att(ftn_mean, xlonID, &
               "add_offset", 0.0))
          call LVT_verify(nf90_put_att(ftn_mean, xlonID, &
               "missing_value", LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn_mean, xlonID, &
               "_FillValue", LVT_rc%udef))

!define time field
          call LVT_verify(nf90_def_var(ftn_mean, 'time', &
               nf90_float, dimids=tdimID, varID=xtimeID))
          write(xtime_units,200) LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
200       format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
               I2.2,':',I2.2)
          write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
          write(xtime_timeInc, fmt='(I20)') &
               LVT_rc%ts

          call LVT_verify(nf90_put_att(ftn_mean, xtimeID,&
               "units", trim(xtime_units)))
          call LVT_verify(nf90_put_att(ftn_mean, xtimeID,&
               "long_name", "time"))
          call LVT_verify(nf90_put_att(ftn_mean, xtimeID,&
               "time_increment", trim(adjustl(xtime_timeInc))))
          call LVT_verify(nf90_put_att(ftn_mean, xtimeID,&
               "begin_date", xtime_begin_date))
          call LVT_verify(nf90_put_att(ftn_mean, xtimeID,&
               "begin_time", xtime_begin_time))

          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "title", &
               "LVT land surface analysis output"))
          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "institution", &
               trim(LVT_rc%institution)))
          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "history", &
               "created on date: "//date(1:4)//"-"//date(5:6)//"-"// &
               date(7:8)//"T"//stime(1:2)//":"//stime(3:4)//":"//stime(5:10)))
          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "references", &
               "Kumar_etal_GMD_2012"))
          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "comment", &
               "website: http://lis.gsfc.nasa.gov/"))

          !grid information
          if (trim(LVT_rc%domain) .eq. "latlon") then !latlon
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL"))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LAT", LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LON", LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DX", &
                  LVT_rc%gridDesc(9)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DY", &
                  LVT_rc%gridDesc(10)))
          elseif (trim(LVT_rc%domain) .eq. "mercator") then
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "MAP_PROJECTION", "MERCATOR"))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LAT", LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LON", LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "STANDARD_LON", LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DY", &
                  LVT_rc%gridDesc(9)))
          elseif (trim(LVT_rc%domain).eq."lambert") then !lambert conformal
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "MAP_PROJECTION", "LAMBERT CONFORMAL"))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LAT", LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LON", LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "TRUELAT1", LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "TRUELAT2", LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "STANDARD_LON", LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DY", &
                  LVT_rc%gridDesc(9)))

          elseif (trim(LVT_rc%domain) .eq. "polar") then ! polar stereographic
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "MAP_PROJECTION", "POLAR STEREOGRAPHIC"))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LAT", LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "SOUTH_WEST_CORNER_LON", LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "ORIENT", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
                  "STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, "DY", &
                  LVT_rc%gridDesc(9)))
          endif

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             !Headers
             call LVT_verify(nf90_def_dim(ftn_ssdev, 'east_west', &
                  LVT_rc%gnc, dimID(1)))
             call LVT_verify(nf90_def_dim(ftn_ssdev, 'north_south', &
                  LVT_rc%gnr,dimID(2)))
             call LVT_verify(nf90_def_dim(ftn_ssdev,'time', 1, tdimID))
             call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                  "missing_value", &
                  LVT_rc%udef))
             call LVT_verify(nf90_def_var(ftn_ssdev, &
                  trim(xlat%short_name), &
                  nf90_float, &
                  dimids=dimID(1:2), varID=xlat_ss_ID))
#if(defined USE_NETCDF4)
             call LVT_verify(nf90_def_var_deflate(ftn_ssdev, &
                  xlat_ss_ID, &
                  shuffle, deflate, deflate_level))
#endif
             call LVT_verify(nf90_def_var(ftn_ssdev, &
                  trim(xlon%short_name), &
                  nf90_float, &
                  dimids=dimID(1:2), varID=xlon_ss_ID))
#if(defined USE_NETCDF4)
             call LVT_verify(nf90_def_var_deflate(ftn_ssdev, &
                  xlon_ss_ID, &
                  shuffle, deflate, deflate_level))
#endif
             call LVT_verify(nf90_put_att(ftn_ssdev, xlat_ss_ID, &
                  "units", trim(xlat%units)))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlat_ss_ID, &
                  "standard_name", trim(xlat%standard_name)))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlat_ss_ID, &
                  "long_name", trim(xlat%long_name)))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlat_ss_ID, &
                  "scale_factor", 1.0))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlat_ss_ID, &
                  "add_offset", 0.0))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlat_ss_ID, &
                  "missing_value", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlat_ss_ID, &
                  "_FillValue", LVT_rc%udef))

             call LVT_verify(nf90_put_att(ftn_ssdev, xlon_ss_ID, &
                  "units", trim(xlon%units)))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlon_ss_ID, &
                  "standard_name", trim(xlon%standard_name)))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlon_ss_ID, &
                  "long_name", trim(xlon%long_name)))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlon_ss_ID, &
                  "scale_factor", 1.0))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlon_ss_ID, &
                  "add_offset", 0.0))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlon_ss_ID, &
                  "missing_value", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_ssdev, xlon_ss_ID, &
                  "_FillValue", LVT_rc%udef))

             !define time field
             call LVT_verify(nf90_def_var(ftn_ssdev,'time', &
                  nf90_float, dimids=tdimID, varID=xtime_ss_ID))
             write(xtime_units,201) LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
                  LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
201          format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
                  I2.2,':',I2.2)
             write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
                  LVT_rc%yr, LVT_rc%mo, LVT_rc%da
             write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
                  LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
             write(xtime_timeInc, fmt='(I20)') &
                  LVT_rc%ts

             call LVT_verify(nf90_put_att(ftn_ssdev, xtime_ss_ID,&
                  "units", trim(xtime_units)))
             call LVT_verify(nf90_put_att(ftn_ssdev, xtime_ss_ID,&
                  "long_name", "time"))
             call LVT_verify(nf90_put_att(ftn_ssdev, xtime_ss_ID,&
                  "time_increment", trim(adjustl(xtime_timeInc))))
             call LVT_verify(nf90_put_att(ftn_ssdev, xtime_ss_ID,&
                  "begin_date", xtime_begin_date))
             call LVT_verify(nf90_put_att(ftn_ssdev, xtime_ss_ID,&
                  "begin_time", xtime_begin_time))

             call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "title", &
                  "LVT land surface analysis output"))
             call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                  "institution", &
                  trim(LVT_rc%institution)))
             call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "history", &
                  "created on date: "//date(1:4)//"-"//date(5:6)//"-"// &
                  date(7:8)//"T"//stime(1:2)//":"//stime(3:4)//":"// &
                  stime(5:10)))
             call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                  "references", "Kumar_etal_GMD_2012"))
             call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "comment", &
                  "website: http://lis.gsfc.nasa.gov/"))

             !grid information
             if (trim(LVT_rc%domain) .eq."latlon" ) then !latlon
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "MAP_PROJECTION", &
                     "EQUIDISTANT CYLINDRICAL"))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LAT", &
                     LVT_rc%gridDesc(4)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LON", &
                     LVT_rc%gridDesc(5)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DX", &
                     LVT_rc%gridDesc(9)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DY", &
                     LVT_rc%gridDesc(10)))
             elseif (trim(LVT_rc%domain) .eq. "mercator") then
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "MAP_PROJECTION", &
                     "MERCATOR"))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LAT", &
                     LVT_rc%gridDesc(4)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LON", &
                     LVT_rc%gridDesc(5)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "TRUELAT1", &
                     LVT_rc%gridDesc(10)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "STANDARD_LON", &
                     LVT_rc%gridDesc(11)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DX", &
                     LVT_rc%gridDesc(8)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DY", &
                     LVT_rc%gridDesc(9)))
             elseif (trim(LVT_rc%domain) .eq. "lambert") then !lambert conformal
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "MAP_PROJECTION", &
                     "LAMBERT CONFORMAL"))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LAT", &
                     LVT_rc%gridDesc(4)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LON", &
                     LVT_rc%gridDesc(5)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "TRUELAT1", &
                     LVT_rc%gridDesc(10)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "TRUELAT2", &
                     LVT_rc%gridDesc(7)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "STANDARD_LON", &
                     LVT_rc%gridDesc(11)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DX", &
                     LVT_rc%gridDesc(8)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DY", &
                     LVT_rc%gridDesc(9)))

             elseif (trim(LVT_rc%domain) .eq. "polar") then ! polar stereographic
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "MAP_PROJECTION", &
                     "POLAR STEREOGRAPHIC"))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LAT", &
                     LVT_rc%gridDesc(4)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "SOUTH_WEST_CORNER_LON", &
                     LVT_rc%gridDesc(5)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "TRUELAT1", &
                     LVT_rc%gridDesc(10)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "ORIENT", &
                     LVT_rc%gridDesc(7)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                     "STANDARD_LON", &
                     LVT_rc%gridDesc(11)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DX", &
                     LVT_rc%gridDesc(8)))
                call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, "DY", &
                     LVT_rc%gridDesc(9)))
             endif
          end if

          dataEntry => LVT_histData%head_ds1_list

          do while(associated(dataEntry))
             !reset the pointers to the head of the linked list
             if (LVT_LIS_rc(1)%anlys_data_class .eq. "LSM") then
                lisdataEntry => LVT_LISoutput(1)%head_lsm_list
             elseif (LVT_LIS_rc(1)%anlys_data_class .eq. "Routing") then
                lisdataEntry => LVT_LISoutput(1)%head_routing_list
             elseif (LVT_LIS_rc(1)%anlys_data_class .eq. "RTM") then
                lisdataEntry => LVT_LISoutput(1)%head_rtm_list
             elseif (LVT_LIS_rc(1)%anlys_data_class .eq. "Irrigation") then
                lisdataEntry => LVT_LISoutput(1)%head_irrig_list
             endif
             do while(associated(lisdataEntry))
                if (lisdataEntry%short_name .eq. dataEntry%short_name) then

                   call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

                   if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
                        LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
                      call defineNETCDFheaderVar_ss(ftn_ssdev, dimID, &
                           lisdataEntry)
                   end if

                endif
                lisdataEntry => lisdataEntry%next
             enddo
             dataEntry => dataEntry%next
          enddo

          ! EMK:  Include number of soil layers and soil layer thicknesses
          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
               "NUM_SOIL_LAYERS", &
               nsoillayers), &
               'nf90_put_att for title failed in LVT_DataStreamsMod')
          call LVT_verify(nf90_put_att(ftn_mean, NF90_GLOBAL, &
               "SOIL_LAYER_THICKNESSES", &
               lyrthk), &
               'nf90_put_att for title failed in LVT_DataStreamsMod')

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call LVT_verify(nf90_put_att(ftn_ssdev, NF90_GLOBAL, &
                  "NUM_SOIL_LAYERS", &
                  nsoillayers), &
                  'nf90_put_att for title failed in LVT_DataStreamsMod')
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL, &
                  "SOIL_LAYER_THICKNESSES", &
                  lyrthk), &
                  'nf90_put_att for title failed in LVT_DataStreamsMod')
          end if

          ! EMK FIXME...Replace HYCOM with NAVGEM
          if (LVT_rc%processHYCOM .eq. 1) then

             ! First, handle water_temp
             call LVT_verify(nf90_def_var(ftn_mean, &
                  trim(LVT_histData%watertemp%short_name), &
                  nf90_float, &
                  dimids=dimID(1:2), &
                  varID=LVT_histData%watertemp%varId_def), &
                  'nf90_def_var for '// &
                  trim(LVT_histData%watertemp%short_name)// &
                  'failed in defineNETCDFheadervar')

#if(defined USE_NETCDF4)
             call LVT_verify(nf90_def_var_deflate(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  shuffle, deflate, deflate_level), &
                  'nf90_def_var_deflate for '// &
                  trim(LVT_histData%watertemp%short_name)// &
                  'failed in defineNETCDFheadervar')
#endif
             !EMK...Add variable attributes
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "units", &
                  trim(LVT_histData%watertemp%units)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "standard_name", &
                  trim(LVT_histData%watertemp%standard_name)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "long_name", &
                  trim(LVT_histData%watertemp%long_name)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "scale_factor", 1.0))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "add_offset", 0.0))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "missing_value", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "_FillValue", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "vmin", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%watertemp%varId_def, &
                  "vmax", LVT_rc%udef))

             ! Next, add aice
             call LVT_verify(nf90_def_var(ftn_mean, &
                  trim(LVT_histData%aice%short_name), &
                  nf90_float, &
                  dimids=dimID(1:2), &
                  varID=LVT_histData%aice%varId_def), &
                  'nf90_def_var for '// &
                  trim(LVT_histData%aice%short_name)// &
                  'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
             call LVT_verify(nf90_def_var_deflate(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  shuffle, deflate, deflate_level), &
                  'nf90_def_var_deflate for '// &
                  trim(LVT_histData%aice%short_name)// &
                  'failed in defineNETCDFheadervar')
#endif
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "units", &
                  trim(LVT_histData%aice%units)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "standard_name", &
                  trim(LVT_histData%aice%standard_name)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "long_name", &
                  trim(LVT_histData%aice%long_name)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "scale_factor", 1.0))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "add_offset", 0.0))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "missing_value", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "_FillValue", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "vmin", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%aice%varId_def, &
                  "vmax", LVT_rc%udef))

             ! Next, add hi
             call LVT_verify(nf90_def_var(ftn_mean, &
                  trim(LVT_histData%hi%short_name), &
                  nf90_float, &
                  dimids=dimID(1:2), &
                  varID=LVT_histData%hi%varId_def), &
                  'nf90_def_var for '// &
                  trim(LVT_histData%hi%short_name)// &
                  'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
             call LVT_verify(nf90_def_var_deflate(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  shuffle, deflate, deflate_level), &
                  'nf90_def_var_deflate for '// &
                  trim(LVT_histData%hi%short_name)// &
                  'failed in defineNETCDFheadervar')
#endif
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "units", &
                  trim(LVT_histData%hi%units)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "standard_name", &
                  trim(LVT_histData%hi%standard_name)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "long_name", &
                  trim(LVT_histData%hi%long_name)))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "scale_factor", 1.0))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "add_offset", 0.0))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "missing_value", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "_FillValue", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "vmin", LVT_rc%udef))
             call LVT_verify(nf90_put_att(ftn_mean, &
                  LVT_histData%hi%varId_def, &
                  "vmax", LVT_rc%udef))
          endif

          ! EMK...Add additional PS41 snow variable headers.
          if (jules_ps41_ens_snow) then

             lisdataEntry => SWE
             call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

             lisdataEntry => SnowDensity
             call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

             lisdataEntry => LayerSnowDensity
             call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

             lisdataEntry => SnowGrain
             call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

             lisdataEntry => SnowDepth
             call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

             lisdataEntry => GrndSnow
             call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

             lisdataEntry => SurftSnow
             call defineNETCDFheaderVar(ftn_mean, dimID, lisdataEntry)

          end if
          ! EMK END PS41 snow headers

          call LVT_verify(nf90_enddef(ftn_mean))
          call LVT_verify(nf90_put_var(ftn_mean, xtimeID, 0.0))

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call LVT_verify(nf90_enddef(ftn_ssdev))
             call LVT_verify(nf90_put_var(ftn_ssdev, xtime_ss_ID, 0.0))
          end if

          ! EMK...lat/lon calculated above for all output file types.
          call LVT_verify(nf90_put_var(ftn_mean, xlatID, &
               lat, (/1, 1/), &
               (/LVT_rc%gnc, LVT_rc%gnr/)), &
               'nf90_put_var failed for lat')
          call LVT_verify(nf90_put_var(ftn_mean, xlonID, &
               lon, (/1, 1/), &
               (/LVT_rc%gnc, LVT_rc%gnr/)), &
               'nf90_put_var failed for lon')

          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call LVT_verify(nf90_put_var(ftn_ssdev, xlat_ss_ID, &
                  lat, (/1, 1/),&
                  (/LVT_rc%gnc, LVT_rc%gnr/)),&
                  'nf90_put_var failed for lat')
             call LVT_verify(nf90_put_var(ftn_ssdev, xlon_ss_ID, &
                  lon, (/1, 1/),&
                  (/LVT_rc%gnc, LVT_rc%gnr/)),&
                  'nf90_put_var failed for lon')
          end if
       endif

       ! EMK Output updated PS41 snow variables not read in from LIS file
       if (jules_ps41_ens_snow) then
          if (LVT_rc%lvt_out_format .eq. "netcdf") then

             gtmp1_1d = 0.0
             call LVT_fetch_jules_ps41_ens_snow_final( &
                  LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d, &
                  1, "SWE_inst", is_ps41_snow_var)
             call writeSingleNetcdfVar(ftn_mean, &
                  gtmp1_1d, &
                  SWE%varid_def, &
                  1)

             gtmp1_1d = 0.0
             call LVT_fetch_jules_ps41_ens_snow_final( &
                  LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d,  &
                  1, "SnowDensity_inst", is_ps41_snow_var)
             call writeSingleNetcdfVar(ftn_mean, &
                  gtmp1_1d, &
                  SnowDensity%varid_def, &
                  1)

             do k = 1, 3
                gtmp1_1d = 0.0
                call LVT_fetch_jules_ps41_ens_snow_final( &
                     LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d, &
                     k, "LayerSnowDensity_inst", is_ps41_snow_var)
                call writeSingleNetcdfVar(ftn_mean, &
                     gtmp1_1d, &
                     LayerSnowDensity%varid_def, &
                     k)
             end do

             gtmp1_1d = 0.0
             call LVT_fetch_jules_ps41_ens_snow_final( &
                  LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d, &
                  1, "SnowGrain_inst", is_ps41_snow_var)
             call writeSingleNetcdfVar(ftn_mean, &
                  gtmp1_1d, &
                  SnowGrain%varid_def, &
                  1)

             gtmp1_1d = 0.0
             call LVT_fetch_jules_ps41_ens_snow_final( &
                  LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d, &
                  1, "SnowDepth_inst", is_ps41_snow_var)
             call writeSingleNetcdfVar(ftn_mean, &
                  gtmp1_1d, &
                  SnowDepth%varid_def, &
                  1)

             gtmp1_1d = 0.0
             call LVT_fetch_jules_ps41_ens_snow_final( &
                  LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d, &
                  1, "GrndSnow_inst", is_ps41_snow_var)
             call writeSingleNetcdfVar(ftn_mean,&
                  gtmp1_1d, &
                  GrndSnow%varid_def, &
                  1)

             gtmp1_1d = 0.0
             call LVT_fetch_jules_ps41_ens_snow_final( &
                  LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d, &
                  1, "SurftSnow_inst", is_ps41_snow_var)
             call writeSingleNetcdfVar(ftn_mean, &
                  gtmp1_1d, &
                  SurftSnow%varid_def, &
                  1)

             ! Cleanup
             call LVT_deallocate_metadata(SWE)
             call LVT_deallocate_metadata(SnowDensity)
             call LVT_deallocate_metadata(LayerSnowDensity)
             call LVT_deallocate_metadata(SnowGrain)
             call LVT_deallocate_metadata(SnowDepth)
             call LVT_deallocate_metadata(GrndSnow)
             call LVT_deallocate_metadata(SurftSnow)

          else
             write(LVT_logunit,*) &
                  '[ERR] GRIB OUTPUT NOT SUPPORTED YET FOR PS41 SNOW'
             flush(LVT_logunit)
             stop
          end if
       end if

       dataEntry => LVT_histData%head_ds1_list

       do while (associated(dataEntry))
!reset the pointers to the head of the linked list
          if (LVT_LIS_rc(1)%anlys_data_class .eq. "LSM") then
             lisdataEntry => LVT_LISoutput(1)%head_lsm_list
          elseif (LVT_LIS_rc(1)%anlys_data_class .eq. "Routing") then
             lisdataEntry => LVT_LISoutput(1)%head_routing_list
          elseif (LVT_LIS_rc(1)%anlys_data_class .eq. "RTM") then
             lisdataEntry => LVT_LISoutput(1)%head_rtm_list
          elseif (LVT_LIS_rc(1)%anlys_data_class .eq. "Irrigation") then
             lisdataEntry => LVT_LISoutput(1)%head_irrig_list
          endif
          do while (associated(lisdataEntry))

             if (lisdataEntry%short_name .eq. dataEntry%short_name) then

                ! Set timerange indicator equal to 133 for AFWA's
                ! specifications for surface runoff, baseflow, and total
                ! precipitation to make the LIS-7 output match the LIS-6
                ! style. - dmm

                ! EMK...Revised settings based on name of variable
                if (index(trim(dataEntry%short_name), "_max") .gt. 0) then
                   stepType = "max"
                   timeRange = 7
                   pdTemplate = 12 ! Derived fcsts from ensemble over time interval
                else if (index(trim(dataEntry%short_name), "_min") .gt. 0) then
                   stepType = "min"
                   timeRange = 7
                   pdTemplate = 12 ! Derived fcsts from ensemble over time interval
                else if (dataEntry%timeAvgOpt .eq. 0) then
                   stepType = "instant"
                   timeRange = 1
                   pdTemplate = 2 ! Derived fcst from ensemble at point in time

                else if (dataEntry%timeAvgOpt .eq. 1 .or. &
                     dataEntry%timeAvgOpt .eq. 2) then
                   stepType = "avg"
                   timeRange = 7
                   pdTemplate = 12 ! Derived fcsts from ensemble over time interval
                else if (dataEntry%timeAvgOpt .eq. 3) then
                   stepType = "accum"
                   timeRange = 7 ! "between first and second"
                   pdTemplate = 12 ! Derived fcsts from ensemble over time interval
                else
                   write(LVT_logunit,*)'[ERR] Cannot handle ', &
                        trim(dataEntry%short_name)
                   call LVT_endrun()
                end if

                if ((lisdataEntry%index .eq. LVT_LIS_MOC_QS(1)) .or.   &
                     (lisdataEntry%index .eq. LVT_LIS_MOC_QSB(1)) .or.   &
                     (lisdataEntry%index .eq. LVT_LIS_MOC_TOTALPRECIP(1))) then
                   ! EMK...GRIB1 only
                   if(LVT_rc%lvt_out_format .ne. "grib2") then
                      timeRange = 133
                   end if
                endif

                !EMK...Special handling for RHMin, which is an extreme
                ! (minimum) value.
                if (trim(dataEntry%short_name) == "RHMin") then
                   stepType = "min"
                   timeRange = 7
                   pdTemplate = 12
                end if

                ! EMK...Reworked ensemble statistics code.  Allow application
                ! of noises smoother to each ensemble member *before*
                ! calculating ensemble mean and spread.
                do k = 1, dataEntry%vlevels
                   gtmp1_1d(:) = 0.0
                   ngtmp1_1d(:) = 0
                   gtmp1_ss(:) = 0.0

                   ! EMK...Special handling for JULES PS41 snow variables.
                   ! In this case, we do not take raw ensemble means, but
                   ! instead apply a JULES-based relayering.  This calculation
                   ! was done higher up; here we pull the requested variable
                   ! for output to file.
                   is_ps41_snow_var = .false.
                   if (jules_ps41_ens_snow) then

                      call LVT_fetch_jules_ps41_ens_snow_final( &
                           LVT_rc%lnc, LVT_rc%lnr, gtmp1_1d, &
                           k, trim(dataEntry%short_name)//"_inst", &
                           is_ps41_snow_var)

                      ! Not all PS41 variables involve snow.  Check to
                      ! see if this did; if it didn't, normal ensemble
                      ! post-processing will occur later down.
                      if (is_ps41_snow_var) then

                         ! Only write ensemble mean for PS41 snow variables
                         if (LVT_rc%lvt_out_format .eq. "grib2") then

                            call writeSingleGrib2Var(ftn_mean, &
                                 gtmp1_1d, &
                                 lisdataentry%varid_def, &
                                 lisdataentry%gribSF, &
                                 lisdataentry%gribSfc, &
                                 lisdataentry%gribLvl, &
                                 lisdataentry%gribDis, &
                                 lisdataentry%gribCat, &
                                 pdTemplate, &
                                 stepType, &
                                 time_unit, &
                                 time_past, &
                                 time_curr, &
                                 timeRange, &
                                 k, &
                                 toplev(k:k), &
                                 botlev(k:k), &
                                 depscale(k:k), &
                                 typeOfGeneratingProcess=4, &
                                 typeOfProcessedData=4)

                         elseif (LVT_rc%lvt_out_format .eq. "grib1") then
                            call writeSingleGrib1Var(ftn_mean, &
                                 gtmp1_1d, &
                                 lisdataentry%varid_def, &
                                 lisdataentry%gribSF, &
                                 lisdataentry%gribSfc, &
                                 lisdataentry%gribLvl, &
                                 stepType, &
                                 time_unit, &
                                 time_past, &
                                 time_curr, &
                                 timeRange, &
                                 k, &
                                 toplev(k:k), &
                                 botlev(k:k))

                         elseif (LVT_rc%lvt_out_format .eq. "netcdf") then

                            call writeSingleNetcdfVar(ftn_mean, &
                                 gtmp1_1d, &
                                 lisdataentry%varid_def, &
                                 k)

                         end if ! output format

                         ! If we processed a PS41 ensemble snow variable, we
                         ! don't need to continue to normal ensemble
                         ! processing. Just go to the next vertical level.
                         cycle

                      end if ! if PS41 snow variable
                   end if ! If processing JULES PS41 snow ensembles.

                   ! Normal ensemble postprocessing starts here.
                   do m = 1, LVT_rc%nensem

                      ! Must initialize ensemble member with "undefined" for
                      ! noise smoother
                      gtmp1_1d_mem(:) = LVT_rc%udef
                      do r = 1, LVT_rc%lnr
                         do c = 1, LVT_rc%lnc
                            if (LVT_domain%gindex(c,r) .ne. -1) then
                               gid = LVT_domain%gindex(c,r)
                               gtmp1_1d_mem(c + (r-1)*LVT_rc%lnc) = &
                                    dataEntry%value(gid,m,k)
                            endif
                         enddo ! c
                      enddo ! r

                      ! Apply the smoother
                      ! EMK...Removed the hardwired exceptions to
                      ! smoothing.  The original exception list did not
                      ! consider forcing perturbations.  It seams best
                      ! to just trust the setting in the lvt.config file.
                      ! EMK...Restored exception list for categorical
                      ! variables, since smoothing makes no physical sense
                      if (.not. ( &
                           (dataEntry%short_name .eq. "Landcover") .or. &
                           (dataEntry%short_name .eq. "Landmask") .or. &
                           (dataEntry%short_name .eq. "Soiltype"))) then
                         if (LVT_rc%applyNoiseReductionFilter .eq. 1) then
                            call applyNoiseReductionFilter(gtmp1_1d_mem)
                         end if
                      end if

                      ! Now provide smoothed field to ensemble mean and
                      ! spread
                      do r = 1,LVT_rc%lnr
                         do c = 1,LVT_rc%lnc
                            if (LVT_domain%gindex(c,r) .ne. -1) then
                               gid = LVT_domain%gindex(c,r)

                               if (LVT_rc%nensem > 1) then
                                  ! Use Welford algorithm to calculate
                                  ! mean and standard deviation
                                  count = ngtmp1_1d(c + (r-1)*LVT_rc%lnc)
                                  mean = gtmp1_1d(c + (r-1)*LVT_rc%lnc)
                                  m2 = gtmp1_ss(c + (r-1)*LVT_rc%lnc)
                                  new_value = &
                                       gtmp1_1d_mem(c + (r-1)*LVT_rc%lnc)
                                  call welford_update(count, mean, m2, &
                                       new_value)
                                  ngtmp1_1d(c + (r-1)*LVT_rc%lnc) = count
                                  gtmp1_1d(c + (r-1)*LVT_rc%lnc) = mean
                                  gtmp1_ss(c + (r-1)*LVT_rc%lnc) = m2
                               else
                                  gtmp1_1d(c + (r-1)*LVT_rc%lnc) = &
                                       gtmp1_1d(c + (r-1)*LVT_rc%lnc) + &
                                       gtmp1_1d_mem(c + (r-1)*LVT_rc%lnc)
                                  ngtmp1_1d(c + (r-1)*LVT_rc%lnc) = &
                                       ngtmp1_1d(c + (r-1)*LVT_rc%lnc) + 1
                               end if
                            endif
                         enddo ! c
                      enddo ! r

                   end do ! m

                   ! Finalize the ensemble mean and spread
                   do r = 1, LVT_rc%lnr
                      do c = 1, LVT_rc%lnc
                         if (LVT_rc%nensem > 1) then
                            ! Use Welford algorithm to calculate mean and
                            ! standard deviation
                            count = ngtmp1_1d(c + (r-1)*LVT_rc%lnc)
                            if (count < 1) then
                               gtmp1_1d(c + (r-1)*LVT_rc%lnc) = LVT_rc%udef
                               gtmp1_ss(c + (r-1)*LVT_rc%lnc) = LVT_rc%udef
                            else
                               mean = gtmp1_1d(c + (r-1)*LVT_rc%lnc)
                               m2 = gtmp1_ss(c + (r-1)*LVT_rc%lnc)
                               call welford_finalize(count, mean, m2, stddev)
                               gtmp1_1d(c + (r-1)*LVT_rc%lnc) = mean
                               gtmp1_ss(c + (r-1)*LVT_rc%lnc) = stddev
                            end if
                         else ! No ensembles, just calculate mean
                            if(ngtmp1_1d(c + (r-1)*LVT_rc%lnc).gt.0) then
                               gtmp1_1d(c + (r-1)*LVT_rc%lnc) = &
                                    gtmp1_1d(c + (r-1)*LVT_rc%lnc)/&
                                    ngtmp1_1d(c + (r-1)*LVT_rc%lnc)
                            else
                               gtmp1_1d(c + (r-1)*LVT_rc%lnc) = LVT_rc%udef
                            end if
                            gtmp1_ss(c + (r-1)*LVT_rc%lnc) = LVT_rc%udef
                         end if
                      enddo ! c
                   enddo ! r
                   ! EMK END...k loop ends further down

                   if (LVT_rc%lvt_out_format .eq. "grib2") then

                      call writeSingleGrib2Var(ftn_mean, &
                           gtmp1_1d, &
                           lisdataentry%varid_def, &
                           lisdataentry%gribSF, &
                           lisdataentry%gribSfc, &
                           lisdataentry%gribLvl, &
                           lisdataentry%gribDis, &
                           lisdataentry%gribCat, &
                           pdTemplate, &
                           stepType, &
                           time_unit, &
                           time_past, &
                           time_curr, &
                           timeRange, &
                           k, &
                           toplev(k:k), &
                           botlev(k:k), &
                           depscale(k:k), &
                           typeOfGeneratingProcess=4, &
                           typeOfProcessedData=4)

                      if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
                           LVT_rc%nensem > 1 &
                           .and. .not. jules_ps41_ens_snow) then
                         call writeSingleGrib2Var(ftn_ssdev, &
                              gtmp1_ss, &
                              lisdataentry%varid_def, &
                              lisdataentry%gribSF, &
                              lisdataentry%gribSfc, &
                              lisdataentry%gribLvl, &
                              lisdataentry%gribDis, &
                              lisdataentry%gribCat, &
                              pdTemplate, &
                              stepType, &
                              time_unit, &
                              time_past, &
                              time_curr, &
                              timeRange, &
                              k, &
                              toplev(k:k), &
                              botlev(k:k), &
                              depscale(k:k), &
                              ensembleSpread=.true., &
                              typeOfGeneratingProcess=4, &
                              typeOfProcessedData=4)
                      end if

                   elseif(LVT_rc%lvt_out_format.eq."grib1") then
                      call writeSingleGrib1Var(ftn_mean, &
                           gtmp1_1d, &
                           lisdataentry%varid_def, &
                           lisdataentry%gribSF, &
                           lisdataentry%gribSfc, &
                           lisdataentry%gribLvl, &
                           stepType, &
                           time_unit, &
                           time_past, &
                           time_curr, &
                           timeRange, &
                           k, &
                           toplev(k:k), &
                           botlev(k:k))

                      if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
                           LVT_rc%nensem > 1 &
                           .and. .not. jules_ps41_ens_snow) then

                         call writeSingleGrib1Var(ftn_ssdev, &
                              gtmp1_ss, &
                              lisdataentry%varid_def, &
                              lisdataentry%gribSF, &
                              lisdataentry%gribSfc, &
                              lisdataentry%gribLvl, &
                              stepType, &
                              time_unit, &
                              time_past, &
                              time_curr, &
                              timeRange, &
                              k, &
                              toplev(k:k), &
                              botlev(k:k))
                      end if

                   elseif (LVT_rc%lvt_out_format .eq. "netcdf") then
                      call writeSingleNetcdfVar(ftn_mean, &
                           gtmp1_1d, &
                           lisdataentry%varid_def, &
                           k)
                      if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
                           LVT_rc%nensem > 1 &
                           .and. .not. jules_ps41_ens_snow) then
                         call writeSingleNetcdfVar(ftn_ssdev, &
                              gtmp1_ss, &
                              lisdataentry%varid_ss, &
                              k)
                      end if
                   endif

                enddo ! k
                exit
             endif
             lisdataEntry => lisdataEntry%next
          enddo
          dataEntry => dataEntry%next
       enddo

       ! Free up memory for PS41 ensemble snow postprocessing.
       if (jules_ps41_ens_snow) then
          call LVT_cleanup_jules_ps41_ens_snow()
       end if

       ! EMK...Use HYCOM for sea ice, and NAVGEM for SST.
       ! EMK 20220519...Reinstate HYCOM SST.
       call LVT_append_HYCOM_fields(ftn_mean, &
          time_unit, &
          time_past, &
          time_curr, &
          timeRange, &
          toplev(1), &
          botlev(1), &
          lat, lon)
       ! call LVT_append_HYCOM_cice_fields(ftn_mean, &
       !    time_unit, &
       !    time_past, &
       !    time_curr, &
       !    timeRange, &
       !    toplev(1), &
       !    botlev(1), &
       !    lat, lon)
       ! call LVT_append_navgem_sst_field(ftn_mean, &
       !       time_unit, &
       !       time_past, &
       !       time_curr, &
       !       timeRange, &
       !       toplev(1), &
       !       botlev(1))
       if (LVT_rc%lvt_out_format .eq. "grib1") then
          call grib_close_file(ftn_mean, iret)
          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call grib_close_file(ftn_ssdev, iret)
          end if
       elseif (LVT_rc%lvt_out_format .eq. "grib2") then
          call grib_close_file(ftn_mean, iret)
          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call grib_close_file(ftn_ssdev, iret)
          end if
       elseif (LVT_rc%lvt_out_format .eq. "netcdf") then
          call LVT_verify(nf90_close(ftn_mean))
          if (LVT_rc%tavgInterval == LVT_rc%ts .and. &
               LVT_rc%nensem > 1 .and. .not. jules_ps41_ens_snow) then
             call LVT_verify(nf90_close(ftn_ssdev))
          end if
       endif
    endif

  end subroutine LVT_writeDataStreams

  ! EMK...Return logical indicating if alarm should ring.
  ! Used by "557 post" runmode.
  logical function alarm_is_on() result(alarmCheck)
     use LVT_timeMgrMod,      only : LVT_get_julhr
     implicit none

     logical, save :: firstTime = .true.
     integer, save :: starttime = 0
     integer :: curtime
     integer :: difftime

     alarmCheck = .false.

     if (firstTime) then
        call LVT_get_julss(LVT_rc%syr, LVT_rc%smo, LVT_rc%sda, &
             LVT_rc%shr, LVT_rc%smn, LVT_rc%sss, &
             starttime)
        firstTime = .false.
     end if
     call LVT_get_julss(LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
          LVT_rc%hr, LVT_rc%mn, LVT_rc%ss, &
          curtime)
     difftime = curtime - starttime
     if (difftime .gt. 0) then
        if (mod(difftime, LVT_rc%statswriteint) .eq. 0) then
           alarmCheck = .true.
        end if
     end if

     return
  end function alarm_is_on

  ! Add NAVGEM fields to output file.
  subroutine LVT_append_navgem_sst_field(ftn_mean, time_unit, time_past, &
       time_curr, timeRange, toplev, botlev)

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: ftn_mean
    integer, intent(in) :: time_unit
    integer, intent(in) :: time_past
    integer, intent(in) :: time_curr
    integer, intent(in) :: timeRange
    real, intent(in) :: toplev(1)
    real, intent(in) :: botlev(1)

    ! Locals
    character(250) :: navgem_sst_fname
    real :: gridDesci(50) ! Full NAVGEM grid
    character(10) :: cdate
    logical :: file_exists
    real, allocatable :: sst(:)
    real, allocatable :: cice(:)
    real, allocatable :: icethick(:)
    integer :: npts
    real :: rlat(LVT_rc%lnc * LVT_rc%lnr)
    real :: rlon(LVT_rc%lnc * LVT_rc%lnr)
    integer :: n11(LVT_rc%lnc * LVT_rc%lnr)
    integer :: n12(LVT_rc%lnc * LVT_rc%lnr)
    integer :: n21(LVT_rc%lnc * LVT_rc%lnr)
    integer :: n22(LVT_rc%lnc * LVT_rc%lnr)
    real :: w11(LVT_rc%lnc * LVT_rc%lnr)
    real :: w12(LVT_rc%lnc * LVT_rc%lnr)
    real :: w21(LVT_rc%lnc * LVT_rc%lnr)
    real :: w22(LVT_rc%lnc * LVT_rc%lnr)
    real :: interp_var(LVT_rc%lnc * LVT_rc%lnr)
    logical*1, allocatable :: li(:)
    logical*1 :: lo(LVT_rc%lnc * LVT_rc%lnr)
    integer :: mi, mo
    integer :: year, month, day, hour, fcst_hr
    real :: udef
    integer :: ivar
    integer :: iret
    integer :: gribSF, gribSfc, gribLvl, gribCat, gribDis
    character*10 :: stepType
    integer :: pdTemplate
    integer :: varid_def
    real :: depscale(1)
    real, allocatable :: thin_latitudes(:,:)

    external :: bilinear_interp_input
    external :: bilinear_interp

    if (LVT_rc%processHYCOM .ne. 1) return

    ! Check for SST GRIB file.  (This actually contains merged sea surface
    ! temperature and land surface temperature; we treat as SST for
    ! simplicity.)
    call LVT_get_navgem_sst_gr1_filename(navgem_sst_fname, &
         year, month, day, hour, fcst_hr)
    !call LVT_get_navgem_sst_bin_filename(navgem_sst_fname, &
    !     year, month, day, hour)
    if (trim(navgem_sst_fname) .eq. "NONE") then
       file_exists = .false.
    else
       file_exists = .true.
    end if
    if (.not. file_exists) then
       write(LVT_logunit,*) '[INFO] No NAVGEM fields to append!'
       return
    end if

    ! Fetch SST from the NAVGEM file.
    call LVT_fetch_navgem_sst_gr1_field(navgem_sst_fname, sst, gridDesci)
    !call LVT_fetch_navgem_sst_bin_field(navgem_sst_fname, sst, gridDesci)

    ! Prepare to interpolate.
    npts = LVT_rc%lnc*LVT_rc%lnr
    call bilinear_interp_input(gridDesci, LVT_rc%gridDesc, npts, &
         rlat, rlon, n11, n12, n21, n22, &
         w11, w12, w21, w22)
    allocate(li(size(sst)))
    li = .true.
    mo = npts
    udef = -9999.
    interp_var = udef
    li = .true.
    lo = .true.

    ! Interpolate the SST
    call bilinear_interp(LVT_rc%gridDesc, li, sst, lo, interp_var, &
         size(sst), size(interp_var), rlat, rlon, &
         w11, w12, w21, w22, &
         n11, n12, n21, n22, udef, iret)

    ! Prepare output field settings.
    gribDis   = 10
    stepType = "instant"
    pdTemplate = 0
    gribCat   = 3
    varid_def = 0
    gribSfc   = 1
    gribSF    = 10
    gribLvl   = 1

    ! Now write the interpolated field to output
    if (LVT_rc%lvt_out_format .eq. "grib2") then
       call writeSingleGrib2Var(ftn_mean, &
            interp_var, &
            varid_def, &
            gribSF, &
            gribSfc, &
            gribLvl, &
            gribDis, &
            gribCat, &
            pdTemplate, &
            stepType, &
            time_unit, &
            time_past, &
            time_curr, &
            timeRange, &
            1, &
            toplev(1), &
            botlev(1), &
            depscale(1), &
            typeOfGeneratingProcess=2, &
            typeOfProcessedData=1, &
            ref_year=year, & ! FIXME
            ref_month=month, & ! FIXME
            ref_day=day, &     ! FIXME
            ref_hour=hour, &   ! FIXME
            ref_fcst_hr=fcst_hr) ! FIXME
    else if (LVT_rc%lvt_out_format .eq. "grib1") then
       call writeSingleGrib1Var(ftn_mean, &
            interp_var, &
            varid_def, &
            gribSF, &
            gribSfc, &
            gribLvl, &
            stepType, &
            time_unit, &
            time_past, &
            time_curr, &
            timeRange, &
            1, &
            toplev(1), &
            botlev(1))
    else if (LVT_rc%lvt_out_format .eq. "netcdf") then
       call writeSingleNetcdfVar(ftn_mean, &
            interp_var, &
            LVT_histData%watertemp%varId_def, &
            1)
    end if

    ! Clean up
    if (allocated(li)) deallocate(li)
    if (allocated(sst)) deallocate(sst)
  end subroutine LVT_append_navgem_sst_field

!BOP
!
! !ROUTINE: LVT_append_HYCOM_cice_fields
! \label{LVT_append_HYCOM_cice_fields}
!
! !INTERFACE:
  !EMK 20220519...Reinstated HYCOM SST.
  !subroutine LVT_append_HYCOM_cice_fields(ftn_mean, time_unit, time_past, &
  !     time_curr, timeRange, toplev, botlev, lat, lon)
  subroutine LVT_append_HYCOM_fields(ftn_mean, time_unit, time_past, &
       time_curr, timeRange, toplev, botlev, lat, lon)

!
! !DESCRIPTION:
!  This subroutine read the water temperature fields from the HYCOM output,
!  reprojects it to the LVT/LIS grid and appends to the grib1 formatted file.
!
!EOP

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    integer                 :: ftn_mean
    integer                 :: time_unit
    integer                 :: time_past
    integer                 :: time_curr
    integer                 :: timeRange
    real                    :: toplev(1)
    real                    :: botlev(1)
    real, intent(in) :: lat(LVT_rc%lnc,LVT_rc%lnr)
    real, intent(in) :: lon(LVT_rc%lnc,LVT_rc%lnr)

    character*100           :: hycom_fname
    character*10            :: cdate
    logical                 :: file_exists
    integer                 :: nid,ios
    integer                 :: c,r,c1,r1,k,cindex,rindex
    integer                 :: watertid
    real                    :: watert_ip(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: lo(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: lb(LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr)
    real                    :: watert(LVT_rc%HYCOM_nc,LVT_rc%HYCOM_nr,1,1)
    real                    :: watert_1d(LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr)

    ! EMK...Support aice_arc
    integer                 :: aice_arc_id
    real                    :: aice_arc_ip(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: &
         aice_arc_lb(LVT_rc%HYCOM_aice_arc_nc*LVT_rc%HYCOM_aice_arc_nr)
    real                    :: &
         aice_arc(LVT_rc%HYCOM_aice_arc_nc,LVT_rc%HYCOM_aice_arc_nr,1,1)
    real                    :: &
         aice_arc_1d(LVT_rc%HYCOM_aice_arc_nc*LVT_rc%HYCOM_aice_arc_nr)

    ! EMK...Support aice_ant
    integer                 :: aice_ant_id
    real                    :: aice_ant_ip(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: &
         aice_ant_lb(LVT_rc%HYCOM_aice_ant_nc*LVT_rc%HYCOM_aice_ant_nr)
    real                    :: &
         aice_ant(LVT_rc%HYCOM_aice_ant_nc,LVT_rc%HYCOM_aice_ant_nr,1,1)
    real                    :: &
         aice_ant_1d(LVT_rc%HYCOM_aice_ant_nc*LVT_rc%HYCOM_aice_ant_nr)

    real                    :: aice_ip(LVT_rc%lnc*LVT_rc%lnr)

    ! EMK...Support hi_arc
    integer                 :: hi_arc_id
    real                    :: hi_arc_ip(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: &
         hi_arc_lb(LVT_rc%HYCOM_hi_arc_nc*LVT_rc%HYCOM_hi_arc_nr)
    real                    :: &
         hi_arc(LVT_rc%HYCOM_hi_arc_nc,LVT_rc%HYCOM_hi_arc_nr,1,1)
    real                    :: &
         hi_arc_1d(LVT_rc%HYCOM_hi_arc_nc*LVT_rc%HYCOM_hi_arc_nr)

    ! EMK...Support hi_ant
    integer                 :: hi_ant_id
    real                    :: hi_ant_ip(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: &
         hi_ant_lb(LVT_rc%HYCOM_hi_ant_nc*LVT_rc%HYCOM_hi_ant_nr)
    real                    :: &
         hi_ant(LVT_rc%HYCOM_hi_ant_nc,LVT_rc%HYCOM_hi_ant_nr,1,1)
    real                    :: &
         hi_ant_1d(LVT_rc%HYCOM_hi_ant_nc*LVT_rc%HYCOM_hi_ant_nr)

    real                    :: hi_ip(LVT_rc%lnc*LVT_rc%lnr)

    character*10            :: stepType
    integer                 :: varid_def
    integer                 :: gribSF, gribSfc,gribLvl,gribCat,gribDis
    real                    :: depscale(1)
    integer                 :: pdTemplate

    integer :: sst_year,sst_month,sst_day, sst_hour, sst_fcst_hr
    integer :: cice_arc_year, cice_arc_month, cice_arc_day, &
         cice_arc_hour, cice_arc_fcst_hr
    integer :: cice_ant_year, cice_ant_month, cice_ant_day, &
         cice_ant_hour, cice_ant_fcst_hr
    integer :: hi_arc_year, hi_arc_month, hi_arc_day, &
         hi_arc_hour, hi_arc_fcst_hr
    integer :: hi_ant_year, hi_ant_month, hi_ant_day, &
         hi_ant_hour, hi_ant_fcst_hr

    integer :: gid

    ! find the filename, open the file, read the field

    if (LVT_rc%processHYCOM .eq. 1) then

       ! *** HANDLE SST ***
      ! write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2)') &
      !      LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr
      ! ! FIXME...Update HYCOM file name convention
      ! hycom_fname = trim(LVT_rc%HYCOMdir)//'/'//&
      !     'hycom_glb_928_'//trim(cdate)//'_t000_ts3z.nc'

      ! watert  = LVT_rc%udef
      ! inquire(file=hycom_fname,exist=file_exists)

      ! if (.not. file_exists) then
      !    write(LVT_logunit,*)'[WARN], missing file ',trim(hycom_fname)
      ! end if

       watert  = LVT_rc%udef
       call get_hycom_sst_filename(hycom_fname, &
            sst_year, sst_month, sst_day, sst_hour, sst_fcst_hr)
       if (trim(hycom_fname) == "NONE") then
          file_exists = .false.
       else
          file_exists = .true.
       end if

       if(file_exists) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          write(LVT_logunit,*) '[INFO] Reading HYCOM data ',trim(hycom_fname)

          ios = nf90_open(path=trim(hycom_fname),mode=NF90_NOWRITE,ncid=nid)
          call LVT_verify(ios, 'Error opening file'//trim(hycom_fname))

          !variable ids
          ios = nf90_inq_varid(nid, 'water_temp',watertid)
          call LVT_verify(ios, 'Error nf90_inq_varid: water_temp')

          !values
          ios = nf90_get_var(nid,watertid, watert,&
               start=(/1,1,1,1/), count=(/LVT_rc%HYCOM_nc,LVT_rc%HYCOM_nr,1,1/))
          call LVT_verify(ios, 'Error nf90_get_var: water_temp')

          ios = nf90_close(nid)
          call LVT_verify(ios, 'Error in nf90_close')
#endif
          watert_1d = -9999.0
          lb = .false.

          do r=1,LVT_rc%HYCOM_nr
             do c=1,LVT_rc%HYCOM_nc
                if(watert(c,r,1,1).ne.-30000) then
                   if(c.gt.2250) then
                      c1 = c-2250
                      r1 = r
                   else
                      c1 = c+2250
                      r1 = r
                   endif
                   !EMK...Change from Celsius to Kelvin
                   !watert_1d(c1+(r1-1)*LVT_rc%HYCOM_nc) = watert(c,r,1,1)*0.001+20.0
                   watert_1d(c1+(r1-1)*LVT_rc%HYCOM_nc) = watert(c,r,1,1)*0.001+20.0+273.15

                   lb(c1+(r1-1)*LVT_rc%HYCOM_nc) = .true.

                endif
             enddo
          enddo

          call upscaleByAveraging(&
               LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
               LVT_rc%HYCOM_n11, lb, &
               watert_1d, lo, watert_ip)

          !EMK: Since SST is missing north of 80N, we need to set water points
          !in this region to a reasonable value.  We follow the typical
          !UKMET SURF value of 271.35K.
          do r = 1, LVT_rc%gnr
             do c = 1, LVT_rc%gnc
                gid = LVT_domain%gindex(c,r)

                if (gid .eq. -1 .and. lat(c,r) >= 80.) then
                   if (watert_ip(c+(r-1)*LVT_rc%gnc) == -9999) then
                      watert_ip(c+(r-1)*LVT_rc%gnc) = 271.35
                   end if
                end if
             end do ! c
          end do ! r

          ! GRIB2 settings...Updated by EMK
          gribDis   = 10
          stepType  = "avg"
          stepType = "instant" ! EMK
          pdTemplate = 0
          gribCat   = 3
          varid_def = 0
          gribSfc   = 1
          gribSF    = 10
          gribLvl   = 1

          if(LVT_rc%lvt_out_format.eq."grib2") then
             !add to the grib file
             call writeSingleGrib2Var(ftn_mean,&
                  watert_ip,&
                  varid_def,&
                  gribSF,&
                  gribSfc,&
                  gribLvl,&
                  gribDis,&
                  gribCat,&
                  pdTemplate,&
                  stepType,&
                  time_unit,&
                  time_past,&
                  time_curr,&
                  timeRange,&
                  1,&
                  toplev(1),&
                  botlev(1),&
                  depscale(1),&
                  typeOfGeneratingProcess=2, &
                  typeOfProcessedData=1, &
                  ref_year=sst_year, &
                  ref_month=sst_month, &
                  ref_day=sst_day, &
                  ref_hour=sst_hour, &
                  ref_fcst_hr=sst_fcst_hr)

          elseif(LVT_rc%lvt_out_format.eq."grib1") then
             call writeSingleGrib1Var(ftn_mean,&
                  watert_ip,&
                  varid_def,&
                  gribSF,&
                  gribSfc,&
                  gribLvl,&
                  stepType,&
                  time_unit,&
                  time_past,&
                  time_curr,&
                  timeRange,&
                  1,&
                  toplev(1),&
                  botlev(1))
          elseif(LVT_rc%lvt_out_format.eq."netcdf") then
             call writeSingleNetcdfVar(ftn_mean,&
                  watert_ip,&
                  LVT_histData%watertemp%varId_def,&
                  1)

          endif
       endif

       ! *** HANDLE AICE_ARC ***

       ! FIXME...Update HYCOM file name convention
!       hycom_fname = trim(LVT_rc%HYCOMdir)//'/'//&
!           'hycom-cice_inst_ARCu0.08_928_'//trim(cdate)//'_t000.nc'
!
!       aice_arc  = LVT_rc%udef
!       inquire(file=hycom_fname,exist=file_exists)
!
!       watert  = LVT_rc%udef

       aice_arc  = LVT_rc%udef
       call get_hycom_cice_filename('ARC', hycom_fname, &
            cice_arc_year, cice_arc_month, cice_arc_day, cice_arc_hour, &
            cice_arc_fcst_hr)
       if (trim(hycom_fname) == "NONE") then
          file_exists = .false.
       else
          file_exists = .true.
       end if

       aice_arc_1d(:) = -9999
       aice_arc_ip(:) = -9999

       if (file_exists) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          write(LVT_logunit,*) '[INFO] Reading HYCOM data ', trim(hycom_fname)

          ios = nf90_open(path=trim(hycom_fname), mode=NF90_NOWRITE, ncid=nid)
          call LVT_verify(ios, 'Error opening file'//trim(hycom_fname))

!variable ids
          ios = nf90_inq_varid(nid, 'aice', aice_arc_id)
          call LVT_verify(ios, 'Error nf90_inq_varid: aice')

!values
          ios = nf90_get_var(nid,aice_arc_id, aice_arc,&
               start=(/1, 1, 1, 1/), &
               count=(/LVT_rc%HYCOM_aice_arc_nc, LVT_rc%HYCOM_aice_arc_nr, &
               1, 1/))
          call LVT_verify(ios, 'Error nf90_get_var: aice_arc')

          ios = nf90_close(nid)
          call LVT_verify(ios, 'Error in nf90_close')
#endif
          aice_arc_1d = -9999.0
          aice_arc_lb = .false.

          do r = 1, LVT_rc%HYCOM_aice_arc_nr
             do c = 1, LVT_rc%HYCOM_aice_arc_nc
                if (aice_arc(c,r,1,1) .ne. -30000) then
                   c1 = c
                   r1 = r
                   aice_arc_1d(c1 + (r1-1)*LVT_rc%HYCOM_aice_arc_nc) = &
                        aice_arc(c,r,1,1)*0.0001

                   aice_arc_lb(c1 + (r1-1)*LVT_rc%HYCOM_aice_arc_nc) = .true.

                endif
             enddo
          enddo

          call upscaleByAveraging( &
               LVT_rc%HYCOM_aice_arc_nc*LVT_rc%HYCOM_aice_arc_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
               LVT_rc%HYCOM_aice_arc_n11, aice_arc_lb, &
               aice_arc_1d, lo, aice_arc_ip)

       end if

       ! FIXME...Update HYCOM file name convention
!       hycom_fname = trim(LVT_rc%HYCOMdir)//'/'//&
!           'hycom-cice_inst_ANTu0.08_928_'//trim(cdate)//'_t000.nc'
!
!       aice_ant  = LVT_rc%udef
!       inquire(file=hycom_fname,exist=file_exists)
!
!       if (.not. file_exists) then
!          write(LVT_logunit,*)'[WARN], missing file ',trim(hycom_fname)
!       end if

       aice_ant = LVT_rc%udef
       call get_hycom_cice_filename('ANT', hycom_fname, &
            cice_ant_year, cice_ant_month, cice_ant_day, &
            cice_ant_hour, cice_ant_fcst_hr)
       if (trim(hycom_fname) == "NONE") then
          file_exists = .false.
       else
          file_exists = .true.
       end if

       aice_ant_1d = -9999.0
       aice_ant_ip(:) = -9999

       if (file_exists) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          write(LVT_logunit,*) '[INFO] Reading HYCOM data ', trim(hycom_fname)

          ios = nf90_open(path=trim(hycom_fname), mode=NF90_NOWRITE, ncid=nid)
          call LVT_verify(ios, 'Error opening file'//trim(hycom_fname))

!variable ids
          ios = nf90_inq_varid(nid, 'aice', aice_ant_id)
          call LVT_verify(ios, 'Error nf90_inq_varid: aice')

!values
          ios = nf90_get_var(nid, aice_ant_id, aice_ant,&
               start=(/1, 1, 1, 1/), &
               count=(/LVT_rc%HYCOM_aice_ant_nc, LVT_rc%HYCOM_aice_ant_nr, &
               1, 1/))
          call LVT_verify(ios, 'Error nf90_get_var: aice_ant')

          ios = nf90_close(nid)
          call LVT_verify(ios, 'Error in nf90_close')
#endif
          aice_ant_1d = -9999.0
          aice_ant_lb = .false.

          do r = 1, LVT_rc%HYCOM_aice_ant_nr
             do c = 1, LVT_rc%HYCOM_aice_ant_nc
                if (aice_ant(c,r,1,1) .ne. -30000) then
                   c1 = c
                   r1 = r
                   aice_ant_1d(c1 + (r1-1)*LVT_rc%HYCOM_aice_ant_nc) = &
                        aice_ant(c,r,1,1)*0.0001

                   aice_ant_lb(c1 + (r1-1)*LVT_rc%HYCOM_aice_ant_nc) = .true.
                endif
             enddo
          enddo

          call upscaleByAveraging( &
               LVT_rc%HYCOM_aice_ant_nc*LVT_rc%HYCOM_aice_ant_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
               LVT_rc%HYCOM_aice_ant_n11, aice_ant_lb, &
               aice_ant_1d, lo, aice_ant_ip)

       end if

       ! Merge the two interpolated aice fields together.
       aice_ip(:) = -9999
       do c = 1,LVT_rc%lnc*LVT_rc%lnr
          if (aice_ant_ip(c) > -9999) then
             aice_ip(c) = aice_ant_ip(c)
          else if (aice_arc_ip(c) > -9999) then
             aice_ip(c) = aice_arc_ip(c)
          end if
       end do ! c

       ! EMK: Since sea ice is missing north of -49.5N and south of 40N, we
       ! need to set water points in this region to a reasonable value.  We
       ! assume sea ice fraction is zero in this region.
       do r = 1, LVT_rc%gnr
          do c = 1, LVT_rc%gnc
             gid = LVT_domain%gindex(c,r)
             if (gid .eq. -1) then
                if (aice_ip(c + (r-1)*LVT_rc%gnc) == -9999) then
                   aice_ip(c + (r-1)*LVT_rc%gnc) = 0
                end if
             end if
          end do ! c
       end do ! r

       ! GRIB2 Settings...updated by EMK
       gribDis   = 10
       !stepType  = "avg"
       stepType = "instant" ! EMK
       pdTemplate = 0
       gribCat   = 2
       varid_def = 0
       gribSfc   = 1
       gribSF    = 100
       gribLvl   = 1

       if (LVT_rc%lvt_out_format .eq. "grib2" ) then
          ! EMK...Use older cice date/time
          if (cice_ant_year .lt. cice_arc_year .or. &
               cice_ant_month .lt. cice_arc_month .or. &
               cice_ant_day .lt. cice_arc_day .or. &
               cice_ant_hour .lt. cice_arc_hour) then
             cice_arc_year = cice_ant_year
             cice_arc_month = cice_ant_month
             cice_arc_day = cice_ant_day
             cice_arc_hour = cice_arc_hour
             cice_arc_fcst_hr = cice_ant_fcst_hr
          end if
          call writeSingleGrib2Var(ftn_mean, &
               aice_ip, &
               varid_def, &
               gribSF, &
               gribSfc, &
               gribLvl, &
               gribDis, &
               gribCat, &
               pdTemplate, &
               stepType, &
               time_unit, &
               time_past, &
               time_curr, &
               timeRange, &
               1, &
               toplev(1), &
               botlev(1), &
               depscale(1), &
               typeOfGeneratingProcess=2, &
               typeOfProcessedData=1, &
               ref_year=cice_arc_year, &
               ref_month=cice_arc_month,&
               ref_day=cice_arc_day, &
               ref_hour=cice_arc_hour, &
               ref_fcst_hr=cice_arc_fcst_hr)

       elseif (LVT_rc%lvt_out_format .eq. "grib1") then
          call writeSingleGrib1Var(ftn_mean, &
               aice_ip, &
               varid_def, &
               gribSF, &
               gribSfc, &
               gribLvl, &
               stepType, &
               time_unit, &
               time_past, &
               time_curr, &
               timeRange, &
               1, &
               toplev(1), &
               botlev(1))
       elseif (LVT_rc%lvt_out_format .eq. "netcdf") then
          call writeSingleNetcdfVar(ftn_mean, &
               aice_ip, &
               LVT_histData%aice%varId_def, &
               1)

       endif

       ! *** HANDLE HI_ARC ***

!        ! FIXME...Update HYCOM file name convention
!        hycom_fname = trim(LVT_rc%HYCOMdir)//'/'//&
!            'hycom-cice_inst_ARCu0.08_928_'//trim(cdate)//'_t000.nc'

!        hi_arc  = LVT_rc%udef
!        inquire(file=hycom_fname,exist=file_exists)

!        if (.not. file_exists) then
!           write(LVT_logunit,*)'[WARN], missing file ',trim(hycom_fname)
!        end if

       hi_arc = LVT_rc%udef
       call get_hycom_cice_filename('ARC', hycom_fname, &
            hi_arc_year, hi_arc_month, hi_arc_day, &
            hi_arc_hour, hi_arc_fcst_hr)

       if (trim(hycom_fname) == "NONE") then
          file_exists = .false.
       else
          file_exists = .true.
       end if

       hi_arc_1d(:) = -9999
       hi_arc_ip(:) = -9999

       if(file_exists) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          write(LVT_logunit,*) '[INFO] Reading HYCOM data ', trim(hycom_fname)

          ios = nf90_open(path=trim(hycom_fname), mode=NF90_NOWRITE, ncid=nid)
          call LVT_verify(ios, 'Error opening file'//trim(hycom_fname))

!variable ids
          ios = nf90_inq_varid(nid, 'hi', hi_arc_id)
          call LVT_verify(ios, 'Error nf90_inq_varid: hi')

!values
          ios = nf90_get_var(nid,hi_arc_id, hi_arc,&
               start=(/1, 1, 1, 1/), &
               count=(/LVT_rc%HYCOM_hi_arc_nc, LVT_rc%HYCOM_hi_arc_nr, 1, 1/))
          call LVT_verify(ios, 'Error nf90_get_var: hi_arc')

          ios = nf90_close(nid)
          call LVT_verify(ios, 'Error in nf90_close')
#endif
          hi_arc_1d = -9999.0
          hi_arc_lb = .false.

          do r = 1, LVT_rc%HYCOM_hi_arc_nr
             do c = 1, LVT_rc%HYCOM_hi_arc_nc
                if (hi_arc(c,r,1,1) .ne. -30000) then
                   c1 = c
                   r1 = r
                   hi_arc_1d(c1 + (r1-1)*LVT_rc%HYCOM_hi_arc_nc) = &
                        hi_arc(c,r,1,1)*0.001

                   hi_arc_lb(c1 + (r1-1)*LVT_rc%HYCOM_hi_arc_nc) = .true.

                endif
             enddo
          enddo

          call upscaleByAveraging( &
               LVT_rc%HYCOM_hi_arc_nc*LVT_rc%HYCOM_hi_arc_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
               LVT_rc%HYCOM_hi_arc_n11, hi_arc_lb, &
               hi_arc_1d, lo, hi_arc_ip)

       end if

       ! FIXME...Update HYCOM file name convention
!        hycom_fname = trim(LVT_rc%HYCOMdir)//'/'//&
!            'hycom-cice_inst_ANTu0.08_928_'//trim(cdate)//'_t000.nc'

!        hi_ant  = LVT_rc%udef
!        inquire(file=hycom_fname,exist=file_exists)

!        if (.not. file_exists) then
!           write(LVT_logunit,*)'[WARN], missing file ',trim(hycom_fname)
!        end if

       hi_ant  = LVT_rc%udef
       call get_hycom_cice_filename('ANT', hycom_fname, &
            hi_ant_year, hi_ant_month, hi_ant_day, &
            hi_ant_hour, hi_ant_fcst_hr)

       if (trim(hycom_fname) == "NONE") then
          file_exists = .false.
       else
          file_exists = .true.
       end if

       hi_ant_1d = -9999.0
       hi_ant_ip(:) = -9999

       if(file_exists) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          write(LVT_logunit,*) '[INFO] Reading HYCOM data ', trim(hycom_fname)

          ios = nf90_open(path=trim(hycom_fname), mode=NF90_NOWRITE, ncid=nid)
          call LVT_verify(ios, 'Error opening file'//trim(hycom_fname))

!variable ids
          ios = nf90_inq_varid(nid, 'hi', hi_ant_id)
          call LVT_verify(ios, 'Error nf90_inq_varid: hi')

!values
          ios = nf90_get_var(nid, hi_ant_id, hi_ant,&
               start=(/1, 1, 1, 1/), &
               count=(/LVT_rc%HYCOM_hi_ant_nc, LVT_rc%HYCOM_hi_ant_nr, 1, 1/))
          call LVT_verify(ios, 'Error nf90_get_var: hi_ant')

          ios = nf90_close(nid)
          call LVT_verify(ios, 'Error in nf90_close')
#endif
          hi_ant_1d = -9999.0
          hi_ant_lb = .false.

          do r = 1, LVT_rc%HYCOM_hi_ant_nr
             do c = 1, LVT_rc%HYCOM_hi_ant_nc
                if (hi_ant(c,r,1,1) .ne. -30000) then
                   c1 = c
                   r1 = r
                   hi_ant_1d(c1 + (r1-1)*LVT_rc%HYCOM_hi_ant_nc) = &
                        hi_ant(c,r,1,1)*0.001

                   hi_ant_lb(c1 + (r1-1)*LVT_rc%HYCOM_hi_ant_nc) = .true.
                endif
             enddo
          enddo

          call upscaleByAveraging( &
               LVT_rc%HYCOM_hi_ant_nc*LVT_rc%HYCOM_hi_ant_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
               LVT_rc%HYCOM_hi_ant_n11, hi_ant_lb, &
               hi_ant_1d, lo, hi_ant_ip)

       end if

       ! Merge the two interpolated hi fields together.
       hi_ip(:) = -9999
       do c = 1, LVT_rc%lnc*LVT_rc%lnr
          if (hi_ant_ip(c) > -9999) then
             hi_ip(c) = hi_ant_ip(c)
          else if (hi_arc_ip(c) > -9999) then
             hi_ip(c) = hi_arc_ip(c)
          end if
       end do ! c

       ! EMK: Since sea ice is missing north of -49.5N and south of 40N, we
       ! need to set water points in this region to a reasonable value.  We
       ! assume sea ice thickness is zero in this region.
       do r = 1, LVT_rc%gnr
          do c = 1, LVT_rc%gnc
             gid = LVT_domain%gindex(c,r)
             if (gid .eq. -1) then
                if (hi_ip(c + (r-1)*LVT_rc%gnc) == -9999) then
                   hi_ip(c + (r-1)*LVT_rc%gnc) = 0
                end if
             end if
          end do ! c
       end do ! r

       ! GRIB2 Settings...updated by EMK
       gribDis   = 10
       !stepType  = "avg"
       stepType = "instant" ! EMK
       pdTemplate = 0
       gribCat   = 2
       varid_def = 1
       gribSfc   = 1
       gribSF    = 10
       gribLvl   = 1

       if (LVT_rc%lvt_out_format .eq. "grib2") then

          ! EMK...Use older hi date/time
          if (hi_ant_year .lt. hi_arc_year .or. &
               hi_ant_month .lt. hi_arc_month .or. &
               hi_ant_day .lt. hi_arc_day .or. &
               hi_ant_hour .lt. hi_arc_hour) then
             hi_arc_year = hi_ant_year
             hi_arc_month = hi_ant_month
             hi_arc_day = hi_ant_day
             hi_arc_hour = hi_arc_hour
             hi_arc_fcst_hr = hi_ant_fcst_hr
          end if

          ! add to the grib file
          call writeSingleGrib2Var(ftn_mean, &
               hi_ip, &
               varid_def, &
               gribSF, &
               gribSfc, &
               gribLvl, &
               gribDis, &
               gribCat, &
               pdTemplate, &
               stepType, &
               time_unit, &
               time_past, &
               time_curr, &
               timeRange, &
               1, &
               toplev(1), &
               botlev(1), &
               depscale(1), &
               typeOfGeneratingProcess=2, &
               typeOfProcessedData=1, &
               ref_year=hi_arc_year, &
               ref_month=hi_arc_month, &
               ref_day=hi_arc_day, &
               ref_hour=hi_arc_hour, &
               ref_fcst_hr=hi_arc_fcst_hr)

       elseif (LVT_rc%lvt_out_format .eq. "grib1") then
          call writeSingleGrib1Var(ftn_mean, &
               hi_ip, &
               varid_def, &
               gribSF, &
               gribSfc, &
               gribLvl, &
               stepType, &
               time_unit, &
               time_past, &
               time_curr, &
               timeRange, &
               1, &
               toplev(1), &
               botlev(1))
       elseif (LVT_rc%lvt_out_format .eq. "netcdf") then
          call writeSingleNetcdfVar(ftn_mean, &
               hi_ip, &
               LVT_histData%hi%varId_def, &
               1)

       endif

    endif

  end subroutine LVT_append_HYCOM_fields


  subroutine applyNoiseReductionFilter(gvar)

    real   :: gvar(LVT_rc%lnc*LVT_rc%lnr)
    real   :: gtmp(LVT_rc%lnc*LVT_rc%lnr)

    integer :: c,r, c1,r1,c_s, c_e, r_s, r_e
    real    :: avg_val
    real    :: navg_val
    real    :: sigma,wt

    gtmp = LVT_rc%udef

    if (LVT_rc%smoothingFilterType .eq. "box filter") then
       do r = 1, LVT_rc%lnr
          do c = 1, LVT_rc%lnc

             c_s = max(1, c-2)
             c_e = min(LVT_rc%lnc, c + 2)
             r_s = max(1, r - 2)
             r_e = min(LVT_rc%lnr, r + 2)

             avg_val = 0.0
             navg_val = 0
             do c1 = c_s, c_e
                do r1 = r_s,r_e
                   if (gvar(c1 + (r1-1)*LVT_rc%lnc) .ne. LVT_rc%udef) then
                      avg_val = avg_val + gvar(c1 + (r1-1)*LVT_rc%lnc)
                      navg_val = navg_val + 1
                   endif
                enddo
             enddo
             if (navg_val .gt. 0) then
                avg_val = avg_val / navg_val
             else
                avg_val = LVT_rc%udef
             endif

             gtmp(c + (r-1)*LVT_rc%lnc) = avg_val

          enddo
       enddo

    elseif (LVT_rc%smoothingFilterType .eq. "gaussian filter") then
       sigma = 1.0
       do r = 1, LVT_rc%lnr
          do c = 1, LVT_rc%lnc

             c_s = max(1, c - 2)
             c_e = min(LVT_rc%lnc, c + 2)
             r_s = max(1, r - 2)
             r_e = min(LVT_rc%lnr, r + 2)

             avg_val = 0.0
             navg_val = 0
             do c1 = c_s, c_e
                do r1 = r_s,r_e
                   if(gvar(c1 + (r1-1)*LVT_rc%lnc) .ne. LVT_rc%udef) then
                      wt = exp(-((c1-c)**2+(r1-r)**2)/(2*sigma**2))/&
                           (2*LVT_CONST_PI*sigma**2)
                      avg_val = avg_val + wt*gvar(c1 + (r1-1)*LVT_rc%lnc)
                      navg_val = navg_val + wt
                   endif
                enddo
             enddo
             if (navg_val .gt. 0) then
                if (gvar(c + (r-1)*LVT_rc%lnc) .ne. LVT_rc%udef) then
                   avg_val = avg_val / navg_val
                else
                   avg_val = LVT_rc%udef
                endif
             else
                avg_val = LVT_rc%udef
             endif

             gtmp(c + (r-1)*LVT_rc%lnc) = avg_val

          enddo
       enddo

    endif

    gvar = gtmp

  end subroutine applyNoiseReductionFilter
!BOP
!
! !ROUTINE: writeSingleGrib1Var
! \label{writeSingleGrib1Var}
!
! !INTERFACE:
  subroutine writeSingleGrib1Var(ftn, gtmp, gribId, gribSF, gribSfc, gribLvl, &
       sType, time_unit, time_p1, time_p2, &
       timeRange, k, toplev, botlev)
!
! !DESCRIPTION:
!  This subroutine writes a single variable to a grib file
!
!EOP

    integer                       :: ftn
    real                          :: gtmp(LVT_rc%lnc*LVT_rc%lnr)
    integer, intent(in)           :: gribid
    integer, intent(in)           :: gribSF
    integer, intent(in)           :: gribSfc
    integer, intent(in)           :: gribLvl
    character(len=*), intent(in)  :: sType
    integer, intent(in)           :: timeRange
    integer, intent(in)           :: time_unit
    integer, intent(in)           :: time_p1
    integer, intent(in)           :: time_p2
    integer, intent(in)           :: k
    real                          :: toplev(1)
    real                          :: botlev(1)


    character*8                   :: date
    integer                       :: yr1, mo1,da1,hr1,mn1
    integer                       :: idate,idate1
    real                          :: lat_ur, lon_ur
    real                          :: lat_ll, lon_ll
    integer                       :: igrib,iret
    integer                       :: decimalPrecision,gribSFtemp
    character*100                 :: message(20)

    ! Note passing string of defined points only to output
    ! because bitmap in GRIB-1 file will fill in the rest

#if (defined USE_ECCODES)
    call grib_new_from_samples(igrib, "GRIB1", iret)
    call LVT_verify(iret, 'grib_new_from_samples failed in LVT_DataStreamsMod')
#else
    call grib_new_from_template(igrib, "GRIB1", iret)
    call LVT_verify(iret, &
         'grib_new_from_template failed in LVT_DataStreamsMod')
#endif

    call grib_set(igrib, 'table2Version', LVT_rc%grib_table, iret)
    call LVT_verify(iret, &
         'grib_set:table2version failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'generatingProcessIdentifier', &
         LVT_rc%grib_process_id, iret)
    call LVT_verify(iret, &
         'grib_set:generatingProcessIdentifier failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'gridDefinition', LVT_rc%grib_grid_id, iret)
    call LVT_verify(iret, 'grib_set:grid ID failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'indicatorOfParameter', gribid, iret)
    call LVT_verify(iret, &
         'grib_set:indicatorOfParameter failed in LVT_DataStreamsMod')

    !    call grib_set(igrib,'paramId',gribid, iret)
    !    call LVT_verify(iret,'grib_set:paramId failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'indicatorOfTypeOfLevel', gribSfc, iret)
    call LVT_verify(iret, &
         'grib_set:indicatorOfTypeOfLevel failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'level', gribLvl, iret)
    call LVT_verify(iret, 'grib_set:level failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'topLevel', toplev(1), iret)
    call LVT_verify(iret, 'grib_set:topLevel failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'bottomLevel', botlev(1), iret)
    call LVT_verify(iret, 'grib_set:bottomLevel failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'stepType', sType, iret)
    call LVT_verify(iret, 'grib_set:stepType failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'stepUnits', time_unit, iret)
    call LVT_verify(iret, 'grib_set:stepUnits failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'startStep', time_p1, iret)
    call LVT_verify(iret, 'grib_set:startStep failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'endStep', time_p2, iret)
    call LVT_verify(iret, 'grib_set:endStep failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'timeRangeIndicator', timeRange, iret)
    call LVT_verify(iret, &
         'grib_set:timeRangeIndicator failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'swapScanningLat', 1, iret)
    call LVT_verify(iret, &
         'grib_set:swapScanningLat failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'Ni', LVT_rc%gnc, iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')

    call grib_set(igrib,'Nj', LVT_rc%gnr, iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')

    call ij_to_latlon(LVT_domain%lvtproj, float(LVT_rc%gnc),&
         float(LVT_rc%gnr), lat_ur, lon_ur)
    call ij_to_latlon(LVT_domain%lvtproj, 1.0, 1.0, &
         lat_ll, lon_ll)

    call grib_set(igrib, 'latitudeOfFirstGridPointInDegrees', lat_ll, iret)
    call LVT_verify(iret, &
         'grib_set:latitudeOfFirstGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    call grib_set(igrib, 'longitudeOfFirstGridPointInDegrees', lon_ll, iret)
    call LVT_verify(iret, &
         'grib_set:longitudeOfFirstGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    call grib_set(igrib, 'latitudeOfLastGridPointInDegrees', lat_ur, iret)
    call LVT_verify(iret, &
         'grib_set:latitudeOfLastGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    call grib_set(igrib, 'longitudeOfLastGridPointInDegrees', lon_ur, iret)
    call LVT_verify(iret, &
         'grib_set:longitudeOfLastGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    call grib_set(igrib, 'missingValue', LVT_rc%udef, iret)
    call LVT_verify(iret, &
         'grib_set:missingValue failed in LVT_DataStreamsMod')

! Should not need to fix the "num bits" value for each parameter
! if the "decimalPrecision" (aka, "DecScale") is set properly. - dmm
!     call grib_set(igrib, 'bitsPerValue',12,iret)
!     call LVT_verify(iret, 'grib_set:bitsPerValue failed in LVT_DataStreamsMod')

! Set the "decimalPrecision" (aka, "DecScale") based on the
! gribSF (grib scale factor) set in the MODEL OUTPUT TBL. - dmm
    gribSFtemp = gribSF
    decimalPrecision = 0
    do while (gribSFtemp.ge.10)
       decimalPrecision = decimalPrecision + 1
       gribSFtemp = gribSFtemp / 10
    enddo
    call grib_set(igrib, 'decimalPrecision', decimalPrecision, iret)
    call LVT_verify(iret, &
         'grib_set:decimalPrecision failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'bitmapPresent', 1, iret)
    call LVT_verify(iret, &
         'grib_set:bitmapPresent failed in LVT_DataStreamsMod')

    if (LVT_rc%domain .eq. "latlon") then
       call grib_set(igrib, 'gridType', 'regular_ll', iret)
       call LVT_verify(iret, 'grib_set: gridType failed in LVT_DataStreamsMod')

       call grib_set(igrib, 'iDirectionIncrementInDegrees', &
            LVT_rc%gridDesc(9), iret)
       call LVT_verify(iret, &
            'grib_set:iDirectionIncrementInDegrees failed in '// &
            'LVT_DataStreamsMod')

       call grib_set(igrib, 'jDirectionIncrementInDegrees', &
            LVT_rc%gridDesc(10), iret)
       call LVT_verify(iret, &
            'grib_set:jDirectionIncrementInDegrees failed in '// &
            'LVT_DataStreamsMod')

    else  !Unsupported Map Projection for GRIB output

       message(1) = 'program:  LVT_DataStreamsMod'
       message(2) = ' subroutine:  writevar_grib1_withstats_real'
       message(3) = '  Unsupported map projection for GRIB1 output!'
       call lvt_abort(message)
       stop

    endif

    da1 = LVT_rc%da
    mo1 = LVT_rc%mo
    yr1 = LVT_rc%yr

    write(unit=date, fmt='(i4.4,i2.2,i2.2)') yr1, mo1, da1
    read(date,'(I8)') idate

    call grib_set(igrib, 'dataDate', idate, iret)
    call LVT_verify(iret, 'grib_set:dataDate failed in LVT_DataStreamsMod')

    hr1 = LVT_rc%hr
    mn1 = LVT_rc%mn

    write(unit=date,fmt='(i2.2,i2.2)') hr1, mn1
    read(date,'(I4)') idate1

    call grib_set(igrib, 'dataTime', idate1, iret)
    call LVT_verify(iret, 'grib_set:dataTime failed in LVT_DataStreamsMod')


    call grib_set(igrib, 'values', gtmp, iret)
    call LVT_verify(iret, 'grib_set:values failed in LVT_DataStreamsMod')

    ! Move setting of centre and subCentre to the end of the settings.
    ! The order these are written is important and will affect output. - dmm
    call grib_set(igrib, 'centre', LVT_rc%grib_center_id, iret)
    call LVT_verify(iret, 'grib_set:centre failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'subCentre', LVT_rc%grib_subcenter_id, iret)
    call LVT_verify(iret, 'grib_set:subCentre failed in LVT_DataStreamsMod')

    call grib_write(igrib, ftn, iret)
    call LVT_verify(iret, 'grib_write failed in LVT_DataStreamsMod')

    call grib_release(igrib, iret)
    call LVT_verify(iret, 'grib_release failed in LVT_DataStreamsMod')


  end subroutine writeSingleGrib1Var

!BOP
!
! !ROUTINE: writeSingleGrib2Var
! \label{writeSingleGrib2Var}
!
! !INTERFACE:
  subroutine writeSingleGrib2Var(ftn, gtmp, gribId, gribSF, gribSfc, gribLvl,&
       gribDis, gribCat, pdTemplate, &
       sType, time_unit, time_p1, time_p2, &
       timeRange, k, toplev, botlev, depscale, &
       ensembleSpread, typeOfGeneratingProcess, &
       typeOfProcessedData, &
       ref_year, ref_month, ref_day, ref_hour, ref_fcst_hr)
!
! !DESCRIPTION:
!  This subroutine writes a single variable to a grib2 file based on
!  the implementation by Hiroko Kato within LIS.
!
!
!EOP

    integer                       :: ftn
    real                          :: gtmp(LVT_rc%lnc*LVT_rc%lnr)
    integer, intent(in)           :: gribid
    integer, intent(in)           :: gribSF
    integer, intent(in)           :: gribSfc
    integer, intent(in)           :: gribLvl
    integer, intent(in)           :: gribDis
    integer, intent(in)           :: gribCat
    integer, intent(in)           :: pdTemplate
    character(len=*), intent(in)  :: sType
    integer, intent(in)           :: timeRange
    integer, intent(in)           :: time_unit
    integer, intent(in)           :: time_p1
    integer, intent(in)           :: time_p2
    integer, intent(in)           :: k
    real                          :: toplev(1)
    real                          :: botlev(1)
    real                          :: depscale(1)
    logical, optional, intent(in) :: ensembleSpread
    integer, optional, intent(in) :: typeOfGeneratingProcess
    integer, optional, intent(in) :: typeOfProcessedData
    integer, optional, intent(in) :: ref_year
    integer, optional, intent(in) :: ref_month
    integer, optional, intent(in) :: ref_day
    integer, optional, intent(in) :: ref_hour
    integer, optional, intent(in) :: ref_fcst_hr

    integer                       :: sType_int
    character*8                   :: date
    integer                       :: yr1, mo1,da1,hr1,mn1
    integer                       :: idate,idate1
    real                          :: lat_ur, lon_ur
    real                          :: lat_ll, lon_ll
    integer                       :: igrib,iret
    integer                       :: decimalPrecision,gribSFtemp
    character*100                 :: message(20)
    logical :: ensembleSpread_local
    integer :: typeOfGeneratingProcess_local
    integer :: typeOfProcessedData_local
    type(ESMF_Time) :: time1, time2
    type(ESMF_TimeInterval) :: timeinterval,timeinterval12

    ! EMK...Handle optional ensemble spread flag
    ensembleSpread_local = .false.
    if (present(ensembleSpread)) then
       ensembleSpread_local = ensembleSpread
    end if

    ! EMK...Handle optional typeOfGeneratingProcess
    typeOfGeneratingProcess_local = 0 ! Analysis
    if (present(typeOfGeneratingProcess)) then
       typeOfGeneratingProcess_local = typeOfGeneratingProcess
    end if

    ! EMK...Handle optional typeOfProcessedData
    typeOfProcessedData_local = 0
    if (present(typeOfProcessedData)) then
       typeOfProcessedData_local = typeOfProcessedData
    end if

    ! Note passing string of defined points only to output
    ! because bitmap in GRIB-1 file will fill in the rest

#if (defined USE_ECCODES)
    call grib_new_from_samples(igrib, "GRIB2", iret)
    call LVT_verify(iret, 'grib_new_from_samples failed in LVT_DataStreamsMod')
#else
    call grib_new_from_template(igrib, "GRIB2", iret)
    call LVT_verify(iret, 'grib_new_from_template failed in LVT_DataStreamsMod')
#endif

    ! Section 0: Indicator
    ! Octet 7
    call grib_set(igrib, 'discipline', gribDis, iret)
    call LVT_verify(iret, 'grib_set: discipline failed in LVT_DataStreamsMod')

    ! Section 1: Identification
    ! Octets 6-7
    call grib_set(igrib, 'centre', LVT_rc%grib_center_id, iret)
    call LVT_verify(iret, 'grib_set:centre failed in LVT_DataStreamsMod')
    ! Octets 8-9
    call grib_set(igrib, 'subCentre', LVT_rc%grib_subcenter_id, iret)
    call LVT_verify(iret, 'grib_set:subCentre failed in LVT_DataStreamsMod')
    ! Octet 10
    call grib_set(igrib,'tablesVersion', LVT_rc%grib_table, iret)
    call LVT_verify(iret, &
         'grib_set:tablesversion failed in LVT_DataStreamsMod')
    ! Octet 11
    call grib_set(igrib, 'localTablesVersion', 1, iret)
    call LVT_verify(iret, &
         'grib_set:localTablesVersion failed in LVT_DataStreamsMod')
    ! Octet 12
    ! EMK 8 May 2018...Reference time will always be start of forecast.
    ! Since this is not available in the LIS history file, we will use
    ! the start day/time specified in the lvt.config file.
    ! Exception is for GOFS analyses.
    if (typeOfGeneratingProcess_local == 0) then ! Analysis
       call grib_set(igrib,'significanceOfReferenceTime', 0, iret)
       call LVT_verify(iret, &
            'grib_set:significanceOfReferenceTime failed in '// &
            'LVT_DataStreamsMod')
    else
       call grib_set(igrib, 'significanceOfReferenceTime', 1, iret)
       call LVT_verify(iret, &
            'grib_set:significanceOfReferenceTime failed in '// &
            'LVT_DataStreamsMod')
    end if

    if (present(ref_year) .and. present(ref_month) .and. present(ref_day) &
         .and. present(ref_hour)) then
       yr1 = ref_year
       mo1 = ref_month
       da1 = ref_day
       hr1 = ref_hour
       mn1 = 0
    else
       yr1 = LVT_rc%syr
       mo1 = LVT_rc%smo
       da1 = LVT_rc%sda
       hr1 = LVT_rc%shr
       mn1 = LVT_rc%smn
    end if

    ! Octets 13-16
    write(unit=date, fmt='(i4.4,i2.2,i2.2)') yr1, mo1, da1
    read(date,'(I8)') idate
    call grib_set(igrib, 'dataDate', idate, iret)
    call LVT_verify(iret, 'grib_set:dataDate failed in LVT_DataStreamsMod')

    ! Octets 17-19
    write(unit=date, fmt='(i2.2,i2.2)') hr1, mn1
    read(date,'(I4)') idate1
    call grib_set(igrib, 'dataTime', idate1, iret)
    call LVT_verify(iret, 'grib_set:dataTime failed in LVT_DataStreamsMod')

    ! Octet 20...Hardcode to operations for 557
    ! FIXME...Set in lvt.config
    call grib_set(igrib, 'productionStatusOfProcessedData', 0, iret)
    call LVT_verify(iret, &
         'grib_set:productionStatusOfProcessedData failed in '// &
         'LVT_DataStreamsMod')
    ! Octet 21
    call grib_set(igrib,'typeOfProcessedData', typeOfProcessedData_local, iret)
    call LVT_verify(iret, &
         'grib_set:typeOfProcessedData failed in LVT_DataStreamsMod')

!    ! ????
!    call grib_set(igrib,'stepType',sType, iret)
!    call LVT_verify(iret,'grib_set:stepType failed in LVT_DataStreamsMod')


    ! Section 2: Local Use Section (Optional) --none for now

    ! Section 3: Grid
    call grib_set(igrib, 'gridDefinitionTemplateNumber', &
         LVT_rc%grib_grid_id, iret)
    call LVT_verify(iret, &
         'grib_set:gridDefinitionTemplateNumber failed in LVT_DataStreamsMod')

    ! Hard-coded: shape of the Earth 0=radius = 6,367,470.0 m; 3.2.table
    call grib_set(igrib, 'shapeOfTheEarth', 0, iret)
    call LVT_verify(iret, &
         'grib_set:shapeOfTheEarth failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'swapScanningLat', 1, iret)
    call LVT_verify(iret, &
         'grib_set:swapScanningLat failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'Ni', LVT_rc%gnc, iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')

    call grib_set(igrib, 'Nj', LVT_rc%gnr,iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')

    call ij_to_latlon(LVT_domain%lvtproj, float(LVT_rc%gnc), &
         float(LVT_rc%gnr), lat_ur, lon_ur)
    call ij_to_latlon(LVT_domain%lvtproj, 1.0, 1.0, &
         lat_ll, lon_ll)

    call grib_set(igrib, 'latitudeOfFirstGridPointInDegrees', lat_ll, iret)
    call LVT_verify(iret, &
         'grib_set:latitudeOfFirstGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    call grib_set(igrib, 'longitudeOfFirstGridPointInDegrees', lon_ll, iret)
    call LVT_verify(iret, &
         'grib_set:longitudeOfFirstGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    call grib_set(igrib, 'latitudeOfLastGridPointInDegrees', lat_ur, iret)
    call LVT_verify(iret, &
         'grib_set:latitudeOfLastGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    call grib_set(igrib, 'longitudeOfLastGridPointInDegrees', lon_ur, iret)
    call LVT_verify(iret, &
         'grib_set:longitudeOfLastGridPointInDegrees failed in '// &
         'LVT_DataStreamsMod')

    if (LVT_rc%domain .eq. "latlon") then
       call grib_set(igrib, 'gridType', 'regular_ll', iret)
       call LVT_verify(iret, 'grib_set: gridType failed in LVT_DataStreamsMod')

       call grib_set(igrib, 'iDirectionIncrementInDegrees', &
            LVT_rc%gridDesc(9), iret)
       call LVT_verify(iret, &
            'grib_set:iDirectionIncrementInDegrees failed in '// &
            'LVT_DataStreamsMod')

       call grib_set(igrib, 'jDirectionIncrementInDegrees', &
            LVT_rc%gridDesc(10), iret)
       call LVT_verify(iret, &
            'grib_set:jDirectionIncrementInDegrees failed in '// &
            'LVT_DataStreamsMod')

    else  !Unsupported Map Projection for GRIB output

       message(1) = 'program:  LVT_DataStreamsMod'
       message(2) = ' subroutine:  writevar_grib1_withstats_real'
       message(3) = '  Unsupported map projection for GRIB1 output!'
       call lvt_abort(message)
       stop

    endif

    ! Section 4: Product Definition Section

    ! Octets 8-9
    call grib_set(igrib, 'productDefinitionTemplateNumber', pdTemplate, iret)
    call LVT_verify(iret, &
         'grib_set:productDefinitionTemplateNumber failed in '// &
         'LVT_DataStreamsMod')

    if (pdTemplate .ne. 0 .and. &
         pdTemplate .ne. 2 .and. &
         pdTemplate .ne. 12) then
       write(LVT_logunit,*) &
            '[ERR] Unsupported Product Definition Template ', pdTemplate
       call LVT_endrun()
    end if

    ! Common settings for Product Definition Templates 4.0, 4.2 and 4.12
    if (pdTemplate == 0 .or. pdTemplate == 2 .or. pdTemplate == 12) then
       ! Octet 10
       call grib_set(igrib, 'parameterCategory', gribCat, iret)
       call LVT_verify(iret, &
            'grib_set:parameterCategory failed in LVT_DataStreamsMod')
       ! Octet 11
       call grib_set(igrib, 'parameterNumber', gribid, iret)
       call LVT_verify(iret, &
            'grib_set:parameterNumber failed in LVT_DataStreamsMod')
       ! Octet 12
       call grib_set(igrib, 'typeOfGeneratingProcess', &
            typeOfGeneratingProcess_local, iret)
       call LVT_verify(iret, &
            'grib_set:typeOfGeneratingProcess failed in LVT_DataStreamsMod')
       ! Octet 13
       call grib_set(igrib, 'backgroundGeneratingProcessIdentifier', &
            LVT_rc%grib_process_id, iret)
       call LVT_verify(iret, &
            'grib_set:backgroundGeneratingProcessIdentifier failed in '// &
            'LVT_DataStreamsMod')
       ! Octet 14
       call grib_set(igrib, 'generatingProcessIdentifier', &
            LVT_rc%grib_process_id, iret)
       call LVT_verify(iret, &
            'grib_set:generatingProcessIdentifier failed in '// &
            'LVT_DataStreamsMod')
       ! Octet 15-17 are skipped

       ! Octet 18...Use hours
       call grib_set(igrib, 'indicatorOfUnitOfTimeRange', 1, iret)
       call LVT_verify(iret, &
            'grib_set:indicatorOfUnitOfTimeRange failed in LVT_DataStreamsMod')

       ! Octets 19-22...Forecast time is in hours.  Must calculate.
       ! For analyses, forecast time is always zero.
       ! In the case of PDT 4.12, the forecast time is for the start of
       ! the processing interval.
       ! If reference time is explicitly passed to routine (like for GOFS
       ! data), use that.
       if (typeOfGeneratingProcess_local == 0) then ! Analysis
          call ESMF_TimeIntervalSet(timeinterval, s=0, rc=iret)

       else if (present(ref_year) .and. present(ref_month) .and. &
            present(ref_day) .and. present(ref_hour) .and. &
            present(ref_fcst_hr)) then
          call ESMF_TimeIntervalSet(timeinterval, h=ref_fcst_hr, rc=iret)

       else
          call ESMF_TimeSet(time1, yy=LVT_rc%syr, mm=LVT_rc%smo, &
               dd=LVT_rc%sda, &
               h=LVT_rc%shr, m=LVT_rc%smn, s=LVT_rc%sss, rc=iret)
          call LVT_verify(iret,&
               'ESMF_TimeSet:time1 failed in LVT_DataStreamsMod')
          call ESMF_TimeSet(time2, yy=LVT_rc%yr, mm=LVT_rc%mo, dd=LVT_rc%da, &
               h=LVT_rc%hr, m=LVT_rc%mn, s=LVT_rc%ss, rc=iret)
          call LVT_verify(iret, &
               'ESMF_TimeSet:time2 failed in LVT_DataStreamsMod')
          if (pdTemplate == 12) then
             call ESMF_TimeIntervalSet(timeinterval12, &
                  s=LVT_rc%statswriteint, rc=iret)
             call LVT_verify(iret,&
                  'ESMF_TimeIntervalSet:timeinterval12 failed in '// &
                  'LVT_DataStreamsMod')
             timeinterval = time2 - time1 - timeinterval12
          else
             timeinterval = time2 - time1
          end if
       end if
       call ESMF_TimeIntervalGet(timeinterval, h=hr1, rc=iret)
       call LVT_verify(iret, &
            'ESMF_TimeIntervalGet:timeinterval failed in LVT_DataStreamsMod')
       call grib_set(igrib, 'forecastTime', hr1, iret)
       call LVT_verify(iret, &
            'grib_set:forecast_time failed in LVT_DataStreamsMod')

       ! Octets 23-34.  Varies by type of level/layer.
       call grib_set(igrib, 'typeOfFirstFixedSurface', gribSfc, iret)
       call LVT_verify(iret, &
            'grib_set:typeOfFirstFixedSurface failed in LVT_DataStreamsMod')

       if ( gribSfc .eq. 106 ) then   ! soil layers

          ! Must set this before scale factor/values of surfaces
          call grib_set(igrib, 'typeOfSecondFixedSurface', gribSfc, iret)
          call LVT_verify(iret, &
               'grib_set:typeOfSecondFixedSurface failed in '// &
               'LVT_DataStreamsMod')

          call grib_set(igrib, 'scaleFactorOfFirstFixedSurface', &
               depscale(1), &
               iret)
          call LVT_verify(iret, &
               'grib_set:scaleFactorOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')
          call grib_set(igrib, 'scaledValueOfFirstFixedSurface', &
               toplev(1), iret)
          call LVT_verify(iret, &
               'grib_set:scaledValueOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')

          call grib_set(igrib, 'scaleFactorOfSecondFixedSurface', &
               depscale(1), &
               iret)
          call LVT_verify(iret, &
               'grib_set:scaledFactorOfSecondFixedSurface failed in '// &
               'LVT_DataStreamsMod')
          call grib_set(igrib, 'scaledValueOfSecondFixedSurface', &
               botlev(1), &
               iret)
          call LVT_verify(iret, &
               'grib_set:scaledValueOfSecondFixedSurface failed in '// &
               'LVT_DataStreamsMod')

       elseif ( gribSfc .eq. 1 ) then    ! surface
          ! Must set this before scale factor/value of surfaces
          call grib_set(igrib, 'typeOfSecondFixedSurface', 255, iret)
          call LVT_verify(iret, &
               'grib_set:typeOfFirstFixedSurface failed in LVT_DataStreamsMod')

          call grib_set(igrib, 'scaleFactorOfFirstFixedSurface', 0, iret)
          call LVT_verify(iret, &
               'grib_set:scaledFactorOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')
          call grib_set(igrib, 'scaledValueOfFirstFixedSurface', &
               toplev(1), iret)
          call LVT_verify(iret, &
               'grib_set:scaledValueOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')

          call grib_set(igrib, 'scaleFactorOfSecondFixedSurface', 255, iret)
          call LVT_verify(iret, &
               'grib_set:scaledFactorOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')
          call grib_set(igrib, 'scaledValueOfSecondFixedSurface', 255, iret)
          call LVT_verify(iret, &
               'grib_set:scaledValueOfSecondFixedSurface failed in '// &
               'LVT_DataStreamsMod')

       else if ( gribSfc .eq. 103 ) then   ! EMK...Meters AGL
          call grib_set(igrib, 'scaleFactorOfFirstFixedSurface', &
               depscale(1), iret)
          call LVT_verify(iret, &
               'grib_set:scaledFactorOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')

          call grib_set(igrib, 'level', gribLvl, iret)
          call LVT_verify(iret, 'grib_set:level failed in LVT_DataStreamsMod')
       else   ! 114 (snow level) or old 112 ??
          write(LVT_logunit,*) 'Warning: special surface type !! '//&
               'verify scale/depth for ', gribSfc

          call grib_set(igrib, 'typeOfSecondFixedSurface', gribSfc, iret)
          call LVT_verify(iret,&
               'grib_set:typeOfFirstFixedSurface failed in LVT_DataStreamsMod')

          call grib_set(igrib, 'scaleFactorOfFirstFixedSurface', 0, iret)
          call LVT_verify(iret,&
               'grib_set:scaledFactorOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')
          call grib_set(igrib, 'scaledValueOfFirstFixedSurface', &
               toplev(1), iret)
          call LVT_verify(iret, &
               'grib_set:scaledValueOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')

          call grib_set(igrib, 'scaleFactorOfSecondFixedSurface', 0, iret)
          call LVT_verify(iret, &
               'grib_set:scaledFactorOfFirstFixedSurface failed in '// &
               'LVT_DataStreamsMod')
          call grib_set(igrib, 'scaledValueOfSecondFixedSurface', &
               botlev(1), iret)
          call LVT_verify(iret, &
               'grib_set:scaledValueOfSecondFixedSurface failed in '// &
               'LVT_DataStreamsMod')
       endif

    end if

    ! Common settings for Product Definition Templates 4.2 and 4.12, but not
    ! 4.0
    if (pdTemplate == 2 .or. pdTemplate == 12) then

       ! Octet 35
       if (ensembleSpread_local) then
          call grib_set(igrib, 'derivedForecast', 4, iret)
          call LVT_verify(iret, &
               'grib_set:derivedForecast failed in LVT_DataStreamsMod')
       else
          call grib_set(igrib, 'derivedForecast', 0, iret)
          call LVT_verify(iret, &
               'grib_set:derivedForecast failed in LVT_DataStreamsMod')
       end if

       ! Octet 36.
       call grib_set(igrib, 'numberOfForecastsInEnsemble', LVT_rc%nensem, iret)
       call LVT_verify(iret, &
            'grib_set:numberOfForecastsInEnsemble failed in '// &
            'LVT_DataStreamsMod')

    end if ! PDT 4.2 or PDT 4.12

    ! Additional entries for Product Definition Template 4.12
    if (pdTemplate == 12) then
       ! Octet 37-38
       call grib_set(igrib, 'yearOfEndOfOverallTimeInterval', LVT_rc%yr, iret)
       call LVT_verify(iret, &
            'grib_set:yearOfEndOfOverallTimeInterval failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 39
       call grib_set(igrib, 'monthOfEndOfOverallTimeInterval', LVT_rc%mo, iret)
       call LVT_verify(iret, &
            'grib_set:monthOfEndOfOverallTimeInterval failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 40
       call grib_set(igrib, 'dayOfEndOfOverallTimeInterval', LVT_rc%da, iret)
       call LVT_verify(iret, &
            'grib_set:dayOfEndOfOverallTimeInterval failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 41
       call grib_set(igrib, 'hourOfEndOfOverallTimeInterval', LVT_rc%hr, iret)
       call LVT_verify(iret, &
            'grib_set:hourOfEndOfOverallTimeInterval failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 42
       call grib_set(igrib, 'minuteOfEndOfOverallTimeInterval', &
            LVT_rc%mn, iret)
       call LVT_verify(iret, &
            'grib_set:minuteOfEndOfOverallTimeInterval failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 43
       call grib_set(igrib, 'secondOfEndOfOverallTimeInterval', &
            LVT_rc%ss, iret)
       call LVT_verify(iret, &
            'grib_set:secondOfEndOfOverallTimeInterval failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 49
       if (trim(sType) .eq. "avg") then
          sType_int = 0
       elseif (trim(sType) .eq. "accum") then
          sType_int = 1
       else if (trim(sType) .eq. "max") then
          sType_int = 2
       else if (trim(sType) .eq. "min") then
          sType_int = 3
       endif
       call grib_set(igrib, 'typeOfStatisticalProcessing', sType_int, iret)
       call LVT_verify(iret, &
            'grib_set:typeOfStatisticalProcessing failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 50
       ! Use 2 -- Successive times processed have same start time of
       ! forecast, forecast time is incremented.
       call grib_set(igrib, 'typeOfTimeIncrement', 2, iret)
       call LVT_verify(iret, &
            'grib_set:typeOfTimeIncrement failed in LVT_DataStreamsMod')

       ! Octet 51...Use hours
       call grib_set(igrib, 'indicatorOfUnitForTimeRange', 1, iret) ! Hour
       call LVT_verify(iret, &
            'grib_set:indicatorOfUnitForTimeRange failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 52-55...Time range for statistical processing
       call ESMF_TimeIntervalGet(timeinterval12, h=hr1, rc=iret)
       call LVT_verify(iret,&
            'ESMF_TimeIntervalGet:timeinterval12 failed in LVT_DataStreamsMod')
       call grib_set(igrib, 'lengthOfTimeRange', hr1, iret)
       call LVT_verify(iret,&
            'grib_set:lengthOfTimeRange failed in LVT_DataStreamsMod')

       ! Octet 56...Use minutes
       call grib_set(igrib, 'indicatorOfUnitForTimeIncrement', 0, iret)
       call LVT_verify(iret, &
            'grib_set:indicatorOfUnitForTimeIncrement failed in '// &
            'LVT_DataStreamsMod')

       ! Octet 57-60...Time increment.  This should be the LIS time step in
       ! minutes
       call grib_set(igrib, 'timeIncrement', LVT_rc%lis_ts/60, iret)
       call LVT_verify(iret, &
            'grib_set:timeIncrement failed in LVT_DataStreamsMod')

    end if ! PDT 4.12


    ! Section 5: Data Representation
    call grib_set(igrib, 'packingType', LVT_rc%grib_packing_type, iret)
    call LVT_verify(iret, 'grib_set:packingType failed in LVT_DataStreamsMod')
    call grib_set(igrib, 'missingValue', LVT_rc%udef, iret)
    call LVT_verify(iret, 'grib_set:missingValue failed in LVT_DataStreamsMod')

    ! Should not need to fix the "num bits" value for each parameter
    ! if the "decimalPrecision" (aka, "DecScale") is set properly. - dmm
    !     call grib_set(igrib, 'bitsPerValue',12,iret)
    !     call LVT_verify(iret, 'grib_set:bitsPerValue failed in LVT_DataStreamsMod')

    ! Set the "decimalPrecision" (aka, "DecScale") based on the
    ! gribSF (grib scale factor) set in the MODEL OUTPUT TBL. - dmm
     gribSFtemp = gribSF
     decimalPrecision = 0
     do while (gribSFtemp .ge. 10)
        decimalPrecision = decimalPrecision + 1
        gribSFtemp = gribSFtemp / 10
     enddo
     call grib_set(igrib, 'decimalPrecision', decimalPrecision, iret)
     call LVT_verify(iret, &
          'grib_set:decimalPrecision failed in LVT_DataStreamsMod')

     ! Section 6: Bit-Map
     call grib_set(igrib, 'bitmapPresent', 1, iret)
     call LVT_verify(iret, &
          'grib_set:bitmapPresent failed in LVT_DataStreamsMod')

     call grib_set(igrib, 'values', gtmp, iret)
     call LVT_verify(iret, 'grib_set:values failed in LVT_DataStreamsMod')

     call grib_write(igrib, ftn, iret)
     call LVT_verify(iret, 'grib_write failed in LVT_DataStreamsMod')

     call grib_release(igrib, iret)
     call LVT_verify(iret,'grib_release failed in LVT_DataStreamsMod')

  end subroutine writeSingleGrib2Var

!BOP
! !ROUTINE: defineNETCDFheaderVar
! \label{defineNETCDFheaderVar}
!
! !INTERFACE:
  subroutine defineNETCDFheaderVar(ftn, dimID, dataEntry)

! !USES:

! !ARGUMENTS:
    integer                                 :: ftn
    type(LVT_lismetadataEntry), pointer     :: dataEntry
    integer                                 :: dimID(4)
!
! !DESCRIPTION:
!    This routine writes the required NETCDF header for a single variable
!
!   The arguments are:
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be
    !    written
!   \end{description}
!
!   The routines invoked are:
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected.
!   \item[LIS\_verify](\ref{LVT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    logical       :: nmodel_status
    integer       :: data_index
    integer       :: shuffle, deflate, deflate_level
    integer       :: fill_value
    character(len=100) :: short_name

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    data_index = dataEntry%index
    fill_value = LVT_rc%udef

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level = NETCDF_deflate_level

    if (dataEntry%selectOpt .eq. 1 )then
       if (dataEntry%vlevels .gt. 1) then
          call LVT_verify(nf90_def_dim(ftn, &
               trim(dataEntry%short_name)//'_profiles', &
               dataEntry%vlevels, dimID(3)), &
               'nf90_def_dim failed (2d gridspace) in LVT_DataStreamsMod')
       endif

       !EMK...Added suffix to clarify if field is instantaneous, time averaged,
       !or an accumulation.
       if (index(trim(dataEntry%short_name),"_max") .gt. 0) then
          short_name = trim(dataEntry%short_name)
       else if (index(trim(dataEntry%short_name),"_min") .gt. 0) then
          short_name = trim(dataEntry%short_name)
       else if (dataEntry%timeAvgOpt.eq.0) then
          short_name = trim(dataEntry%short_name)//"_inst"
       else if (dataEntry%timeAvgOpt.eq.1 .or. &
            dataEntry%timeAvgOpt.eq.2) then
          short_name = trim(dataEntry%short_name)//"_tavg"
       else if (dataEntry%timeAvgOpt.eq.3) then
          short_name = trim(dataEntry%short_name)//"_acc"
       else
          write(LVT_logunit,*) '[ERR] Cannot handle ', &
               trim(dataEntry%short_name)
          call LVT_endrun()
       end if

       if (dataEntry%vlevels .gt. 1) then
          call LVT_verify(nf90_def_var(ftn,trim(short_name), &
               nf90_float, &
               dimids = dimID(1:3), varID=dataEntry%varId_def), &
               'nf90_def_var for '//trim(short_name)// &
               'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn, &
               dataEntry%varId_def, &
               1,fill_value), 'nf90_def_var_fill failed for '// &
               dataEntry%short_name)

          call LVT_verify(nf90_def_var_deflate(ftn, &
               dataEntry%varId_def, &
               shuffle, deflate, deflate_level), &
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)// &
               'failed in defineNETCDFheadervar')
#endif
       else
          call LVT_verify(nf90_def_var(ftn,trim(short_name), &
               nf90_float, &
               dimids = dimID(1:2), varID=dataEntry%varId_def), &
               'nf90_def_var for '//trim(short_name)// &
               'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn, &
               dataEntry%varId_def, &
               1,fill_value), 'nf90_def_var_fill failed for '// &
               dataEntry%short_name)

          call LVT_verify(nf90_def_var_deflate(ftn, &
               dataEntry%varId_def, &
               shuffle, deflate, deflate_level), &
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)// &
               'failed in defineNETCDFheadervar')
#endif
       endif


       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def, &
            "units", trim(dataEntry%units)), &
            'nf90_put_att for units failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def, &
            "standard_name", trim(dataEntry%standard_name)), &
            'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def, &
            "long_name",trim(dataEntry%long_name)), &
            'nf90_put_att for long_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def, &
            "scale_factor", 1.0), &
            'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def, &
            "add_offset", 0.0), &
            'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def, &
            "missing_value", LVT_rc%udef), &
            'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def, &
            "_FillValue", LVT_rc%udef), &
            'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def,&
            "vmin", dataEntry%valid_min),&
            'nf90_put_att for vmin failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varId_def,&
            "vmax", dataEntry%valid_max),&
            'nf90_put_att for vmax failed in defineNETCDFheaderVar')

    endif
#endif
  end subroutine defineNETCDFheaderVar


!BOP
! !ROUTINE: defineNETCDFheaderVar_SS
! \label{defineNETCDFheaderVar_SS}
!
! !INTERFACE:
  subroutine defineNETCDFheaderVar_SS(ftn, dimID, dataEntry)

! !USES:

! !ARGUMENTS:
    integer                                 :: ftn
    type(LVT_lismetadataEntry), pointer     :: dataEntry
    integer                                 :: dimID(4)
!
! !DESCRIPTION:
!    This routine writes the required NETCDF header for a single variable
!
!   The arguments are:
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be
    !    written
!   \end{description}
!
!   The routines invoked are:
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected.
!   \item[LIS\_verify](\ref{LVT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    logical       :: nmodel_status
    integer       :: data_index
    integer       :: shuffle, deflate, deflate_level
    integer       :: fill_value
    character(len=100) :: short_name

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    data_index = dataEntry%index
    fill_value = LVT_rc%udef

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level = NETCDF_deflate_level

    if (dataEntry%selectOpt .eq. 1) then
       if (dataEntry%vlevels .gt. 1) then
          call LVT_verify(nf90_def_dim(ftn, &
               trim(dataEntry%short_name)//'_profiles', &
               dataEntry%vlevels, dimID(3)), &
               'nf90_def_dim failed (2d gridspace) in LVT_DataStreamsMod')
       endif

       !EMK...Added suffix to clarify if field is instantaneous, time averaged,
       !or an accumulation.
       if (index(trim(dataEntry%short_name),"_max") .gt. 0) then
          short_name = trim(dataEntry%short_name)
       else if (index(trim(dataEntry%short_name),"_min") .gt. 0) then
          short_name = trim(dataEntry%short_name)
       else if (dataEntry%timeAvgOpt.eq.0) then
          short_name = trim(dataEntry%short_name)//"_inst"
       else if (dataEntry%timeAvgOpt.eq.1 .or. &
            dataEntry%timeAvgOpt.eq.2) then
          short_name = trim(dataEntry%short_name)//"_tavg"
       else if (dataEntry%timeAvgOpt.eq.3) then
          short_name = trim(dataEntry%short_name)//"_acc"
       else
          write(LVT_logunit,*)'[ERR] Cannot handle ',trim(dataEntry%short_name)
          call LVT_endrun()
       end if

       if(dataEntry%vlevels.gt.1) then
          call LVT_verify(nf90_def_var(ftn, trim(short_name),&
               nf90_float,&
               dimids=dimID(1:3), varID=dataEntry%varid_ss),&
               'nf90_def_var for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn, &
               dataEntry%varid_ss, &
               1,fill_value), 'nf90_def_var_fill failed for '// &
               dataEntry%short_name)

          call LVT_verify(nf90_def_var_deflate(ftn, &
               dataEntry%varid_ss, &
               shuffle, deflate, deflate_level), &
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)// &
               'failed in defineNETCDFheadervar')
#endif
       else
          call LVT_verify(nf90_def_var(ftn,trim(short_name), &
               nf90_float, &
               dimids = dimID(1:2), varID=dataEntry%varid_ss), &
               'nf90_def_var for '//trim(dataEntry%short_name)// &
               'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn, &
               dataEntry%varid_ss, &
               1,fill_value), 'nf90_def_var_fill failed for '// &
               dataEntry%short_name)

          call LVT_verify(nf90_def_var_deflate(ftn, &
               dataEntry%varid_ss, &
               shuffle, deflate, deflate_level), &
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)// &
               'failed in defineNETCDFheadervar')
#endif
       endif

       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss, &
            "units",trim(dataEntry%units)), &
            'nf90_put_att for units failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss, &
            "standard_name", trim(dataEntry%standard_name)), &
            'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss, &
            "long_name", trim(dataEntry%long_name)), &
            'nf90_put_att for long_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss, &
            "scale_factor", 1.0), &
            'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss,&
            "add_offset", 0.0),&
            'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss,&
            "missing_value", LVT_rc%udef),&
            'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss,&
            "_FillValue", LVT_rc%udef),&
            'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss,&
            "vmin", dataEntry%valid_min),&
            'nf90_put_att for vmin failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn, dataEntry%varid_ss,&
            "vmax" ,dataEntry%valid_max),&
            'nf90_put_att for vmax failed in defineNETCDFheaderVar')

    endif
#endif
  end subroutine defineNETCDFheaderVar_SS

!BOP
!
! !ROUTINE: writeSingleNetcdfVar
! \label{writeSingleNetcdfVar}
!
! !INTERFACE:
  subroutine writeSingleNetcdfVar(ftn, gtmp, varID, k)
!
! !DESCRIPTION:
!  This subroutine writes a single variable to a grib file
!
!EOP

! !ARGUMENTS:
    integer                       :: ftn
    real                          :: gtmp(LVT_rc%lnc*LVT_rc%lnr)
    integer, intent(in)           :: varid
    integer, intent(in)           :: k

    integer                       :: c,r
    real                          :: gtmp2d(LVT_rc%lnc,LVT_rc%lnr)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    do r = 1, LVT_rc%lnr
       do c = 1, LVT_rc%lnc
          gtmp2d(c,r) = gtmp(c + (r-1)*LVT_rc%lnc)
       enddo
    enddo

    call LVT_verify(nf90_put_var(ftn, varID, gtmp2d, (/1, 1, k/),&
         (/LVT_rc%gnc, LVT_rc%gnr, 1/)),&
         'nf90_put_var failed for in LVT_DataStreamsMod')
#endif

  end subroutine writeSingleNetcdfVar

!BOP
!
! !ROUTINE: LVT_tavgDataStreams
! \label{LVT_tavgDataStreams}
!
! !INTERFACE:
  subroutine LVT_tavgDataStreams
!
! !USES:
    use LVT_statsDataMod

    implicit none
!
!
! !DESCRIPTION:
!   This routine invokes the calls to compute temporal averages of
!   desired set of variables, based on the specified
!   temporal averaging frequency.
!
!   The routines invoked are:
!   \begin{description}
!    \item[tavgSingleDataStream](\ref{tavgSingleDataStream})
!     computes the temporal average for a single variable
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP

    integer      :: kk
    type(LVT_metadataEntry), pointer :: dataEntry
    type(LVT_metadataEntry), pointer :: ds1, ds2
    logical :: local_computeFlag

    ! EMK...557 post runmode has different requirements for when to set
    ! the computeFlag than the normal applications of LVT.  To prevent
    ! surprises for normal LVT users, we will use a local computeFlag
    ! variable here.
    local_computeFlag = LVT_rc%computeFlag
    if (LVT_rc%runmode.eq."557 post") then
       local_computeFlag = LVT_557post_alarm_is_on()
    end if
    !if(LVT_rc%computeFlag) then
    if (local_computeFlag) then

!data stream 1
       do kk = 1, LVT_rc%nDataStreams
          if (kk .eq. 1) then
             dataEntry => LVT_histData%head_ds1_list
          elseif (kk .eq. 2) then
             dataEntry => LVT_histData%head_ds2_list
          elseif (kk .eq. 3) then
             dataEntry => LVT_histData%head_ds3_list
          endif

          do while (associated(dataEntry))
             call tavgSingleDataStream(dataEntry)
             dataEntry => dataEntry%next
          enddo

! copy duplicate entries
! Note that this check is not enabled for three datastrems.
! The responsibility of ensuring non-duplicate entries is
! on the user.
          if (LVT_rc%ds1_dup) then
             ds1 => LVT_histData%head_ds1_list
             do while (associated(ds1))
                ds2 => ds1%next
                do while (associated(ds2))
                   if (ds2%index .ne. ds1%index .and.&
                        ds1%short_name .eq. ds2%short_name) then
                      ds2%value = ds1%value
                      ds2%count = ds1%count
                   endif
                   ds2 => ds2%next
                enddo
                ds1 => ds1%next
             enddo

          endif

          if (LVT_rc%ds2_dup) then
             ds1 => LVT_histData%head_ds2_list
             do while (associated(ds1))
                ds2 => ds1%next
                do while (associated(ds2))
                   if(ds2%index .ne. ds1%index .and.&
                        ds1%short_name .eq. ds2%short_name) then
                      ds2%value = ds1%value
                      ds2%count = ds1%count
                   endif
                   ds2 => ds2%next
                enddo
                ds1 => ds1%next
             enddo

          endif

          if(LVT_rc%var_based_strat .gt. 0) then
             call tavgSingleDataStream(LVT_histData%strat_varEntry)
             LVT_stats%strat_var(:,:,:) = &
                  LVT_histData%strat_varEntry%value(:,:,:)
          endif
       enddo
    endif
  end subroutine LVT_tavgDataStreams

!BOP
!
! !ROUTINE: tavgSingleDataStream
! \label{tavgSingleDataStream}
!
! !INTERFACE:
  subroutine tavgSingleDataStream (dataEntry)
!
! !USES:

    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This routine temporally averages the accumulated data in a
!   given datastream
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
!BOP
! !ARGUMENTS:
    type(LVT_metadataEntry) :: dataEntry
!EOP
    integer :: k,t,c,r,m,gid

    ! EMK...Special rules for accumulations, max, and min values are now
    ! handled by incrementing count variable.
!    ! EMK Do not average the data if the raw input are accumulations.
!    ! This is indicated in the lvt.config file
!    if (dataEntry%timeAvgOpt .eq. 3) return
!    ! EMK....Do not average data if Tair_f_max or Tair_f_min
!    if (trim(dataEntry%short_name) == "Tair_f_max") return
!    if (trim(dataEntry%short_name) == "Tair_f_min") return

    if (dataEntry%selectNlevs .ge. 1) then
       if (LVT_rc%computeEnsMetrics .eq. 1) then
          do t = 1, LVT_LIS_rc(1)%ntiles
             do k = 1, dataEntry%vlevels
                c = LVT_LIS_domain(1)%tile(t)%col
                r = LVT_LIS_domain(1)%tile(t)%row
                if (LVT_LIS_domain(1)%gindex(c,r) .ne. -1) then
                   gid = LVT_LIS_domain(1)%gindex(c,r)
                   do m = 1,LVT_rc%nensem
                      if (dataEntry%count(t,m,k) .ne. 0) then
                         dataEntry%value(t,m,k) = &
                              dataEntry%value(t,m,k)/dataEntry%count(t,m,k)

                      endif
                   enddo
                endif
             enddo
          enddo
       else
          do r = 1,LVT_rc%lnr
             do c = 1,LVT_rc%lnc
                do k = 1,dataEntry%vlevels
                   if (LVT_domain%gindex(c,r) .ne. -1) then
                      gid = LVT_domain%gindex(c,r)
                      do m = 1,LVT_rc%nensem
                         if (dataEntry%count(gid,m,k).ne.0) then
                            dataEntry%value(gid,m,k) = &
                                 dataEntry%value(gid,m,k) / &
                                 dataEntry%count(gid,m,k)
                         endif
                      enddo
                   endif
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine tavgSingleDataStream


!BOP
!
! !ROUTINE: LVT_resetDataStreams
! \label{LVT_resetDataStreams}
!
! !INTERFACE:
  subroutine LVT_resetDataStreams
!
! !USES:
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This routine reinitializes the data structures that hold the observational
!   data
!
!   The routines invoked are:
!   \begin{description}
!    \item[resetSingleDataStream2](\ref{resetSingleDataStream2})
!     resets the datastructures for a single variable
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP

    type(LVT_metadataEntry), pointer :: ds1
    type(LVT_metadataEntry), pointer :: ds2
    type(LVT_metadataEntry), pointer :: ds3
    logical :: local_computeFlag

    ! EMK..."557 post" runmode has different requirements for setting the
    ! compute flag compared to normal LVT users.  So we use a local variable
    ! to accomodate both.
    local_computeFlag = LVT_rc%computeFlag
    if (LVT_rc%runmode .eq. "557 post") then
       local_computeFlag = LVT_557post_alarm_is_on()
    end if

!    if(LVT_rc%computeFlag) then
    if (local_computeFlag) then
!data stream 1
       ds1 => LVT_histData%head_ds1_list

       do while (associated(ds1))
          call resetSingleDataStream(ds1)
          ds1 => ds1%next
       enddo

!data stream 2
       ds2 => LVT_histData%head_ds2_list

       do while (associated(ds2))
          call resetSingleDataStream(ds2)
          ds2 => ds2%next
       enddo

       if (LVT_rc%nDataStreams .gt. 2) then

!data stream 3
          ds3 => LVT_histData%head_ds3_list

          do while (associated(ds3))
             call resetSingleDataStream(ds3)
             ds3 => ds3%next
          enddo
       endif
!need special handler for LIS output
       if (LVT_rc%lis_output_obs) then
          if (LVT_rc%obssource(1) .eq. "LIS output") then
             call LVT_resetLISoutputContainers(1)
          endif
          if (LVT_rc%obssource(2) .eq. "LIS output") then
             call LVT_resetLISoutputContainers(2)
          endif
          if (LVT_rc%nDataStreams .gt. 2) then
             if (LVT_rc%obssource(3) .eq. "LIS output") then
                call LVT_resetLISoutputContainers(3)
             endif
          endif
       endif
    endif
  end subroutine LVT_resetDataStreams


!BOP
!
! !ROUTINE: resetSingleDataStream
! \label{resetSingleDataStream}
!
! !INTERFACE:
  subroutine resetSingleDataStream(dataEntry)
!
! !USES:
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine resets the data structures that hold the observational
!  data and the temporal averaging counters
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
!BOP
!
! !ARGUMENTS:
    type(LVT_metadataEntry) :: dataEntry
!
!EOP

    integer                 :: k

    if (dataEntry%selectNlevs .ge. 1) then
       do k = 1, dataEntry%vlevels
          dataEntry%value(:,:,k) = 0
          dataEntry%count(:,:,k) = 0
          dataEntry%count_status(:,:,k) = 0
          if(dataEntry%stdev_flag) then
             dataEntry%count_stdev(:,k)= 0
             dataEntry%stdev(:,k) = 0
          endif
       enddo
    endif
  end subroutine resetSingleDataStream

  ! EMK...Construct filename for HYCOM SST.  Targets GOFS 93.0
  subroutine get_hycom_sst_filename(sst_filename, sst_year, sst_month, &
       sst_day, sst_hour, sst_fcst_hr)

     implicit none

     ! Arguments
     character(len=100), intent(inout) :: sst_filename
     integer, intent(out) :: sst_year
     integer, intent(out) :: sst_month
     integer, intent(out) :: sst_day
     integer, intent(out) :: sst_hour
     integer, intent(out) :: sst_fcst_hr

     ! GOFS SST fields are generated from a single 00Z cycle, with output
     ! 6-hrly from 0 to 24 hours.
     integer :: sst_julhr, lvt_julhr
     logical :: file_exists

     ! First guess for SST file
     sst_year = LVT_rc%yr
     sst_month = LVT_rc%mo
     sst_day = LVT_rc%da
     sst_hour = 0
     if (LVT_rc%hr .lt. 6) then
        sst_fcst_hr = 0
     else if (LVT_rc%hr .lt. 12) then
        sst_fcst_hr = 6
     else if (LVT_rc%hr .lt. 18) then
        sst_fcst_hr = 12
     else
        sst_fcst_hr = 18
     end if
     call construct_hycom_sst_filename(LVT_rc%HYCOMdir, &
          sst_year, sst_month, sst_day, &
          sst_hour, sst_fcst_hr, sst_filename)

     write(LVT_logunit,*)'----------------------------------------------------'
     write(LVT_logunit,*)'[INFO] *** SEARCHING FOR GOFS SST ***'
     write(LVT_logunit,*)'[INFO] Trying ', trim(sst_filename)
     inquire(file=trim(sst_filename), exist=file_exists)
     if (file_exists) then
        write(LVT_logunit,*)'[INFO] Will use ', trim(sst_filename)
        return
     end if

     ! At this point, we are rolling back to earlier SST file
     call LVT_get_julhr(LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
          0, 0, 0, lvt_julhr)
     sst_julhr = lvt_julhr

     ! Start looping for earlier files
     do
        write(LVT_logunit,*) '[WARN] Cannot find ', trim(sst_filename)
        sst_fcst_hr = sst_fcst_hr - 6
        if (sst_fcst_hr < 0) then
           sst_fcst_hr = 24
           sst_julhr = sst_julhr - 24 ! Roll back to previous 00Z cycle
           ! Give up after 5 days
           if ((lvt_julhr - sst_julhr) > 24*5) then
              write(LVT_logunit,*)'[WARN] *** GIVING UP ON GOFS SST! ***'
              write(LVT_logunit,*)'[WARN] *** NO GOFS SST AVAILABLE!!! ***'
              sst_filename = 'NONE'
              return
           end if
           call LVT_julhr_date(sst_julhr, sst_year, sst_month, sst_day, &
                sst_hour)
        end if

        call construct_hycom_sst_filename(LVT_rc%HYCOMdir, &
             sst_year, sst_month, sst_day, &
             sst_hour, sst_fcst_hr, sst_filename)
        write(LVT_logunit,*) '[INFO] Trying ', trim(sst_filename)
        inquire(file=trim(sst_filename), exist=file_exists)
        if (file_exists) then
           write(LVT_logunit,*) '[INFO] Will use ', trim(sst_filename)
           return
        end if
     end do
     return
  end subroutine get_hycom_sst_filename

  ! EMK
  subroutine construct_hycom_sst_filename(rootdir, &
       yr, mo, da, hr, fcst_hr, sst_filename)

     implicit none

     ! Arguments
     character(len=*), intent(in)  :: rootdir
     integer, intent(in) :: yr
     integer, intent(in) :: mo
     integer, intent(in) :: da
     integer, intent(in) :: hr
     integer, intent(in) :: fcst_hr
     character(len=*), intent(out) :: sst_filename

     ! Local variables
     character(len=10) :: yyyymmddhh
     character(len=4)  :: thhh

     write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yr, mo, da, hr
     write(thhh,'(a1,i3.3)') 't', fcst_hr

     sst_filename = trim(rootdir) // '/hycom_glb_sfc_u_' //yyyymmddhh// &
          '_'//thhh//'.nc'
     return
  end subroutine construct_hycom_sst_filename

  ! EMK...Construct filename for HYCOM sea ice.  Targets GOFS 93.0
  subroutine get_hycom_cice_filename(region, cice_filename, &
       cice_year, cice_month, cice_day, cice_hour, cice_fcst_hr)

     implicit none

     ! Arguments
     character(len=3), intent(in) :: region
     character(len=100), intent(inout) :: cice_filename
     integer, intent(out) :: cice_year
     integer, intent(out) :: cice_month
     integer, intent(out) :: cice_day
     integer, intent(out) :: cice_hour
     integer, intent(out) :: cice_fcst_hr

     ! GOFS CICE fields are generated from a single 12Z cycle, with output
     ! 12-hrly from 0 to 180 hours.
     integer :: cice_julhr, lvt_julhr
     logical :: file_exists

     ! First guess
     call LVT_get_julhr(LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
          12,0,0, lvt_julhr)
     if (LVT_rc%hr .ge. 12) then
        cice_julhr = lvt_julhr
        cice_fcst_hr = 0
     else
        cice_julhr = lvt_julhr - 24 ! Must be previous day
        cice_fcst_hr = 12
     end if
     call LVT_julhr_date(cice_julhr, cice_year, cice_month, cice_day, &
          cice_hour)
     call construct_hycom_cice_filename(LVT_rc%HYCOMdir, &
          region, &
          cice_year, cice_month, cice_day, &
          cice_hour, cice_fcst_hr, cice_filename)

     write(LVT_logunit,*)'----------------------------------------------------'
     write(LVT_logunit,*)'[INFO] *** SEARCHING FOR GOFS CICE FOR ',&
          trim(region), ' REGION ***'
     write(LVT_logunit,*)'[INFO] Trying ', trim(cice_filename)
     inquire(file=trim(cice_filename), exist=file_exists)
     if (file_exists) then
        write(LVT_logunit,*)'[INFO] Will use ', trim(cice_filename)
        return
     end if

     ! At this point, we are rolling back to earlier CICE file
     ! Start looping for earlier files
     do
        write(LVT_logunit,*) '[WARN] Cannot find ', trim(cice_filename)
        cice_fcst_hr = cice_fcst_hr + 24
        cice_julhr = cice_julhr - 24
        if ((lvt_julhr - cice_julhr) > 24*5) then
           write(LVT_logunit,*) &
                '[WARN] *** GIVING UP ON GOFS CICE FOR ', trim(region), ' ***'
           write(LVT_logunit,*) &
                '[WARN] *** NO GOFS CICE DATA FOR ', trim(region), &
                ' AVAILABLE!!! ***'
           cice_filename = 'NONE'
           return
        end if
        call LVT_julhr_date(cice_julhr, cice_year, cice_month, cice_day, &
             cice_hour)

        call construct_hycom_cice_filename(LVT_rc%HYCOMdir, &
             region, &
             cice_year, cice_month, cice_day, &
             cice_hour, cice_fcst_hr, cice_filename)
        write(LVT_logunit,*)'[INFO] Trying ', trim(cice_filename)
        inquire(file=trim(cice_filename), exist=file_exists)
        if (file_exists) then
           write(LVT_logunit,*)'[INFO] Will use ', trim(cice_filename)
           return
        end if
     end do

     return
  end subroutine get_hycom_cice_filename

  ! EMK
  subroutine construct_hycom_cice_filename(rootdir, &
       region, yr, mo, da, hr, fcst_hr, cice_filename)

     implicit none

     ! Arguments
     character(len=*), intent(in)  :: rootdir
     character(len=3), intent(in) ::  region
     integer, intent(in) :: yr
     integer, intent(in) :: mo
     integer, intent(in) :: da
     integer, intent(in) :: hr
     integer, intent(in) :: fcst_hr
     character(len=*), intent(out) :: cice_filename

     ! Local variables
     character(len=10) :: yyyymmddhh
     character(len=4)  :: thhh

     write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yr, mo, da, hr
     write(thhh,'(a1,i3.3)') 't', fcst_hr

     cice_filename = trim(rootdir) // '/hycom-cice_inst_' // trim(region) &
          // 'u0.08_930_' // yyyymmddhh // '_'//thhh//'.nc'
     return
  end subroutine construct_hycom_cice_filename

  ! EMK...Calculate mean and standard deviation using Welford algorithm.
  ! See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  ! Intent is to avoid catastrophic cancellation when calculating variance,
  ! which can result in negative variance and imaginary standard deviation.
  ! REFERENCES:
  ! Chan, T F, G H Golub, and R J LeVeque, 1983:  Algorithms for computing
  !   the sample variance:  Analysis and recommendations.  The American
  !   Statistician, 37, 242-247.  doi:10.1080/00031305.1983.10483115.
  ! Knuth, D E, 1998:  The Art of Computer Programming, volume 2:
  !   Seminumerical Algorithms.  Third Edition, p 232.  Boston:
  !   Addison-Wesley.
  ! Ling, R F, 1974:  Comparisons of several algorithms for computing
  !   sample means and variances.  Journal of the American Statistical
  !   Association, 69, 859-866.  doi:10.2307/2286154.
  ! Welford, B P, 1962:  Note on a method for calculating corrected sums of
  !   squares and products.  Technometrics, 4, 419-420. doi:10.2307/1266577.

  subroutine welford_update(count, mean, m2, new_value)
     implicit none
     integer, intent(inout) :: count
     real, intent(inout) :: mean
     real, intent(inout) :: m2
     real, intent(in) :: new_value
     real :: delta, delta2
     count = count + 1
     delta = new_value - mean
     mean = mean + (delta / real(count))
     delta2 = new_value - mean
     m2 = m2 + (delta * delta2)
  end subroutine welford_update
  subroutine welford_finalize(count, mean, m2, stddev)
     implicit none
     integer, intent(in) :: count
     real, intent(inout) :: mean
     real, intent(in) :: m2
     real, intent(out) :: stddev
     if (count < 1) then
        mean = 0
        stddev = 0
     else
        stddev = sqrt(m2 / real(count))
     end if
  end subroutine welford_finalize

end module LVT_DataStreamsMod
