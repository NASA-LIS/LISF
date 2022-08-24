! micromet_code.f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "LIS_misc.h"

      subroutine MICROMET_CODE(nx,ny,xmn,ymn,deltax,deltay,&
     &  iyear_init,imonth_init,iday_init,xhour_init,dt,undef,&
     &  ifill,iobsint,dn,iter,curve_len_scale,slopewt,curvewt,&
     &  topo,curvature,terrain_slope,slope_az,Tair_grid,&
     &  rh_grid,uwind_grid,vwind_grid,Qsi_grid,prec_grid,&
     &  i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,&
     &  i_prec_flag,isingle_stn_flag,igrads_metfile,&
     &  windspd_grid,winddir_grid,windspd_flag,winddir_flag,&
     &  sprec,windspd_min,Qli_grid,i_longwave_flag,vegtype,&
     &  forest_LAI,iyear,imonth,iday,xhour,corr_factor,&
     &  icorr_factor_index,lapse_rate_user_flag,&
     &  iprecip_lapse_rate_user_flag,use_shortwave_obs,&
     &  use_longwave_obs,use_sfc_pressure_obs,sfc_pressure,&
     &  run_enbal,run_snowpack,calc_subcanopy_met,vegsnowd_xy,&
     &  gap_frac,cloud_frac_factor,barnes_lg_domain,n_stns_used,&
     &  k_stn,xlat_grid,xlon_grid,UTC_flag,icorr_factor_loop,&
     &  snowmodel_line_flag,xg_line,yg_line,irun_data_assim,&
     &  wind_lapse_rate,iprecip_scheme,cf_precip_flag,cf_precip,&
     &  cloud_frac_grid,snowfall_frac,seaice_run,metforce_opt)

      use snowmodel_inc
!KRA
      use LIS_logMod,   only : LIS_logunit
      use LIS_coreMod
!KRA

      implicit none

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell
      integer nstns_orig ! number of input values

      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      double precision xstn_orig(nstns_max)     ! input stn x coords
      double precision ystn_orig(nstns_max)     ! input stn y coords
      real Tair_orig(nstns_max)     ! input values
      real rh_orig(nstns_max)       ! input values
      real winddir_orig(nstns_max)  ! input values
      real windspd_orig(nstns_max)  ! input values
      real prec_orig(nstns_max)     ! input values
      real elev_orig(nstns_max)     ! station elevation
      real dn               ! average observation spacing
      real topo(nx,ny)      ! grid topography
      real xlat_grid(nx,ny) ! lat (dec deg) of cell centers
      real xlon_grid(nx,ny) ! lon (dec deg) of cell centers

      real Tair_grid(nx,ny)   ! output values
      real rh_grid(nx,ny)     ! output values
      real uwind_grid(nx,ny)  ! output, E-W wind component
      real vwind_grid(nx,ny)  ! output, N-S wind component
      real windspd_grid(nx,ny)
      real winddir_grid(nx,ny)
      real Qsi_grid(nx,ny)    ! output
      real Qli_grid(nx,ny)    ! output
      real prec_grid(nx,ny)   ! output
      real sprec(nx,ny)   ! output
      real sfc_pressure(nx,ny)

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour
      real dt                    ! model time step, in seconds
      integer iter               ! model iteration
      integer iyear_init     ! model start year
      integer imonth_init    ! model start month
      integer iday_init      ! model start day
      real xhour_init        ! model start hour
      integer J_day          ! model Julian day, actually day-of-year

      real undef       ! undefined value
      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      real curvature(nx,ny)     ! topographic curvature
      real slope_az(nx,ny)      ! azimuth of topographic slope
      real terrain_slope(nx,ny) ! terrain slope
      real vegtype(nx,ny)
      real vegsnowd_xy(nx,ny)
      real, save, allocatable :: topo_ref_grid(:,:)! reference surface

      real curve_len_scale   ! length scale for curvature calculation
      real slopewt           ! wind model slope weight
      real curvewt           ! wind model curvature weight

      integer i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,&
     &  i_prec_flag,i_longwave_flag,isingle_stn_flag,igrads_metfile,&
     &  lapse_rate_user_flag,iprecip_lapse_rate_user_flag,n_stns_used,&
     &  icorr_factor_loop,irun_data_assim,iprecip_scheme,metforce_opt

      real windspd_flag,winddir_flag,windspd_min,calc_subcanopy_met,&
     &  T_lapse_rate,Td_lapse_rate,precip_lapse_rate,&
     &  use_shortwave_obs,use_longwave_obs,use_sfc_pressure_obs,&
     &  run_enbal,run_snowpack,gap_frac,cloud_frac_factor,&
     &  barnes_lg_domain,UTC_flag,wind_lapse_rate

      integer, parameter :: nftypes = 5
      real forest_LAI(nftypes)

      real corr_factor(nx_max,ny_max,max_obs_dates+1)
      integer icorr_factor_index(max_time_steps)
      integer k_stn(nx,ny,9)

      real cf_precip(nx,ny)
      real cf_precip_flag

      real cloud_frac_grid(nx,ny)
      real snowfall_frac

      real seaice_run
      integer i,j,irec,irec_day

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(.not.allocated(topo_ref_grid)) then
        allocate(topo_ref_grid(nx,ny))
        topo_ref_grid = 0.0
      end if

! Calculate what the current simulation date should be.
      call get_model_time(iyear_init,imonth_init,iday_init,&
     &  xhour_init,iter,dt,iyear,imonth,iday,xhour,J_day)

      if (irun_data_assim.eq.1) then
        write (*,150) icorr_factor_loop,iyear,imonth,iday,xhour
  150   format('In Assim Loop #',i1,';  WORKING ON MODEL TIME =', &
     &    i5,2i4,f6.1)
      else
        write (*,151) iyear,imonth,iday,xhour
  151   format('                   WORKING ON MODEL TIME =',&
     &    i5,2i4,f6.1)
      endif

! Read in the observations for this time step, and build an array of
!   valid observations to be interpolated.
      call get_obs_data(nstns_orig,Tair_orig,rh_orig,xstn_orig,&
     &  ystn_orig,elev_orig,iyear,imonth,iday,xhour,undef,&
     &  windspd_orig,winddir_orig,prec_orig,isingle_stn_flag,&
     &  igrads_metfile,iter)

! Make the topographic calculations required by the wind and solar
!   radiation models.  These calculations are not fixed in time
!   because as the snow depth evolves it modifies the "topography".
      call topo_data(nx,ny,deltax,deltay,topo,&
     &  curvature,terrain_slope,slope_az,curve_len_scale)

! Calculate the temperature and dew-point lapse rates to be used in
!   the interpolations.
      call get_lapse_rates(imonth,iday,T_lapse_rate,&
     &  Td_lapse_rate,xlat_grid(1,1),lapse_rate_user_flag,&
     &  precip_lapse_rate,iprecip_lapse_rate_user_flag)

! Calculate the forest lai for each of the five forest types, and
!   for this day of the simulation (the lai varies seasonally for
!   the case of deciduous trees).
      call get_lai(J_day,forest_LAI)

! If this is a Lagrangian sea ice parcel trajectory simulation,
!   extract the lat-lon of the parcels at this time step.  Note
!   that these files have daily data in them, so the record has
!   to be adjusted to account for this.
      if (seaice_run.eq.4.0) then 

        if (iter.eq.1) then 
          open (91,file='extra_met/grid_lat_time.gdat',&
     &      form='unformatted',access='direct',recl=4*nx*ny)
          open (92,file='extra_met/grid_lon_time.gdat',&
     &      form='unformatted',access='direct',recl=4*nx*ny)
        endif

        if (dt.eq.86400.0) then 
          irec = iter 
        elseif (dt.eq.10800.0) then 
          call get_daily_irec (iter,dt,irec_day)
          irec = irec_day
        else
          print *,'This has not been set up to work on dt'
          print *,'  values other than 1 day and 3-hours.'
          print *,'  dt =',dt
          stop
        endif

        read (91,rec=irec) ((xlat_grid(i,j),i=1,nx),j=1,ny)
        read (92,rec=irec) ((xlon_grid(i,j),i=1,nx),j=1,ny)

      endif

! TEMPERATURE.
      if (i_tair_flag.eq.1) then
!       print *,'   solving for temperature'
        call temperature(nx,ny,deltax,deltay,xmn,ymn,&
     &    nstns_orig,xstn_orig,ystn_orig,Tair_orig,dn,Tair_grid,&
     &    undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,&
     &    topo,T_lapse_rate,barnes_lg_domain,n_stns_used,k_stn,&
     &    snowmodel_line_flag,xg_line,yg_line,seaice_run)
      endif

! RELATIVE HUMIDITY.
      if (i_rh_flag.eq.1) then
!       print *,'   solving for relative humidity'
        call relative_humidity(nx,ny,deltax,deltay,xmn,ymn,&
     &    nstns_orig,xstn_orig,ystn_orig,rh_orig,dn,rh_grid,undef,&
     &    ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,topo,&
     &    Tair_orig,Tair_grid,Td_lapse_rate,barnes_lg_domain,&
     &    n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line,&
     &    seaice_run)
      endif

! WIND SPEED AND DIRECTION.
      if (i_wind_flag.eq.1) then
!       print *,'   solving for wind speed and direction'
        call wind(nx,ny,deltax,deltay,xmn,ymn,windspd_orig,&
     &    nstns_orig,xstn_orig,ystn_orig,dn,undef,ifill,&
     &    iobsint,iyear,imonth,iday,xhour,elev_orig,&
     &    winddir_orig,uwind_grid,vwind_grid,slopewt,curvewt,&
     &    curvature,slope_az,terrain_slope,windspd_grid,&
     &    winddir_grid,windspd_flag,winddir_flag,windspd_min,&
     &    vegtype,forest_LAI,calc_subcanopy_met,vegsnowd_xy,&
     &    barnes_lg_domain,n_stns_used,k_stn,snowmodel_line_flag,&
     &    xg_line,yg_line,topo_ref_grid,topo,wind_lapse_rate,&
     &    curve_len_scale,seaice_run)

! Provide the ability to read in an alternate wind dataset that
!   has been previously generated with another program, like
!   NUATMOS.  The following assumes you are reading in a single
!   file with u and v values.  The file name is hardcoded here.
      elseif (i_wind_flag.eq.-1) then
        call read_wind_file(nx,ny,iter,uwind_grid,vwind_grid,&
     &    windspd_grid,winddir_grid,windspd_flag,winddir_flag,&
     &    windspd_min)
      endif

! SOLAR RADIATION.
      if (i_solar_flag.eq.1) then
!       print *,'   solving for solar radiation'
        call solar(nx,ny,xhour,J_day,topo,rh_grid,Tair_grid,&
     &    xlat_grid,Qsi_grid,slope_az,terrain_slope,dt,vegtype,&
     &    forest_LAI,T_lapse_rate,Td_lapse_rate,&
     &    calc_subcanopy_met,gap_frac,cloud_frac_factor,UTC_flag,&
     &    xlon_grid,cloud_frac_grid)

! If requested, modify the model output to account for shortwave
!   radiation observations.
        if (use_shortwave_obs.eq.1.0) then
          if (barnes_lg_domain.eq.1.0) then
            print *,'The model is not configured to assimilate'
            print *,'  solar data with barnes_lg_domain = 1.0.'
            stop
          endif
          call shortwave_data(nx,ny,deltax,deltay,xmn,ymn,&
     &      iyear,imonth,iday,xhour,undef,Qsi_grid,iter)
        endif
      endif

! INCOMING LONGWAVE RADIATION.
      if (i_longwave_flag.eq.1) then
!       print *,'   solving for incoming longwave radiation'
        call longwave(nx,ny,rh_grid,Tair_grid,Qli_grid,topo,&
     &    vegtype,forest_LAI,T_lapse_rate,Td_lapse_rate,&
     &    calc_subcanopy_met,cloud_frac_factor)

! If requested, modify the model output to account for longwave
!   radiation observations.
        if (use_longwave_obs.eq.1.0) then
          if (barnes_lg_domain.eq.1.0) then
            print *,'The model is not configured to assimilate'
            print *,'  longwave data with barnes_lg_domain = 1.0.'
            stop
          endif
          call longwave_data(nx,ny,deltax,deltay,xmn,ymn,&
     &      iyear,imonth,iday,xhour,undef,Qli_grid,iter)
        endif
      endif

! PRECIPITATION.
      if (i_prec_flag.eq.1) then
!       print *,'   solving for precipitation'
        call precipitation(nx,ny,deltax,deltay,xmn,ymn,&
     &    nstns_orig,xstn_orig,ystn_orig,prec_orig,dn,prec_grid,&
     &    undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,&
     &    topo,Tair_grid,sprec,corr_factor,icorr_factor_index,iter,&
     &    precip_lapse_rate,barnes_lg_domain,n_stns_used,k_stn,&
     &    snowmodel_line_flag,xg_line,yg_line,topo_ref_grid,&
     &    iprecip_scheme,cf_precip_flag,cf_precip,snowfall_frac,&
     &    seaice_run)
      endif

! SURFACE PRESSURE.
! Surface pressure is used in EnBal and SnowMass.  If needed for
!   this SnowModel simulation, calculate the distribution here.
      if (run_enbal.eq.1.0 .or. run_snowpack.eq.1.0) then
        call pressure(nx,ny,topo,sfc_pressure)

! If requested, modify the model output to account for surface
!   pressure observations.
        if (use_sfc_pressure_obs.eq.1.0) then
          if (barnes_lg_domain.eq.1.0) then
            print *,'The model is not configured to assimilate'
            print *,'  pressure data with barnes_lg_domain = 1.0.'
            stop
          endif
          call sfc_pressure_data(nx,ny,deltax,deltay,xmn,ymn,&
     &      iyear,imonth,iday,xhour,undef,sfc_pressure,iter)
        endif
      endif

      return
      end subroutine MICROMET_CODE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine precipitation(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns_orig,xstn_orig,ystn_orig,prec_orig,dn,prec_grid,&
     &  undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,&
     &  topo,Tair_grid,sprec,corr_factor,icorr_factor_index,iter,&
     &  precip_lapse_rate,barnes_lg_domain,n_stns_used,k_stn,&
     &  snowmodel_line_flag,xg_line,yg_line,topo_ref_grid,&
     &  iprecip_scheme,cf_precip_flag,cf_precip,snowfall_frac,&
     &  seaice_run)

! Interpolate the observed precipitation values to the grid.  Also
!   interpolate the station elevations to a reference surface.  Use
!   a precipitation "lapse rate", or adjustment factor to define
!   the precipitation on the actual elevation grid.  The reason the
!   interpolated station elevations are used as the topographic
!   reference surface (instead of something like sea level), is
!   because the precipitation adjustment factor is a non-linear 
!   function of elevation difference.
 
! The adjustment factor that is used comes from: Thornton, P. E.,
!   S. W. Running, and M. A. White, 1997: Generating surfaces of
!   daily meteorological variables over large regions of complex
!   terrain.  J. Hydrology, 190, 214-251.

      use snowmodel_inc
      implicit none

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real prec(nstns_max) ! input values
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation
      real prec_orig(nstns_max) ! input values

      real dn                   ! average observation spacing
      real topo(nx,ny)          ! grid topography
      real prec_grid(nx,ny)     ! output values
      real Tair_grid(nx,ny)     ! input values
      real sprec(nx,ny)         ! output values
      real topo_ref_grid(nx,ny) ! reference surface

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real delta_topo,alfa,Tf,precip_lapse_rate_m,precip_lapse_rate,&
     &  barnes_lg_domain
      integer i,j,iter,n_stns_used
      integer k_stn(nx,ny,9)

      real corr_factor(nx_max,ny_max,max_obs_dates+1)
      integer icorr_factor_index(max_time_steps)
      integer iprecip_scheme

      real cf_precip(nx,ny)
      real cf_precip_flag

      real Tair_C,Tair_C_center,slope,b,snowfall_frac
      real snowfall_frac_1,snowfall_frac_2,snowfall_frac_3

      real seaice_run

! Filter through the original input data, and eliminate any
!   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,prec_orig,prec)

! Use the barnes oi scheme to interpolate the station elevation data
!   to the grid, so that it can be used as a topographic reference
!   surface.
      call interpolate(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,elev,dn,topo_ref_grid,undef,ifill,iobsint,&
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,&
     &  k_stn,snowmodel_line_flag,xg_line,yg_line,seaice_run)

! Use the barnes oi scheme to interpolate the station data to
!   the grid.
      call interpolate(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,prec,dn,prec_grid,undef,ifill,iobsint,&
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,&
     &  k_stn,snowmodel_line_flag,xg_line,yg_line,seaice_run)

! Convert the precipitation "lapse rate" (km-1) to m-1.
      precip_lapse_rate_m = precip_lapse_rate / 1000.0

! Choose between Glen's original precipitation increase with
!   elevation scheme, and Ward van Pelt's scheme used in our
!   Svalbard simulations (see van Pelt et al. 2016).

! This is Glen's original MicroMet precipitation adjustment
!   scheme.
      if (iprecip_scheme.eq.1) then

        do j=1,ny
          do i=1,nx

! Convert the gridded station data to the SnowModel-grid elevations.
            delta_topo = topo(i,j) - topo_ref_grid(i,j)

! Don't let the elevation difference be greater than some number
!  (like 1800 meters gives a factor of 4.4).  If it is too large
!   you get huge precipitation adjustment, a divide by zero, or
!   even negative adjustments for high elevations).
            delta_topo = min(delta_topo,1800.0)
            alfa = precip_lapse_rate_m * delta_topo
            prec_grid(i,j) = prec_grid(i,j) * (1.0 + alfa)/(1.0 - alfa)

          enddo
        enddo

! This is van Pelt's precipitation adjustment scheme.
      elseif (iprecip_scheme.eq.2) then

        do j=1,ny
          do i=1,nx
     
! Don't correct precipitation above 1000 m a.s.l..
            if (topo_ref_grid(i,j).le.1000.0) then
              if (topo(i,j).gt.1000.0) then
                delta_topo = 1000.0 - topo_ref_grid(i,j)
              else
                delta_topo = topo(i,j) - topo_ref_grid(i,j)
              endif
            endif
            if (topo_ref_grid(i,j).gt.1000.0) then
              if (topo(i,j).gt.1000.0) then
                delta_topo = 0.0 
              else
                 delta_topo = topo(i,j) - 1000.0
              endif
            endif
                
! Don't let the elevation difference be greater than some number
!  (like 1800 meters gives a factor of 4.4).  If it is too large
!   you get huge precipitation adjustment, a divide by zero, or
!   even negative adjustments for high elevations).
            delta_topo = min(delta_topo,1800.0)
            alfa = 1.75 * delta_topo / 1000.0
            prec_grid(i,j) = prec_grid(i,j) * max(1.0+alfa,0.1)

          enddo
        enddo

      endif

! Convert the precipitation values from mm to m swe.  Also, make
!   sure the interpolation has not created any negetive
!   precipitation values.
      do j=1,ny
        do i=1,nx
          prec_grid(i,j) = prec_grid(i,j) / 1000.0
          prec_grid(i,j) = max(0.0,prec_grid(i,j))
        enddo
      enddo

! This is my original code:
! Use the temperature distribution to define whether this
!   precipitation is falling as rain or snow (generate a
!   snow-precipitation array following the air temperature
!   threshold defined by Auer (1974) = 2.0 C).  Note here that,
!   if you ever want it, rain_prec = prec_grid - sprec.  This
!   snow precipitation (sprec) is assumed to be in meters
!   snow-water-equivalent per time step.
!     Tf = 273.15
!     do j=1,ny
!       do i=1,nx
!         if (Tair_grid(i,j).lt.2.0+Tf) then
!           if (icorr_factor_index(iter).gt.0) then
!             prec_grid(i,j) =
!    &          corr_factor(i,j,icorr_factor_index(iter)) *
!    &          prec_grid(i,j)
!           else
!             prec_grid(i,j) = prec_grid(i,j)
!           endif
!           sprec(i,j) = prec_grid(i,j)
!         else
!           sprec(i,j) = 0.0
!         endif
!       enddo
!     enddo

! First apply any precipitation correction factors calculated as
!   part of any data assimilation.
      do j=1,ny
        do i=1,nx
          if (icorr_factor_index(iter).gt.0) then
            prec_grid(i,j) = &
     &        corr_factor(i,j,icorr_factor_index(iter)) * &
     &        prec_grid(i,j)
          endif
        enddo
      enddo

! Apply the user-defined precipitation correction factor, if one
!   exists.
      if (cf_precip_flag.ne.0.0) then
        do j=1,ny
          do i=1,nx
            prec_grid(i,j) = cf_precip(i,j) * prec_grid(i,j)
          enddo
        enddo
      endif

! Now calculate whether the precipiation is falling as rain or
!   snow, and how much of each.
      Tf = 273.15

! Auer (1974); a rain-snow threshold at +2.0 C.
      if (snowfall_frac.eq.1.0) then

        do j=1,ny
          do i=1,nx
            Tair_C = Tair_grid(i,j) - Tf
            if (Tair_C.lt.2.0) then
              snowfall_frac_1 = 1.0
            else
              snowfall_frac_1 = 0.0
            endif
            sprec(i,j) = snowfall_frac_1 * prec_grid(i,j)
          enddo
        enddo

! Dai, A. (2008): Temperature and pressure dependence of the
!   rain-snow phase transition over land and ocean, Geophys.
!   Res. Lett., 35, L12802, doi:10.1029/2008GL033295.
! In this implementation I have clipped Dai's values to 0
!   and 1.
      elseif (snowfall_frac.eq.2.0) then

        do j=1,ny
          do i=1,nx
            Tair_C = Tair_grid(i,j) - Tf
            snowfall_frac_2 = - 0.482292 * &
     &        (tanh(0.7205 * (Tair_C - 1.1662)) - 1.0223) 
            if (Tair_C.lt.-4.0) then
              snowfall_frac_2 = 1.0
            elseif (Tair_C.gt.6.0) then
              snowfall_frac_2 = 0.0
            endif
            sprec(i,j) = snowfall_frac_2 * prec_grid(i,j)
          enddo
        enddo

! Glen's linear approximation to Dai (2008).  This plots right
!   over the top of Dai, between frac = 0.1 and 0.9.
      elseif (snowfall_frac.eq.3.0) then

! Define where you want the center temperature to be when frac = 0.5.
        Tair_C_center = 1.1
! Define the slope of the line.
        slope = -0.30
! Calculate the intercept (the 0.5 is the frac for Tair_C_center).
        b = 0.5 - slope * Tair_C_center

        do j=1,ny
          do i=1,nx
            Tair_C = Tair_grid(i,j) - Tf

! Solve the equation in the form y = m*x + b
            snowfall_frac_3 = slope * Tair_C + b
            snowfall_frac_3 = max(0.0,snowfall_frac_3)
            snowfall_frac_3 = min(1.0,snowfall_frac_3)

            sprec(i,j) = snowfall_frac_3 * prec_grid(i,j)
          enddo
        enddo

      endif

      return
      end subroutine precipitation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine longwave(nx,ny,rh_grid,Tair_grid,Qli_grid,topo,&
     &  vegtype,forest_LAI,T_lapse_rate,Td_lapse_rate,&
     &  calc_subcanopy_met,cloud_frac_factor)

      use snowmodel_inc
      implicit none

      integer nx       ! number of x output values
      integer ny       ! number of y output values

      real Tair_grid(nx,ny)
      real rh_grid(nx,ny)
      real Qli_grid(nx,ny)
      real topo(nx,ny)
      real vegtype(nx,ny)

      real T_lapse_rate,Td_lapse_rate,es,e,emiss_cloud,&
     &  A,B,C,Tf,Stef_Boltz,cloud_frac,E1,X1,Y1,Z1,E2,X2,Y2,Z2,&
     &  Xs,Ys,Zs,forest_frac,E3,X3,Y3,Z3,alfa,calc_subcanopy_met,&
     &  cloud_frac_factor

      integer, parameter :: nftypes = 5

      real forest_LAI(nftypes)

      integer i,j,nveg

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: temperatures for Buck`s equations are in deg C, and
!   vapor pressures are in mb.  Do the adjustments so that the
!   calculations are done with temperatures in K, and vapor
!   pressures in Pa.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

! Over ice.
!     A = 6.1115 * 100.0
!     B = 22.452
!     C = 272.55

! Define the freezing temperature to be used to convert from C to K.
      Tf = 273.15

! Define the Stefan Boltzmann constant.
      Stef_Boltz = 5.6696e-8

! Constants required for Iziomon et al. (2003).
      E1 = 200.0
      X1 = 0.35
      Y1 = 0.100
      Z1 = 0.224

      E2 = 1500.0
      X2 = 0.43
      Y2 = 0.115
      Z2 = 0.320

! Assume the X and Y coefficients increase linearly to 3000 m,
!   and adjust Z to create a best fit to the CLPX data.
      E3 = 3000.0
      X3 = 0.51
      Y3 = 0.130
      Z3 = 1.100

      do j=1,ny
        do i=1,nx

! Compute the cloud fraction.
          call get_cloudfrac(Tair_grid(i,j),rh_grid(i,j),topo(i,j),&
     &      cloud_frac,T_lapse_rate,Td_lapse_rate,cloud_frac_factor)

! Calculate the vapor pressure.
          es = A * exp((B * (Tair_grid(i,j) - Tf))/ &
     &      (C + (Tair_grid(i,j) - Tf)))
          e = es * rh_grid(i,j) / 100.0

! Compute Qli following Iziomon et al. (2003).
          if (topo(i,j).lt.E1) then
            Xs = X1
            Ys = Y1
            Zs = Z1
          elseif (topo(i,j).gt.E2) then
            Xs = X3
            Ys = Y3
            Zs = Z3
          else
            Xs = X1 + (topo(i,j) - E1) * (X3 - X1)/(E3 - E1)
            Ys = Y1 + (topo(i,j) - E1) * (Y3 - Y1)/(E3 - E1)
            Zs = Z1 + (topo(i,j) - E1) * (Z3 - Z1)/(E3 - E1)
          endif

          alfa = 1.083
          emiss_cloud = alfa * &
     &      (1.0 - Xs * exp((- Ys) * e/Tair_grid(i,j))) * &
     &      (1.0 + Zs * cloud_frac**2)
          emiss_cloud = min(1.0,emiss_cloud)

          Qli_grid(i,j) = emiss_cloud * Stef_Boltz * Tair_grid(i,j)**4

! Modify the incoming longwave radiation for the forest canopy.
          if (vegtype(i,j).le.5.0) then
            if (calc_subcanopy_met.eq.1.0) then

! Define the forest-canopy parameters.
              nveg = nint(vegtype(i,j))
              if (forest_LAI(nveg).lt.0.2) then
                forest_frac = 0.5 * forest_LAI(nveg)
              else
                forest_frac = &
     &            min(1.0,0.55 + 0.29 * log(forest_LAI(nveg)))
              endif
              Qli_grid(i,j) = Qli_grid(i,j) * (1.0 - forest_frac) + &
     &          (Stef_Boltz * Tair_grid(i,j)**4) * forest_frac
            endif
          endif

        enddo
      enddo

      return
      end subroutine longwave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine solar(nx,ny,xhour,J_day,topo,rh_grid,Tair_grid,&
     &  xlat_grid,Qsi_grid,slope_az,terrain_slope,dt,vegtype,&
     &  forest_LAI,T_lapse_rate,Td_lapse_rate,&
     &  calc_subcanopy_met,gap_frac,cloud_frac_factor,UTC_flag,&
     &  xlon_grid,cloud_frac_grid)

! First take the surface gridded fields of Tair and RH, and
!   calculate Td for the topographic surface.  Then use Tair and
!   Td, and the associated lapse rates, and calculate Tair and Td
!   for the 700 mb level.  Use these surfaces to calculate RH at
!   700 mb, and convert these values to cloud fraction following
!   Walcek, C. J., 1994: Cloud cover and its relationship to
!   relative humidity during a spring midlatitude cyclone.  Mon.
!   Wea. Rev., 122, 1021-1035.

      use snowmodel_inc
      implicit none

      integer nx       ! number of x output values
      integer ny       ! number of y output values

      real topo(nx,ny)
      real Tair_grid(nx,ny)
      real rh_grid(nx,ny)
      real Qsi_grid(nx,ny)
      real slope_az(nx,ny)
      real terrain_slope(nx,ny)
      real vegtype(nx,ny)
      real xlat_grid(nx,ny)
      real xlon_grid(nx,ny)
      real cloud_frac_grid(nx,ny)

      real xhour                 ! model decimal hour
      real xxhour
      integer J_day              ! model day of year
      integer i,j,ihrs_day,ihour,nveg
      real dt,cloud_frac,Qsi_tmp,Qsi_sum,trans_veg,&
     &  T_lapse_rate,Td_lapse_rate,calc_subcanopy_met,gap_frac,&
     &  cloud_frac_factor,UTC_flag

      integer, parameter :: nftypes = 5
      real forest_LAI(nftypes)

      ihrs_day = 24

      do j=1,ny
        do i=1,nx

! Compute the cloud fraction.
          call get_cloudfrac(Tair_grid(i,j),rh_grid(i,j),topo(i,j),&
     &      cloud_frac,T_lapse_rate,Td_lapse_rate,cloud_frac_factor)

! Save a copy of the cloud fraction distribution.
          cloud_frac_grid(i,j) = cloud_frac

! Compute the incoming solar radiation.  The solar_rad subroutine
!   calculates the instantaneous incoming solar radiation, so
!   if the time step is very long, account for this by calculating
!   the incoming solar radiation every 3 hours and then taking the
!   average.
          if (dt.le.10800.0) then
            call solar_rad(Qsi_grid(i,j),J_day,xlat_grid(i,j),&
     &        cloud_frac,xhour,slope_az(i,j),terrain_slope(i,j),&
     &        UTC_flag,xlon_grid(i,j))
          elseif (dt.eq.86400.0) then
            Qsi_sum = 0.0
            do ihour=3,ihrs_day,3
              xxhour = real(ihour)
              call solar_rad(Qsi_tmp,J_day,xlat_grid(i,j),&
     &          cloud_frac,xxhour,slope_az(i,j),terrain_slope(i,j),&
     &          UTC_flag,xlon_grid(i,j))
                Qsi_sum = Qsi_sum + Qsi_tmp
            enddo
            Qsi_grid(i,j) = Qsi_sum / (real(ihrs_day)/3.0)
          else
            print *,'The model may not do what you want with this dt'
            stop
          endif

! Modify the incoming solar radiation for the forest canopy.
          if (vegtype(i,j).le.5.0) then
            if (calc_subcanopy_met.eq.1.0) then

! Define the forest-canopy transmissivity.  0.71 provided a
!   best-fit to the observations, when averaged over the two years
!   of hourly data.
              nveg = nint(vegtype(i,j))
              trans_veg = exp((- 0.71) * forest_LAI(nveg))

! Account for any gaps in the forest canopy that will allow
!   direct incoming solar radiation to reach the snow surface.
              trans_veg = gap_frac * (1.0 - trans_veg) + trans_veg

              Qsi_grid(i,j) = trans_veg * Qsi_grid(i,j)

            endif
          endif

        enddo
      enddo

      return
      end subroutine solar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_cloudfrac(Tair_grid,rh_grid,topo,&
     &  cloud_frac,T_lapse_rate,Td_lapse_rate,cloud_frac_factor)

      implicit none

      real Td_lapse_rate,topo_ref,delta_topo,A,B,C,e,es,dx,&
     &  T_lapse_rate,press_ratio,f_max,one_minus_RHe,f_1,&
     &  Td_700,Tair_700,rh_700,cloud_frac,Tair_grid,rh_grid,&
     &  topo,Td_grid,Tf,cloud_frac_factor

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: temperatures for Buck`s equations are in deg C, and
!   vapor pressures are in mb.  Do the adjustments so that the
!   calculations are done with temperatures in K, and vapor
!   pressures in Pa.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

! Over ice.
!     A = 6.1115 * 100.0
!     B = 22.452
!     C = 272.55

! Define the freezing temperature to be used to convert from C to K.
      Tf = 273.15

! Assume that 700 mb is equivalent to 3000 m in a standard
!   atmosphere.
      topo_ref = 3000.0

! Define the ratio of 700 mb level pressure to the surface pressure
!   (~1000 mb).
      press_ratio = 0.7

! Assume dx = 80.0 km, for Walcek (1994).
      dx = 80.0

! Walcek coefficients.
      f_max = 78.0 + 80.0/15.5
      one_minus_RHe = 0.196 + (0.76-80.0/2834.0) * (1.0 - press_ratio)
      f_1 = f_max * (press_ratio - 0.1) / 0.6 / 100.0

! Convert the gridded topo-surface RH to Td.
      es = A * exp((B * (Tair_grid - Tf))/(C + (Tair_grid - Tf)))
      e = es * max(10.0,rh_grid) / 100.0
      Td_grid = C * log(e/A) / (B - log(e/A)) + Tf

! Convert the topo-surface temperature values to 700 mb values.
!     delta_topo = topo - topo_ref
      delta_topo = topo_ref - topo
      Td_700 = Td_grid + Td_lapse_rate * delta_topo
      Tair_700 = Tair_grid + T_lapse_rate * delta_topo

! Convert each Td to a gridded relative humidity (0-1).
      e = A * exp((B * (Td_700 - Tf))/(C + (Td_700 - Tf)))
      es = A * exp((B * (Tair_700 - Tf))/(C + (Tair_700 - Tf)))
      rh_700 = e/es
      rh_700 = min(1.0,rh_700)
      rh_700 = max(0.0,rh_700)

! Use this RH at 700 mb to define the cloud fraction (0-1).
      cloud_frac = f_1 * exp((rh_700 - 1.0)/one_minus_RHe)

! If the user wants to, reduce the calculate cloud fraction by the
!   cloud_frac_factor.
      cloud_frac = cloud_frac_factor * cloud_frac

      return
      end subroutine get_cloudfrac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine solar_rad(Qsi,J_day,xlat,&
     &  cloud_frac,xxhour,slope_az,terrain_slope,&
     &  UTC_flag,xlon)

      implicit none

      integer J_day

      real solar_const,days_yr,Trop_Can,solstice,pi,deg2rad,&
     &  Qsi_direct,Qsi_diffuse,cos_i,cos_Z,Qsi,xlat,sin_z,xxhour,&
     &  cloud_frac,slope_az,terrain_slope,sol_dec,hr_angl,&
     &  trans_direct,trans_diffuse,Qsi_trans_dir,Qsi_trans_dif,&
     &  sun_azimuth,slope_az_S0,xxxhour,UTC_flag,xlon

! Required constants.
      solar_const = 1370.
      days_yr = 365.25
      Trop_Can = 0.41
      solstice = 173.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0

! COMPUTE THE BASIC SOLAR RADIATION PARAMETERS.

! Compute the solar declination angle (radians).
      sol_dec = Trop_Can * &
     &  cos(2.*pi * (real(J_day) - solstice)/days_yr)
      
! For the case of running UTC time and a latitudinal variation
!   in solar radiation, adjust the time to correspond to the
!   local time at this longitude position.
      if (UTC_flag.ne.0.0) then
        xxxhour = xxhour + xlon / 15.0
        if (xxxhour.ge.24.0) xxxhour = xxxhour - 24.0
        if (xxxhour.lt.0.0) xxxhour = xxxhour + 24.0
      else
        xxxhour = xxhour
      endif

! Compute the sun's hour angle (radians).
      hr_angl = (xxxhour * 15.0 - 180.0) * deg2rad

! Compute cos_Z.  Note that the sin of the solar elevation angle,
!   sin_alfa, is equal to the cosine of the solar zenith angle,
!   cos_Z.
      cos_Z = sin(sol_dec) * sin(xlat * deg2rad) + &
     &  cos(sol_dec) * cos(xlat * deg2rad) * cos(hr_angl)
      cos_Z = max(0.0,cos_Z)

! Account for clouds, water vapor, pollution, etc.
      trans_direct = (0.6 + 0.2 * cos_Z) * (1.0 - cloud_frac)
      trans_diffuse = (0.3 + 0.1 * cos_Z) * cloud_frac

! Compute the solar radiation transmitted through the atmosphere.
      Qsi_trans_dir = solar_const * trans_direct
      Qsi_trans_dif = solar_const * trans_diffuse

! COMPUTE THE CORRECTIONS TO ALLOW FOR TOPOGRAPHIC SLOPE AND ASPECT.

! The sine of the solar zenith angle.
      sin_Z = sqrt(1.0 - cos_Z*cos_Z)

! Azimuth of the sun, with south having zero azimuth for the
!   northern hemisphere.
      sun_azimuth = &
     &  asin(max(-1.0,min(1.0,cos(sol_dec)*sin(hr_angl)/sin_Z)))
      if (xlat.lt.0.0) then
        sun_azimuth = - sun_azimuth
      endif

! Make the corrections so that the angles below the local horizon
!   are still measured from the normal to the slope.
      if (xlat.ge.0.0) then
        if (hr_angl.lt.0.0) then
          if (hr_angl.lt.sun_azimuth) sun_azimuth = - pi - sun_azimuth
        elseif (hr_angl.gt.0.0) then
          if (hr_angl.gt.sun_azimuth) sun_azimuth = pi - sun_azimuth
        endif
      else
!       if (hr_angl.lt.0.0) then
!         if (hr_angl.lt.sun_azimuth) sun_azimuth = - pi - sun_azimuth
!       elseif (hr_angl.gt.0.0) then
!         if (hr_angl.gt.sun_azimuth) sun_azimuth = pi - sun_azimuth
!       endif
      endif

! Build, from the variable with north having zero azimuth, a 
!   slope_azimuth value with south having zero azimuth.  Also
!   make north have zero azimuth if in the southern hemsisphere.
      if (xlat.ge.0.0) then
        if (slope_az.ge.180.0) then
          slope_az_S0 = slope_az - 180.0
        else
          slope_az_S0 = slope_az + 180.0
        endif
      else
        slope_az_S0 = slope_az
      endif

! Compute the angle between the normal to the slope and the angle
!   at which the direct solar radiation impinges on the sloping
!   terrain (radians).
      cos_i = cos(terrain_slope * deg2rad) * cos_Z + &
     &  sin(terrain_slope * deg2rad) * sin_Z * &
     &  cos(sun_azimuth - slope_az_S0 * deg2rad)

! Adjust the topographic correction due to local slope so that
!   the correction is zero if the sun is below the local horizon 
!   (i.e., the slope is in the shade) or if the sun is below the
!   global horizon.
      if (cos_i.lt.0.0) cos_i = 0.0
      if (cos_Z.le.0.0) cos_i = 0.0

! Adjust the solar radiation for slope, etc.
      Qsi_direct = cos_i * Qsi_trans_dir
      Qsi_diffuse = cos_Z * Qsi_trans_dif

! Combine the direct and diffuse solar components.
      Qsi = Qsi_direct + Qsi_diffuse

      return
      end subroutine solar_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine wind(nx,ny,deltax,deltay,xmn,ymn,windspd_orig,&
     &  nstns_orig,xstn_orig,ystn_orig,dn,undef,ifill,&
     &  iobsint,iyear,imonth,iday,xhour,elev_orig,&
     &  winddir_orig,uwind_grid,vwind_grid,slopewt,curvewt,&
     &  curvature,slope_az,terrain_slope,windspd_grid,&
     &  winddir_grid,windspd_flag,winddir_flag,windspd_min,&
     &  vegtype,forest_LAI,calc_subcanopy_met,vegsnowd_xy,&
     &  barnes_lg_domain,n_stns_used,k_stn,snowmodel_line_flag,&
     &  xg_line,yg_line,topo_ref_grid,topo,wind_lapse_rate,&
     &  curve_len_scale,seaice_run)

! This program takes the station wind speed and direction, converts
!   them to u and v components, interpolates u and v to a grid,
!   converts the gridded values to speed and direction, and then
!   runs a simple wind model that adjusts those speeds and
!   directions according to topographic slope and curvature
!   relationships.  The resulting speeds and directions are
!   converted to u and v components and passed back to the main
!   program to be written to a file.  (All of the conversion
!   between u-v and speed-dir is done because of the problems
!   with interpolating over the 360/0 direction line.)

      use snowmodel_inc
      ! KRA
      use LIS_coreMod
      use LIS_mpiMod
      use snowmodel_lsmMod, only : snowmodel_struc
      ! KRA
      implicit none

      integer :: ierr  ! KRA

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation

      real windspd_orig(nstns_max) ! input values
      real winddir_orig(nstns_max) ! input values

      real speed(nstns_max)  ! input values
      real dir(nstns_max)    ! input values
      real u(nstns_max)      ! u component of wind
      real v(nstns_max)      ! v component of wind

      real dn                  ! average observation spacing
      real uwind_grid(nx,ny)  ! output values
      real vwind_grid(nx,ny)  ! output values
      real u_grid(nx,ny)  ! temporary u wind component
      real v_grid(nx,ny)  ! temporary v wind component
      real winddir_grid(nx,ny)  ! temporary wind direction
      real windspd_grid(nx,ny)  ! temporary wind speed

      real curvature(nx,ny)     ! topographic curvature
      real slope_az(nx,ny)      ! azimuth of topographic slope
      real terrain_slope(nx,ny) ! terrain slope
      real vegtype(nx,ny)
      real vegsnowd_xy(nx,ny)
      real topo_ref_grid(nx,ny) ! reference surface
      real topo(nx,ny)

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real pi,deg2rad,rad2deg,slopewt,curvewt,curve_len_scale
      integer i,j,k,n_stns_used
      integer k_stn(nx,ny,9)
      real windspd_flag,winddir_flag,u_sum,v_sum,windspd_min,&
     &  calc_subcanopy_met,barnes_lg_domain,wind_lapse_rate,&
     &  delta_topo,alfa1,alfa2

      integer, parameter ::  nftypes = 5
      real forest_LAI(nftypes)
      real seaice_run

! Define the required constants.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0
      rad2deg = 180.0 / pi

! Filter through the original input data, and eliminate any
!   missing values (i.e., make sure each wind direction is paired
!   up with a wind speed.
      call get_good_values2(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,windspd_orig,winddir_orig,&
     &  speed,dir)

! Convert these station data to u and v wind components.
      do k=1,nstns
        speed(k) = max(windspd_min,speed(k))
        u(k) = (- speed(k)) * sin(deg2rad * dir(k))
        v(k) = (- speed(k)) * cos(deg2rad * dir(k))
      enddo

! Use the barnes oi scheme to interpolate the station data to
!   the grid.
! U component.
      call interpolate(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,u,dn,u_grid,undef,ifill,iobsint,&
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,&
     &  k_stn,snowmodel_line_flag,xg_line,yg_line,seaice_run)

! V component.
      call interpolate(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,v,dn,v_grid,undef,ifill,iobsint,&
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,&
     &  k_stn,snowmodel_line_flag,xg_line,yg_line,seaice_run)

! If desired, impose a wind speed increase with elevation.  Here
!   the wind_lapse_rate = the wind speed increase per 1-km elevation
!   gain.  The adjustment is patterned after the precipitation-
!   elevation adjustment.  topo_ref_grid here comes from the
!   precipitation adjustment.
      if (wind_lapse_rate.ne.0.0) then
        alfa1 = (wind_lapse_rate - 1.0) / (1.0 + wind_lapse_rate)
! Convert to m-1.
        alfa1 = alfa1 / 1000.0
        do j=1,ny
          do i=1,nx
            delta_topo = topo(i,j) - topo_ref_grid(i,j)
! Impose some limits to the adjustment.
            delta_topo = min(delta_topo,1800.0)
            alfa2 = alfa1 * delta_topo
            u_grid(i,j) = u_grid(i,j) * (1.0 + alfa2)/(1.0 - alfa2)
            v_grid(i,j) = v_grid(i,j) * (1.0 + alfa2)/(1.0 - alfa2)
          enddo
        enddo
      endif

! Convert these u and v components to speed and directions.
      do j=1,ny
        do i=1,nx

! Some compilers do not allow both u and v to be 0.0 in
!   the atan2 computation.
          if (abs(u_grid(i,j)).lt.1e-10) u_grid(i,j) = 1e-10

          winddir_grid(i,j) = rad2deg * atan2(u_grid(i,j),v_grid(i,j))
          if (winddir_grid(i,j).ge.180.0) then
            winddir_grid(i,j) = winddir_grid(i,j) - 180.0
          else
            winddir_grid(i,j) = winddir_grid(i,j) + 180.0
          endif

!         winddir_grid(i,j) = 270.0 -
!    &      rad2deg*atan2(v_grid(i,j),u_grid(i,j))
!         if (winddir_grid(i,j).ge.360.0)
!    &      winddir_grid(i,j) = winddir_grid(i,j)-360.0

          windspd_grid(i,j) = sqrt(u_grid(i,j)**2 + v_grid(i,j)**2)
        enddo
      enddo

! Modify the wind speed and direction according to simple
!   wind-topography relationships.
      call topo_mod_winds(nx,ny,winddir_grid,slopewt,curvewt,&
     &  windspd_grid,uwind_grid,vwind_grid,curvature,slope_az,&
     &  terrain_slope,vegtype,forest_LAI,calc_subcanopy_met,&
     &  vegsnowd_xy,curve_len_scale,deltax,deltay)

! Avoid problems of zero (low) winds (for example, turbulence
!   theory, log wind profile, etc., says that we must have some
!   wind.  Thus, some equations blow up when the wind speed gets
!   very small).
      do j=1,ny
        do i=1,nx
          if (windspd_grid(i,j).lt.windspd_min) then
            windspd_grid(i,j) = windspd_min
            uwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &        sin(deg2rad*winddir_grid(i,j))
            vwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &        cos(deg2rad*winddir_grid(i,j))
          endif
        enddo
      enddo

! Find the maximum wind speed in the domain, and the
!   domain-averaged wind direction.
      windspd_flag = 0.0
      u_sum = 0.0
      v_sum = 0.0
      do j=1,ny
        do i=1,nx
          windspd_flag = max(windspd_flag,windspd_grid(i,j))
          u_sum = u_sum + uwind_grid(i,j)
          v_sum = v_sum + vwind_grid(i,j)
        enddo
      enddo

! KRA
#if (defined SPMD)
! CAF calls for locating totals ... (Alessandro)
!  call co_sum(u_sum)
!  call co_sum(v_sum)
!  call co_max(windspd_flag)

!  print *, "wind: ",LIS_localPet+1, u_sum, v_sum, windspd_flag
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call MPI_ALLREDUCE(u_sum, snowmodel_struc(1)%usum_glb, 1,&
           MPI_REAL, MPI_SUM,&
           LIS_mpi_comm, ierr)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call MPI_ALLREDUCE(v_sum, snowmodel_struc(1)%vsum_glb, 1,&
           MPI_REAL, MPI_SUM,&
           LIS_mpi_comm, ierr)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call MPI_ALLREDUCE(windspd_flag, snowmodel_struc(1)%windspdflg_glb, 1,&
           MPI_REAL, MPI_MAX,&
           LIS_mpi_comm, ierr)

      u_sum = snowmodel_struc(1)%usum_glb
      v_sum = snowmodel_struc(1)%vsum_glb
      windspd_flag = snowmodel_struc(1)%windspdflg_glb

!  print *, " final u_sum: ",u_sum
!  print *, " final v_sum: ",v_sum
!  print *, " final wspdflag: ",windspd_flag

#endif
! KRA

      u_sum = u_sum / real(nx*ny)
      v_sum = v_sum / real(nx*ny)

! Some compilers do not allow both u and v to be 0.0 in
!   the atan2 computation.
      if (abs(u_sum).lt.1e-10) u_sum = 1e-10

      winddir_flag = rad2deg * atan2(u_sum,v_sum)
      if (winddir_flag.ge.180.0) then
        winddir_flag = winddir_flag - 180.0
      else
        winddir_flag = winddir_flag + 180.0
      endif

!     winddir_flag = 270.0 - rad2deg*atan2(v_sum,u_sum)
!     if (winddir_flag.ge.360.0) winddir_flag = winddir_flag-360.0

      return
      end subroutine wind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine relative_humidity(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns_orig,xstn_orig,ystn_orig,rh_orig,dn,rh_grid,undef,&
     &  ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,topo,&
     &  Tair_orig,Tair_grid,Td_lapse_rate,barnes_lg_domain,&
     &  n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line,&
     &  seaice_run)

! This procedure follows: Kunkel, K. E., 1989: Simple procedures for
!   extrapolation of humidity variables in the mountainous Western
!   United States. J. Climate, 2, 656-669.

! First convert stn relative humidity to dew-point temperature.  Use
!   the Td lapse rate to take the stn Td to sea level.  Interpolate
!   the stn Td to the grid.  Use the Td lapse rate to take the sea
!   level grid to the actual elevations.  Convert each Td to
!   relative humidity.

      use snowmodel_inc
      implicit none

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation

      real Tair_orig(nstns_max)  ! input values
      real rh_orig(nstns_max)    ! input values

      real Tair(nstns_max)  ! input values
      real Td(nstns_max)    ! input values
      real rh(nstns_max)    ! input values

      real dn                  ! average observation spacing
      real topo(nx,ny) ! grid topography
      real Tair_grid(nx,ny)  ! output values
      real Td_grid(nx,ny)    ! output values
      real rh_grid(nx,ny)    ! output values

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real Td_lapse_rate,topo_ref,delta_topo,A,B,C,e,es,Tf,&
     &  barnes_lg_domain
      integer i,j,k,n_stns_used
      integer k_stn(nx,ny,9)
      real seaice_run

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: temperatures for Buck`s equations are in deg C.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

! Over ice.
!     A = 6.1115 * 100.0
!     B = 22.452
!     C = 272.55

! Define the freezing temperature to be used to convert from C to K.
      Tf = 273.15

! Filter through the original input data, and eliminate any
!   missing values.
      call get_good_values2(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,Tair_orig,rh_orig,&
     &  Tair,rh)

! Convert the stn relative humidity to Td.
      do k=1,nstns

! Saturation vapor pressure at temperature, T.
        es = A * exp((B * (Tair(k) - Tf))/(C + (Tair(k) - Tf)))

! Dew point temperature for a given temperature and relative humidity.
        e = es * max(10.0,rh(k)) / 100.0
        Td(k) = C * log(e/A) / (B - log(e/A)) + Tf

      enddo

! Define the topographic reference surface.
      topo_ref = 0.0

! Convert the station data to sea level values.
      do k=1,nstns
        delta_topo = topo_ref - elev(k)
        Td(k) = Td(k) + Td_lapse_rate * delta_topo
      enddo

! Use the barnes oi scheme to interpolate the station data to
!   the grid.
      call interpolate(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,Td,dn,Td_grid,undef,ifill,iobsint,&
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,&
     &  k_stn,snowmodel_line_flag,xg_line,yg_line,seaice_run)

! Convert these grid values back to the actual gridded elevations.
      do j=1,ny
        do i=1,nx
          delta_topo = topo(i,j) - topo_ref
          Td_grid(i,j) = Td_grid(i,j) + Td_lapse_rate * delta_topo
        enddo
      enddo

! Convert each Td to a gridded relative humidity.
      do j=1,ny
        do i=1,nx
          e = A * exp((B * (Td_grid(i,j) - Tf)) / &
     &      (C + (Td_grid(i,j) - Tf)))
          es = A * exp((B * (Tair_grid(i,j) - Tf)) /&
     &     (C + (Tair_grid(i,j) - Tf)))
          rh_grid(i,j) = 100.0 * e/es

! Make sure the interpolation processes has not created any values
!   above 100 and below 0.
          rh_grid(i,j) = min(100.0,rh_grid(i,j))
          rh_grid(i,j) = max(0.0,rh_grid(i,j))
        enddo
      enddo

      return
      end subroutine relative_humidity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine temperature(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns_orig,xstn_orig,ystn_orig,Tair_orig,dn,Tair_grid,&
     &  undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,&
     &  topo,T_lapse_rate,barnes_lg_domain,n_stns_used,k_stn,&
     &  snowmodel_line_flag,xg_line,yg_line,seaice_run)

! The lapse rate used depends on the month of the year, and is
!   defined by: Kunkel, K. E., 1989: Simple procedures for
!   extrapolation of humidity variables in the mountainous Western
!   United States. J. Climate, 2, 656-669.

! First adjust the stn temperatures to a common level (sea level),
!   assuming this lapse rate.  Then interpolate the temperatures
!   to the model grid.  Then use the topography data and lapse
!   rate to adjust the gridded temperatures values back to the
!   actual elevation.

! Contact Glen if you are interested in the temperature inversion
!   code presented in: Mernild, S. H., and G. E. Liston, 2010:
!   The influence of air temperature inversions on snowmelt and
!   glacier mass-balance simulations, Ammassalik Island, SE
!   Greenland. J. Applied Meteorology and Climatology, 49, 47-67.

      use snowmodel_inc
      use LIS_coreMod   ! KRA

      implicit none

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real Tair(nstns_max) ! input values
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation
      real Tair_orig(nstns_max) ! input values

      real dn                  ! average observation spacing
      real topo(nx,ny) ! grid topography
      real Tair_grid(nx,ny) ! output values

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real T_lapse_rate,topo_ref,delta_topo,barnes_lg_domain
      integer i,j,k,n_stns_used
      integer k_stn(nx,ny,9)
      real seaice_run

! Filter through the original input data, and eliminate any
!   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,Tair_orig,Tair)

! Define the topographic reference surface.
      topo_ref = 0.0

! Convert the station data to sea level values.
      do k=1,nstns
        delta_topo = topo_ref - elev(k)
        Tair(k) = Tair(k) + T_lapse_rate * delta_topo
      enddo

! Use the barnes oi scheme to interpolate the station data to
!   the grid.
      call interpolate(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,Tair,dn,Tair_grid,undef,ifill,iobsint,&
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,&
     &  k_stn,snowmodel_line_flag,xg_line,yg_line,seaice_run)

! Convert these grid values back to the actual gridded elevations.
      do j=1,ny
        do i=1,nx
          delta_topo = topo(i,j) - topo_ref
          Tair_grid(i,j) = Tair_grid(i,j) + T_lapse_rate * delta_topo
        enddo
      enddo

      return
      end subroutine temperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine topo_mod_winds(nx,ny,winddir_grid,slopewt,curvewt,&
     &  windspd_grid,uwind_grid,vwind_grid,curvature,slope_az,&
     &  terrain_slope,vegtype,forest_LAI,calc_subcanopy_met,&
     &  vegsnowd_xy,curve_len_scale,deltax,deltay)

      use snowmodel_inc
! KRA
      use LIS_coreMod
      use LIS_mpiMod
      use snowmodel_lsmMod, only : snowmodel_struc
! KRA
      implicit none

      integer :: ierr   ! KRA

      integer i,j,nx,ny,nveg,k

      real pi,deg2rad,rad2deg,slopewt,curvewt,dirdiff,curve_len_scale,&
     &  wslope_max,beta,veg_ht,a,canopy_windwt,calc_subcanopy_met,&
     &  deltax,deltay,xmult
      integer, parameter :: nftypes = 5
      integer :: loops_windwt_smoother

      real forest_LAI(nftypes)

      real curvature(nx,ny)
      real windspd_grid(nx,ny)
      real winddir_grid(nx,ny)
      real uwind_grid(nx,ny)
      real vwind_grid(nx,ny)
      real wind_slope(nx,ny)
      real slope_az(nx,ny)
      real terrain_slope(nx,ny)
      real vegtype(nx,ny)
      real vegsnowd_xy(nx,ny)
      real windwt(nx,ny)

! Compute the wind modification factor which is a function of
!   topography and wind direction following Liston and Sturm (1998).

! Define the required constants.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0
      rad2deg = 180.0 / pi

! Compute the slope in the direction of the wind.
      do i=1,nx
        do j=1,ny
          wind_slope(i,j) = deg2rad * terrain_slope(i,j) * &
     &      cos(deg2rad * (winddir_grid(i,j) - slope_az(i,j)))
        enddo
      enddo

! Scale the wind slope such that the max abs(wind slope) has a value
!   of abs(0.5).  Include a 1 mm slope in slope_max to prevent
!   divisions by zero in flat terrain where the slope is zero.
      wslope_max = 0.0 + 0.001
      do j=1,ny
        do i=1,nx
          wslope_max = max(wslope_max,abs(wind_slope(i,j)))
        enddo
      enddo

! KRA
#if (defined SPMD)
! call co_max(wslope_max)
!  print *, "wslope_max: ",LIS_localPet+1, wslope_max
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call MPI_ALLREDUCE(wslope_max, snowmodel_struc(1)%wslopemax_glb, 1,&
           MPI_REAL, MPI_MAX,&
           LIS_mpi_comm, ierr)
      wslope_max = snowmodel_struc(1)%wslopemax_glb
!  print *, "final wslope_max: ",wslope_max
#endif
! KRA

      do j=1,ny
        do i=1,nx
          wind_slope(i,j) = wind_slope(i,j) / (2.0 * wslope_max)
        enddo
      enddo

! Calculate the wind speed and direction adjustments.  The
!   curvature and wind_slope values range between -0.5 and +0.5.
!   Valid slopewt and curvewt values are between 0 and 1, with
!   values of 0.5 giving approximately equal weight to slope and
!   curvature.  I suggest that slopewt and curvewt be set such
!   that slopewt + curvewt = 1.0.  This will limit the total
!   wind weight to between 0.5 and 1.5 (but this is not required).

! Compute the wind weighting factor.
      do i=1,nx
        do j=1,ny
          windwt(i,j) = 1.0 + slopewt * wind_slope(i,j) + &
     &      curvewt * curvature(i,j)
        enddo
      enddo

! Smooth the wind weighting factor to eliminate any sharp speed
!   increases resulting from the merging of the curve wt and the
!   slope wt.  Define the number of times this is done to be a
!   function of the curvature length scale and the grid increment.
!   The first 0.5 here just means that half of the caclulated
!   number appears to be about right (additional loops with this
!   smoother does not change the results much).  If there are
!   unwanted wave features in the snow distribution, this 0.5
!   factor can be increased to 1.0 or more, to get rid of these
!   waves.  Also see "loops_snowd_smoother" in snowtran_code.f.
!     xmult = 0.5
      xmult = 1.0
!     xmult = 1.5
      loops_windwt_smoother = nint(xmult * curve_len_scale / &
     &  (0.5 * (deltax + deltay)))

!     print *
!     print *, 'loops_windwt_smoother ',loops_windwt_smoother
!     print *

! Don't do this smoothing if the domain is arbitrarily small.
      if (nx.gt.100 .and. ny.gt.100) then
        do k=1,loops_windwt_smoother
          call smoother9(nx,ny,windwt)
        enddo
      endif

! Continue with the wind calculations.
      do i=1,nx
        do j=1,ny

! Generate the terrain-modified wind speed.
          windspd_grid(i,j) = windwt(i,j) * windspd_grid(i,j)

! Further modify the wind speed to account for forest canopies.
          if (vegtype(i,j).le.5.0) then
            if (calc_subcanopy_met.eq.1.0) then
              nveg = nint(vegtype(i,j))

! Define the canopy wind-weighting factor.  Assume z=0.6*canopy_ht,
!   and the canopy_ht equals the vegetation snow-holding depth.
              beta = 0.9
              veg_ht = vegsnowd_xy(i,j)
              a = beta * forest_LAI(nveg)
              canopy_windwt = exp((- a)*(1.0 - (0.6*veg_ht)/veg_ht))
              windspd_grid(i,j) = canopy_windwt * windspd_grid(i,j)

            endif
          endif

! Modify the wind direction according to Ryan (1977).  Note that it
!   is critical that "dirdiff" handles the cases where the slope
!   azimuth and the wind direction are on different sides of the
!   360-0 line.
          if (slope_az(i,j).gt.270.0.and. &
     &      winddir_grid(i,j).lt.90.0) then
            dirdiff = slope_az(i,j) - winddir_grid(i,j) - 360.0
          elseif (slope_az(i,j).lt.90.0.and. &
     &      winddir_grid(i,j).gt.270.0) then
            dirdiff = slope_az(i,j) - winddir_grid(i,j) + 360.0
          else
            dirdiff = slope_az(i,j) - winddir_grid(i,j)
          endif
          if (abs(dirdiff).le.90.0) then
            winddir_grid(i,j) = winddir_grid(i,j) - 0.5 * &
     &        min(wind_slope(i,j)*rad2deg,45.0) * &
     &        sin(deg2rad * (2.0 * dirdiff))
            if (winddir_grid(i,j).gt.360.0) then
              winddir_grid(i,j) = winddir_grid(i,j) - 360.0
            elseif (winddir_grid(i,j).lt.0.0) then
              winddir_grid(i,j) = winddir_grid(i,j) + 360.0
            endif
          endif

! Extract the u and v wind components.
          uwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &      sin(deg2rad*winddir_grid(i,j))
          vwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &      cos(deg2rad*winddir_grid(i,j))
        enddo
      enddo

      return
      end subroutine topo_mod_winds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine topo_data(nx,ny,deltax,deltay,topo,&
     &  curvature,terrain_slope,slope_az,curve_len_scale)

      use snowmodel_inc
! KRA
      use LIS_coreMod
      use LIS_mpiMod
      use snowmodel_lsmMod, only : snowmodel_struc
! KRA
      implicit none

      integer :: ierr  ! KRA

      integer i,j,nx,ny,inc

      real pi,rad2deg,deltax,deltay,deltaxy,curve_len_scale,curve_max

      real topo(nx,ny)
      real dzdx(nx,ny)
      real dzdy(nx,ny)
      real curvature(nx,ny)
      real slope_az(nx,ny)
      real terrain_slope(nx,ny)

! Compute the topographic information required to run the wind
!   model.

! Deal with the model running at a point, or along single or double
!   lines.
      if (nx.le.2  .or.  ny.le.2) then
        do i=1,nx
          do j=1,ny
            curvature(i,j) = 0.0
            terrain_slope(i,j) = 0.0
            slope_az(i,j) = 0.0
          enddo
        enddo
      else

! Define the required constants.
        pi = 2.0 * acos(0.0)
        rad2deg = 180.0 / pi

! CURVATURE CALCULATIONS.

! Compute the average grid increment.
        deltaxy = 0.5 * (deltax + deltay)

! Convert the length scale to an appropriate grid increment.
        inc = max(1,nint(curve_len_scale/deltaxy))

! Compute the curvature.
        do i=1,nx
          do j=1,ny
            curvature(i,j) = (4.0 * topo(i,j) - &
     &        topo(max(1,i-inc),max(1,j-inc)) - &
     &        topo(min(nx,i+inc),min(ny,j+inc)) - &
     &        topo(min(nx,i+inc),max(1,j-inc)) - &
     &        topo(max(1,i-inc),min(ny,j+inc))) / &
     &        (sqrt(2.0) * 16.0 * real(inc) * deltaxy) + &
     &        (4.0 * topo(i,j) - &
     &        topo(min(nx,i+inc),j) - topo(max(1,i-inc),j) - &
     &        topo(i,min(ny,j+inc)) - topo(i,max(1,j-inc))) / &
     &        (16.0 * real(inc) * deltaxy)
          enddo
        enddo

! Scale the curvature such that the max abs(curvature) has a value
!   of abs(0.5).  Include a 1 mm curvature in curve_max to prevent
!   divisions by zero in flat terrain where the curvature is zero.
        curve_max = 0.0 + 0.001
        do j=1,ny
          do i=1,nx
            curve_max = max(curve_max,abs(curvature(i,j)))
          enddo
        enddo
! KRA
#if (defined SPMD)
! call co_max(curve_max)
        call MPI_Barrier(LIS_MPI_COMM, ierr)
        call MPI_ALLREDUCE(curve_max, snowmodel_struc(1)%curvemax_glb, 1,&
             MPI_REAL, MPI_MAX,&
             LIS_mpi_comm, ierr)
        curve_max = snowmodel_struc(1)%curvemax_glb
#endif
! KRA
        do j=1,ny
          do i=1,nx
            curvature(i,j) = curvature(i,j) / (2.0 * curve_max)
          enddo
        enddo

! SLOPE CALCULATIONS.

! Find dzdx.
        do j=1,ny
          dzdx(1,j) = (topo(2,j) - topo(1,j)) / deltax
          do i=2,nx-1
            dzdx(i,j) = (topo(i+1,j) - topo(i-1,j)) / (2.0 * deltax)
          enddo
          dzdx(nx,j) = (topo(nx,j) - topo(nx-1,j)) / deltax
        enddo

! Find dzdy.
        do i=1,nx
          dzdy(i,1) = (topo(i,2) - topo(i,1)) / deltay
          do j=2,ny-1
            dzdy(i,j) = (topo(i,j+1) - topo(i,j-1)) / (2.0 * deltay)
          enddo
          dzdy(i,ny) = (topo(i,ny) - topo(i,ny-1)) / deltay
        enddo

! Calculate the terrain slope and azimuth.
        do i=1,nx
          do j=1,ny

! Some compilers will not allow dzdx and dzdy to both be 0.0 in
!   the atan2 computation.
!           if (abs(dzdx(i,j)).lt.1e-10) dzdx(i,j) = 1e-10
            if (abs(dzdy(i,j)).lt.1e-10) dzdy(i,j) = 1e-10

! Compute the slope azimuth, making sure that north has zero
!   azimuth.  Also note that for the Ryan wind rotation, the
!   azimuth values must range from 0 to 360.
            slope_az(i,j) = rad2deg * &
     &        (3.0 / 2.0 * pi - atan2(dzdy(i,j),dzdx(i,j)))
            if (slope_az(i,j).ge.360.0) slope_az(i,j) = &
     &        slope_az(i,j) - 360.0

! Compute the slope of the terrain.
            terrain_slope(i,j) = rad2deg * &
     &        atan(sqrt(dzdx(i,j)*dzdx(i,j) + dzdy(i,j)*dzdy(i,j)))

          enddo
        enddo

      endif

      return
      end subroutine topo_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interpolate(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,var,dn,grid,undef,ifill,iobsint,&
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,&
     &  k_stn,snowmodel_line_flag,xg_line,yg_line,seaice_run)

      use snowmodel_inc
      use LIS_coreMod    ! KRA

      implicit none

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real var(nstns_max)  ! input values

      double precision xstn_tmp(nstns_max) ! input stn x coords
      double precision ystn_tmp(nstns_max) ! input stn y coords
      real var_tmp(nstns_max)  ! input values

      integer nstns        ! number of input values, all good
      real undef           ! undefined value
      real dn              ! average observation spacing
      real grid(nx,ny)     ! output values

      integer i,j      ! col, row counters
      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      integer k_stn(nx,ny,9)
      integer k,n_stns_used
      real barnes_lg_domain,seaice_run

! Use the barnes oi scheme to grid the station data. If there is
!   only a single station, distribute those data uniformly over
!   the grid.  In the event that there are no valid observations,
!   send an error message and stop (although this should have been
!   caught as part of a preprocessor step).

! The interpolation can be done two different ways:
!   First, barnes_oi does the interpolation by processing all of
!     the available station data for each model grid cell.
!   Second, barnes_oi_ij does the interpolation by processing only
!     the "n_stns_used" number of stations for each model grid cell.
!   For small domains, with relatively few met stations (100's),
!   the first way is best.  For large domains (like the
!   United States, Globe, Pan-Arctic, North America, Greenland)
!   and many met stations (like 1000's), the second approach is the
!   most efficient.  But, the second approach carries the following
!   restrictions: 1) there can be no missing data for the fields of
!   interest; 2) there can be no missing stations (all stations
!   must exist throughout the simulation period); and 3) the
!   station met file must list the stations in the same order for
!   all time steps.  In addition, the code limits the number of
!   nearest stations used to be 5 or less.

! In this first case you are doing a Lagrangian sea ice parcel
!   simulation.
      if (seaice_run.eq.4.0) then
        do j=1,ny
          do i=1,nx
            grid(i,j) = var(i)
          enddo
        enddo
      else
        if (nstns.ge.2) then

!          call get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)  !KRA
          call get_dn(LIS_rc%gnc(1),LIS_rc%gnr(1),&
                      deltax,deltay,nstns,dn,iobsint)  !KRA

          if (barnes_lg_domain.eq.1.0 .or. barnes_lg_domain.eq.2.0) then

            if (barnes_lg_domain.eq.2.0) then
              call get_nearest_stns_2(nx,ny,xmn,ymn,deltax,deltay,&
     &          n_stns_used,k_stn,nstns,xstn,ystn,xg_line,yg_line,&
     &          snowmodel_line_flag)
            endif

            do j=1,ny
              do i=1,nx

! Use that nearest station list to extract the station information
!   to be used in the interpolation.
                do k=1,n_stns_used
                  xstn_tmp(k) = xstn(k_stn(i,j,k))
                  ystn_tmp(k) = ystn(k_stn(i,j,k))
                  var_tmp(k) = var(k_stn(i,j,k))
                enddo

! Do the interpolation for this model grid cell.
                call barnes_oi_ij(nx,ny,deltax,deltay,xmn,ymn,&
     &            n_stns_used,xstn_tmp,ystn_tmp,var_tmp,dn,grid,&
     &            undef,ifill,i,j,snowmodel_line_flag,xg_line,yg_line)

              enddo
            enddo

          else

            call barnes_oi(nx,ny,deltax,deltay,xmn,ymn,&
     &        nstns,xstn,ystn,var,dn,grid,undef,ifill)

          endif

        elseif (nstns.eq.1) then

          call single_stn(nx,ny,nstns,var,grid)

        else

          print *,'found no valid obs data at this time step'
          print *,'  model time =', iyear,imonth,iday,xhour
          stop

        endif

      endif

      return
      end subroutine interpolate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_nearest_stns_2(nx,ny,xmn,ymn,deltax,deltay,&
     &  n_stns_used,k_stn,nstns,xstn,ystn,xg_line,yg_line,&
     &  snowmodel_line_flag)

      use snowmodel_inc
      implicit none

      double precision xstn(nstns_max)
      double precision ystn(nstns_max)
      double precision dsq(nstns_max)
      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      double precision xg,yg,xmn,ymn,dist_min
      real deltax,deltay

      integer i,j,k,kk,nstns,n_stns_used,nx,ny
      integer k_stn(nx,ny,9)

      do j=1,ny
        do i=1,nx

! xcoords of grid nodes at index i,j
! ycoords of grid nodes at index i,j
          if (snowmodel_line_flag.eq.1.0) then
            xg = xg_line(i,j)
            yg = yg_line(i,j)
          elseif (snowmodel_line_flag.eq.0.0) then
            print *, 'This can be very slow.  I suggest you study'
            print *, 'and implement some version of the code below'
            print *, 'this subroutine before continuing along this'
            print *, 'simulation approach.'
            stop
            xg = xmn + deltax * (real(i) - 1.0)
            yg = ymn + deltay * (real(j) - 1.0)
          endif

! Loop through all of the stations, calculating the distance
!   between the current grid point and each of the stations.
          do k=1,nstns
            dsq(k) = (xg - xstn(k))**2 + (yg - ystn(k))**2
          enddo

! Loop through each of the station distances and find the
!   stations closest to the grid point in question.
          do kk=1,n_stns_used
            dist_min = 1.0e30
            do k=1,nstns
              if (dsq(k).le.dist_min) then
                k_stn(i,j,kk) = k
                dist_min = dsq(k)
              endif
            enddo

! Eliminate the last found minimum from the next search by making
!   its distance a big number.
            dsq(k_stn(i,j,kk)) = 1.0e30
          enddo

        enddo
      enddo

      return
      end subroutine get_nearest_stns_2

! get_nearest.f

! Find the ii,jj locations of the nearest observation grid
!   point (veg grid) closest to each of the output grid cells
!   (the topo grid).

! I'm starting out with original veg data on a UTM grid, and
!   interpolating to the UTM coords of the topo data grid.

! nnx, nny is the original veg data grid.
! nx, ny is the topo data grid.

!!    implicit none

!!    integer nx,ny,nstns,nnx,nny,n,icount

!!    parameter (nx=23401,ny=27502,nnx=14288,nny=16176)
!!    parameter (nstns=1)

!!    real dsq(nnx,nny)
!!    real xk_min(nx,ny,nstns)

!!    real del_x,del_y,gx_ll,gy_ll
!!    real xmn,ymn,deltax,deltay,xg,yg,xs,ys,dist_min
!!    integer i,j,k,kk,ii,jj,iii,iiii,jjj,jjjj,inc

!!    real gx_stn(nnx,nny)
!!    real gy_stn(nnx,nny)

!!    print *,'preforming the station indexing calculations'

! Grid increment for the veg grid, in m.
!!    del_x = 30.0
!!    del_y = 30.0

! Coordinates of center of lower left grid cell of veg (input) grid,
!   and grid increments.
!!    gx_ll = 398580.0 + del_x/2.0
!!    gy_ll = 8495405.0 + del_y/2.0

! UTM coordinates of the 'stations' or input grid-cell centers.
!!    do j=1,nny
!!      do i=1,nnx
!!        gx_stn(i,j) = gx_ll + del_x * (real(i) - 1.0)
!!        gy_stn(i,j) = gy_ll + del_y * (real(j) - 1.0)
!!      enddo
!!    enddo

! Grid increment for the topo grid, in m.
!!    deltax = 20.0
!!    deltay = 20.0

! Coordinates of center of lower left grid cell of topo (output) grid,
!   and grid increments.
!!    xmn = 381990.0 + deltax/2.0
!!    ymn = 8449790.0 + deltay/2.0

! Search for the stations nearest to the grid point of interest.
!!    icount = 0
!!    do j=1,ny
!!      if (mod(j,100).eq.0.0) print *, 'j = ',j
!!      do i=1,nx

! xcoords of topo grid nodes at index i,j
! ycoords of topo grid nodes at index i,j
!!        xg = xmn + deltax * (real(i) - 1.0)
!!        yg = ymn + deltay * (real(j) - 1.0)

! Find the corresponding veg grid cell coordinate to this grid cell.
!!        iii = nint((xg - gx_ll)/del_x + 1.0)
!!        jjj = nint((yg - gy_ll)/del_y + 1.0)

! Don't let things go outside the array bounds.
!!        inc = 5
!!        if (iii.le.inc) iii = 1 + inc
!!        if (iii.ge.nnx-inc+1) iii = nnx - inc
!!        if (jjj.le.inc) jjj = 1 + inc
!!        if (jjj.ge.nny-inc+1) jjj = nny - inc

! Loop through all of the stations, calculating the distance
!   between the current grid point and each of the stations.
!!        do jj=jjj-inc,jjj+inc
!!          do ii=iii-inc,iii+inc
!!            xs = gx_stn(ii,jj)
!!            ys = gy_stn(ii,jj)
!!            dsq(ii,jj) = (xg - xs)**2 + (yg - ys)**2
!!          enddo
!!        enddo

! Loop through each of the station distances and find the
!   stations closest to the grid point in question.
!!        do kk=1,nstns
!!          dist_min = 1.0e30
!!          do jj=jjj-inc,jjj+inc
!!            do ii=iii-inc,iii+inc
!!              if (dsq(ii,jj).le.dist_min) then
!!                k = ii + (jj - 1) * nnx
!!                xk_min(i,j,kk) = real(k)
!!                dist_min = dsq(ii,jj)
!!                iiii = ii
!!                jjjj = jj
!!              endif
!!            enddo
!!          enddo

! Eliminate the last found minimum from the next search by making
!   its distance a big number.
!!          dsq(iiii,jjjj) = 1.0e30
!!        enddo

!!      enddo
!!    enddo

! Save xk_min.
!!    print *,'--------------------------------------'
!!    print *,'--------------------------------------'
!!    print *, 'saving data'

!!    open (unit=41,file='nearest_stns_1.gdat',
!!   &  form='unformatted',access='direct',recl=4*nx)

!!    print *, '1'
!!    do j=1,ny
!!      write(41,rec=j) (xk_min(i,j,1),i=1,nx)
!!    enddo

!!    end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_good_values1(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var)

      use snowmodel_inc
      implicit none

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! input stn elevation
      real var(nstns_max)  ! input values
      real undef           ! undefined value
      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! input stn elevation
      real var_orig(nstns_max)  ! input values

      integer k

! Before performing the interpolation, sort through the data and
!   toss out any undefined values.
      nstns = 0
      do k=1,nstns_orig
        if (var_orig(k).ne.undef) then
          nstns = nstns + 1
          xstn(nstns) = xstn_orig(k)
          ystn(nstns) = ystn_orig(k)
          var(nstns) = var_orig(k)
          elev(nstns) = elev_orig(k)
        endif
      enddo

      return
      end subroutine get_good_values1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_good_values2(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,var1_orig,var2_orig,&
     &  var1,var2)

! Account for the special case where you must have two coincident
!   values to do the interpolation, like Tair and rh to interpolate
!   rh, and wind speed and dir to interpolate the winds.

      use snowmodel_inc
      implicit none

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! input stn elevation
      real undef           ! undefined value
      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! input stn elevation

      real var1_orig(nstns_max)  ! input values
      real var2_orig(nstns_max)  ! input values
      real var1(nstns_max)  ! input values
      real var2(nstns_max)  ! input values

      integer k

! Before performing the interpolation, sort through the data and
!   toss out any undefined values.
      nstns = 0
      do k=1,nstns_orig
        if (var1_orig(k).ne.undef .and. var2_orig(k).ne.undef) then
          nstns = nstns + 1
          xstn(nstns) = xstn_orig(k)
          ystn(nstns) = ystn_orig(k)
          var1(nstns) = var1_orig(k)
          var2(nstns) = var2_orig(k)
          elev(nstns) = elev_orig(k)
        endif
      enddo

      return
      end subroutine get_good_values2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_obs_data(nstns_orig,Tair_orig,rh_orig,xstn_orig,&
     &  ystn_orig,elev_orig,iyear,imonth,iday,xhour,undef,&
     &  windspd_orig,winddir_orig,prec_orig,isingle_stn_flag,&
     &  igrads_metfile,iter)

      use snowmodel_inc
      implicit none

      integer iyr,imo,idy      ! year, month, and day of data
      real xhr                 ! decimal hour
      integer idstn            ! station id number

      integer k,nstns_orig,isingle_stn_flag,igrads_metfile,iter
      integer iyear,imonth,iday

      real Tair_orig(nstns_max),rh_orig(nstns_max)
      real winddir_orig(nstns_max),windspd_orig(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)
      real elev_orig(nstns_max),xhour,prec_orig(nstns_max)
      real undef               ! undefined value

      if (isingle_stn_flag.eq.1) then
        nstns_orig = 1
      else
        read(20,*) nstns_orig
      endif

      if (nstns_orig.gt.nstns_max) then
        print *
        print *, 'You must increase the value of nstns_max in'
        print *, '  snowmodel.inc to be greater than nstns_orig.'
        print *, 'nstns_max = ',nstns_max
        print *, 'nstns_orig = ',nstns_orig
        print *
        stop
      endif

      do k=1,nstns_orig
        if (igrads_metfile.eq.1) then
          read(20,rec=iter) iyr,imo,idy,xhr,idstn,xstn_orig(k),&
     &      ystn_orig(k),elev_orig(k),Tair_orig(k),rh_orig(k),&
     &      windspd_orig(k),winddir_orig(k),prec_orig(k)
        else
          read(20,*) iyr,imo,idy,xhr,idstn,xstn_orig(k),&
     &      ystn_orig(k),elev_orig(k),Tair_orig(k),rh_orig(k),&
     &      windspd_orig(k),winddir_orig(k),prec_orig(k)
        endif

! MicroMet assumes that the air temperature comes in as deg C, but
!   all computations must be done in K.  Check for this and make
!   an appropriate adjustment.
        if (Tair_orig(k).lt.150.0 .and. Tair_orig(k).ne.undef) &
     &    Tair_orig(k) = Tair_orig(k) + 273.15

! Do some error checking.  Check for the presence of data from
!  the same date.
        if (iyr.ne.iyear .or. imo.ne.imonth .or. idy.ne.iday &
     &    .or. xhr.ne.xhour) then
          print *,'model time does not match data input time'
          print *,'  model =', iyear,imonth,iday,xhour
          print *,'  obs   =', iyr,imo,idy,xhr
          stop
        endif

      enddo

      return
      end subroutine get_obs_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_model_time(iyear_init,imonth_init,iday_init,&
     &  xhour_init,iter,dt,iyear,imonth,iday,xhour,J_day)

      implicit none

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour
      real dt                    ! model time step, in seconds
      integer iter               ! model iteration
      integer iyear_init         ! model start year
      integer imonth_init        ! model start month
      integer iday_init          ! model start day
      real xhour_init            ! model decimal start hour
      integer J_day              ! model day of year

! Misc. variables.
      real xmin,xhour_frac,xhour_tmp,xday
      integer ihour,last

      integer lastday(12)
      data lastday/31,28,31,30,31,30,31,31,30,31,30,31/

! Convert the simulation time to the exact year, month, day, and
!   decimal hour.

! Number of minutes into the simulation.  Here I have assumed that
!   dt will never be in fractions of minutes.
!     xmin = ((real(iter) - 1.0) * dt) / 60.0
      xmin = (real(iter) - 1.0) * (dt / 60.0)

! Model integration time in decimal hours.  The xhour_frac variable
!   needs to be fixed for dt < 3600 sec.
      xhour_tmp = xhour_init + xmin / 60.0
      ihour = mod(int(xhour_tmp),24)
      xhour_frac = 0.0
      xhour = real(ihour) + xhour_frac

! Number of days.
      xday = xhour_tmp / 24.0
      iday = iday_init + int(xday)

! Month and year, while accounting for leap-years.
      imonth = imonth_init
      iyear = iyear_init
 20   continue
      last = lastday(imonth)
      if (imonth.eq.2 .and. mod(iyear,4).eq.0 &
     &  .and. (mod(iyear,100).ne.0 .or. mod(iyear,1000).eq.0)) then
        last = last + 1
      endif
      if (iday.gt.last) then
        iday = iday - last
        imonth = imonth + 1
        if (imonth.gt.12) then
          imonth = 1
          iyear = iyear + 1
        endif
        go to 20
      endif

! Calculate the day of year (1...365,366) corresponding to the date
!   iyear-imonth-iday.
      J_day = iday &
     &  + min(1,max(0,imonth-1))*31 &
     &  + min(1,max(0,imonth-2))*(28+(1-min(1,mod(iyear,4)))) &
     &  + min(1,max(0,imonth-3))*31 &
     &  + min(1,max(0,imonth-4))*30 &
     &  + min(1,max(0,imonth-5))*31 &
     &  + min(1,max(0,imonth-6))*30 &
     &  + min(1,max(0,imonth-7))*31 &
     &  + min(1,max(0,imonth-8))*31 &
     &  + min(1,max(0,imonth-9))*30 &
     &  + min(1,max(0,imonth-10))*31 &
     &  + min(1,max(0,imonth-11))*30 &
     &  + min(1,max(0,imonth-12))*31

!     print *, iyear,imonth,iday,xhour,J_day

      return
      end subroutine get_model_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_lapse_rates(imonth,iday,T_lapse_rate,&
     &  Td_lapse_rate,xlat,lapse_rate_user_flag,&
     &  precip_lapse_rate,iprecip_lapse_rate_user_flag)

      implicit none

      integer imonth,iday,mbefore,mafter,k,&
     &  lapse_rate_user_flag,iprecip_lapse_rate_user_flag
      real weight,T_lapse_rate,Td_lapse_rate,A,B,C,xlat,&
     &  precip_lapse_rate
      integer, parameter :: months = 12
      integer lastday(months)
      data lastday/31,28,31,30,31,30,31,31,30,31,30,31/

! The lapse rate units are in deg_C km-1.  They are converted to
!   negative_deg_C m-1 below.
      real lapse_rate(months)
      real lapse_rate_nohem(months)
      real lapse_rate_sohem(months)
      real lapse_rate_user(months)
      data lapse_rate_nohem /4.4,5.9,7.1,7.8,8.1,8.2,&
     &                       8.1,8.1,7.7,6.8,5.5,4.7/
      data lapse_rate_sohem /8.1,8.1,7.7,6.8,5.5,4.7,&
     &                       4.4,5.9,7.1,7.8,8.1,8.2/

! If you want to use the 'user' array, put your monthly values in
!   here and set lapse_rate_user_flag = 1 in the .par file.
      data lapse_rate_user /4.4,5.9,7.1,7.8,8.1,8.2, &
     &                      8.1,8.1,7.7,6.8,5.5,4.7/

! The initial vapor pressure coeffs are in units of km-1.
      real am(months)
      real am_nohem(months)
      real am_sohem(months)
      real am_user(months)
      data am_nohem /0.41,0.42,0.40,0.39,0.38,0.36, &
     &               0.33,0.33,0.36,0.37,0.40,0.40/
      data am_sohem /0.33,0.33,0.36,0.37,0.40,0.40, &
     &               0.41,0.42,0.40,0.39,0.38,0.36/

! If you want to use the 'user' array, put your monthly values in
!   here and set lapse_rate_user_flag = 1 in the .par file.
      data am_user /0.41,0.42,0.40,0.39,0.38,0.36, &
     &              0.33,0.33,0.36,0.37,0.40,0.40/

! The precipitation lapse rate units are in km-1.
      real prec_lapse_rate(months)
      real precip_lapse_rate_nohem(months)
      real precip_lapse_rate_sohem(months)
      real precip_lapse_rate_user(months)
      data precip_lapse_rate_nohem /0.35,0.35,0.35,0.30,0.25,0.20,&
     &                              0.20,0.20,0.20,0.25,0.30,0.35/
      data precip_lapse_rate_sohem /0.20,0.20,0.20,0.25,0.30,0.35,&
     &                              0.35,0.35,0.35,0.30,0.25,0.20/

! If you want to use the 'user' array, put your monthly values in
!   here and set iprecip_lapse_rate_user_flag = 1 in the .par file.
      data precip_lapse_rate_user /0.35,0.35,0.35,0.30,0.25,0.20,&
     &                             0.20,0.20,0.20,0.25,0.30,0.35/

! Air and dewpoint temperature.
      do k=1,months
        if (lapse_rate_user_flag.eq.0) then
          if (xlat.lt.0.0) then
            lapse_rate(k) = lapse_rate_sohem(k)
            am(k) = am_sohem(k)
          else
            lapse_rate(k) = lapse_rate_nohem(k)
            am(k) = am_nohem(k)
          endif
        elseif (lapse_rate_user_flag.eq.1) then
          lapse_rate(k) = lapse_rate_user(k)
          am(k) = am_user(k)
        endif
      enddo

! Precipitation.
      do k=1,months
        if (iprecip_lapse_rate_user_flag.eq.0) then
          if (xlat.lt.0.0) then
            prec_lapse_rate(k) = precip_lapse_rate_sohem(k)
          else
            prec_lapse_rate(k) = precip_lapse_rate_nohem(k)
          endif
        elseif (iprecip_lapse_rate_user_flag.eq.1) then
          prec_lapse_rate(k) = precip_lapse_rate_user(k)
        endif
      enddo

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: temperatures for Buck`s equations are in deg C, and
!   vapor pressures are in mb.  Do the adjustments so that the
!   calculations are done with temperatures in K, and vapor
!   pressures in Pa.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

! Over ice.
!     A = 6.1115 * 100.0
!     B = 22.452
!     C = 272.55

! Find the month before and after the day in question.
      if (iday.le.15) then
        mbefore = imonth - 1
        if (mbefore.eq.0) mbefore = 12
        mafter = imonth
        weight = (real(lastday(mbefore)) - 15. + real(iday)) / &
     &    real(lastday(mbefore))
      else
        mbefore = imonth
        mafter = imonth + 1
        if (mafter.eq.13) mafter = 1
        weight = (real(iday) - 15.) / real(lastday(mbefore))
      endif

! Define the temperature lapse rate (deg C/m).
      T_lapse_rate = (- (weight * lapse_rate(mafter) + &
     &  (1. - weight) * lapse_rate(mbefore))) / 1000.0

! Define the dew-point temperature lapse rate (deg C/m).
      Td_lapse_rate = (- ((weight * am(mafter) + &
     &  (1. - weight) * am(mbefore)) * C)) / (B * 1000.0)

! Define the precipitation lapse rate (km-1).
      precip_lapse_rate = weight * prec_lapse_rate(mafter) + &
     &  (1. - weight) * prec_lapse_rate(mbefore)

      return
      end subroutine get_lapse_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)

      implicit none

      integer nx,ny,nstns
      real deltax,deltay,dn
      real dn_max           ! the max obs spacing, dn_r
      real dn_min           ! dn_r, for large n
      integer iobsint       ! flag (=1) use dn value from .par file

! Calculate an appropriate filtered wavelength value.  First
!   calculate dn for the case of severely nonuniform data, and
!   then for the case where there is just about a station for
!   every grid cell.  Then assume that the average of these two
!   is a reasonable value to use in the interpolation.
        dn_max = sqrt(deltax*real(nx) * deltay*real(ny)) * &
     &    ((1.0 + sqrt(real(nstns))) / (real(nstns) - 1.0))
        dn_min = sqrt((deltax*real(nx) * deltay*real(ny)) / &
     &    real(nstns))

        if (iobsint.eq.1) then
!         dn = dn
        else
          dn = 0.5 * (dn_min + dn_max)
        endif

!       print *,'You are using an average obs spacing of',dn
!       print *,'  the program indicates a min, max range of',
!    &    dn_min,dn_max

      return
      end subroutine get_dn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine barnes_oi(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,var,dn,grid,undef,ifill)

! This is an implementation of the Barnes objective analysis scheme
!   as described in:
!
!   Koch, S. E., M. DesJardins, and P. J. Kocin, 1983: An
!   interactive Barnes objective map analysis scheme for use with
!   satellite and conventional data. J. Climate and Applied
!   Meteorology, 22(9), 1487-1503.

      use snowmodel_inc
      implicit none

      real, parameter :: gamma = 0.2

      real pi

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn !center x coords of lower left grid cell
      double precision ymn !center y coords of lower left grid cell

      integer nstns        ! number of input values, all good
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real var(nstns_max)  ! input values
      integer nflag        ! determines if output will be undef value
      real undef           ! undefined value

      real dn                  ! average observation spacing
      real grid(nx,ny) ! output values

      integer i,j      ! col, row counters
      integer mm,nn    ! station counters
      integer ifill    ! flag (=1) forces a value in every cell

      double precision xg,yg !temporary x and y coords of current cell
      real w1,w2       ! weights for Gauss-weighted average
      real wtot1,wtot2 ! sum of weights
      real ftot1,ftot2 ! accumulators for values, corrections
      real dsq         ! delx**2 + dely**2
      double precision xa,ya       ! x, y coords of current station
      double precision xb,yb       ! x, y coords of current station
      real dvar(nstns_max)   ! estimated error

      real xkappa_1    ! Gauss constant for first pass
      real xkappa_2    ! Gauss constant for second pass
      real rmax_1      ! maximum scanning radii, for first
      real rmax_2      ! and second passes
      real anum_1      ! numerator, beyond scanning radius,
      real anum_2      ! for first and second passes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the first and second pass values of the scaling parameter
!   and the maximum scanning radius used by the Barnes scheme.
!   Values above this maximum will use a 1/r**2 weighting.  Here I
!   have assumed a gamma value of 0.2.

! First-round values, Eqn (13).
      pi = 2.0 * acos(0.0)
      xkappa_1 = 5.052 * (2.0*dn/pi)**2

! Define the maximum scanning radius to have weight defined by
!   wt = 1.0 x 10**(-30) = exp(-rmax_1/xkappa_1)
! Also scale the 1/r**2 wts so that when dsq = rmax, the wts match.
      rmax_1 = xkappa_1 * 30.0 * log(10.0)
      anum_1 = 1.0e-30 * rmax_1

! Second-round values, Eqn (4).
      xkappa_2 = gamma * xkappa_1
      rmax_2 = rmax_1 * gamma
      anum_2 = 1.0e-30 * rmax_2

! Scan each input data point and construct estimated error, dvar, at
!   that point.
      do nn=1,nstns

        xa = xstn(nn)
        ya = ystn(nn)
        wtot1 = 0.0
        ftot1 = 0.0

        do mm=1,nstns

          xb = xstn(mm)
          yb = ystn(mm)
          dsq = (xb - xa)**2 + (yb - ya)**2

          if (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)

          else

! Assume a 1/r**2 weight.
            w1 = anum_1/dsq

          endif

          wtot1 = wtot1 + w1
          ftot1 = ftot1 + w1 * var(mm)

        end do    ! end loop on sites m

        if (wtot1.eq.0.0) print *,'stn wt totals zero'

        dvar(nn) = var(nn) - ftot1/wtot1

      end do        ! end prediction loop on sites nn

! Grid-prediction loop.  Generate the estimate using first set of
!   weights, and correct using error estimates, dvar, and second
!   set of weights.

      do j=1,ny
      do i=1,nx

! xcoords of grid nodes at index i,j
! ycoords of grid nodes at index i,j
        xg = xmn + deltax * (real(i) - 1.0)
        yg = ymn + deltay * (real(j) - 1.0)

! Scan each input data point.
        ftot1 = 0.0
        wtot1 = 0.0
        ftot2 = 0.0
        wtot2 = 0.0
        nflag = 0

        do nn=1,nstns
           
          xa = xstn(nn)
          ya = ystn(nn)
          dsq = (xg - xa)**2 + (yg - ya)**2

          if (dsq.le.rmax_2) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = exp((- dsq)/xkappa_2)

          elseif (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = anum_2/dsq

          else

! Assume a 1/r**2 weight.
            w1 = anum_1/dsq
            nflag = nflag + 1
! With anum_2/dsq.
            w2 = gamma * w1

          endif

          wtot1 = wtot1 + w1
          wtot2 = wtot2 + w2
          ftot1 = ftot1 + w1 * var(nn)
          ftot2 = ftot2 + w2 * dvar(nn)
           
        end do    ! end loop on data sites nn

        if (wtot1.eq.0.0 .or. wtot2.eq.0.0) print *,'wts total zero'

        if (ifill.eq.1) then
          grid(i,j) = ftot1/wtot1 + ftot2/wtot2
        else
          if (nflag.lt.nstns) then
            grid(i,j) = ftot1/wtot1 + ftot2/wtot2
          else
            grid(i,j) = undef
          endif
        endif

      end do         ! end loop on cols i
      end do         ! end loop on rows j

      return
      end subroutine barnes_oi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine barnes_oi_ij(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns,xstn,ystn,var,dn,grid,&
     &  undef,ifill,i,j,snowmodel_line_flag,xg_line,yg_line)

! This is an implementation of the Barnes objective analysis scheme
!   as described in:
!
!   Koch, S. E., M. DesJardins, and P. J. Kocin, 1983: An
!   interactive Barnes objective map analysis scheme for use with
!   satellite and conventional data. J. Climate and Applied
!   Meteorology, 22(9), 1487-1503.

      use snowmodel_inc
      implicit none

      real, parameter ::  gamma = 0.2

      real pi

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn !center x coords of lower left grid cell
      double precision ymn !center y coords of lower left grid cell

      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      integer nstns        ! number of input values, all good
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real var(nstns_max)  ! input values
      integer nflag        ! determines if output will be undef value
      real undef           ! undefined value

      real dn                  ! average observation spacing
      real grid(nx,ny) ! output values

      integer i,j      ! col, row counters
      integer mm,nn    ! station counters
      integer ifill    ! flag (=1) forces a value in every cell

      double precision xg,yg !temporary x and y coords of current cell
      real w1,w2       ! weights for Gauss-weighted average
      real wtot1,wtot2 ! sum of weights
      real ftot1,ftot2 ! accumulators for values, corrections
      real dsq         ! delx**2 + dely**2
      double precision xa,ya       ! x, y coords of current station
      double precision xb,yb       ! x, y coords of current station
      real dvar(nstns_max)   ! estimated error

      real xkappa_1    ! Gauss constant for first pass
      real xkappa_2    ! Gauss constant for second pass
      real rmax_1      ! maximum scanning radii, for first
      real rmax_2      ! and second passes
      real anum_1      ! numerator, beyond scanning radius,
      real anum_2      ! for first and second passes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the first and second pass values of the scaling parameter
!   and the maximum scanning radius used by the Barnes scheme.
!   Values above this maximum will use a 1/r**2 weighting.  Here I
!   have assumed a gamma value of 0.2.

! First-round values, Eqn (13).
      pi = 2.0 * acos(0.0)
      xkappa_1 = 5.052 * (2.0*dn/pi)**2

! Define the maximum scanning radius to have weight defined by
!   wt = 1.0 x 10**(-30) = exp(-rmax_1/xkappa_1)
! Also scale the 1/r**2 wts so that when dsq = rmax, the wts match.
      rmax_1 = xkappa_1 * 30.0 * log(10.0)
      anum_1 = 1.0e-30 * rmax_1

! Second-round values, Eqn (4).
      xkappa_2 = gamma * xkappa_1
      rmax_2 = rmax_1 * gamma
      anum_2 = 1.0e-30 * rmax_2

! Scan each input data point and construct estimated error, dvar, at
!   that point.
      do nn=1,nstns

        xa = xstn(nn)
        ya = ystn(nn)
        wtot1 = 0.0
        ftot1 = 0.0

        do mm=1,nstns

          xb = xstn(mm)
          yb = ystn(mm)
          dsq = (xb - xa)**2 + (yb - ya)**2

          if (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)

          else

! Assume a 1/r**2 weight.
            w1 = anum_1/dsq

          endif

          wtot1 = wtot1 + w1
          ftot1 = ftot1 + w1 * var(mm)

        end do    ! end loop on sites m

        if (wtot1.eq.0.0) print *,'stn wt totals zero'

        dvar(nn) = var(nn) - ftot1/wtot1

      end do        ! end prediction loop on sites nn

! Grid-prediction loop.  Generate the estimate using first set of
!   weights, and correct using error estimates, dvar, and second
!   set of weights.

!     do 666 j=1,ny
!     do 555 i=1,nx

! xcoords of grid nodes at index i,j
! ycoords of grid nodes at index i,j
        if (snowmodel_line_flag.eq.1.0) then
          xg = xg_line(i,j)
          yg = yg_line(i,j)
        else
          xg = xmn + deltax * (real(i) - 1.0)
          yg = ymn + deltay * (real(j) - 1.0)
        endif

! Scan each input data point.
        ftot1 = 0.0
        wtot1 = 0.0
        ftot2 = 0.0
        wtot2 = 0.0
        nflag = 0

        do nn=1,nstns
           
          xa = xstn(nn)
          ya = ystn(nn)
          dsq = (xg - xa)**2 + (yg - ya)**2

          if (dsq.le.rmax_2) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = exp((- dsq)/xkappa_2)

          elseif (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = anum_2/dsq

          else

! Assume a 1/r**2 weight.
            w1 = anum_1/dsq
            nflag = nflag + 1
! With anum_2/dsq.
            w2 = gamma * w1

          endif

          wtot1 = wtot1 + w1
          wtot2 = wtot2 + w2
          ftot1 = ftot1 + w1 * var(nn)
          ftot2 = ftot2 + w2 * dvar(nn)
           
        end do   ! end loop on data sites nn

        if (wtot1.eq.0.0 .or. wtot2.eq.0.0) print *,'wts total zero'

        if (ifill.eq.1) then
          grid(i,j) = ftot1/wtot1 + ftot2/wtot2
        else
          if (nflag.lt.nstns) then
            grid(i,j) = ftot1/wtot1 + ftot2/wtot2
          else
            grid(i,j) = undef
          endif
        endif

!c 555 continue         ! end loop on cols i
!c 666 continue         ! end loop on rows j

      return
      end subroutine barnes_oi_ij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine single_stn(nx,ny,nstns,var,grid)

      use snowmodel_inc
      implicit none

      integer nstns    ! number of input values, all good
      integer nx       ! number of x output values
      integer ny       ! number of y output values

      real var(nstns_max)      ! input values
      real grid(nx,ny) ! output values
      integer i,j              ! col, row counters

!c Assign the station value to every grid cell.
      do j=1,ny
        do i=1,nx
          grid(i,j) = var(nstns)
        enddo
      enddo

      return
      end subroutine single_stn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine pressure(nx,ny,topo,sfc_pressure)

      use snowmodel_inc
      implicit none

      integer i,j,nx,ny

      real topo(nx,ny),sfc_pressure(nx,ny)
      real one_atmos,scale_ht

      one_atmos = 101300.0
!     scale_ht = 8500.0
      scale_ht = 8000.0

! Compute the average station pressure (in Pa).
      do j=1,ny
        do i=1,nx
          sfc_pressure(i,j) = one_atmos * exp((- topo(i,j))/scale_ht)
        enddo
      enddo

      return
      end subroutine pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine shortwave_data(nx,ny,deltax,deltay,xmn,ymn,&
     &  iyear,imonth,iday,xhour,undef,var_grid,iter)

! This program takes observations as discrete points, compares
!   those observations with a gridded model representation of those
!   observations, at the corresponding grid cells, computes a
!   difference between the observations and modeled values, fits
!   a gridded surface through those differences, and adds the
!   difference grid to the model grid.  Thus, correcting the model
!   outputs with the observations.

      use snowmodel_inc
      implicit none

      real deltax,deltay,xhour,xhr,elev(nstns_max),undef
      real var_grid(nx,ny),delta_var_grid(nx,ny)
      real var_obs(nstns_max),elev_orig(nstns_max),var_orig(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)

      integer nx,ny,k,nstns,iter,iyr,imo,idy,iyear,imonth,iday,&
     &  idstn_orig,nstns_orig,i,j

! Open the observation data file.
      if (iter.eq.1) open (unit=71,file='extra_met/shortwave.dat')

! Read the data describing the time, location, and variable values
!   for each station, at this time step.  Here I have assumed that
!   the data file is in the 'non-single-station' format (with a
!   station count listed at the begining at each new time step).
      read(71,*) nstns_orig
      do k=1,nstns_orig
        read(71,*) iyr,imo,idy,xhr,idstn_orig,xstn_orig(k),&
     &    ystn_orig(k),elev_orig(k),var_orig(k)
      enddo

! Compare the observation time with the model time.
      if (iyr.ne.iyear .or. imo.ne.imonth .or. idy.ne.iday &
     &  .or. xhr.ne.xhour) then
        print *,'model time does not match obs data input time'
        print *,'  model =', iyear,imonth,iday,xhour
        print *,'  obs   =', iyr,imo,idy,xhr
        stop
      endif

! Filter through the original input data, and eliminate any
!   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var_obs)

! If there are no observational data at this time step, use the
!   modeled values without any modification.  If there are some
!   good data, do the correction/data assimilation.
      if (nstns.gt.0) then
        call DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,&
     &    nstns,var_obs,delta_var_grid,var_grid)
      endif

! For incoming shortwave, incoming longwave, and surface pressure,
!   make sure no negetive numbers have been produced.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = max(0.0,var_grid(i,j))
        enddo
      enddo

      open (74,file='extra_met/shortwave_grid.gdat',&
     &  form='unformatted',access='direct',recl=4*nx*ny)
      write (74,rec=iter) ((delta_var_grid(i,j),i=1,nx),j=1,ny)
      close(74)

      return
      end subroutine shortwave_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine longwave_data(nx,ny,deltax,deltay,xmn,ymn,&
     &  iyear,imonth,iday,xhour,undef,var_grid,iter)

! This program takes observations as discrete points, compares
!   those observations with a gridded model representation of those
!   observations, at the corresponding grid cells, computes a
!   difference between the observations and modeled values, fits
!   a gridded surface through those differences, and adds the
!   difference grid to the model grid.  Thus, correcting the model
!   outputs with the observations.

      use snowmodel_inc
      implicit none

      real deltax,deltay,xhour,xhr,elev(nstns_max),undef
      real var_grid(nx,ny),delta_var_grid(nx,ny)
      real var_obs(nstns_max),elev_orig(nstns_max),var_orig(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)

      integer nx,ny,k,nstns,iter,iyr,imo,idy,iyear,imonth,iday,&
     &  idstn_orig,nstns_orig,i,j

! Open the observation data file.
      if (iter.eq.1) open (unit=72,file='extra_met/longwave.dat')

! Read the data describing the time, location, and variable values
!   for each station, at this time step.  Here I have assumed that
!   the data file is in the 'non-single-station' format (with a
!   station count listed at the begining at each new time step).
      read(72,*) nstns_orig
      do k=1,nstns_orig
        read(72,*) iyr,imo,idy,xhr,idstn_orig,xstn_orig(k),&
     &    ystn_orig(k),elev_orig(k),var_orig(k)
      enddo

! Compare the observation time with the model time.
      if (iyr.ne.iyear .or. imo.ne.imonth .or. idy.ne.iday &
     &  .or. xhr.ne.xhour) then
        print *,'model time does not match obs data input time'
        print *,'  model =', iyear,imonth,iday,xhour
        print *,'  obs   =', iyr,imo,idy,xhr
        stop
      endif

! Filter through the original input data, and eliminate any
!   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var_obs)

! If there are no observational data at this time step, use the
!   modeled values without any modification.  If there are some
!   good data, do the correction/data assimilation.
      if (nstns.gt.0) then
        call DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,&
     &    nstns,var_obs,delta_var_grid,var_grid)
      endif

! For incoming shortwave, incoming longwave, and surface pressure,
!   make sure no negetive numbers have been produced.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = max(0.0,var_grid(i,j))
        enddo
      enddo

!     open (75,file='extra_met/longwave_grid.gdat',
!    &  form='unformatted',access='direct',recl=4*nx*ny)
!     write (75,rec=iter) ((delta_var_grid(i,j),i=1,nx),j=1,ny)
!     close(75)

      return
      end subroutine longwave_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sfc_pressure_data(nx,ny,deltax,deltay,xmn,ymn,&
     &  iyear,imonth,iday,xhour,undef,var_grid,iter)

! This program takes observations as discrete points, compares
!   those observations with a gridded model representation of those
!   observations, at the corresponding grid cells, computes a
!   difference between the observations and modeled values, fits
!   a gridded surface through those differences, and adds the
!   difference grid to the model grid.  Thus, correcting the model
!   outputs with the observations.

      use snowmodel_inc
      implicit none

      real deltax,deltay,xhour,xhr,elev(nstns_max),undef
      real var_grid(nx,ny),delta_var_grid(nx,ny)
      real var_obs(nstns_max),elev_orig(nstns_max),var_orig(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)

      integer nx,ny,k,nstns,iter,iyr,imo,idy,iyear,imonth,iday,&
     &  idstn_orig,nstns_orig,i,j

! Open the observation data file.
      if (iter.eq.1) open (unit=73,file='extra_met/sfc_pressure.dat')

! Read the data describing the time, location, and variable values
!   for each station, at this time step.  Here I have assumed that
!   the data file is in the 'non-single-station' format (with a
!   station count listed at the begining at each new time step).
      read(73,*) nstns_orig
      do k=1,nstns_orig
        read(73,*) iyr,imo,idy,xhr,idstn_orig,xstn_orig(k),&
     &    ystn_orig(k),elev_orig(k),var_orig(k)
      enddo

! Compare the observation time with the model time.
      if (iyr.ne.iyear .or. imo.ne.imonth .or. idy.ne.iday &
     &  .or. xhr.ne.xhour) then
        print *,'model time does not match obs data input time'
        print *,'  model =', iyear,imonth,iday,xhour
        print *,'  obs   =', iyr,imo,idy,xhr
        stop
      endif

! Filter through the original input data, and eliminate any
!   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,&
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var_obs)

! If there are no observational data at this time step, use the
!   modeled values without any modification.  If there are some
!   good data, do the correction/data assimilation.
      if (nstns.gt.0) then
        call DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,&
     &    nstns,var_obs,delta_var_grid,var_grid)
      endif

! For incoming shortwave, incoming longwave, and surface pressure,
!   make sure no negetive numbers have been produced.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = max(0.0,var_grid(i,j))
        enddo
      enddo

!     open (76,file='extra_met/sfc_pressure_grid.gdat',
!    &  form='unformatted',access='direct',recl=4*nx*ny)
!     write (76,rec=iter) ((delta_var_grid(i,j),i=1,nx),j=1,ny)
!     close(76)

      return
      end subroutine sfc_pressure_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,&
     &  nstns,var_obs,delta_var_grid,var_grid)

      use snowmodel_inc
      implicit none

      real deltax,deltay,undef,dn
      real var_grid(nx,ny),delta_var_grid(nx,ny)
      real var_model(nstns_max),var_obs(nstns_max),&
     &  delta_var(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)

      integer ii(nstns_max),jj(nstns_max)
      integer nx,ny,i,j,ifill,iobsint,k,nstns

! Convert the x and y locations to (ii,jj) locations.
      do k=1,nstns
        ii(k) = 1 + nint((xstn(k) - xmn) / deltax)
        jj(k) = 1 + nint((ystn(k) - ymn) / deltay)
      enddo

! Extract the modeled data at the appropriate grid cells.
      do k=1,nstns
        var_model(k) = var_grid(ii(k),jj(k))
      enddo

! Calculate the difference between the modeled variable and the
!   observation at each point/grid cell.
      do k=1,nstns
        delta_var(k) = var_obs(k) - var_model(k)
      enddo

! Now that I have the differences calculated at each observation
!   point, interpolate them over the simulation domain.  Use the
!   barnes oi scheme to create the distribution. If there is
!   only a single station, distribute those data uniformly over
!   the domain.  Make sure that ifill=1, and then undef is not
!   really used (so it does not have to be the same as defined in
!   the .par file).
      undef = -9999.0
      ifill = 1
      iobsint = 0

      if (nstns.ge.2) then
        call get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)
        call barnes_oi(nx,ny,deltax,deltay,xmn,ymn,&
     &    nstns,xstn,ystn,delta_var,dn,delta_var_grid,undef,ifill)
      elseif (nstns.eq.1) then
        call single_stn(nx,ny,nstns,delta_var,delta_var_grid)
      endif

! Use the gridded delta surface to correct the modeled variable.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = var_grid(i,j) + delta_var_grid(i,j)
        enddo
      enddo

      return
      end subroutine DATA_ASSIM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_lai(J_day,forest_LAI)

      implicit none

      integer J_day,n
      integer, parameter :: nftypes = 5

      real vlai_summer(nftypes),vlai_winter(nftypes),&
     &  forest_LAI(nftypes)

      real pi,daysinyr,tmax,tmin,peak_jday,dtseason,vtrans,tseason,&
     &  fseason

! The five forest types in MicroMet/SnowModel:
!  1  coniferous forest
!  2  deciduous forest
!  3  mixed forest
!  4  scattered short-conifer
!  5  clearcut conifer

      data vlai_summer /2.5, 2.5, 2.5, 1.5, 1.0/
      data vlai_winter /2.5, 0.5, 1.5, 1.5, 1.0/

! Note: A maximum forest LAI of 5.0 will give almost zero (like
!   10 W m^2) incoming solar under the canopy.  Values for Fraser
!   Experimental Forest in Colorado are 2-3 (LSOS site = 1.8,
!   Kelly's/Gus' site = 2.3).

! Calculate a seasonally varying temperature, assuming a max and
!   min temperature and a cos distribution peaking in mid July
!   (J_day = 200).  Then use this to define the seasonal lai
!   variation.
      pi = 2.0 * acos(0.0)
      daysinyr = 366.0
      tmax = 298.0
      tmin = 273.0
      peak_jday = 200.0

      dtseason = tmax - tmin
      vtrans = tmin + dtseason / 2.0

      tseason = vtrans + dtseason / 2.0 * &
     &  cos(2.0 * pi / daysinyr * (real(J_day) - peak_jday))

      fseason = 0.0016 * (tmax - tseason)**2

      do n=1,nftypes
        forest_LAI(n) = (1.0 - fseason) * vlai_summer(n) + &
     &    fseason * vlai_winter(n)
      enddo

!     print *,J_day,(forest_LAI(n),n=1,nftypes)

      return
      end subroutine get_lai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_wind_file(nx,ny,iter,uwind_grid,vwind_grid,&
     &  windspd_grid,winddir_grid,windspd_flag,winddir_flag,&
     &  windspd_min)

      use snowmodel_inc
      implicit none

      integer nx
      integer ny

      real uwind_grid(nx,ny)
      real vwind_grid(nx,ny)
      real winddir_grid(nx,ny)
      real windspd_grid(nx,ny)

      real pi,deg2rad,rad2deg
      integer i,j,iter,irec
      real windspd_flag,winddir_flag,u_sum,v_sum,windspd_min

! Define the required constants.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0
      rad2deg = 180.0 / pi

! Read in the u and v arrays for this time step.
      if (iter.eq.1) then
        open (66,file='/data10/baffin/nuatmos/wind_spd_dir.gdat',&
     &    form='unformatted',access='direct',recl=4*nx*ny)
      endif

      irec = (iter - 1) * 2
      read (66,rec=irec+1) ((uwind_grid(i,j),i=1,nx),j=1,ny)
      read (66,rec=irec+2) ((vwind_grid(i,j),i=1,nx),j=1,ny)

! Convert these u and v components to speed and directions.
      do j=1,ny
        do i=1,nx

! Some compilers do not allow both u and v to be 0.0 in
!   the atan2 computation.
          if (abs(uwind_grid(i,j)).lt.1e-10) uwind_grid(i,j) = 1e-10

          winddir_grid(i,j) = rad2deg * &
     &      atan2(uwind_grid(i,j),vwind_grid(i,j))
          if (winddir_grid(i,j).ge.180.0) then
            winddir_grid(i,j) = winddir_grid(i,j) - 180.0
          else
            winddir_grid(i,j) = winddir_grid(i,j) + 180.0
          endif
          windspd_grid(i,j) = &
     &      sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
        enddo
      enddo

! Avoid problems of zero (low) winds (for example, turbulence
!   theory, log wind profile, etc., says that we must have some
!   wind.  Thus, some equations blow up when the wind speed gets
!   very small).
      do j=1,ny
        do i=1,nx
          if (windspd_grid(i,j).lt.windspd_min) then
            windspd_grid(i,j) = windspd_min
            uwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &        sin(deg2rad*winddir_grid(i,j))
            vwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &        cos(deg2rad*winddir_grid(i,j))
          endif
        enddo
      enddo

! Find the maximum wind speed in the domain, and the
!   domain-averaged wind direction.
      windspd_flag = 0.0
      u_sum = 0.0
      v_sum = 0.0
      do j=1,ny
        do i=1,nx
          windspd_flag = max(windspd_flag,windspd_grid(i,j))
          u_sum = u_sum + uwind_grid(i,j)
          v_sum = v_sum + vwind_grid(i,j)
        enddo
      enddo
      u_sum = u_sum / real(nx*ny)
      v_sum = v_sum / real(nx*ny)

! Some compilers do not allow both u and v to be 0.0 in
!   the atan2 computation.
      if (abs(u_sum).lt.1e-10) u_sum = 1e-10

      winddir_flag = rad2deg * atan2(u_sum,v_sum)
      if (winddir_flag.ge.180.0) then
        winddir_flag = winddir_flag - 180.0
      else
        winddir_flag = winddir_flag + 180.0
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_daily_irec (iter,dt,irec_day)

      implicit none

      integer iter,irec_day
      real dt,secs_in_day,secs_in_sim

      secs_in_day = 60.0 * 60.0 * 24.0

      secs_in_sim = dt * real(iter - 1)
      irec_day = int(secs_in_sim / secs_in_day) + 1

      return
      end subroutine get_daily_irec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
