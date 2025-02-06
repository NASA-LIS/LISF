! snowpack_code.f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SNOWPACK_CODE(nx,ny,Tair_grid,rh_grid,ro_nsnow,&
     &  dt,swe_depth,Tsfc,snow_d,prec_grid,runoff,Qm,rain,&
     &  sprec,iter,w_balance,sum_prec,sum_runoff,xro_snow,&
     &  undef,ro_snow,ro_snow_grid,soft_snow_d,sum_sprec,&
     &  snow_depth,windspd_grid,Qsi_grid,sum_Qcs,canopy_int,&
     &  Qcs,vegtype,forest_LAI,albedo,glacier_melt,&
     &  canopy_unload,sum_unload,sum_glacmelt,run_snowtran,&
     &  swemelt,d_canopy_int,sum_d_canopy_int,snow_d_init,&
     &  sfc_pressure,Qe,sfc_sublim_flag,sum_sfcsublim,&
     &  sum_swemelt,corr_factor,icorr_factor_index,swesublim,&
     &  swe_depth_old,canopy_int_old,KK,max_layers,melt_flag,&
     &  ro_snowmax,tsls_threshold,dz_snow_min,tslsnowfall,&
     &  change_layer,snod_layer,swed_layer,ro_layer,T_old,gamma,&
     &  multilayer_snowpack,seaice_run,seaice_conc,ht_windobs,&
     &  windspd_2m_grid,diam_layer,flux_layer,sum_trans)

      use snowmodel_inc
      implicit none

      integer nx,ny,iter,i,j

      integer max_layers,multilayer_snowpack,k,n_tsteps_in_day,irec
      integer KK(nx,ny)
      integer melt_flag(nx,ny,nz_max)

      real ro_snowmax,tsls_threshold,dz_snow_min,Cp_snow
      real tslsnowfall(nx,ny)
      real change_layer(nx,ny)
      real snod_layer(nx,ny,nz_max)
      real swed_layer(nx,ny,nz_max)
      real ro_layer(nx,ny,nz_max)
      real T_old(nx,ny,nz_max)
      real gamma(nx,ny,nz_max)
      real diam_layer(nx,ny,nz_max)
      real flux_layer(nx,ny,nz_max)

      integer melt_flag_z(nz_max)
      real snod_layer_z(nz_max)
      real swed_layer_z(nz_max)
      real ro_layer_z(nz_max)
      real T_old_z(nz_max)
      real gamma_z(nz_max)
      real diam_z(nz_max)
      real flux_z(nz_max)

      real Tair_grid(nx,ny)
      real rh_grid(nx,ny)
      real prec_grid(nx,ny)
      real windspd_grid(nx,ny)
      real windspd_2m_grid(nx,ny)
      real Qsi_grid(nx,ny)
      real vegtype(nx,ny)
      real albedo(nx,ny)
      real glacier_melt(nx,ny)
      real canopy_unload(nx,ny)
      real sum_unload(nx,ny)
      real sum_glacmelt(nx,ny)
      real sum_swemelt(nx,ny)
      real swemelt(nx,ny)
      real swesublim(nx,ny)
      real snow_d_init(nx,ny)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)
      real seaice_conc(nx,ny)
      real sum_trans(nx,ny)
      real ro_snow

      real, dimension(nx,ny) :: ro_nsnow,snow_d,&
     &  runoff,rain, &
     &  sprec,w_balance,&
     &  sum_prec,sum_runoff,&
     &  xro_snow,sfc_pressure,&
     &  ro_snow_grid,swe_depth,&
     &  Tsfc,Qm,&
     &  soft_snow_d,sum_sprec,&
     &  snow_depth,sum_Qcs,&
     &  canopy_int,Qcs,&
     &  d_canopy_int,sum_d_canopy_int,&
     &  Qe,sum_sfcsublim

      real dt,undef,Cp,xLf,Tf,A1,A2,ro_water,xLs,ro_ice,Twb,&
     &  run_snowtran,sfc_sublim_flag,seaice_run,ht_windobs

      real corr_factor(nx_max,ny_max,max_obs_dates+1)
      real corr_factor_ij
      integer icorr_factor_index(max_time_steps)

      integer, parameter :: nftypes = 5
      real forest_LAI(nftypes)

!     print *,'   solving the snow-cover evolution'

! Define the constants used in the computations.
      CALL CONSTS_SNOWPACK(Cp,xLs,ro_ice,xLf,Tf,A1,A2,ro_water,&
     &  Cp_snow,ro_snowmax)

! Run the snowpack evolution sub-model.
      do j=1,ny
        do i=1,nx

! Extract the vertical column for this i,j point, and send it
!   to the subroutine. *** Note that I should use f95, then I would
!   not have to do this (I could pass in subsections of the arrays).
          if (multilayer_snowpack.eq.1) then
            do k=1,nz_max
              melt_flag_z(k) = melt_flag(i,j,k)
              snod_layer_z(k) = snod_layer(i,j,k)
              swed_layer_z(k) = swed_layer(i,j,k)
              ro_layer_z(k) = ro_layer(i,j,k)
              T_old_z(k) = T_old(i,j,k)
              gamma_z(k) = gamma(i,j,k)
              diam_z(k) = diam_layer(i,j,k)
            enddo
          endif

! Extract the correction factor from the data assimilation array
!   so it can be passed into the SnowPack routines without the
!   negative array index.
          if (icorr_factor_index(iter).lt.0) then
            k = -icorr_factor_index(iter)
            corr_factor_ij = corr_factor(i,j,k)
          else
            corr_factor_ij = undef
          endif

          CALL SNOWPACK_CORE(Twb,Tf,Tair_grid(i,j),rh_grid(i,j),xLs,&
     &      Cp,sfc_pressure(i,j),ro_nsnow(i,j),dt,ro_snow,&
     &      swe_depth(i,j),Tsfc(i,j),A1,A2,snow_d(i,j),ro_water,&
     &      ro_ice,prec_grid(i,j),runoff(i,j),Qm(i,j),xLf,rain(i,j),&
     &      sprec(i,j),iter,w_balance(i,j),sum_prec(i,j),&
     &      sum_runoff(i,j),xro_snow(i,j),undef,&
     &      soft_snow_d(i,j),sum_sprec(i,j),ro_snow_grid(i,j),&
     &      snow_depth(i,j),windspd_grid(i,j),Qsi_grid(i,j),&
     &      sum_Qcs(i,j),canopy_int(i,j),Qcs(i,j),vegtype(i,j),&
     &      forest_LAI,albedo(i,j),canopy_unload(i,j),&
     &      sum_unload(i,j),sum_glacmelt(i,j),run_snowtran,&
     &      swemelt(i,j),d_canopy_int(i,j),sum_d_canopy_int(i,j),&
     &      snow_d_init(i,j),Qe(i,j),glacier_melt(i,j),&
     &      sfc_sublim_flag,sum_sfcsublim(i,j),sum_swemelt(i,j),&
     &      corr_factor_ij,icorr_factor_index(iter),swesublim(i,j),&
     &      swe_depth_old(i,j),canopy_int_old(i,j),KK(i,j),&
     &      max_layers,melt_flag_z,ro_snowmax,tsls_threshold,&
     &      dz_snow_min,tslsnowfall(i,j),change_layer(i,j),snod_layer_z,&
     &      swed_layer_z,ro_layer_z,T_old_z,gamma_z,multilayer_snowpack,&
     &      Cp_snow,seaice_run,ht_windobs,windspd_2m_grid(i,j),&
     &      diam_z,flux_z)

! Re-build the 3-D arrays.  See note above about using f95 to avoid this.
          if (multilayer_snowpack.eq.1) then
            do k=1,nz_max
              melt_flag(i,j,k) = melt_flag_z(k)
              snod_layer(i,j,k) = snod_layer_z(k)
              swed_layer(i,j,k) = swed_layer_z(k)
              ro_layer(i,j,k) = ro_layer_z(k)
              T_old(i,j,k) = T_old_z(k)
              gamma(i,j,k) = gamma_z(k)
              diam_layer(i,j,k) = diam_z(k)
              flux_layer(i,j,k) = flux_z(k)
            enddo
          endif

        enddo
      enddo

      if (run_snowtran.eq.0.0) then
        do j=1,ny
          do i=1,nx
          swe_depth_old(i,j) = swe_depth(i,j)
          canopy_int_old(i,j) = canopy_int(i,j)
          enddo
        enddo
      endif

! Read in the sea ice concentration.  These are daily data, so
!   first calculate which record in the data file this time step
!   corresponds to.
      if (seaice_run.ne.0.0) then
        n_tsteps_in_day = nint(86400.0 / dt)
        if (mod(iter-1,n_tsteps_in_day).eq.0) then
          irec = int((real(iter) - 0.5) * dt / 86400.0) + 1
          print *,'sea ice irec =',irec
          read (445,rec=irec) ((seaice_conc(i,j),i=1,nx),j=1,ny)
        endif
      endif

! If this simulation is not running SnowTran-3D, then zero out
!   the ocean grid cells that have no sea ice here.  If it is
!   running with SnowTran-3D, then do this in the SnowTran-3D
!   subroutine.
      if (run_snowtran.eq.0.0) then
        if (seaice_run.ne.0.0) then
          CALL ZERO_SEAICE_SNOW(nx,ny,snow_depth,ro_snow_grid,&
     &      ro_snow,swe_depth,swe_depth_old,canopy_int_old,KK,&
     &      tslsnowfall,snod_layer,swed_layer,ro_layer,T_old,&
     &      multilayer_snowpack,tsls_threshold,seaice_conc,&
     &      sum_sprec,sum_trans)
        endif
      endif

      return
      end SUBROUTINE SNOWPACK_CODE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SNOWPACK_CORE(Twb,Tf,Tair,rh,xLs,&
     &  Cp,sfc_pressure,ro_nsnow,dt,ro_snow,&
     &  swe_depth,Tsfc,A1,A2,snow_d,ro_water,&
     &  ro_ice,prec,runoff,Qm,xLf,rain,&
     &  sprec,iter,w_balance,sum_prec,&
     &  sum_runoff,xro_snow,undef,&
     &  soft_snow_d,sum_sprec,ro_snow_grid,&
     &  snow_depth,windspd,Qsi,&
     &  sum_Qcs,canopy_int,Qcs,vegtype,&
     &  forest_LAI,albedo,canopy_unload,&
     &  sum_unload,sum_glacmelt,run_snowtran,&
     &  swemelt,d_canopy_int,sum_d_canopy_int,&
     &  snow_d_init,Qe,glacier_melt,&
     &  sfc_sublim_flag,sum_sfcsublim,sum_swemelt,&
     &  corr_factor,icorr_factor_index,swesublim,&
     &  swe_depth_old,canopy_int_old,KK,&
     &  max_layers,melt_flag,ro_snowmax,tsls_threshold,&
     &  dz_snow_min,tslsnowfall,change_layer,snod_layer,&
     &  swed_layer,ro_layer,T_old,gamma,multilayer_snowpack,&
     &  Cp_snow,seaice_run,ht_windobs,windspd_2m,&
     &  diam_layer,flux_layer)

      use snowmodel_inc
      implicit none

      integer iter,icorr_factor_index

      integer KK,max_layers,multilayer_snowpack
      
      real ro_snowmax,tsls_threshold,dz_snow_min,tslsnowfall,Cp_snow

      integer melt_flag(nz_max)
      real change_layer
      real snod_layer(nz_max)
      real swed_layer(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      real gamma(nz_max)
      real diam_layer(nz_max)
      real flux_layer(nz_max)

      real Twb,Tf,Tair,rh,xLs,Cp,ro_nsnow,dt,ro_snow,swe_depth,&
     &  Tsfc,A1,A2,snow_d,ro_water,ro_ice,prec,runoff,Qm,xLf,rain,&
     &  sprec,w_balance,sum_prec,sum_runoff,xro_snow,undef,&
     &  soft_snow_d,sum_sprec,ro_snow_grid,snow_depth,sprec_grnd,&
     &  windspd,Qsi,sum_Qcs,canopy_int,Qcs,canopy_unload,&
     &  vegtype,albedo,glacier_melt,sum_unload,sum_glacmelt,&
     &  run_snowtran,swemelt,d_canopy_int,sfc_pressure,&
     &  sum_d_canopy_int,snow_d_init,Qe,sfc_sublim_flag,&
     &  sum_sfcsublim,sum_swemelt,corr_factor,swesublim,&
     &  swe_depth_old,canopy_int_old,sprec_grnd_ml,seaice_run,&
     &  ro_nsnow_wind,ht_windobs,windspd_2m

      integer,parameter ::  nftypes = 5
      real forest_LAI(nftypes)

! Calculate the canopy sublimation, loading and unloading.  Note
!   that here I have assumed that trees are type 1-5.
      if (vegtype.le.5.0) then
        CALL CANOPY_SNOW(rh,Tair,windspd,Qsi,sum_Qcs,albedo,&
     &    canopy_int,sprec,Qcs,dt,canopy_unload,&
     &    forest_LAI(nint(vegtype)),sum_unload,d_canopy_int,&
     &    sum_d_canopy_int)
        sprec_grnd = sprec + canopy_unload - d_canopy_int
        sprec_grnd_ml = sprec - d_canopy_int
      else
        Qcs = 0.0
        sprec_grnd = sprec
        sprec_grnd_ml = sprec
      endif

! Calculate the wind speed at 2 meters.
      CALL WINDSPEED_2M(windspd,ht_windobs,windspd_2m)

! Solve for the wet bulb temperature.
      CALL SOLVEWB(Twb,Tf,Tair,rh,xLs,Cp,sfc_pressure)

! Compute the new snow density.
      CALL NSNOWDEN(ro_nsnow,Twb,Tf,dt)

! Call the subroutine that increases the snow density due to
!   blowing snow wind speeds during a snowfall event.  Here we
!   are only calculating a increment to the previous NSNOWDEN
!   new snow density calculation.
      CALL NSNOW_DENSITY_FROM_BLOWING_SNOW(windspd_2m,sprec,dt,&
     &  ro_nsnow_wind)

! Update the new snow density with the wind contribution.
      ro_nsnow = ro_nsnow + ro_nsnow_wind

! Make sure the snow density falls within reasonable limits.
      ro_nsnow = min(ro_nsnow,ro_snowmax)

! Call the multi-layer snowpack model.
      if (multilayer_snowpack.eq.1) then

        CALL MULTI_LAYER_SNOW(KK,ro_layer,Tf,dt,ro_water,&
     &    ro_ice,T_old,snod_layer,swed_layer,Qm,ro_snowmax,rain,&
     &    xLf,Cp_snow,melt_flag,runoff,tslsnowfall,ro_nsnow,&
     &    sprec,Tsfc,tsls_threshold,gamma,max_layers,change_layer,&
     &    dz_snow_min,snow_depth,swe_depth,undef,canopy_unload,&
     &    vegtype,glacier_melt,sum_glacmelt,sum_swemelt,snow_d,&
     &    Qe,sfc_sublim_flag,sum_sfcsublim,soft_snow_d,ro_snow,&
     &    sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,&
     &    ro_snow_grid,xro_snow,swesublim,A1,A2,windspd_2m,&
     &    sfc_pressure,diam_layer,flux_layer,corr_factor,&
     &    icorr_factor_index)

! Call the original single-layer snowpack model.
      else

! Compute the snow density change due to settling.
        CALL DDENSITY(ro_snow_grid,swe_depth,Tf,Tsfc,dt,A1,A2,&
     &    snow_depth,ro_water,ro_snowmax)

! Compute the melt, rain, and snow contributions to modifying
!   the snowpack depth, density, and snow water equivalent.
        CALL SNOWPACK(swe_depth,snow_d,ro_snow_grid,&
     &    prec,ro_water,ro_nsnow,runoff,Qm,xLf,dt,rain,sprec,&
     &    sum_prec,sum_runoff,soft_snow_d,sum_sprec,ro_snow,&
     &    snow_depth,sprec_grnd,vegtype,glacier_melt,sum_glacmelt,&
     &    swemelt,canopy_unload,Qe,sfc_sublim_flag,sum_sfcsublim,&
     &    sum_swemelt,corr_factor,icorr_factor_index,swesublim,&
     &    ro_snowmax)

! Post process the data for output.
        CALL POSTPROC(ro_snow_grid,xro_snow,snow_depth,undef)

      endif

! Perform a water balance check (see notes in this subroutine).
      if (seaice_run.eq.0.0) then
        if (run_snowtran.eq.0.0) then
          CALL WATERBAL_SNOWPACK(w_balance,prec,Qcs,runoff,&
     &    d_canopy_int,swe_depth,glacier_melt,swe_depth_old,iter,&
     &    swesublim,canopy_unload,canopy_int_old,canopy_int)
        endif
      endif

      return
      end SUBROUTINE SNOWPACK_CORE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CANOPY_SNOW(rh,Tair,windspd,Qsi,sum_Qcs,albedo,&
     &  canopy_int,sprec,Qcs,dt,canopy_unload,&
     &  forest_LAI,sum_unload,d_canopy_int,&
     &  sum_d_canopy_int)

      implicit none

      real rh,Tair,windspd,V_s,Qsi,forest_LAI,dt,xImax,canopy_int,&
     &  d_canopy_int,Qcs,Ce,sprec,C_0,unload_melt,canopy_unload,&
     &  sum_Qcs,albedo,sum_unload,sum_d_canopy_int

! Note that all of this must deal with the (kg/m2)=(mm), => (m)
!   issues.  Precip is in (m), all of these equations are in
!   (kg/m2), and I want the outputs to be in (m).

! Compute the sublimation loss rate coefficient for canopy snow.
      CALL SUBLIM_COEF(rh,Tair,windspd,V_s,Qsi,albedo)

! Maximum interception storage.
      xImax = 4.4 * forest_LAI

! Change in canopy load due to snow precipitation during this time
!   step.  Convert the canopy interception to mm.
      canopy_int = 1000.0 * canopy_int
      d_canopy_int = 0.7 * (xImax - canopy_int) * &
     &  ( 1.0 - exp((- sprec)*1000.0/xImax))

! Update the interception load.
      canopy_int = canopy_int + d_canopy_int

! Canopy exposure coefficient.
      if (canopy_int.eq.0.0) then
        Ce = 0.0
      else
! Pomeroy's k_c value
!       Ce = 0.0114 * (canopy_int/xImax)**(-0.4)
! My k_c value.
        Ce = 0.00995 * (canopy_int/xImax)**(-0.4)
      endif

! Canopy sublimation (kg/m2), (a negative mumber).  Make sure that
!   you don't sublimate more than is available.
      Qcs = Ce * canopy_int * V_s * dt
      Qcs = -min(canopy_int,-Qcs)

! Remove the sublimated moisture from the canopy store.
      canopy_int = canopy_int + Qcs

! Save the sublimation in (m).
      Qcs = Qcs / 1000.0
      sum_Qcs = sum_Qcs + Qcs

! Perform a second unloading due to melt.  Assume an unloading rate
!   of 5.0 mm/day/C.
      C_0 = 5.0 / 86400.0
      unload_melt = C_0 * max(0.0,Tair-273.15) * dt
      unload_melt = min(canopy_int,unload_melt)
      canopy_int = canopy_int - unload_melt

! Keep track of the unloaded snow that reached the ground during
!   this time step (m) (this will add to the snow depth).
      canopy_unload = unload_melt / 1000.0
      d_canopy_int = d_canopy_int / 1000.0

! Save a summing array of this unloaded snow.
      sum_unload = sum_unload + canopy_unload

! Save a summing array of the change in canopy load.
      sum_d_canopy_int = sum_d_canopy_int + d_canopy_int

! Save the interception load for the next time step.  Convert to m.
      canopy_int = canopy_int / 1000.0

      return
      end SUBROUTINE CANOPY_SNOW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SUBLIM_COEF(rh,Tair,windspd,V_s,Qsi,albedo)

! Compute the sublimation loss rate coefficient for canopy snow.

      implicit none

      real pi,ro_ice,xM,R,R_dryair,vonKarman,visc_air,h_s,xlamdaT,&
     &  D,ro_sat,sigma,V_s,radius,xmass,windspd,rh,Tair,Qsi,Sp,&
     &  xN_r,xNu,xSh,top,bottom,omega,albedo

! Constants.
      pi = 2.0 * acos(0.0)
      ro_ice = 917.0
      xM = 18.01
      R = 8313.
      R_dryair = 287.
      vonKarman = 0.4
      visc_air = 13.e-6
      h_s = 2.838e6
      xlamdaT = 0.024

! Particle radius.
      radius = 5.0e-4

! Particle mass.
      xmass = 4.0/3.0 * pi * ro_ice * radius**3

! Diffusivity of water vapor in the atmosphere.
      D = 2.06e-5 * (Tair/273.0)**(1.75)

! Saturation density of water vapor.
      ro_sat = 0.622 / (R_dryair * Tair) * &
     &  611.15 * exp(22.452 * (Tair - 273.15) / (Tair - 0.61))

! Humidity deficit.
      sigma = 0.01 * rh - 1.0
      sigma = min(0.0,sigma)
      sigma = max(-1.0,sigma)

! Reynolds, Nusselt, and Sherwood numbers.
      xN_r = 2.0 * radius * windspd / visc_air
      xNu = 1.79 + 0.606 * xN_r**(0.5)
      xSh = xNu

! Solar radiation absorbed by the snow particle.  Here assume that
!   the general snow albedo is the same as the snow particle albedo.
      Sp = pi * radius**2 * (1.0 - albedo) * Qsi

! Sublimation-loss rate coefficient for an ice sphere.
      omega = ((h_s * xM)/(R * Tair) - 1.0) / (xlamdaT * Tair * xNu)
      top = 2.0 * pi * radius * sigma - Sp * omega
      bottom = h_s * omega + 1.0/(D * ro_sat * xSh)
      V_s = (top/bottom)/xmass

      return
      end SUBROUTINE SUBLIM_COEF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE WATERBAL_SNOWPACK(w_balance,prec,Qcs,runoff,&
     &  d_canopy_int,swe_depth,glacier_melt,swe_depth_old,iter,&
     &  swesublim,canopy_unload,canopy_int_old,canopy_int)

      implicit none

      integer iter

      real w_balance,prec,Qcs,runoff,d_canopy_int,swe_depth_old,&
     &  swe_depth,glacier_melt,swesublim,canopy_unload,canopy_int_old,&
     &  canopy_int

! Note that the following balances should hold.  These aren't quite
!   right, but it is a place to start.
!   Canopy Balance (forest):
!     canopy = sprec - unload + Qcs ==> unload = sprec - canopy + Qcs
!
!   Snowpack Balance (forest):
!     swe_d = unload + rain - runoff ==>
!       canopy + swe_d = sprec + rain + Qcs - runoff
!     prec = sprec + rain
!     sum_rain  = sum_sprec - sum_prec
!
!   Snowpack Balance (non-forest):
!     swe_d = sprec + rain - runoff + subl + salt + susp + subgrid +
!       glaciermelt
!
!   Everywhere:
!     w_balance = sum_prec + sum_Qcs - sum_runoff + sum_subl +
!       sum_trans - canopy_int - swe_depth + sum_glacmelt
!
!   The related variables that would need to be brought in are:
!      d_canopy_int,sum_d_canopy_int,sum_unload

! This subroutine is called for the case where SnowTran-3D is not
!   run.  The subroutine WATERBAL_SNOWTRAN is used if the model
!   simulation includes SnowTran-3D.
!     w_balance = swe_depth_old - swe_depth + prec - runoff +
!    &  glacier_melt - swesublim + canopy_int_old - canopy_int -
!    &  d_canopy_int + Qcs + canopy_unload

! Do the snowpack.
!     w_balance = swe_depth_old - swe_depth + prec - runoff -
!    &  glacier_melt - swesublim

! Do the canopy.
!     w_balance = canopy_int_old - canopy_int + d_canopy_int +
!    &  Qcs - canopy_unload

! Do the snowpack and canopy store.
      w_balance = swe_depth_old - swe_depth + prec - runoff + &
     &  glacier_melt - swesublim + canopy_int_old - canopy_int + &
     &  Qcs

      if (abs(w_balance).gt.1.0e-5) &
     &  print*,'water imbalance found, iter =',iter,' ',w_balance

      return
      end SUBROUTINE WATERBAL_SNOWPACK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SNOWPACK(swe_depth,snow_d,ro_snow_grid,&
     &  prec,ro_water,ro_nsnow,runoff,Qm,xLf,dt,rain,sprec,&
     &  sum_prec,sum_runoff,soft_snow_d,sum_sprec,ro_snow,&
     &  snow_depth,sprec_grnd,vegtype,glacier_melt,sum_glacmelt,&
     &  swemelt,canopy_unload,Qe,sfc_sublim_flag,sum_sfcsublim,&
     &  sum_swemelt,corr_factor,icorr_factor_index,swesublim,&
     &  ro_snowmax)

      implicit none

      real ro_snowmax,runoff,Qm,swe_depth,potmelt,swemelt,dt,&
     &  ro_water,xLf,snow_depth,ro_snow_grid,snow_d_melt,dz_water,&
     &  soft_snow_d,prec,rain,snow_d,sum_sprec,sum_prec,&
     &  sum_runoff,ro_nsnow,sprec,ro_snow,snow_d_new,sprec_grnd,&
     &  vegtype,glacier_melt,sum_glacmelt,canopy_unload,Qe,&
     &  xLsublim,potsublim,swesublim,snow_d_sublim,sfc_sublim_flag,&
     &  sum_sfcsublim,sum_swemelt,corr_factor,potmelt_tmp
      integer icorr_factor_index

      runoff = 0.0

! SURFACE SUBLIMATION.

! Whether static-surface (non-blowing snow) sublimation is included
!   in the model calculations is controlled by the sfc_sublim_flag.
!   I am waiting for the flux-tower data Matthew and I are collecting
!   in Alaska, to compare with the model simulations, before
!   including this part of the model in all simulations.

! If the sfc_sublim_flag is turned on, the latent heat flux (Qe)
!   calculated in ENBAL is used to add/remove snow from the snowpack.
!   xLsublim = xLf + xLv = 2.5104x10^6 J/kg + 3.334x10^5 J/kg, and
!   potsublim is in m swe.

      if (swe_depth.gt.0.0  .and.  sfc_sublim_flag.eq.1.0) then
        if (Qe.lt.0.0) then

! Compute the snow-surface sublimation (m, swe).
          xLsublim = 2.844e6
          potsublim = (- dt) * Qe / (ro_water * xLsublim)
          swesublim = min(potsublim,swe_depth)

! Save a summing array of the static surface snow sublimation.
          sum_sfcsublim = sum_sfcsublim + swesublim

! Compute the change in snow depth.  Assume that this sublimated
!   snow does not change the snow density and does not change the
!   soft snow depth.  It only reduces the snow depth and the
!   associated swe depth.
          swe_depth = swe_depth - swesublim
          if (swe_depth.eq.0.0) then
            snow_depth = 0.0
          else
            snow_d_sublim = swesublim * ro_water / ro_snow_grid
            snow_depth = snow_depth - snow_d_sublim
          endif
        else
          swesublim = 0.0
        endif
      else
        swesublim = 0.0
      endif

! MELTING.

! If melting occurs, decrease the snow depth, and place the melt
!   water in the 'runoff' variable.  Keep track of the liquid water
!   produced.

      if (Qm.gt.0.0) then

! Compute the snow melt (m).
        potmelt = dt * Qm / (ro_water * xLf)

! Account for any snowmelt data assimilation.
        if (icorr_factor_index.lt.0) then
          potmelt_tmp = potmelt * corr_factor
          swemelt = min(potmelt_tmp,swe_depth)
! Handle the case of no snowmelt data assimilation.
        else
          swemelt = min(potmelt,swe_depth)
        endif

! Compute any glacier or permanent snow-field melt (m water equiv.).
        if (vegtype.eq.20.0) then
          glacier_melt = potmelt - swemelt
        else
          glacier_melt = 0.0
        endif

! Save a summing array of the glacier melt.
        sum_glacmelt = sum_glacmelt + glacier_melt

! Save the runoff contribution.
        runoff = runoff + glacier_melt

! Save a summing array of the snow melt.
        sum_swemelt = sum_swemelt + swemelt

! Compute the change in snow depth.
        snow_d_melt = swemelt * ro_water / ro_snow_grid
        snow_depth = snow_depth - snow_d_melt
        snow_depth = max(0.0,snow_depth)

! Clip snow depth values that are within machine precision of
!   zero.  Add the clipped moisture to the swe_melt value to
!   satisfy the moisture budget.
        if (snow_depth.lt.1.0e-6) then
          swemelt = swemelt + snow_depth
          snow_depth = 0.0
        endif

! Compute the changes in snow density resulting from the melt.
!   Assume that the melted snow is redistributed through the new
!   snow depth up to a maximum density.  Any additional melt water
!   is added to the runoff.
        if (snow_depth.eq.0.0) then
          ro_snow_grid = ro_snowmax
          runoff = runoff + swemelt
        else
          ro_snow_grid = swe_depth * ro_water / snow_depth
        endif

        if (ro_snow_grid.gt.ro_snowmax) then
          dz_water = snow_depth * &
     &      (ro_snow_grid - ro_snowmax) / ro_water
          ro_snow_grid = ro_snowmax
          swe_depth = snow_depth * ro_snow_grid / ro_water
          runoff = runoff + dz_water
        else
          swe_depth = snow_depth * ro_snow_grid / ro_water
        endif

        soft_snow_d = 0.0

      else

! These prevent values from the previous time step from being
!   carried through to the next time step.
        swemelt = 0.0
        glacier_melt = 0.0

      endif

! PRECIPITATION.

! Precipitation falling as rain on snow contributes to a snow
!   density increase, precipitation falling as snow adds to the
!   snow depth, and rain falling on bare ground contributes to the
!   runoff.  This latest version of this section follows Justin
!   Pflug's code updates.

! We have precipitation.
      if (prec.gt.0.0) then
        rain = prec - sprec

! If there is snow on the ground, all of the precipitation (rain
!   and snowfall) is added to the snowpack.  If there is no
!   snowpack, then snowfall builds a new snowpack, and rain goes
!   into runoff.
        if (snow_depth.gt.0.0) then
          swe_depth = swe_depth + rain + sprec_grnd
        else
          swe_depth = sprec_grnd
          runoff = runoff + rain
        endif

! Update the new snow depth, the total snow depth, and the snow
!   density.
        snow_d_new = ro_water / ro_nsnow * sprec_grnd
        snow_depth = snow_depth + snow_d_new
        if (snow_depth.gt.0.0) then
          ro_snow_grid = ro_water * swe_depth / snow_depth
        endif

! If the density threshold is exceeded, adjust that and place the
!   excess moisture in the runoff array.
        if (ro_snow_grid.gt.ro_snowmax) then
          dz_water = snow_depth * (ro_snow_grid - ro_snowmax) / ro_water
          ro_snow_grid = ro_snowmax
          swe_depth = snow_depth * ro_snow_grid / ro_water
          runoff = runoff + dz_water
        endif

! Here we handle the case where there is no precipitation, but
!   there is snow falling from the canopy to the snowpack.
      else
        rain = 0.0
        if (sprec_grnd.gt.0.0) then
          swe_depth = swe_depth + sprec_grnd
          snow_d_new = ro_water / ro_snow * sprec_grnd
          snow_depth = snow_depth + snow_d_new
          ro_snow_grid = ro_water * swe_depth / snow_depth
        endif
      endif

! The following are set up to be compatible with SnowTran-3D, and
!   are in snow-depth units.  The sum_sprec corrections are done
!   in the SnowTran-3D code.

! Assume any rain sets the soft snow depth to zero.
      if (rain.eq.0.0) then
        soft_snow_d = soft_snow_d + sprec_grnd * ro_water / ro_snow
      else
        soft_snow_d = 0.0
      endif

      snow_d = swe_depth * ro_water / ro_snow
      sum_sprec = sum_sprec + sprec_grnd

! The following are in swe-depth units.
      sum_prec = sum_prec + prec
      sum_runoff = sum_runoff + runoff

      return
      end SUBROUTINE SNOWPACK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE DDENSITY(ro_snow_grid,swe_depth,Tf,Tsfc,dt,A1,A2,&
     &  snow_depth,ro_water,ro_snowmax)

      implicit none

      real snow_depth,Tsg,Tf,Tsnow,Tsfc,ro_snow_grid,dt,A1,A2,&
     &  swe_depth_star,ro_snowmax,ro_water,swe_depth,ro_adjust

! ro_adjust is a snow density rate adjustment factor that can be
!   used to make the snow density increase faster (ro_adjust > 1.0)
!   or slower (ro_adjust < 1.0).
      ro_adjust = 5.0

      if (snow_depth.gt.0.0) then

! Assume that the snow-ground interface temperature is -1.0 C.
        Tsg = Tf - 1.0
        Tsnow = 0.5 * (Tsg + Tsfc)
        swe_depth_star= 0.5 * swe_depth
        ro_snow_grid = ro_snow_grid + ro_adjust * dt * &
     &    (A1 * swe_depth_star * ro_snow_grid * &
     &    exp((- 0.08)*(Tf-Tsnow)) * exp((- A2)*ro_snow_grid))
        ro_snow_grid = min(ro_snowmax,ro_snow_grid)
        snow_depth = ro_water * swe_depth / ro_snow_grid

      endif

      return
      end SUBROUTINE DDENSITY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SOLVEWB(xnew,Tf,Tair,rh,xLs,Cp,sfc_pressure)

      implicit none

      real A,B,C,ea,rh,Tair,Tf,tol,old,fprime,xLs,Cp,funct,xnew,&
     &  sfc_pressure

      integer maxiter,i

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: temperatures for Buck`s equations are in deg C, and
!   vapor pressures are in mb.  Do the adjustments so that the
!   calculations are done with temperatures in K, and vapor
!   pressures in Pa.

! Over water.
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
! Over ice.
!       A = 6.1115 * 100.0
!       B = 22.452
!       C = 272.55

! Atmospheric vapor pressure from relative humidity data.
      ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)))

! Solve for the wet bulb temperature.
      tol = 1.0e-2
      maxiter = 20
      old = Tair

      do i=1,maxiter
        fprime = 1.0 + xLs/Cp * 0.622/sfc_pressure * log(10.0) * &
     &    2353. * (10.0**(11.40 - 2353./old)) / old**2
        funct = old - Tair + xLs/Cp * 0.622/sfc_pressure * &
     &    (10.0**(11.40-2353./old) - ea)
        xnew = old - funct/fprime
        if (abs(xnew - old).lt.tol) return
        old = xnew
      end do

! If the maximum iterations are exceeded, send a message and set
!   the wet bulb temperature to the air temperature.
      write (*,102)
  102 format('max iteration exceeded when solving for Twb')
      xnew = Tair

      return
      end SUBROUTINE SOLVEWB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE NSNOWDEN(ro_nsnow,Twb,Tf,dt)

      implicit none

      real Twgmax,Tf,Twb,ro_nsnow,scalefact,dt,wt

      Twgmax = Tf + 1.0
      if (Twb.ge.258.15 .and. Twb.le.Twgmax) then
        ro_nsnow = 50. + 1.7 * (Twb - 258.15)**1.5
      elseif (Twb.lt.258.15) then
        ro_nsnow = 50.0
      else
        ro_nsnow = 158.8
      endif

! For one day time steps, this equation gives a new snow density at
!   the end of the 24 hour period which is too low, by an approximate
!   factor of X.  Thus, for a daily time step, I scale the density by
!   X before returning it to the main program.

      scalefact = 1.0
      if (dt.eq.86400.0) then
        if (ro_nsnow.le.158.8) then
          wt = 1.0 + (50.0 - ro_nsnow) / 108.8
          ro_nsnow = wt * scalefact * ro_nsnow + ro_nsnow
          ro_nsnow = min(158.8,ro_nsnow)
        endif
      endif

      return
      end SUBROUTINE NSNOWDEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE POSTPROC(ro_snow_grid,xro_snow,snow_depth,undef)

      implicit none

      real snow_depth,xro_snow,undef,ro_snow_grid

      if (snow_depth.eq.0.0) then
        xro_snow = undef
      else
        xro_snow = ro_snow_grid
      endif

      return
      end SUBROUTINE POSTPROC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CONSTS_SNOWPACK(Cp,xLs,ro_ice,xLf,Tf,A1,A2,ro_water,&
     &  Cp_snow,ro_snowmax)

      implicit none

      real Cp,xLs,ro_ice,xLf,Tf,A1,A2,ro_water,Cp_snow,ro_snowmax

      Cp = 1004.
      xLs = 2.500e6
      ro_ice = 917.0
      xLf = 3.34e5
      Tf = 273.15
      A1 = 0.0013
      A2 = 0.021
      ro_water = 1000.0
      Cp_snow = 2106.
      ro_snowmax = 550.0

      return
      end SUBROUTINE CONSTS_SNOWPACK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MULTI_LAYER_SNOW(KK,ro_layer,Tf,dt,ro_water,&
     &  ro_ice,T_old,snod_layer,swed_layer,Qm,ro_snowmax,rain,&
     &  xLf,Cp_snow,melt_flag,runoff,tslsnowfall,ro_nsnow,&
     &  sprec,Tsfc,tsls_threshold,gamma,max_layers,change_layer,&
     &  dz_snow_min,snow_depth,swe_depth,undef,canopy_unload,&
     &  vegtype,glacier_melt,sum_glacmelt,sum_swemelt,snow_d,&
     &  Qe,sfc_sublim_flag,sum_sfcsublim,soft_snow_d,ro_snow,&
     &  sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,&
     &  ro_snow_grid,xro_snow,swesublim,A1,A2,windspd_2m,&
     &  sfc_pressure,diam_layer,flux_layer,corr_factor,&
     &  icorr_factor_index)

      use snowmodel_inc
      implicit none

      integer KK,max_layers,k
      integer melt_flag(nz_max)

      real snod_layer(nz_max)
      real swed_layer(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      real gamma(nz_max)
      real diam_layer(nz_max)
      real flux_layer(nz_max)
      real frac_liq(nz_max)

      real Tf,dt,ro_water,ro_ice,Qm,ro_snowmax,rain,xLf,Cp_snow,&
     &  runoff,tslsnowfall,ro_nsnow,sprec,Tsfc,tsls_threshold,&
     &  dz_snow_min,snow_depth,swe_depth,undef,change_layer,&
     &  canopy_unload,vegtype,glacier_melt,sum_glacmelt,sum_swemelt,&
     &  soft_snow_d,Qe,sfc_sublim_flag,sum_sfcsublim,snow_d,&
     &  ro_snow,sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,&
     &  ro_snow_grid,xro_snow,swesublim,A1,A2,windspd_2m,&
     &  sfc_pressure,corr_factor

      integer icorr_factor_index

! THIS IS THE MULTI-LAYER SNOWPACK MODEL.

! Note there is some confusion with the dy - dz notation used here.
!   In the multi-layer code 'y' is the vertical coordinate.  This is
!   a hold-over from a time when my temperature solution code had
!   y going up-down.

! Compute the snow density change due to compaction and the impact
!   of wind on the top snow layer.
      CALL DDENSITY_ML(ro_layer,Tf,dt,ro_water,ro_snowmax,&
     &  T_old,KK,snod_layer,A1,A2,windspd_2m)

! Calculate the rainfall from prec and sprec.
      if (prec.gt.0.0) then
        rain = prec - sprec
      else
        rain = 0.0
      endif

! Distribute surface melt and rain precipitation through the snowpack.
      CALL MELT_SNOW_ML(KK,swed_layer,ro_water,ro_layer,Qm,dt,&
     &  snod_layer,ro_snowmax,rain,xLf,Cp_snow,Tf,T_old,melt_flag,&
     &  runoff,canopy_unload,swe_depth,snow_depth,vegtype,&
     &  glacier_melt,sum_glacmelt,sum_swemelt,soft_snow_d,Qe,&
     &  sfc_sublim_flag,sum_sfcsublim,swesublim,ro_snow,&
     &  corr_factor,icorr_factor_index)

! Account for the accumulation of snow precipitation on the snowpack.
      CALL PRECIP_ML(KK,ro_layer,snod_layer,ro_water,tslsnowfall,&
     &  swed_layer,ro_nsnow,T_old,Tsfc,tsls_threshold,dt,&
     &  melt_flag,soft_snow_d,ro_snow,sum_sprec,sprec_grnd_ml,&
     &  sum_prec,prec,sum_runoff,runoff,snow_d,snow_depth,swe_depth,&
     &  diam_layer)

! Merge layers if the number of layers exceeds some maximum number of
!   layers or if a layer gets thinner than some minimum thickness.
      CALL MERGE_LAYERS_ML(KK,ro_layer,snod_layer,swed_layer,&
     &  T_old,ro_water,max_layers,change_layer,dz_snow_min,melt_flag,&
     &  diam_layer)

! The grain-growth model can account for the liquid fraction in the
!   snow.  The model deals with three different wetness conditions
!   differently:
!     frac_liq < 1.0e-4 is considered dry
!     1.0e-4 <= frac_liq < 0.09 is considered wet
!     frac_liq >= 0.09 is considered very wet
! Note that there is also a different vapor diffusion threshold
!   that is also used (0.02).
! I have not implemented these yet, and for now assume that the
!   snow is dry.
      do k=1,KK
        frac_liq(k) = 0.0
      enddo

! Update the grain size using the SNTHERM equations.
      CALL GET_GRAIN_SIZE_SNTHERM(KK,dt,ro_layer,undef,&
     &  diam_layer,frac_liq,T_old,snod_layer,Tf,sfc_pressure,&
     &  flux_layer)

! Calculate the temperature of each snow layer.
      CALL SNOWTEMP_ML(gamma,T_old,Tsfc,KK,dt,ro_layer,Cp_snow,&
     &  Tf,snod_layer,melt_flag,diam_layer)

! Postprocess the data.
      CALL POST_PROC_ML(KK,snod_layer,snow_depth,swe_depth,undef,&
     &  swed_layer,gamma,ro_layer,melt_flag,T_old,Tf,ro_snow_grid,&
     &  ro_water,xro_snow,ro_snow)

      return
      end SUBROUTINE MULTI_LAYER_SNOW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GET_GRAIN_SIZE_SNTHERM(nlayers,dt,ro_layer,undef,&
     &  diam_layer,frac_liq,T_layer,dz,Tf,press,&
     &  flux_layer)

! Someday I need to make the above terminology consistent with the
!   rest of the code (like that below).
!     SUBROUTINE GET_GRAIN_SIZE_SNTHERM(KK,dt,ro_layer,diam_layer,
!    &  frac_liq,T_old,snod_layer,Tf,sfc_pressure)

! This code follows Jordan 1991 (Tech Doc for SNTHERM.89; this is
!   referenced as "SNTH89" in the comments below).  The equations
!   that are solved are those in the original SNTHERM89 code.

! De0 : Effective diffusion coefficient for water vapor in
!       snow (m^2/s) at 100000 Pa and 273.15K
! press : Barometric pressure (Pa)
! dc(nlayers) : Vapor diffusion coefficient for = de*change in
!               saturation vapor pressure with temperature
! vapor_flux(nlayers) : Average diffusive flux at nodal centroid
!                       [(upper+lower)/2] (kg/m^2 s)
! T_layer(nlayers) : Temperature in center (node) of layer (K)
! g1 : Grain growth constant for dry snow (5.0e-7)
! g2 : Grain growth constant for wet snow (4.0e-12)
! diam : Nodal effective grain diameter (m)
! frac_liq : Volume fraction of liquid water
! dt : Time step (s)
! e0 : Saturation vapor pressure at 0 degrees C (mb)
! Rw : Gas constant for water vapor (461.5) (j/kg-K)
! ro_ice : Ice density (917.0) (kg/m^3)
! ro_snow(nlayers) : Nodal bulk snow density (kg/m^3)

      use snowmodel_inc
      implicit none

!     integer k,KK
      integer nlayers,k,nzm

      real diam_layer(nz_max)
      real flux_layer(nz_max)
      real ro_layer(nz_max)
      real T_layer(nz_max)
      real frac_liq(nz_max)
      real dz(nz_max)
      real porosity(nz_max)
      real diff_coef(nz_max)
      real C_kT(nz_max)
      real vapor_flux(nz_max)

! These are the values at the control volume boundaries.
      real T_layer_bndry(nz_max+1)
      real dz_bndry(nz_max+1)
      real diff_coef_bndry(nz_max+1)

      real g1,g2,ro_ice,De0,Des,dt,Tf,press,grain_scale,pi,&
     &  xLvi_o_Rw,c1i,xLvw_o_Rw,c1w,vaporvol,Uv_bot,Uv_top,&
     &  diam_mm,undef

!     real e0,Rw,xLvw,xLvi

      data g1 /5.0e-7/
      data g2 /4.0e-12/
      data ro_ice /917.0/
      data De0 /9.2e-5/

! Calculate the vapor diffusion constants for water(w) and ice(i).
!   These are constants, so I just save the values as data values.
!     data e0 /613.6/
!     data Rw /461.5/
!     data xLvw /2.505e6/
!     data xLvi /2.838e6/
!     xLvw_o_Rw = xLvw / Rw
!     c1w = e0 * exp(xLvw_o_Rw / Tf) / Rw
!     xLvi_o_Rw = xLvi / Rw
!     c1i = e0 * exp(xLvi_o_Rw / Tf) / Rw
!     print *,xLvw_o_Rw
!     print *,c1w
!     print *,xLvi_o_Rw
!     print *,c1i

      data xLvw_o_Rw /5427.952/
      data c1w /5.6739e8/
      data xLvi_o_Rw /6149.513/
      data c1i /7.9639e9/

      pi = 2.0 * acos(0.0)

! This 'grain_scale' controls allows the user additional control
!   on the grain-growth calculations.  Values greater than 1.0
!   grows grains faster, and values less than 1.0 grows grains
!   slower.
!     data grain_scale /1.0/
!     data grain_scale /0.5/
      data grain_scale /2.5/

! Initialize the flux output variable so the no-layer outputs make
!   sense.
      do k=1,nz_max
        flux_layer(k) = undef
      enddo

! Define the snow porosity.
      do k=1,nlayers
        porosity(k) = 1.0 - (ro_layer(k) / ro_ice)

! The only place this value should have problems is if ro_layer is
!   undef.  It should be well-constrained before getting to this
!   subroutine.
        if (porosity(k) .lt. 0.0) porosity(k) = 0.0
      enddo

! Diffusion coefficient for each snow layer.
      do k=1,nlayers

! Eqn. 20 of SNTH89.
        if (frac_liq(k) .gt. 0.02) then
          C_kT(k) = c1w * exp(- xLvw_o_Rw / T_layer(k)) * &
     &      (xLvw_o_Rw / T_layer(k) - 1.0) / T_layer(k)**2
        else
          C_kT(k) = c1i * exp(- xLvi_o_Rw / T_layer(k)) * &
     &      (xLvi_o_Rw / T_layer(k) - 1.0) / T_layer(k)**2
        endif

! SNTH89 assumed the porosity does not affect the fluxes.  SNTH89
!   had a note that said "The diffusive vapor flux in snow is
!   customarily taken as independent of porosity, which is
!   generally a consequence of the 'hand-to-hand' process of
!   vapor diffusion."  Conceptually, it seems like if the porosity
!   goes to zero the fluxes should stop.  I do that here.  If you
!   don't like this idea, you can just comment out these 3 lines.
! Scale the flux by the available vapor.
        vaporvol = porosity(k) - frac_liq(k)
        vaporvol = max(0.0,vaporvol)
        C_kT(k) = vaporvol * C_kT(k)

! Left part of Eqn. 21 of SNTH89.
        Des = De0 * (100000.0 / press) * (T_layer(k) / Tf)**6

! Eqn. 21 of SNTH89.  The vapor diffusion coefficents.
        diff_coef(k) = Des * C_kT(k)

        if (diff_coef(k).le.0.0) then
          print *,'diff_coef(k) <= 0.0',k,diff_coef(k)
          print *,porosity(k),ro_layer(k)
        endif

      enddo

! Include the boundary information in the required arrays.  This
!   information, and the flux calcuations below, follow Glen's 3D
!   thermal sea ice growth model code (see /nice/heat/temp_3d.f,
!   or SeaIce-3D).

! Number of shifted layers.
      nzm = nlayers + 1

! Control volume size.
      dz_bndry(1) = 0.0
      do k=2,nzm
        dz_bndry(k) = dz(k-1)
      enddo
      dz_bndry(nzm+1) = 0.0

!     do k=1,nzm+1
!       print *, k,dz_bndry(k)
!     enddo

! Temperature.
      T_layer_bndry(1) = T_layer(1)
      do k=2,nzm
        T_layer_bndry(k) = T_layer(k-1)
      enddo
      T_layer_bndry(nzm+1) = T_layer(nzm-1)

!     do k=1,nzm+1
!       print *, k,T_layer_bndry(k)
!     enddo

! Diffusion coefficients.
      diff_coef_bndry(1) = diff_coef(1)
      do k=2,nzm
        diff_coef_bndry(k) = diff_coef(k-1)
      enddo
      diff_coef_bndry(nzm+1) = diff_coef(nzm-1)

!     do k=1,nzm+1
!       print *, k,diff_coef_bndry(k)
!     enddo

! Calculate the vapor flux across the control volume walls (Eqn.
!   21 of SNTH89).  This is: flux = - diff_coef * dT/dz, where 
!   diff_coef = Des * C_kT.  See page 45 of Patankar (1980) for
!   a description of what is being done with the harmonic mean
!   calculations.  Here I am solving Patankar's Eqn. 4.8, with
!   delx_e- = 0.5*dx_P, and delx_e+ = 0.5*dx_E.
      do k=2,nzm

        if (dz_bndry(k-1).eq.0.0 .and. dz_bndry(k).eq.0.0) then
          Uv_bot = 0.0
        elseif (dz_bndry(k).eq.0.0 .and. dz_bndry(k+1).eq.0.0) then
          Uv_top = 0.0
        else
          Uv_bot = - 2.0 * (T_layer_bndry(k) - T_layer_bndry(k-1)) / &
     &      (dz_bndry(k-1) / diff_coef_bndry(k-1) + &
     &      dz_bndry(k) / diff_coef_bndry(k))
          Uv_top = - 2.0 * (T_layer_bndry(k+1) - T_layer_bndry(k)) / &
     &      (dz_bndry(k) / diff_coef_bndry(k) + &
     &      dz_bndry(k+1) / diff_coef_bndry(k+1))
        endif

! Assume the flux at the center of the control volume is the average
!   of the fluxes at the control volume walls.  Also note that we
!   don't care about the direction of the fluxes; we are just
!   assuming that the vapor transport is growing grains, regardless
!   of the direction.  This is used in the grain growth algorithm.
        vapor_flux(k-1) = (abs(Uv_bot) + abs(Uv_top)) / 2.0

! Save a record of the layer fluxes, including the direction of the
!   flow.  Here I am saving the flux across the top of each layer.
!   Under the assumption that the flux across the bottom of the
!   bottom layer is zero, this should be enough information to
!   calculate d_flux/d_z and get the mass loss and/or gain in each
!   layer.
        flux_layer(k-1) = Uv_top

      enddo

! Because the zero temperature gradient assumed at the boundaries
!   is not realistic, set the boundary fluxes equal to those just
!   inside the boundaries. 
      vapor_flux(1) = vapor_flux(2)
      vapor_flux(nlayers) = vapor_flux(nlayers-1)

! Because below we don't allow the abs(fluxes) to be over 1.0e-6,
!   do the same thing here.
      do k=1,nlayers

        if (flux_layer(k).lt.-1.0e-6) then
          flux_layer(k) = -1.0e-6
        elseif (flux_layer(k).gt.1.0e-6) then
          flux_layer(k) = 1.0e-6
        endif

! Convert these to fluxes per dt (instead of per sec).  Without
!   this the values are something like 10^-11.
          flux_layer(k) = dt * flux_layer(k)

      enddo

! Update the snow grain diameter.
      do k=1,nlayers

        if (diam_layer(k) .le. 0.0) then
          print *, 'execution halted because diam_layer <= 0.0'
          print *, 'layer = ',k,'  diam_layer(k) =',diam_layer(k)
          stop
        endif

! Dry snow: The cut-off bewteen dry and wet snow is arbitrarily
!   set at 0.0001.
        if (frac_liq(k) .lt. 1.0e-4) then

! The max vapor flux available for growth is arbitrarily set at
!   1.0e-6.  This can be increased if you want to allow larger
!   grains to grow, if the temperature gradients are available to
!   drive greater fluxes.  Here I have also included a scaling
!   term that can be used to increase or decrease the growth rate
!   of the grains to better match any observational datasets you
!   might have.  Eqn. 33 of SNTH89.
          if (abs(vapor_flux(k)) .lt. 1.0e-6) then
            diam_layer(k) = diam_layer(k) + grain_scale * dt * g1 * &
     &        abs(vapor_flux(k)) / diam_layer(k)
          else
            diam_layer(k) = diam_layer(k) + grain_scale * dt * g1 * &
     &        1.0e-6 / diam_layer(k)
          endif

! Wet snow: Different equations for liquid volume fraction
!   above and below 0.09.
        else
          if (frac_liq(k) .lt. 0.09) then 
! Eqn. 34a of SNTH89.
            diam_layer(k) = diam_layer(k) + grain_scale * dt * g2 * &
     &        (frac_liq(k)+0.05) / diam_layer(k)
          else
! Eqn. 34b of SNTH89.
            diam_layer(k) = diam_layer(k) + grain_scale * dt * g2 * &
     &        0.14 / diam_layer(k)
          endif
        endif

! Max grain size set at 5mm, based on Arctic Alaska depth hoar
!   observations (larger sizes are possible but not common).
        if (diam_layer(k) .gt. 5.0e-3) diam_layer(k) = 5.0e-3

      enddo

      return
      end SUBROUTINE GET_GRAIN_SIZE_SNTHERM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GETGAMMA(KK,ro_layer,gamma)

      use snowmodel_inc
      implicit none

      integer k,KK
      real ro_layer(nz_max)
      real gamma(nz_max)

! Compute the snow thermal conductivity (gamma) from the snow density.
      do k=1,KK
        if (ro_layer(k).lt.156.0) then
          gamma(k) = 0.023 + 0.234 * (ro_layer(k)/1000.0)
        else
          gamma(k) = 0.138 - 1.01 * (ro_layer(k)/1000.0) + 3.233 * &
     &      (ro_layer(k)/1000.0)**2
        endif
      enddo

      return
      end SUBROUTINE GETGAMMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
      SUBROUTINE GET_THERMAL_CONDUCTIVITY(KK,ro_layer,diam_layer,gamma)

! This program calculates the thermal conductivity for a given
!   grain size.

! The grain size here is provided in mm, and is assumed to range
!   from 0.5 mm (wind slab) to 5.0 mm (depth hoar).

      use snowmodel_inc
      implicit none

      integer k,KK
      real ro_layer(nz_max)
      real gamma(nz_max)
      real diam_layer(nz_max)

      real cond_slab,cond_hoar,wt,x1,x2,xx

! Convert the diameter limits from mm to m.
      x1 = 0.5 / 1000.0
      x2 = 5.0 / 1000.0

      do k=1,KK
        if (ro_layer(k).lt.180.0) then
          cond_slab = 3.94e-2 + 2.00e-4*ro_layer(k)
        else
          cond_slab = 1.55e-1 - 1.02e-3*ro_layer(k) + &
     &      3.21e-6*ro_layer(k)**2
        endif

        cond_hoar = 3.00e-2 + 2.00e-4*ro_layer(k)

! Calculate the weighting factors.
        xx = diam_layer(k)
        wt = (xx - x1) / (x2 - x1)

! Calculate the thermal conductivity.
        gamma(k) = wt * cond_hoar + (1.0 - wt) * cond_slab

      enddo

      return
      end SUBROUTINE GET_THERMAL_CONDUCTIVITY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE POST_PROC_ML(KK,snod_layer,snow_depth,swe_depth,undef,&
     &  swed_layer,gamma,ro_layer,melt_flag,T_old,Tf,ro_snow_grid,&
     &  ro_water,xro_snow,ro_snow)

      use snowmodel_inc
      implicit none

      integer k,KK

      real snod_layer(nz_max)
      real swed_layer(nz_max)
      real gamma(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      integer melt_flag(nz_max)

      real snow_depth,swe_depth,undef,Tf,ro_snow_grid,ro_water,&
     &  xro_snow,ro_snow

! Calculate the total snow and swe depth, and the bulk snow density.
      snow_depth = 0.0
      swe_depth = 0.0
      do k=1,KK
        snow_depth = snow_depth + snod_layer(k)
        swe_depth = swe_depth + swed_layer(k)
      enddo
      if (snow_depth.le.0.0) then
        ro_snow_grid = ro_snow
      else
        ro_snow_grid = swe_depth * ro_water / snow_depth
      endif

! Set any areas outside the snowpack to undef.
      do k=KK+1,nz_max
        gamma(k) = undef
        ro_layer(k) = undef
        T_old(k) = undef + Tf
        melt_flag(k) = nint(undef)
        snod_layer(k) = undef
        swed_layer(k) = undef
      enddo

! Clean up the snow density array so there are no values when
!   there is no snow.
      if (snow_depth.eq.0.0) then
        xro_snow = undef
      else
        xro_snow = ro_snow_grid
      endif

      return
      end SUBROUTINE POST_PROC_ML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MERGE_LAYERS_ML(KK,ro_layer,snod_layer,swed_layer,&
     &  T_old,ro_water,max_layers,change_layer,dz_snow_min,melt_flag,&
     &  diam_layer)

      use snowmodel_inc
      implicit none

      integer k,KK,kkk,k_small,max_layers,icount,kkkk

      real swed_layer(nz_max)
      real ro_layer(nz_max)
      real snod_layer(nz_max)
      real T_old(nz_max)
      real diam_layer(nz_max)
      integer melt_flag(nz_max)

      real snod_layer_small,ro_water,dz_snow_min,change_layer

! Merge layers if the number of layers exceeds some maximum number of
!   layers or if a layer gets thinner than some minimum thickness.
!   Do this in snow_depth space because that is the grid the snow
!   temperatures are being calculated on.

! If any layer is thinner than the minimum layer thickness, merge it
!   with the layer below.  If that layer is layer 1, merge it with
!   layer 2.  If there is only one layer left, let it be smaller than
!   the minimum thickness.  Don't do any of this if a new layer is
!   being built; only do it for layers below the top layer.
      change_layer = 0.0

! Count how many layers are less than the minimum thickness, excluding
!   the case where there is only one layer.
      icount = 0
      if (KK.gt.1) then
!       do k=1,KK
        do k=1,KK-1
          if (snod_layer(k).lt.dz_snow_min) then
            icount = icount + 1
          endif
        enddo
      endif

! Note that if two thin layers are together, the merge may take
!   out the other one.
      do k=1,icount
        change_layer = 1.0

! This gets and processes the last occurance.
!       do kkkk=1,KK
        do kkkk=1,KK-1
          if (snod_layer(kkkk).lt.dz_snow_min) then
            k_small = kkkk
          endif
        enddo

        if (k_small.eq.1) then
          snod_layer(1) = snod_layer(1) + snod_layer(2)
          swed_layer(1) = swed_layer(1) + swed_layer(2)
          ro_layer(1) = swed_layer(1) * ro_water / snod_layer(1)
          T_old(1) = T_old(2)
          diam_layer(1) = diam_layer(2)
          melt_flag(1) = melt_flag(2)
          KK = KK - 1
          do kkk=2,KK
            snod_layer(kkk) = snod_layer(kkk+1)
              swed_layer(kkk) = swed_layer(kkk+1)
            ro_layer(kkk) = swed_layer(kkk) * ro_water / snod_layer(kkk)
            T_old(kkk) = T_old(kkk+1)
            diam_layer(kkk) = diam_layer(kkk+1)
            melt_flag(kkk) = melt_flag(kkk+1)
          enddo
        else
          snod_layer(k_small-1) = snod_layer(k_small-1) + &
     &      snod_layer(k_small)
          swed_layer(k_small-1) = swed_layer(k_small-1) + &
     &      swed_layer(k_small)
          ro_layer(k_small-1) = swed_layer(k_small-1) * ro_water / &
     &      snod_layer(k_small-1)
          T_old(k_small-1) = T_old(k_small)
          diam_layer(k_small-1) = diam_layer(k_small)
          melt_flag(k_small-1) = melt_flag(k_small)
          KK = KK - 1
          do kkk=k_small,KK
            snod_layer(kkk) = snod_layer(kkk+1)
            swed_layer(kkk) = swed_layer(kkk+1)
            ro_layer(kkk) = swed_layer(kkk) * ro_water / snod_layer(kkk)
            T_old(kkk) = T_old(kkk+1)
            diam_layer(kkk) = diam_layer(kkk+1)
            melt_flag(kkk) = melt_flag(kkk+1)
          enddo
        endif
      enddo

! Where the number of layers exceeds some maximum number of layers,
!   find the thinnest layer and merge it with the one below.  For the
!   case where the thinnest layer is the bottom layer, merge it with
!   layer 2.
      if (KK.eq.max_layers+1) then
        change_layer = 1.0
! Find the thinnest layer.
        snod_layer_small = 1000.0
        do k=1,KK
          if (snod_layer(k).lt.snod_layer_small) then
            snod_layer_small = snod_layer(k)
            k_small = k
          endif
        enddo

! Adjust accordingly.  Note that layers below the thin layer do not
!   change, unless the thin layer is layer 1.  Also, since the layer
!   is thin, assign the new layer the thick layer temperature.
        if (k_small.eq.1) then
          snod_layer(1) = snod_layer(1) + snod_layer(2)
          swed_layer(1) = swed_layer(1) + swed_layer(2)
          ro_layer(1) = swed_layer(1) * ro_water / snod_layer(1)
          T_old(1) = T_old(2)
          diam_layer(1) = diam_layer(2)
          melt_flag(1) = melt_flag(2)
          KK = KK - 1
          do kkk=2,KK
            snod_layer(kkk) = snod_layer(kkk+1)
              swed_layer(kkk) = swed_layer(kkk+1)
            ro_layer(kkk) = swed_layer(kkk) * ro_water / snod_layer(kkk)
            T_old(kkk) = T_old(kkk+1)
            diam_layer(kkk) = diam_layer(kkk+1)
            melt_flag(kkk) = melt_flag(kkk+1)
          enddo
        else
          snod_layer(k_small-1) = snod_layer(k_small-1) + &
     &      snod_layer(k_small)
          swed_layer(k_small-1) = swed_layer(k_small-1) + &
     &      swed_layer(k_small)
          ro_layer(k_small-1) = swed_layer(k_small-1) * ro_water / &
     &      snod_layer(k_small-1)
          T_old(k_small-1) = T_old(k_small)
          diam_layer(k_small-1) = diam_layer(k_small)
          melt_flag(k_small-1) = melt_flag(k_small)
          KK = KK - 1
          do kkk=k_small,KK
            snod_layer(kkk) = snod_layer(kkk+1)
            swed_layer(kkk) = swed_layer(kkk+1)
            ro_layer(kkk) = swed_layer(kkk) * ro_water / snod_layer(kkk)
            T_old(kkk) = T_old(kkk+1)
            diam_layer(kkk) = diam_layer(kkk+1)
            melt_flag(kkk) = melt_flag(kkk+1)
          enddo
        endif
      endif

! Now that we are done with change_layer, set it equal to k_small,
!   the position of the change.
      if (change_layer.eq.1.0) change_layer = real(k_small)

      return
      end SUBROUTINE MERGE_LAYERS_ML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MELT_SNOW_ML(KK,swed_layer,ro_water,ro_layer,Qm,dt,&
     &  snod_layer,ro_snowmax,rain,xLf,Cp_snow,Tf,T_old,melt_flag,&
     &  runoff,canopy_unload,swe_depth,snow_depth,vegtype,&
     &  glacier_melt,sum_glacmelt,sum_swemelt,soft_snow_d,Qe,&
     &  sfc_sublim_flag,sum_sfcsublim,swesublim,ro_snow,&
     &  corr_factor,icorr_factor_index)

      use snowmodel_inc
      implicit none

      integer k,KK

      real swed_layer(nz_max)
      real ro_layer(nz_max)
      real snod_layer(nz_max)
      real T_old(nz_max)
      real swemelt,extra,ro_water,swe_space,add,runoff,ro_snowmax,&
     &  rain,delta_T,xLf,Cp_snow,extra_delta_T,Tf,dt,Qm,canopy_unload,&
     &  potmelt,swe_depth,snow_depth,vegtype,glacier_melt,&
     &  sum_glacmelt,sum_swemelt,soft_snow_d,Qe,sfc_sublim_flag,&
     &  xLsublim,potsublim,swesublim,sum_sfcsublim,ro_snow,&
     &  potmelt_tmp,corr_factor

      integer icorr_factor_index

      integer melt_flag(nz_max)

! Initialize the runoff array.
      runoff = 0.0

! SURFACE SUBLIMATION.

! Whether static-surface (non-blowing snow) sublimation is included
!   in the model calculations is controlled by the sfc_sublim_flag.
!   I am waiting for the flux-tower data Matthew and I are collecting
!   in Alaska, to compare with the model simulations, before
!   including this part of the model in all simulations.

! If the sfc_sublim_flag is turned on, the latent heat flux (Qe)
!   calculated in ENBAL is used to add/remove snow from the snowpack.
!   xLsublim = xLf + xLv = 2.5104x10^6 J/kg + 3.334x10^5 J/kg, and
!   potsublim is in m swe.

      if (swe_depth.gt.0.0  .and.  sfc_sublim_flag.eq.1.0) then
        if (Qe.lt.0.0) then
! Compute the snow-surface sublimation (m, swe).
          xLsublim = 2.844e6
          potsublim = (- dt) * Qe / (ro_water * xLsublim)
          swesublim = min(potsublim,swe_depth)
! Save a summing array of the static surface snow sublimation.
          sum_sfcsublim = sum_sfcsublim + swesublim
        else
          swesublim = 0.0
        endif
      else
        swesublim = 0.0
      endif

! Modify the swe layer thicknesses, and reduce the number of layers
!   if needed.
      if (swesublim.gt.0.0) then
! Check to see whether this sublimation requires a layer reduction.
        CALL REDUCE_LAYERS(swesublim,swed_layer,KK)

! Build the new snow layer thicknesses, and recalculate the total
!   snow and swe depths.  Assume this sublimated snow does not
!   change the snow density and does not change the soft snow depth.
!   It only reduces the snow depth and the associated swe depth.
        snow_depth = 0.0
        swe_depth = 0.0
        do k=1,KK
          snod_layer(k) = swed_layer(k) * ro_water / ro_layer(k)
          snow_depth = snow_depth + snod_layer(k)
          swe_depth = swe_depth + swed_layer(k)
        enddo
      endif

! MELTING.

      if (Qm.gt.0.0) then

! Convert the melt energy to water equivalent melt depth (m).
        potmelt = dt * Qm / (ro_water * xLf)

! Account for any snowmelt data assimilation.
        if (icorr_factor_index.lt.0) then
          potmelt_tmp = potmelt * corr_factor
          swemelt = min(potmelt_tmp,swe_depth)
! Handle the case of no snowmelt data assimilation.
        else
          swemelt = min(potmelt,swe_depth)
        endif

! Compute any glacier or permanent snow-field melt (m water equiv.).
        if (vegtype.eq.20.0) then
          glacier_melt = potmelt - swemelt
        else
          glacier_melt = 0.0
        endif

! Save a summing array of the glacier melt.
        sum_glacmelt = sum_glacmelt + glacier_melt

! Save the runoff contribution.
        runoff = runoff + glacier_melt

! Save a summing array of the snow melt.
        sum_swemelt = sum_swemelt + swemelt

! In the presence of melt, zero out the soft snow layer.
        soft_snow_d = 0.0

      else

! These prevent values from the previous time step from being
!   carried through to the next time step.
        swemelt = 0.0
        glacier_melt = 0.0

      endif

! Handle the case where rain and canopy_unload fall on snow-free
!   ground (this is not included in the code below, nor in the
!   PRECIP_ML subroutine, so I include it here).
      if (swe_depth.eq.0.0) then
        runoff = runoff + rain + canopy_unload
      endif

! Deal with melting snow.

      if (swemelt.gt.0.0) then
! Check to see whether this melt leads to a reduction in layers.
        CALL REDUCE_LAYERS(swemelt,swed_layer,KK)

! Build the new snow layer thicknesses, and initiate the melt_flag.
        do k=1,KK
          snod_layer(k) = swed_layer(k) * ro_water / ro_layer(k)
          melt_flag(k) = 0
        enddo
      endif

! Add the melt, rain, and canopy unloading (assumed to be wet as rain)
!   to the remaining snow layer thicknesses, up to the maximum snow
!   density, and let the rest of the melt drain out the snowpack bottom
!   as runoff.
      extra = swemelt + rain + canopy_unload
      if (extra.gt.0.0) then
        do k=KK,1,-1
          if (extra.gt.0.0) then
            swe_space = snod_layer(k) * (ro_snowmax - ro_layer(k)) / &
     &        ro_water
            add = min(swe_space,extra)
            swed_layer(k) = swed_layer(k) + add
            extra = extra - add
          endif
          if (snod_layer(k).le.0.0) then
            ro_layer(k) = ro_snow
          else
            ro_layer(k) = swed_layer(k) * ro_water / snod_layer(k)
          endif
        enddo
        runoff = extra
      endif

! Also take into account the refreezing of this liquid in a cold
!   snowpack.  Assume that the liquid will fully warm each layer before
!   moving on to the next layer.
      extra = swemelt + rain + canopy_unload
      if (extra.gt.0.0) then
        do k=KK,1,-1
          if (extra.gt.0.0) then

! Compute the change in temperature that would result if this liquid
!   was used to freeze and raise the snow temperature.
            if (snod_layer(k).le.0.0) then
              delta_T = 0.0
            else
              delta_T = (extra * xLf) / (Cp_snow * snod_layer(k))
            endif

! Use this potential temperature change to adjust the snow
!   temperature in the presence of the liquid.
            T_old(k) = T_old(k) + delta_T
            extra_delta_T = max(0.0,T_old(k)-Tf)
            T_old(k) = min(T_old(k),Tf)

! Keep track of which layers have been pushed to Tf.  This will be
!   used in the temperature solution subroutine to fix the snow
!   temperature at Tf (if melt_flag = 1).
            if (T_old(k).eq.Tf) then
              melt_flag(k) = 1
            else
              melt_flag(k) = 0
            endif

! Define the amount of liquid this remaining temperature change
!   represents, so it can be used in the layer below (that may be
!   a different size).
            extra = Cp_snow * snod_layer(k) * extra_delta_T / xLf

          endif
        enddo
      endif

      return
      end SUBROUTINE MELT_SNOW_ML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE REDUCE_LAYERS(swemelt,swed_layer,KK)

      use snowmodel_inc
      implicit none

      integer k,KK
      real swed_layer(nz_max)
      real eps,swemelt_tmp,swemelt,excess

      eps = 1e-6
      swemelt_tmp = swemelt

! The use of eps here does not allow the vertical grid increment to
!   be less that eps.

      do k=KK,1,-1
        excess = swed_layer(k) - swemelt_tmp

!       if (excess.gt.0.0) then
        if (excess.gt.eps) then
          swed_layer(k) = excess
          KK = k
          return
!       elseif (excess.eq.0.0) then
        elseif (excess.ge.0.0 .and. excess.le.eps) then
          KK = k - 1
          return
        else
          swemelt_tmp = - excess
        endif
      enddo

! If there is no snow left.
      KK = 0

      return
      end SUBROUTINE REDUCE_LAYERS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE PRECIP_ML(KK,ro_layer,snod_layer,ro_water,tslsnowfall,&
     &  swed_layer,ro_nsnow,T_old,Tsfc,tsls_threshold,dt,&
     &  melt_flag,soft_snow_d,ro_snow,sum_sprec,sprec_grnd_ml,&
     &  sum_prec,prec,sum_runoff,runoff,snow_d,snow_depth,swe_depth,&
     &  diam_layer)

      use snowmodel_inc
      implicit none

      integer k,KK

      real snod_layer(nz_max)
      real ro_layer(nz_max)
      real swed_layer(nz_max)
      real T_old(nz_max)
      real diam_layer(nz_max)
      real ro_nsnow,ro_water,z_nsnow,tsls_threshold,&
     &  z_snowtopl,sweq_topl,Tsfc,tslsnowfall,dt,soft_snow_d,ro_snow,&
     &  sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,runoff,snow_d,&
     &  snow_depth,swe_depth

      integer melt_flag(nz_max)

! If the melt from the previous subroutine reduced the snowpack
!   to no snow, reset the time since last snowfall to the threshold,
!   otherwise you will not build a new layer on the bare ground.
      if (KK.eq.0) tslsnowfall = tsls_threshold

! Create and/or modify the snow c.v.'s to account for new snowfall.
      if (sprec_grnd_ml.gt.0.0) then
        if (tslsnowfall.ge.tsls_threshold) then
! Create a new layer if snowfall has stopped for a period of time
!   greater or equal to the defined threshold.
          KK = KK + 1
          z_nsnow = ro_water / ro_nsnow * sprec_grnd_ml
          snod_layer(KK) = z_nsnow
          ro_layer(KK) = ro_nsnow
          swed_layer(KK) =  ro_layer(KK) * snod_layer(KK) / ro_water
! Define this new snow layer to have the surface temperature.
          T_old(KK) = Tsfc
! Define this new layer to have the initial grain size (0.5 mm).
          diam_layer(KK) = 0.5 / 1000.0

          melt_flag(KK) = 0
        else
! Add to the existing top layer.
          z_nsnow = ro_water / ro_nsnow * sprec_grnd_ml
          z_snowtopl = snod_layer(KK) + z_nsnow
          sweq_topl = sprec_grnd_ml + snod_layer(KK) * ro_layer(KK) / &
     &      ro_water
          snod_layer(KK) = snod_layer(KK) + z_nsnow
          ro_layer(KK) = ro_water * sweq_topl / z_snowtopl
          swed_layer(KK) = ro_layer(KK) * snod_layer(KK) / ro_water
        endif

! Update the total swe and snow depths.
        snow_depth = 0.0
        swe_depth = 0.0
        do k=1,KK
          snow_depth = snow_depth + snod_layer(k)
          swe_depth = swe_depth + swed_layer(k)
        enddo
      endif

! Define the time since last snowfall, in hours.  Handle the case
!   where there is no snow on the ground.
      if (sprec_grnd_ml.gt.0.0) then
        tslsnowfall = 0.0
      else
        tslsnowfall = tslsnowfall + dt / 3600.0
      endif
      if (KK.eq.0) tslsnowfall = tsls_threshold

! The following are set up to be compatible with SnowTran-3D, and
!   are in snow-depth units.  The sum_sprec corrections are done
!   in the SnowTran-3D code.
      soft_snow_d = soft_snow_d + sprec_grnd_ml * ro_water / ro_snow
      snow_d = swe_depth * ro_water / ro_snow
!     sum_sprec = sum_sprec + sprec_grnd_ml * ro_water / ro_snow
      sum_sprec = sum_sprec + sprec_grnd_ml

! The following are in swe-depth units.
      sum_prec = sum_prec + prec
      sum_runoff = sum_runoff + runoff

      return
      end SUBROUTINE PRECIP_ML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE DDENSITY_ML(ro_layer,Tf,dt,ro_water,ro_snowmax,&
     &  T_old,KK,snod_layer,A1,A2,windspd_2m)

      use snowmodel_inc
      implicit none

      integer k,KK,kkk

      real snod_layer(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      real sweqstar(nz_max)
      real sweql(nz_max)

      real A1,A2,ro_water,ro_snowmax,dt,Tf,ro_adjust,C,U,alfa,&
     &  windspd_2m

! Define the density rate coefficients.
!     C = 1.00
      C = 0.50
!     C = 0.10

! Define alfa. 
      alfa = 0.2

! ro_adjust is a snow density rate adjustment factor that can be
!   used to make the snow density increase faster (ro_adjust > 1.0)
!   or slower (ro_adjust < 1.0).
      ro_adjust = 5.0

      if (KK.gt.0) then

        do k=1,KK
          sweql(k) = ro_layer(k) / ro_water * snod_layer(k)
        enddo

        do kkk=1,KK
          sweqstar(kkk) = sweql(kkk) / 2.0
          do k=kkk+1,KK
            sweqstar(kkk) = sweqstar(kkk) + sweql(k)
          enddo
        enddo

! Pre-wind-adjustment code.
!       do k=1,KK
!         ro_layer(k) = ro_layer(k) + ro_adjust * dt *
!    &      (A1 * sweqstar(k) * ro_layer(k) *
!    &      exp(-0.08*(Tf-T_old(k))) * exp(-A2*ro_layer(k)))
!         ro_layer(k) = min(ro_snowmax,ro_layer(k))
!         snod_layer(k) = sweql(k) * ro_water / ro_layer(k)
!       enddo

! Update the density of all the layers except the top layer.
        do k=1,KK-1
          ro_layer(k) = ro_layer(k) + ro_adjust * dt * &
     &      (A1 * sweqstar(k) * ro_layer(k) * &
     &      exp(-0.08*(Tf-T_old(k))) * exp(-A2*ro_layer(k)))
          ro_layer(k) = min(ro_snowmax,ro_layer(k))
          snod_layer(k) = sweql(k) * ro_water / ro_layer(k)
        enddo

! Update the density of the top layer while including the influence
!   of wind in this layer only.  This is Equation (18) for U in
!   Liston et al. (2007).  C here has just been an assigned an
!   arbitrary value; it controls the impact of U on the density
!   evolution.
        k=KK

        if (windspd_2m.ge.5.0) then
          U = C * &
     &      (5.0 + 15.0 * (1.0 - exp(-(alfa*(windspd_2m - 5.0)))))
        else
          U = 1.0
        endif
! If you want to turn this off, you can uncomment the following
!   line.
!       U = 1.0

! The equation below is the same as that for the below-surface
!   layers, except for the addition of the U term.
        ro_layer(k) = ro_layer(k) + ro_adjust * dt * &
     &    (U * A1 * sweqstar(k) * ro_layer(k) * &
     &    exp(-0.08*(Tf-T_old(k))) * exp(-A2*ro_layer(k)))
        ro_layer(k) = min(ro_snowmax,ro_layer(k))
        snod_layer(k) = sweql(k) * ro_water / ro_layer(k)

      endif

      return
      end SUBROUTINE DDENSITY_ML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SNOWTEMP_ML(gamma,T_old,Tsfc,KK,dt,ro_layer,Cp_snow,&
     &  Tf,snod_layer,melt_flag,diam_layer)

      use snowmodel_inc
      implicit none

      real gamma(nz_max)
      real diam_layer(nz_max)
      real ro_layer(nz_max)
      real snod_layer(nz_max)
      real g_b_ns(nz_max+1)
      real f_n(nz_max+1)
      real aN(nz_max)
      real aP0(nz_max)
      real aS(nz_max)
      real dely_p(nz_max+1)
      real dy_p(nz_max)
      real y_crds(nz_max+2)
      real y_wall(nz_max+1)
      real A_sub(nz_max)
      real A_super(nz_max)
      real A_main(nz_max)
      real b_vector(nz_max)
      real T_old(nz_max)
      real Sc(nz_max)
      real Sp(nz_max)

      integer melt_flag(nz_max)

      integer k,KK
      real Tsfc,T_N,bc_N,bc_S,Cp_snow,Tf,dt,Tsg

! Define the snow thermal conductivity (gamma) for each layer.
!     CALL GETGAMMA(KK,ro_layer,gamma)

! This is the thermal conductivity that depends on grain size.
      CALL GET_THERMAL_CONDUCTIVITY(KK,ro_layer,diam_layer,gamma)

      if (KK.gt.1) then

! Update the control volume information.
        CALL GETCV(KK,dy_p,snod_layer)
        CALL CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,KK)

! Compute the general equation coefficients.
        CALL GAMMA1(g_b_ns,gamma,f_n,KK)
        CALL GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,KK,&
     &    ro_layer,Cp_snow)

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Account for the boundary conditions.
!   South boundary condition:
!     For T_S = known, define 
!       bc_S = aS(1) * T_S;         where T_S = known
!     For dT_S/dn = 0, define
!       bc_S = 0.0
!       aS(1) = 0.0
!   North boundary condition:
!     For T_N = known, define 
!       bc_N = aN(KK) * T_N;        where T_N = known
!     For dT_N/dn = 0, define
!       bc_N = 0.0
!       aN(KK) = 0.0
!---------------------------------------------------------------------
!---------------------------------------------------------------------

! Define the upper and lower boundary conditions.
        T_N = Tsfc
        bc_N = aN(KK) * T_N
        bc_S = 0.0
        aS(1) = 0.0

! Provide the source terms.

! Force the source terms to produce Tf at the positions where melting
!   occurred during this time step.
        do k=1,KK
          if (melt_flag(k).eq.1) then
            Sc(k) = 10e30 * Tf
            Sp(k) = -10e30
          else
            Sc(k) = 0.0
            Sp(k) = 0.0
          endif
        enddo

! Configure the information for the matrix solver.
        CALL PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old,&
     &    dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,KK)

! Solve the system of equations.
        CALL TRISOLVE(T_old,A_sub,A_main,A_super,b_vector,KK)

      elseif (KK.eq.1) then
! Assume that the snow-ground interface temperature is -1.0 C.
        Tsg = Tf - 1.0
        T_old(1) = 0.5 * (Tsg + Tsfc)
      endif

      return
      end SUBROUTINE SNOWTEMP_ML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GETCV(KK,dy_p,snod_layer)

      use snowmodel_inc
      implicit none

      real dy_p(nz_max)
      real snod_layer(nz_max)

      integer k,KK

! Provide values of Control Volume size in the y direction.
      do k=1,KK
        dy_p(k) = snod_layer(k)
      enddo

      return
      end SUBROUTINE GETCV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old,&
     &  dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,KK)

      use snowmodel_inc
      implicit none

      real aP(nz_max)
      real aN(nz_max)
      real aS(nz_max)
      real Sp(nz_max)
      real Sc(nz_max)
      real aP0(nz_max)
      real dy_p(nz_max)
      real T_old(nz_max)
      real b_vector(nz_max)
      real A_sub(nz_max)
      real A_super(nz_max)
      real A_main(nz_max)

      integer k,KK
      real bc_S,bc_N

! Compute matrix diagonal and b coeffs.
      do k=1,KK
        aP(k) = aN(k) + aS(k) + aP0(k) - Sp(k) * dy_p(k)
        b_vector(k) = Sc(k) * dy_p(k) + aP0(k) * T_old(k)
      enddo

! Modify b to account for dirichlet boundary conditions.
      b_vector(1) = b_vector(1) + bc_S
      b_vector(KK) = b_vector(KK) + bc_N

! Prepare to call the tridiagonal solver.
      do k=1,KK-1
        A_sub(k) = - aS(k+1)
        A_super(k) = - aN(k)
      enddo

      do k=1,KK
        A_main(k) = aP(k)
      enddo

      return
      end SUBROUTINE PREPSOLVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,KK)

      use snowmodel_inc
      implicit none

      real dy_pbc(nz_max+2)
      real dely_p(nz_max+1)
      real f_n(nz_max+1)
      real dy_p(nz_max)
      real y_crds(nz_max+2)
      real y_wall(nz_max+1)

      integer k,KK
      real temp

! PRESSURE CONTROL VOLUME SIZE AND POSITION INFORMATION

! Include exterior boundary pressure grid points.
      dy_pbc(1) = 0.0
      do k=2,KK+1
        dy_pbc(k) = dy_p(k-1)
      enddo
      dy_pbc(KK+2) = 0.0

! Compute the distance between pressure grid points.
      do k=1,KK+1
        dely_p(k) = .5 * (dy_pbc(k) + dy_pbc(k+1))
      enddo

! Compute the distance between the pressure grid points and the control
!   volume wall.  (The following is true because the grid points do
!   pressure are defined to be in the center of the control volume.)
!   And then compute f_e and f_n.  These two steps are combined below.
      do k=1,KK+1
        f_n(k) = .5 * dy_pbc(k+1) / dely_p(k)
      enddo

! Compute the x and y coordinates of the pressure c.v. grid points,
!   including boundaries.
      temp = 0.0
      do k=1,KK+2
        y_crds(k) = temp + .5 * dy_pbc(k)
        temp = temp + dy_pbc(k)
      enddo

! Compute the x and y coordinates of the pressure c.v. walls.
      y_wall(1) = 0.0
      do k=2,KK+1
        y_wall(k) = y_wall(k-1) + dy_p(k-1)
      enddo

      return
      end SUBROUTINE CV_INFO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GAMMA1(g_b_ns,gamma,f_n,KK)

      use snowmodel_inc
      implicit none

      real g_b_ns(nz_max+1)
      real gamma(nz_max)
      real g_ns(nz_max+2)
      real f_n(nz_max+1)

      integer k,KK

! This provides gamma information on c.v. walls.

! Include gamma just outside of n, s boundaries.
      g_ns(1) = gamma(1)
      do k=2,KK+1
        g_ns(k) = gamma(k-1)
      enddo
      g_ns(KK+2) = gamma(KK)

! Compute gamma (diffusion coefficient) at the n, s control
!   volume boundaries using equation 4.9, p. 45.
      do k=1,KK+1
        g_b_ns(k) = 1.0/((1.0 - f_n(k))/g_ns(k) + f_n(k)/g_ns(k+1))
      enddo

      return
      end SUBROUTINE GAMMA1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,KK,&
     &  ro_layer,Cp_snow)

      use snowmodel_inc
      implicit none

      real aN(nz_max)
      real aS(nz_max)
      real aP0(nz_max)
      real dely_p(nz_max+1)
      real g_b_ns(nz_max+1)
      real dy_p(nz_max)
      real ro_layer(nz_max)

      integer k,KK
      real Cp_snow,dt

! CALCULATE THE COEFFICIENTS aP, for the general phi equation.
      do k=2,KK+1
        aN(k-1) = g_b_ns(k)   / dely_p(k)
        aS(k-1) = g_b_ns(k-1) / dely_p(k-1)
      enddo

      do k=1,KK
        aP0(k) = ro_layer(k) * Cp_snow * dy_p(k) / dt
      enddo

      return
      end SUBROUTINE GE_COEF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE TRISOLVE(x,asub,amain,asuper,b,KK)

      use snowmodel_inc
      implicit none

      real asub(nz_max)
      real asuper(nz_max)
      real amain(nz_max)
      real b(nz_max)
      real x(nz_max)
      real z(nz_max)
      real lmain(nz_max)
      real lsub(nz_max)
      real usuper(nz_max)

      integer k,KK

      lmain(1) = amain(1)
      usuper(1) = asuper(1)/lmain(1)

      do k=2,KK-1
        lsub(k-1) = asub(k-1)
        lmain(k) = amain(k) - lsub(k-1) * usuper(k-1)
        usuper(k) = asuper(k) / lmain(k)
      enddo

      lsub(KK-1) = asub(KK-1)
      lmain(KK) = amain(KK) - lsub(KK-1) * usuper(KK-1)
      z(1) = b(1) / lmain(1)

      do k=2,KK
        z(k) = 1.0 / lmain(k) * (b(k) - lsub(k-1) * z(k-1))
      enddo

      x(KK) = z(KK)

      do k=KK-1,1,-1
        x(k) = z(k) - usuper(k) * x(k+1)
      enddo

      return
      end SUBROUTINE TRISOLVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ZERO_SNOW(nx,ny,snow_depth,ro_snow_grid,ro_snow,&
     &  swe_depth,swe_depth_old,canopy_int_old,KK,sum_swemelt,&
     &  tslsnowfall,snod_layer,swed_layer,ro_layer,T_old,&
     &  sum_sprec,multilayer_snowpack,tsls_threshold,&
     &  sum_trans)

      use snowmodel_inc
      implicit none

      integer nx,ny,i,j,k

      integer multilayer_snowpack
      integer KK(nx,ny)

      real tsls_threshold,ro_snow
      real tslsnowfall(nx,ny)
      real snod_layer(nx,ny,nz_max)
      real swed_layer(nx,ny,nz_max)
      real ro_layer(nx,ny,nz_max)
      real T_old(nx,ny,nz_max)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)
      real swe_depth(nx,ny)
      real snow_depth(nx,ny)
      real ro_snow_grid(nx,ny)
      real sum_sprec(nx,ny)
      real sum_swemelt(nx,ny)
      real sum_trans(nx,ny)

      print *,'ZEROING OUT THE SNOW ARRAYS'
      print *,'ZEROING OUT THE SNOW ARRAYS'
      print *,'ZEROING OUT THE SNOW ARRAYS'

      do j=1,ny
        do i=1,nx
          canopy_int_old(i,j) = 0.0
          swe_depth_old(i,j) = 0.0
          snow_depth(i,j) = 0.0
          ro_snow_grid(i,j) = ro_snow
          swe_depth(i,j) = 0.0
          sum_sprec(i,j) = 0.0
          sum_swemelt(i,j) = 0.0
          sum_trans(i,j) = 0.0
        enddo
      enddo

      if (multilayer_snowpack.eq.1) then
        do j=1,ny
          do i=1,nx
            tslsnowfall(i,j) = tsls_threshold
            do k=1,KK(i,j)
              snod_layer(i,j,k) = 0.0
              swed_layer(i,j,k) = 0.0
              ro_layer(i,j,k) = ro_snow
              T_old(i,j,k) = 273.15
            enddo
            KK(i,j) = 0
          enddo
        enddo
      endif

      return
      end SUBROUTINE ZERO_SNOW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ZERO_SNOW_SEAICE_4(nx,ny,snow_depth,ro_snow_grid,&
     &  ro_snow,swe_depth,swe_depth_old,canopy_int_old,KK,&
     &  sum_swemelt,tslsnowfall,snod_layer,swed_layer,&
     &  ro_layer,T_old,sum_sprec,multilayer_snowpack,&
     &  tsls_threshold,iyear,diam_layer,output_path_wo_assim)

! NOTE: This program is VERY specific to the original Lagrangian
!   simulations.  If you are doing something else, you will have to
!   make numerous changes to this code and the subroutines that are
!   called from it.

      use snowmodel_inc
      implicit none

      integer nx,ny,i,j,k,iyear,icount

      parameter (icount=70000)

      integer multilayer_snowpack
      integer KK(nx,ny)

      real tsls_threshold,ro_snow
      real tslsnowfall(nx,ny)
      real snod_layer(nx,ny,nz_max)
      real swed_layer(nx,ny,nz_max)
      real ro_layer(nx,ny,nz_max)
      real T_old(nx,ny,nz_max)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)
      real swe_depth(nx,ny)
      real snow_depth(nx,ny)
      real ro_snow_grid(nx,ny)
      real sum_sprec(nx,ny)
      real sum_swemelt(nx,ny)
      real diam_layer(nx,ny,nz_max)

      real snod_init(icount)
      real swed_init(icount)
      real sden_init(icount)

      character*80 output_path_wo_assim

      print *,'TRANSFERRING THE SNOW DATA TO THE NEXT YEAR'
      print *,'TRANSFERRING THE SNOW DATA TO THE NEXT YEAR'
      print *,'TRANSFERRING THE SNOW DATA TO THE NEXT YEAR'

! Calculate the snow data initial conditions that will be used
!   in the next year.
      CALL SUPERIMPOSED_ICE(iyear,output_path_wo_assim,&
     &  snod_init,swed_init,sden_init,ro_snow)

      do j=1,ny
        do i=1,nx

          canopy_int_old(i,j) = 0.0
          swe_depth_old(i,j) = 0.0
          snow_depth(i,j) = 0.0
          ro_snow_grid(i,j) = ro_snow
          swe_depth(i,j) = 0.0
          sum_sprec(i,j) = 0.0
          sum_swemelt(i,j) = 0.0

          if (snod_init(i).ge.0.0) then
            swe_depth_old(i,j) = swed_init(i)
            snow_depth(i,j) = snod_init(i)
            ro_snow_grid(i,j) = sden_init(i)
            swe_depth(i,j) = swed_init(i)
          endif

        enddo
      enddo

      if (multilayer_snowpack.eq.1) then
        do j=1,ny
          do i=1,nx
            tslsnowfall(i,j) = tsls_threshold
            do k=1,KK(i,j)
              snod_layer(i,j,k) = 0.0
              swed_layer(i,j,k) = 0.0
              ro_layer(i,j,k) = ro_snow
              T_old(i,j,k) = 273.15
              diam_layer(i,j,k) = 0.5 / 1000.0
            enddo

            if (snod_init(i).ge.0.0) then
              KK(i,j) = 1
              snod_layer(i,j,1) = snod_init(i)
              swed_layer(i,j,1) = swed_init(i)
              ro_layer(i,j,1) = sden_init(i)
              T_old(i,j,1) = 273.15
              diam_layer(i,j,1) = 0.5 / 1000.0
            else
              KK(i,j) = 0
            endif

          enddo
        enddo
      endif

      return
      end SUBROUTINE ZERO_SNOW_SEAICE_4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ZERO_SEAICE_SNOW(nx,ny,snow_depth,ro_snow_grid,&
     &  ro_snow,swe_depth,swe_depth_old,canopy_int_old,KK,&
     &  tslsnowfall,snod_layer,swed_layer,ro_layer,T_old,&
     &  multilayer_snowpack,tsls_threshold,seaice_conc,&
     &  sum_sprec,sum_trans)

      use snowmodel_inc
      implicit none

      integer nx,ny,i,j,k

      integer multilayer_snowpack
      integer KK(nx,ny)

      real tsls_threshold,ro_snow
      real tslsnowfall(nx,ny)
      real snod_layer(nx,ny,nz_max)
      real swed_layer(nx,ny,nz_max)
      real ro_layer(nx,ny,nz_max)
      real T_old(nx,ny,nz_max)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)
      real swe_depth(nx,ny)
      real snow_depth(nx,ny)
      real ro_snow_grid(nx,ny)
      real seaice_conc(nx,ny)
      real sum_sprec(nx,ny)
      real sum_trans(nx,ny)

      do j=1,ny
        do i=1,nx
          if (seaice_conc(i,j).eq.0.0) then
            canopy_int_old(i,j) = 0.0
            swe_depth_old(i,j) = 0.0
            snow_depth(i,j) = 0.0
            ro_snow_grid(i,j) = ro_snow
            swe_depth(i,j) = 0.0
            sum_sprec(i,j) = 0.0
            sum_trans(i,j) = 0.0
          endif
        enddo
      enddo

      if (multilayer_snowpack.eq.1) then
        do j=1,ny
          do i=1,nx
            if (seaice_conc(i,j).eq.0.0) then
              tslsnowfall(i,j) = tsls_threshold
              do k=1,KK(i,j)
                snod_layer(i,j,k) = 0.0
                swed_layer(i,j,k) = 0.0
                ro_layer(i,j,k) = ro_snow
                T_old(i,j,k) = 273.15
              enddo
              KK(i,j) = 0
            endif
          enddo
        enddo
      endif

      return
      end SUBROUTINE ZERO_SEAICE_SNOW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE NSNOW_DENSITY_FROM_BLOWING_SNOW(windspd_2m,sprec,dt,&
     &  ro_nsnow_wind)

      implicit none

      real alfa,windspd_2m,ro_nsnow_wind,sprec,dt_24hours,dt

! Required constants.
      alfa = 0.2
      dt_24hours = 86400.0

! This new snow density increase due to wind, is only for the single-
!   layer snowpack.
      if (sprec.gt.0.0) then

! To define the offset I assumed under cold conditions (ro_nsnow =
!   50.0 kg/m3):
!   1) 24-hour ave wind speed >= 20 m/s gives ro_nsnow = 350 kg/m3.
!   2) 24-hour ave wind speed = 5 m/s gives ro_nsnow = 150 kg/m3,
!      (actually what really matters here is the density difference
!      of 200 kg/m3 for the wind speed difference from 5 to 20 m/s).
!   3) 24-hour ave wind speed < 5 m/s gives ro_nsnow = ro_nsnow.
!   4) It is appropriate to use an exponential decay function
!      between the 5 and 20 m/s values.
        if (windspd_2m.lt.5.0) then
          ro_nsnow_wind = 0.0
        else
          ro_nsnow_wind = 25.0 + &
     &      250.0 * (1.0 - exp(-(alfa*(windspd_2m - 5.0))))

! Now scale this 24-hour value by dt.  I checked, and this is the
!   perfect way to do this.  It makes the 3-hour and 24-hour dt
!   cases identical!
! In the end I decided to not do this scaling because the resulting
!   snow depth is scaled at this time step by the precipitation
!   inputs, and because I think the wind speed is the main
!   controlling factor, not the duration of the wind at this speed.
!         ro_nsnow_wind = dt / dt_24hours * ro_nsnow_wind
        endif

      else

        ro_nsnow_wind = 0.0

      endif

! If you want to turn this off, you can uncomment the following
!   line.
!     ro_nsnow_wind = 0.0

      return
      end SUBROUTINE NSNOW_DENSITY_FROM_BLOWING_SNOW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE WINDSPEED_2M(windspd,ht_windobs,windspd_2m)

      implicit none

      real ht_windobs,snow_z0,windspd,windspd_2m

! Required constants.
      snow_z0 = 0.001

! Calculate the 2-m wind speed.
      windspd_2m = windspd * &
     &  log(2.0/snow_z0)/log(ht_windobs/snow_z0)

      return
      end SUBROUTINE WINDSPEED_2M

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SUPERIMPOSED_ICE(iyear,output_path_wo_assim,&
     &  snod_col_init2,swed_col_init2,sden_col_init1,ro_snow)

! Here I am going to sweep through June and July and determine
!   whether the snowpack ever became isothermal during this
!   period.

      implicit none

      integer k,iend_year,iyear,i,j,ii,jj

! Number of tracked parcels.
      integer, parameter :: icount = 70000

      integer, parameter :: nnx=361,nny=361

      character*80 output_path_wo_assim

! Parcels.
      real roff_col(icount)
      real snod_col(icount)
      real swed_col(icount)
      real conc_col(icount)
      real simp_col(icount)
      real snod_col_init1(icount)
      real snod_col_init2(icount)
      real swed_col_init1(icount)
      real swed_col_init2(icount)
      real sden_col_init1(icount)

      real zi(icount)
      real zj(icount)
      real conc(icount)

! EASE grid.
      real snod1(nnx,nny)
      real snod2(nnx,nny)
      real swed1(nnx,nny)

      integer iyr_start,imo_start,idy_start,imo_end,idy_end,iiyear,&
     &  ioptn,julian_end,julian_start,irec,imo_beg,idy_beg,irec1,&
     &  irec2,irec3,julian_beg,i_len_wo,trailing_blanks,istart_year

      real undef,ro_snow

! Parcel tracks and concentration from NSIDC.
      character path3*(*) 
      parameter (path3= &
     &  '../../parcel_tracks/3_mk_daily/')

! SnowModel run outputs on original parcels.
!     character path1*(*) 
!     parameter (path1 = '/data4/lagrangian/sm_38yrs/merra/outputs/')
      i_len_wo = 80 - trailing_blanks(output_path_wo_assim)

      undef = -9999.0

      iend_year = 2018
      istart_year = 1981

      iiyear = iyear - istart_year + 1

! Open the parcel track data.
      if (iyear.eq.istart_year) then
        open (7021,file=path3//'ij_parcels_1980-2018.gdat', &
     &    form='unformatted',access='direct',recl=4*2*icount)

        open (7022,file=path3//'conc_parcels_1980-2018.gdat',&
     &    form='unformatted',access='direct',recl=4*1*icount)
      endif

! Input files.
!     open (7031,file=output_path_wo_assim(1:i_len_wo)//'roff.gdat',
!    &  form='unformatted',access='direct',recl=4*icount)

!     open (7032,file=output_path_wo_assim(1:i_len_wo)//'snod.gdat',
!    &  form='unformatted',access='direct',recl=4*icount)

!     open (7033,file=output_path_wo_assim(1:i_len_wo)//'swed.gdat',
!    &  form='unformatted',access='direct',recl=4*icount)

! These have to be closed and reopened so the data are completely
!   available for the following read statements.
      close (234)
      close (236)
      close (238)

      open (234,file=output_path_wo_assim(1:i_len_wo)//'roff.gdat',&
     &  form='unformatted',access='direct',recl=4*icount)

      open (236,file=output_path_wo_assim(1:i_len_wo)//'snod.gdat',&
     &  form='unformatted',access='direct',recl=4*icount)

      open (238,file=output_path_wo_assim(1:i_len_wo)//'swed.gdat',&
     &  form='unformatted',access='direct',recl=4*icount)

! Output files.
      if (iyear.eq.istart_year) then
        open (7051,file='seaice/superimposed_ice_flag_parcel.gdat',&
     &    form='unformatted',access='direct',recl=4*3*icount,&
     &    status='replace')

        open (7061,file='seaice/superimposed_ice_flag_ease.gdat',&
     &    form='unformatted',access='direct',recl=4*2*nnx*nny,&
     &    status='replace')

        open (7071,file='seaice/snod_init.gdat',&
     &    form='unformatted',access='direct',recl=4*3*icount,&
     &    status='replace')
      endif

! Find the irecs between 1 June and 31 July of each year.

! This is the begining of the simulation.
      iyr_start = 1980
      imo_start = 8
      idy_start = 1

! Julian start.
      ioptn = 3
      call calndr (ioptn,idy_start,imo_start,iyr_start,julian_start)

! This is the begining of the summer period.
      imo_beg = 6
      idy_beg = 1

! This is the end of the simulation year.
      imo_end = 7
      idy_end = 31

! Extracting the last day of the simulation year.

! Initialize the array.
      do k=1,icount
        simp_col(k) = undef
      enddo

! Find the irecs for 1 June and 31 July.
      call calndr (ioptn,idy_beg,imo_beg,iyear,julian_beg)
      irec1 = julian_beg - julian_start + 1
      call calndr (ioptn,idy_end,imo_end,iyear,julian_end)
      irec2 = julian_end - julian_start + 1

      print *,iiyear,iyear,irec1,irec2
      print *,iiyear,iyear,irec1,irec2
      print *,iiyear,iyear,irec1,irec2

      do irec=irec1,irec2

! Read in the parcel data.
!       read (7031,rec=irec) (roff_col(k),k=1,icount)
        read (234,rec=irec) (roff_col(k),k=1,icount)
        read (7022,rec=irec) (conc_col(k),k=1,icount)

! Perform some cleanup duties.  This should be done in the
!   SnowModel code before roff is written out.  But right now
!   it is not.
        do k=1,icount
          if (conc_col(k).eq.0.0) roff_col(k) = undef
        enddo

! These runoff data = 0.0 for parcels with snow (or no snow) and
!   zero runoff at this time step, and = undef if there is no
!   parcel (i.e., if conc = 0.0).

! If there is any runoff out the base of the snowpack during
!   June or July, then the snowpack is isothermal and we assume
!   that any remaining snow on this parcel on 31 July is so wet
!   that if and when it eventually re-freezes it will create
!   superimposed ice.
        do k=1,icount
          if (roff_col(k).gt.0.0) then
            simp_col(k) = 1.0
          endif
        enddo

      enddo

! Apply the 'ZERO THE SNOW DEPTH' requirement to the 31 July
!   snow depth.
!     read (7032,rec=irec2) (snod_col(k),k=1,icount)
!     read (7033,rec=irec2) (swed_col(k),k=1,icount)
      read (236,rec=irec2) (snod_col(k),k=1,icount)
      read (238,rec=irec2) (swed_col(k),k=1,icount)

      do k=1,icount
        snod_col_init1(k) = snod_col(k)
        swed_col_init1(k) = swed_col(k)
        if (simp_col(k).eq.1.0) then
          snod_col_init1(k) = 0.0
          swed_col_init1(k) = 0.0
        endif
      enddo

      write (7051,rec=iiyear) &
     &  (snod_col_init1(k),k=1,icount),&
     &  (snod_col(k),k=1,icount),&
     &  (simp_col(k),k=1,icount)

! Grid these data to the EASE grid.
      call parcel_to_ease(nnx,nny,icount,irec2,snod_col_init1,snod1,&
     &  undef)
      call parcel_to_ease(nnx,nny,icount,irec2,snod_col,snod2,&
     &  undef)
      call parcel_to_ease(nnx,nny,icount,irec2,swed_col_init1,swed1,&
     &  undef)

! Save the EASE grid data.
      write (7061,rec=iiyear)&
     &  ((snod1(i,j),i=1,nnx),j=1,nny),&
     &  ((snod2(i,j),i=1,nnx),j=1,nny)

! Prepare the initial condition snow depth array for 1 Aug.
!   snod1 is my initial condition on the EASE grid.  Use this to
!   define an initial condition on the 1 Aug parcels.  Do this
!   by finding the EASE grid (i,j) position that corresponds to
!   each 1 Aug parcel location.  And then assigning this initial
!   snow depth (and the other corresponding variables) to that
!   parcel.

! NOTE: this is going to crash on the last year!  I need to fix
!   that!
! The +1 here is what is converting the end of year data to a
!   begining of the year IC on the new parcels.
! This is a stupid way to fix that problem.
      if (iyear.ne.iend_year) then
        irec3 = irec2 + 1
      else
        irec3 = irec2
      endif

      read (7021,rec=irec3) (zi(k),k=1,icount),(zj(k),k=1,icount)
      read (7022,rec=irec3) (conc(k),k=1,icount)

      do k=1,icount
        snod_col_init2(k) = 0.0
        swed_col_init2(k) = 0.0
        if (conc(k).ne.0.0) then
          ii = nint(zi(k))
          jj = nint(zj(k))
          if (snod1(ii,jj).ne.undef) then
            snod_col_init2(k) = snod1(ii,jj)
            swed_col_init2(k) = swed1(ii,jj)
          else
            snod_col_init2(k) = 0.0
            swed_col_init2(k) = 0.0
          endif
        endif

        if (snod_col_init2(k).eq.undef) then
          print *,'snod_col_init2(k)=undef',k,conc(k)
        endif
        if (swed_col_init2(k).eq.undef) then
          print *,'swed_col_init2(k)=undef',k,conc(k)
        endif

      enddo

! Calculate the snow density IC by using the swed and snod values.
!   This cleans up the problem of having undef sden when the
!   depths are 0.0.
      do k=1,icount
        if (snod_col_init2(k).ne.0.0 .and. &
     &    snod_col_init2(k).ne.undef) then
          sden_col_init1(k) = 1000.0 * swed_col_init2(k) / &
     &      snod_col_init2(k)
        else
          sden_col_init1(k) = ro_snow
        endif
      enddo

! Save the snow depth initial condition data.
      write (7071,rec=iiyear) (snod_col_init2(k),k=1,icount),&
     &                        (swed_col_init2(k),k=1,icount),&
     &                        (sden_col_init1(k),k=1,icount)

!     close (7021)
!     close (7022)
!     close (7031)
!     close (7032)
!     close (7033)
!     close (7051)
!     close (7061)
!     close (7071)

      return
      end SUBROUTINE SUPERIMPOSED_ICE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine parcel_to_ease(nnx,nny,icount,iter,snow_col,snow,&
     &  undef)

      implicit none

      integer ii,jj,nnx,nny,i,j,k,iter,icount

      real zi(icount),zj(icount),conc(icount)

      real snow_col(icount)
      real snow(nnx,nny)

      real cnt(nnx,nny)
      real conc_area_sum(nnx,nny)

      real undef,zii,zjj
      real center_x1,center_y1,center_x2,center_y2,area

! Read in the parcel track data.
      read (7021,rec=iter) (zi(k),k=1,icount),(zj(k),k=1,icount)
      read (7022,rec=iter) (conc(k),k=1,icount)

! Make any 0.0 concentrations undefined, and convert them
!   from 0-100 to 0-1.
      do k=1,icount
        conc(k) = conc(k) / 100.0
        if (conc(k).eq.0.0) conc(k) = undef
      enddo

! Initialize the working arrays.
      do j=1,nny
        do i=1,nnx

! Snow depth.
          snow(i,j) = 0.0

! Counting array.
          cnt(i,j) = 0.0

! Concentration-area product summing array.
          conc_area_sum(i,j) = 0.0

        enddo
      enddo

      do k=1,icount
        if (conc(k).ne.undef) then

! Decimal coordinate of the center of the parcel.
          center_x2 = zi(k)
          center_y2 = zj(k)

! Coordinate of the EASE grid cell center that this parcel
!   center sits in, integer and real versions.
          ii = nint(zi(k))
          jj = nint(zj(k))

! This eliminates any problems at the boundaries, if you have
!   a situation where there is non-zero conc on the boundary
!   grid cells.
          if (ii.eq.1) ii = 2
          if (ii.eq.nnx) ii = nnx-1
          if (jj.eq.1) jj = 2
          if (jj.eq.nny) jj = nny-1

          zii = real(ii)
          zjj = real(jj)

! Loop through the 3x3 grid cells surrounding the EASE grid cell.
          do j=jj-1,jj+1
            do i=ii-1,ii+1

! Define the center coords of the surrounding EASE grid cells.
              center_x1 = real(i)
              center_y1 = real(j)

! Find the overlap in fractional area of each EASE grid cell.  So,
!   the "area" coming out of this is: 0-1 = the fraction of the
!   25-km by 25-km EASE grid cell that this parcel covers.
              call overlap_area (center_x1,center_y1,center_x2,&
     &          center_y2,area)

! Use the concentration to weight the snow depth under the
!   assumption that the concentration is proportional to the
!   contributing area.
              area = area * conc(k)

! Count how many parcels contribute to each EASE grid cell.
              if (area.gt.0.0) then
                cnt(i,j) = cnt(i,j) + 1.0
              endif

! Keep track of the concentrations that sum to over 1.0.  They
!   will be used to scale the result at the end of the averaging
!   calculations.
              conc_area_sum(i,j) = conc_area_sum(i,j) + area

! Sum the snow depth contributions from each parcel in each EASE
!   grid cell, by scaling the depth values by the area and conc.
              snow(i,j) = snow(i,j) + area * snow_col(k)

            enddo
          enddo

        endif
      enddo

! Take care of the problem where the area-concentrations summed
!   to over 1.0.
      do j=1,nny
        do i=1,nnx

! Clip the values below 1.0 so when you do the divide, only
!   concentration sums above 1.0 scale the result (values below
!   1.0 are already grid-averaged values).
          conc_area_sum(i,j) = max(1.0,conc_area_sum(i,j))

! Scale the snow depth averages for cases where the concentration
!   summed to greater than 1.0.
          snow(i,j) = snow(i,j) / conc_area_sum(i,j)

        enddo
      enddo

! Clean up any grid cells that never had any ice in them.
      do j=1,nny
        do i=1,nnx
          if (cnt(i,j).eq.0.0) then
            snow(i,j) = undef
          endif
        enddo
      enddo

      return
      end SUBROUTINE parcel_to_ease

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine overlap_area (center_x1,center_y1,center_x2,&
     & center_y2,area)

! Here center_x1,y1 is the EASE grid cell, and center_x2,y2 is
!   the parcel cell.

      implicit none

      real center_x1,center_y1
      real center_x2,center_y2
      real x1_min,x1_max,y1_min,y1_max
      real x2_min,x2_max,y2_min,y2_max
      real x_overlap,y_overlap,area

! These are the boundary coords of the EASE grid cell.
      x1_min = center_x1 - 0.5
      x1_max = center_x1 + 0.5
      y1_min = center_y1 - 0.5
      y1_max = center_y1 + 0.5

! These are the boundary coords of the parcel.
      x2_min = center_x2 - 0.5
      x2_max = center_x2 + 0.5
      y2_min = center_y2 - 0.5
      y2_max = center_y2 + 0.5

! Define the overlap in x and y.
      x_overlap = max(0.0,min(x1_max,x2_max) - max(x1_min,x2_min))
      y_overlap = max(0.0,min(y1_max,y2_max) - max(y1_min,y2_min))

! Calculate the area overlap.
      area = x_overlap * y_overlap

!     print 98, x_overlap,y_overlap,area
!  98 format (3f10.3)

      return
      end subroutine overlap_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
