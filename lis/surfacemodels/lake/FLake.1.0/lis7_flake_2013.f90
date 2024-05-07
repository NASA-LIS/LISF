!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!BOP
!
! !ROUTINE: LIS7_FLake_2013
! \label{LIS7_FLake_2013}
!
! !REVISION HISTORY:
!  Apr 2013: Shugong Wang intial implementation
!
! !INTERFACE:
subroutine LIS7_FLake_2013(t, swd, lwd, wind_e, wind_n, tmp, q2, psurf,    & ! forcing data 
                           height_wind, height_tq,                      & ! forcing forcing parameter
                           flake_dt, lon, lat,                          & ! model time step, longitude and latitude
                           depth_w, fetch, depth_bs, T_bs,              & ! lake parameters
                           T_snow, T_ice, T_mnw, T_wML, T_bot, T_b1,    & ! state variables 
                           C_T, H_snow, H_ice, H_ML, H_B1, T_sfc,       & ! state variables 
                           albedo_water, albedo_ice, albedo_snow,       & ! state variables (model-driven time-varying parameters)
                           ufr_a, ufr_w, Wconv, Q_se, Q_la, I_w,        & ! output variables 
                           Q_lwa, Q_lww, Q_bot)                           ! output variables 


! ! USES:

!
! ! DESCRIPTION:
!   A standard interface defined for FLake model
!
!EOP
    use data_parameters, only    : ireals, iintegers
    use SfcFlx, only             : u_star_a_sf,     &  ! friction velocity [m/s]
                                   Q_sens_a_sf,     &  ! sensible heat flux [W/m2]
                                   Q_lat_a_sf          ! latent heat flux [W/m2]
    use flake, only              : u_star_w_flk,    &  ! friction velocity in the surface layer of later water [m/s]
                                   T_bot_n_flk,     &  ! temperature at the water-bottom sediment interface [K]             
                                   I_w_flk,         &  ! radiation flux through the ice-water or air-water interface [W m^{-2}]
                                   Q_bot_flk,       &  ! heat flux through the water-bottom sediment interface [W m^{-2}]
                                   w_star_sfc_flk      !  Convective velocity scale, using a generalized heat flux scale [m s^{-1}]
    use SfcFlx, only             : SfcFlx_lwradwsfc    ! function, returns the surface longwave radiation from water surface to atmosphere
    use flake_derivedtypes
                               
    implicit none
    integer, intent(in) :: t
    !!! interface variables 
    real, intent(in)    :: swd             ! incoming solar radiation at the surface [W/m2]
    real, intent(in)    :: lwd             ! incoming longwave radiation at the surface [W/m2]
    real, intent(in)    :: wind_e           ! zonal wind speed (x-coordinate, along latitude) at height_wind [m/s] 
    real, intent(in)    :: wind_n           ! meridional wind speed (y-coordinate, along longitude) at height_wind [m/s] 
    real, intent(in)    :: tmp             ! air temperature at height_tq [K]
    real, intent(in)    :: q2              ! specific air humidity at height_tq [kg/kg]
    real, intent(in)    :: psurf           ! surface air pressure [Pa]
    real, intent(in)    :: height_wind     ! height where wind speed is measured [m]
    real, intent(in)    :: height_tq       ! height where temperature and humidity are measured
    real, intent(in)    :: flake_dt        ! model time step [s]
    real, intent(in)    :: lon             ! longitude [-]
    real, intent(in)    :: lat             ! latitude [-]
    real, intent(in)    :: depth_w         ! lake depth [m]
    real, intent(in)    :: fetch           ! typical wind fetch [m]
    real, intent(in)    :: depth_bs        ! depth of the thermally active layer of the bottom sediments [m]
    real, intent(in)    :: T_bs            ! temperature at the outer edge of the thermally active layer of the bottom sediments [K]
    real, intent(inout) :: T_snow          ! temperature at the air-snow interface [K] 
    real, intent(inout) :: T_ice           ! temperature at the snow-ice interface [K]
    real, intent(inout) :: T_mnw           ! mean temperature of the water column [K]
    real, intent(inout) :: T_wML           ! temperature of mixed layer [K]
    real, intent(inout) :: T_bot           ! temperature at the water-bottom sediment interface [K]
    real, intent(inout) :: T_b1            ! temperature at the bottom of the upper layer of the sediments [K]
    real, intent(inout) :: C_T             ! shape factor (thermocline) [-]
    real, intent(inout) :: H_snow          ! snow thickness [m]
    real, intent(inout) :: H_ice           ! ice thickness [m] 
    real, intent(inout) :: H_ML            ! thickness of mixed layer [m] 
    real, intent(inout) :: H_B1            ! thickness of the upper layer of bottom sediments [m]
    real, intent(inout) :: T_sfc           ! surface temperature [K]
    real, intent(inout) :: albedo_water    ! water surface albedo with resect to solar radiation [-] 
    real, intent(inout) :: albedo_ice      ! ice surface albedo with respect to the solar radiation [-]
    real, intent(inout) :: albedo_snow     ! snow surface albedo with respect to the solar radiation [-]
    real, intent(out)   :: ufr_a           ! friction velocity in air [m/s]
    real, intent(out)   :: ufr_w           ! friction velocity in surface water [m/s]
    real, intent(out)   :: Wconv           ! convective velocity scale [m/s]
    real, intent(out)   :: Q_se            ! sensible surface heat flux [W/m2]
    real, intent(out)   :: Q_la            ! latent surface heat flux [W/m2] 
    real, intent(out)   :: I_w             ! shortwave radiation [W/m2]
    real, intent(out)   :: Q_lwa           ! longwave radiation flux from atmosphere [W/m2]
    real, intent(out)   :: Q_lww           ! longwave radiation flux from water [W/m2]
    real, intent(out)   :: Q_bot           ! heat flux across water-sediments boundary [W/m2]
    !!! local variables
    type(opticpar_medium) :: opticpar_water_inout  ! Optical characteristics of water, dummy parameter 
    type(opticpar_medium) :: opticpar_ice_inout    ! Optical characteristics of ice,   dummy parameter
    type(opticpar_medium) :: opticpar_snow_inout   ! Optical characteristics of snow,  dummy parameter
    real(kind=ireals)   :: dMsnowdt_in     , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
                           I_atm_in        , & ! Solar radiation flux at the surface [W m^{-2}]
                           Q_atm_lw_in     , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
                           height_u_in     , & ! Height above the lake surface where the wind speed is measured [m]
                           height_tq_in    , & ! Height where temperature and humidity are measured [m]
                           U_a_in          , & ! Wind speed at z=height_u_in [m s^{-1}]
                           T_a_in          , & ! Air temperature at z=height_tq_in [K]
                           q_a_in          , & ! Air specific humidity at z=height_tq_in
                           P_a_in          , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]
                           depth_w_in      , & ! The lake depth [m]
                           fetch_in        , & ! Typical wind fetch [m]
                           depth_bs_in     , & ! Depth of the thermally active layer of the bottom sediments [m]
                           T_bs_in         , & ! Temperature at the outer edge of the thermally active layer of the bottom sediments [K]
                           par_Coriolis_in , & ! The Coriolis parameter [s^{-1}]
                           del_time_in     , & ! The model time step [s]
                           T_snow_in       , & ! Temperature at the air-snow interface [K] 
                           T_ice_in        , & ! Temperature at the snow-ice or air-ice interface [K]
                           T_mnw_in        , & ! Mean temperature of the water column [K]
                           T_wML_in        , & ! Mixed-layer temperature [K]
                           T_bot_in        , & ! Temperature at the water-bottom sediment interface [K]
                           T_B1_in         , & ! Temperature at the bottom of the upper layer of the sediments [K]
                           C_T_in          , & ! Shape factor (thermocline)
                           H_snow_in       , & ! Snow thickness [m]
                           H_ice_in        , & ! Ice thickness [m]
                           H_ML_in         , & ! Thickness of the mixed-layer [m]
                           H_B1_in         , & ! Thickness of the upper layer of bottom sediments [m]
                           T_sfc_p_in      , & ! Surface temperature at the previous time step [K] 
                           albedo_water_inout    , & ! Water surface albedo with respect to the solar radiation
                           albedo_ice_inout      , & ! Ice surface albedo with respect to the solar radiation
                           albedo_snow_inout     , & ! Snow surface albedo with respect to the solar radiation
                           T_snow_out      , & ! Temperature at the air-snow interface [K] 
                           T_ice_out       , & ! Temperature at the snow-ice or air-ice interface [K]
                           T_mnw_out       , & ! Mean temperature of the water column [K]
                           T_wML_out       , & ! Mixed-layer temperature [K]
                           T_bot_out       , & ! Temperature at the water-bottom sediment interface [K]
                           T_B1_out        , & ! Temperature at the bottom of the upper layer of the sediments [K]
                           C_T_out         , & ! Shape factor (thermocline)
                           H_snow_out      , & ! Snow thickness [m]
                           H_ice_out       , & ! Ice thickness [m]
                           H_ML_out        , & ! Thickness of the mixed-layer [m]
                           H_B1_out        , & ! Thickness of the upper layer of bottom sediments [m]
                           T_sfc_n_out         ! Updated surface temperature [K] 
    !!! 
    real, parameter :: Omega  = 7.27e-5   ! angular speed of earth rad/s 
    real :: PI



    !!! Set the rate of snow accumulation to be 0 because the snow module is not used in this version. 
    !!! http://www.flake.igb-berlin.de/docs.shtml : The snow module of FLake is not used at the time being. 
    !!! No switch is required to deactivate the snow module. All one has to do is to pass a zero value of 
    !!! dMsnowdt_in to FLAKE_INTERFACE. This choice can be hard-coded by explicitly setting dMsnowdt_flk to 
    !!! zero in FLAKE_INTERFACE.
    dMsnowdt_in = 0.0        
    
    !!! solar radiation flux at the surface [W/m2]
    I_atm_in     = real(swd, kind=ireals)

    !!! longwave radiation flux from the atmosphere [W/m2]
    Q_atm_lw_in  = real(lwd, kind=ireals)

    !!! wind speed at z=height_u_in [m/s]

    U_a_in       = real(sqrt(wind_e*wind_e + wind_n*wind_n), kind=ireals)

    !!! air temperature at z=height_tq_in [K]
    T_a_in       = real(tmp, kind=ireals)

    !!! air specific humidity at height_tq_in [kg/kg]
    q_a_in       = real(q2, kind=ireals)

    !!! surface air pressure [Pa=N/m2]
    P_a_in       = real(psurf, kind=ireals)

    !!! height above the lake surface where the wind speed is measured [m]
    height_u_in  = real(height_wind, kind=ireals)

    !!! height above the lake surface where the temperature and humidity are measured [m]
    height_tq_in = real(height_tq, kind=ireals)

    !!!!! Lake Parameters 
    !!! lake depth [m]
    depth_w_in   = real(depth_w, kind=ireals)
    
    !!! typical wind fetch [m]
    fetch_in     = real(fetch, kind=ireals)
    
    !!! depth of the thermally active layer of the bottom sediments [m]
    depth_bs_in  = real(depth_bs, kind=ireals)

    !!! temperature at the outer edge of the thermally active layer of the bottom sediments [K]
    T_bs_in      = real(T_bs, kind=ireals)

    !!! The Coriolis parameter [s^-1]
    PI = 4.0 * atan(1.0)
    par_Coriolis_in = 2 * Omega * sin(PI * lat/180.0)
    
    !!! model time step [s]
    del_time_in   = real(flake_dt, kind=ireals)

    !!!! state variables 
    T_snow_in     = real(T_snow, kind=ireals)
    T_ice_in      = real(T_ice,  kind=ireals)
    T_mnw_in      = real(T_mnw,  kind=ireals)
    T_wML_in      = real(T_wML,  kind=ireals)
    T_bot_in      = real(T_bot,  kind=ireals)
    T_B1_in       = real(T_B1,   kind=ireals)
    C_T_in        = real(C_T,    kind=ireals)
    H_snow_in     = real(H_snow, kind=ireals)
    H_ice_in      = real(H_ice,  kind=ireals)
    H_ML_in       = real(H_ML,   kind=ireals)
    H_B1_in       = real(h_B1,   kind=ireals)
    T_sfc_p_in    = real(T_sfc,  kind=ireals)
    !!! model-driven dynamic parameters 
    albedo_water_inout    = real(albedo_water, kind=ireals)
    albedo_ice_inout      = real(albedo_ice,   kind=ireals)
    albedo_snow_inout     = real(albedo_snow,  kind=ireals)
   
    !!! call FLake model physics 
    call flake_interface ( lat, dMsnowdt_in, I_atm_in, Q_atm_lw_in, height_u_in, height_tq_in,     &  ! forcing 
                           U_a_in, T_a_in, q_a_in, P_a_in,                                    &  ! forcing  
                           depth_w_in, fetch_in, depth_bs_in, T_bs_in,                        &  ! parameters
                           par_Coriolis_in, del_time_in,                                      &  ! parameters 
                           T_snow_in,  T_ice_in,  T_mnw_in,  T_wML_in,  T_bot_in,  T_B1_in,   &  ! state variables  
                           C_T_in,  H_snow_in,  H_ice_in,  H_ML_in,  H_B1_in, T_sfc_p_in,     &  ! state variables  
                           albedo_water_inout, albedo_ice_inout, albedo_snow_inout,           &  ! model-driven dynamic parameters  
                           opticpar_water_inout, opticpar_ice_inout, opticpar_snow_inout,     &  ! model-driven dynamic parameters 
                           T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out, T_B1_out,  &  ! model output variables 
                           C_T_out, H_snow_out, H_ice_out, H_ML_out, H_B1_out, T_sfc_n_out )     ! model output variables 
        
    !!! save state variables
    T_snow  = real(T_snow_out)
    T_ice   = real(T_ice_out)
    T_mnw   = real(T_mnw_out)
    T_wML   = real(T_wML_out)
    T_bot   = real(T_bot_out)
    T_B1    = real(T_B1_out)
    C_T     = real(C_T_out)
    H_snow  = real(H_snow_out)
    H_ice   = real(H_ice_out)
    H_ML    = real(H_ML_out)
    H_B1    = real(H_B1_out)
    T_sfc   = real(T_sfc_n_out)

    !!! save model-driven dynamic parameters 
    albedo_water    = real(albedo_water_inout)
    albedo_ice      = real(albedo_ice_inout)
    albedo_snow     = real(albedo_snow_inout)

    !!! other output variables
    ufr_a   = real(u_star_a_sf)
    ufr_w   = real(u_star_w_flk)
    Wconv   = real(w_star_sfc_flk)
    Q_se    = real(Q_sens_a_sf)
    Q_la    = real(Q_lat_a_sf)
    I_w     = real(I_w_flk)
    Q_bot   = real(Q_bot_flk)
    Q_lwa   = lwd                  ! forcing 
    Q_lww   = real(SfcFlx_lwradwsfc(T_sfc_n_out))

end subroutine LIS7_FLake_2013    

