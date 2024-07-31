!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!#if !defined(LIS_JULES)
SUBROUTINE lis_control ()

!-------------------------------------------------------------------------------
! Control level routine, to call the main parts of the model.
!-------------------------------------------------------------------------------

!Import the main subroutines we will use
USE surf_couple_radiation_mod, ONLY: surf_couple_radiation
USE surf_couple_explicit_mod,  ONLY: surf_couple_explicit
USE surf_couple_implicit_mod,  ONLY: surf_couple_implicit
USE surf_couple_extra_mod,     ONLY: surf_couple_extra
USE update_mod,                ONLY: calc_downward_rad, & 
                                     update_derived_variables
use model_time_mod,            only: timestep
USE forcing,                  ONLY:                                         &
  con_rain_ij, con_snow_ij, ls_rain_ij, ls_snow_ij, lw_down_ij, pstar_ij,   &
  qw_1_ij, sw_down_ij, tl_1_ij, u_0_ij, u_1_ij, v_0_ij, v_1_ij, diff_rad_ij
!Variables- modules in alphabetical order

USE aero, ONLY:                                                               &
  co2_3d_ij, rho_aresist_ij, aresist_ij, resist_b_ij, rho_aresist_surft,      &
  aresist_surft, resist_b_surft, r_b_dust_ij, cd_std_dust_ij, u_s_std_surft

USE ancil_info,               ONLY:                                           &
  ice_fract_ij, ice_fract_ncat_sicat, pond_frac_cat_sicat,                    &
  pond_depth_cat_sicat, land_pts, land_index, nsurft, surft_pts,              &
  surft_index, frac_surft, row_length, rows, sstfrz_ij,  z1_tq_ij,            &
  dim_cs1, dim_cs2, z1_uv_ij, ti_cat_sicat, land_pts_trif, lice_index,        &
  lice_pts, soil_index, soil_pts, nsoilt

USE atm_fields_bounds_mod,    ONLY:                                           &
  tdims, pdims, udims, vdims, tdims_s, pdims_s, udims_s, vdims_s

USE bl_option_mod,            ONLY: on

USE boundary_layer_var_mod,   ONLY:                                           &
  zh

USE c_elevate,                ONLY: z_land_ij

USE coastal,                  ONLY:                                           &
  flandg, taux_land_ij, taux_ssi_ij, tauy_land_ij, tauy_ssi_ij,               &
  tstar_sice_sicat, fland, tstar_sea_ij,  tstar_ssi_ij, surf_ht_flux_land_ij, &
  surf_ht_flux_sice_sicat, tstar_land_ij, tstar_sice_ij,                      &
  taux_land_star, tauy_land_star, taux_ssi_star, tauy_ssi_star, vshr_land_ij, &
  vshr_ssi_ij

USE datetime_mod,             ONLY: l_360, l_leap

USE datetime_utils_mod,       ONLY: day_of_year

USE diag_swchs,               ONLY: stf_sub_surf_roff

USE fluxes,                   ONLY:                                           &
  alb_surft, e_sea_ij, ecan_ij, ecan_surft, ei_ij, ei_surft, esoil_ij_soilt,  &
  esoil_surft, ext_soilt, fqw_surft, fqw_sicat, fsmc_pft,                     &
  ftl_sicat, ftl_surft, h_sea_ij, hf_snow_melt_gb, land_albedo_ij,            &
  le_surft, melt_surft, sea_ice_htf_sicat, snomlt_sub_htf_gb, snow_melt_gb,   &
  snowmelt_ij, sub_surf_roff_gb, surf_ht_flux_ij, surf_htf_surft,             &
  surf_roff_gb, radnet_surft, tot_tfall_gb, tstar_ij,                         &
  sw_surft, emis_surft, rflow_gb, rrun_gb, snow_soil_htf, z0m_surft, z0h_surft
USE gridmean_fluxes,          ONLY : fqw_1_ij, ftl_1_ij,taux_1_ij, tauy_1_ij 
USE jules_radiation_mod,      ONLY: i_sea_alb_method, l_cosz

USE jules_sea_seaice_mod,     ONLY: l_ctile, buddy_sea, nice_use

USE jules_soil_mod,           ONLY: sm_levels

USE jules_surface_types_mod,  ONLY: ntype

USE model_time_mod,           ONLY: current_time

USE orog,                     ONLY:                                           &
  ho2r2_orog_gb, sil_orog_land_gb, h_blend_orog_ij, z0m_eff_ij

USE planet_constants_mod,     ONLY: c_virtual

USE prognostics,              ONLY:                                           &
  snow_mass_ij, snow_mass_sea_sicat, di_ncat_sicat, lai_pft, canht_pft,       &
  rgrain_surft, snow_surft, soot_ij, tstar_surft, canopy_surft, smc_soilt,    &
  k_sice_sicat, t_soil_soilt, ti_sicat, tsurf_elev_surft, z0msea_ij, gs_gb,   &
  gc_surft, canopy_gb, smcl_soilt, snow_grnd_surft

USE p_s_parms,                ONLY:                                           &
  cosz_ij, albsoil_soilt, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb,         &
  z0_surft, catch_surft, catch_snow_surft, hcon_soilt, smvccl_soilt,          &
  smvcst_soilt, smvcwt_soilt, sthf_soilt, sthu_soilt, z0h_bare_surft,         &
  z0m_soil_gb, infil_surft

USE sf_diags_mod,             ONLY: sf_diag

USE switches,                 ONLY:                                           &
  lq_mix_bl, l_co2_interactive, l_spec_z0

USE theta_field_sizes,        ONLY: t_i_length, t_j_length

USE top_pdm,                  ONLY:                                           &
  a_fsat_soilt, c_fsat_soilt, a_fwet_soilt, c_fwet_soilt, fexp_soilt,         &
  fsat_soilt, fwetl_soilt, gamtot_soilt, sthzw_soilt, ti_mean_soilt,          &
  ti_sig_soilt, zw_soilt

USE trifctl,                  ONLY:                                           &
  asteps_since_triffid, g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,            &
  gpp_gb, npp_gb, resp_p_gb, g_leaf_pft, gpp_pft, npp_pft, resp_p_pft,        &
  resp_s_soilt, resp_w_pft, g_leaf_phen_acc_pft, frac_agr_gb

USE u_v_grid,                 ONLY:                                           &
   dtrdz_charney_grid_1_ij, u_0_p_ij, u_1_p_ij, v_0_p_ij, v_1_p_ij

USE zenith_mod,               ONLY: zenith

!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
!INTEGER :: timestep    ! IN Atmospheric timestep number.

!!-------------------------------------------------------------------------------
!! Arguments with intent(in)
!!-------------------------------------------------------------------------------
!!Forcing INTENT(IN)
!REAL, INTENT(IN) ::                                                           &
!  qw_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
!!    Total water content (Kg/Kg)
!  tl_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
!!    Ice/liquid water temperature (k)
!  u_0_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!!    W'ly component of surface current (m/s)
!  v_0_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!!    S'ly component of surface current (m/s)
!  u_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!!    W'ly wind component (m/s)
!  v_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!!    S'ly wind component (m/s)
!  pstar_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
!!    Surface pressure (Pascals)
!  ls_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
!!    Large-scale rain (kg/m2/s)
!  con_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
!!     Convective rain (kg/m2/s)
!  ls_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
!!    Large-scale snowfall (kg/m2/s)
!  con_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
!!    Convective snowfall (kg/m2/s)
!  sw_down_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
!!    Surface downward SW radiation (W/m2)
!  lw_down_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!!    Surface downward LW radiation (W/m2)


!!-------------------------------------------------------------------------------
!! Arguments with intent(out)
!!-------------------------------------------------------------------------------
!!Fluxes INTENT(OUT)
!REAL, INTENT(OUT) ::                                                          &
!  fqw_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
!!    Moisture flux between layers (kg per square metre per sec).
!  ftl_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
!!    FTL(,K) contains net turbulent sensible heat flux into layer K from below;
!!    so FTL(,1) is the surface sensible heat, H.(W/m2)
!  taux_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
!!    W'ly component of surface wind stress (N/sq m)
!  tauy_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!!    S'ly component of surface wind stress (N/sq m)
!!    On V-grid; comments as per TAUX

!-------------------------------------------------------------------------------
!Parameters
!-------------------------------------------------------------------------------

LOGICAL, PARAMETER ::                                                         &
  l_aero_classic = .FALSE.
    !switch for CLASSIC aerosol- NEVER USED!

INTEGER, PARAMETER ::                                                         &
  n_swbands = 1,                                                              &
    !  Nnumber of SW bands
  numcycles = 1,                                                              &
    !Number of cycles (iterations) for iterative SISL.
  cycleno = 1,                                                                &
    !Iteration no
  n_proc = 2,                                                                 &
  river_row_length = 1,                                                       &
  river_rows = 1,                                                             &
  aocpl_row_length = 1,                                                       &
  aocpl_p_rows = 1

REAL, PARAMETER ::                                                            &
  r_gamma = 2.0,                                                              &
    ! Implicit weighting coefficient
  spec_band_bb(n_swbands) = -1.0
    !spectral band boundary (bb=-1)

!-------------------------------------------------------------------------------
! Local variables
! In the vast majority of cases, these are arguments only required when 
! coupled to the UM.
!-------------------------------------------------------------------------------

LOGICAL ::                                                                    &
  l_correct,                                                                  &
  at_extremity(4),                                                            &
  land_sea_mask(row_length,rows)

INTEGER ::                                                                    &
  a_steps_since_riv,                                                          &
  g_p_field,                                                                  &
  g_r_field,                                                                  &
  global_row_length,                                                          &
  global_rows,                                                                &
  global_river_row_length,                                                    &
  global_river_rows,                                                          &
  halo_i,                                                                     &
  halo_j,                                                                     &
  model_levels,                                                               &
  n_rows,                                                                     &
  offx,                                                                       &
  offy,                                                                       &
  n_procx,                                                                    &
  n_procy,                                                                    &
  g_rows (0:n_proc-1),                                                        &
  g_row_length (0:n_proc-1),                                                  &
  error,                                                                      &
    ! OUT 0 - AOK; 1 to 7  - bad grid definition detected
  curr_year,                                                                  &
  curr_day_number,                                                            &
  curr_hour,                                                                  &
  curr_minute,                                                                &
  curr_second

REAL ::                                                                       &
  olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !TOA-surface upward LW on last radiation timestep Corrected TOA outward LW
  sea_ice_albedo(row_length,rows,4),                                          &
    ! Sea ice albedo
  ws10m(t_i_length,t_j_length),                                               &
    ! 10m wind speed (m s-1)
  photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),    &
    ! Net downward shortwave radiation in band 1 (w/m2).
  bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    ! A buoyancy parameter (beta T tilde).
  bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    ! A buoyancy parameter (beta q tilde).
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    ! Surface net radiation on sea-ice (W/m2)
  rhokm_land(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),    &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),     &
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  alpha1(land_pts,nsurft),                                                    &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    ! ALPHA1 for open sea.
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    ! ALPHA1 for sea-ice.
  ashtf_prime(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    ! Adjusted SEB coefficient for sea-ice
  ashtf_prime_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
    ! Adjusted SEB coefficient for open sea
  ashtf_prime_surft(land_pts,nsurft),                                         &
    ! Adjusted SEB coefficient for land tiles
  rhokm_1(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
    ! Exchange coefficients for momentum on P-grid
  epot_surft(land_pts,nsurft),                                                &
    ! Local EPOT for land tiles.
  fraca(land_pts,nsurft),                                                     &
    !Fraction of surface moisture flux with only aerodynamic resistance for
    !snow-free land tiles.
  resfs(land_pts,nsurft),                                                     &
    ! Combined soil, stomatal and aerodynamic resistance factor for fraction
    !(1-FRACA) of snow-free land tiles.
  resft(land_pts,nsurft),                                                     &
    !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1
    !for snow.
  rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    ! Grid-box surface exchange coefficients
  rhokh_surft(land_pts,nsurft),                                               &
    ! Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),   &
    ! Surface exchange coefficients for sea-ice
  rhokh_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    ! Surface exchange coefficients for sea
  dtstar_ij_surft(land_pts,nsurft),                                           &
    ! Change in tstar_ij over timestep for land tiles
  dtstar_ij_sea(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end),                                   &
    ! Change is tstar_ij over timestep for open sea
  dtstar_ij_sice(tdims%i_start:tdims%i_end,                                   &
                tdims%j_start:tdims%j_end,nice_use),                          &
    ! Change is tstar_ij over timestep for sea-ice
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    ! Roughness length for heat and moisture over sea (m).
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    ! Roughness length for momentum over sea (m).
  chr1p5m(land_pts,nsurft),                                                   &
    ! Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    ! CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                               &
    ! Areal heat capacity of canopy for land tiles (J/K/m2).
  wt_ext_surft(land_pts,sm_levels,nsurft),                                    &
    !Fraction of evapotranspiration which is extracted from each soil layer
    !by each tile.
  flake(land_pts,nsurft),                                                     &
    !Lake fraction.
  tile_frac(land_pts,nsurft),                                                 &
    !Tile fractions including snow cover in the ice tile.
  hcons_soilt(land_pts,nsoilt),                                               &
    !Thermal conductivity of top soil layer, including water and ice (W/m/K)
  flandfac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),            &
  flandfac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),            &
  fseafac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),             &
  fseafac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),             &
  du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),            &
    !Level 1 increment to u wind field
  dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),            &
    !Level 1 increment to v wind field
  r_gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
    !weights for new BL solver
  r_gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
  chloro(row_length,rows),                                                    &
    !nr surface chlorophyll content

!Radiation
  albobs_sc_ij(t_i_length,t_j_length,ntype,2),                                &
    !albedo scaling factors to obs
  open_sea_albedo(row_length,rows,2,n_swbands),                               &
    !Surface albedo for Open Sea (direct and diffuse components, for each
    !band, with zeros for safety where no value applies)

!Implicit
  ctctq1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  dqw1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  dtl1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  du_star1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),      &
  dv_star1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),      &
  cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end),             &
    ! Coefficient in U tri-diagonal implicit matrix
  cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),             &
    ! Coefficient in V tri-diagonal implicit matrix
  flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),              &
    !Land frac (on U-grid, with 1st and last rows undefined or, at present,
    !set to "missing data")
  flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),              &
    !Land frac (on V-grid, with 1st and last rows undefined or, at present,
    !set to "missing data")
  rho1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Density on lowest level
  f3_at_p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Coriolis parameter
  uStargbm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
    ! BM surface friction velocity
  tscrndcl_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !Decoupled screen-level temperature over sea or sea-ice
  tscrndcl_surft(land_pts,nsurft),                                            &
    !Decoupled screen-level temperature over land tiles
  tstbtrans(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Time since the transition
  ei_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Sea ice sublimation
  rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Exchange coeffs for moisture.
  ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &

!Explicit
  !These variables are INTENT(IN) to sf_expl, but not used with the
  !current configuration of standalone JULES (initialised to 0 below)
  z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest uv-layer
  z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest Tq-layer
  ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    ! Convective downdraught mass-flux at cloud base
  cs_pool_gb_um(land_pts,dim_cs1),                                            &
  resp_s_acc_gb_um(land_pts_trif,dim_cs1),                                    &
  resp_s_gb_um(land_pts,dim_cs1),                                             &
  soil_clay_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  ! Charnock parameter from the wave model
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &

  !These variables are required for prescribed roughness lengths in
  !SCM mode in UM - not used standalone
  z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for momentum (m).(SCM)
  z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for heat (m). (SCM)
  recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
    !Reciprocal of the surface Obukhov  length at sea points. (m-1).
  rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Mean bulk Richardson number for lowest layer.
  flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
  fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface flux buoyancy over density (m^2/s^3)
  u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Surface friction velocity (m/s)
  t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent fluctuations of layer 1 temp; used in
    !initiating convection.
  q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent flux of layer 1 humidity; used in
    !initiating convection.
  rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface air density
  vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Magnitude of surface-to-lowest atm level wind shear (m per s).
  resp_s_tot_soilt(dim_cs2,nsoilt),                                           &
    !Total soil respiration (kg C/m2/s).
  emis_soil(land_pts),                                                        &
  cca_2d(row_length,rows),                                                    &
  cs_ch4_soilt(land_pts,nsoilt),                                              &
  dhf_surf_minus_soil(land_pts),                                              &
  flash_rate_ancil(row_length,rows),                                          &
  pop_den_ancil(row_length,rows),                                             &
  acc_lake_evap(row_length,rows),                                             &
  inlandout_atm_gb(land_pts),                                                 &
  delta_lambda,                                                               &
  delta_phi,                                                                  &
  xx_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                        &
                        tdims_s%j_start:tdims_s%j_end),                       &
  xpa(aocpl_row_length+1),                                                    &
  xua(0:aocpl_row_length),                                                    &
  xva(aocpl_row_length+1),                                                    &
  ypa(aocpl_p_rows),                                                          &
  yua(aocpl_p_rows),                                                          &
  yva(0:aocpl_p_rows),                                                        &
  trivdir(river_row_length, river_rows),                                      &
  trivseq(river_row_length, river_rows),                                      &
  r_area(row_length, rows),                                                   &
  slope(row_length, rows),                                                    &
  flowobs1(row_length, rows),                                                 &
  r_inext(row_length, rows),                                                  &
  r_jnext(row_length, rows),                                                  &
  r_land(row_length, rows),                                                   &
  substore(row_length, rows),                                                 &
  surfstore(row_length, rows),                                                &
  flowin(row_length, rows),                                                   &
  bflowin(row_length, rows),                                                  &
  smvcst(land_pts),                                                           &
  smvcwt(land_pts),                                                           &
  satcon(land_pts),                                                           &
  snow_depth(row_length,rows),                                                &
  surf_ht_flux_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
  surf_ht_flux_ld(land_pts),                                                  &
  tot_sub_runoff(land_pts),                                                   &
  tot_surf_runoff(land_pts),                                                  &
  twatstor(river_row_length, river_rows),                                     &
  ls_graup(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
  ls_rainfrac_land(land_pts)

!------------------------------------------------------------------------------
!End of header

!!! added by Shugong Wang
asteps_since_triffid = timestep 

!Initialise the olr diagnostic
olr(:,:) = 0.0

CALL update_derived_variables()  

!-----------------------------------------------------------------------
!   Conditional allocation for calculation of diagnostics
!-----------------------------------------------------------------------

IF (sf_diag%suv10m_n) THEN
  ALLOCATE(sf_diag%cdr10m_n(pdims_s%i_start:pdims_s%i_end,                    &
                            pdims_s%j_start:pdims_s%j_end))
  ALLOCATE(sf_diag%cdr10m_n_u(udims%i_start:udims%i_end,                      &
                              udims%j_start:udims%j_end))
  ALLOCATE(sf_diag%cdr10m_n_v(vdims%i_start:vdims%i_end,                      &
                              vdims%j_start:vdims%j_end))
  ALLOCATE(sf_diag%cd10m_n(pdims_s%i_start:pdims_s%i_end,                     &
                           pdims_s%j_start:pdims_s%j_end))
ELSE
  ALLOCATE(sf_diag%cdr10m_n(1,1))
  ALLOCATE(sf_diag%cdr10m_n_u(1,1))
  ALLOCATE(sf_diag%cdr10m_n_v(1,1))
  ALLOCATE(sf_diag%cd10m_n(1,1))
END IF

!Calculate the cosine of the zenith angle
IF ( l_cosz ) THEN
  CALL zenith(cosz_ij)
ELSE
  !     Set cosz to a default of 1.0
  cosz_ij(:,:) = 1.0
END IF

IF ( i_sea_alb_method == 3 ) THEN
  ! Calculate the 10m wind speed.
  IF (timestep == 1) THEN
    ! the u/v at 10m not calculated yet, so make a first guess:
    ws10m(:,:) = 2.0
  ELSE
    ws10m(:,:) = SQRT(sf_diag%u10m(:,:)**2 + sf_diag%v10m(:,:)**2)
  END IF
END IF


! put something in ocean nr. surface chlorophyll - a global mean value:
chloro(:,:) = 5.0e-7

CALL surf_couple_radiation(                                                   &
  !Fluxes INTENT(IN)
  tstar_ij,                                                                   &
  !Misc INTENT(IN)
  ws10m, chloro,                                                              &
  n_swbands, n_swbands, spec_band_bb, spec_band_bb,                           &
  !Misc INTENT(OUT)
  sea_ice_albedo,                                                             &
  !Fluxes INTENT(OUT)
  alb_surft, land_albedo_ij,                                                  &
  !UM-only args: INTENT(IN)
  pond_frac_cat_sicat, pond_depth_cat_sicat,                                  &
  !(ancil_info mod)
  nsurft, land_pts, land_index, surft_pts, surft_index,                       &
  row_length, rows, ice_fract_ij, ice_fract_ncat_sicat, frac_surft,           &
  !(p_s_parms mod)
  cosz_ij, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb,                        &
  z0_surft, albsoil_soilt,                                                    &
  !(coastal mod)
  flandg, tstar_sice_sicat,                                                   &
  !(prognostics mod)
  snow_mass_ij, snow_mass_sea_sicat, di_ncat_sicat, lai_pft, canht_pft,       &
  rgrain_surft, snow_surft, soot_ij, tstar_surft, ho2r2_orog_gb,              &
  !UM-only args: INTENT(OUT)
  albobs_sc_ij, open_sea_albedo)

!All the calculations for the radiation downward components have been put into
!a subroutine. For standalone only.
CALL calc_downward_rad(sea_ice_albedo, photosynth_act_rad)

!Calculate buoyancy parameters bt and bq. set ct_ctq_1, cq_cm_u_1,
!cq_cm_v_1, dtl_1, dqw_1, du , dv all to zero (explicit coupling)
!and boundary-layer depth.
bt_1(:,:) = 1.0 / tl_1_ij(:,:)
bq_1(:,:) = c_virtual / (1.0 + c_virtual * qw_1_ij(:,:))

!-----------------------------------------------------------------------
! Set date-related information for generate_anthropogenic_heat in sf_expl
!-----------------------------------------------------------------------

curr_year       = current_time%year
curr_day_number = day_of_year(                                                &
            current_time%year, current_time%month, current_time%day,          &
            l_360, l_leap)
curr_hour       = 0 !These could be correctly calculated by taking the
curr_minute     = 0 !appropriate modulus of current_time%time which is
curr_second     = 0 !the number of seconds into the current day

z0m_scm(:,:)       = 0.0
z0h_scm(:,:)       = 0.0
z1_uv_top          = 0.0
z1_tq_top          = 0.0
ddmfx(:,:)         = 0.0
charnock_w(:,:)    = 0.0

!-----------------------------------------------------------------------
!   Explicit calculations.
!-----------------------------------------------------------------------
CALL surf_couple_explicit(                                                    &
  !Misc INTENT(IN)
  bq_1, bt_1, zh, photosynth_act_rad,                                         &
  curr_year, curr_day_number, curr_hour, curr_minute, curr_second,            &
  !Forcing INTENT(IN)
  qw_1_ij, tl_1_ij, pstar_ij, lw_down_ij,                                     &
  !Fluxes INTENT(IN)
  sw_surft, tstar_ij,                                                         &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(OUT)
  fqw_1_ij,ftl_1_ij, ftl_surft, fqw_surft, fqw_sicat, ftl_sicat, fsmc_pft,    &
  emis_surft,                                                                 &
  !Misc INTENT(OUT)
  radnet_sice, rhokm_1, rhokm_land, rhokm_ssi,                                &
  !Out of explicit and into implicit only INTENT(OUT)
  cdr10m,                                                                     &
  alpha1, alpha1_sea, alpha1_sice, ashtf_prime, ashtf_prime_sea,              &
  ashtf_prime_surft, epot_surft,                                              &
  fraca, resfs, resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea_ij,          &
  dtstar_ij_surft, dtstar_ij_sea, dtstar_ij_sice,                             &
  z0hssi, z0h_surft, z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft,   &
  wt_ext_surft, flake,                                                        &
  !Out of explicit and into extra only INTENT(OUT)
  hcons_soilt,                                                                &
  !Out of explicit and into implicit and extra INTENT(OUT)
  tile_frac,                                                                  &
  !Additional arguments for the UM-----------------------------------------
  !JULES prognostics module
  canopy_surft, snow_surft, k_sice_sicat, cs_pool_gb_um, canht_pft, lai_pft,  &
  t_soil_soilt, tsurf_elev_surft, ti_sicat, ti_cat_sicat,                     &
  tstar_surft, z0msea_ij, smc_soilt, gc_surft, gs_gb,                         &
  !JULES ancil_info module
  land_pts, z1_uv_ij, z1_tq_ij, land_index, nsurft, ice_fract_ncat_sicat,     &
  frac_surft, surft_index, surft_pts,                                         &
  ! IN input data from the wave model
  charnock_w,                                                                 &
  !JULES coastal module
  fland, flandg, tstar_sea_ij, tstar_sice_sicat, vshr_land_ij, vshr_ssi_ij,   &
  !JULES aero module
  co2_3d_ij, rho_aresist_ij, aresist_ij, resist_b_ij, rho_aresist_surft,      &
  aresist_surft, resist_b_surft, r_b_dust_ij, cd_std_dust_ij, u_s_std_surft,  &
  !JULES trifctl module
  asteps_since_triffid, g_leaf_acc_pft, npp_acc_pft, resp_w_acc_pft,          &
  resp_s_acc_gb_um, gpp_gb, npp_gb, resp_p_gb, g_leaf_pft, gpp_pft, npp_pft,  &
  resp_p_pft, resp_s_gb_um, resp_w_pft,                                       &
  !JULES p_s_parms module
  catch_surft, catch_snow_surft, hcon_soilt, smvccl_soilt, smvcst_soilt,      &
  smvcwt_soilt, sthf_soilt, sthu_soilt, z0_surft,                             &
  z0h_bare_surft, z0m_soil_gb, albsoil_soilt, cosz_ij, soil_clay_ij,          &
  !JULES orog module
  ho2r2_orog_gb, sil_orog_land_gb, h_blend_orog_ij, z0m_eff_ij,               &
  !JULES u_v_grid module
  u_1_p_ij, v_1_p_ij, u_0_p_ij, v_0_p_ij,                                     &
  !JULES switches module
  l_spec_z0,                                                                  &
  !JULES c_elevate module
  z_land_ij,                                                                  &
  !Not in a JULES module
  numcycles, cycleno, z1_uv_top, z1_tq_top, ddmfx,                            &
  l_aero_classic, z0m_scm, z0h_scm, recip_l_mo_sea, rib,                      &
  flandfac, fseafac, fb_surf, u_s, t1_sd, q1_sd, rhostar,                     &
  vshr, resp_s_tot_soilt, emis_soil)

!  Items that are only passed from explicit to implicit.
!  Cross-check for exlusive use in the UM too
!  rhokm_1, cdr10m, alpha1, alpha1_sice, ashtf_prime, ashtf_prime_surft,
!  epot_surft, fraca, rhokh, rhokh_surft, rhokh_sice, dtstar_ij_surft,
!  dtstar_ij_sea, dtstar_ij_sice, z0hssi, z0h_surft,
!  z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, wt_ext_surft, flake

!  Items that are passed from explicit to implicit and/or extra.
!  Cross-check in UM as above.
!  hcons_soilt, tile_frac

!Calculate variables required by sf_impl2
!In the UM, these are message passing variables and variables required by
!the implicit solver

!(i.e. no points are part land/part sea)
flandfac_u(:,:)   = 1.0
flandfac_v(:,:)   = 1.0
fseafac_u(:,:)    = 1.0
fseafac_v(:,:)    = 1.0

IF (l_ctile .AND. buddy_sea == on) THEN
  taux_land_ij(:,:) =                                                         &
              rhokm_land(:,:) *  u_1_ij(:,:)                * flandfac_u(:,:)
  taux_ssi_ij(:,:)  =                                                         &
              rhokm_ssi(:,:)  * (u_1_ij(:,:) - u_0_ij(:,:)) * fseafac_u(:,:)
  tauy_land_ij(:,:) =                                                         &
              rhokm_land(:,:) *  v_1_ij(:,:)                * flandfac_v(:,:)
  tauy_ssi_ij(:,:)  =                                                         &
              rhokm_ssi(:,:)  * (v_1_ij(:,:) - v_0_ij(:,:)) * fseafac_v(:,:)
ELSE   ! Standard code
  taux_land_ij(:,:) = rhokm_land(:,:) *  u_1_ij(:,:)
  taux_ssi_ij(:,:)  = rhokm_ssi(:,:)  * (u_1_ij(:,:) - u_0_ij(:,:))
  tauy_land_ij(:,:) = rhokm_land(:,:) *  v_1_ij(:,:)
  tauy_ssi_ij(:,:)  = rhokm_ssi(:,:)  * (v_1_ij(:,:) - v_0_ij(:,:))
END IF

taux_1_ij(:,:) = flandg(:,:) * taux_land_ij(:,:) +                            &
                 (1.0 - flandg(:,:)) * taux_ssi_ij(:,:)
tauy_1_ij(:,:) = flandg(:,:) * tauy_land_ij(:,:) +                            &
                 (1.0 - flandg(:,:)) * tauy_ssi_ij(:,:)

IF (sf_diag%suv10m_n) THEN
  sf_diag%cdr10m_n_u(:,:) = sf_diag%cdr10m_n(:,:)
  sf_diag%cdr10m_n_v(:,:) = sf_diag%cdr10m_n(:,:)
END IF

!-----------------------------------------------------------------------
! Implicit calculations.
!
! In the new boundary layer implicit solver in the UM, sf_impl2 is
! called twice - once with l_correct = .FALSE. and once with
! l_correct = .TRUE.
! The variables r_gamma1 and r_gamma2 that determine weights for this new
! solver are set so that the new scheme is the same as the old scheme
! (i.e. fully explicit coupling)
!-----------------------------------------------------------------------

! Set values of inputs for the implicit solver so that we get
! explicit coupling
du(:,:)   = 0.0
dv(:,:)   = 0.0
r_gamma1(:,:) = r_gamma

! To get explicit coupling with the new scheme, r_gamma2 is 0 everywhere
! for the 1st call, and 1 everywhere for the 2nd call
r_gamma2(:,:) = 0.0
! Call sf_impl2 with
l_correct = .FALSE.

ctctq1(:,:)     = 0.0
dqw1_1(:,:)     = 0.0
dtl1_1(:,:)     = 0.0
du_star1(:,:)   = 0.0
dv_star1(:,:)   = 0.0
cq_cm_u_1(:,:)  = 0.0
cq_cm_v_1(:,:)  = 0.0
rho1(:,:)       = 0.0
f3_at_p(:,:)    = 0.0
uStargbm(:,:)   = 0.0
flandg_u(:,:)   = flandg(:,:)
flandg_v(:,:)   = flandg(:,:)

CALL surf_couple_implicit(                                                    &
  !Important switch
  l_correct,                                                                  &
  !Forcing INTENT(IN)
  pstar_ij, lw_down_ij, qw_1_ij, tl_1_ij, u_1_ij, v_1_ij, u_0_ij, v_0_ij,     &
  !Fluxes INTENT(IN)
  sw_surft, emis_surft,                                                       &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  rhokm_1, rhokm_1, r_gamma1, r_gamma2, alpha1, alpha1_sea, alpha1_sice,      &
  ashtf_prime, ashtf_prime_sea, ashtf_prime_surft, du, dv, fraca, resfs,      &
  resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea_ij, z0hssi, z0mssi,        &
  z0h_surft, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, flake, tile_frac, &
  wt_ext_surft, cdr10m, cdr10m, r_gamma,                                      &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(INOUT)
  fqw_sicat, ftl_sicat, fqw_surft, fqw_1_ij, ftl_1_ij, ftl_surft,             &
  !Misc INTENT(INOUT)
  epot_surft, dtstar_ij_surft, dtstar_ij_sea, dtstar_ij_sice, radnet_sice,    &
  olr,                                                                        &
  !Fluxes INTENT(OUT)
  tstar_ij, le_surft, radnet_surft, e_sea_ij, h_sea_ij, taux_1_ij, tauy_1_ij, &
  ecan_surft, ei_ij,                                                          &
  esoil_ij_soilt, ext_soilt, snowmelt_ij, melt_surft,                         &
  ecan_ij, ei_surft, esoil_surft, sea_ice_htf_sicat, surf_ht_flux_ij,         &
  surf_htf_surft,                                                             &
  !Misc INTENT(OUT)
  error,                                                                      &
  !UM-only arguments
  !JULES ancil_info module
  nsurft, land_pts, land_index, surft_index, surft_pts, ice_fract_ij,         &
  sstfrz_ij, ice_fract_ncat_sicat, z1_tq_ij,                                  &
  !JULES prognostics module
  canopy_surft, smc_soilt, k_sice_sicat, t_soil_soilt, ti_sicat, snow_surft,  &
  di_ncat_sicat, tstar_surft,                                                 &
  !JULES coastal module
  fland, flandg, tstar_sea_ij, tstar_sice_sicat, tstar_ssi_ij,                &
  taux_land_ij, tauy_land_ij, taux_ssi_ij, tauy_ssi_ij,                       &
  surf_ht_flux_land_ij, surf_ht_flux_sice_sicat, tstar_land_ij, tstar_sice_ij,&
  !JULES u_v_grid module
  dtrdz_charney_grid_1_ij,                                                    &
  !JULES switches module
  !==l_mr_physics in UM
  l_co2_interactive, lq_mix_bl,                                               &
  !JULES aero module
  co2_3d_ij,                                                                  &
  !Arguments without a JULES module
  ctctq1,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,flandg_u,flandg_v, &
  rho1, f3_at_p, uStarGBM,tscrndcl_ssi,tscrndcl_surft,tStbTrans,              &
  taux_land_star,tauy_land_star,taux_ssi_star,                                &
  tauy_ssi_star,ei_sice,rhokh_mix, ti_gb)


! Adjust the value of r_gamma2 to ensure explicit coupling
r_gamma2(:,:) = 1.0
! Call sf_impl2 again with
l_correct = .TRUE.

ctctq1(:,:)     = 0.0
dqw1_1(:,:)     = 0.0
dtl1_1(:,:)     = 0.0
du_star1(:,:)   = 0.0
dv_star1(:,:)   = 0.0
cq_cm_u_1(:,:)  = 0.0
cq_cm_v_1(:,:)  = 0.0
rho1(:,:)       = 0.0
f3_at_p(:,:)    = 0.0
uStargbm(:,:)   = 0.0
flandg_u(:,:)   = flandg(:,:)
flandg_v(:,:)   = flandg(:,:)


CALL surf_couple_implicit(                                                    &
  !Important switch
  l_correct,                                                                  &
  !Forcing INTENT(IN)
  pstar_ij, lw_down_ij, qw_1_ij, tl_1_ij, u_1_ij, v_1_ij, u_0_ij, v_0_ij,     &
  !Fluxes INTENT(IN)
  sw_surft, emis_surft,                                                       &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  rhokm_1, rhokm_1, r_gamma1, r_gamma2, alpha1, alpha1_sea, alpha1_sice,      &
  ashtf_prime, ashtf_prime_sea, ashtf_prime_surft, du, dv, fraca, resfs,      &
  resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea_ij, z0hssi, z0mssi,        &
  z0h_surft, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, flake, tile_frac, &
  wt_ext_surft, cdr10m, cdr10m, r_gamma,                                      &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(INOUT)
  fqw_sicat, ftl_sicat, fqw_surft, fqw_1_ij, ftl_1_ij, ftl_surft,             &
  !Misc INTENT(INOUT)
  epot_surft, dtstar_ij_surft, dtstar_ij_sea, dtstar_ij_sice, radnet_sice,    &
  olr,                                                                        &
  !Fluxes INTENT(OUT)
  tstar_ij, le_surft, radnet_surft, e_sea_ij, h_sea_ij, taux_1_ij, tauy_1_ij, &
  ecan_surft, ei_ij,                                                          &
  esoil_ij_soilt, ext_soilt, snowmelt_ij, melt_surft,                         &
  ecan_ij, ei_surft, esoil_surft, sea_ice_htf_sicat, surf_ht_flux_ij,         &
  surf_htf_surft,                                                             &
  !Misc INTENT(OUT)
  error,                                                                      &
  !UM-only arguments
  !JULES ancil_info module
  nsurft, land_pts, land_index, surft_index, surft_pts, ice_fract_ij,         &
  sstfrz_ij, ice_fract_ncat_sicat, z1_tq_ij,                                  &
  !JULES prognostics module
  canopy_surft, smc_soilt, k_sice_sicat, t_soil_soilt, ti_sicat, snow_surft,  &
  di_ncat_sicat, tstar_surft,                                                 &
  !JULES coastal module
  fland, flandg, tstar_sea_ij, tstar_sice_sicat, tstar_ssi_ij,                &
  taux_land_ij, tauy_land_ij, taux_ssi_ij, tauy_ssi_ij,                       &
  surf_ht_flux_land_ij, surf_ht_flux_sice_sicat, tstar_land_ij, tstar_sice_ij,&
  !JULES u_v_grid module
  dtrdz_charney_grid_1_ij,                                                    &
  !JULES switches module
  !==l_mr_physics in UM
  l_co2_interactive, lq_mix_bl,                                               &
  !JULES aero module
  co2_3d_ij,                                                                  &
  !Arguments without a JULES module
  ctctq1,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,flandg_u,flandg_v, &
  rho1, f3_at_p, uStarGBM,tscrndcl_ssi,tscrndcl_surft,tStbTrans,              &
  taux_land_star,tauy_land_star,taux_ssi_star,                                &
  tauy_ssi_star,ei_sice,rhokh_mix, ti_gb)

!-------------------------------------------------------------------------------
! Calculation of pseudostresses.
!-------------------------------------------------------------------------------
IF (sf_diag%suv10m_n) THEN
  sf_diag%mu10m_n(:,:) = (taux_1_ij(:,:)) / sf_diag%cd10m_n(:,:)
  sf_diag%mv10m_n(:,:) = (tauy_1_ij(:,:)) / sf_diag%cd10m_n(:,:)
END IF
DEALLOCATE(sf_diag%cdr10m_n)
DEALLOCATE(sf_diag%cdr10m_n_u)
DEALLOCATE(sf_diag%cdr10m_n_v)
DEALLOCATE(sf_diag%cd10m_n)


!Arguments to surf_couple_extra that are not actively used in standalone JULES
!However, we expect a sensible value to be stored in them
ls_graup(:,:) = 0.0
cca_2d(:,:) = 0.0
ls_rainfrac_land(:) = 0.0

  !Note that the commented intents may not be correct. The declarations in 
  !surf_couple_extra will be more correct

CALL surf_couple_extra(                                                       &
  !Driving data and associated INTENT(IN)
  ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij, tl_1_ij, lw_down_ij,      &
  qw_1_ij, u_1_ij, v_1_ij, pstar_ij,                                          &
  !Fluxes INTENT(IN)
  ei_surft, surf_htf_surft, ecan_surft, ext_soilt, sw_surft,                  &
  !Misc INTENT(IN)
  timestep, sf_diag%smlt, tile_frac, hcons_soilt,                             &
  !Fluxes INTENT(INOUT)
  melt_surft,                                                                 &
  !Fluxes INTENT(OUT)
  hf_snow_melt_gb, snowmelt_ij, snomlt_sub_htf_gb, sub_surf_roff_gb,          &
  surf_roff_gb, tot_tfall_gb, snow_melt_gb, rrun_gb, rflow_gb,                &
  snow_soil_htf,                                                              &
  !Arguments for the UM-----------------------------------------
  !IN
  land_pts, row_length, rows, river_row_length, river_rows, land_index,       &
  ls_graup,                                                                   &
  cca_2d, nsurft, surft_pts, surft_index,                                     &
  tstar_surft,                                                                &
  lice_pts, lice_index, soil_pts, soil_index,                                 &
  stf_sub_surf_roff,                                                          &
  fexp_soilt, gamtot_soilt, ti_mean_soilt, ti_sig_soilt,                      &
  cs_ch4_soilt, flash_rate_ancil, pop_den_ancil,                              &
  a_fsat_soilt, c_fsat_soilt, a_fwet_soilt, c_fwet_soilt,                     &
  ntype, fqw_surft,                                                           &
  halo_i, halo_j, model_levels,                                               &
  delta_lambda, delta_phi, xx_cos_theta_latitude,                             &
  aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa, yua, yva,               &
  g_p_field, g_r_field, n_proc, global_row_length, global_rows,               &
  global_river_row_length, global_river_rows, flandg,                         &
  trivdir, trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,        &
  n_rows, offx, offy, n_procx, n_procy, g_rows, g_row_length,                 &
  at_extremity, frac_agr_gb, soil_clay_ij, resp_s_soilt, npp_gb,              &
  z0m_soil_gb,                                                                &
  !INOUT
  a_steps_since_riv,  t_soil_soilt, tsurf_elev_surft,                         &
  rgrain_surft, snow_grnd_surft, snow_surft,                                  &
  smcl_soilt, sthf_soilt, sthu_soilt, canopy_surft, fsat_soilt, fwetl_soilt,  &
  zw_soilt, sthzw_soilt,                                                      &
  snow_mass_ij,  ls_rainfrac_land,                                            &
  substore, surfstore, flowin, bflowin,                                       &
  tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,                   &
  asteps_since_triffid, g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,     &
  resp_s_acc_gb_um,                                                           &
  resp_w_acc_pft, cs_pool_gb_um, frac_surft, lai_pft, canht_pft,              &
  catch_snow_surft, catch_surft, infil_surft,                                 &
  inlandout_atm_gb,                                                           &
  !OUT
  dhf_surf_minus_soil,                                                        &
  canopy_gb,  smc_soilt,                                                      &
  z0_surft, z0h_bare_surft,                                                   &
  land_sea_mask                                                               &
  )

END SUBROUTINE lis_control
!#endif
