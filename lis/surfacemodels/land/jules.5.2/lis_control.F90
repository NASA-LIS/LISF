!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#if defined(LIS_JULES)
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

!Imports for setting array dimensions
USE ancil_info,               ONLY:                                         &
  land_pts, nsurft, row_length, rows, nsoilt

!Imports for driving and flux variables
USE forcing,                  ONLY:                                         &
  con_rain_ij, con_snow_ij, ls_rain_ij, ls_snow_ij, lw_down_ij, pstar_ij,   &
  qw_1_ij, sw_down_ij, tl_1_ij, u_0_ij, u_1_ij, v_0_ij, v_1_ij, diff_rad_ij

USE fluxes,                   ONLY:                                         &
  alb_surft, e_sea_ij, ecan_ij, ecan_surft, ei_ij, ei_surft, esoil_ij_soilt,&
  esoil_surft, ext_soilt, fqw_1_ij, fqw_surft, fqw_sicat, fsmc_pft, ftl_1_ij,&
  ftl_sicat, ftl_surft, h_sea_ij, hf_snow_melt_gb, land_albedo_ij,          &
  le_surft, melt_surft, sea_ice_htf_sicat, snomlt_sub_htf_gb, snow_melt_gb, &
  snowmelt_ij, sub_surf_roff_gb, surf_ht_flux_ij, surf_htf_surft,           &
  surf_roff_gb, radnet_surft, taux_1_ij, tauy_1_ij, tot_tfall_gb, tstar_ij, &
  sw_surft, emis_surft, alb_sicat, rflow_gb, rrun_gb,                       &
  snow_soil_htf, z0m_surft, z0h_surft

!Imports for the remaining calculations in lis_control
USE bl_option_mod,            ONLY:                                         &
  on
USE sf_diags_mod, ONLY: sf_diag
USE planet_constants_mod,     ONLY:                                         &
  c_virtual
USE coastal,                  ONLY:                                         &
  flandg,                                                                   &
  taux_land_ij, taux_ssi_ij, tauy_land_ij, tauy_ssi_ij
USE p_s_parms,                ONLY:                                         &
  cosz_ij
USE jules_radiation_mod,      ONLY:                                         &
  i_sea_alb_method, l_cosz
USE jules_sea_seaice_mod,     ONLY:                                         &
  l_ctile, buddy_sea, nice_use
USE jules_soil_mod,           ONLY:                                         &
  sm_levels
USE zenith_mod,               ONLY:                                         &
  zenith
USE atm_fields_bounds_mod,    ONLY:                                         &
  tdims, pdims, udims, vdims, pdims_s, udims_s, vdims_s
USE theta_field_sizes,        ONLY:                                         &
  t_i_length, t_j_length
USE datetime_mod,             ONLY:                                         &
  secs_in_day, l_360, l_leap
USE datetime_utils_mod,       ONLY:                                         &
  day_of_year
USE model_time_mod,           ONLY:                                         &
  current_time
USE trifctl,                  ONLY: asteps_since_triffid

!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
INTEGER :: a_step     ! Atmospheric timestep number.

! Local variables
REAL ::                                                                     &
  olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !TOA-surface upward LW on last radiation timestep Corrected TOA outward LW
  sea_ice_albedo(row_length,rows,4),                                        &
    ! Sea ice albedo
  ws10m(t_i_length,t_j_length),                                             &
    ! 10m wind speed (m s-1)
  photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),  &
    ! Net downward shortwave radiation in band 1 (w/m2).
  bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    ! A buoyancy parameter (beta T tilde).
  bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    ! A buoyancy parameter (beta q tilde).
  radnet_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    ! Surface net radiation on open sea (W/m2)
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),&
    ! Surface net radiation on sea-ice (W/m2)
  rhokm_land(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),   &
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  alpha1(land_pts,nsurft),                                                  &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    ! ALPHA1 for open sea.
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),&
    ! ALPHA1 for sea-ice.
  ashtf_prime(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),&
    ! Adjusted SEB coefficient for sea-ice
  ashtf_prime_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
    ! Adjusted SEB coefficient for open sea
  ashtf_prime_surft(land_pts,nsurft),                                       &
    ! Adjusted SEB coefficient for land tiles
  zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    ! Height above surface of top of boundary layer (metres).
  rhokm_1(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),     &
    ! Exchange coefficients for momentum on P-grid
  epot_surft(land_pts,nsurft),                                              &
    ! Local EPOT for land tiles.
  fraca(land_pts,nsurft),                                                   &
    !Fraction of surface moisture flux with only aerodynamic resistance for
    !snow-free land tiles.
  resfs(land_pts,nsurft),                                                   &
    ! Combined soil, stomatal and aerodynamic resistance factor for fraction
    !(1-FRACA) of snow-free land tiles.
  resft(land_pts,nsurft),                                                   &
    !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1
    !for snow.
  rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    ! Grid-box surface exchange coefficients
  rhokh_surft(land_pts,nsurft),                                             &
    ! Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use), &
    ! Surface exchange coefficients for sea-ice
  rhokh_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
    ! Surface exchange coefficients for sea
  dtstar_ij_surft(land_pts,nsurft),                                         &
    ! Change in tstar_ij over timestep for land tiles
  dtstar_ij_sea(tdims%i_start:tdims%i_end,                                  &
                tdims%j_start:tdims%j_end),                                 &
    ! Change is tstar_ij over timestep for open sea
  dtstar_ij_sice(tdims%i_start:tdims%i_end,                                 &
                tdims%j_start:tdims%j_end,nice_use),                        &
    ! Change is tstar_ij over timestep for sea-ice
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
    ! Roughness length for heat and moisture over sea (m).
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
    ! Roughness length for momentum over sea (m).
  chr1p5m(land_pts,nsurft),                                                 &
    ! Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
    ! CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                             &
    ! Areal heat capacity of canopy for land tiles (J/K/m2).
  wt_ext_surft(land_pts,sm_levels,nsurft),                                  &
    !Fraction of evapotranspiration which is extracted from each soil layer
    !by each tile.
  flake(land_pts,nsurft),                                                   &
    !Lake fraction.
  tile_frac(land_pts,nsurft),                                               &
    !Tile fractions including snow cover in the ice tile.
  hcons_soilt(land_pts,nsoilt),                                             &
    !Thermal conductivity of top soil layer, including water and ice (W/m/K)
  flandfac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),          &
  flandfac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),          &
  fseafac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),           &
  fseafac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),           &
  du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),          &
    !Level 1 increment to u wind field
  dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),          &
    !Level 1 increment to v wind field
  r_gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
    !weights for new BL solver
  r_gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
  chloro(row_length,rows)
    !nr surface chlorophyll content

!-------------------------------------------------------------------------------
! Variables required for 10 m equivalent neutral winds.
!-------------------------------------------------------------------------------
REAL, ALLOCATABLE :: cdr10m_n(:,:)
    ! Interpolation coefficients for equivalent neutral 10m winds
REAL, ALLOCATABLE :: cd10m_n(:,:)
    ! Neutral drag coefficient for calculation of pseudostress

REAL, PARAMETER :: r_gamma = 2.0
    ! Implicit weighting coefficient

INTEGER :: error, &     ! OUT 0 - AOK; 1 to 7  - bad grid definition detected
           curr_year, curr_day_number, curr_hour, curr_minute, curr_second

LOGICAL :: l_correct

INTEGER, PARAMETER :: n_swbands = 1  !  Nnumber of SW bands
REAL, PARAMETER :: spec_band_bb(n_swbands) = -1.0
  !spectral band boundary (bb=-1)

! reset time steps 
a_step = timestep
asteps_since_triffid = timestep 


!------------------------------------------------------------------------------
!End of header

!Initialise the olr diagnostic
olr(:,:) = 0.0

CALL update_derived_variables()

!-----------------------------------------------------------------------
!   Conditional allocation for calculation of diagnostics
!-----------------------------------------------------------------------

IF (sf_diag%suv10m_n) THEN
  ALLOCATE(cdr10m_n(pdims_s%i_start:pdims_s%i_end,                            &
                    pdims_s%j_start:pdims_s%j_end))
  ALLOCATE(cd10m_n(pdims_s%i_start:pdims_s%i_end,                             &
                   pdims_s%j_start:pdims_s%j_end))
ELSE
  ALLOCATE(cdr10m_n(1,1))
  ALLOCATE(cd10m_n(1,1))
END IF

!Calculate the cosine of the zenith angle
IF ( l_cosz ) THEN
  CALL zenith(cosz_ij)
ELSE
  !     Set cosz to a default of 1.0
  cosz_ij(:,:) = 1.0
END IF
!!! Note
!!! Since we already has 10 m wind, we should use the 10 m wind directly
!!! u_1_ij and v_1_ij are 10m wind from AGRMET and NLDAS forcing. If the height
!!! of wind is not at 10m, the wind speed should be scaled to 10m height. 
!!! Shugong Wang 
ws10m(:,:) = SQRT(u_1_ij(:,:)**2 + v_1_ij(:,:)**2)
!IF ( i_sea_alb_method == 3 ) THEN
!  ! Calculate the 10m wind speed.
!  IF (a_step == 1) THEN
!    ! the u/v at 10m not calculated yet, so make a first guess:
!    ws10m(:,:) = 2.0
!  ELSE
!    ws10m(:,:) = SQRT(sf_diag%u10m(:,:)**2 + sf_diag%v10m(:,:)**2)
!  END IF
!END IF

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
  alb_surft, land_albedo_ij                                                   &
)

!All the calculations for the radiation downward components have been put into
!a subroutine. For standalone only.
CALL calc_downward_rad(sea_ice_albedo, photosynth_act_rad)

!Calculate buoyancy parameters bt and bq. set ct_ctq_1, cq_cm_u_1,
!cq_cm_v_1, dtl_1, dqw_1, du , dv all to zero (explicit coupling)
!and boundary-layer depth.
bt_1(:,:) = 1.0 / tl_1_ij(:,:)
bq_1(:,:) = c_virtual / (1.0 + c_virtual * qw_1_ij(:,:))
zh(:,:)   = 1000.0

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
  radnet_sea, radnet_sice, rhokm_1, rhokm_land, rhokm_ssi,                    &
  !Out of explicit and into implicit only INTENT(OUT)
  cdr10m, cdr10m_n, cd10m_n,                                                  &
  alpha1, alpha1_sea, alpha1_sice, ashtf_prime, ashtf_prime_sea,              &
  ashtf_prime_surft, epot_surft,                                              &
  fraca, resfs, resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea_ij,          &
  dtstar_ij_surft, dtstar_ij_sea, dtstar_ij_sice,                             &
  z0hssi, z0h_surft, z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft,   &
  wt_ext_surft, flake,                                                        &
  !Out of explicit and into extra only INTENT(OUT)
  hcons_soilt,                                                                &
  !Out of explicit and into implicit and extra INTENT(OUT)
  tile_frac                                                                   &
)

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
  wt_ext_surft, cdr10m, cdr10m, cdr10m_n, cdr10m_n, r_gamma,                  &
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
  error                                                                       &
)


! Adjust the value of r_gamma2 to ensure explicit coupling
r_gamma2(:,:) = 1.0
! Call sf_impl2 again with
l_correct = .TRUE.

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
  wt_ext_surft, cdr10m, cdr10m, cdr10m_n, cdr10m_n, r_gamma,                  &
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
  error                                                                       &
)

!-------------------------------------------------------------------------------
! Calculation of pseudostresses.
!-------------------------------------------------------------------------------
IF (sf_diag%suv10m_n) THEN
  sf_diag%mu10m_n(:,:) = (taux_1_ij(:,:)) / cd10m_n(:,:)
  sf_diag%mv10m_n(:,:) = (tauy_1_ij(:,:)) / cd10m_n(:,:)
END IF
IF ( ALLOCATED(cdr10m_n) ) DEALLOCATE(cdr10m_n)
IF ( ALLOCATED(cd10m_n) )  DEALLOCATE(cd10m_n)

CALL surf_couple_extra(                                                       &
  !Driving data and associated INTENT(IN)
  ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij, tl_1_ij, lw_down_ij,      &
  qw_1_ij, u_1_ij, v_1_ij, pstar_ij,                                          &
  !Fluxes INTENT(IN)
  ei_surft, surf_htf_surft, ecan_surft, ext_soilt, sw_surft,                  &
  !Misc INTENT(IN)
  a_step, sf_diag%smlt, tile_frac, hcons_soilt,                               &
  !Fluxes INTENT(INOUT)
  melt_surft,                                                                 &
  !Fluxes INTENT(OUT)
  hf_snow_melt_gb, snowmelt_ij, snomlt_sub_htf_gb, sub_surf_roff_gb,          &
  surf_roff_gb, tot_tfall_gb, snow_melt_gb, rrun_gb, rflow_gb,                &
  snow_soil_htf                                                               &
)

END SUBROUTINE lis_control
#endif
