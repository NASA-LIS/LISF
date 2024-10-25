!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine prog_to_tile(n,t, pft)
  use jules5x_lsmMod
  use prognostics
  use jules_snow_mod, only : nsmax
  use jules_vegetation_mod, only:  l_triffid, l_phenol, irr_crop, l_nitrogen, l_irrig_dmd
  use jules_soil_mod, only: l_bedrock
  use jules_surface_types_mod,  only: npft, nnvg, ntype
  implicit none 
  integer, intent(in) :: n, t, pft


  jules5x_struc(n)%jules5x(t)%nsnow(:)              = nsnow_surft(1, :)                  ! ntiles  number of snow layers on ground on tiles
  if( l_bedrock ) then
    jules5x_struc(n)%jules5x(t)%tsoil_deep(:)         = tsoil_deep_gb(1, :)             ! ns_deep  deep soil temperatures (k)
  endif 
  if (nsmax > 0)then
    jules5x_struc(n)%jules5x(t)%sice(pft,:)             = sice_surft(1, pft, :)                 ! temp_tiles,temp_layers, Snow layer ice mass on tiles (Kg/m2)
    jules5x_struc(n)%jules5x(t)%sliq(pft,:)             = sliq_surft(1, pft, :)                 ! temp_tiles,temp_layers, Snow layer liquid mass on tiles (Kg/m2)
    jules5x_struc(n)%jules5x(t)%tsnow(pft,:)            = tsnow_surft(1, pft,:)                ! temp_tiles,temp_layers,Snow layer temperature (K)
    jules5x_struc(n)%jules5x(t)%rgrainl(pft,:)          = rgrainl_surft(1, pft,:)              ! temp_tiles,temp_layers, Snow layer grain size on tiles (microns)
    jules5x_struc(n)%jules5x(t)%ds(pft,:)               = ds_surft(1, pft,:)                   ! temp_tiles,temp_layers, Snow layer thickness (m)
  else
    jules5x_struc(n)%jules5x(t)%sice(:,:)             = sice_surft(1, :, :)                 ! temp_tiles,temp_layers, Snow layer ice mass on tiles (Kg/m2)
    jules5x_struc(n)%jules5x(t)%sliq(:,:)             = sliq_surft(1, :, :)                 ! temp_tiles,temp_layers, Snow layer liquid mass on tiles (Kg/m2)
    jules5x_struc(n)%jules5x(t)%tsnow(:,:)            = tsnow_surft(1, :,:)                ! temp_tiles,temp_layers,Snow layer temperature (K)
    jules5x_struc(n)%jules5x(t)%rgrainl(:,:)          = rgrainl_surft(1, :,:)              ! temp_tiles,temp_layers, Snow layer grain size on tiles (microns)
    jules5x_struc(n)%jules5x(t)%ds(:,:)               = ds_surft(1, :,:)                   ! temp_tiles,temp_layers, Snow layer thickness (m)
  endif

  jules5x_struc(n)%jules5x(t)%rho_snow_grnd(:)      = rho_snow_grnd_surft(1, :)          ! ntiles, Snowpack bulk density (kg/m3)
  if(nsmax>0) then
    jules5x_struc(n)%jules5x(t)%rho_snow(pft,:)         = rho_snow_surft(1, pft,:)             ! ntiles, nsmax, Snow layer densities (m)
  endif 
  
  if ( l_triffid .or. l_phenol ) then
    jules5x_struc(n)%jules5x(t)%wood_prod_fast        = wood_prod_fast_gb(1)           ! Fast-turnover wood product C pool.
    jules5x_struc(n)%jules5x(t)%wood_prod_med         = wood_prod_med_gb(1)             ! Medium-turnover wood product C pool. 
    jules5x_struc(n)%jules5x(t)%wood_prod_slow        = wood_prod_slow_gb(1)            ! Slow-turnover wood product C pool.
    jules5x_struc(n)%jules5x(t)%frac_agr_prev         = frac_agr_prev_gb(1)             ! Agricultural fraction from previous TRIFFID call
    jules5x_struc(n)%jules5x(t)%frac_past_prev        = frac_past_prev_gb(1)            ! 
  endif 
  
  jules5x_struc(n)%jules5x(t)%n_inorg_soilt_lyrs(:)   = n_inorg_soilt_lyrs(1,1,:)   
  jules5x_struc(n)%jules5x(t)%n_inorg                 = n_inorg_gb(1)                   ! Gridbox Inorganic N pool (kg N m-2)
  jules5x_struc(n)%jules5x(t)%n_inorg_avail_pft(:,:)  = n_inorg_avail_pft(1,:,:) 
  jules5x_struc(n)%jules5x(t)%ns(:,:)                 = ns_pool_gb(1, :,:)                     ! dim_cs1, Soil Organic Nitrogen (kg N m-2)
  
  if(l_triffid) then 
    jules5x_struc(n)%jules5x(t)%triffid_co2_gb          = triffid_co2_gb(1) 
  endif 

  jules5x_struc(n)%jules5x(t)%canht_ft(1:npft)           = canht_pft(1, 1:npft)               ! npft, Canopy height (m)
  jules5x_struc(n)%jules5x(t)%lai(1:npft)                = lai_pft(1, 1:npft)                    ! npft LAI of plant functional types
  jules5x_struc(n)%jules5x(t)%canopy(pft)             = canopy_surft(1, pft)                 ! ntiles, Surface/canopy water for snow-free land tiles (kg m-2)
  jules5x_struc(n)%jules5x(t)%canopy_gb             = canopy_gb(1)                 ! Gridbox canopy water content (kg m-2)

  jules5x_struc(n)%jules5x(t)%cs_pool_soilt(:,:)   = cs_pool_soilt(1, 1,:,:)                     ! dim_cs1, Soil carbon (kg C/m2)
  jules5x_struc(n)%jules5x(t)%di_ncat(:)            = di_ncat_sicat(1,1, :)                ! nice, "Equivalent thickness" of sea-ice catagories (m)
  jules5x_struc(n)%jules5x(t)%k_sice(:)             = k_sice_sicat(1,1, :)                 ! nice, Sea ice effective conductivity (2*kappai/de)
  
  jules5x_struc(n)%jules5x(t)%gc_surft(pft)                 = gc_surft(1, pft)                     ! ntiles, Stomatal" conductance to evaporation for land tiles(m s-1)
  jules5x_struc(n)%jules5x(t)%gs_gb                    = gs_gb(1)                        !  "Stomatal" conductance to evaporation (m s-1)
  
  jules5x_struc(n)%jules5x(t)%rgrain(pft)             = rgrain_surft(1, pft)                 ! ntiles, Snow surface grain size on tiles (microns)
  jules5x_struc(n)%jules5x(t)%smc_soilt                   = smc_soilt(1,1)                       !  Soil moisture in a layer at the surface (kg m-2).
  jules5x_struc(n)%jules5x(t)%smcl_soilt(:)               = smcl_soilt(1,1, :)                   ! sm_levels, Soil moisture content of layers (kg m-2)
  jules5x_struc(n)%jules5x(t)%snowdepth(pft)          = snowdepth_surft(1, pft)              ! ntiles, Snow depth on ground on tiles (m)
  jules5x_struc(n)%jules5x(t)%snow_tile(pft)          = snow_surft(1, pft)              ! ntiles, Lying snow on tiles (kg m-2)
  jules5x_struc(n)%jules5x(t)%snow_grnd(pft)          = snow_grnd_surft(1, pft)              ! ntiles, Snow on the ground (kg m-2)
  jules5x_struc(n)%jules5x(t)%snow_mass_ij             = snow_mass_ij(1,1)                 !  Gridbox snowmass (kg m-2)
  jules5x_struc(n)%jules5x(t)%snow_mass_sea_sicat(:) = snow_mass_sea_sicat(1, 1, :)     !nice_use,  Snow on category sea-ice (Kg/m2)
  jules5x_struc(n)%jules5x(t)%soot_ij                  = soot_ij(1,1)                      !  Snow soot content (kg kg-1)
  jules5x_struc(n)%jules5x(t)%t_soil(:)             = t_soil_soilt(1,1, :)                 ! sm_levels, Sub-surface temperatures (K)
  jules5x_struc(n)%jules5x(t)%ti_sicat(:)                    = ti_sicat(1,1,:)                        !  Sea-ice surface layer
  jules5x_struc(n)%jules5x(t)%tstar_tile(:)         = tstar_surft(1, :)             ! ntiles, Tile surface temperatures (K)
  jules5x_struc(n)%jules5x(t)%tsurf_elev_surft(pft)   = tsurf_elev_surft(1, pft)
  jules5x_struc(n)%jules5x(t)%z0msea                  = z0msea_ij(1,1)                    !  Sea-surface roughness length for momentum (m).
end subroutine prog_to_tile
