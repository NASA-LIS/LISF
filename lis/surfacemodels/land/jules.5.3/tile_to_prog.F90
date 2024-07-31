!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine tile_to_prog(n, t, pft)
  use jules53_lsmMod
  use prognostics
  use jules_snow_mod, only : nsmax
  use jules_vegetation_mod, only:  l_triffid, l_phenol, irr_crop, l_nitrogen, l_irrig_dmd
  use jules_soil_mod, only: l_bedrock
  use jules_surface_types_mod,  only: npft, nnvg, ntype
  implicit none 
  integer, intent(in) :: n, t, pft ! pft here is surface type 1-9, not necessary to be plant 


   nsnow_surft(1, pft)             = jules53_struc(n)%jules53(t)%nsnow(pft)                  ! ntiles  number of snow layers on ground on tiles
   if( l_bedrock ) then
    tsoil_deep_gb(1, :)        = jules53_struc(n)%jules53(t)%tsoil_deep(:)             ! ns_deep  deep soil temperatures (k)
   endif
   if (nsmax > 0)then
     sice_surft(1, pft, :)           = jules53_struc(n)%jules53(t)%sice(pft,:)                  ! temp_tiles,temp_layers, Snow layer ice mass on tiles (Kg/m2)
     sliq_surft(1, pft, :)           = jules53_struc(n)%jules53(t)%sliq(pft,:)                  ! temp_tiles,temp_layers, Snow layer liquid mass on tiles (Kg/m2)
     tsnow_surft(1, pft,:)           = jules53_struc(n)%jules53(t)%tsnow(pft,:)                ! temp_tiles,temp_layers,Snow layer temperature (K)
     rgrainl_surft(1, pft,:)         = jules53_struc(n)%jules53(t)%rgrainl(pft,:)              ! temp_tiles,temp_layers, Snow layer grain size on tiles (microns)
     ds_surft(1, pft,:)              = jules53_struc(n)%jules53(t)%ds(pft,:)                   ! temp_tiles,temp_layers, Snow layer thickness (m)
   else
     sice_surft(1, :, :)           = jules53_struc(n)%jules53(t)%sice(:,:)                  ! temp_tiles,temp_layers, Snow layer ice mass on tiles (Kg/m2)
     sliq_surft(1, :, :)           = jules53_struc(n)%jules53(t)%sliq(:,:)                  ! temp_tiles,temp_layers, Snow layer liquid mass on tiles (Kg/m2)
     tsnow_surft(1, :,:)           = jules53_struc(n)%jules53(t)%tsnow(:,:)                ! temp_tiles,temp_layers,Snow layer temperature (K)
     rgrainl_surft(1, :,:)         = jules53_struc(n)%jules53(t)%rgrainl(:,:)              ! temp_tiles,temp_layers, Snow layer grain size on tiles (microns)
     ds_surft(1, :,:)              = jules53_struc(n)%jules53(t)%ds(:,:)                   ! temp_tiles,temp_layers, Snow layer thickness (m)
   endif
   
   rho_snow_grnd_surft(1, :)     = jules53_struc(n)%jules53(t)%rho_snow_grnd(:)          ! ntiles, Snowpack bulk density (kg/m3)
   if(nsmax>0) then
      rho_snow_surft(1, pft,:)        = jules53_struc(n)%jules53(t)%rho_snow(pft,:)             ! ntiles, nsmax, Snow layer densities (m)
   endif

   if ( l_triffid .or. l_phenol ) then
     wood_prod_fast_gb(1)       = jules53_struc(n)%jules53(t)%wood_prod_fast           ! Fast-turnover wood product C pool.
     wood_prod_med_gb(1)        = jules53_struc(n)%jules53(t)%wood_prod_med             ! Medium-turnover wood product C pool. 
     wood_prod_slow_gb(1)       = jules53_struc(n)%jules53(t)%wood_prod_slow            ! Slow-turnover wood product C pool.
     frac_agr_prev_gb(1)        = jules53_struc(n)%jules53(t)%frac_agr_prev             ! Agricultural fraction from previous TRIFFID call
     frac_past_prev_gb(1)       = jules53_struc(n)%jules53(t)%frac_past_prev
   endif
   
   n_inorg_soilt_lyrs(1,1,:) = jules53_struc(n)%jules53(t)%n_inorg_soilt_lyrs(:)        ! Gridbox Inorganic N pool on soil levels (kg N/m2)
   n_inorg_gb(1)              = jules53_struc(n)%jules53(t)%n_inorg                        ! Gridbox Inorganic N pool (kg N m-2)
   n_inorg_avail_pft(1,:,:)  = jules53_struc(n)%jules53(t)%n_inorg_avail_pft(:,:) 
   ns_pool_gb(1,:, :)                = jules53_struc(n)%jules53(t)%ns(:,:)                     ! **** dim_cslayer, dim_cs1, Soil Organic Nitrogen (kg N m-2)
   if(l_triffid) then
     triffid_co2_gb(1)             =jules53_struc(n)%jules53(t)%triffid_co2_gb !Atmospheric CO2 fluxes from TRIFFID (kgC/m2/yr)
   endif  
   
   canht_pft(1, 1:npft)          = jules53_struc(n)%jules53(t)%canht_ft(1:npft)               ! npft, Canopy height (m)
   lai_pft(1, 1:npft)            = jules53_struc(n)%jules53(t)%lai(1:npft)                    ! npft LAI of plant functional types
   canopy_surft(1, pft)            = jules53_struc(n)%jules53(t)%canopy(pft)                 ! ntiles, Surface/canopy water for snow-free land tiles (kg m-2)
   canopy_gb(1)            = jules53_struc(n)%jules53(t)%canopy_gb                 ! Gridbox canopy water content (kg m-2)
   
   cs_pool_soilt(1, 1,:,:)  = jules53_struc(n)%jules53(t)%cs_pool_soilt(:,:)                     ! dim_cs1, Soil carbon (kg C/m2)
   di_ncat_sicat(1,1, :)           = jules53_struc(n)%jules53(t)%di_ncat(:)                ! nice, "Equivalent thickness" of sea-ice catagories (m)
   k_sice_sicat(1,1, :)            = jules53_struc(n)%jules53(t)%k_sice(:)                 ! nice, Sea ice effective conductivity (2*kappai/de)
   
   gc_surft(1, pft)                = jules53_struc(n)%jules53(t)%gc_surft(pft)                     ! ntiles, Stomatal" conductance to evaporation for land tiles(m s-1)
   gs_gb(1)                   = jules53_struc(n)%jules53(t)%gs_gb                        !  "Stomatal" conductance to evaporation (m s-1)
   
   rgrain_surft(1, pft)            = jules53_struc(n)%jules53(t)%rgrain(pft)                 ! ntiles, Snow surface grain size on tiles (microns)
   smc_soilt(1,1)                  = jules53_struc(n)%jules53(t)%smc_soilt                       !  Soil moisture in a layer at the surface (kg m-2).
   smcl_soilt(1,1, :)              = jules53_struc(n)%jules53(t)%smcl_soilt(:)                   ! sm_levels, Soil moisture content of layers (kg m-2)
   snowdepth_surft(1, pft)         = jules53_struc(n)%jules53(t)%snowdepth(pft)              ! ntiles, Snow depth on ground on tiles (m)
   snow_surft(1, pft)         = jules53_struc(n)%jules53(t)%snow_tile(pft)              ! ntiles, Lying snow on tiles (kg m-2)
   snow_grnd_surft(1, pft)         = jules53_struc(n)%jules53(t)%snow_grnd(pft)              ! ntiles, Snow on the ground (kg m-2)
   snow_mass_ij(1,1)            = jules53_struc(n)%jules53(t)%snow_mass_ij                 !  Gridbox snowmass (kg m-2)
   snow_mass_sea_sicat(1,1,:)= jules53_struc(n)%jules53(t)%snow_mass_sea_sicat(:)     !nice_use,  Snow on category sea-ice (Kg/m2)
   soot_ij(1,1)                 = jules53_struc(n)%jules53(t)%soot_ij                      !  Snow soot content (kg kg-1)
   t_soil_soilt(1,1, :)            = jules53_struc(n)%jules53(t)%t_soil(:)                 ! sm_levels, Sub-surface temperatures (K)
   ti_sicat(1,1,:)                   = jules53_struc(n)%jules53(t)%ti_sicat(:)                        !  Sea-ice surface layer
   tstar_surft(1, :)        = jules53_struc(n)%jules53(t)%tstar_tile(:)             ! ntiles, Tile surface temperatures (K)
   tsurf_elev_surft(1, pft)   = jules53_struc(n)%jules53(t)%tsurf_elev_surft(pft)       ! Tiled land-ice bedrock subsurface temperatures (K)
   z0msea_ij(1,1)               = jules53_struc(n)%jules53(t)%z0msea                    !  Sea-surface roughness length for momentum (m).
end subroutine tile_to_prog

