!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine ps_to_tile(n, t, pft)
  use jules53_lsmMod
  use p_s_parms
  use jules_vegetation_mod,     only: l_use_pft_psi
  use jules_soil_biogeochem_mod, only: soil_model_ecosse, soil_bgc_model
  implicit none 
  integer, intent(in) :: n, t, pft

  jules53_struc(n)%jules53(t)%p_s_b(:)             = bexp_soilt(1,1, :)                 !  Exponent for soil moisture characteristic function Clapp-Hornberger model: b is the Clapp-Hornberger exponent,  van Genuchten model: b=1/(n-1)  (metres)
  jules53_struc(n)%jules53(t)%p_s_sathh(:)         = sathh_soilt(1,1, :)             !  Parameter for soil moisture characteristic functions Clapp-Hornberger model: sathh is the saturated soil water pressure (m), van Genuchten model: sathh=1/alpha
  jules53_struc(n)%jules53(t)%p_s_hcap(:)          = hcap_soilt(1,1, :)              !  Soil heat capacity (J/K/m3)
  jules53_struc(n)%jules53(t)%p_s_hcon(:)          = hcon_soilt(1,1, :)              !  Soil thermal conductivity (W/m/K)
  jules53_struc(n)%jules53(t)%p_s_satcon(:)        = satcon_soilt(1,1, :)            !  Saturated hydraulic conductivity (kg m-2/s)
  jules53_struc(n)%jules53(t)%p_s_smvccl(:)        = smvccl_soilt(1,1, :)            !  Critical volumetric SMC (cubic m per cubic m of soil)    
  jules53_struc(n)%jules53(t)%p_s_smvcst(:)        = smvcst_soilt(1,1, :)            !  Volumetric saturation point (m^3 m-3 of soil)
  jules53_struc(n)%jules53(t)%p_s_smvcwt(:)        = smvcwt_soilt(1,1, :)            !  Volumetric wilting point (cubic m per cubic m of soil)

  if(l_use_pft_psi .and. pft .le. jules53_struc(n)%npft) then
    jules53_struc(n)%jules53(t)%p_s_v_close_pft(:,:) = v_close_pft(1,:,:)
    jules53_struc(n)%jules53(t)%p_s_v_open_pft(:,:)  = v_open_pft(1,:,:) 
  endif
  jules53_struc(n)%jules53(t)%p_s_clay_soilt(:)    = clay_soilt(1,1, :)              ! modified for JULES 5.2

  jules53_struc(n)%jules53(t)%p_s_albsoil          = albsoil_soilt(1,1)              !  Soil albedo
  jules53_struc(n)%jules53(t)%p_s_albobs_sw        = albobs_sw_gb(1)            !  Obs SW albedo    
  jules53_struc(n)%jules53(t)%p_s_albobs_vis       = albobs_vis_gb(1)           !  Obs VIS albedo    
  jules53_struc(n)%jules53(t)%p_s_albobs_nir       = albobs_nir_gb(1)           !  Obs NIR albedo
  
  jules53_struc(n)%jules53(t)%p_s_catch(pft)       = catch_surft(1, pft)             !  Surface/canopy water capacity of snow-free land tiles (kg m-2)
  jules53_struc(n)%jules53(t)%p_s_catch_snow(pft)  = catch_snow_surft(1, pft)        !  Snow interception capacity (kg m-2)
  jules53_struc(n)%jules53(t)%p_s_cosz             = cosz_ij(1,1)                 !  Cosine of the zenith angle    
 
  jules53_struc(n)%jules53(t)%p_s_infil_tile(pft)  = infil_surft(1, pft)        !  Maximum possible surface infiltration for tiles (kg m-2/s)
  jules53_struc(n)%jules53(t)%p_s_z0_tile(pft)       = z0_surft(1, pft)           !  Surface roughness on tiles (m).
  jules53_struc(n)%jules53(t)%p_s_z0h_tile_bare(pft) = z0h_bare_surft(1, pft)     !  Surface thermal roughness on tiles before allowance for snow cover (m).
 
  jules53_struc(n)%jules53(t)%p_s_sthu(:)          = sthu_soilt(1,1, :)              !  Unfrozen soil moisture content of the layers as a fraction of saturation.
  jules53_struc(n)%jules53(t)%p_s_sthf(:)          = sthf_soilt(1,1, :)              !  Frozen soil moisture content of the layers as a fraction of saturation.
  
  if ( soil_bgc_model == soil_model_ecosse ) then
    !jules53_struc(n)%jules53(t)%p_s_clay_frac_soilt(:)     = clay_frac_soilt(1,1,:)          !  Soil clay fraction
    jules53_struc(n)%jules53(t)%p_s_soil_ph_soilt(:)     = soil_ph_soilt(1,1,:)            ! Soil pH, defined on soil layers.  
  endif 
  jules53_struc(n)%jules53(t)%p_s_sthu_min(:)      = sthu_min_soilt(1,1, :)          ! Minimum unfrozen water content for each layer. Used to normalise thaw depth calculation based on unfrozen water content fraction.
  jules53_struc(n)%jules53(t)%p_s_z0m_soil_gb          = z0m_soil_gb(1)  ! Bare soil roughness, for momentum (m). 
end subroutine ps_to_tile 
