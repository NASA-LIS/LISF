!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine tile_to_ps(n, t, pft)
  use jules54_lsmMod
  use p_s_parms
  use jules_vegetation_mod,     only: l_use_pft_psi
  use jules_soil_biogeochem_mod, only: soil_model_ecosse, soil_bgc_model
  implicit none 
  integer, intent(in) :: n, t, pft

  bexp_soilt(1,1, :)           =  jules54_struc(n)%jules54(t)%p_s_b(:)                !  Exponent for soil moisture characteristic function Clapp-Hornberger model: b is the Clapp-Hornberger exponent,  van Genuchten model: b=1/(n-1)  (metres)
  sathh_soilt(1,1, :)          =  jules54_struc(n)%jules54(t)%p_s_sathh(:)            !  Parameter for soil moisture characteristic functions Clapp-Hornberger model: sathh is the saturated soil water pressure (m), van Genuchten model: sathh=1/alpha
  hcap_soilt(1,1, :)           =  jules54_struc(n)%jules54(t)%p_s_hcap(:)             !  Soil heat capacity (J/K/m3)
  hcon_soilt(1,1, :)           =  jules54_struc(n)%jules54(t)%p_s_hcon(:)             !  Soil thermal conductivity (W/m/K)
  satcon_soilt(1,1, :)         =  jules54_struc(n)%jules54(t)%p_s_satcon(:)           !  Saturated hydraulic conductivity (kg m-2/s)
  smvccl_soilt(1,1, :)         =  jules54_struc(n)%jules54(t)%p_s_smvccl(:)           !  Critical volumetric SMC (cubic m per cubic m of soil)    
  smvcst_soilt(1,1, :)         =  jules54_struc(n)%jules54(t)%p_s_smvcst(:)           !  Volumetric saturation point (m^3 m-3 of soil)
  smvcwt_soilt(1,1, :)         =  jules54_struc(n)%jules54(t)%p_s_smvcwt(:)           !  Volumetric wilting point (cubic m per cubic m of soil)

  if(l_use_pft_psi .and. pft .le. jules54_struc(n)%npft) then
    v_close_pft(1,:,pft)           =  jules54_struc(n)%jules54(t)%p_s_v_close_pft(:,pft)
    v_open_pft(1,:,pft)            =  jules54_struc(n)%jules54(t)%p_s_v_open_pft(:,pft)
  endif
  clay_soilt(1,1,:)            =  jules54_struc(n)%jules54(t)%p_s_clay_soilt(:)       ! modified for JULES 5.2 

  albsoil_soilt(1,1)           =  jules54_struc(n)%jules54(t)%p_s_albsoil             !  Soil albedo
  albobs_sw_gb(1)              =  jules54_struc(n)%jules54(t)%p_s_albobs_sw           !  Obs SW albedo    
  albobs_vis_gb(1)             =  jules54_struc(n)%jules54(t)%p_s_albobs_vis          !  Obs VIS albedo    
  albobs_nir_gb(1)             =  jules54_struc(n)%jules54(t)%p_s_albobs_nir          !  Obs NIR albedo
  
  catch_surft(1, pft)          =  jules54_struc(n)%jules54(t)%p_s_catch(pft)            !  Surface/canopy water capacity of snow-free land tiles (kg m-2)
  catch_snow_surft(1, pft)     =  jules54_struc(n)%jules54(t)%p_s_catch_snow(pft)       !  Snow interception capacity (kg m-2)
  cosz_ij(1,1)                 =  jules54_struc(n)%jules54(t)%p_s_cosz                !  Cosine of the zenith angle    
  
  infil_surft(1, pft)          =  jules54_struc(n)%jules54(t)%p_s_infil_tile(pft)       !  Maximum possible surface infiltration for tiles (kg m-2/s)
  z0_surft(1, pft)             =  jules54_struc(n)%jules54(t)%p_s_z0_tile(pft)          !  Surface roughness on tiles (m).
  z0h_bare_surft(1, pft)       =  jules54_struc(n)%jules54(t)%p_s_z0h_tile_bare(pft)    !  Surface thermal roughness on tiles before allowance for snow cover (m).
  
  sthu_soilt(1, 1,:)           =  jules54_struc(n)%jules54(t)%p_s_sthu(:)             !  Unfrozen soil moisture content of the layers as a fraction of saturation.
  sthf_soilt(1, 1,:)           =  jules54_struc(n)%jules54(t)%p_s_sthf(:)             !  Frozen soil moisture content of the layers as a fraction of saturation.
  
  if ( soil_bgc_model == soil_model_ecosse ) then
    !clay_frac_soilt(1,1,:)       =  jules54_struc(n)%jules54(t)%p_s_clay_frac_soilt(:)
    soil_ph_soilt(1,1,:)         =  jules54_struc(n)%jules54(t)%p_s_soil_ph_soilt(:)
  endif 
  sthu_min_soilt(1,1, :)       =  jules54_struc(n)%jules54(t)%p_s_sthu_min(:)         ! Minimum unfrozen water content for each layer. Used to normalise thaw depth calculation based on unfrozen water content fraction.
  z0m_soil_gb(1)               =  jules54_struc(n)%jules54(t)%p_s_z0m_soil_gb
end subroutine tile_to_ps 

