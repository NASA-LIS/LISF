!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: NoahMP50_setgeosexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from NoahMP to GEOS in coupled mode
!
! !INTERFACE:
subroutine NoahMP50_setgeosexport(n)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_logMod
  use LISGEOSGridCompMod, only : LISGEOS_export
  use NoahMP50_lsmMod
  
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
!EOP
  integer               :: i
  real, allocatable     :: temp(:)

  allocate(temp(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  ! Ensure the GEOS-coupling derived-variable store is allocated, then
  ! refresh it from the latest NoahMP state (idempotent allocation; safe
  ! to call every coupling step).
  !!call NoahMP50_struc_init(n)
  call NoahMP50_convert_to_GEOSexport(n)

  ! =====================================================================
  ! 1. BASELINE VARIABLES (Carried over from original WRF coupling)
  ! =====================================================================

  ! Surface Temperature
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%tsfcmean
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%avgsurft_t,temp)

  ! Emissivity
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%emiss
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%emiss_t,temp)

  ! Latent Heat Flux
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%lh
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%qle_t,temp)

  ! Surface Albedo
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%albedo
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%albedo_t,temp)

  ! Soil Moisture Content (Layers 1-4)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%smc(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%smc1_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%smc(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%smc2_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%smc(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%smc3_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%smc(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%smc4_t,temp)

  ! Soil Temperature (Layers 1-4)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%tslb(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%stc1_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%tslb(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%stc2_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%tslb(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%stc3_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%tslb(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%stc4_t,temp)

  ! Volumetric Liquid Soil Moisture (Layers 1-4)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%sh2o(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%sh2o1_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%sh2o(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%sh2o2_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%sh2o(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%sh2o3_t,temp)
  
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%sh2o(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%sh2o4_t,temp)

  ! Snow Water Equivalent 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%sneqv 
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%snow_t,temp)

  ! Snow Height
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%snowh
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%snowh_t,temp)

#ifdef WRF_HYDRO
  ! Infiltration Excess
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%infxs1rt
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%infxsrt_t,temp)

  ! Soil Drainage
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%soldrain1rt
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%soldrain_t,temp)
#endif


  ! =====================================================================
  ! 2. DERIVED VARIABLES 
  ! =====================================================================

  ! Sensible Heat (qh_shout) - native NoahMP output, not a weighted/derived
  ! quantity; sourced directly from the struct (NoahmpIO%hfx).
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%hfx
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%qh_t,temp)

  ! Exchange Coefficient (ch_f)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%ch_f
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%ch_f_t,temp)

  ! Momentum Exchange Coefficient (cm_f)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%cm_f
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%cm_f_t,temp)

  ! 2m Specific Humidity (qair_t2m) - native NoahMP output, not a weighted/
  ! derived quantity; sourced directly from the struct (NoahmpIO%qsfc).
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%qsfcmean
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%qair_t2m_t,temp)

  ! Snow Depth - Snow-covered-area form (snowdepth_sca)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%snowdepth_sca
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%snowdepth_sca_t,temp)

  ! Total Runoff (totalrunoff)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%totalrunoff
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%totalrunoff_t,temp)

  ! Root-zone Soil Moisture (soilmoist_rz)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%soilmoist_rz
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%rootmoist_t,temp)

  ! Albedo Visible Direct (band 1)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%albedo_visdir
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%albedovisdir_t,temp)

  ! Albedo Near-Infrared Direct (band 2)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%albedo_nirdir
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%albedonirdir_t,temp)

  ! Albedo Visible Diffuse (band 1)      
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%albedo_visdif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%albedovisdif_t,temp)

   ! Albedo Near-Infrared Diffuse (band 2) 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%albedo_nirdif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%albedonirdif_t,temp)

  ! 2m Air Temperature (tair_t2m) - native NoahMP output, not a weighted/
  ! derived quantity; sourced directly from the struct (NoahmpIO%t2mxy).
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%t2mean
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%tair_t2m_t,temp)

  ! Canopy Water Content (canopint)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%canopint
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%canopint_t,temp)


  ! =====================================================================
  ! 3. NATIVE VARIABLES 
  ! =====================================================================

  ! Surface Runoff (Qs)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%sfcrunoff
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%qs_t,temp)

  ! Subsurface Runoff (Qsb)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%udrrunoff
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%qsb_t,temp)

  ! Surface Air Density (rhoair)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%rhoair
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%rhos_t,temp)

  ! Snow Cover Fraction
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%snowc
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%snocvr_t,temp)

  ! ZPD (Zero Plane Displacement)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%zpd
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%zpd_t,temp)

  ! Z0 (Surface Roughness Length)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%z0
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%z0_t,temp)

  ! Taux
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%taux
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%taux_t,temp)

  ! Tauy
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%tauy
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%tauy_t,temp)

  ! Friction Velocity - Vegetated (UstarVeg)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%ustarveg
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%ustarveg_t,temp)

  ! Friction Velocity - Bare Soil (UstarBare)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%ustarbare
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%ustarbare_t,temp)

  ! Frication Velocity - grid mean (Ustar)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%ustar
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%ustar_t,temp)

  ! Root Mass (rtmass)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%rtmass
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%rtmass_t,temp)

  ! Stem Mass (stmass)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%stmass
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%stmass_t,temp)

  ! Wood (wood)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%wood
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%wood_t,temp)

  ! Stable Carbon Pool (stblcp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%stblcp
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%stblcp_t,temp)

  ! Fast Carbon Pool (fastcp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%fastcp
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%fastcp_t,temp)

  ! Leaf Area Index (lai)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%lai
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%lai_t,temp)

  ! Stem Area Index (sai)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%sai
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%sai_t,temp)

  ! Vegetation Area Index (vai = lai + sai)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%vai
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%vai_t,temp)

  ! Total Evapotranspiration (evap)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%qfx
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%evap_t,temp)

  ! Greenness Fraction (greenness)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%fveg
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%greenness_t,temp)

  ! Vegetation/Land Cover Type (landcover)
  ! Raw NoahMP vegetation-type passthrough; no NoahMP->CLSM land-cover
  ! category remap table exists yet, so this is exported unconverted.
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = real(NoahMP50_struc(n)%noahmp50(i)%vegetype)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%landcover_t,temp)

  ! Net Ecosystem Exchange (nee)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NoahMP50_struc(n)%noahmp50(i)%nee
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISGEOS_export(n)%nee_t,temp)

  deallocate(temp)

end subroutine NoahMP50_setgeosexport