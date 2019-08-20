module cable_arrays
  
  use cable_types
  
  implicit none
  
  type(air_type)             :: air             ! air property variables
  type(bgc_pool_type)        :: bgc             ! carbon pool variables
  type(canopy_type)          :: canopy          ! vegetation variables
  type(met_type)             :: met             ! met input variables
  type(balances_type)        :: bal             ! energy and water balance variables
  type(radiation_type)       :: rad             ! radiation variables
  type(roughness_type)       :: rough           ! roughness varibles
  type(soil_parameter_type)  :: soil            ! soil parameters
  type(soil_snow_type)       :: ssoil           ! soil and snow variables
  type(sum_flux_type)        :: sum_flux        ! cumulative flux variables
  type(veg_parameter_type)   :: veg             ! vegetation parameters
  type(model_structure_type) :: model_structure ! which parametrisations to use
  
end module cable_arrays
