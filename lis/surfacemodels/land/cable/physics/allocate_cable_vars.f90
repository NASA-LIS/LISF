  SUBROUTINE allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssoil, &
       sum_flux,veg,model_structure,arraysize)
    use cable_types
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(OUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)    :: canopy
    TYPE (roughness_type), INTENT(OUT) :: rough
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal
    TYPE (model_structure_type), INTENT(INOUT) :: model_structure
    INTEGER, INTENT(IN) :: arraysize

    ! First nullify all pointers:
    CALL dealloc_cbm_var(air, arraysize, model_structure)
    CALL dealloc_cbm_var(bgc, arraysize, model_structure)
    CALL dealloc_cbm_var(canopy, arraysize, model_structure)
    CALL dealloc_cbm_var(met, arraysize, model_structure)
    CALL dealloc_cbm_var(bal, arraysize, model_structure)
    CALL dealloc_cbm_var(rad, arraysize, model_structure)
    CALL dealloc_cbm_var(rough, arraysize, model_structure)
    CALL dealloc_cbm_var(soil, arraysize, model_structure)
    CALL dealloc_cbm_var(ssoil, arraysize, model_structure)
    CALL dealloc_cbm_var(sum_flux, arraysize, model_structure)
    CALL dealloc_cbm_var(veg, arraysize, model_structure)

    ! Allocate CABLE's main variables:
    CALL alloc_cbm_var(air, arraysize, model_structure)
    CALL alloc_cbm_var(bgc, arraysize, model_structure)
    CALL alloc_cbm_var(canopy, arraysize, model_structure)
    CALL alloc_cbm_var(met, arraysize, model_structure)
    CALL alloc_cbm_var(bal, arraysize, model_structure)
    CALL alloc_cbm_var(rad, arraysize, model_structure)
    CALL alloc_cbm_var(rough, arraysize, model_structure)
    CALL alloc_cbm_var(soil, arraysize, model_structure)
    CALL alloc_cbm_var(ssoil, arraysize, model_structure)
    CALL alloc_cbm_var(sum_flux, arraysize, model_structure)
    CALL alloc_cbm_var(veg, arraysize, model_structure)

    ! Allocate patch fraction variable:
!    ALLOCATE(patch(arraysize))

  END SUBROUTINE allocate_cable_vars
