module ac_defaultcropsoil

use ac_global , only:   DetermineParametersCR, &
                        GetCropFileFull, &
                        GetCrop_CCo, &
                        GetCrop_PlantingDens, &
                        GetCrop_RootMin, &
                        GetCrop_SizeSeedling, &
                        GetPathNameSimul, &
                        GetProfFileFull, &
                        GetSoilLayer_CRa, &
                        GetSoilLayer_CRb, &
                        GetSoilLayer_FC, &
                        GetSoilLayer_InfRate, &
                        GetSoilLayer_SAT, &
                        GetSoilLayer_SoilClass, &
                        GetSoilLayer_WP, &
                        ModeCycle_CalendarDays, &
                        NumberSoilClass, &
                        plant_Seed, &
                        pMethod_FAOCorrection, &
                        SaveCrop, &
                        SaveProfile, &
                        SetCropDescription, &
                        SetCropFileFull, &
                        SetCrop_aCoeff, &
                        SetCrop_AdaptedToCO2, &
                        SetCrop_AnaeroPoint, &
                        SetCrop_Assimilates_Mobilized, &
                        SetCrop_Assimilates_On, &
                        SetCrop_Assimilates_Period, &
                        SetCrop_Assimilates_Stored, &
                        SetCrop_bCoeff, &
                        SetCrop_CCEffectEvapLate, &
                        SetCrop_CCini, &
                        SetCrop_CCo, &
                        SetCrop_CCsaltDistortion, &
                        SetCrop_CCx, &
                        SetCrop_CCxRoot, &
                        SetCrop_CDC, &
                        SetCrop_CGC, &
                        SetCropFileFull, &
                        SetCrop_DaysToCCini, &
                        SetCrop_DaysToFlowering, &
                        SetCrop_DaysToGermination, &
                        SetCrop_DaysToHarvest, &
                        SetCrop_DaysToHIo, &
                        SetCrop_DaysToMaxRooting, &
                        SetCrop_DaysToSenescence, &
                        SetCrop_DeterminancyLinked, &
                        SetCrop_dHIdt, &
                        SetCrop_DHImax, &
                        SetCrop_DryMatter, &
                        SetCrop_ECemax, &
                        SetCrop_ECemin, &
                        SetCrop_fExcess, &
                        SetCrop_GDDaysToCCini, &
                        SetCrop_GDDCDC, &
                        SetCrop_GDDCGC, &
                        SetCrop_GDDaysToFlowering, &
                        SetCrop_GDDaysToGermination, &
                        SetCrop_GDDaysToHarvest, &
                        SetCrop_GDDaysToHIo, &
                        SetCrop_GDDaysToMaxRooting, &
                        SetCrop_GDDaysToSenescence, &
                        SetCrop_GDDLengthFlowering, &
                        SetCrop_GDtranspLow, &
                        SetCrop_HI, &
                        SetCrop_HIincrease, &
                        SetCrop_KcDecline, &
                        SetCrop_KcTop, &
                        SetCrop_KsShapeFactorLeaf, &
                        SetCrop_KsSHapeFactorSenescence, &
                        SetCrop_KsShapeFactorStomata, &
                        SetCrop_LengthFlowering, &
                        SetCrop_ModeCycle, &
                        SetCrop_Planting, &
                        SetCrop_pdef, &
                        SetCrop_PlantingDens, &
                        SetCrop_pLeafDefLL, &
                        SetCrop_pLeafDefUL, &
                        SetCrop_pMethod, &
                        SetCrop_pPollination, &
                        SetCrop_pSenescence, &
                        SetCrop_ResponseECsw, &
                        SetCrop_RootMax, &
                        SetCrop_RootMin, &
                        SetCrop_RootMinYear1, &
                        SetCrop_RootShape, &
                        SetCrop_SizePlant, &
                        SetCrop_SizeSeedling, &
                        SetCrop_SmaxBotQuarter, &
                        SetCrop_SmaxTopQuarter, &
                        SetCrop_SownYear1, &
                        SetCrop_StressResponse_Calibrated, &
                        SetCrop_StressResponse_ShapeCCX, &
                        SetCrop_StressResponse_ShapeCDecline, &
                        SetCrop_StressResponse_ShapeCGC, &
                        SetCrop_StressResponse_ShapeWP, &
                        SetCrop_StressResponse_Stress, &
                        SetCrop_subkind, &
                        SetCrop_SumEToDelaySenescence, &
                        SetCrop_Tbase, &
                        SetCrop_Tcold, &
                        SetCrop_Theat, &
                        SetCrop_Tupper, &
                        SetCrop_WP, &
                        SetCrop_WPy, &
                        SetCrop_YearCCx, &
                        SetProfFileFull, &
                        SetProfDescription, &
                        SetSoil_CNvalue, &
                        SetSoil_NrSoilLayers, &
                        SetSoil_REW, &
                        SetSoilLayer_CRa, &
                        SetSoilLayer_CRb, &
                        SetSoilLayer_FC, &
                        SetSoilLayer_Description, &
                        SetSoilLayer_GravelMass, &
                        SetSoilLayer_GravelVol, &
                        SetSoilLayer_InfRate, &
                        SetSoilLayer_Penetrability, &
                        SetSoilLayer_SAT, &
                        SetSoilLayer_SoilClass, &
                        SetSoilLayer_Thickness, &
                        SetSoilLayer_WP, &
                        subkind_Grain, &
                        undef_double, &
                        undef_int
use ac_kinds, only: int8, &
                    int16, &
                    int32, &
                    dp
implicit none


contains


subroutine ResetDefaultCrop(use_default_crop_file)
    logical, intent(in) :: use_default_crop_file
        !! Whether to write a 'DEFAULT.CRO' file.

    call SetCropDescription('a generic crop')
    call SetCrop_subkind(subkind_Grain)
    call SetCrop_Planting(plant_Seed)
    call SetCrop_SownYear1(.true.) ! for perennials
    call SetCrop_ModeCycle(ModeCycle_CalendarDays)
    call SetCrop_pMethod(pMethod_FAOCorrection)
    call SetCrop_Tbase(5.5_dp)  ! Basal temperature (degC)
    call SetCrop_Tupper(30.0_dp) ! Cut-off temperature (degC)
    call SetCrop_pLeafDefUL(0.25_dp) ! Soil water depletion factor for leaf
                                     ! expansion (p-leaf) - Upper Limit
    call SetCrop_pLeafDefLL(0.60_dp) ! Soil water depletion factor for leaf
                                     ! expansion (p-leaf) - Lower Limit
    call SetCrop_KsShapeFactorLeaf(3._dp) ! Shape factor for Water stress
                                          ! coefficient Leaf expansion
                                          ! (0 = straight line)
    call SetCrop_pdef(0.50_dp) ! Soil water depletion fraction for stomatal
                               ! control (p - stomatal)
    call SetCrop_KsShapeFactorStomata(3._dp) ! Shape factor for Water stress
                                        ! coefficient Stomatal Control
                                        ! (0 = straight line)
    call SetCrop_pSenescence(0.85_dp) ! Soil water depletion factor for
                                      ! canopy senescence (p-senescence)
    call SetCrop_KsShapeFactorSenescence(3._dp) ! Shape factor for Water
                                        ! stress coefficient Canopy
                                        ! Senescence (0 = straight line)
    call SetCrop_SumEToDelaySenescence(50) ! Sum(ETo) during stress period to be
                                      ! exceeded before senescence is triggered
    call SetCrop_pPollination(0.90_dp)
    call SetCrop_AnaeroPoint(5) ! Vol% for Anaerobiotic point (* (SAT - [vol%])
                                ! at which deficient aeration occurs

    call SetCrop_StressResponse_Stress(50_int8)  ! Soil fertility stress at
                                                 ! calibration (%)
    call SetCrop_StressResponse_ShapeCGC(2.16_dp)  ! Shape factor for response
                                                   ! of Canopy Growth Coefficient
                                                   ! to soil fertility stress
    call SetCrop_StressResponse_ShapeCCX(0.79_dp) ! Shape factor for response of
                              ! Maximum Canopy Coefficient to soil fertility stress
    call SetCrop_StressResponse_ShapeWP(1.67_dp)  ! Shape factor for response of
                              ! Crop water productivity to soil fertility stress
    call SetCrop_StressResponse_ShapeCDecline(1.67_dp) ! Shape factor for response
                              ! of Canopy cover decline to soil fertility stress
    call SetCrop_StressResponse_Calibrated(.true.)

    call SetCrop_ECemin(2_int8) ! Electrical Conductivity of soil saturation
                                ! extract at which crop starts to be affected by
                                ! soil salinity (dS/m)
    call SetCrop_ECemax(12_int8) ! Electrical Conductivity of soil saturation
                                 ! extract at which crop can no longer grow (dS/m)
    call SetCrop_CCsaltDistortion(25_int8) ! Calibrated distortion canopy cover
                                           ! for calibration of simulation of
                                           ! effect of salinity stress (%)
    call SetCrop_ResponseECsw(100) ! Calibrated response of Ks stomata to ECsw for
                                   ! calibration: From 0 (none) to +125 (very strong)
    call SetCrop_Tcold(8_int8) ! Minimum air temperature below which pollination
                               ! starts to fail (cold stress) (degC)
    call SetCrop_Theat(40_int8) ! Maximum air temperature above which pollination
                                ! starts to fail (heat stress) (degC)
    call SetCrop_GDtranspLow(11.1_dp) ! Minimum growing degrees required for full
                                      ! crop transpiration (degC - day)
    call SetCrop_KcTop(1.10_dp) ! Crop coefficient when complete cover and prior to
                                ! senescence (Kc,top)
    call SetCrop_KcDecline(0.150_dp) ! Decline crop coefficient (%/day) as a result
                                     ! of ageing, nitrogen defficiency, etc.
    call SetCrop_RootMin(0.30_dp) ! Minimum rooting depth (m)
    call SetCrop_RootMax(1.00_dp) ! Maximum rooting depth (m)
    call SetCrop_RootMinYear1(GetCrop_RootMin()) ! Minimum rooting depth in first
                                                 ! year in meter (for perennials)
    call SetCrop_RootShape(15_int8) ! Shape factor describing root zone expansion
    call SetCrop_SmaxTopQuarter(0.048_dp) ! Maximum root water extraction (m3water/m3soil.day)
                                          ! in top quarter of root zone
    call SetCrop_SmaxBotQuarter(0.012_dp) ! Maximum root water extraction (m3water/m3soil.day)
                                          ! in bottom quarter of root zone
    call SetCrop_CCEffectEvapLate(50) ! Effect of canopy cover on reduction soil evap
                                      ! in late season stage
    call SetCrop_SizeSeedling(6.50_dp) ! Canopy cover per seedling (cm2)
    call SetCrop_SizePlant(GetCrop_SizeSeedling()) ! Canopy cover when regrowth (cm2)
    call SetCrop_PlantingDens(185000) ! Number of plants per hectare
    call SetCrop_CCo((GetCrop_SizeSeedling()/10000) &
                        * (GetCrop_PlantingDens()/10000._dp)) ! Starting canopy size
                                                              ! (CCo) in fraction
    call SetCrop_CCini(GetCrop_CCo())
    call SetCrop_CGC(0.15_dp) ! Canopy growth coefficient (CGC): Increase in canopy
                              ! cover (in fraction) per day
    call SetCrop_YearCCx(int(undef_int, int8)) ! the number of years at which CCx declines to 90%
                                    ! of its value due to self-thining - Perennials
    call SetCrop_CCxRoot(real(undef_int, kind=dp)) ! shape factor of the decline of CCx over the
                                    ! years due to self-thinning - Perennials
    call SetCrop_CCx(0.80_dp) ! Maximum canopy cover (CCx) in fraction
    call SetCrop_CDC(0.1275_dp) ! Canopy decline coefficient (CDC): Decrease in
                                ! canopy cover (in fraction) per day
    call SetCrop_DaysToCCini(0)
    call SetCrop_DaysToGermination(5) ! Calendar Days: from sowing to germination
    call SetCrop_DaysToMaxRooting(100) ! Calendar Days: from sowing to maximum rooting depth
    call SetCrop_DaysToSenescence(110) ! Calendar Days: from sowing to start senescence
    call SetCrop_DaysToHarvest(125) ! Calendar Days: from sowing to maturity
    call SetCrop_DaysToFlowering(70) ! Calendar Days: from sowing to flowering
    call SetCrop_LengthFlowering(10) ! Length of the flowering stage (days)
    call SetCrop_DaysToHIo(50)
    call SetCrop_DeterminancyLinked(.true.)
    call SetCrop_fExcess(50_int16) ! Potential excess of fruits (%)
    call SetCrop_WP(17.0_dp) ! (normalized) Water productivity (gram/m2)
    call SetCrop_WPy(100) ! (normalized) Water productivity during yield formation
                          ! (Percent of WP)
    call SetCrop_AdaptedToCO2(100_int8) ! Percentage adapted to elevated atmospheric
                                       ! CO2 concentration
    call SetCrop_HI(50) ! HI harvest index (percentage)
    call SetCrop_DryMatter(25_int8) ! dry matter content (%) of fresh yield
    call SetCrop_HIincrease(5_int8) ! Possible increase (%) of HI due to water
                                    ! stress before flowering
    call SetCrop_aCoeff(10._dp) ! Coefficient describing positive impact of
                                ! restricted vegetative growth at flowering on HI
    call SetCrop_bCoeff(8._dp) ! Coefficient describing reduction of impact of
                               ! stomatal closure at flowering on HI
    call SetCrop_DHImax(15_int8) ! Allowable maximum increase (%) of specified HI
    call SetCrop_dHIdt(-9._dp) ! Caluclated as Crop_HI/Crop_DaysToHIo
    call SetCrop_GDDaysToCCini(-9)
    call SetCrop_GDDaysToGermination(-9) ! GDDays: from sowing to germination
    call SetCrop_GDDaysToMaxRooting(-9) ! GDDays: from sowing to maximum rooting depth
    call SetCrop_GDDaysToSenescence(-9) ! GDDays: from sowing to start senescence
    call SetCrop_GDDaysToHarvest(-9) ! GDDays: from sowing to harvest
    call SetCrop_GDDaysToFlowering(-9) ! GDDays: from sowing to flowering
    call SetCrop_GDDLengthFlowering(-9) ! Length of the flowering stage
                                        ! (growing degree days)
    call SetCrop_GDDaysToHIo(-9)
    call SetCrop_GDDCGC(-9.000000_dp) ! CGC for GGDays: Increase in canopy cover
                                      ! (in fraction) per growing-degree day
    call SetCrop_GDDCDC(-9.000000_dp) ! CDC for GGDays: Decrease in canopy cover
                                      ! (in fraction) growing-degree day

    call SetCrop_Assimilates_On(.false.)  ! transfer of assimilates from above ground
                                          ! parts to root system is NOT considered
    call SetCrop_Assimilates_Period(0)  ! Number of days before end of season at
                                        ! which storage starts
    call SetCrop_Assimilates_Stored(0_int8)  ! Percentage of assimilates transferred
                                             ! to root system at end of season
    call SetCrop_Assimilates_Mobilized(0_int8) ! Percentage stored assimilates,
                                    ! transferred to above ground parts in next season
    if (use_default_crop_file) then
        call SetCropFilefull(GetPathNameSimul() // 'DEFAULT.CRO')
        call SaveCrop(GetCropFilefull())
    end if
end subroutine ResetDefaultCrop


subroutine ResetDefaultSoil(use_default_soil_file)
    logical, intent(in) :: use_default_soil_file
        !! Whether to write a 'DEFAULT.SOL' file.

    real(dp) :: cra_temp, crb_temp
    character(len=25) :: TempString

    call SetProfDescription('deep loamy soil profile')
    call SetSoil_CNvalue(61_int8) ! for an initial abstraction of 0.05 S
    call SetSoil_REW(9_int8)
    call SetSoil_NrSoilLayers(1_int8)
    call SetSoilLayer_Thickness(1, 4._dp)
    call SetSoilLayer_SAT(1, 50.0_dp)
    call SetSoilLayer_FC(1, 30.0_dp)
    call SetSoilLayer_WP(1, 10.0_dp)
    call SetSoilLayer_InfRate(1, 500.0_dp)
    call SetSoilLayer_Penetrability(1, 100_int8)
    call SetSoilLayer_GravelMass(1, 0_int8)
    call SetSoilLayer_GravelVol(1, 0._dp)
    TempString = 'Loamy soil horizon'
    call SetSoilLayer_Description(1, TempString)
    call SetSoilLayer_SoilClass(1, &
                NumberSoilClass(GetSoilLayer_SAT(1), GetSoilLayer_FC(1), &
                                GetSoilLayer_WP(1), GetSoilLayer_InfRate(1)))
    cra_temp = GetSoilLayer_CRa(1)
    crb_temp = GetSoilLayer_CRb(1)
    call DetermineParametersCR(GetSoilLayer_SoilClass(1), &
                                GetSoilLayer_InfRate(1), cra_temp, crb_temp)
    call SetSoilLayer_CRa(1, cra_temp)
    call SetSoilLayer_CRb(1, crb_temp)

    if (use_default_soil_file) then
        call SetProfFilefull(GetPathNameSimul() // 'DEFAULT.SOL')
        call SaveProfile(GetProfFilefull())
    end if
end subroutine ResetDefaultSoil

end module ac_defaultcropsoil
