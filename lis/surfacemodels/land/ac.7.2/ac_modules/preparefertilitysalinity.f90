module ac_preparefertilitysalinity

use ac_global, only: CO2ref,& 
                     CropStressParametersSoilFertility, & 
                     datatype_daily, & 
                     datatype_decadely, & 
                     DaysInMonth, &  
                     DaysToReachCCwithGivenCGC, & 
                     DegreesDay, & 
                     DetermineDate,& 
                     DetermineDayNr, & 
                     DetermineDayNr, & 
                     FileExists, & 
                     GetCO2FileFull, & 
                     GetDaySwitchToLinear, & 
                     GetPathNameSimul, & 
                     GetSimulParam_GDDMethod, & 
                     GetTemperatureFile, & 
                     GetTemperatureFilefull, & 
                     GetTemperatureRecord_DataType, & 
                     GetTemperatureRecord_FromDayNr, & 
                     GetTemperatureRecord_FromM, & 
                     GetTemperatureRecord_FromY, & 
                     GetTemperatureRecord_ToDayNr, & 
                     GetTemperatureRecord_ToM, & 
                     GetTemperatureRecord_ToY, & 
                     GetTmaxCropReferenceRun, & 
                     GetTmaxCropReferenceRun_i, & 
                     GetTmaxTnxReference365DaysRun,& 
                     GetTmaxTnxReference365DaysRun_i,& 
                     GetTminCropReferenceRun, & 
                     GetTminCropReferenceRun_i, & 
                     GetTminTnxReference365DaysRun,& 
                     GetTminTnxReference365DaysRun_i,& 
                     GetTnxReferenceFile, & 
                     GetTnxReferenceFileFull, & 
                     GetTnxReferenceYear,& 
                     HarvestIndexGrowthCoefficient, & 
                     LeapYear, & 
                     modeCycle_CalendarDays, & 
                     modeCycle_GDDays, & 
                     rep_DayEventDbl,& 
                     rep_EffectStress,& 
                     rep_Shapes,& 
                     SeasonalSumOfKcPot, & 
                     SetSimulation_DelayedDays, & 
                     SetTmaxCropReferenceRun, & 
                     SetTmaxCropReferenceRun_i, & 
                     SetTmaxTnxReference365DaysRun,& 
                     SetTmaxTnxReference365DaysRun_i,& 
                     SetTminCropReferenceRun, & 
                     SetTminCropReferenceRun_i, & 
                     SetTminTnxReference365DaysRun,& 
                     SetTminTnxReference365DaysRun_i,& 
                     SetTnxReferenceFile, & 
                     SetTnxReferenceFileFull, & 
                     SplitStringInTwoParams, & 
                     Subkind_Forage, & 
                     subkind_Grain,  & 
                     subkind_Grain, & 
                     subkind_Tuber, & 
                     subkind_Vegetative, &
                     SumCalendarDaysReferenceTnx, &
                     TimeToMaxCanopySF, & 
                     undef_double, & 
                     undef_int
use ac_kinds, only: int8, &
                    int16, &
                    int32, &
                    intEnum, &
                    sp
use ac_tempprocessing, only: Bnormalized, &
                        CropStressParametersSoilSalinity, &
                        fTnxReference_open, &
                        fTnxReference_write, & 
                        fTnxReference_close,&
                        GDDCDCToCDC, &
                        GrowingDegreeDays
use ac_project_input, only: GetNumberSimulationRuns, &
                            ProjectInput
use ac_utils, only: GetReleaseDate, &
                    GetVersionString, &
                    roundc, &
                    trunc
use iso_fortran_env, only: iostat_end
implicit none


contains



subroutine AdjustCalendarDaysReferenceTnx(PlantDayNr, TheCropType, &
                    Tbase, Tupper, TDayMin, TDayMax, &
                    GDDL0, GDDL12, GDDFlor, GDDLengthFlor, GDDL123, &
                    GDDL1234, GDDHImax, GDDCGC, GDDCDC, &
                    CCo, CCx, RefHI, TheDaysToCCini, TheGDDaysToCCini, &
                    ThePlanting, L0, L12, LFlor, LengthFlor, &
                    L123, L1234, LHImax, CGC, CDC, RatedHIdt)
    integer(int32), intent(in) :: PlantDayNr
    integer(intEnum), intent(in) :: TheCropType
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDFlor
    integer(int32), intent(in) :: GDDLengthFlor
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    integer(int32), intent(in) :: GDDHImax
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    integer(int32), intent(in) :: RefHI
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(intEnum), intent(in) :: ThePlanting
    integer(int32), intent(inout) :: L0
    integer(int32), intent(inout) :: L12
    integer(int32), intent(inout) :: LFlor
    integer(int32), intent(inout) :: LengthFlor
    integer(int32), intent(inout) :: L123
    integer(int32), intent(inout) :: L1234
    integer(int32), intent(inout) :: LHImax
    real(sp), intent(inout) :: CGC
    real(sp), intent(inout) :: CDC
    real(sp), intent(inout) :: RatedHIdt

    integer(int32) :: ExtraGDDays, ExtraDays

    if (TheDaysToCCini == 0) then
        ! planting/sowing
        L0 = SumCalendarDaysReferenceTnx(GDDL0, PlantDayNr, PlantDayNr, &
            Tbase, Tupper, TDayMin, TDayMax)
        L12 = SumCalendarDaysReferenceTnx(GDDL12, PlantDayNr, PlantDayNr, &
            Tbase, Tupper, TDayMin, TDayMax)
    else
        ! regrowth
        if (TheDaysToCCini > 0) then
            ! CCini < CCx
            ExtraGDDays = GDDL12 - GDDL0 - TheGDDaysToCCini
            ExtraDays = SumCalendarDaysReferenceTnx(ExtraGDDays, PlantDayNr, &
                PlantDayNr, Tbase, Tupper, TDayMin, TDayMax)
            L12 = L0 + TheDaysToCCini + ExtraDays
        end if
    end if
    if (TheCropType /= subkind_Forage) then
        L123 = SumCalendarDaysReferenceTnx(GDDL123, PlantDayNr, PlantDayNr, &
            Tbase, Tupper, TDayMin, TDayMax)
        L1234 = SumCalendarDaysReferenceTnx(GDDL1234, PlantDayNr, PlantDayNr, &
            Tbase, Tupper, TDayMin, TDayMax)
    end if

    select case (TheCropType)
    case (subkind_Grain, subkind_Tuber)
        LFlor = SumCalendarDaysReferenceTnx(GDDFlor, PlantDayNr, PlantDayNr, &
            Tbase, Tupper, TDayMin, TDayMax)
        if (TheCropType == subkind_Grain) then
            LengthFlor = SumCalendarDaysReferenceTnx(GDDLengthFlor, PlantDayNr, &
                (PlantDayNr+LFlor), Tbase, Tupper, TDayMin, TDayMax)
        else
            LengthFlor = 0
        end if
        LHImax = SumCalendarDaysReferenceTnx(GDDHImax, PlantDayNr, &
            (PlantDayNr+LFlor), Tbase, Tupper, TDayMin, TDayMax)
    case (subkind_Vegetative, subkind_Forage)
        LHImax = SumCalendarDaysReferenceTnx(GDDHImax, PlantDayNr, &
            PlantDayNr, Tbase, Tupper, TDayMin, TDayMax)
    end select

    CGC = (real(GDDL12, kind=sp)/real(L12, kind=sp)) * GDDCGC
    call GDDCDCToCDC(PlantDayNr, L123, GDDL123, GDDL1234, CCx, GDDCDC, &
        Tbase, Tupper, TDayMin, TDayMax, CDC, .true.)
    if ((TheCropType == subkind_Grain) .or. (TheCropType == subkind_Tuber)) then
        RatedHIdt = real(RefHI, kind=sp)/real(LHImax, kind=sp)
    end if
    if ((TheCropType == subkind_Vegetative) .or. (TheCropType == subkind_Forage)) then
        if (LHImax > 0) then
            if (LHImax > L1234) then
                RatedHIdt = real(RefHI, kind=sp)/real(L1234, kind=sp)
            else
                RatedHIdt = real(RefHI, kind=sp)/real(LHImax, kind=sp)
            end if
            if (RatedHIdt > 100._sp) then
                RatedHIdt = 100._sp ! 100 is maximum TempdHIdt (See SetdHIdt)
                LHImax = 0
            end if
        else
            RatedHIdt = 100._sp ! 100 is maximum TempdHIdt (See SetdHIdt)
            LHImax = 0
        end if
    end if
end subroutine AdjustCalendarDaysReferenceTnx


subroutine DailyTnxReferenceFileCoveringCropPeriod(CropFirstDay)
    integer(int32), intent(in) :: CropFirstDay

    integer(int32) :: DayNr1
    integer(int32) :: Dayi, Monthi, Yeari, i
    real(sp) :: Tlow, Thigh
    character(len=1025) :: TempString

    if (FileExists(GetTnxReferenceFileFull()) &
        .or.(GetTnxReferenceFile() == '(External)')) then
        ! CropFirstDay = DayNr1 in undefined year
        call DetermineDate(CropFirstDay, Dayi, Monthi, Yeari)
        call DetermineDayNr(Dayi, Monthi, (1901), DayNr1)
        
        if (GetTnxReferenceFile() /= '(External)') then
            ! create SIM file
            call fTnxReference_open(trim(GetPathNameSimul()) // 'TCropReference.SIM','w')
        end if
        
        i = 0
        do Dayi = DayNr1, 365
            i=i+1
            Tlow = GetTminTnxReference365DaysRun_i(Dayi)
            Thigh = GetTmaxTnxReference365DaysRun_i(Dayi)
            call SetTminCropReferenceRun_i(i,Tlow)
            call SetTmaxCropReferenceRun_i(i,Thigh)
            if (GetTnxReferenceFile() /= '(External)') then
                write(TempString, '(f10.2, f10.2)') Tlow,Thigh
                call fTnxReference_write(trim(TempString))
            end if
        end do
        do Dayi = 1, (DayNr1-1)
            i=i+1
            Tlow = GetTminTnxReference365DaysRun_i(Dayi)
            Thigh = GetTmaxTnxReference365DaysRun_i(Dayi)
            call SetTminCropReferenceRun_i(i,Tlow)
            call SetTmaxCropReferenceRun_i(i,Thigh)
            if (GetTnxReferenceFile() /= '(External)') then
                write(TempString, '(f10.2, f10.2)') Tlow,Thigh
                call fTnxReference_write(trim(TempString))
            end if
        end do
        
        if (GetTnxReferenceFile() /= '(External)') then
            ! Close files
            call fTnxReference_close()
        end if
    end if
end subroutine DailyTnxReferenceFileCoveringCropPeriod


real(sp) function CO2ForTnxReferenceYear(TnxReferenceYear)
    integer(int32), intent(in) :: TnxReferenceYear

    integer(int32) :: i
    integer(int32) :: fhandle, rc
    character(len=1025) :: TempString
    real(sp) :: TheCO2, CO2a, CO2b, YearA, YearB

    if (TnxReferenceYear == 2000) then
        TheCO2 = CO2Ref
    else
        open(newunit=fhandle, file=trim(GetCO2FileFull()), status='old', &
                                                    action='read',iostat=rc)
        do i= 1, 3
            read(fhandle, *, iostat=rc) ! Description and Title
        end do
        ! from year
        read(fhandle, '(a)', iostat=rc) TempString
        call SplitStringInTwoParams(trim(TempString), YearB, CO2b)
        if (roundc(YearB, mold=1) >= TnxReferenceYear) then
            TheCO2 = CO2b
        else
            loop: do
                YearA = YearB
                CO2a = CO2b
                read(fhandle, '(a)', iostat=rc) TempString
                call SplitStringInTwoParams(trim(TempString), YearB, CO2b)
                if ((roundc(YearB, mold=1) >= TnxReferenceYear) .or. (rc == iostat_end)) exit loop
            end do loop
            if (TnxReferenceYear > roundc(YearB, mold=1)) then
                TheCO2 = CO2b
            else
                TheCO2 = CO2a + (CO2b-CO2a)*(TnxReferenceYear - &
                    roundc(YearA, mold=1))/(roundc(YearB, mold=1) - &
                    roundc(YearA, mold=1))
            end if
        end if
        Close(fhandle)
    end if
    CO2ForTnxReferenceYear = TheCO2
end function CO2ForTnxReferenceYear


subroutine StressBiomassRelationshipForTnxReference(TheDaysToCCini, TheGDDaysToCCini,&
            L0, L12, L123, L1234, LFlor, LengthFlor, GDDL0, GDDL12,&
            GDDL123, GDDL1234, WPyield, RefHI, CCo, CCx, CGC, GDDCGC,&
            CDC, GDDCDC, KcTop, KcDeclAgeing, CCeffectProcent,&
            Tbase, Tupper, TDayMin, TDayMax, GDtranspLow, WPveg, RatedHIdt,&
            CO2TnxReferenceYear, RefCropDay1, CropDeterm, CropSResp, TheCropType,&
            TheModeCycle, b0, b1, b2, &
            BM10, BM20, BM30, BM40, BM50, BM60, BM70)
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: LFlor
    integer(int32), intent(in) :: LengthFlor
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    integer(int32), intent(in) :: WPyield
    integer(int32), intent(in) :: RefHI
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: KcTop
    real(sp), intent(in) :: KcDeclAgeing
    real(sp), intent(in) :: CCeffectProcent
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax
    real(sp), intent(in) :: GDtranspLow
    real(sp), intent(in) :: WPveg
    real(sp), intent(in) :: RatedHIdt
    real(sp), intent(in) :: CO2TnxReferenceYear
    integer(int32), intent(in) :: RefCropDay1
    logical, intent(in) :: CropDeterm
    type(rep_Shapes), intent(in) :: CropSResp
    integer(intEnum), intent(in) :: TheCropType
    integer(intEnum), intent(in) :: TheModeCycle
    real(sp), intent(inout) :: b0
    real(sp), intent(inout) :: b1
    real(sp), intent(inout) :: b2
    real(sp), intent(inout) :: BM10
    real(sp), intent(inout) :: BM20
    real(sp), intent(inout) :: BM30
    real(sp), intent(inout) :: BM40
    real(sp), intent(inout) :: BM50
    real(sp), intent(inout) :: BM60
    real(sp), intent(inout) :: BM70

    type StressIndexes
        integer(int32) :: StressProc
            !! Undocumented
        real(sp) :: BioMProc
            !! Undocumented
        real(sp) :: BioMSquare
            !! Undocumented
    end type StressIndexes

    type(StressIndexes), dimension(8) :: StressMatrix
    integer(int8) :: Si
    integer(int32) :: L12SF, GDDL12SF
    type(rep_EffectStress) :: StressResponse
    real(sp) :: RatDGDD, BNor, BNor100, Yavg, X1avg, X2avg,&
                y, x1, x2, x1y, x2y, x1Sq, x2Sq, x1x2, &
                SUMx1y, SUMx2y, SUMx1Sq, SUMx2Sq, SUMx1x2
    integer(int32) :: SiPr
    real(sp) :: SumKcTop, HIGC, HIGClinear
    integer(int32) :: DaysYieldFormation, tSwitch
    real(sp) :: TDayMax_temp, TDayMin_temp

    ! 1. initialize
    call SetSimulation_DelayedDays(0) ! required for CalculateETpot
    L12SF = L12 ! to calculate SumKcTop (no stress)
    GDDL12SF = GDDL12 ! to calculate SumKcTop (no stress)
    ! Maximum sum Kc (no stress)
    SumKcTop = SeasonalSumOfKcPot(TheDaysToCCini, TheGDDaysToCCini,&
        L0, L12, L123, L1234, GDDL0, GDDL12, GDDL123, GDDL1234,&
        CCo, CCx, CGC, GDDCGC, CDC, GDDCDC, KcTop, KcDeclAgeing,&
        CCeffectProcent, Tbase, Tupper, TDayMin, TDayMax, &
        GDtranspLow, CO2TnxReferenceYear, TheModeCycle, .true.)

    ! Get PercentLagPhase (for estimate WPi during yield formation)
    if ((TheCropType == subkind_Tuber) .or. (TheCropType == subkind_grain)) then
        ! DaysToFlowering corresponds with Tuberformation
        DaysYieldFormation = roundc(RefHI/RatedHIdt, mold=1)
        if (CropDeterm) then
            HIGC = HarvestIndexGrowthCoefficient(real(RefHI, kind=sp), RatedHIdt)
            call GetDaySwitchToLinear(RefHI, RatedHIdt, HIGC, tSwitch,&
                  HIGClinear)
        else
            tSwitch = roundc(DaysYieldFormation/3._sp, mold=1)
        end if
    end if

    ! 2. Biomass production for various stress levels
    do Si = 1, 8
        ! various stress levels
        ! stress effect
        SiPr = int(10*(Si-1), kind=int32)
        StressMatrix(Si)%StressProc = SiPr
        call CropStressParametersSoilFertility(CropSResp, SiPr, StressResponse)
        ! adjusted length of Max canopy cover
        RatDGDD = 1
        if ((StressResponse%RedCCX == 0) .and. &
            (StressResponse%RedCGC == 0))then
            L12SF = L12
            GDDL12SF = GDDL12
        else
            call TimeToMaxCanopySF(CCo, CGC, CCx, L0, L12, L123, LFlor,&
                   LengthFlor, CropDeterm, L12SF, StressResponse%RedCGC,&
                   StressResponse%RedCCX, SiPr)
            if (TheModeCycle == modeCycle_GDDays) then
                TDayMin_temp = TDayMin
                TDayMax_temp = TDayMax
                GDDL12SF = SumCalendarDaysReferenceTnx(L12SF, RefCropDay1, RefCropDay1, Tbase, Tupper,&
                                 TDayMin_temp, TDayMax_temp)
            end if
            if ((TheModeCycle == modeCycle_GDDays) .and. (GDDL12SF < GDDL123)) then
                RatDGDD = (L123-L12SF)*1._sp/(GDDL123-GDDL12SF)
            end if
        end if
        ! biomass production
        BNor = Bnormalized(TheDaysToCCini, TheGDDaysToCCini,&
                L0, L12, L12SF, L123, L1234, LFlor,&
                GDDL0, GDDL12, GDDL12SF, GDDL123, GDDL1234, WPyield, &
                DaysYieldFormation, tSwitch, CCo, CCx, CGC, GDDCGC, CDC,&
                GDDCDC, KcTop, KcDeclAgeing, CCeffectProcent, WPveg, CO2TnxReferenceYear,&
                Tbase, Tupper, TDayMin, TDayMax, GDtranspLow, RatDGDD,&
                SumKcTop, SiPr, StressResponse%RedCGC, StressResponse%RedCCX,&
                StressResponse%RedWP, StressResponse%RedKsSto, 0_int8, 0 ,&
                StressResponse%CDecline, -0.01_sp, TheModeCycle, .true.,&
                .true.)
        if (Si == 1) then
            BNor100 = BNor
            StressMatrix(1)%BioMProc = 100._sp
        else
            if (BNor100 > 0.00001_sp) then
                StressMatrix(Si)%BioMProc = 100._sp * BNor/BNor100
            else
                StressMatrix(Si)%BioMProc = 100._sp
            end if
        end if
        StressMatrix(Si)%BioMSquare =&
             StressMatrix(Si)%BioMProc *&
             StressMatrix(Si)%BioMProc
        ! end stress level
    end do

    ! 5. Stress - Biomass relationship
    Yavg = 0._sp
    X1avg = 0._sp
    X2avg = 0._sp
    do Si = 1, 8
        ! various stress levels
        Yavg = Yavg + StressMatrix(Si)%StressProc
        X1avg = X1avg + StressMatrix(Si)%BioMProc
        X2avg = X2avg + StressMatrix(Si)%BioMSquare
    end do
    Yavg  = Yavg/8._sp
    X1avg = X1avg/8._sp
    X2avg = X2avg/8._sp
    SUMx1y  = 0._sp
    SUMx2y  = 0._sp
    SUMx1Sq = 0._sp
    SUMx2Sq = 0._sp
    SUMx1x2 = 0._sp
    do Si = 1, 8
        ! various stress levels
        y     = StressMatrix(Si)%StressProc - Yavg
        x1    = StressMatrix(Si)%BioMProc - X1avg
        x2    = StressMatrix(Si)%BioMSquare - X2avg
        x1y   = x1 * y
        x2y   = x2 * y
        x1Sq  = x1 * x1
        x2Sq  = x2 * x2
        x1x2  = x1 * x2
        SUMx1y  = SUMx1y + x1y
        SUMx2y  = SUMx2y + x2y
        SUMx1Sq = SUMx1Sq + x1Sq
        SUMx2Sq = SUMx2Sq + x2Sq
        SUMx1x2 = SUMx1x2 + x1x2
    end do

    if (abs(roundc(SUMx1x2*1000._sp, mold=1)) /= 0) then
        b2 = (SUMx1y - (SUMx2y * SUMx1Sq)/SUMx1x2)/&
             (SUMx1x2 - (SUMx1Sq * SUMx2Sq)/SUMx1x2)
        b1 = (SUMx1y - b2 * SUMx1x2)/SUMx1Sq
        b0 = Yavg - b1*X1avg - b2*X2avg

        BM10 =  StressMatrix(2)%BioMProc
        BM20 =  StressMatrix(3)%BioMProc
        BM30 =  StressMatrix(4)%BioMProc
        BM40 =  StressMatrix(5)%BioMProc
        BM50 =  StressMatrix(6)%BioMProc
        BM60 =  StressMatrix(7)%BioMProc
        BM70 =  StressMatrix(8)%BioMProc
    else
        b2 = real(undef_int, kind=sp)
        b1 = real(undef_int, kind=sp)
        b0 = real(undef_int, kind=sp)
    end if
end subroutine StressBiomassRelationshipForTnxReference


subroutine CCxSaltStressRelationshipForTnxReference(TheDaysToCCini, TheGDDaysToCCini,&
       L0, L12, L123, L1234, LFlor, LengthFlor, GDDFlor, GDDLengthFlor,&
       GDDL0, GDDL12, GDDL123, GDDL1234, WPyield, RefHI, CCo, CCx, CGC,&
       GDDCGC, CDC, GDDCDC, KcTop, KcDeclAgeing, CCeffectProcent, Tbase,&
       Tupper, TDayMin, TDayMax, GDbioLow, WPveg, RatedHIdt, CO2TnxReferenceYear,&
       CropDNr1, CropDeterm, TheCropType, TheModeCycle, TheCCsaltDistortion,&
       Coeffb0Salt, Coeffb1Salt, Coeffb2Salt, Salt10, Salt20, Salt30,&
       Salt40, Salt50, Salt60, Salt70, Salt80, Salt90)
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: LFlor
    integer(int32), intent(in) :: LengthFlor
    integer(int32), intent(in) :: GDDFlor
    integer(int32), intent(in) :: GDDLengthFlor
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    integer(int32), intent(in) :: WPyield
    integer(int32), intent(in) :: RefHI
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: KcTop
    real(sp), intent(in) :: KcDeclAgeing
    real(sp), intent(in) :: CCeffectProcent
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax
    real(sp), intent(in) :: GDbioLow
    real(sp), intent(in) :: WPveg
    real(sp), intent(in) :: RatedHIdt
    real(sp), intent(in) :: CO2TnxReferenceYear
    integer(int32), intent(in) :: CropDNr1
    logical, intent(in) :: CropDeterm
    integer(intEnum), intent(in) :: TheCropType
    integer(intEnum), intent(in) :: TheModeCycle
    integer(int8), intent(in) :: TheCCsaltDistortion
    real(sp), intent(inout) :: Coeffb0Salt
    real(sp), intent(inout) :: Coeffb1Salt
    real(sp), intent(inout) :: Coeffb2Salt
    real(sp), intent(inout) :: Salt10
    real(sp), intent(inout) :: Salt20
    real(sp), intent(inout) :: Salt30
    real(sp), intent(inout) :: Salt40
    real(sp), intent(inout) :: Salt50
    real(sp), intent(inout) :: Salt60
    real(sp), intent(inout) :: Salt70
    real(sp), intent(inout) :: Salt80
    real(sp), intent(inout) :: Salt90

    type StressIndexes
        integer(int32) :: CCxReduction
            !! Undocumented
        real(sp) :: SaltProc
            !! Undocumented
        real(sp) :: SaltSquare
            !! Undocumented
    end type StressIndexes

    integer(int32) :: L12SS, GDDL12SS, DaysYieldFormation, tSwitch
    real(sp) :: SumKcTop, HIGC, HIGClinear, CCToReach
    integer(int32) :: Si, SiPr
    type(StressIndexes), dimension(10) :: StressMatrix
    type(rep_EffectStress) :: StressResponse
    real(sp) :: RatDGDD, BNor, BNor100, BioMProc
    real(sp) :: Yavg, X1avg, X2avg, SUMx1y, SUMx2y, SUMx1Sq, &
         SUMx2Sq, SUMx1x2, y, x1, x2, x1y, x2y, x1Sq, x2Sq, x1x2
    real(sp) :: TDayMax_temp, TDayMin_temp

    ! 1. initialize
    call SetSimulation_DelayedDays(0) ! required for CalculateETpot
    GDDL12SS = GDDL12 ! to calculate SumKcTop (no stress)
    BNor100 = real(undef_int, kind=sp)
    ! Maximum sum Kc (no stress)
    SumKcTop = SeasonalSumOfKcPot(TheDaysToCCini, TheGDDaysToCCini,&
        L0, L12, L123, L1234, GDDL0, GDDL12, GDDL123, GDDL1234,&
        CCo, CCx, CGC, GDDCGC, CDC, GDDCDC, KcTop, KcDeclAgeing, &
        CCeffectProcent,Tbase, Tupper, TDayMin, TDayMax, GDbioLow, &
        CO2TnxReferenceYear, TheModeCycle, .true.)
    ! Get PercentLagPhase (for estimate WPi during yield formation)
    if ((TheCropType == subkind_Tuber) .or. (TheCropType == subkind_grain)) then
        ! DaysToFlowering corresponds with Tuberformation
        DaysYieldFormation = roundc(RefHI/RatedHIdt, mold=1)
        if (CropDeterm) then
            HIGC = HarvestIndexGrowthCoefficient(real(RefHI, kind=sp), RatedHIdt)
            call GetDaySwitchToLinear(RefHI, RatedHIdt, &
                    HIGC, tSwitch, HIGClinear)
        else
            tSwitch = roundc(DaysYieldFormation/3._sp, mold=1)
        end if
    end if

    ! 2. Biomass production (or Salt stress) for various CCx reductions
    do Si = 1, 10
        ! various CCx reduction
        ! CCx reduction
        SiPr = int(10*(Si-1), kind=int32)
        StressMatrix(Si)%CCxReduction = int(SiPr, kind=int8)
        ! adjustment CC
        call CropStressParametersSoilSalinity(int(SiPr, kind=int8), TheCCsaltDistortion, &
            CCo, CCx, CGC, GDDCGC, CropDeterm, L12, LFlor, LengthFlor, L123,&
            GDDL12, GDDFlor, GDDLengthFlor, GDDL123, TheModeCycle,&
            StressResponse)
        ! adjusted length of Max canopy cover
        RatDGDD = 1
        if ((StressResponse%RedCCX == 0) .and.&
            (StressResponse%RedCGC == 0)) then
            L12SS = L12
            GDDL12SS = GDDL12
        else
            CCToReach = 0.98_sp*(1._sp-StressResponse%RedCCX/100._sp)*CCx
            L12SS = DaysToReachCCwithGivenCGC(CCToReach, CCo, &
                 (1._sp-StressResponse%RedCCX/100._sp)*CCx,&
                 CGC*(1._sp-StressResponse%RedCGC/100._sp), L0)
            if (TheModeCycle == modeCycle_GDDays) then
                TDayMax_temp = TDayMax
                TDayMin_temp = TDayMin
                GDDL12SS = GrowingDegreeDays(L12SS, CropDNr1, Tbase, &
                           Tupper, TDayMin_temp, TDayMax_temp)
            end if
            if ((TheModeCycle == modeCycle_GDDays) .and.&
                (GDDL12SS < GDDL123)) then
                RatDGDD = (L123-L12SS)*1._sp/(GDDL123-GDDL12SS)
            end if
        end if

        ! biomass production
        BNor = Bnormalized(TheDaysToCCini, TheGDDaysToCCini,&
                L0, L12, L12SS, L123, L1234, LFlor,&
                GDDL0, GDDL12, GDDL12SS, GDDL123, GDDL1234,&
                WPyield, DaysYieldFormation, tSwitch,&
                CCo, CCx, CGC, GDDCGC, CDC, GDDCDC,&
                KcTop, KcDeclAgeing, CCeffectProcent, WPveg, CO2TnxReferenceYear,&
                Tbase, Tupper, TDayMin, TDayMax, GDbioLow, RatDGDD, SumKcTop,&
                SiPr, StressResponse%RedCGC, StressResponse%RedCCX,&
                StressResponse%RedWP, StressResponse%RedKsSto, &
                0_int8, 0, StressResponse%CDecline, -0.01_sp,&
                TheModeCycle, .false., .true.)
        if (Si == 1) then
            BNor100 = BNor
            BioMProc = 100._sp
            StressMatrix(1)%SaltProc = 0._sp
        else
            if (BNor100 > 0.00001_sp) then
                BioMProc = 100._sp * BNor/BNor100
                StressMatrix(Si)%SaltProc = 100._sp - BioMProc
            else
                StressMatrix(Si)%SaltProc = 0._sp
            end if
        end if
        StressMatrix(Si)%SaltSquare = &
             StressMatrix(Si)%SaltProc *&
             StressMatrix(Si)%SaltProc
        ! end stress level
    end do

    ! 3. CCx - Salt stress relationship
    Yavg = 0._sp
    X1avg = 0._sp
    X2avg = 0._sp
    do Si = 1, 10
        ! various CCx reduction
        Yavg = Yavg + StressMatrix(Si)%CCxReduction
        X1avg = X1avg + StressMatrix(Si)%SaltProc
        X2avg = X2avg + StressMatrix(Si)%SaltSquare
    end do
    Yavg  = Yavg/10._sp
    X1avg = X1avg/10._sp
    X2avg = X2avg/10._sp
    SUMx1y  = 0._sp
    SUMx2y  = 0._sp
    SUMx1Sq = 0._sp
    SUMx2Sq = 0._sp
    SUMx1x2 = 0._sp
    do Si = 1, 10
        ! various CCx reduction
        y     = StressMatrix(Si)%CCxReduction - Yavg
        x1    = StressMatrix(Si)%SaltProc - X1avg
        x2    = StressMatrix(Si)%SaltSquare - X2avg
        x1y   = x1 * y
        x2y   = x2 * y
        x1Sq  = x1 * x1
        x2Sq  = x2 * x2
        x1x2  = x1 * x2
        SUMx1y  = SUMx1y + x1y
        SUMx2y  = SUMx2y + x2y
        SUMx1Sq = SUMx1Sq + x1Sq
        SUMx2Sq = SUMx2Sq + x2Sq
        SUMx1x2 = SUMx1x2 + x1x2
    end do

    if (abs(roundc(SUMx1x2*1000._sp, mold=1)) /= 0) then
        Coeffb2Salt = (SUMx1y - (SUMx2y * SUMx1Sq)/SUMx1x2)/&
                      (SUMx1x2 - (SUMx1Sq * SUMx2Sq)/SUMx1x2)
        Coeffb1Salt = (SUMx1y - Coeffb2Salt * SUMx1x2)/SUMx1Sq
        Coeffb0Salt = Yavg - Coeffb1Salt*X1avg - Coeffb2Salt*X2avg

        Salt10 =  StressMatrix(2)%SaltProc
        Salt20 =  StressMatrix(3)%SaltProc
        Salt30 =  StressMatrix(4)%SaltProc
        Salt40 =  StressMatrix(5)%SaltProc
        Salt50 =  StressMatrix(5)%SaltProc
        Salt60 =  StressMatrix(7)%SaltProc
        Salt70 =  StressMatrix(8)%SaltProc
        Salt80 =  StressMatrix(9)%SaltProc
        Salt90 =  StressMatrix(10)%SaltProc
    else
        Coeffb2Salt = real(undef_int, kind=sp)
        Coeffb1Salt = real(undef_int, kind=sp)
        Coeffb0Salt = real(undef_int, kind=sp)
    end if
end subroutine CCxSaltStressRelationshipForTnxReference


subroutine ReferenceStressBiomassRelationship(TheDaysToCCini, &
        TheGDDaysToCCini, &
        L0, L12, L123, L1234, LFlor, LengthFlor, GDDL0, GDDL12, &
        GDDL123, GDDL1234, WPyield, RefHI, CCo, CCx, CGC, GDDCGC, &
        CDC, GDDCDC, KcTop, KcDeclAgeing, CCeffectProcent, &
        Tbase, Tupper, TDayMin, TDayMax, GDtranspLow, WPveg, RatedHIdt, &
        CropDNr1, CropDeterm, CropSResp, TheCropType, &
        TheModeCycle, b0, b1, b2, &
        BM10, BM20, BM30, BM40, BM50, BM60, BM70, &
        GDDFlor, GDDLengthFlor, GDDHImax, ThePlanting, LHImax)
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: LFlor
    integer(int32), intent(in) :: LengthFlor
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    integer(int32), intent(in) :: WPyield
    integer(int32), intent(in) :: RefHI
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: KcTop
    real(sp), intent(in) :: KcDeclAgeing
    real(sp), intent(in) :: CCeffectProcent
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax
    real(sp), intent(in) :: GDtranspLow
    real(sp), intent(in) :: WPveg
    real(sp), intent(in) :: RatedHIdt
    integer(int32), intent(in) :: CropDNr1
    logical, intent(in) :: CropDeterm
    type(rep_Shapes), intent(in) :: CropSResp
    integer(intEnum), intent(in) :: TheCropType
    integer(intEnum), intent(in) :: TheModeCycle
    real(sp), intent(inout) :: b0
    real(sp), intent(inout) :: b1
    real(sp), intent(inout) :: b2
    real(sp), intent(inout) :: BM10
    real(sp), intent(inout) :: BM20
    real(sp), intent(inout) :: BM30
    real(sp), intent(inout) :: BM40
    real(sp), intent(inout) :: BM50
    real(sp), intent(inout) :: BM60
    real(sp), intent(inout) :: BM70
    integer(int32), intent(in) :: GDDFlor
    integer(int32), intent(in) :: GDDLengthFlor
    integer(int32), intent(in) :: GDDHImax
    integer(intEnum), intent(in) :: ThePlanting
    integer(int32), intent(in) :: LHImax

    integer(int32) :: RefCropDay1
    integer(int32) :: Dayi, Monthi, Yeari
    real(sp) :: CO2TnxReferenceYear
    integer(int32) :: L0_loc
    integer(int32) :: L12_loc
    integer(int32) :: L123_loc
    integer(int32) :: L1234_loc
    integer(int32) :: LFlor_loc
    integer(int32) :: LengthFlor_loc
    integer(int32) :: LHImax_loc
    real(sp) :: CGC_loc
    real(sp) :: CDC_loc
    real(sp) :: RatedHIdt_loc

    L0_loc = L0
    L12_loc = L12
    LFlor_loc = LFlor
    LengthFlor_loc = LengthFlor
    L123_loc = L123
    L1234_loc = L1234
    LHImax_loc = LHImax
    CGC_loc = CGC
    CDC_loc =  CDC 
    RatedHIdt_loc = RatedHIdt

    ! 1. Day 1 of the GrowingCycle
    call DetermineDate(CropDNr1, Dayi, Monthi, Yeari)
    call DetermineDayNr(Dayi, Monthi, (1901), RefCropDay1)  ! not linked to a specific year

    ! 2. Create TCropReference.SIM (i.e. daily mean Tnx for 365 days from Onset onwards)
    if (GetTnxReferenceFile() /= '(None)') then
        call DailyTnxReferenceFileCoveringCropPeriod(RefCropDay1)
    end if

    ! 3. Determine coresponding calendar days if crop cycle is defined in GDDays
    if (TheModeCycle == modeCycle_GDDays) then
        call AdjustCalendarDaysReferenceTnx(RefCropDay1,&
        TheCropType,&
        Tbase, Tupper, TDayMin, TDayMax,&
        GDDL0, GDDL12, GDDFlor, GDDLengthFlor, GDDL123, GDDL1234,&
        GDDHImax,&
        GDDCGC, GDDCDC, CCo, CCx,&
        RefHI,&
        TheDaysToCCini, TheGDDaysToCCini,&
        ThePlanting,&
        L0_loc, L12_loc, LFlor_loc, LengthFlor_loc, L123_loc, &
        L1234_loc, LHImax_loc,&
        CGC_loc, CDC_loc, RatedHIdt_loc)
    end if

    ! 4. CO2 concentration for TnxReferenceYear
    if (GetTnxReferenceYear() == 2000) then
        CO2TnxReferenceYear = CO2Ref
    else
        CO2TnxReferenceYear = CO2ForTnxReferenceYear(GetTnxReferenceYear())
    end if

    ! 5. Stress Biomass relationship
    call StressBiomassRelationshipForTnxReference(TheDaysToCCini, TheGDDaysToCCini,&
    L0_loc, L12_loc, L123_loc, L1234_loc,&
    LFlor_loc, LengthFlor_loc,&
    GDDL0, GDDL12, GDDL123, GDDL1234, WPyield, RefHI,&
    CCo, CCx, CGC_loc, GDDCGC, CDC_loc, GDDCDC,&
    KcTop, KcDeclAgeing, CCeffectProcent,&
    Tbase, Tupper, TDayMin, TDayMax, GDtranspLow,&
    WPveg, RatedHIdt_loc, CO2TnxReferenceYear,&
    RefCropDay1,&
    CropDeterm,&
    CropSResp,&
    TheCropType,&
    TheModeCycle,&
    b0, b1, b2,&
    BM10, BM20, BM30, BM40, BM50, BM60, BM70)
end subroutine ReferenceStressBiomassRelationship


subroutine ReferenceCCxSaltStressRelationship(TheDaysToCCini, &
        TheGDDaysToCCini, L0, L12, L123, L1234, LFlor, LengthFlor, &
        GDDFlor, GDDLengthFlor, GDDL0, GDDL12, GDDL123, GDDL1234, &
        WPyield, RefHI, CCo, CCx, CGC, GDDCGC, CDC, GDDCDC, KcTop, &
        KcDeclAgeing, CCeffectProcent, Tbase, Tupper, TDayMin, &
        TDayMax, GDbioLow, WPveg, RatedHIdt, CropDNr1, CropDeterm, &
        TheCropType, TheModeCycle, TheCCsaltDistortion, Coeffb0Salt, &
        Coeffb1Salt, Coeffb2Salt, Salt10, Salt20, Salt30, Salt40, &
        Salt50, Salt60, Salt70, Salt80, Salt90, GDDHImax, &
        ThePlanting, LHImax)
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: LFlor
    integer(int32), intent(in) :: LengthFlor
    integer(int32), intent(in) :: GDDFlor
    integer(int32), intent(in) :: GDDLengthFlor
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    integer(int32), intent(in) :: WPyield
    integer(int32), intent(in) :: RefHI
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: KcTop
    real(sp), intent(in) :: KcDeclAgeing
    real(sp), intent(in) :: CCeffectProcent
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax
    real(sp), intent(in) :: GDbioLow
    real(sp), intent(in) :: WPveg
    real(sp), intent(in) :: RatedHIdt
    integer(int32), intent(in) :: CropDNr1
    logical, intent(in) :: CropDeterm
    integer(intEnum), intent(in) :: TheCropType
    integer(intEnum), intent(in) :: TheModeCycle
    integer(int8), intent(in) :: TheCCsaltDistortion
    real(sp), intent(inout) :: Coeffb0Salt
    real(sp), intent(inout) :: Coeffb1Salt
    real(sp), intent(inout) :: Coeffb2Salt
    real(sp), intent(inout) :: Salt10
    real(sp), intent(inout) :: Salt20
    real(sp), intent(inout) :: Salt30
    real(sp), intent(inout) :: Salt40
    real(sp), intent(inout) :: Salt50
    real(sp), intent(inout) :: Salt60
    real(sp), intent(inout) :: Salt70
    real(sp), intent(inout) :: Salt80
    real(sp), intent(inout) :: Salt90
    integer(int32), intent(in) :: GDDHImax
    integer(intEnum), intent(in) :: ThePlanting
    integer(int32), intent(in) :: LHImax

    integer(int32) :: RefCropDay1
    integer(int32) :: Dayi, Monthi, Yeari
    real(sp) :: CO2TnxReferenceYear
    integer(int32) :: L0_loc
    integer(int32) :: L12_loc
    integer(int32) :: L123_loc
    integer(int32) :: L1234_loc
    integer(int32) :: LFlor_loc
    integer(int32) :: LengthFlor_loc
    integer(int32) :: LHImax_loc
    real(sp) :: CGC_loc
    real(sp) :: CDC_loc
    real(sp) :: RatedHIdt_loc

    L0_loc = L0
    L12_loc = L12
    L123_loc = L123
    L1234_loc = L1234
    LFlor_loc = LFlor
    LengthFlor_loc = LengthFlor
    LHImax_loc = LHImax
    CGC_loc = CGC
    CDC_loc = CDC
    RatedHIdt_loc = RatedHIdt

    ! 1. Day 1 of the GrowingCycle
    call DetermineDate(CropDNr1, Dayi, Monthi, Yeari)
    call DetermineDayNr(Dayi, Monthi, (1901), RefCropDay1)  ! not linked to a specific year

    ! 2. Create TCropReference.SIM (i.e. daily mean Tnx for 365 days from Onset onwards)
    if (GetTnxReferenceFile() /= '(None)') then
        call DailyTnxReferenceFileCoveringCropPeriod(RefCropDay1)
    end if

    ! 3. Determine coresponding calendar days if crop cycle is defined in GDDays
    if (TheModeCycle == modeCycle_GDDays) then
        call AdjustCalendarDaysReferenceTnx(RefCropDay1,&
        TheCropType,&
        Tbase, Tupper, TDayMin, TDayMax,&
        GDDL0, GDDL12, GDDFlor, GDDLengthFlor, GDDL123, GDDL1234,&
        GDDHImax,&
        GDDCGC, GDDCDC, CCo, CCx,&
        RefHI,&
        TheDaysToCCini, TheGDDaysToCCini,&
        ThePlanting,&
        L0_loc, L12_loc, LFlor_loc, LengthFlor_loc, L123_loc, &
        L1234_loc, LHImax_loc,&
        CGC_loc, CDC_loc, RatedHIdt_loc)
    end if

    ! 4. CO2 concentration for TnxReferenceYear
    if (GetTnxReferenceYear() == 2000) then
        CO2TnxReferenceYear = CO2Ref
    else
        CO2TnxReferenceYear = CO2ForTnxReferenceYear(GetTnxReferenceYear())
    end if

    ! 5. Stress Biomass relationship for salinity
    call CCxSaltStressRelationshipForTnxReference(TheDaysToCCini, TheGDDaysToCCini,&
    L0_loc, L12_loc, L123_loc, L1234_loc,&
    LFlor_loc, LengthFlor_loc, GDDFlor, GDDLengthFlor,&
    GDDL0, GDDL12, GDDL123, GDDL1234, WPyield, RefHI,&
    CCo, CCx, CGC_loc, GDDCGC, CDC_loc, GDDCDC,&
    KcTop, KcDeclAgeing, CCeffectProcent,&
    Tbase, Tupper, TDayMin, TDayMax, GDbioLow, WPveg, RatedHIdt_loc, CO2TnxReferenceYear,&
    CropDNr1,&
    CropDeterm,&
    TheCropType,&
    TheModeCycle,&
    TheCCsaltDistortion,&
    Coeffb0Salt, Coeffb1Salt, Coeffb2Salt,&
    Salt10, Salt20, Salt30, Salt40, Salt50, Salt60, Salt70, Salt80, Salt90)
end subroutine ReferenceCCxSaltStressRelationship

end module ac_preparefertilitysalinity
