module ac_tempprocessing

use ac_global , only: undef_int, &
                      GetTmaxTnxReference12MonthsRun_i,&
                      GetTminTnxReference12MonthsRun_i,&
                      GetTmaxTnxReference12MonthsRun,&
                      GetTminTnxReference12MonthsRun,&
                      SetTmaxTnxReference12MonthsRun_i,&
                      SetTminTnxReference12MonthsRun_i,&
                      SetTmaxTnxReference12MonthsRun,&
                      SetTminTnxReference12MonthsRun,&
                      fTnxReference,&
                      fTnxReference_iostat,&
                      fTnxReference365Days,& 
                      fTnxReference365Days_iostat,&
                      GetTnxReference365DaysFileFull,&
                      modeCycle_GDDays, &
                      modeCycle_CalendarDays, &
                      DaysinMonth, &
                      rep_DayEventDbl, &
                      rep_CropFileSet, &
                      rep_EffectStress, &
                      rep_Shapes, &
                      rep_clim, &
                      max_No_compartments, &
                      CompartmentIndividual, &
                      subkind_Vegetative, &
                      subkind_Grain, &
                      subkind_Tuber, &
                      subkind_Forage, &
                      datatype_daily, &
                      datatype_decadely, &
                      datatype_monthly, &
                      CalculateETpot, &
                      CCiNoWaterStressSF, &
                      CanopyCoverNoStressSF, &
                      DetermineDayNr, &
                      DetermineDate, &
                      ECeComp, &
                      FileExists, &
                      GetDaySwitchToLinear, &
                      HarvestIndexGrowthCoefficient, &
                      LeapYear, &
                      SeasonalSumOfKcPot, &
                      SplitStringInTwoParams, &
                      KsTemperature, &
                      DegreesDay, &
                      LengthCanopyDecline, &
                      DetermineLengthGrowthStages, &
                      AdjustSizeCompartments, &
                      AdjustOnsetSearchPeriod, &
                      adjustcropyeartoclimfile, &
                      GenerateCO2Description, &
                      LoadProfile, &
                      LoadProfileProcessing, &
                      LoadClim, &
                      LoadIrriScheduleInfo,&
                      LoadManagement, &
                      GetECiAqua,&
                      SetECiAqua,&
                      GetECStorage,&
                      GetSurfaceStorage,&
                      SetSurfaceStorage, &
                      GetNrCompartments, &
                      GetManagement_FertilityStress, &
                      CanopyCoverNoStressSF, &
                      CCmultiplierWeed, &
                      CCiNoWaterStressSF, &
                      CalculateETpot, &
                      CalculateAdjustedFC, &
                      GetPathNameSimul, &
                      GetTemperatureFile, &
                      GetTemperatureFilefull, &
                      GetTemperatureDescription, &
                      GetEToFile, GetEToFilefull, &
                      GetEToDescription,&
                      GetEToRecord, &
                      GetRainFile, GetRainFileFull,&
                      Getraindescription,&
                      GetCO2File, GetCO2FileFull,&
                      GetCO2Description,&
                      GetClimFile,&
                      GetIrriFile, GetIrriFileFull,&
                      GetManFile, GetManFilefull,&
                      GetProfFile, GetProfFilefull,&
                      GetOffSeasonFile, GetOffSeasonFilefull,&
                      GetGroundWaterFile, GetGroundWaterFilefull,&
                      GetObservationsFile,&
                      GetObservationsDescription,&
                      GetObservationsFileFull,&
                      GetCompartment_Thickness, &
                      GetCompartment_Theta, &
                      GetCompartment_i, &
                      GetZiAqua, &
                      GetManagement_BundHeight, &
                      GetRainRecord, &
                      GetClimRecord_NrObs, &
                      GetClimRecord_FromY, &
                      GetTemperatureRecord, &
                      GetTemperatureRecord_FromD, &
                      GetTemperatureRecord_FromM, &
                      GetTemperatureRecord_FromY, &
                      GetTemperatureRecord_NrObs, &
                      GetTemperatureRecord_ToD, &
                      GetTemperatureRecord_ToM, &
                      GetTemperatureRecord_ToY, &
                      GetTemperatureRecord_DataType, &
                      GetTemperatureRecord_FromDayNr, &
                      GetTemperatureRecord_ToDayNr, &
                      GetSimulParam_GDDMethod, &
                      GetCalendarFile, GetCalendarFileFull,&
                      GetCalendarDescription, &
                      GetCropFilefull,&
                      GetCrop_RootMin,&
                      GetCrop_SizePlant,&
                      SetRainFile,SetManFile,&
                      SetIrriFileFull,&
                      FullUndefinedRecord, &
                      SetEToFile, SetEToFilefull, &
                      SetObservationsFileFull, &
                      SetProfFilefull,&
                      SetGroundWaterFilefull,&
                      SetTemperatureFileFull,&
                      SetCO2FileFull,&
                      SetSimulation_ECStorageIni,&
                      SetCrop_Planting,&
                      GetCrop_SownYear1, &
                      SetCrop_DayN, &
                      SetCrop_GDDaysToFullCanopy, &
                      SetCrop_GDDaysToHIo, &
                      SetCrop_DaysToHarvest, &
                      SetCrop_DaysToSenescence, &
                      SetCrop_CDC, &
                      SetCrop_DaysToHIo, &
                      SetCrop_DaysToGermination, &
                      SetCrop_LengthFlowering, &
                      SetCrop_Length, &
                      SetCrop_DaysToFullCanopy, &
                      SetCrop_dHIdt, &
                      SetCrop_DaysToMaxRooting, &
                      SetCrop_CGC, &
                      SetCrop_DaysToFlowering, &
                      SetCrop_GDDaysToHarvest, &
                      SetSWCIniFile, SetSWCIniFileFull,&
                      SetTemperatureDescription,&
                      SetZiAqua, &
                      SetGroundWaterDescription, &
                      SetClimateFileFull, &
                      SetSimulation_LinkCropToSimPeriod, &
                      SetRainFilefull, &
                      SetCalendarFileFull, &
                      SetManFileFull, &
                      GetCropFileSet, &
                      GetSimulation_EffectStress_RedCGC, &
                      GetSimulation_EffectStress_RedCCX, &
                      GetCrop_PlantingDens, &
                      GetCrop_SizeSeedling, &
                      GetCrop_DaysToFullCanopySF, &
                      GetCrop_DayN,GetCrop_Day1, &
                      SetSimulation_DelayedDays, &
                      SetSimulation_DelayedDays, &
                      GetCrop_GDDaysToGermination, &
                      GetCrop_GDDaysToFullCanopy, &
                      GetCrop_GDDaysToHarvest,   GetCrop_GDDaysToHIo, &
                      GetCrop_GDDaysToFlowering, GetCrop_GDDLengthFlowering, &
                      GetCrop_GDDaysToSenescence,GetCrop_GDDaysToHarvest, &
                      GetCrop_GDDaysToMaxRooting, &
                      GetCrop_DaysToGermination, &
                      GetCrop_DaysToFullCanopy, &
                      GetCrop_DaysToFlowering, &
                      GetCrop_LengthFlowering, &
                      GetCrop_DaysToSenescence, &
                      GetCrop_DaysToHarvest, &
                      GetCrop_DaysToMaxRooting, &
                      GetCrop_DaysToHIo, &
                      GetCrop_Length, &
                      GetCrop_GDDCGC, GetCrop_GDDCDC, GetCrop_CCo, &
                      GetCrop_CCx, GetCrop_HI, &
                      GetCrop_CGC, &
                      GetCrop_CDC, &
                      GetCrop_dHIdt, &
                      GetCrop_DaysToCCini, GetCrop_Planting, &
                      GetCrop_subkind, &
                      GetCrop_ModeCycle, &
                      GetCrop_Tbase, &
                      GetCrop_Tupper, &
                      GetCrop_Planting, &
                      GetCrop_CCini, &
                      GetCrop_GDDaysToCCini, &
                      GetCrop_DeterminancyLinked, &
                      GetCrop_RootMax, &
                      GetSoil_RootMax, &
                      GetCrop_RootMinYear1, &
                      GetCompartment,&
                      GetSimulParam_Tmin, GetSimulParam_Tmax,&
                      GetSimulation_ResetIniSWC,&
                      GetSimulation_MultipleRunWithKeepSWC,&
                      GetSimulation_MultipleRunConstZrx,&
                      GetSimulation_ToDayNr, &
                      GetSimulation_FromDayNr,&
                      GetSimulation_YearSeason, &
                      GetSimulation_IniSWC_NrLoc, &
                      GetSimulation_IniSWC_AtDepths, &
                      GetSimulation_IniSWC_Loc, &
                      GetSimulation_IniSWC_VolProc, &
                      GetSimulation_IniSWC_SaltECe, &
                      GetSimulation_IniSWC_AtFC, &
                      GetSWCIniFile, GetSWCiniFileFull, &
                      GetClimateFile, GetClimateFileFull, &
                      GetSimulParam_Tmin, GetSimulParam_Tmax, &
                      SetSimulation_DelayedDays,&
                      GetManagement_WeedAdj, &
                      GetSimulParam_Tmin, GetSimulParam_Tmax,&
                      GetWeedRC, &
                      DaysToReachCCwithGivenCGC, &
                      timetomaxcanopysf, &
                      cropstressparameterssoilfertility,&
                      GetCropFile, &
                      setclimatedescription,&
                      setoffseasonfilefull,&
                      settemperaturefile,&
                      setirrifile,&
                      loadoffseason,&
                      setco2description,&
                      setco2file,&
                      setoffseasonfile,&
                      completecropdescription,&
                      completeprofiledescription,&
                      setcrop_rootmin,&
                      setetodescription,&
                      setcrop_gddaystoccini,&
                      loadinitialconditions,&
                      timetomaxcanopysf,&
                      setsimulation_fromdaynr,&
                      resetswctofc,&
                      setsimulparam_constgwt,&
                      setclimdata,&
                      setcrop_gddaystosenescence,&
                      setobservationsfile,&
                      setmanagement_fertilitystress,&
                      setcalendardescription,&
                      setcropfile,&
                      adjustyearperennials,&
                      completeclimatedescription,&
                      settemperaturerecord,&
                      setcrop_daystoccini,&
                      setetorecord,&
                      setcrop_day1,&
                      setcrop_sizeplant,&
                      adjustsimperiod,&
                      setswcinidescription,&
                      setsimulation_yearseason,&
                      setsimulation_effectstress_redcgc,&
                      setproffile,&
                      loadgroundwater,&
                      setcrop_daystofullcanopysf,&
                      setcropfilefull,&
                      setoffseasondescription,&
                      setcompartment,&
                      setmandescription,&
                      setraindescription,&
                      setsimulation_eceini_i,&
                      setclimatefile,&
                      setrainrecord,&
                      setobservationsdescription,&
                      adjustclimrecordto,&
                      loadcrop,&
                      setsimulation_todaynr,&
                      noirrigation,&
                      setsimulation_thetaini_i,&
                      translateinilayerstoswprofile,&
                      setsimulation_effectstress_redccx,&
                      setcalendarfile,&
                      setcrop_ccini,&
                      setgroundwaterfile,&
                      setsimulation_surfacestorageini,&
                      translateinipointstoswprofile,&
                      KsAny, &
                      SetTminRun_i, &
                      SetTminRun, &
                      SetTmaxRun_i, &
                      SetTmaxRun, &
                      GetTminRun_i, &
                      GetTminRun, &
                      GetTmaxRun_i, &
                      GetTmaxRun, &
                      SetTminCropReferenceRun_i, &
                      SetTminCropReferenceRun, &
                      SetTmaxCropReferenceRun_i, &
                      SetTmaxCropReferenceRun, &
                      GetTminCropReferenceRun_i, &
                      GetTminCropReferenceRun, &
                      GetTmaxCropReferenceRun_i, &
                      GetTmaxCropReferenceRun, &
                      SetTminTnxReference365DaysRun_i, &
                      SetTminTnxReference365DaysRun, &
                      SetTmaxTnxReference365DaysRun_i, &
                      SetTmaxTnxReference365DaysRun, &
                      GetTminTnxReference365DaysRun_i, &
                      GetTminTnxReference365DaysRun, &
                      GetTmaxTnxReference365DaysRun_i, &
                      GetTmaxTnxReference365DaysRun, &
                      GetTnxReferenceYear, &
                      SetTnxReferenceYear, &
                      GetTnxReferenceFile, &
                      SetTnxReferenceFile, &
                      GetTnxReferenceFileFull, &
                      SetTnxReferenceFileFull, &
                      TminTnxReference365DaysRun, &
                      TmaxTnxReference365DaysRun, &
                      TminCropReferenceRun, &
                      TmaxCropReferenceRun
use ac_kinds,  only: sp,&
                     dp, &
                     int8, &
                     int16, &
                     int32, &
                     intEnum
use ac_project_input, only: ProjectInput
use ac_utils, only: roundc, &
                    write_file, &
                    open_file
use iso_fortran_env, only: iostat_end
implicit none


logical :: TemperatureFilefull_exists
real(dp), dimension(:), allocatable :: Tmin, Tmax  !! (daily) temperature data
type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet


contains


subroutine AdjustMONTHandYEAR(MFile, Yfile)
    integer(int32), intent(inout) :: MFile
    integer(int32), intent(inout) :: Yfile

    Mfile = Mfile - 12
    YFile = Yfile + 1
end subroutine AdjustMONTHandYEAR


subroutine AdjustDecadeMONTHandYEAR(DecFile, Mfile, Yfile)
    integer(int32), intent(inout) :: DecFile
    integer(int32), intent(inout) :: Mfile
    integer(int32), intent(inout) :: Yfile

    DecFile = 1
    Mfile = Mfile + 1
    if (Mfile > 12) then
        Mfile = 1
        YFile = Yfile + 1
    end if
end subroutine AdjustDecadeMONTHandYEAR


subroutine ReadTemperatureFilefull()
    !! Reads the contents of the TemperatureFilefull file,
    !! storing the temperatures in the Tmin and Tmax arrays.

    character(len=:), allocatable :: filename
    integer :: fhandle, i, nlines, nrows, rc

    filename = GetTemperatureFilefull()
    open(newunit=fhandle, file=trim(filename), status='old', action='read')

    ! Count the number of lines
    nlines = 0
    do
        read(fhandle, '(a)', iostat=rc)
        if (rc == iostat_end) then
            exit
        else
            nlines = nlines + 1
        end if
    end do

    ! Now read in the actual content
    nrows = nlines - 8
    if (allocated(Tmin)) deallocate(Tmin)
    if (allocated(Tmax)) deallocate(Tmax)
    allocate(Tmin(nrows), Tmax(nrows))

    rewind(fhandle)
    read(fhandle, *) ! description
    read(fhandle, *) ! time step
    read(fhandle, *) ! day
    read(fhandle, *) ! month
    read(fhandle, *) ! year
    read(fhandle, *)
    read(fhandle, *)
    read(fhandle, *)
    do i = 1, nrows
        read(fhandle, *) Tmin(i), Tmax(i)
    end do

    close(fhandle)
end subroutine ReadTemperatureFilefull


subroutine SetDayNrToYundef(DayNri)
    integer(int32), intent(inout) :: DayNri

    integer(int32) :: Dayi, Monthi, Yeari

    call DetermineDate(DayNri, Dayi, Monthi, Yeari)
    Yeari = 1901
    call DetermineDayNr(Dayi, Monthi, Yeari, DayNri)
end subroutine SetDayNrToYundef


subroutine GetDecadeTemperatureDataSet(DayNri, TminDataSet, TmaxDataSet)
    integer(int32), intent(in) :: DayNri
    type(rep_DayEventDbl), dimension(31) , intent(inout) :: TminDataSet
    type(rep_DayEventDbl), dimension(31) , intent(inout) :: TmaxDataSet

    integer(int32) :: Nri, ni, Dayi, Deci, Monthi, Yeari, DayN
    integer(int32) :: DNR
    real(dp) :: C1Min, C1Max, C2Min, C2Max, C3Min, C3Max
    real(dp) :: UlMin, LLMin, MidMin, UlMax, LLMax, MidMax

    call DetermineDate(DayNri, Dayi, Monthi, Yeari)
    if (Dayi > 20) then
        Deci = 3
        Dayi = 21
        DayN = DaysInMonth(Monthi)
        if ((Monthi == 2) .and. LeapYear(Yeari)) then
            DayN = DayN + 1
        end if
        ni = DayN - Dayi + 1
    elseif (Dayi > 10) then
        Deci = 2
        Dayi = 11
        DayN = 20
        ni = 10
    else
        Deci = 1
        Dayi = 1
        DayN = 10
        ni = 10
    end if
    call GetSetofThree(DayN, Deci, Monthi, Yeari, &
               C1Min, C1Max, C2Min, C2Max, C3Min, C3Max)
    call DetermineDayNr(Dayi, Monthi, Yeari, DNR)

    call GetParameters(C1Min, C2Min, C3Min, ULMin, LLMin, MidMin)
    do Nri = 1, ni
        TMinDataSet(Nri)%DayNr = DNR+Nri-1
        if (Nri <= (ni/2._dp+0.01_dp)) then
            TMinDataSet(Nri)%Param = (2._dp*ULMin + &
                          (MidMin-ULMin)*(2._dp*Nri-1._dp)/(ni/2._dp))/2._dp
        else
            if (((ni == 11) .or. (ni == 9)) .and. (Nri < (ni+1.01_dp)/2._dp)) then
                TminDataSet(Nri)%Param = MidMin
            else
                TminDataSet(Nri)%Param = (2._dp*MidMin + &
                          (LLMin-MidMin)*(2._dp*Nri-(ni+1))/(ni/2._dp))/2._dp
            end if
        end if
    end do

    call GetParameters(C1Max, C2Max, C3Max, ULMax, LLMax, MidMax)
    do Nri = 1, ni
        TMaxDataSet(Nri)%DayNr = DNR+Nri-1
        if (Nri <= (ni/2._dp+0.01_dp)) then
            TMaxDataSet(Nri)%Param = (2._dp*ULMax + &
                          (MidMax-ULMax)*(2._dp*Nri-1)/(ni/2._dp))/2._dp
        else
            if (((ni == 11) .or. (ni == 9)) .and. (Nri < (ni+1.01_dp)/2._dp)) then
                 TmaxDataSet(Nri)%Param = MidMax
            else
                TmaxDataSet(Nri)%Param = (2._dp*MidMax + &
                          (LLMax-MidMax)*(2._dp*Nri-(ni+1))/(ni/2._dp))/2._dp
            end if
        end if
    end do

    do Nri = (ni+1), 31
        TminDataSet(Nri)%DayNr = DNR+ni-1
        TminDataSet(Nri)%Param = 0._dp
        TmaxDataSet(Nri)%DayNr = DNR+ni-1
        TmaxDataSet(Nri)%Param = 0._dp
    end do


    contains


    subroutine GetSetofThree(DayN, Deci, Monthi, Yeari, &
                     C1Min, C1Max, C2Min, C2Max, C3Min, C3Max)
        integer(int32), intent(in) :: DayN
        integer(int32), intent(in) :: Deci
        integer(int32), intent(in) :: Monthi
        integer(int32), intent(in) :: Yeari
        real(dp), intent(inout) :: C1Min
        real(dp), intent(inout) :: C1Max
        real(dp), intent(inout) :: C2Min
        real(dp), intent(inout) :: C2Max
        real(dp), intent(inout) :: C3Min
        real(dp), intent(inout) :: C3Max

        integer(int32) :: fhandle
        integer(int32) :: DecFile, Mfile, Yfile, Nri, Obsi, rc
        logical :: OK3
        character(len=255) :: StringREAD

        !! 1 = previous decade, 2 = Actual decade, 3 = Next decade;
        open(newunit=fhandle, file=trim(GetTemperatureFilefull()), &
                     status='old', action='read', iostat=rc)
        read(fhandle, *, iostat=rc) ! description
        read(fhandle, *, iostat=rc) ! time step
        read(fhandle, *, iostat=rc) ! day
        read(fhandle, *, iostat=rc) ! month
        read(fhandle, *, iostat=rc) ! year
        read(fhandle, *, iostat=rc)
        read(fhandle, *, iostat=rc)
        read(fhandle, *, iostat=rc)

        if (GetTemperatureRecord_FromD() > 20) then
            DecFile = 3
        elseif (GetTemperatureRecord_FromD() > 10) then
            DecFile = 2
        else
            DecFile = 1
        end if
        Mfile = GetTemperatureRecord_FromM()
        if (GetTemperatureRecord_FromY() == 1901) then
            Yfile = Yeari
        else
            Yfile = GetTemperatureRecord_FromY()
        end if
        OK3 = .false.

        if (GetTemperatureRecord_NrObs() <= 2) then
            read(fhandle, '(a)', iostat=rc) StringREAD
            call SplitStringInTwoParams(StringREAD, C1Min, C1Max)
            select case (GetTemperatureRecord_NrObs())
            case (0)
                C2Min = C1Min
                C2Max = C2Max
                C3Min = C1Min
                C3Max = C1Max
            case (1)
                DecFile = DecFile + 1
                if (DecFile > 3) then
                    call AdjustDecadeMONTHandYEAR(DecFile, Mfile, Yfile)
                end if
                read(fhandle, '(a)', iostat=rc) StringREAD
                call SplitStringInTwoParams(StringREAD, C3Min, C3Max)
                if (Deci == DecFile) then
                    C2Min = C3Min
                    C2Max = C3Max
                    C3Min = C2Min+(C2Min-C1Min)/4._dp
                    C3Max = C2Max+(C2Max-C1Max)/4._dp
                else
                    C2Min = C1Min
                    C2Max = C1Max
                    C1Min = C2Min + (C2Min-C3Min)/4._dp
                    C1Max = C2Max + (C2Max-C3Max)/4._dp
                end if
            end select
            OK3 = .true.
        end if

       if ((.not. OK3) .and. ((Deci == DecFile) .and. (Monthi == Mfile) &
            .and. (Yeari == Yfile))) then
            read(fhandle, '(a)', iostat=rc) StringREAD
            call SplitStringInTwoParams(StringREAD, C1Min, C1Max)
            C2Min = C1Min
            C2Max = C1Max
            read(fhandle, '(a)', iostat=rc) StringREAD
            call SplitStringInTwoParams(StringREAD, C3Min, C3Max)
            C1Min = C2Min + (C2Min-C3Min)/4._dp
            C1Max = C2Max + (C2Max-C3Max)/4._dp
            OK3 = .true.
        end if

        if ((.not. OK3) .and. ((DayN == GetTemperatureRecord_ToD()) &
             .and. (Monthi == GetTemperatureRecord_ToM()))) then
            if ((GetTemperatureRecord_FromY() == 1901) .or. &
                (Yeari == GetTemperatureRecord_ToY())) then
                do Nri = 1, (GetTemperatureRecord_NrObs()-2)
                     read(fhandle, *, iostat=rc)
                end do
                read(fhandle, '(a)', iostat=rc) StringREAD
                call SplitStringInTwoParams(StringREAD, C1Min, C1Max)
                read(fhandle, '(a)', iostat=rc) StringREAD
                call SplitStringInTwoParams(StringREAD, C2Min, C2Max)
                C3Min = C2Min+(C2Min-C1Min)/4._dp
                C3Max = C2Max+(C2Max-C1Max)/4._dp
                OK3 = .true.
            end if
        end if

        if (.not. OK3) then
            Obsi = 1
            do while (.not. OK3)
                if ((Deci == DecFile) .and. (Monthi == Mfile) &
                    .and. (Yeari == Yfile)) then
                    OK3 = .true.
                else
                    DecFile = DecFile + 1
                    if (DecFile > 3) then
                        call AdjustDecadeMONTHandYEAR(DecFile, Mfile, Yfile)
                    end if
                    Obsi = Obsi + 1
                end if
            end do
            if (GetTemperatureRecord_FromD() > 20) then
                DecFile = 3
            elseif (GetTemperatureRecord_FromD() > 10) then
                DecFile = 2
            else
                DecFile = 1
            end if
            do Nri = 1, (Obsi-2)
                read(fhandle, *, iostat=rc)
            end do
            read(fhandle, '(a)', iostat=rc) StringREAD
            call SplitStringInTwoParams(StringREAD, C1Min, C1Max)
            read(fhandle, '(a)', iostat=rc) StringREAD
            call SplitStringInTwoParams(StringREAD, C2Min, C2Max)
            read(fhandle, '(a)', iostat=rc) StringREAD
            call SplitStringInTwoParams(StringREAD, C3Min, C3Max)
        end if
        close(fhandle)
    end subroutine GetSetofThree


    subroutine GetParameters(C1, C2, C3, UL, LL, Mid)
        real(dp), intent(in) :: C1
        real(dp), intent(in) :: C2
        real(dp), intent(in) :: C3
        real(dp), intent(inout) :: UL
        real(dp), intent(inout) :: LL
        real(dp), intent(inout) :: Mid

        UL = (C1+C2)/2._dp
        LL = (C2+C3)/2._dp
        Mid = 2._dp*C2 - (UL+LL)/2._dp
        ! --previous decade-->/UL/....... Mid ......../LL/<--next decade--
    end subroutine GetParameters
end subroutine GetDecadeTemperatureDataSet


subroutine GetMonthlyTemperatureDataSet(DayNri, TminDataSet, TmaxDataSet)
    integer(int32), intent(in) :: DayNri
    type(rep_DayEventDbl), dimension(31) , intent(inout) :: TminDataSet
    type(rep_DayEventDbl), dimension(31) , intent(inout) :: TmaxDataSet

    integer(int32) :: Dayi, Monthi, Yeari, DayN
    integer(int32) :: DNR
    integer(int32) :: X1, X2, X3, t1, t2
    real(dp) :: C1Min, C2Min, C3Min
    real(dp) :: C1Max, C2Max, C3Max
    real(dp) :: aOver3Min, bOver2Min, cMin
    real(dp) :: aOver3Max, bOver2Max, cMax

    call DetermineDate(DayNri, Dayi, Monthi, Yeari)
    call GetSetofThreeMonths(Monthi, Yeari, &
         C1Min, C2Min, C3Min, C1Max, C2Max, C3Max, X1, X2, X3, t1)

    Dayi = 1
    call DetermineDayNr(Dayi, Monthi, Yeari, DNR)
    DayN = DaysInMonth(Monthi)
    if ((Monthi == 2) .and. LeapYear(Yeari)) then
        DayN = DayN + 1
    end if

    call GetInterpolationParameters(C1Min, C2Min, C3Min, &
                   aOver3Min, bOver2Min, cMin)
    call GetInterpolationParameters(C1Max, C2Max, C3Max, &
                   aOver3Max, bOver2Max, cMax)
    do Dayi = 1, DayN
        t2 = t1 + 1
        TminDataSet(Dayi)%DayNr = DNR+Dayi-1
        TmaxDataSet(Dayi)%DayNr = DNR+Dayi-1
        TminDataSet(Dayi)%Param = aOver3Min*(t2*t2*t2-t1*t1*t1) &
            + bOver2Min*(t2*t2-t1*t1) + cMin*(t2-t1)
        TmaxDataSet(Dayi)%Param = aOver3Max*(t2*t2*t2-t1*t1*t1) &
            + bOver2Max*(t2*t2-t1*t1) + cMax*(t2-t1)
        t1 = t2
    end do
    do Dayi = (DayN+1), 31
        TminDataSet(Dayi)%DayNr = DNR+DayN-1
        TmaxDataSet(Dayi)%DayNr = DNR+DayN-1
        TminDataSet(Dayi)%Param = 0._dp
        TmaxDataSet(Dayi)%Param = 0._dp
    end do


    contains


    subroutine GetSetofThreeMonths(Monthi, Yeari, &
            C1Min, C2Min, C3Min, C1Max, C2Max, C3Max, X1, X2, X3, t1)
        integer(int32), intent(in) :: Monthi
        integer(int32), intent(in) :: Yeari
        real(dp), intent(inout) :: C1Min
        real(dp), intent(inout) :: C2Min
        real(dp), intent(inout) :: C3Min
        real(dp), intent(inout) :: C1Max
        real(dp), intent(inout) :: C2Max
        real(dp), intent(inout) :: C3Max
        integer(int32), intent(inout) :: X1
        integer(int32), intent(inout) :: X2
        integer(int32), intent(inout) :: X3
        integer(int32), intent(inout) :: t1

        integer(int32), parameter :: n1 = 30
        integer(int32), parameter :: n2 = 30
        integer(int32), parameter :: n3 = 30
        integer(int32) :: fhandle
        integer(int32) :: Mfile, Yfile, Nri, Obsi, rc
        logical :: OK3

        ! 1. Prepare record
        open(newunit=fhandle, file=trim(GetTemperatureFilefull()), &
                     status='old', action='read', iostat=rc)
        read(fhandle, *, iostat=rc) ! description
        read(fhandle, *, iostat=rc) ! time step
        read(fhandle, *, iostat=rc) ! day
        read(fhandle, *, iostat=rc) ! month
        read(fhandle, *, iostat=rc) ! year
        read(fhandle, *, iostat=rc)
        read(fhandle, *, iostat=rc)
        read(fhandle, *, iostat=rc)

        Mfile = GetTemperatureRecord_FromM()
        if (GetTemperatureRecord_FromY() == 1901) then
            Yfile = Yeari
        else
            Yfile = GetTemperatureRecord_FromY()
        end if
        OK3 = .false.

        ! 2. IF 3 or less records
        if (GetTemperatureRecord_NrObs() <= 3) then
            call ReadMonth(C1Min, C1Max, fhandle, rc)
            X1 = n1
            select case (GetTemperatureRecord_NrObs())
            case (0)
                t1 = X1
                X2 = X1 + n1
                C2Min = C1Min
                C2Max = C1Max
                X3 = X2 + n1
                C3Min = C1Min
                C3Max = C1Max
            case (1)
                t1 = X1
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
                call ReadMonth(C3Min, C3Max, fhandle, rc)
                if (Monthi == Mfile) then
                    C2Min = C3Min
                    C2Max = C3Max
                    X2 = X1 + n3
                    X3 = X2 + n3
                else
                    C2Min = C1Min
                    C2Max = C1Max
                    X2 = X1 + n1
                    X3 = X2 + n3
               end if
            case (2)
               if (Monthi == Mfile) then
                   t1 = 0
               end if
               Mfile = Mfile + 1
               if (Mfile > 12) then
                   call AdjustMONTHandYEAR(Mfile, Yfile)
               end if
               call ReadMonth(C2Min, C2Max, fhandle, rc)
               X2 = X1 + n2
               if (Monthi == Mfile) then
                   t1 = X1
               end if
               Mfile = Mfile + 1
               if (Mfile > 12) then
                   call AdjustMONTHandYEAR(Mfile, Yfile)
               end if
               call ReadMonth(C3Min, C3Max, fhandle, rc)
               X3 = X2 + n3
               if (Monthi == Mfile) then
                   t1 = X2
               end if
           end select
           OK3 = .true.
        end if

        ! 3. If first observation
        if ((.not. OK3) .and. ((Monthi == Mfile) .and. (Yeari == Yfile))) then
            t1 = 0
            call ReadMonth(C1Min, C1Max, fhandle, rc)
            X1 = n1
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            call ReadMonth(C2Min, C2Max, fhandle, rc)
            X2 = X1 + n2
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            call ReadMonth(C3Min, C3Max, fhandle, rc)
            X3 = X2 + n3
            OK3 = .true.
        end if

        ! 4. If last observation
        if ((.not. OK3) .and. (Monthi == GetTemperatureRecord_ToM())) then
            if ((GetTemperatureRecord_FromY() == 1901) &
                .or. (Yeari == GetTemperatureRecord_ToY())) then
                do Nri = 1, (GetTemperatureRecord_NrObs()-3)
                    read(fhandle, *, iostat=rc)
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                end do
                call ReadMonth(C1Min, C1Max, fhandle, rc)
                X1 = n1
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
                call ReadMonth(C2Min, C2Max, fhandle, rc)
                X2 = X1 + n2
                t1 = X2
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
                call ReadMonth(C3Min, C3Max, fhandle, rc)
                X3 = X2 + n3
                OK3 = .true.
            end if
        end if

        ! 5. IF not previous cases
        if (.not. OK3) then
            Obsi = 1
            do while (.not. OK3)
                if ((Monthi == Mfile) .and. (Yeari == Yfile)) then
                   OK3 = .true.
                else
                   Mfile = Mfile + 1
                   if (Mfile > 12) then
                       call AdjustMONTHandYEAR(Mfile, Yfile)
                   end if
                  Obsi = Obsi + 1
                end if
            end do
            Mfile = GetTemperatureRecord_FromM()
            do Nri = 1, (Obsi-2)
                read(fhandle, *, iostat=rc)
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
            end do
            call ReadMonth(C1Min, C1Max, fhandle, rc)
            X1 = n1
            t1 = X1
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            call ReadMonth(C2Min, C2Max, fhandle, rc)
            X2 = X1 + n2
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            call ReadMonth(C3Min, C3Max, fhandle, rc)
            X3 = X2 + n3
        end if

        close(fhandle)
    end subroutine GetSetofThreeMonths


    subroutine ReadMonth(CiMin, CiMax, fhandle, rc)
        real(dp), intent(inout) :: CiMin
        real(dp), intent(inout) :: CiMax
        integer(int32), intent(in) :: fhandle
        integer(int32), intent(inout) :: rc

        integer(int32), parameter :: ni = 30
        character(len=255) :: StringREAD

        read(fhandle, '(a)', iostat=rc) StringREAD
        call SplitStringInTwoParams(StringREAD, CiMin, CiMax)
        ! simplification give better results for all cases
        CiMin = CiMin * ni
        CiMax = CiMax * ni
    end subroutine ReadMonth


    subroutine GetInterpolationParameters(C1, C2, C3, &
                          aOver3, bOver2, c)
        real(dp), intent(in) :: C1
        real(dp), intent(in) :: C2
        real(dp), intent(in) :: C3
        real(dp), intent(inout) :: aOver3
        real(dp), intent(inout) :: bOver2
        real(dp), intent(inout) :: c

        ! n1=n2=n3=30 --> better parabola
        aOver3 = (C1-2._dp*C2+C3)/(6._dp*30._dp*30._dp*30._dp)
        bOver2 = (-6._dp*C1+9._dp*C2-3._dp*C3)/(6._dp*30._dp*30._dp)
        c = (11._dp*C1-7._dp*C2+2._dp*C3)/(6._dp*30._dp)
    end subroutine GetInterpolationParameters
end subroutine GetMonthlyTemperatureDataSet


integer(int32) function GrowingDegreeDays(ValPeriod, FirstDayPeriod, Tbase, &
                                          Tupper, TDayMin, TDayMax)
    integer(int32), intent(in) :: ValPeriod
    integer(int32), intent(in) :: FirstDayPeriod
    real(dp), intent(in) :: Tbase
    real(dp), intent(in) :: Tupper
    real(dp), intent(in) :: TDayMin
    real(dp), intent(in) :: TDayMax

    integer(int32) :: i, RemainingDays
    integer(int32) :: DayNri
    real(dp)       :: GDDays, DayGDD
    type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet
    logical :: AdjustDayNri
    real(dp) :: TDayMin_local, TDayMax_local

    TDayMin_local = TDayMin
    TDayMax_local = TDayMax
    GDDays = 0._dp

    if (ValPeriod > 0) then
        if (GetTemperatureFile() == '(None)') then
            ! given average Tmin and Tmax
            DayGDD = DegreesDay(Tbase, Tupper, &
                     TDayMin_local, TDayMax_local, GetSimulParam_GDDMethod())
            GDDays = roundc(ValPeriod * DayGDD, mold=1_int32)
        else if (GetTemperatureFile() == '(External)') then
            DayNri = FirstDayPeriod
            RemainingDays = ValPeriod
            i = DayNri - GetSimulation_FromDayNr() + 1
            TDayMin_local = real(GetTminRun_i(i),kind=dp)
            TDayMax_local = real(GetTmaxRun_i(i),kind=dp)
            DayGDD = DegreesDay(Tbase, Tupper, TDayMin_local, &
                                        TDayMax_local, &
                                        GetSimulParam_GDDMethod())
            GDDays = GDDays + DayGDD
            RemainingDays = RemainingDays - 1

            do while ((RemainingDays > 0) &
                        .and. (i<(GetSimulation_ToDayNr()-GetSimulation_FromDayNr()+1)))
                        i = i + 1
                        TDayMin_local = real(GetTminRun_i(i),kind=dp)
                        TDayMax_local = real(GetTmaxRun_i(i),kind=dp)
                        DayGDD = DegreesDay(Tbase, Tupper, TDayMin_local, &
                                            TDayMax_local, &
                                            GetSimulParam_GDDMethod())
                        GDDays = GDDays + DayGDD
                        RemainingDays = RemainingDays - 1
            end do

           if (RemainingDays > 0) then
                 GDDays = undef_int
           end if
        else
            ! temperature file
            DayNri = FirstDayPeriod
            if (FullUndefinedRecord(GetTemperatureRecord_FromY(),&
                  GetTemperatureRecord_FromD(), GetTemperatureRecord_FromM(),&
                  GetTemperatureRecord_ToD(), GetTemperatureRecord_ToM())) then
                AdjustDayNri = .true.
                call SetDayNrToYundef(DayNri)
            else
                AdjustDayNri = .false.
            end if

            if (TemperatureFilefull_exists .and. &
                (GetTemperatureRecord_ToDayNr() > DayNri) .and. &
                (GetTemperatureRecord_FromDayNr() <= DayNri)) then
                RemainingDays = ValPeriod

                select case (GetTemperatureRecord_DataType())
                case (datatype_daily)
                    ! Tmin and Tmax arrays contain the TemperatureFilefull data
                    i = DayNri - GetTemperatureRecord_FromDayNr() + 1
                    TDayMin_local = Tmin(i)
                    TDayMax_local = Tmax(i)

                    DayGDD = DegreesDay(Tbase, Tupper, TDayMin_local, &
                                        TDayMax_local, &
                                        GetSimulParam_GDDMethod())
                    GDDays = GDDays + DayGDD
                    RemainingDays = RemainingDays - 1
                    DayNri = DayNri + 1

                    do while ((RemainingDays > 0) &
                        .and. ((DayNri < GetTemperatureRecord_ToDayNr()) &
                        .or. AdjustDayNri))

                        i = i + 1
                        if (i == size(Tmin)) then
                            i = 1
                        end if
                        TDayMin_local = Tmin(i)
                        TDayMax_local = Tmax(i)

                        DayGDD = DegreesDay(Tbase, Tupper, TDayMin_local, &
                                            TDayMax_local, &
                                            GetSimulParam_GDDMethod())
                        GDDays = GDDays + DayGDD
                        RemainingDays = RemainingDays - 1
                        DayNri = DayNri + 1
                    end do

                    if (RemainingDays > 0) then
                        GDDays = undef_int
                    end if

                case(datatype_decadely)
                    call GetDecadeTemperatureDataSet(DayNri, TminDataSet,&
                                TmaxDataSet)
                    i = 1
                    do while (TminDataSet(i)%DayNr /= DayNri)
                        i = i+1
                    end do
                    TDayMin_local = TminDataSet(i)%Param
                    TDayMax_local = TmaxDataSet(i)%Param
                    DayGDD = DegreesDay(Tbase, Tupper, &
                                 TDayMin_local, TDayMax_local, GetSimulParam_GDDMethod())
                    GDDays = GDDays + DayGDD
                    RemainingDays = RemainingDays - 1
                    DayNri = DayNri + 1
                    do while ((RemainingDays > 0) &
                        .and. ((DayNri < GetTemperatureRecord_ToDayNr()) &
                               .or. AdjustDayNri))
                        if (DayNri > TminDataSet(31)%DayNr) then
                            call GetDecadeTemperatureDataSet(DayNri, &
                                    TminDataSet, TmaxDataSet)
                        end if
                        i = 1
                        do while (TminDataSet(i)%DayNr /= DayNri)
                            i = i+1
                        end do
                        TDayMin_local = TminDataSet(i)%Param
                        TDayMax_local = TmaxDataSet(i)%Param
                        DayGDD = DegreesDay(Tbase, Tupper, &
                                     TDayMin_local, TDayMax_local,&
                                     GetSimulParam_GDDMethod())
                        GDDays = GDDays + DayGDD
                        RemainingDays = RemainingDays - 1
                        DayNri = DayNri + 1
                    end do
                    if (RemainingDays > 0) then
                        GDDays = undef_int
                    end if

                case(datatype_monthly)
                    call GetMonthlyTemperatureDataSet(DayNri, &
                            TminDataSet, TmaxDataSet)
                    i = 1
                    do while (TminDataSet(i)%DayNr /= DayNri)
                        i = i+1
                    end do
                    TDayMin_local = TminDataSet(i)%Param
                    TDayMax_local = TmaxDataSet(i)%Param
                    DayGDD = DegreesDay(Tbase, Tupper, &
                                 TDayMin_local, TDayMax_local, GetSimulParam_GDDMethod())
                    GDDays = GDDays + DayGDD
                    RemainingDays = RemainingDays - 1
                    DayNri = DayNri + 1
                    do while((RemainingDays > 0) &
                        .and. ((DayNri < GetTemperatureRecord_ToDayNr()) &
                        .or. AdjustDayNri))
                        if (DayNri > TminDataSet(31)%DayNr) then
                            call GetMonthlyTemperatureDataSet(DayNri, &
                                 TminDataSet, TmaxDataSet)
                        end if
                        i = 1
                        do while (TminDataSet(i)%DayNr /= DayNri)
                            i = i+1
                        end do
                        TDayMin_local = TminDataSet(i)%Param
                        TDayMax_local = TmaxDataSet(i)%Param
                        DayGDD = DegreesDay(Tbase, Tupper, &
                              TDayMin_local, TDayMax_local, GetSimulParam_GDDMethod())
                        GDDays = GDDays + DayGDD
                        RemainingDays = RemainingDays - 1
                        DayNri = DayNri + 1
                    end do
                    if (RemainingDays > 0) then
                        GDDays = undef_int
                    end if
                end select
            end if !if temperaturefull file exists
        end if !if temperature file
    else
        GDDays = undef_int
    endif !end valperiod>0
    GrowingDegreeDays = roundc(GDDays, mold=1_int32)
end function GrowingDegreeDays


integer(int32) function SumCalendarDays(ValGDDays, FirstDayCrop, Tbase, Tupper,&
                                        TDayMin, TDayMax)
    integer(int32), intent(in) :: ValGDDays
    integer(int32), intent(in) :: FirstDayCrop
    real(dp), intent(in) :: Tbase
    real(dp), intent(in) :: Tupper
    real(dp), intent(in) :: TDayMin
    real(dp), intent(in) :: TDayMax

    integer(int32) :: i
    integer(int32) :: NrCDays
    real(dp) :: RemainingGDDays, DayGDD
    integer(int32) :: DayNri
    type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet
    logical :: AdjustDayNri
    real(dp) :: TDayMin_loc, TDayMax_loc

    TDayMin_loc = TDayMin
    TDayMax_loc = TDayMax


    NrCdays = 0
    if (ValGDDays > 0) then
        if (GetTemperatureFile() == '(None)') then
            ! given average Tmin and Tmax
            DayGDD = DegreesDay(Tbase, Tupper, &
                       TDayMin_loc, TDayMax_loc, GetSimulParam_GDDMethod())
            if (abs(DayGDD) < epsilon(1._dp)) then
                NrCDays = -9
            else
                NrCDays = roundc(ValGDDays/DayGDD, mold=1_int32)
            end if
        else if (GetTemperatureFile() == '(External)') then
            DayNri = FirstDayCrop
            RemainingGDDays = ValGDDays
            i = DayNri-GetSimulation_FromDayNr()+1
            TDayMin_loc = real(GetTminRun_i(i),kind=dp)
            TDayMax_loc = real(GetTmaxRun_i(i),kind=dp)
            DayGDD = DegreesDay(Tbase, Tupper, TDayMin_loc, &
                                        TDayMax_loc, &
                                        GetSimulParam_GDDMethod())
            NrCDays = NrCDays + 1
            RemainingGDDays = RemainingGDDays - DayGDD

            do while ((RemainingGDDays > 0) &
                           .and. (i < (GetSimulation_ToDayNr()-GetSimulation_FromDayNr()+1)))
                  i = i + 1
                  TDayMin_loc = real(GetTminRun_i(i),kind=dp)
                  TDayMax_loc = real(GetTmaxRun_i(i),kind=dp)

                  DayGDD = DegreesDay(Tbase, Tupper, TDayMin_loc, &
                                       TDayMax_loc, &
                                       GetSimulParam_GDDMethod())
                  NrCDays = NrCDays + 1
                  RemainingGDDays = RemainingGDDays - DayGDD
            end do

            if (RemainingGDDays > 0) then
                NrCDays = undef_int
            end if
        else
            DayNri = FirstDayCrop
            if (FullUndefinedRecord(GetTemperatureRecord_FromY(), &
                  GetTemperatureRecord_FromD(), GetTemperatureRecord_FromM(), &
                  GetTemperatureRecord_ToD(), GetTemperatureRecord_ToM())) then
                AdjustDayNri = .true.
                call SetDayNrToYundef(DayNri)
            else
                AdjustDayNri = .false.
            end if

            if (TemperatureFilefull_exists .and. &
                (GetTemperatureRecord_ToDayNr() > DayNri) .and. &
                (GetTemperatureRecord_FromDayNr() <= DayNri)) then
                RemainingGDDays = ValGDDays

                select case (GetTemperatureRecord_DataType())
                case (datatype_daily)
                    ! Tmin and Tmax arrays contain the TemperatureFilefull data
                    i = DayNri - GetTemperatureRecord_FromDayNr() + 1
                    TDayMin_loc = Tmin(i)
                    TDayMax_loc = Tmax(i)

                    DayGDD = DegreesDay(Tbase, Tupper, TDayMin_loc, &
                                        TDayMax_loc, &
                                        GetSimulParam_GDDMethod())
                    NrCDays = NrCDays + 1
                    RemainingGDDays = RemainingGDDays - DayGDD
                    DayNri = DayNri + 1

                    do while ((RemainingGDDays > 0) &
                        .and. ((DayNri < GetTemperatureRecord_ToDayNr()) &
                        .or. AdjustDayNri))

                        i = i + 1
                        if (i == size(Tmin)) then
                            i = 1
                        end if
                        TDayMin_loc = Tmin(i)
                        TDayMax_loc = Tmax(i)

                        DayGDD = DegreesDay(Tbase, Tupper, TDayMin_loc, &
                                            TDayMax_loc, &
                                            GetSimulParam_GDDMethod())
                        NrCDays = NrCDays + 1
                        RemainingGDDays = RemainingGDDays - DayGDD
                        DayNri = DayNri + 1
                    end do

                    if (RemainingGDDays > 0) then
                        NrCDays = undef_int
                    end if

                case(datatype_decadely)
                    call GetDecadeTemperatureDataSet(DayNri, &
                      TminDataSet, TmaxDataSet)
                    i = 1
                    do while (TminDataSet(i)%DayNr /= DayNri)
                        i = i+1
                    end do
                    TDayMin_loc = TminDataSet(i)%Param
                    TDayMax_loc = TmaxDataSet(i)%Param
                    DayGDD = DegreesDay(Tbase, Tupper, &
                               TDayMin_loc, TDayMax_loc, GetSimulParam_GDDMethod())
                    NrCDays = NrCDays + 1
                    RemainingGDDays = RemainingGDDays - DayGDD
                    DayNri = DayNri + 1
                    do while ((RemainingGDDays > 0) &
                        .and. ((DayNri < GetTemperatureRecord_ToDayNr()) &
                         .or. AdjustDayNri))
                        if (DayNri > TminDataSet(31)%DayNr) then
                            call GetDecadeTemperatureDataSet(DayNri, &
                              TminDataSet, TmaxDataSet)
                        end if
                        i = 1
                        do while (TminDataSet(i)%DayNr /= DayNri)
                            i = i+1
                        end do
                        TDayMin_loc = TminDataSet(i)%Param
                        TDayMax_loc = TmaxDataSet(i)%Param
                        DayGDD = DegreesDay(Tbase, Tupper, &
                             TDayMin_loc, TDayMax_loc, GetSimulParam_GDDMethod())
                        NrCDays = NrCDays + 1
                        RemainingGDDays = RemainingGDDays - DayGDD
                        DayNri = DayNri + 1
                    end do
                    if (RemainingGDDays > 0) then
                        NrCDays = undef_int
                    end if

                case(datatype_monthly)
                    call GetMonthlyTemperatureDataSet(DayNri, &
                           TminDataSet, TmaxDataSet)
                    i = 1
                    do while (TminDataSet(i)%DayNr /= DayNri)
                        i = i+1
                    end do
                    TDayMin_loc = TminDataSet(i)%Param
                    TDayMax_loc = TmaxDataSet(i)%Param
                    DayGDD = DegreesDay(Tbase, Tupper, &
                               TDayMin_loc, TDayMax_loc, GetSimulParam_GDDMethod())
                    NrCDays = NrCDays + 1
                    RemainingGDDays = RemainingGDDays - DayGDD
                    DayNri = DayNri + 1
                    do while ((RemainingGDDays > 0) &
                        .and. ((DayNri < GetTemperatureRecord_ToDayNr()) &
                         .or. AdjustDayNri))
                        if (DayNri > TminDataSet(31)%DayNr) then
                            call GetMonthlyTemperatureDataSet(DayNri, &
                                   TminDataSet, TmaxDataSet)
                        end if
                        i = 1
                        do while (TminDataSet(i)%DayNr /= DayNri)
                            i = i+1
                        end do
                        TDayMin_loc = TminDataSet(i)%Param
                        TDayMax_loc = TmaxDataSet(i)%Param
                        DayGDD = DegreesDay(Tbase, Tupper, &
                             TDayMin_loc, TDayMax_loc, GetSimulParam_GDDMethod())
                        NrCDays = NrCDays + 1
                        RemainingGDDays = RemainingGDDays - DayGDD
                        DayNri = DayNri + 1
                    end do
                    if (RemainingGDDays > 0) then
                        NrCDays = undef_int
                    end if
                end select
            else
                NrCDays = undef_int
            endif
        endif
    endif
    SumCalendarDays = NrCDays
end function SumCalendarDays


real(dp) function MaxAvailableGDD(FromDayNr, Tbase, Tupper, TDayMin, TDayMax)
    integer(int32), intent(inout) :: FromDayNr
    real(dp), intent(in) :: Tbase
    real(dp), intent(in) :: Tupper
    real(dp), intent(inout) :: TDayMin
    real(dp), intent(inout) :: TDayMax

    integer(int32) :: i
    real(dp) :: MaxGDDays, DayGDD
    integer(int32) :: DayNri
    type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet

    MaxGDDays = 100000._dp
    if (GetTemperatureFile() == '(None)') then
        DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                       GetSimulParam_GDDMethod())
        if (DayGDD <= epsilon(1._dp)) then
            MaxGDDays = 0._dp
        end if
    else if (GetTemperatureFile() == '(External)') then
        DayNri = FromDayNr
        MaxGDDays = 0._dp
        i = DayNri - GetSimulation_FromDayNr() + 1
        TDayMin = real(GetTminRun_i(i),kind=dp)
        TDayMax = real(GetTmaxRun_i(i),kind=dp)
        DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                                    GetSimulParam_GDDMethod())
        MaxGDDays = MaxGDDays + DayGDD
        do while (i < (GetSimulation_ToDayNr()-GetSimulation_FromDayNr()+1))
            i = i + 1
            TDayMin = real(GetTminRun_i(i),kind=dp)
            TDayMax = real(GetTmaxRun_i(i),kind=dp)
            DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                                GetSimulParam_GDDMethod())
            MaxGDDays = MaxGDDays + DayGDD
        end do
    else
        MaxGDDays = 0._dp
        if (FullUndefinedRecord(GetTemperatureRecord_FromY(),&
               GetTemperatureRecord_FromD(), GetTemperatureRecord_FromM(),&
               GetTemperatureRecord_ToD(), GetTemperatureRecord_ToM())) then
            FromDayNr = GetTemperatureRecord_FromDayNr()  ! since we have 365 days anyway
        end if
        DayNri = FromDayNr

        if (TemperatureFilefull_exists .and. &
            (GetTemperatureRecord_ToDayNr() > FromDayNr) .and. &
            (GetTemperatureRecord_FromDayNr() <= FromDayNr)) then

            select case (GetTemperatureRecord_DataType())
            case (datatype_daily)
                ! Tmin and Tmax arrays contain the TemperatureFilefull data
                i = DayNri - GetTemperatureRecord_FromDayNr() + 1
                TDayMin = Tmin(i)
                TDayMax = Tmax(i)

                DayNri = DayNri + 1
                DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                                    GetSimulParam_GDDMethod())
                MaxGDDays = MaxGDDays + DayGDD

                do while (DayNri < GetTemperatureRecord_ToDayNr())
                    i = i + 1
                    if (i == size(Tmin)) then
                        i = 1
                    end if
                    TDayMin = Tmin(i)
                    TDayMax = Tmax(i)

                    DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                                        GetSimulParam_GDDMethod())
                    MaxGDDays = MaxGDDays + DayGDD
                    DayNri = DayNri + 1
                end do

            case (datatype_decadely)
                call GetDecadeTemperatureDataSet(DayNri, TminDataSet,&
                         TmaxDataSet)
                i = 1
                do while (TminDataSet(i)%DayNr /= DayNri)
                    i = i+1
                end do
                TDaymin = TminDataSet(i)%Param
                TDaymax = TmaxDataSet(i)%Param
                DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                              GetSimulParam_GDDMethod())
                MaxGDDays = MaxGDDays + DayGDD
                DayNri = DayNri + 1
                do while(DayNri < GetTemperatureRecord_ToDayNr())
                    if (DayNri > TminDataSet(31)%DayNr) then
                        call GetDecadeTemperatureDataSet(DayNri, TminDataSet,&
                                TmaxDataSet)
                    end if
                    i = 1
                    do while (TminDataSet(i)%DayNr /= DayNri)
                        i = i+1
                    end do
                    TDayMin = TminDataSet(i)%Param
                    TDayMax = TmaxDataSet(i)%Param
                    DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax,&
                                 GetSimulParam_GDDMethod())
                    MaxGDDays = MaxGDDays + DayGDD
                    DayNri = DayNri + 1
                end do

            case (datatype_monthly)
                call GetMonthlyTemperatureDataSet(DayNri, TminDataSet,&
                           TmaxDataSet)
                i = 1
                do while (TminDataSet(i)%DayNr /= DayNri)
                    i = i+1
                end do
                TDayMin = TminDataSet(i)%Param
                TDayMax = TmaxDataSet(i)%Param
                DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax,&
                             GetSimulParam_GDDMethod())
                MaxGDDays = MaxGDDays + DayGDD
                DayNri = DayNri + 1
                do while (DayNri < GetTemperatureRecord_ToDayNr())
                    if (DayNri > TminDataSet(31)%DayNr) then
                        call GetMonthlyTemperatureDataSet(DayNri, TminDataSet,&
                                  TmaxDataSet)
                    end if
                    i = 1
                    do while (TminDataSet(i)%DayNr /= DayNri)
                        i = i+1
                    end do
                    TDayMin = TminDataSet(i)%Param
                    TDayMax = TmaxDataSet(i)%Param
                    DayGDD = DegreesDay(Tbase, Tupper, TDayMin, TDayMax,&
                                 GetSimulParam_GDDMethod())
                    MaxGDDays = MaxGDDays + DayGDD
                    DayNri = DayNri + 1
                end do
            end select
        end if
    end if
    MaxAvailableGDD = MaxGDDays
end function MaxAvailableGDD


subroutine AdjustCalendarDays(PlantDayNr, InfoCropType,&
              Tbase, Tupper, NoTempFileTMin, NoTempFileTMax,&
              GDDL0, GDDL12, GDDFlor, GDDLengthFlor, GDDL123,&
              GDDHarvest, GDDLZmax, GDDHImax, GDDCGC, GDDCDC,&
              CCo, CCx, IsCGCGiven, HIndex, TheDaysToCCini, TheGDDaysToCCini,&
              ThePlanting, D0, D12, DFlor, LengthFlor,&
              D123, DHarvest, DLZmax, LHImax, StLength,&
              CGC, CDC, dHIdt, Succes)
    integer(int32), intent(in) :: PlantDayNr
    integer(intEnum), intent(in) :: InfoCropType
    real(dp), intent(in) :: Tbase
    real(dp), intent(in) :: Tupper
    real(dp), intent(in) :: NoTempFileTMin
    real(dp), intent(in) :: NoTempFileTMax
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDFlor
    integer(int32), intent(in) :: GDDLengthFlor
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDHarvest
    integer(int32), intent(in) :: GDDLZmax
    integer(int32), intent(inout) :: GDDHImax
    real(dp), intent(in) :: GDDCGC
    real(dp), intent(in) :: GDDCDC
    real(dp), intent(in) :: CCo
    real(dp), intent(in) :: CCx
    logical, intent(in) :: IsCGCGiven
    integer(int32), intent(in) :: HIndex
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(intEnum), intent(in) :: ThePlanting
    integer(int32), intent(inout) :: D0
    integer(int32), intent(inout) :: D12
    integer(int32), intent(inout) :: DFlor
    integer(int32), intent(inout) :: LengthFlor
    integer(int32), intent(inout) :: D123
    integer(int32), intent(inout) :: DHarvest
    integer(int32), intent(inout) :: DLZmax
    integer(int32), intent(inout) :: LHImax
    integer(int32), dimension(4), intent(inout) :: StLength
    real(dp), intent(inout) :: CGC
    real(dp), intent(inout) :: CDC
    real(dp), intent(inout) :: dHIdt
    logical, intent(inout) :: Succes

    real(dp) :: tmp_NoTempFileTMin, tmp_NoTempFileTMax
    integer :: ExtraDays, ExtraGDDays

    tmp_NoTempFileTMin = NoTempFileTMin
    tmp_NoTempFileTMax = NoTempFileTMax

    Succes = .true.
    if (TheDaysToCCini == 0) then
        ! planting/sowing
        D0 = SumCalendarDays(GDDL0, PlantDayNr, Tbase, Tupper, &
                             NoTempFileTMin, NoTempFileTMax)
        D12 = SumCalendarDays(GDDL12, PlantDayNr, Tbase, Tupper, &
                              NoTempFileTMin, NoTempFileTMax)
    else
        ! regrowth
        if (TheDaysToCCini > 0) THEN
           ! CCini < CCx
           ExtraGDDays = GDDL12 - GDDL0 - TheGDDaysToCCini
           ExtraDays = SumCalendarDays(ExtraGDDays, PlantDayNr, Tbase, &
                                       Tupper, NoTempFileTMin, NoTempFileTMax)
           D12 = D0 + TheDaysToCCini + ExtraDays
        end if
    end if

    if (InfoCropType /= subkind_Forage) then
        D123 = SumCalendarDays(GDDL123, PlantDayNr,&
                 Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax)
        DHarvest = SumCalendarDays(GDDHarvest, PlantDayNr,&
                     Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax)
    end if

    DLZmax = SumCalendarDays(GDDLZmax, PlantDayNr,&
               Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax)
    select case (InfoCropType)
    case (subkind_Grain, subkind_Tuber)
        DFlor = SumCalendarDays(GDDFlor, PlantDayNr,&
                  Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax)
        if (DFlor /= undef_int) then
            if (InfoCropType == subkind_Grain) then
                LengthFlor = SumCalendarDays(GDDLengthFlor, (PlantDayNr+DFlor),&
                   Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax)
            else
                LengthFlor = 0
            end if
            LHImax = SumCalendarDays(GDDHImax, (PlantDayNr+DFlor),&
                       Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax)
            if ((LengthFlor == undef_int) .or. (LHImax == undef_int)) then
                Succes = .false.
            end if
        else
            LengthFlor = undef_int
            LHImax = undef_int
            Succes = .false.
        end if
    case (subkind_Vegetative, subkind_Forage)
        LHImax = SumCalendarDays(GDDHImax, PlantDayNr,&
                   Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax)
    end select
    if ((D0 == undef_int) .or. (D12 == undef_int) .or. &
        (D123 == undef_int) .or. (DHarvest == undef_int) .or. &
        (DLZmax == undef_int)) then
        Succes = .false.
    end if

    if (Succes) then
        CGC = (real(GDDL12, kind=dp)/real(D12, kind=dp)) * GDDCGC
        call GDDCDCToCDC(PlantDayNr, D123, GDDL123, GDDHarvest,&
               CCx, GDDCDC, Tbase, Tupper, tmp_NoTempFileTMin, tmp_NoTempFileTMax, CDC)
        call DetermineLengthGrowthStages(CCo, CCx, CDC, D0, DHarvest,&
               IsCGCGiven, TheDaysToCCini, &
               ThePlanting, D123, StLength, D12, CGC)
        if ((InfoCropType == subkind_Grain) .or. (InfoCropType == subkind_Tuber)) then
            dHIdt = real(HIndex, kind=dp)/real(LHImax, kind=dp)
        end if
        if ((InfoCropType == subkind_Vegetative) &
            .or. (InfoCropType == subkind_Forage)) then
            if (LHImax > 0) then
                if (LHImax > DHarvest) then
                    dHIdt = real(HIndex, kind=dp)/real(DHarvest, kind=dp)
                else
                    dHIdt = real(HIndex, kind=dp)/real(LHImax, kind=dp)
                end if
                if (dHIdt > 100) then
                    dHIdt = 100 ! 100 is maximum TempdHIdt (See SetdHIdt)
                    LHImax = 0
                end if
            else
                dHIdt = 100 ! 100 is maximum TempdHIdt (See SetdHIdt)
                LHImax = 0
            end if
        end if
    end if
end subroutine AdjustCalendarDays


subroutine AdjustCalendarCrop(FirstCropDay)
    integer(int32), intent(in) :: FirstCropDay

    logical :: succes
    logical :: CGCisGiven
    integer(int32) :: Crop_GDDaysToHIo_temp
    integer(int32) :: Crop_DaysToGermination_temp
    integer(int32) :: Crop_DaysToFullCanopy_temp
    integer(int32) :: Crop_DaysToFlowering_temp
    integer(int32) :: Crop_LengthFlowering_temp
    integer(int32) :: Crop_DaysToSenescence_temp
    integer(int32) :: Crop_DaysToHarvest_temp
    integer(int32) :: Crop_DaysToMaxRooting_temp
    integer(int32) :: Crop_DaysToHIo_temp
    integer(int32), dimension(4) :: Crop_Length_temp
    real(dp) :: Crop_CGC_temp
    real(dp) :: Crop_CDC_temp
    real(dp) :: Crop_dHIdt_temp

    CGCisGiven = .true.

    select case (GetCrop_ModeCycle())
    case (modeCycle_GDDays)
        call SetCrop_GDDaysToFullCanopy(GetCrop_GDDaysToGermination() &
           + roundc(log((0.25_dp*GetCrop_CCx()*GetCrop_CCx()/GetCrop_CCo()) &
               /(GetCrop_CCx()-(0.98_dp*GetCrop_CCx())))/GetCrop_GDDCGC(), &
               mold=int32))
        if (GetCrop_GDDaysToFullCanopy() > GetCrop_GDDaysToHarvest()) then
            call SetCrop_GDDaysToFullCanopy(GetCrop_GDDaysToHarvest())
        end if
        Crop_GDDaysToHIo_temp = GetCrop_GDDaysToHIo()
        Crop_DaysToGermination_temp = GetCrop_DaysToGermination()
        Crop_DaysToFullCanopy_temp = GetCrop_DaysToFullCanopy()
        Crop_DaysToFlowering_temp = GetCrop_DaysToFlowering()
        Crop_LengthFlowering_temp = GetCrop_LengthFlowering()
        Crop_DaysToSenescence_temp = GetCrop_DaysToSenescence()
        Crop_DaysToHarvest_temp = GetCrop_DaysToHarvest()
        Crop_DaysToMaxRooting_temp = GetCrop_DaysToMaxRooting()
        Crop_DaysToHIo_temp = GetCrop_DaysToHIo()
        Crop_Length_temp = GetCrop_Length()
        Crop_CGC_temp = GetCrop_CGC()
        Crop_CDC_temp = GetCrop_CDC()
        Crop_dHIdt_temp = GetCrop_dHIdt()
        call AdjustCalendarDays(FirstCropDay, GetCrop_subkind(), &
          GetCrop_Tbase(), GetCrop_Tupper(), &
          GetSimulParam_Tmin(), GetSimulParam_Tmax(), &
          GetCrop_GDDaysToGermination(), GetCrop_GDDaysToFullCanopy(), &
          GetCrop_GDDaysToFlowering(), GetCrop_GDDLengthFlowering(), &
          GetCrop_GDDaysToSenescence(), GetCrop_GDDaysToHarvest(), &
          GetCrop_GDDaysToMaxRooting(), Crop_GDDaysToHIo_temp, &
          GetCrop_GDDCGC(), GetCrop_GDDCDC(), GetCrop_CCo(), &
          GetCrop_CCx(), CGCisGiven, GetCrop_HI(), &
          GetCrop_DaysToCCini(), GetCrop_GDDaysToCCini(), GetCrop_Planting(), &
          Crop_DaysToGermination_temp, Crop_DaysToFullCanopy_temp,&
          Crop_DaysToFlowering_temp, Crop_LengthFlowering_temp, &
          Crop_DaysToSenescence_temp, Crop_DaysToHarvest_temp, &
          Crop_DaysToMaxRooting_temp, Crop_DaysToHIo_temp,&
          Crop_Length_temp, Crop_CGC_temp, &
          Crop_CDC_temp, Crop_dHIdt_temp, Succes)
        call SetCrop_GDDaysToHIo(Crop_GDDaysToHIo_temp)
        call SetCrop_DaysToGermination(Crop_DaysToGermination_temp)
        call SetCrop_DaysToFullCanopy(Crop_DaysToFullCanopy_temp)
        call SetCrop_DaysToFlowering(Crop_DaysToFlowering_temp)
        call SetCrop_LengthFlowering(Crop_LengthFlowering_temp)
        call SetCrop_DaysToSenescence(Crop_DaysToSenescence_temp)
        call SetCrop_DaysToHarvest(Crop_DaysToHarvest_temp)
        call SetCrop_DaysToMaxRooting(Crop_DaysToMaxRooting_temp)
        call SetCrop_DaysToHIo(Crop_DaysToHIo_temp)
        call SetCrop_Length(Crop_Length_temp)
        call SetCrop_CGC(Crop_CGC_temp)
        call SetCrop_CDC(Crop_CDC_temp)
        call SetCrop_dHIdt(Crop_dHIdt_temp)
    case default
        Succes = .true.
    end select
end subroutine AdjustCalendarCrop


subroutine GDDCDCToCDC(PlantDayNr, D123, GDDL123, &
                       GDDHarvest, CCx, GDDCDC, Tbase, Tupper, &
                       NoTempFileTMin, NoTempFileTMax, CDC)
    integer(int32), intent(in) :: PlantDayNr
    integer(int32), intent(in) :: D123
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDHarvest
    real(dp), intent(in) :: CCx
    real(dp), intent(in) :: GDDCDC
    real(dp), intent(in) :: Tbase
    real(dp), intent(in) :: Tupper
    real(dp), intent(in) :: NoTempFileTMin
    real(dp), intent(in) :: NoTempFileTMax
    real(dp), intent(inout) :: CDC

    integer(int32) :: ti, GDDi
    real(dp) :: CCi

    GDDi = LengthCanopyDecline(CCx, GDDCDC)
    if ((GDDL123+GDDi) <= GDDHarvest) then
        CCi = 0._dp ! full decline
    else
        ! partly decline
        if (GDDL123 < GDDHarvest) then
            GDDi = GDDHarvest - GDDL123
        else
            GDDi = 5._dp
        end if
        CCi = CCx * (1._dp - 0.05_dp  &
                 * (exp(real(GDDi,kind=dp)*(GDDCDC*3.33_dp)/(CCx+2.29_dp))-1._dp) )
       ! CC at time ti
    end if
    ti = SumCalendarDays(GDDi, (PlantDayNr+D123),&
              Tbase, Tupper, NoTempFileTMin, NoTempFileTMax)
    if (ti > 0) then
        CDC = (((CCx+2.29_dp)/real(ti, kind=dp)) &
                * log(1._dp + ((1._dp-CCi/CCx)/0.05_dp)))/3.33_dp
    else
        CDC = undef_int
    end if
end subroutine GDDCDCToCDC


integer(int32) function RoundedOffGDD(PeriodGDD, PeriodDay,&
           FirstDayPeriod, TempTbase, TempTupper, TempTmin, TempTmax)
    integer(int32), intent(in) :: PeriodGDD
    integer(int32), intent(in) :: PeriodDay
    integer(int32), intent(in) :: FirstDayPeriod
    real(dp), intent(in) :: TempTbase
    real(dp), intent(in) :: TempTupper
    real(dp), intent(in) :: TempTmin
    real(dp), intent(in) :: TempTmax

    integer(int32) :: DayMatch, PeriodUpdatedGDD
    real(dp) :: TempTmin_t, TempTmax_t

    TempTmin_t = TempTmin
    TempTmax_t = TempTmax

    if (PeriodGDD > 0) then
        DayMatch = SumCalendarDays(PeriodGDD, FirstDayPeriod, &
                     TempTbase, TempTupper, TempTmin_t, TempTmax_t)
        PeriodUpdatedGDD = GrowingDegreeDays(PeriodDay, FirstDayPeriod, &
                     TempTbase, TempTupper, TempTmin_t, TempTmax_t)
        if (PeriodDay == DayMatch) then
            RoundedOffGDD = PeriodGDD
        else
            RoundedOffGDD = PeriodUpdatedGDD
        end if
    else
        RoundedOffGDD = GrowingDegreeDays(PeriodDay, FirstDayPeriod,&
                     TempTbase, TempTupper, TempTmin_t, TempTmax_t)
    end if
end function RoundedOffGDD

integer(int32) function ResetCropDay1(CropDay1IN, SwitchToYear1)
    integer(int32), intent(in) :: CropDay1IN
    logical, intent(in) :: SwitchToYear1

    integer(int32) :: CropDay1OUT
    integer(int32) :: dayi, monthi, yeari

    call DetermineDate(CropDay1IN, dayi, monthi, yeari)
    if (GetTemperatureRecord_FromY() == 1901) then
        yeari = 1901
        call DetermineDayNr(Dayi, Monthi, Yeari, CropDay1OUT)
    else
        if (SwitchToYear1) then
            call DetermineDayNr(Dayi, Monthi, &
                 GetTemperatureRecord_FromY(), CropDay1OUT)
        else
            CropDay1OUT = CropDay1IN
        end if
    end if
    ResetCropDay1 = CropDay1OUT
end function ResetCropDay1


subroutine CropStressParametersSoilSalinity(CCxRed, CCdistortion, &
             CCo, CCx, CGC, GDDCGC, CropDeterm, L12, LFlor, &
             LengthFlor, L123, GDDL12, GDDLFlor, GDDLengthFlor, &
             GDDL123, TheModeCycle, StressResponse)
    integer(int8), intent(in) :: CCxRed
    integer(int8), intent(in) :: CCdistortion
    real(dp), intent(in) :: CCo
    real(dp), intent(in) :: CCx
    real(dp), intent(in) :: CGC
    real(dp), intent(in) :: GDDCGC
    logical, intent(in) :: CropDeterm
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: LFlor
    integer(int32), intent(in) :: LengthFlor
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDLFlor
    integer(int32), intent(in) :: GDDLengthFlor
    integer(int32), intent(in) :: GDDL123
    integer(intEnum), intent(in) :: TheModeCycle
    type(rep_EffectStress), intent(inout) :: StressResponse

    real(dp) :: CCToReach, CCxAdj, L12Double, L12SS, &
                CGCadjMax, CGCAdjMin, L12SSmax, CGCadj, CCxFinal, &
                GDDL12Double, GDDCGCadjMax, GDDL12SSmax, GDDCGCAdjMin, GDDCGCadj

    ! initialize
    StressResponse%RedCCX = CCxRed
    StressResponse%RedWP = 0_int8
    L12Double = L12
    L12SSmax = L12
    GDDL12Double = GDDL12

    ! CGC reduction
    CCToReach = 0.98_dp * CCx
    if ((CCo > CCToReach) .or. (CCo >= CCx) .or. (CCxRed == 0)) then
        StressResponse%RedCGC = 0_int8
    else
        StressResponse%RedCGC = undef_int
        ! reference for no salinity stress
        if (TheModeCycle == modeCycle_CalendarDays) then
            L12Double = log((0.25_dp*CCx*CCx/CCo)/(CCx-CCToReach))/CGC
            if (L12Double <= epsilon(1._dp)) then
                StressResponse%RedCGC = 0_int8
            end if
        else
            GDDL12Double = log((0.25_dp*CCx*CCx/CCo)/(CCx-CCToReach))/GDDCGC
            if (GDDL12Double <= epsilon(1._dp)) then
                StressResponse%RedCGC = 0_int8
            end if
        end if
        ! with salinity stress
        CCxAdj = 0.90_dp * CCx * (1._dp - CCxRed/100._dp)
        CCToReach = 0.98_dp * CCxAdj
        if ((StressResponse%RedCGC /= 0) .and. &
            ((CCxAdj-CCToReach) >= 0.0001_dp)) then
            if (TheModeCycle == modeCycle_CalendarDays) then
                CGCadjMax = log((0.25_dp*CCxAdj*CCxAdj/CCo)&
                                /(CCxAdj-CCToReach))/L12Double
                L12SSmax = L12 + (L123 - L12)/2._dp
                if (CropDeterm .and. (L12SSmax > &
                       (LFlor + roundc(LengthFlor/2._dp, mold=1_int32)))) then
                    L12SSmax = LFlor + roundc(LengthFlor/2._dp, mold=1_int32)
                end if
                if (L12SSmax > L12Double) then
                    CGCAdjMin = log((0.25_dp*CCxAdj*CCxAdj/CCo)&
                                    /(CCxAdj-CCToReach))/L12SSmax
                else
                    CGCAdjMin = CGCadjMax
                end if
                if (CCxRed < 10) then ! smooth start required
                    CGCadj = CGCadjMax - (CGCadjMax-CGCAdjMin)&
                     *(exp(CCxRed*log(1.5_dp))/exp(10*log(1.5_dp)))&
                     *(CCdistortion/100._dp)
                else
                    CGCadj = CGCadjMax - (CGCadjMax-CGCAdjMin)&
                                         *(CCdistortion/100._dp)
                end if
                StressResponse%RedCGC = &
                  roundc(100._dp*(CGC-CGCadj)/CGC, mold=1_int8)
            else
                GDDCGCadjMax = log((0.25_dp*CCxAdj*CCxAdj/CCo) &
                                   /(CCxAdj-CCToReach))/GDDL12Double
                GDDL12SSmax = GDDL12 + (GDDL123 - GDDL12)/2._dp
                if (CropDeterm .and. (GDDL12SSmax > &
                      (GDDLFlor + roundc(LengthFlor/2._dp, mold=1_int32)))) then
                    GDDL12SSmax = GDDLFlor + &
                       roundc(GDDLengthFlor/2._dp, mold=1_int32)
                end if
                if (GDDL12SSmax > GDDL12Double) then
                    GDDCGCAdjMin = log((0.25_dp*CCxAdj*CCxAdj/CCo) &
                                       /(CCxAdj-CCToReach))/GDDL12SSmax
                else
                    GDDCGCAdjMin = GDDCGCadjMax
                end if
                if (CCxRed < 10) then ! smooth start required
                    GDDCGCadj = GDDCGCadjMax - (GDDCGCadjMax-GDDCGCAdjMin)*&
                                   (exp(real(CCxRed, kind=dp))/&
                                    exp(10._dp))*(CCdistortion/100._dp)
                else
                    GDDCGCadj = GDDCGCadjMax - &
                      (GDDCGCadjMax-GDDCGCAdjMin)*(CCdistortion/100._dp)
                end if
                StressResponse%RedCGC = &
                      roundc(100._dp*(GDDCGC-GDDCGCadj)/GDDCGC, mold=1_int8)
           end if
        else
            StressResponse%RedCGC = 0_int8
        end if
    end if

    ! Canopy decline
    if (CCxRed == 0) then
        StressResponse%CDecline = 0._dp
    else
        CCxAdj = 0.98_dp*CCx*(1._dp - CCxRed/100._dp)
        L12SS = L12SSmax - (L12SSmax-L12Double) * (CCdistortion/100._dp)
        if ((L123 > L12SS) .and. (CCdistortion > 0)) then
            if (CCxRed < 10) then ! smooth start required
                CCxFinal = CCxAdj - &
                  (exp(CCxRed*log(1.5_dp))/exp(10._dp*log(1.5_dp)))* &
                  (0.5_dp*CCdistortion/100._dp)*(CCxAdj - CCo)
            else
                CCxFinal = CCxAdj - (0.5_dp*CCdistortion/100._dp)*(CCxAdj - CCo)
            end if
            if (CCxFinal < CCo) then
                CCxFinal = CCo
            end if
            StressResponse%CDecline = &
                100._dp*(CCxAdj - CCxFinal)/real(L123 - L12SS, kind=dp)
            if (StressResponse%CDecline > 1) then
                StressResponse%CDecline = 1.0_dp
            end if
            if (StressResponse%CDecline <= epsilon(1._dp)) then
                StressResponse%CDecline = 0.001_dp
            end if
        else
            StressResponse%CDecline = 0.001_dp ! no shift of maturity
        end if
    end if

    ! Stomata closure
    StressResponse%RedKsSto = CCxRed
end subroutine CropStressParametersSoilSalinity


subroutine TemperatureFileCoveringCropPeriod(CropFirstDay, CropLastDay)
    integer(int32), intent(in) :: CropFirstDay
    integer(int32), intent(in) :: CropLastDay

    character(len=:), allocatable :: totalnameOUT
    integer(int32) :: fhandle
    integer(int32) :: i, RunningDay
    type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet
    real(dp) :: Tlow, Thigh

    if (TemperatureFilefull_exists) then
        ! open file and find first day of cropping period
        select case (GetTemperatureRecord_DataType())
        case (datatype_daily)
            ! Tmin and Tmax arrays contain the TemperatureFilefull data
            i = CropFirstDay - GetTemperatureRecord_FromDayNr() + 1
            Tlow = Tmin(i)
            Thigh = Tmax(i)

        case (datatype_decadely)
            call GetDecadeTemperatureDataSet(CropFirstDay, TminDataSet, &
                        TmaxDataSet)
            i = 1
            do while (TminDataSet(i)%DayNr /= CropFirstDay)
                i = i+1
            end do
            Tlow = TminDataSet(i)%Param
            Thigh = TmaxDataSet(i)%Param

        case (datatype_monthly)
            call GetMonthlyTemperatureDataSet(CropFirstDay, TminDataSet, TmaxDataSet)
            i = 1
            do while (TminDataSet(i)%DayNr /= CropFirstDay)
                i = i+1
            end do
            Tlow = TminDataSet(i)%Param
            Thigh = TmaxDataSet(i)%Param
        end select

        ! create SIM file and record first day
        totalnameOUT = trim(GetPathNameSimul()) // 'TCrop.SIM'
        open(newunit=fhandle, file=trim(totalnameOUT), &
             action='write')
        write(fhandle, '(2f10.4)') Tlow, Thigh

        ! next days of simulation period
        do RunningDay = (CropFirstDay + 1), CropLastDay
            select case (GetTemperatureRecord_DataType())
            case (datatype_daily)
                i = i + 1
                if (i == size(Tmin)) then
                    i = 1
                end if
                Tlow = Tmin(i)
                Thigh = Tmax(i)

            case (datatype_decadely)
                if (RunningDay > TminDataSet(31)%DayNr) then
                    call GetDecadeTemperatureDataSet(RunningDay, TminDataSet,&
                        TmaxDataSet)
                end if
                i = 1
                do while (TminDataSet(i)%DayNr /= RunningDay)
                    i = i+1
                end do
                Tlow = TminDataSet(i)%Param
                Thigh = TmaxDataSet(i)%Param

            case (datatype_monthly)
               if (RunningDay > TminDataSet(31)%DayNr) then
                    call GetMonthlyTemperatureDataSet(RunningDay, TminDataSet,&
                        TmaxDataSet)
               end if
               i = 1
               do while (TminDataSet(i)%DayNr /= RunningDay)
                   i = i+1
               end do
               Tlow = TminDataSet(i)%Param
               Thigh = TmaxDataSet(i)%Param
            end select

            write(fhandle, '(2f10.4)') Tlow, Thigh
        end do

        close(fhandle)
    else
        if (GetTemperatureFile() /= '(External)') then
           write(*,*) 'ERROR: no valid air temperature file'
           return
           ! fatal error if no air temperature file
        end if
    endif
end subroutine TemperatureFileCoveringCropPeriod


subroutine AdjustCropFileParameters(TheCropFileSet, LseasonDays,&
                TheCropDay1, TheModeCycle, TheTbase, TheTupper,&
                L123, L1234, GDD123, GDD1234)
    type(rep_CropFileSet), intent(in) :: TheCropFileSet
    integer(int32), intent(in) :: LseasonDays
    integer(int32), intent(in) :: TheCropDay1
    integer(intEnum), intent(in) :: TheModeCycle
    real(dp), intent(in) :: TheTbase
    real(dp), intent(in) :: TheTupper
    integer(int32), intent(inout) :: L123
    integer(int32), intent(inout) :: L1234
    integer(int32), intent(inout) :: GDD123
    integer(int32), intent(inout) :: GDD1234

    real(dp) :: Tmin_tmp, Tmax_tmp

    ! Adjust some crop parameters (CROP.*) as specified by the generated length
    ! season (LseasonDays)
    ! time to maturity
    L1234 = LseasonDays ! days
    if (TheModeCycle == modeCycle_GDDays) then
        Tmin_tmp = GetSimulParam_Tmin()
        Tmax_tmp = GetSimulParam_Tmax()
        GDD1234 = GrowingDegreeDays(LseasonDays, TheCropDay1,&
                       TheTbase, TheTupper, &
                       Tmin_tmp, Tmax_tmp)
    else
        GDD1234 = undef_int
    end if

    ! time to senescence  (reference is given in TheCropFileSet
    if (TheModeCycle == modeCycle_GDDays) then
        GDD123 = GDD1234 - TheCropFileSet%GDDaysFromSenescenceToEnd
        if (GDD123 >= GDD1234) then
            GDD123 = GDD1234
            L123 = LseasonDays
        else
            Tmin_tmp = GetSimulParam_Tmin()
            Tmax_tmp = GetSimulParam_Tmax()
            L123 = SumCalendarDays(GDD123, TheCropDay1, TheTbase, TheTupper, &
                                   Tmin_tmp, Tmax_tmp)
        end if
    else
        L123 = L1234 - TheCropFileSet%DaysFromSenescenceToEnd
        if (L123 >= L1234) L123 = LseasonDays
        GDD123 = undef_int
    end if
end subroutine AdjustCropFileParameters


subroutine LoadSimulationRunProject(NrRun)
    integer(int32), intent(in) :: NrRun

    integer(int32) :: fClim, i, rc
    character(len=1025) :: TempString, TempString1, TempString2
    character(len=1025) :: observations_descr, eto_descr
    character(len=1025) :: CO2descr, rain_descr
    character(len=1025) :: CalendarDescriptionLocal
    character(len=1025) :: TemperatureDescriptionLocal

    real(dp) :: TotDepth
    integer(int32)  :: FertStress
    type(rep_clim) :: temperature_record
    integer(int8)  :: RedCGC_temp, RedCCX_temp
    type(CompartmentIndividual), dimension(max_No_compartments) :: &
                      Compartment_temp
    integer(intEnum) :: Crop_Planting_temp
    real(dp) :: Crop_CCini_temp, Crop_RootMin_temp, Crop_SizePlant_temp
    integer(int32) :: Crop_DaysToCCini_temp, Crop_GDDaysToCCini_temp
    integer(int32) :: Crop_DaysToSenescence_temp, Crop_DaysToHarvest_temp
    integer(int32) :: Crop_GDDaysToSenescence_temp, Crop_GDDaysToHarvest_temp
    integer(int32) :: Crop_Day1_temp
    integer(int32) :: Crop_DayN_temp
    integer(int32) :: Crop_DaysToFullCanopySF_temp
    integer(int32) :: ZiAqua_temp
    type(rep_clim) :: etorecord_tmp, rainrecord_tmp
    real(dp)       :: ECiAqua_temp, SurfaceStorage_temp

    ! 0. Year of cultivation and Simulation and Cropping period
    call SetSimulation_YearSeason(ProjectInput(NrRun)%Simulation_YearSeason)
    call SetCrop_Day1(ProjectInput(NrRun)%Crop_Day1)
    call SetCrop_DayN(ProjectInput(NrRun)%Crop_DayN)
    call SetSimulation_FromDayNr(ProjectInput(NrRun)%Simulation_DayNr1)
    call SetSimulation_ToDayNr(ProjectInput(NrRun)%Simulation_DayNrN)

    ! 1. Climate
    call SetClimateFile(ProjectInput(NrRun)%Climate_Filename)
    if ((GetClimateFile() == '(None)') .or. (GetClimateFile() == '(External)')) then
        call SetClimateFileFull(GetClimateFile())
    else
        call SetClimateFileFull(ProjectInput(NrRun)%Climate_Directory &
                                // GetClimateFile())
        open(newunit=fClim, file=trim(GetClimateFileFull()), &
             status='old', action='read', iostat=rc)
        ! 1.0 Description
        read(fClim, '(a)', iostat=rc) TempString
        call SetClimateDescription(trim(TempString))
        close(fClim)
    end if

    ! 1.1 Temperature
    call SetTemperatureFile(ProjectInput(NrRun)%Temperature_Filename)

    if ((GetTemperatureFile() == '(None)') .or. &
        (GetTemperatureFile() == '(External)')) then
        call SetTemperatureFilefull(GetTemperatureFile())  ! no file
        write(TempString1,'(f8.1)') GetSimulParam_Tmin()
        write(TempString2,'(f8.1)') GetSimulParam_Tmax()
        call SetTemperatureDescription(('Default temperature data: Tmin = '// &
          trim(TempString1)// ' and Tmax = '// trim(TempString2) // ' deg'))
    else
        call SetTemperatureFilefull(ProjectInput(NrRun)%Temperature_Directory &
                                    // GetTemperatureFile())
        TemperatureFilefull_exists = FileExists(GetTemperatureFilefull())
        if (TemperatureFilefull_exists) call ReadTemperatureFilefull()

        temperature_record = GetTemperatureRecord()
        TemperatureDescriptionLocal = GetTemperatureDescription()
        call LoadClim(GetTemperatureFileFull(), TemperatureDescriptionLocal,&
                       temperature_record)
        call SetTemperatureDescription(TemperatureDescriptionLocal)
        call CompleteClimateDescription(temperature_record)
        call SetTemperatureRecord(temperature_record)
    end if

    ! Create Temperature Reference file 
    if (GetTemperatureFile() /= '(External)') then
        call CreateTnxReferenceFile(GetTemperatureFile(),GetTnxReferenceFile(),GetTnxReferenceYear());
    end if
    call CreateTnxReference365Days();

    ! 1.2 ETo
    call SetEToFile(ProjectInput(NrRun)%ETo_Filename)
    if ((GetEToFile() == '(None)') .or. &
        (GetEToFile() == '(External)')) then
        call SetEToFilefull(GetEToFile())  ! no file
        call SetEToDescription('Specify ETo data when Running AquaCrop')
    else
        call SetEToFilefull(ProjectInput(NrRun)%ETo_Directory // GetEToFile())
        eto_descr = GetEToDescription()
        etorecord_tmp = GetEToRecord()
        call LoadClim(GetEToFilefull(), eto_descr, etorecord_tmp)
        call SetEToDescription(eto_descr)
        call CompleteClimateDescription(etorecord_tmp)
        call SetEToRecord(etorecord_tmp)
    end if

    ! 1.3 Rain
    call SetRainFile(ProjectInput(NrRun)%Rain_Filename)
    if ((GetRainFile() == '(None)') .or. &
        (GetRainFile() == '(External)')) then
        call SetRainFilefull(GetRainFile())  ! no file
        call SetRainDescription('Specify Rain data when Running AquaCrop')
    else
        call SetRainFileFull(ProjectInput(NrRun)%Rain_Directory &
                             // GetRainFile())
        rain_descr = Getraindescription()
        rainrecord_tmp = GetRainRecord()
        call LoadClim(GetRainFilefull(), rain_descr, rainrecord_tmp)
        call SetRainDescription(rain_descr)
        call CompleteClimateDescription(rainrecord_tmp)
        call SetRainRecord(rainrecord_tmp)
    end if

    ! 1.4 CO2
    call SetCO2File(ProjectInput(NrRun)%CO2_Filename)
    if (GetCO2File() /= '(None)') then
        call SetCO2FileFull(ProjectInput(NrRun)%CO2_Directory &
                            // GetCO2File())
        CO2descr =  GetCO2Description()
        call GenerateCO2Description(GetCO2FileFull(), CO2descr)
        call SetCO2Description(CO2descr)
    end if
    if (GetClimateFile() /= '(External)') then
        call SetClimData()
    end if
    call AdjustOnsetSearchPeriod() ! Set initial StartSearch and StopSearchDayNr

    ! 2. Calendar
    call SetCalendarFile(trim(ProjectInput(NrRun)%Calendar_Filename))
    if (GetCalendarFile() == '(None)') then
        call SetCalendarDescription('No calendar for the Seeding/Planting year')
    else
        call SetCalendarFileFull(ProjectInput(NrRun)%Calendar_Directory &
                                 // GetCalendarFile())
        CalendarDescriptionLocal = GetCalendarDescription()
        call GetFileDescription(GetCalendarFileFull(), CalendarDescriptionLocal)
        call SetCalendarDescription(CalendarDescriptionLocal)
    end if

    ! 3. Crop
    call SetSimulation_LinkCropToSimPeriod(.true.)
    call SetCropFile(ProjectInput(NrRun)%Crop_Filename)
    call SetCropFilefull(ProjectInput(NrRun)%Crop_Directory // GetCropFile())
    call LoadCrop(GetCropFilefull())

    ! Adjust crop parameters of Perennials
    if (GetCrop_subkind() == subkind_Forage) then
        ! adjust crop characteristics to the Year (Seeding/Planting or
        ! Non-seesing/Planting year)
        Crop_Planting_temp = GetCrop_Planting()
        Crop_RootMin_temp = GetCrop_RootMin()
        Crop_SizePlant_temp = GetCrop_SizePlant()
        Crop_CCini_temp = GetCrop_CCini()
        Crop_DaysToCCini_temp = GetCrop_DaysToCCini()
        Crop_GDDaysToCCini_temp = GetCrop_GDDaysToCCini()
        call AdjustYearPerennials(GetSimulation_YearSeason(),&
              GetCrop_SownYear1(), GetCrop_ModeCycle(), &
              GetCrop_RootMax(), GetCrop_RootMinYear1(), &
              GetCrop_CCo(), GetCrop_SizeSeedling(), GetCrop_CGC(),&
              GetCrop_CCx(), GetCrop_GDDCGC(), GetCrop_PlantingDens(), &
              Crop_Planting_temp, Crop_RootMin_temp, Crop_SizePlant_temp,&
              Crop_CCini_temp, Crop_DaysToCCini_temp, Crop_GDDaysToCCini_temp)
        call SetCrop_Planting(Crop_Planting_temp)
        call SetCrop_RootMin(Crop_RootMin_temp)
        call SetCrop_SizePlant(Crop_SizePlant_temp)
        call SetCrop_CCini(Crop_CCini_temp)
        call SetCrop_DaysToCCini(Crop_DaysToCCini_temp)
        call SetCrop_GDDaysToCCini(Crop_GDDaysToCCini_temp)
        ! adjust length of season
        call  SetCrop_DaysToHarvest(GetCrop_DayN() - GetCrop_Day1() + 1)
        Crop_DaysToSenescence_temp = GetCrop_DaysToSenescence()
        Crop_DaysToHarvest_temp = GetCrop_DaysToHarvest()
        Crop_GDDaysToSenescence_temp = GetCrop_GDDaysToSenescence()
        Crop_GDDaysToHarvest_temp = GetCrop_GDDaysToHarvest()
        call AdjustCropFileParameters(GetCropFileSet(),&
              GetCrop_DaysToHarvest(), GetCrop_Day1(), &
              GetCrop_ModeCycle(), GetCrop_Tbase(), GetCrop_Tupper(),&
              Crop_DaysToSenescence_temp, Crop_DaysToHarvest_temp,&
              Crop_GDDaysToSenescence_temp, Crop_GDDaysToHarvest_temp)
        call SetCrop_DaysToSenescence(Crop_DaysToSenescence_temp)
        call SetCrop_DaysToHarvest(Crop_DaysToHarvest_temp)
        call SetCrop_GDDaysToSenescence(Crop_GDDaysToSenescence_temp)
        call SetCrop_GDDaysToHarvest(Crop_GDDaysToHarvest_temp)
    end if

    call AdjustCalendarCrop(GetCrop_Day1())
    call CompleteCropDescription
    ! Onset.Off := true;
    if (GetClimFile() == '(None)') then
        Crop_Day1_temp = GetCrop_Day1()
        Crop_DayN_temp = GetCrop_DayN()
        call AdjustCropYearToClimFile(Crop_Day1_temp, Crop_DayN_temp)
        ! adjusting Crop.Day1 and Crop.DayN to ClimFile
        call SetCrop_Day1(Crop_Day1_temp)
        call SetCrop_DayN(Crop_DayN_temp)
    else
        call SetCrop_DayN(GetCrop_Day1() + GetCrop_DaysToHarvest() - 1)
    end if

    ! adjusting ClimRecord.'TO' for undefined year with 365 days
    if ((GetClimFile() /= '(None)') .and. (GetClimRecord_FromY() == 1901) &
        .and. (GetClimRecord_NrObs() == 365)) then
        call AdjustClimRecordTo(GetCrop_DayN())
    end if
    ! adjusting simulation period
    call AdjustSimPeriod

    ! 4. Irrigation
    call SetIrriFile(ProjectInput(NrRun)%Irrigation_Filename)
    if (GetIrriFile() == '(None)') then
        call SetIrriFileFull(GetIrriFile())
        call NoIrrigation
        ! IrriDescription := 'Rainfed cropping';
    else
        call SetIrriFileFull(ProjectInput(NrRun)%Irrigation_Directory &
                             // GetIrriFile())
        call LoadIrriScheduleInfo(GetIrriFileFull())
    end if

    ! 5. Field Management
    call SetManFile(ProjectInput(NrRun)%Management_Filename)
    if (GetManFile() == '(None)') then
        call SetManFileFull(GetManFile())
        call SetManDescription('No specific field management')
    else
        call SetManFileFull(ProjectInput(NrRun)%Management_Directory &
                            // GetManFile())
        call LoadManagement(GetManFilefull())
        ! reset canopy development to soil fertility
        FertStress = GetManagement_FertilityStress()
        Crop_DaysToFullCanopySF_temp = GetCrop_DaysToFullCanopySF()
        RedCGC_temp = GetSimulation_EffectStress_RedCGC()
        RedCCX_temp = GetSimulation_EffectStress_RedCCX()
        call TimeToMaxCanopySF(GetCrop_CCo(), GetCrop_CGC(), GetCrop_CCx(),&
               GetCrop_DaysToGermination(), GetCrop_DaysToFullCanopy(),&
               GetCrop_DaysToSenescence(), GetCrop_DaysToFlowering(),&
               GetCrop_LengthFlowering(), GetCrop_DeterminancyLinked(),&
               Crop_DaysToFullCanopySF_temp, RedCGC_temp,&
               RedCCX_temp, FertStress)
        call SetCrop_DaysToFullCanopySF(Crop_DaysToFullCanopySF_temp)
        call SetManagement_FertilityStress(FertStress)
        call SetSimulation_EffectStress_RedCGC(RedCGC_temp)
        call SetSimulation_EffectStress_RedCCX(RedCCX_temp)
    end if

    ! 6. Soil Profile
    call SetProfFile(ProjectInput(NrRun)%Soil_Filename)
    if (GetProfFile() == '(External)') then
        call SetProfFilefull(GetProfFile())
    elseif (GetProfFile() == '(None)') then
        call SetProfFilefull(GetPathNameSimul() // 'DEFAULT.SOL')
    else
        call SetProfFilefull(ProjectInput(NrRun)%Soil_Directory &
                             // GetProfFile())
    end if

    ! The load of profile is delayed to check if soil water profile need to be
    ! reset (see 8.)

    ! 7. Groundwater
    call SetGroundWaterFile(ProjectInput(NrRun)%GroundWater_Filename)
    if (GetGroundWaterFile() == '(None)') then
        call SetGroundWaterFilefull(GetGroundWaterFile())
        call SetGroundWaterDescription('no shallow groundwater table')
    else
        call SetGroundWaterFilefull(ProjectInput(NrRun)%GroundWater_Directory &
                                    // GetGroundWaterFile())
        ! Loading the groundwater is done after loading the soil profile (see
        ! 9.)
    end if

    ! 8. Set simulation period
    call SetSimulation_FromDayNr(ProjectInput(NrRun)%Simulation_DayNr1)
    call SetSimulation_ToDayNr(ProjectInput(NrRun)%Simulation_DayNrN)
    if ((GetCrop_Day1() /= GetSimulation_FromDayNr()) .or. &
        (GetCrop_DayN() /= GetSimulation_ToDayNr())) then
        call SetSimulation_LinkCropToSimPeriod(.false.)
    end if

    ! 9. Initial conditions
    if (ProjectInput(NrRun)%SWCIni_Filename == 'KeepSWC') then
        ! No load of soil file (which reset thickness compartments and Soil
        ! water content to FC)
        call SetSWCIniFile('KeepSWC')
        call SetSWCIniDescription('Keep soil water profile of previous run')
    else
        ! start with load and complete profile description (see 5.) which reset
        ! SWC to FC by default
        if (GetProfFile() == '(External)') then
            call LoadProfileProcessing(ProjectInput(NrRun)%VersionNr)
        else
            call LoadProfile(GetProfFilefull())
        end if
        call CompleteProfileDescription

        ! Adjust size of compartments if required
        TotDepth = 0._dp
        do i = 1, GetNrCompartments()
            TotDepth = TotDepth + GetCompartment_Thickness(i)
        end do
        if (GetSimulation_MultipleRunWithKeepSWC()) then
        ! Project with a sequence of simulation runs and KeepSWC
            if (roundc(GetSimulation_MultipleRunConstZrx()*1000._dp, mold=1) > &
                roundc(TotDepth*1000._dp, mold=1)) then
                call AdjustSizeCompartments(GetSimulation_MultipleRunConstZrx())
            end if
        else
            if (roundc(GetCrop_RootMax()*1000._dp, mold=1) > &
                roundc(TotDepth*1000._dp, mold=1)) then
                if (roundc(GetSoil_RootMax()*1000._dp, mold=1) == &
                    roundc(GetCrop_RootMax()*1000._dp, mold=1)) then
                    call AdjustSizeCompartments(&
                            real(GetCrop_RootMax(), kind=dp))
                    ! no restrictive soil layer
                else
                    ! restrictive soil layer
                    if (roundc(GetSoil_RootMax()*1000._dp, mold=1) > &
                        roundc(TotDepth*1000._dp, mold=1)) then
                        call AdjustSizeCompartments(&
                            real(GetSoil_RootMax(), kind=dp))
                    end if
                end if
            end if
        end if

        call SetSWCIniFile(ProjectInput(NrRun)%SWCIni_Filename)
        if (GetSWCIniFile() == '(None)') then
            call SetSWCiniFileFull(GetSWCiniFile()) ! no file
            call SetSWCiniDescription(&
                     'Soil water profile at Field Capacity')
        else
            call SetSWCiniFileFull(ProjectInput(NrRun)%SWCIni_Directory &
                                   // GetSWCIniFile())
            SurfaceStorage_temp = GetSurfaceStorage()
            call LoadInitialConditions(GetSWCiniFileFull(),&
                  SurfaceStorage_temp)
            call SetSurfaceStorage(SurfaceStorage_temp)
        end if

        Compartment_temp = GetCompartment()

        select case (GetSimulation_IniSWC_AtDepths())
        case (.true.)
            call TranslateIniPointsToSWProfile(&
               GetSimulation_IniSWC_NrLoc(), &
               GetSimulation_IniSWC_Loc(), GetSimulation_IniSWC_VolProc(), &
               GetSimulation_IniSWC_SaltECe(), GetNrCompartments(), &
               Compartment_temp)
        case default
            call TranslateIniLayersToSWProfile(&
               GetSimulation_IniSWC_NrLoc(),&
               GetSimulation_IniSWC_Loc(), GetSimulation_IniSWC_VolProc(), &
               GetSimulation_IniSWC_SaltECe(), GetNrCompartments(),&
               Compartment_temp)
        end select
        call SetCompartment(Compartment_temp)

        if (GetSimulation_ResetIniSWC()) then
             ! to reset SWC and SALT at end of simulation run
            do i = 1, GetNrCompartments()
                 call SetSimulation_ThetaIni_i(i, GetCompartment_Theta(i))
                 call SetSimulation_ECeIni_i(i, &
                          ECeComp(GetCompartment_i(i)))
            end do
            ! ADDED WHEN DESINGNING 4.0 BECAUSE BELIEVED TO HAVE FORGOTTEN -
            ! CHECK LATER
            if (GetManagement_BundHeight() >= 0.01_dp) then
                 call SetSimulation_SurfaceStorageIni(GetSurfaceStorage())
                 call SetSimulation_ECStorageIni(GetECStorage())
             end if
        end if
    end if

    ! 10. load the groundwater file if it exists (only possible for Version 4.0
    ! and higher)
    if ((roundc(10*ProjectInput(NrRun)%VersionNr, mold=1) >= 40) .and. &
        (GetGroundWaterFile() /= '(None)')) then
          ! the groundwater file is only available in Version 4.0 or higher
        ZiAqua_temp = GetZiAqua()
        ECiAqua_temp = GetECiAqua()
        call LoadGroundWater(GetGroundWaterFilefull(),&
                GetSimulation_FromDayNr(), ZiAqua_temp, ECiAqua_temp)
        call SetZiAqua(ZiAqua_temp)
        call SetECiAqua(ECiAqua_temp)
    else
        call SetZiAqua(undef_int)
        call SetECiAqua(real(undef_int, kind=dp))
        call SetSimulParam_ConstGwt(.true.)
    end if
    Compartment_temp = GetCompartment()
    call CalculateAdjustedFC((GetZiAqua()/100._dp), Compartment_temp)
    call SetCompartment(Compartment_temp)
    if (GetSimulation_IniSWC_AtFC() .and. (GetSWCIniFile() /= 'KeepSWC')) then
        call ResetSWCToFC()
    end if

    ! 11. Off-season conditions
    call SetOffSeasonFile(ProjectInput(NrRun)%OffSeason_Filename)
    if (GetOffSeasonFile() == '(None)') then
        call SetOffSeasonFileFull(GetOffSeasonFile())
        call SetOffSeasonDescription('No specific off-season conditions')
    else
        call SetOffSeasonFileFull(ProjectInput(NrRun)%OffSeason_Directory &
                                  // GetOffSeasonFile())
        call LoadOffSeason(GetOffSeasonFilefull())
    end if

    ! 12. Field data
    call SetObservationsFile(ProjectInput(NrRun)%Observations_Filename)
    if (GetObservationsFile() == '(None)') then
        call SetObservationsFileFull(GetObservationsFile())
        call SetObservationsDescription('No field observations')
    else
        call SetObservationsFileFull(ProjectInput(NrRun)%Observations_Directory &
                                     // GetObservationsFile())
        observations_descr = GetObservationsDescription()
        call GetFileDescription(GetObservationsFileFull(), observations_descr)
        call SetObservationsDescription(observations_descr)
    end if


    contains


    subroutine GetFileDescription(TheFileFullName, TheDescription)
        character(len=*), intent(in) :: TheFileFullName
        character(len=*), intent(inout) :: TheDescription

        integer(int32) :: f0, rc

        open(newunit=f0, file=trim(TheFileFullName), &
                     status='old', action='read', iostat=rc)
        read(f0, '(a)', iostat=rc) TheDescription
        close(f0)
    end subroutine GetFileDescription
end subroutine LoadSimulationRunProject


subroutine BTransferPeriod(TheDaysToCCini, TheGDDaysToCCini,&
              L0, L12, L123, L1234, GDDL0, GDDL12, GDDL123, GDDL1234,&
              CCo, CCx, CGC, GDDCGC, CDC, GDDCDC, KcTop, &
              KcDeclAgeing, CCeffectProcent, WPbio, TheCO2,&
              Tbase, Tupper, TDayMin, TDayMax, GDtranspLow, RatDGDD,&
              TheModeCycle, TempAssimPeriod, TempAssimStored,&
              SumBtot, SumBstored)
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    real(dp), intent(in) :: CCo
    real(dp), intent(in) :: CCx
    real(dp), intent(in) :: CGC
    real(dp), intent(in) :: GDDCGC
    real(dp), intent(in) :: CDC
    real(dp), intent(in) :: GDDCDC
    real(dp), intent(in) :: KcTop
    real(dp), intent(in) :: KcDeclAgeing
    real(dp), intent(in) :: CCeffectProcent
    real(dp), intent(in) :: WPbio
    real(dp), intent(in) :: TheCO2
    real(dp), intent(in) :: Tbase
    real(dp), intent(in) :: Tupper
    real(dp), intent(in) :: TDayMin
    real(dp), intent(in) :: TDayMax
    real(dp), intent(in) :: GDtranspLow
    real(dp), intent(in) :: RatDGDD
    integer(intEnum), intent(in) :: TheModeCycle
    integer(int32), intent(in) :: TempAssimPeriod
    integer(int8), intent(in) :: TempAssimStored
    real(dp), intent(inout) :: SumBtot
    real(dp), intent(inout) :: SumBstored

    real(dp), parameter :: EToStandard = 5._dp

    integer(int32) :: fTemp, rc
    real(dp) :: SumGDDfromDay1, SumGDDforPlot, SumGDD, DayFraction, &
                GDDayFraction, CCinitial, Tndayi, Txdayi, GDDi, CCi, &
                CCxWitheredForB, TpotForB, EpotTotForB
    logical :: GrowthON
    integer(int32) :: GDDTadj, Tadj, DayCC, Dayi, StartStorage

    ! 1. Open Temperature file
    if ((GetTemperatureFile() /= '(None)') .and. &
        (GetTemperatureFile() /= '(External)')) then
        open(newunit=fTemp, file=trim(GetPathNameSimul()//'TCrop.SIM'), &
             status='old', action='read', iostat=rc)
    end if
     ! 2. initialize
    call SetSimulation_DelayedDays(0) ! required for CalculateETpot
    SumBtot = 0._dp
    SumBstored = 0._dp
    SumGDDforPlot = undef_int
    SumGDD = undef_int
    SumGDDfromDay1 = 0._dp
    GrowthON = .false.
    GDDTadj = undef_int
    DayFraction = undef_int
    GDDayFraction = undef_int
    StartStorage = L1234 - TempAssimPeriod + 1
    CCxWitheredForB = 0._dp

    ! 3. Initialise 1st day
    if (TheDaysToCCini /= 0) then
       ! regrowth which starts on 1st day
        GrowthON = .true.
        if (TheDaysToCCini == undef_int) then
            ! CCx on 1st day
            Tadj = L12 - L0
            if (TheModeCycle == modeCycle_GDDays) then
                GDDTadj = GDDL12 - GDDL0
                SumGDD = GDDL12
            end if
            CCinitial = CCx
        else
            ! CC on 1st day is < CCx
            Tadj = TheDaysToCCini
            DayCC = Tadj + L0
            if (TheModeCycle == modeCycle_GDDays) then
                GDDTadj = TheGDDaysToCCini
                SumGDD = GDDL0 + TheGDDaysToCCini
                SumGDDforPlot = SumGDD
            end if
            CCinitial = CanopyCoverNoStressSF(DayCC, L0, L123, L1234,&
                GDDL0, GDDL123, GDDL1234, CCo, CCx, CGC, CDC,&
                GDDCGC, GDDCDC, SumGDDforPlot, TheModeCycle, 0_int8, 0_int8)
        end if
        ! Time reduction for days between L12 and L123
        DayFraction = (L123-L12) *1._dp/ &
                      real(Tadj + L0 + (L123-L12), kind=dp)
        if (TheModeCycle == modeCycle_GDDays) then
            GDDayFraction = (GDDL123-GDDL12) *1._dp/&
                            real(GDDTadj + GDDL0 + (GDDL123-GDDL12), kind=dp)
        end if
    else
        ! growth starts after germination/recover
        Tadj = 0
        if (TheModeCycle == modeCycle_GDDays) then
            GDDTadj = 0._dp
            SumGDD = 0._dp
        end if
        CCinitial = CCo
    end if

    ! 4. Calculate Biomass
    do Dayi = 1, L1234
        ! 4.1 growing degrees for dayi
        if (GetTemperatureFile() == '(None)') then
            GDDi = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                              GetSimulParam_GDDMethod())
        elseif (GetTemperatureFile() == '(External)') then 
            Tndayi = real(GetTminRun_i(GetCrop_Day1()-GetSimulation_FromDayNr()+Dayi),kind=dp)
            Txdayi = real(GetTmaxRun_i(GetCrop_Day1()-GetSimulation_FromDayNr()+Dayi),kind=dp)
            GDDi = DegreesDay(Tbase, Tupper, Tndayi, Txdayi, &
                                    GetSimulParam_GDDMethod())
        else
            read(fTemp, *, iostat=rc) Tndayi, Txdayi
            GDDi = DegreesDay(Tbase, Tupper, Tndayi, Txdayi, &
                              GetSimulParam_GDDMethod())
        end if
        if (TheModeCycle == modeCycle_GDDays) then
            SumGDD = SumGDD + GDDi
            SumGDDfromDay1 = SumGDDfromDay1 + GDDi
        end if

        ! 4.2 green Canopy Cover (CC)
        DayCC = Dayi
        if (GrowthON .eqv. .false.) then
            ! not yet canopy development
            CCi = 0._dp
            if (TheDaysToCCini /= 0) then
                ! regrowth
                CCi = CCinitial
                GrowthON = .true.
            else
                ! sowing or transplanting
                if (TheModeCycle == modeCycle_CalendarDays) then
                    if (Dayi == (L0+1)) then
                        CCi = CCinitial
                        GrowthON = .true.
                    end if
                else
                    if (SumGDD > GDDL0) then
                        CCi = CCinitial
                        GrowthON = .true.
                    end if
                end if
            end if
        else
            if (TheDaysToCCini == 0) then
                DayCC = Dayi
            else
                DayCC = Dayi + Tadj + L0 ! adjusted time scale
                if (DayCC > L1234) then
                    DayCC = L1234 ! special case where L123 > L1234
                end if
                if (DayCC > L12) then
                    if (Dayi <= L123) then
                        DayCC = L12 + roundc(DayFraction *&
                             real(Dayi+Tadj+L0 - L12, kind=dp),mold=1) ! slow down
                    else
                        DayCC = Dayi ! switch time scale
                    end if
                end if
            end if
            if (TheModeCycle == modeCycle_GDDays) then
                if (TheGDDaysToCCini == 0) then
                    SumGDDforPlot = SumGDDfromDay1
                else
                    SumGDDforPlot = SumGDD
                    if (SumGDDforPlot > GDDL1234) then
                        SumGDDforPlot = GDDL1234 ! special case where L123 > L1234
                    end if
                    if (SumGDDforPlot > GDDL12) then
                        if (SumGDDfromDay1 <= GDDL123) then
                            SumGDDforPlot = GDDL12 + real(GDDayFraction * &
                              real(SumGDDfromDay1+GDDTadj+GDDL0 - GDDL12,&
                                   kind=dp)) ! slow down
                        else
                            SumGDDforPlot = SumGDDfromDay1 ! switch time scale
                        end if
                    end if
                    CCi = CCiNoWaterStressSF(DayCC, L0, L12, L123, L1234,&
                        GDDL0, GDDL12, GDDL123, GDDL1234,&
                        CCo, CCx, CGC, GDDCGC, CDC, GDDCDC, SumGDDforPlot,&
                        RatDGDD, 0_int8, 0_int8, 0._dp, TheModeCycle)
                end if
                if (CCi > CCxWitheredForB) then
                     CCxWitheredForB = CCi
                end if

                ! 4.3 potential transpiration (TpotForB)
                if (CCi > 0.0001_dp) then
                    ! 5.3 potential transpiration of total canopy cover
                    call CalculateETpot(DayCC, L0, L12, L123, L1234, (0), CCi,&
                         EToStandard, KcTop, KcDeclAgeing,&
                         CCx, CCxWitheredForB, CCeffectProcent, TheCO2, GDDi, &
                         GDtranspLow, TpotForB, EpotTotForB)
                else
                    TpotForB = 0._dp
                end if

                ! 4.4 Biomass (B)
                if (Dayi >= StartStorage) then
                    SumBtot = SumBtot +  WPbio * (TpotForB/EToStandard)
                    SumBstored = SumBstored + WPbio*(TpotForB/EToStandard)*&
                            (0.01_dp*TempAssimStored)*&
                            (1-KsAny(((Dayi-StartStorage+1._dp)/&
                               real(TempAssimPeriod, kind=dp)),&
                               0._dp,1._dp,-5._dp));
               end if
           end if

           ! 5. Close Temperature file
           if ((GetTemperatureFile() /= '(None)') .and. &
               (GetTemperatureFile() /= '(External)')) then
               close(fTemp)
           end if
       end if
    end do
end subroutine BTransferPeriod


real(dp) function Bnormalized(TheDaysToCCini, TheGDDaysToCCini,&
            L0, L12, L12SF, L123, L1234, LFlor, &
            GDDL0, GDDL12, GDDL12SF, GDDL123, GDDL1234, &
            WPyield, DaysYieldFormation, tSwitch, CCo, CCx, &
            CGC, GDDCGC, CDC, GDDCDC, KcTop, KcDeclAgeing, &
            CCeffectProcent, WPbio, TheCO2, Tbase, Tupper, &
            TDayMin, TDayMax, GDtranspLow, RatDGDD, SumKcTop, &
            StressInPercent, StrResRedCGC, StrResRedCCx, StrResRedWP, &
            StrResRedKsSto, WeedStress, DeltaWeedStress, StrResCDecline, &
            ShapeFweed, TheModeCycle, FertilityStressOn, ReferenceClimate)
     integer(int32), intent(in) :: TheDaysToCCini
     integer(int32), intent(in) :: TheGDDaysToCCini
     integer(int32), intent(in) :: L0
     integer(int32), intent(in) :: L12
     integer(int32), intent(in) :: L12SF
     integer(int32), intent(in) :: L123
     integer(int32), intent(in) :: L1234
     integer(int32), intent(in) :: LFlor
     integer(int32), intent(in) :: GDDL0
     integer(int32), intent(in) :: GDDL12
     integer(int32), intent(in) :: GDDL12SF
     integer(int32), intent(in) :: GDDL123
     integer(int32), intent(in) :: GDDL1234
     integer(int32), intent(in) :: WPyield
     integer(int32), intent(in) :: DaysYieldFormation
     integer(int32), intent(in) :: tSwitch
     real(dp), intent(in) :: CCo
     real(dp), intent(in) :: CCx
     real(dp), intent(in) :: CGC
     real(dp), intent(in) :: GDDCGC
     real(dp), intent(in) :: CDC
     real(dp), intent(in) :: GDDCDC
     real(dp), intent(in) :: KcTop
     real(dp), intent(in) :: KcDeclAgeing
     real(dp), intent(in) :: CCeffectProcent
     real(dp), intent(in) :: WPbio
     real(dp), intent(in) :: TheCO2
     real(dp), intent(in) :: Tbase
     real(dp), intent(in) :: Tupper
     real(dp), intent(in) :: TDayMin
     real(dp), intent(in) :: TDayMax
     real(dp), intent(in) :: GDtranspLow
     real(dp), intent(in) :: RatDGDD
     real(dp), intent(in) :: SumKcTop
     integer(int32), intent(in) :: StressInPercent
     integer(int8), intent(in) :: StrResRedCGC
     integer(int8), intent(in) :: StrResRedCCx
     integer(int8), intent(in) :: StrResRedWP
     integer(int8), intent(in) :: StrResRedKsSto
     integer(int8), intent(in) :: WeedStress
     integer(int32), intent(in) :: DeltaWeedStress
     real(dp), intent(in) :: StrResCDecline
     real(dp), intent(in) :: ShapeFweed
     integer(intEnum), intent(in) :: TheModeCycle
     logical, intent(in) :: FertilityStressOn
     logical, intent(in) :: ReferenceClimate

     real(dp), parameter :: EToStandard = 5._dp
     integer(int32), parameter :: k = 2

     integer(int32) ::  fTemp, rc
     real(dp) :: SumGDD, Tndayi, Txdayi, GDDi, CCi,&
                 CCxWitheredForB, TpotForB, EpotTotForB, SumKCi,&
                 fSwitch, WPi, SumBnor, SumKcTopSF, fCCx
     integer(int32) :: Dayi, DayCC, Tadj, GDDTadj
     real(dp) :: CCoadj, CCxadj, CDCadj, GDDCDCadj, CCw, CCtotStar, CCwStar
     real(dp) :: SumGDDfromDay1, SumGDDforPlot, CCinitial,&
                 DayFraction, GDDayFraction, fWeed, WeedCorrection
     logical :: GrowthON
     integer(int32) :: DeltaWeedStress_local

     ! 1. Adjustment for weed infestation
     if (WeedStress > 0) then
         if (StressInPercent > 0) then ! soil fertility stress
             fWeed = 1._dp  ! no expansion of canopy cover possible
         else
             fWeed = CCmultiplierWeed(WeedStress, CCx, ShapeFweed)
         end if
         CCoadj = CCo*fWeed
         CCxadj = CCx*fWeed
         CDCadj = CDC*(fWeed*CCx + 2.29_dp)/(CCx + 2.29_dp)
         GDDCDCadj = GDDCDC*(fWeed*CCx + 2.29_dp)/(CCx + 2.29_dp)
     else
         CCoadj = CCo
         CCxadj = CCx
         CDCadj = CDC
         GDDCDCadj = GDDCDC
     end if

     ! 2. Open Temperature file
     if ((GetTemperatureFile() /= '(None)') .and. &
         (GetTemperatureFile() /= '(External)')) then
        if (ReferenceClimate .eqv. .true.) then
            open(newunit=fTemp, file=trim(GetPathNameSimul()//'TCropReference.SIM'), &
                 status='old', action='read', iostat=rc)
        else
            open(newunit=fTemp, file=trim(GetPathNameSimul()//'TCrop.SIM'), &
                 status='old', action='read', iostat=rc)
        end if
     end if

     ! 3. Initialize
     SumKcTopSF = (1._dp - real(StressInPercent, kind=dp)/100._dp) * SumKcTop
     !! only required for soil fertility stress

     call SetSimulation_DelayedDays(0) ! required for CalculateETpot
     SumKci = 0._dp
     SumBnor = 0._dp
     SumGDDforPlot = undef_int
     SumGDD = undef_int
     SumGDDfromDay1 = 0._dp
     GrowthON = .false.
     GDDTadj = undef_int
     DayFraction = undef_int
     GDDayFraction = undef_int
     CCxWitheredForB = 0._dp

     ! 4. Initialise 1st day
     if (TheDaysToCCini /= 0) then
         ! regrowth which starts on 1st day
         GrowthON = .true.
         if (TheDaysToCCini == undef_int) then
             ! CCx on 1st day
             Tadj = L12 - L0
             if (TheModeCycle == modeCycle_GDDays) then
                 GDDTadj = GDDL12 - GDDL0
                 SumGDD = GDDL12
             end if
             CCinitial = CCxadj * (1._dp-StrResRedCCX/100._dp)
         else
         ! CC on 1st day is < CCx
             Tadj = TheDaysToCCini
             DayCC = Tadj + L0
             if (TheModeCycle == modeCycle_GDDays) then
                 GDDTadj = TheGDDaysToCCini
                 SumGDD = GDDL0 + TheGDDaysToCCini
                 SumGDDforPlot = SumGDD
             end if
             CCinitial = CanopyCoverNoStressSF(DayCC, L0, L123, L1234,&
                 GDDL0, GDDL123, GDDL1234, CCoadj, CCxadj, CGC, CDCadj,&
                 GDDCGC, GDDCDCadj, SumGDDforPlot, TheModeCycle, &
                 StrResRedCGC, StrResRedCCX)
         end if
         ! Time reduction for days between L12 and L123
         DayFraction = (L123-L12) * 1._dp/&
                       real(Tadj + L0 + (L123-L12) ,kind=dp)
         if (TheModeCycle == modeCycle_GDDays) then
             GDDayFraction = (GDDL123-GDDL12) * 1._dp/&
                             (GDDTadj + GDDL0 + (GDDL123-GDDL12))
         end if
     else
         ! growth starts after germination/recover
         Tadj = 0
         if (TheModeCycle == modeCycle_GDDays) then
             GDDTadj = 0
             SumGDD = 0._dp
         end if
         CCinitial = CCoadj
     end if

     ! 5. Calculate Bnormalized
     do Dayi = 1, L1234
         ! 5.1 growing degrees for dayi
         if (GetTemperatureFile() == '(None)') then
             GDDi = DegreesDay(Tbase, Tupper, TDayMin, TDayMax,&
                               GetSimulParam_GDDMethod())
         elseif (GetTemperatureFile() == '(External)') then 
             Tndayi = real(GetTminRun_i(GetCrop_Day1()-GetSimulation_FromDayNr()+Dayi),kind=dp)
             Txdayi = real(GetTmaxRun_i(GetCrop_Day1()-GetSimulation_FromDayNr()+Dayi),kind=dp)
             GDDi = DegreesDay(Tbase, Tupper, Tndayi, Txdayi, &
                                    GetSimulParam_GDDMethod())
         else
            read(fTemp, *, iostat=rc) Tndayi, Txdayi
            if ((rc == iostat_end) .and. (ReferenceClimate .eqv. .true.)) then
                rewind(fTemp)
                read(fTemp, *, iostat=rc) Tndayi, Txdayi
            end if
            GDDi = DegreesDay(Tbase, Tupper, Tndayi, Txdayi,&
                              GetSimulParam_GDDMethod())
         end if
         if (TheModeCycle == modeCycle_GDDays) then
             SumGDD = SumGDD + GDDi
             SumGDDfromDay1 = SumGDDfromDay1 + GDDi
         end if

         ! 5.2 green Canopy Cover (CC)
         DayCC = Dayi
         if (GrowthON .eqv. .false.) then
             ! not yet canopy development
             CCi = 0._dp
             if (TheDaysToCCini /= 0) then
                 ! regrowth
                 CCi = CCinitial
                 GrowthON = .true.
             else
                 ! sowing or transplanting
                 if (TheModeCycle == modeCycle_CalendarDays) then
                     if (Dayi == (L0+1)) then
                         CCi = CCinitial
                         GrowthON = .true.
                     end if
                 else
                     if (SumGDD > GDDL0) then
                         CCi = CCinitial
                         GrowthON = .true.
                     end if
                 end if
             end if
         else
             if (TheDaysToCCini == 0) then
                 DayCC = Dayi
             else
                 DayCC = Dayi + Tadj + L0 ! adjusted time scale
                 if (DayCC > L1234) then
                     DayCC = L1234 ! special case where L123 > L1234
                 end if
                 if (DayCC > L12) then
                     if (Dayi <= L123) then
                          DayCC = L12 + roundc(DayFraction * &
                                     (Dayi+Tadj+L0 - L12), mold=1) ! slow down
                     else
                         DayCC = Dayi ! switch time scale
                     end if
                 end if
             end if

             if (TheModeCycle == modeCycle_GDDays) then
                 if (TheGDDaysToCCini == 0) then
                     SumGDDforPlot = SumGDDfromDay1
                 else
                     SumGDDforPlot = SumGDD
                     if (SumGDDforPlot > GDDL1234) then
                         SumGDDforPlot = GDDL1234
                         ! special case where L123 > L1234
                     end if
                     if (SumGDDforPlot > GDDL12) then
                         if (SumGDDfromDay1 <= GDDL123) then
                             SumGDDforPlot = GDDL12 + roundc(GDDayFraction * &
                                 (SumGDDfromDay1+GDDTadj+GDDL0 - &
                                 GDDL12),mold=1) ! slow down
                         else
                             SumGDDforPlot = SumGDDfromDay1 ! switch time scale
                         end if
                     end if
                 end if
             endif
             CCi = CCiNoWaterStressSF(DayCC, L0, L12SF, L123, L1234,&
                         GDDL0, GDDL12SF, GDDL123, GDDL1234,&
                         CCoadj, CCxadj, CGC, GDDCGC, CDCadj, GDDCDCadj, &
                         SumGDDforPlot, RatDGDD,&
                         StrResRedCGC, StrResRedCCX, StrResCDecline,&
                         TheModeCycle)
         end if

         if (CCi > CCxWitheredForB) then
             CCxWitheredForB = CCi
         end if
         if (DayCC >= L12SF) then
             CCxWitheredForB = CCxadj*(1._dp-StrResRedCCX/100._dp)
         end if
         CCw = CCi

         if (CCi > 0.0001_dp) then
             ! 5.3 potential transpiration of total canopy cover (crop and weed)
             call CalculateETpot(DayCC, L0, L12, L123, L1234, (0), CCi, &
                            EToStandard, KcTop, KcDeclAgeing,&
                            CCxadj, CCxWitheredForB, CCeffectProcent, TheCO2,&
                            GDDi, GDtranspLow, TpotForB, EpotTotForB)

             ! 5.4 Sum of Kc (only required for soil fertility stress)
             SumKci = SumKci + (TpotForB/EToStandard)

             ! 5.5 potential transpiration of crop canopy cover (without weed)
             if (WeedStress > 0) then
                 ! green canopy cover of the crop (CCw) in weed-infested field
                 ! (CCi is CC of crop and weeds)
                 fCCx = 1.0_dp ! only for non perennials (no self-thinning)
                 if (DeltaWeedStress /= 0) then
                     DeltaWeedStress_local = DeltaWeedStress
                     WeedCorrection = GetWeedRC(DayCC, SumGDDforPlot, fCCx,&
                            WeedStress, GetManagement_WeedAdj(),&
                            DeltaWeedStress_local, L12SF, L123, &
                            GDDL12SF, GDDL123, TheModeCycle)
                 else
                     WeedCorrection = WeedStress
                 end if
                 CCw = CCi * (1._dp - WeedCorrection/100._dp)
                 ! correction for micro-advection
                 CCtotStar = 1.72_dp*CCi - 1._dp*(CCi*CCi) + &
                                 0.30_dp*(CCi*CCi*CCi)
                 if (CCtotStar < 0._dp) then
                     CCtotStar = 0._dp
                 end if
                 if (CCtotStar > 1._dp) then
                     CCtotStar = 1._dp
                 end if
                 if (CCw > 0.0001_dp) then
                     CCwStar = CCw + (CCtotStar - CCi)
                 else
                     CCwStar = 0._dp
                 end if
                 ! crop transpiration in weed-infested field
                 if (CCtotStar <= 0.0001_dp) then
                     TpotForB = 0._dp
                 else
                     TpotForB = TpotForB * (CCwStar/CCtotStar)
                 end if
             end if
         else
             TpotForB = 0._dp
         end if

         ! 5.6 biomass water productivity (WP)
         WPi = WPbio ! vegetative stage
         ! 5.6a. vegetative versus yield formation stage
         if (((GetCrop_subkind() == subkind_Tuber) &
             .or. (GetCrop_subkind() == subkind_Grain)) .and.&
             (WPyield < 100) .and. (Dayi > LFlor)) then
             ! yield formation stage
             fSwitch = 1._dp
             if ((DaysYieldFormation > 0) .and. (tSwitch > 0)) then
                 fSwitch = (Dayi-LFlor) * 1._dp/real(tSwitch, kind=dp)
                 if (fSwitch > 1) then
                     fSwitch = 1._dp
                 end if
             end if
             WPi = WPi * (1._dp - (1._dp - WPyield/100._dp)*fSwitch)
         end if

         ! 5.7 Biomass (B)
         if (FertilityStressOn) then
             ! 5.7a - reduction for soil fertiltiy
             if ((StrResRedWP > 0) .and. (SumKci > 0._dp) &
                 .and. (SumKcTopSF > epsilon(1.0))) then
                 if (SumKci < SumKcTopSF) then
                     if (SumKci > 0) then
                         WPi = WPi * (1._dp - (StrResRedWP/100._dp) *&
                                 exp(k*log(SumKci/SumKcTopSF)))
                     end if
                 else
                     WPi = WPi * (1._dp - StrResRedWP/100._dp)
                 end if
             end if
             ! 5.7b - Biomass (B)
             SumBnor = SumBnor +  WPi * (TpotForB/EToStandard)
         else
             SumBnor = SumBnor +  WPi * (1._dp - StrResRedKsSto/100._dp) *&
                           (TpotForB/EToStandard) ! for salinity stress
         end if

     enddo

     ! 4. Close Temperature file
     if ((GetTemperatureFile() /= '(None)') .and. &
         (GetTemperatureFile() /= '(External)')) then
         close(fTemp)
     end if

     ! 5. Export
     Bnormalized = SumBnor
end function Bnormalized


subroutine CreateTnxReference365Days()

    character(len=:), allocatable :: totalnameOUT
    character(len=1025) :: TempString
    integer(int32) :: DayNri
    integer(int32) :: i, Monthi
    type(rep_DayEventDbl), dimension(31) :: TminDataSet_temp,TmaxDataSet_temp
    real(sp) :: Tlow, Thigh

    if ((FileExists(GetTnxReferenceFileFull())) .OR. (GetTnxReferenceFile() == '(External)')) then
        ! create SIM file
        if (GetTnxReferenceFile() /= '(External)') then
            totalnameOUT = trim(GetPathNameSimul()) // 'TnxReference365Days.SIM'
            call fTnxReference365Days_open(totalnameOUT, 'w')
        end if
        
        ! get data set for 1st month
        Monthi = 1
        call GetMonthlyTemperatureDataSetFromTnxReferenceFile(Monthi, TminDataSet_temp, TmaxDataSet_temp)
        
        ! Tnx for Day 1 to 365
        do DayNri = 1, 365
            if (DayNri > TminDataSet_temp(31)%DayNr) then
                ! next month
                Monthi = Monthi + 1
                call GetMonthlyTemperatureDataSetFromTnxReferenceFile(Monthi, TminDataSet_temp, TmaxDataSet_temp)
            end if
            i = 1
            do while (TminDataSet_temp(i)%DayNr /= DayNri)
                i = i+1
            end do
            Tlow = real(roundc(100*TminDataSet_temp(i)%Param, mold=int32),kind=sp)/100._sp
            Thigh = real(roundc(100*TmaxDataSet_temp(i)%Param, mold=int32),kind=sp)/100._sp
            call SetTminTnxReference365DaysRun_i(DayNri,Tlow)
            call SetTmaxTnxReference365DaysRun_i(DayNri,Thigh)
            if (GetTnxReferenceFile() /= '(External)') then
                write(TempString, '(2f10.2)') Tlow, Thigh
                call fTnxReference365Days_write(trim(TempString))
            end if
        end do
        if (GetTnxReferenceFile() /= '(External)') then
            ! Close file
            call fTnxReference365Days_close()
        end if
    end if
end subroutine CreateTnxReference365Days


subroutine CreateTnxReferenceFile(TemperatureFile, TnxReferenceFile, TnxReferenceYear)
    character(len=*), intent(in) :: TemperatureFile
    character(len=*), intent(in) :: TnxReferenceFile
    integer(int32), intent(in) :: TnxReferenceYear

    real(dp), dimension(12) :: MonthVal1, MonthVal2
    character(len=:), allocatable :: FullName
    integer(int32)   :: fhandle, rc
    integer(int32) :: i, NrYears
    integer(int32) :: Yeari, Monthi, MonthDays, MonthDecs, Deci, Dayi
    real(dp) :: SUM1, SUM2, Val1, Val2
    integer(int32) :: DayNri, EndDayNr
    logical :: EndMonth
    integer(int32), dimension(12) :: MonthNrYears
    character(len=:), allocatable :: TempString
    character(len=1025) :: TempString2

    ! 1a. Delete existing TnxReferenceFile and adjust TnxReferenceYear
    FullName = trim(GetPathNameSimul()) // TnxReferenceFile
    if (FileExists(FullName)) then
        call unlink(FullName)
        call SetTnxReferenceYear(2000)
    end if
    ! 1b. Delete existing TnxReference365Days.SIM
    FullName = trim(GetPathNameSimul()) // 'TnxReference365Days.SIM'
    if (FileExists(FullName)) then
        call unlink(FullName)
    end if
    ! 1c. Delete existing TCropReference.SIM
    FullName = trim(GetPathNameSimul()) // 'TCropReference.SIM'
    if (FileExists(FullName)) then
        call unlink(FullName)
    end if


    ! 2. Get Mean monthly data
    if (TemperatureFile /= '(None)') then
        ! 2.a Preparation
        ! Get Number of Years
        NrYears = GetTemperatureRecord_ToY() - GetTemperatureRecord_FromY() + 1
        ! Get TnxReferenceYear
        call SetTnxReferenceYear(roundc((GetTemperatureRecord_FromY()+GetTemperatureRecord_ToY())/2._dp,mold=1_int32))
        if (GetTnxReferenceYear() == 1901) then
            call SetTnxReferenceYear(2000)
        end if
        ! check number of years for each month
        do Monthi = 1, 12
            MonthNrYears(Monthi) = NrYears
        end do
        if (NrYears > 1) then
            if (GetTemperatureRecord_FromM() > 1) then
                do Monthi = 1, (GetTemperatureRecord_FromM()-1)
                    MonthNrYears(Monthi) = MonthNrYears(Monthi) - 1
                end do
            end if
            if (GetTemperatureRecord_ToM() < 12) then
                do Monthi = 12, (GetTemperatureRecord_ToM()+1), -1
                    MonthNrYears(Monthi) = MonthNrYears(Monthi) - 1
                end do
            end if
        end if
                
        ! initialize
        do Monthi = 1, 12
            MonthVal1(Monthi) = 0._dp ! Tmin data
            MonthVal2(Monthi) = 0._dp ! Tmax data
        end do
        Yeari = GetTemperatureRecord_FromY()
        Monthi = GetTemperatureRecord_FromM()
        DayNri = GetTemperatureRecord_FromDayNr()
        EndMonth = .false.
        EndDayNr = undef_int
        Deci = undef_int
        SUM1 = 0._dp
        SUM2 = 0._dp
        MonthDays = 0
        MonthDecs = 0
        select case (GetTemperatureRecord_DataType())
        case (datatype_daily)
            ! EndDayNr = DayNr of last day in month
            Dayi = DaysInMonth(Monthi)
            if ((LeapYear(Yeari) .eqv. .true.) .and. (Monthi == 2)) then
                Dayi = Dayi + 1
            end if
            call DetermineDayNr(Dayi, Monthi, Yeari, EndDayNr)
        end select
            
        ! open temperature file
        open(newunit=fhandle, file=trim(GetTemperatureFilefull()), &
                     status='old', action='read', iostat=rc)
        read(fhandle, *, iostat=rc) ! Description
        read(fhandle, *, iostat=rc) ! Data type
        read(fhandle, *, iostat=rc) ! Day
        read(fhandle, *, iostat=rc) ! Month
        read(fhandle, *, iostat=rc) ! Year
        read(fhandle, *, iostat=rc) ! Title
        read(fhandle, *, iostat=rc) ! Title
        read(fhandle, *, iostat=rc) ! Title

        ! 2.b Determine mean monthly data
        do while (rc /= iostat_end)
            ! read data
            read(fhandle, *, iostat=rc) Val1, Val2 ! Tmin, Tmax
            if (rc == iostat_end) exit
            SUM1 = SUM1 + Val1 ! minimum temperature
            SUM2 = SUM2 + Val2 ! maximum temperature
            ! check end of month
            select case (GetTemperatureRecord_DataType())
            case (datatype_daily)
                MonthDays = MonthDays + 1
                if (DayNri == EndDayNr) then
                    EndMonth = .true.
                end if
            case (datatype_decadely)
                MonthDecs = MonthDecs + 1
                if (Deci == 3) then
                    EndMonth = .true.
                end if
            case (datatype_monthly)
                EndMonth = .true.
            end select
            if (rc == iostat_end) then
                EndMonth = .true.
            end if
            ! End of Month
            if (EndMonth .eqv. .true.) then
                ! mean monthly values Tmin and Tmax
                select case (GetTemperatureRecord_DataType())
                case (datatype_daily)
                    SUM1 = SUM1/real(MonthDays, kind=dp)
                    SUM2 = SUM2/real(MonthDays, kind=dp)
                case (datatype_decadely)
                    SUM1 = SUM1/real(MonthDecs, kind=dp)
                    SUM2 = SUM2/real(MonthDecs, kind=dp)
                end select
                MonthVal1(Monthi) = MonthVal1(Monthi) + SUM1/real(MonthNrYears(Monthi),kind=dp) ! Tmin
                MonthVal2(Monthi) = MonthVal2(Monthi) + SUM2/real(MonthNrYears(Monthi),kind=dp) ! Tmax
                ! Next month
                SUM1 = 0._dp
                SUM2 = 0._dp
                if (Monthi == 12) then
                    Monthi = 1
                    Yeari = Yeari + 1
                else
                    Monthi = Monthi + 1
                end if
                EndMonth = .false.
                select case (GetTemperatureRecord_DataType())
                case (datatype_daily)
                    ! EndDayNr = DayNr of last day in month
                    MonthDays = 0
                    Dayi = DaysInMonth(Monthi)
                    if ((LeapYear(Yeari) .eqv. .true.) .and. (Monthi == 2)) then
                        Dayi = Dayi + 1
                    end if
                    call DetermineDayNr(Dayi, Monthi, Yeari, EndDayNr)
                case (datatype_decadely)
                    MonthDecs = 0
                    Deci = 0
                end select
            end if
            ! Next day, decade, month
            select case (GetTemperatureRecord_DataType())
            case (datatype_daily)
                DayNri = DayNri + 1
            end select
        end do
        
        ! close temperature file
        close(fhandle)
    end if
    ! 3. Create and Save TnxReference File
    if (TemperatureFile /= '(None)') then
        ! Determine Name of TnxReference File
        i = len(trim(TemperatureFile))
        TempString = trim(TemperatureFile)
        call SetTnxReferenceFile(TempString(:i-4) // 'Reference.Tnx')
        FullName = trim(GetPathNameSimul()) // GetTnxReferenceFile()
        call SetTnxReferenceFileFull(FullName)
        ! Save mean monhtly data data
        call fTnxReference_open(FullName, 'w')
        call fTnxReference_write('Reference : ' // trim(GetTemperatureDescription()))
        call fTnxReference_write('     3  : Monthly records (1=daily, 2=10-daily and 3=monthly data)')
        call fTnxReference_write('     1  : First day of record (1, 11 or 21 for 10-day or 1 for months)')
        call fTnxReference_write('     1  : First month of record')
        call fTnxReference_write('  1901  : First year of record (1901 if not linked to a specific year)')
        call fTnxReference_write('')
        call fTnxReference_write('  Tmin (C)   TMax (C) ')
        call fTnxReference_write('  =======================')
        do Monthi = 1, 12
            write(TempString2, '(2f10.2)') MonthVal1(Monthi), MonthVal2(Monthi)
            call fTnxReference_write(trim(TempString2))
            call SetTminTnxReference12MonthsRun_i(Monthi,real(roundc(100*MonthVal1(Monthi),mold=int32),kind=sp)/100._sp)
            call SetTmaxTnxReference12MonthsRun_i(Monthi,real(roundc(100*MonthVal2(Monthi),mold=int32),kind=sp)/100._sp)
        end do
        call fTnxReference_close()
    end if
end subroutine CreateTnxReferenceFile

!TminDatSet

function GetTminDataSet() result(TminDataSet_out)
    !! Getter for the "TminDataSet" global variable.
    type(rep_DayEventDbl), dimension(31) :: TminDataSet_out

    TminDataSet_out = TminDataSet
end function GetTminDataSet


function GetTminDataSet_i(i) result(TminDataSet_i)
    !! Getter for individual elements of the "TminDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl) :: TminDataSet_i

    TminDataSet_i = TminDataSet(i)
end function GetTminDataSet_i


integer(int32) function GetTminDataSet_DayNr(i)
    integer(int32), intent(in) :: i

    GetTminDataSet_DayNr = TminDataSet(i)%DayNr
end function GetTminDataSet_DayNr


real(dp) function GetTminDataSet_Param(i)
    integer(int32), intent(in) :: i

    GetTminDataSet_Param = TminDataSet(i)%Param
end function GetTminDataSet_Param


subroutine SetTminDataSet(TminDataSet_in)
    !! Setter for the "TminDatSet" global variable.
    type(rep_DayEventDbl), dimension(31), intent(in) :: TminDataSet_in

    TminDataSet = TminDataSet_in
end subroutine SetTminDataSet


subroutine SetTminDataSet_i(i, TminDataSet_i)
    !! Setter for individual element for the "TminDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl), intent(in) :: TminDataSet_i

    TminDataSet(i) = TminDataSet_i
end subroutine SetTminDataSet_i


subroutine SetTminDataSet_DayNr(i, DayNr_in)
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr_in

    TminDataSet(i)%DayNr = DayNr_in
end subroutine SetTminDataSet_DayNr


subroutine SetTminDataSet_Param(i, Param_in)
    integer(int32), intent(in) :: i
    real(dp), intent(in) :: Param_in

    TminDataSet(i)%Param = Param_in
end subroutine SetTminDataSet_Param

! TmaxDataSet

function GetTmaxDataSet() result(TmaxDataSet_out)
    !! Getter for the "TmaxDataSet" global variable.
    type(rep_DayEventDbl), dimension(31) :: TmaxDataSet_out

    TmaxDataSet_out = TmaxDataSet
end function GetTmaxDataSet


function GetTmaxDataSet_i(i) result(TmaxDataSet_i)
    !! Getter for individual elements of the "TmaxDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl) :: TmaxDataSet_i

    TmaxDataSet_i = TmaxDataSet(i)
end function GetTmaxDataSet_i


integer(int32) function GetTmaxDataSet_DayNr(i)
    integer(int32), intent(in) :: i

    GetTmaxDataSet_DayNr = TmaxDataSet(i)%DayNr
end function GetTmaxDataSet_DayNr


real(dp) function GetTmaxDataSet_Param(i)
    integer(int32), intent(in) :: i

    GetTmaxDataSet_Param = TmaxDataSet(i)%Param
end function GetTmaxDataSet_Param


subroutine SetTmaxDataSet(TmaxDataSet_in)
    !! Setter for the "TmaxDatSet" global variable.
    type(rep_DayEventDbl), dimension(31), intent(in) :: TmaxDataSet_in

    TmaxDataSet = TmaxDataSet_in
end subroutine SetTmaxDataSet


subroutine SetTmaxDataSet_i(i, TmaxDataSet_i)
    !! Setter for individual element for the "TmaxDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl), intent(in) :: TmaxDataSet_i

    TmaxDataSet(i) = TmaxDataSet_i
end subroutine SetTmaxDataSet_i


subroutine SetTmaxDataSet_DayNr(i, DayNr_in)
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr_in

    TmaxDataSet(i)%DayNr = DayNr_in
end subroutine SetTmaxDataSet_DayNr


subroutine SetTmaxDataSet_Param(i, Param_in)
    integer(int32), intent(in) :: i
    real(dp), intent(in) :: Param_in

    TmaxDataSet(i)%Param = Param_in
end subroutine SetTmaxDataSet_Param


! fTnxReference

subroutine fTnxReference_open(filename, mode)
    !! Opens the given file, assigning it to the 'fTnxReference' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fTnxReference, filename, mode, fTnxReference_iostat)
end subroutine fTnxReference_open


subroutine fTnxReference_write(line, advance_in)
    !! Writes the given line to the fTnxReference file.
    character(len=*), intent(in) :: line
        !! line to write
    logical, intent(in), optional :: advance_in
        !! whether or not to append a newline character

    logical :: advance

    if (present(advance_in)) then
        advance = advance_in
    else
        advance = .true.
    end if
    call write_file(fTnxReference, line, advance, fTnxReference_iostat)
end subroutine fTnxReference_write


subroutine fTnxReference_close()
    close(fTnxReference)
end subroutine fTnxReference_close


subroutine fTnxReference_erase()
    call unlink(GetTnxReferenceFileFull())
end subroutine fTnxReference_erase


! fTnxReference365Days

subroutine fTnxReference365Days_open(filename, mode)
    !! Opens the given file, assigning it to the 'fTnxReference365Days' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fTnxReference365Days, filename, mode, fTnxReference365Days_iostat)
end subroutine fTnxReference365Days_open


subroutine fTnxReference365Days_write(line, advance_in)
    !! Writes the given line to the fTnxReference365Days file.
    character(len=*), intent(in) :: line
        !! line to write
    logical, intent(in), optional :: advance_in
        !! whether or not to append a newline character

    logical :: advance

    if (present(advance_in)) then
        advance = advance_in
    else
        advance = .true.
    end if
    call write_file(fTnxReference365Days, line, advance, fTnxReference365Days_iostat)
end subroutine fTnxReference365Days_write


subroutine fTnxReference365Days_close()
    close(fTnxReference365Days)
end subroutine fTnxReference365Days_close


subroutine fTnxReference365Days_erase()
    call unlink(GetTnxReference365DaysFileFull())
end subroutine fTnxReference365Days_erase


subroutine GetMonthlyTemperatureDataSetFromTnxReferenceFile(Monthi, TminDataSet, TmaxDataSet)
    integer(int32), intent(in) :: Monthi
    type(rep_DayEventDbl), dimension(31) , intent(inout) :: TminDataSet
    type(rep_DayEventDbl), dimension(31) , intent(inout) :: TmaxDataSet

    integer(int32) :: Dayi, DayN
    integer(int32) :: DNR
    integer(int32) :: t1, t2, ni
    real(dp) :: C1Min, C2Min, C3Min
    real(dp) :: C1Max, C2Max, C3Max
    real(dp) :: aOver3Min, bOver2Min, cMin
    real(dp) :: aOver3Max, bOver2Max, cMax

    ni=30

    if (Monthi == 1) then
        C1Min = GetTminTnxReference12MonthsRun_i(12)
        C1Max = GetTmaxTnxReference12MonthsRun_i(12)
    else
        C1Min = GetTminTnxReference12MonthsRun_i(Monthi-1)
        C1Max = GetTmaxTnxReference12MonthsRun_i(Monthi-1)
    end if
    C2Min = GetTminTnxReference12MonthsRun_i(Monthi)
    C2Max = GetTmaxTnxReference12MonthsRun_i(Monthi)
    if (Monthi == 12) then
        C3Min = GetTminTnxReference12MonthsRun_i(1)
        C3Max = GetTmaxTnxReference12MonthsRun_i(1)
    else
        C3Min = GetTminTnxReference12MonthsRun_i(Monthi+1)
        C3Max = GetTmaxTnxReference12MonthsRun_i(Monthi+1)
    end if

    C1Min = C1Min*ni
    C1Max = C1Max*ni
    C2Min = C2Min*ni
    C2Max = C2Max*ni
    C3Min = C3Min*ni
    C3Max = C3Max*ni

    call DetermineDayNr(1, Monthi, 1901, DNR)
    DayN = DaysInMonth(Monthi)

    call GetInterpolationParameters(C1Min, C2Min, C3Min, &
                   aOver3Min, bOver2Min, cMin)
    call GetInterpolationParameters(C1Max, C2Max, C3Max, &
                   aOver3Max, bOver2Max, cMax)

    t1 = 30
    do Dayi = 1, DayN
        t2 = t1 + 1
        TminDataSet(Dayi)%DayNr = DNR+Dayi-1
        TmaxDataSet(Dayi)%DayNr = DNR+Dayi-1
        TminDataSet(Dayi)%Param = aOver3Min*(t2*t2*t2-t1*t1*t1) &
            + bOver2Min*(t2*t2-t1*t1) + cMin*(t2-t1)
        TmaxDataSet(Dayi)%Param = aOver3Max*(t2*t2*t2-t1*t1*t1) &
            + bOver2Max*(t2*t2-t1*t1) + cMax*(t2-t1)
        t1 = t2
    end do
    do Dayi = (DayN+1), 31 ! Give to remaining days (day 29,30 or 31) the DayNr of the last day of the month
        TminDataSet(Dayi)%DayNr = DNR+DayN-1
        TmaxDataSet(Dayi)%DayNr = DNR+DayN-1
        TminDataSet(Dayi)%Param = 0._dp
        TmaxDataSet(Dayi)%Param = 0._dp
    end do

    contains


    subroutine GetInterpolationParameters(C1, C2, C3, &
                          aOver3, bOver2, c)
        real(dp), intent(in) :: C1
        real(dp), intent(in) :: C2
        real(dp), intent(in) :: C3
        real(dp), intent(inout) :: aOver3
        real(dp), intent(inout) :: bOver2
        real(dp), intent(inout) :: c

        ! n1=n2=n3=30 --> better parabola
        aOver3 = (C1-2._dp*C2+C3)/(6._dp*30._dp*30._dp*30._dp)
        bOver2 = (-6._dp*C1+9._dp*C2-3._dp*C3)/(6._dp*30._dp*30._dp)
        c = (11._dp*C1-7._dp*C2+2._dp*C3)/(6._dp*30._dp)
    end subroutine GetInterpolationParameters
end subroutine GetMonthlyTemperatureDataSetFromTnxReferenceFile

end module ac_tempprocessing
