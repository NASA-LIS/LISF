!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: AC72_setup
! \label{AC72_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   04 NOV 2024; Louise Busschaert, Initial implementation
!
! !INTERFACE:
subroutine AC72_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use LIS_timeMgrMod
    use LIS_mpiMod,    only: LIS_mpi_comm
    use LIS_constantsMod

    use module_sf_aclsm_72, only: &
           OC, WP, SAT, FC, INFRATE, SD, CL, SI 
    use AC72_lsmMod
    use ac72_prep_f,    only: ac72_read_Trecord

    use ac_global, only: &
            CheckForKeepSWC, &
            DegreesDay,&
            FromGravelMassToGravelVolume,&
            GetCCiActual,&
            GetCCiprev,&
            GetCCiTopEarlySen,&
            GetClimRecord,&
            GetClimRecord_FromDayNr,&
            GetClimRecord_FromString,&
            GetClimRecord_FromY,&
            GetClimRecord_ToDayNr,&
            GetClimRecord_ToString,&
            GetCompartment,&
            GetCompartment_theta,&
            GetCrop,&
            GetCrop_Day1,&
            GetCrop_DaysToHarvest,&
            GetCrop_ModeCycle,&
            GetCRsalt,&  
            GetCRwater,& 
            GetDaySubmerged,&
            GetDrain,&  
            GetEact,& 
            GetECdrain,& 
            GetECiAqua,& 
            GetECstorage,& 
            Geteffectiverain,&
            GetEndSeason,&
            GetEpot,& 
            GetETo,&
            GetEToRecord,&
            GetEvapoEntireSoilSurface,&
            GetFullFileNameProgramParameters, &
            GetGenerateDepthMode,&
            GetGenerateTimeMode,&
            GetInfiltrated,&
            GetIniPercTAW,&
            GetIrriAfterSeason,&
            GetIrriBeforeSeason,&
            GetIrriECw,&
            GetIrriFirstDayNr,&
            GetIrrigation,&
            GetIrriMethod,&
            GetIrriMode,&
            GetManagement,&
            GetManagement_Cuttings,&
            GetMaxPlotNew,&
            GetMaxPlotTr,&
            GetMultipleProjectFile, &
            GetMultipleProjectFileFull, &
            GetNrCompartments,&
            GetNumberSimulationRuns,&
            GetOnset,&
            GetOut1Wabal,&
            GetOut2Crop,&
            GetOut3Prof,&
            GetOut4Salt,&
            GetOut5CompWC,&
            GetOut6CompEC,&
            GetOut7Clim,&
            GetOut8Irri,&
            GetOutDaily,&
            GetOutputAggregate,&
            GetPart1Mult,&
            GetPart2Eval,&
            GetPathNameList,&
            GetPerennialPeriod,&
            GetPreDay,&
            GetRain,& 
            GetRainRecord,&
            GetRootingDepth,&
            GetRootZoneSalt,&
            GetRootZoneWC_Actual,&
            GetRootZoneWC_FC,&
            GetRootZoneWC_Leaf,&
            GetRootZoneWC_SAT,&
            GetRootZoneWC_Sen,&
            GetRootZoneWC_Thresh,&
            GetRootZoneWC_WP,&
            GetRootZoneWC_ZtopAct,&
            GetRootZoneWC_ZtopFC,&
            GetRootZoneWC_ZtopThresh,&
            GetRootZoneWC_ZtopWP,&
            GetRunoff,& 
            GetSaltInfiltr,&
            GetSimulation,&
            GetSimulation_MultipleRunConstZrx,&
            GetSimulation_MultipleRunConstZrx,&
            GetSimulation_MultipleRunWithKeepSWC,&
            GetSimulation_MultipleRunWithKeepSWC,&
            GetSimulation_SumGDD,&
            GetSimulation_SumGDDfromDay1,&
            GetSimulation_ToDayNr,&
            GetSimulation_ToDayNr,&
            GetSimulParam,&
            GetSimulParam_GDDMethod,&
            GetSimulParam_ThicknessTopSWC,&
            GetSoil,&
            GetSoilLayer,&
            GetSumWaBal,&
            GetSurf0,& 
            GetSurfaceStorage,&
            GetTact,&
            GetTactWeedInfested,&
            GetTemperatureRecord,&
            GetTmax,& 
            GetTmin,& 
            GetTotalSaltContent,&
            GetTotalWaterContent,&
            GetTpot,&
            GetZiAqua,&
            InitializeGlobalStrings,&
            IrriMode_Generate,&
            IrriMode_Manual, &
            ModeCycle_GDDays, &
            SetCCiActual,&
            SetCCiprev,&
            SetCCiTopEarlySen,&
            SetClimDescription,&
            SetClimFile,&
            SetClimRecord,&
            SetClimRecord_DataType,&
            SetClimRecord_fromd,&
            SetClimRecord_FromDayNr,&
            SetClimRecord_fromm,&
            SetClimRecord_FromString,&
            SetClimRecord_fromy,&
            SetClimRecord_FromY,&
            SetClimRecord_NrObs,&
            SetClimRecord_tod,&
            SetClimRecord_ToDayNr,&
            SetClimRecord_tom, &
            SetClimRecord_ToString,&
            SetClimRecord_toy,&
            SetCO2Description,&
            SetCompartment,&
            SetCompartment_theta,&
            SetCrop,&
            SetCRsalt,&  
            SetCRwater,& 
            SetDaySubmerged,&
            SetDrain,&  
            SetEact,& 
            SetECdrain,& 
            SetECiAqua,& 
            SetECstorage,& 
            Seteffectiverain,&
            SetEndSeason,&
            SetEpot,& 
            SetETo,&
            SetEToRecord,&
            SetEToRecord,&
            SetEvapoEntireSoilSurface,&
            SetFullFileNameProgramParameters,&
            SetGenerateDepthMode,&
            SetGenerateTimeMode,&
            SetInfiltrated,&
            SetIniPercTAW,&
            SetIrriAfterSeason,&
            SetIrriBeforeSeason,&
            SetIrriECw,&
            SetIrriFirstDayNr,&
            SetIrrigation,&
            SetIrriMethod,&
            SetIrriMode,&
            SetManagement,&
            SetManagement_Cuttings,&
            SetMaxPlotNew,&
            SetMaxPlotTr,&
            SetMultipleProjectDescription,&
            SetMultipleProjectFile, &
            SetMultipleProjectFileFull,&
            SetNrCompartments,&
            SetOnset,&
            SetOut1Wabal,&
            SetOut2Crop,&
            SetOut3Prof,&
            SetOut4Salt,&
            SetOut5CompWC,&
            SetOut6CompEC,&
            SetOut7Clim,&
            SetOutDaily,&
            SetOutputAggregate,&
            SetPart1Mult,&
            SetPart2Eval,&
            SetPathNameList,&
            SetPathNameOutp,&
            SetPathNameParam,&
            SetPathNameProg,&
            SetPathNameSimul,&
            SetPerennialPeriod,&
            SetPreDay,&
            SetProfFilefull,&
            SetRain,& 
            SetRainRecord,&
            SetRainRecord,&
            SetRootingDepth,&
            SetRootZoneSalt,&
            SetRootZoneWC_Actual,&
            SetRootZoneWC_FC,&
            SetRootZoneWC_Leaf,&
            SetRootZoneWC_SAT,&
            SetRootZoneWC_Sen,&
            SetRootZoneWC_Thresh,&
            SetRootZoneWC_WP,&
            SetRootZoneWC_ZtopAct,&
            SetRootZoneWC_ZtopFC,&
            SetRootZoneWC_ZtopThresh,&
            SetRootZoneWC_ZtopWP,&
            SetRunoff,& 
            SetSaltInfiltr,&
            SetSimulation,&
            SetSimulation_MultipleRun,&
            SetSimulation_MultipleRunConstZrx,&
            SetSimulation_MultipleRunWithKeepSWC,&
            SetSimulation_NrRuns, &
            SetSimulation_Storage_CropString,&
            SetSimulation_SumGDD, &
            SetSimulation_SumGDDfromDay1,&
            SetSimulation_ToDayNr, &
            SetSimulParam,&
            SetSimulParam_ThicknessTopSWC, &
            SetSoil,&
            SetSoilLayer,&
            SetSumWaBal,&
            SetSurf0,& 
            SetSurfaceStorage,&
            SetTact,&
            SetTactWeedInfested,&
            SetTemperatureRecord,&
            SetTmax,& 
            SetTmaxRun,& 
            SetTmaxTnxReference12MonthsRun,&
            SetTmin,&
            SetTminRun,&
            SetTminTnxReference12MonthsRun, &
            SetTnxReferenceFile,&
            SetTnxReferenceYear,&
            SetTotalSaltContent,&
            SetTotalWaterContent,&
            SetTpot,&
            SetZiAqua,&
            typeproject_typeprm,&
            typeproject_typepro,&
            undef_int
    use ac_project_input, only: ProjectInput, allocate_project_input, set_project_input

    use ac_run, only: &
            AdvanceOneTimeStep,&
            FinalizeRun1,&
            FinalizeRun2,&
            FinalizeSimulation,&
            fIrri_close,&
            GetalfaHI,&
            GetalfaHIAdj,&
            GetBin,&
            GetBout,&
            GetBprevSum,&
            GetCCiActualWeedInfested,&
            GetCCoTotal,&
            GetCCxCropWeedsNoSFstress,&
            GetCCxTotal,&
            GetCCxWitheredTpotNoS,&
            GetCDCTotal,&
            GetCGCref,&
            GetCO2i,&
            GetCoeffb0,&
            GetCoeffb0Salt,&
            GetCoeffb1,&
            GetCoeffb1Salt,&
            GetCoeffb2,&
            GetCoeffb2Salt,&
            GetCrop_Tbase, &
            GetCrop_Tupper, &
            GetCutInfoRecord1,&
            GetCutInfoRecord2,&
            GetDayFraction,&
            GetDayLastCut,&
            GetDayNr1Eval,&
            GetDayNrEval,&
            GetDayNri,&
            GetFracBiomassPotSF,&
            GetfWeedNoS,&
            GetGDDayFraction,&
            GetGDDayi,&
            GetGDDCDCTotal,&
            GetGDDCGCref ,&
            GetGDDTadj,&
            GetGwTable,&
            GetHItimesAT,&
            GetHItimesAT1,&
            GetHItimesAT2,&
            GetHItimesBEF,&
            GetIrriInfoRecord1,&
            GetIrriInfoRecord2,&
            GetIrriInterval,&
            GetLineNrEval,&
            GetNextSimFromDayNr,&
            GetNoMoreCrop,&
            GetNoYear,&
            GetNrCut,&
            GetPlotVarCrop,&
            GetPreviousBmob,&
            GetPreviousBsto,&
            GetPreviousDayNr,&
            GetPreviousStressLevel,&
            GetPreviousSum,&
            GetPreviousSumETo,&
            GetPreviousSumGDD,&
            GetScorAT1,&
            GetScorAT2,&
            GetSimulation_ToDayNr,&
            GetStageCode,&
            GetStartMode,&
            GetStressLeaf,&
            GetStressSenescence,&
            GetStressSFadjNEW,&
            GetStressTot,&
            GetSumETo,&
            GetSumGDD,&
            GetSumGDDcuts,&
            GetSumGDDPrev,&
            GetSumInterval,&
            GetSumKci,&
            GetSumKcTop,&
            GetSumKcTopStress,&
            GetTadj,&
            GetTheProjectFile,&
            GetTimeSenescence ,&
            GetTransfer,&
            GetWaterTableInProfile,&
            GetWeedRCi,&
            GetYprevSum,&
            GetZeval,&
            GetZiprev,&
            InitializeClimate,&
            InitializeClimate,&
            InitializeRunPart1,&
            InitializeRunPart2,&
            InitializeSimulation,&
            InitializeSimulationRunPart2,&
            ReadClimateNextDay,&
            SetalfaHI,&
            SetalfaHIAdj,&
            SetBin,&
            SetBout,&
            SetBprevSum,&
            SetCCiActualWeedInfested,&
            SetCCoTotal,&
            SetCCxCropWeedsNoSFstress,&
            SetCCxTotal,&
            SetCCxWitheredTpotNoS,&
            SetCDCTotal,&
            SetCGCref,&
            SetCO2i,&
            SetCoeffb0,&
            SetCoeffb0Salt,&
            SetCoeffb1,&
            SetCoeffb1Salt,&
            SetCoeffb2,&
            SetCoeffb2Salt,&
            SetCutInfoRecord1,&
            SetCutInfoRecord2,&
            SetDayFraction,&
            SetDayLastCut,&
            SetDayNr1Eval,&
            SetDayNrEval,&
            SetDayNri,&
            SetETo,&
            SetFracBiomassPotSF,&
            SetfWeedNoS,&
            SetGDDayFraction,&
            SetGDDayi,&
            SetGDDCDCTotal,&
            SetGDDCGCref,&
            SetGDDTadj,&
            SetGDDVariablesNextDay,&
            SetGwTable,&
            SetHItimesAT,&
            SetHItimesAT1,&
            SetHItimesAT2,&
            SetHItimesBEF,&
            SetIrriInfoRecord1,&
            SetIrriInfoRecord2,&
            SetIrriInterval,&
            SetLineNrEval,&
            SetNextSimFromDayNr,&
            SetNoMoreCrop,&
            SetNoYear,&
            SetNrCut,&
            SetPlotVarCrop,&
            SetPreviousBmob,&
            SetPreviousBsto,&
            SetPreviousDayNr,&
            SetPreviousStressLevel,&
            SetPreviousSum,&
            SetPreviousSumETo,&
            SetPreviousSumGDD,&
            SetRain,& 
            SetScorAT1,&
            SetScorAT2,&
            SetStageCode,&
            SetStartMode,&
            SetStressLeaf,&
            SetStressSenescence ,&
            SetStressSFadjNEW,&
            SetStressTot,&
            SetSumETo,&
            SetSumGDD,&
            SetSumGDDcuts,&
            SetSumGDDPrev,&
            SetSumInterval,&
            SetSumKci,&
            SetSumKcTop,&
            SetSumKcTopStress,&
            SetTadj,&
            SetTheProjectFile,&
            SetTimeSenescence,&
            SetTmax,& 
            SetTmin,&
            SetTransfer,&
            SetWaterTableInProfile,&
            SetWeedRCi,&
            SetYprevSum,&
            SetZeval,&
            SetZiprev

        use ac_kinds, only: intEnum, &
                            int32, &
                            int8, &
                            sp
                             
        use ac_startunit, only: &
            ComposeFileForProgramParameters, &
            GetListProjectsFile, &
            GetNumberOfProjects, &
            GetProjectFileName, &
            GetProjectType, &
            GetRequestDailyResults, &
            GetRequestParticularResults, &
            GetSimulation_NrRuns, &
            GetTimeAggregationResults, &
            InitializeProject, &
            InitializeTheProgram, &
            LoadProgramParametersProjectPlugIn, &
            WriteProjectsInfo

        use ac_initialsettings, only: InitializeSettings

    !
    ! !DESCRIPTION:
    !
    !  This routine is the entry point to set up the parameters
    !  required for AC72.  
    !  The routines invoked are:
    !  \begin{description}
    !  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
    !    retrieves LIS parameter data from NetCDF file
    !  \item[AC72\_read\_MULTILEVEL\_param](\ref{AC72_read_MULTILEVEL_param}) \newline
    !    retrieves MULTILEVEL spatial parameter from NetCDF file
    !  \item[AC72\_read\_croptype](\ref{AC72_read_croptype}) \newline
    !    retrieves the crop type from LIS_param and table, then assigns
    !    a crop type to each tile.
    !  \end{description}
    !EOP
        implicit none
        integer           :: mtype
        integer           :: t, k, n, l
        integer           :: col, row
        real, allocatable :: placeholder(:,:)

        real              :: Z_surf, cl_tmp, si_tmp, sd_tmp, InfRate_tmp
        integer           :: REW, descr, ierr
        integer           :: time1julhours, timerefjulhours
        integer           :: time1days, time2days
        integer           :: yr2, mo2, da2, hr2, mn2
        
        integer :: daynr, todaynr, iproject, nprojects, NrRuns
        integer(intEnum) :: TheProjectType
        logical :: ListProjectFileExist
        character(len=:), allocatable :: ListProjectsFile, TheProjectFile

        logical ::  ProgramParametersAvailable 
        integer(int32) :: TotalSimRuns

        logical :: MultipleRunWithKeepSWC_temp    
        real    :: MultipleRunConstZrx_temp

        mtype = LIS_rc%lsm_index

        
        do n=1, LIS_rc%nnest
            ! allocate memory for place holder for #n nest
            allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))

            ! Check for AC simulation period start at LIS start for coldstart
            if (trim(LIS_rc%startcode) .eq. "coldstart") then
                if ((AC72_struc(n)%Sim_AnnualStartMonth.eq.LIS_rc%smo) &
                    .and.(AC72_struc(n)%Sim_AnnualStartDay.eq.LIS_rc%sda)) then
                    write(LIS_logunit, *) 'AC coldstart: simulation period start matches LIS start'
                else
                    write(LIS_logunit, *) 'AC coldstart: simulation period start does not match LIS start'
                    write(LIS_logunit, *) 'program stopping ...'
                    call LIS_endrun
                endif
            endif
            
            !------------------------------------!
            ! reading spatial spatial parameters !
            !------------------------------------!
            ! croptype
            ! Read crop type:

            call ac72_read_croptype(n)

            ! soiltype takes value from the LIS built-in parameter soilt
            if(LIS_rc%usetexturemap(n) .ne. 'none') then
                write(LIS_logunit,*) "AC72: retrieve parameter SOILTYPE from LIS"
                do t=1, LIS_rc%npatch(n, mtype)
                    AC72_struc(n)%ac72(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
                enddo
            else 
                ! read: soiltype
                write(LIS_logunit,*) "AC72: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
                call LIS_read_param(n, trim(AC72_struc(n)%LDT_ncvar_soiltype), placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    AC72_struc(n)%ac72(t)%soiltype = placeholder(col, row)
                enddo 
            endif

            write(LIS_logunit,*) "[INFO] AC72: reading parameter AC_Tmin_clim from ",&
                trim(LIS_rc%paramfile(n))
            do k = 1, 12
                call AC72_read_MULTILEVEL_param(n, AC72_struc(n)%LDT_ncvar_tmincli_monthly, k, placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    ! Climatology is rounded to 2 decimals in AquaCrop
                    AC72_struc(n)%ac72(t)%tmincli_monthly(k) = anint(placeholder(col, row)*100)/100 - LIS_CONST_TKFRZ
                enddo 
            enddo

            write(LIS_logunit,*) "[INFO] AC72: reading parameter AC_Tmax_clim from ",&
                trim(LIS_rc%paramfile(n))
            do k = 1, 12
                call AC72_read_MULTILEVEL_param(n, AC72_struc(n)%LDT_ncvar_tmaxcli_monthly, k, placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    ! Climatology is rounded to 2 decimals in AquaCrop
                    AC72_struc(n)%ac72(t)%tmaxcli_monthly(k) = anint(placeholder(col, row)*100)/100 - LIS_CONST_TKFRZ
                enddo 
            enddo 
            deallocate(placeholder)
            ! Read soil table
            call SOIL_PARM_AC72(AC72_struc(n)%soil_tbl_name)
            TheProjectType = typeproject_typeprm

            ! Create AquaCrop 'Run' Structure
            call LIS_get_julhr(1901,1,1,0,0,0,timerefjulhours)
            TotalSimRuns = LIS_rc%eyr - LIS_rc%syr + 1
            call allocate_project_input(TotalSimRuns)
            do l=1, TotalSimRuns  ! TotalSimRuns
                call set_project_input(l, 'Simulation_YearSeason', 1_int8)
                ! Simulation
                call LIS_get_julhr(LIS_rc%syr+(l-1), AC72_struc(n)%Sim_AnnualStartMonth, &
                                   AC72_struc(n)%Sim_AnnualStartDay,0,0,0,time1julhours)
                time1days = (time1julhours - timerefjulhours)/24 + 1
                ! Find last day of simulation (check for leap year)
                if ((((mod(LIS_rc%syr+(l-1),4).eq.0.and.mod(LIS_rc%syr+(l-1),100).ne.0) &
                    .or.(mod(LIS_rc%syr+(l-1),400).eq.0)).and.(LIS_rc%smo.le.2)) &
                    .or.(((mod(LIS_rc%syr+(l),4).eq.0.and.mod(LIS_rc%syr+(l),100).ne.0) &
                    .or.(mod(LIS_rc%syr+(l),400).eq.0)).and.(LIS_rc%smo.gt.2))) then ! leap year sim period
                    time2days = time1days + 365
                else ! no leap year
                    time2days = time1days + 364
                endif
                call set_project_input(l, 'Simulation_DayNr1', time1days)
                call set_project_input(l, 'Simulation_DayNrN', time2days)
                ! Crop
                call LIS_get_julhr(LIS_rc%syr+(l-1),AC72_struc(n)%Crop_AnnualStartMonth, &
                                   AC72_struc(n)%Crop_AnnualStartDay,0,0,0,time1julhours)
                time1days = (time1julhours - timerefjulhours)/24 + 1
                !Note: end of cropping period is defined by crop params
                call set_project_input(l, 'Crop_Day1', time1days)
                call set_project_input(l, 'Crop_DayN', time2days)
                ! Note: '(External)' input sources can vary spatially (e.g. soil, crop, meteo)
                call set_project_input(l, 'Description', ' LIS ')
                call set_project_input(l, 'Climate_Info', '(External)')
                call set_project_input(l, 'Climate_Filename', '(External)')
                call set_project_input(l, 'Climate_Directory', '(None)')
                call set_project_input(l, 'VersionNr', 7.2)
                call set_project_input(l, 'Temperature_Info', '(None)')
                call set_project_input(l, 'Temperature_Filename', '(External)')
                call set_project_input(l, 'Temperature_Directory', '(None)')
                call set_project_input(l, 'ETo_Info', '(None)')
                call set_project_input(l, 'ETo_Filename', '(External)')
                call set_project_input(l, 'ETo_Directory', '(None)')
                call set_project_input(l, 'Rain_Info', ' LIS ')
                call set_project_input(l, 'Rain_Filename', '(External)')
                call set_project_input(l, 'Rain_Directory', '(None)')
                call set_project_input(l, 'CO2_Info', ' LIS ')
                call set_project_input(l, 'CO2_Filename', trim(AC72_struc(n)%CO2_Filename))
                call set_project_input(l, 'CO2_Directory', trim(AC72_struc(n)%PathNameSimul))
                call set_project_input(l, 'Calendar_Info', '(None)')
                call set_project_input(l, 'Calendar_Filename', '(None)')
                call set_project_input(l, 'Calendar_Directory', '(None)')
                call set_project_input(l, 'Crop_Info', '(None)')
                call set_project_input(l, 'Crop_Filename', '(External)')
                call set_project_input(l, 'Crop_Directory', trim(AC72_struc(n)%PathCropFiles))
                call set_project_input(l, 'Irrigation_Info', ' LIS ')
                call set_project_input(l, 'Irrigation_Filename', trim(AC72_struc(n)%Irrigation_Filename))
                call set_project_input(l, 'Irrigation_Directory', trim(AC72_struc(n)%PathNameSimul))
                call set_project_input(l, 'Management_Info', ' LIS ')
                call set_project_input(l, 'Management_Filename',  trim(AC72_struc(n)%Management_Filename))
                call set_project_input(l, 'Management_Directory', trim(AC72_struc(n)%PathNameSimul))
                call set_project_input(l, 'GroundWater_Info', '(None)')
                call set_project_input(l, 'GroundWater_Filename', '(None)')
                call set_project_input(l, 'GroundWater_Directory', '(None)')
                call set_project_input(l, 'Soil_Info', '(External)')
                call set_project_input(l, 'Soil_Filename', '(External)')
                call set_project_input(l, 'Soil_Directory', '(External)')
                call set_project_input(l, 'SWCIni_Info', '(None)')
                if (l == 1) then
                    call set_project_input(l, 'SWCIni_Filename', '(None)') 
                    !Initial conditions are deifned in config
                else
                    call set_project_input(l, 'SWCIni_Filename', 'KeepSWC')
                end if
                call set_project_input(l, 'SWCIni_Directory', '(None)')
                call set_project_input(l, 'OffSeason_Info', '(None)')
                call set_project_input(l, 'OffSeason_Filename', '(None)')
                call set_project_input(l, 'OffSeason_Directory', '(None)')
                call set_project_input(l, 'Observations_Info', '(None)')
                call set_project_input(l, 'Observations_Filename', '(None)')
                call set_project_input(l, 'Observations_Directory', '(None)')
            end do

            ! Read annual temperature record
            call ac72_read_Trecord(n)

            do t = 1, LIS_rc%npatch(n, mtype)
                
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row

                ! Initialize global strings
                call InitializeGlobalStrings()

                call SetPathNameSimul(trim(AC72_struc(n)%PathNameSimul)) ! needed for InitializeSettings
                ! Initializes all settings to default, later overwritten in main
                call InitializeSettings(use_default_soil_file=.false.,use_default_crop_file=.false.)
                ! It uses the defaults anyway
                ! Set crop file (crop parameters are read when calling InitializeRunPart1)
                ! Set the crop type cst in time
                do l=1, TotalSimRuns
                    call set_project_input(l, 'Crop_Filename', &
                                        trim(AC72_struc(n)%ac72(t)%cropt)//'.CRO')
                enddo
                call CheckForKeepSWC(MultipleRunWithKeepSWC_temp, &
                                     MultipleRunConstZrx_temp)
                call SetSimulation_MultipleRunWithKeepSWC(MultipleRunWithKeepSWC_temp) 
                call SetSimulation_MultipleRunConstZrx(MultipleRunConstZrx_temp)  
                !!! Start the Program AC72
                call SetPathNameOutp('')
                call SetPathNameList('')
                call SetPathNameParam('')
                call SetPathNameProg('')


                call GetTimeAggregationResults()
                call GetRequestDailyResults()
                call GetRequestParticularResults()

                call SetClimRecord_DataType(0_int8)
                call SetClimRecord_fromd(LIS_rc%sda)
                call SetClimRecord_fromdaynr(ProjectInput(1)%Simulation_DayNr1)
                call SetClimRecord_fromm(LIS_rc%smo)
                call SetClimRecord_fromstring("")
                call SetClimRecord_fromy(LIS_rc%syr)
                call SetClimRecord_NrObs(999)
                call SetClimRecord_tod(LIS_rc%eda)
                call SetClimRecord_todaynr(ProjectInput(GetSimulation_NrRuns())%Simulation_DayNrN)
                call SetClimRecord_tom(LIS_rc%emo)
                call SetClimRecord_tostring("")
                call SetClimRecord_toy(LIS_rc%eyr)
                call SetClimFile('(External)')
                call SetCO2Description('')
                call SetProfFilefull('External')
                call SetTemperatureRecord(GetClimRecord())
                call SetEToRecord(GetClimRecord())
                call SetRainRecord(GetClimRecord())
                call SetSimulation_Storage_CropString('None')

                call SetMultipleProjectDescription('undefined')
                ProgramParametersAvailable = .false.
                call SetFullFileNameProgramParameters("(None)")
                call LoadProgramParametersProjectPlugIn(&
                      GetFullFileNameProgramParameters(), &
                             ProgramParametersAvailable)
                call SetSimulation_MultipleRun(.true.)
                call SetSimulation_NrRuns(TotalSimRuns) 
                MultipleRunWithKeepSWC_temp = GetSimulation_MultipleRunWithKeepSWC() 
                MultipleRunConstZrx_temp = GetSimulation_MultipleRunConstZrx()

                ! Set AC72 soil parameters based on soiltype
                ! Note: up to 5 soil layers with the same texture
                AC72_struc(n)%ac72(t)%SoilLayer(1)%wp = WP(AC72_struc(n)%ac72(t)%soiltype) * 100
                AC72_struc(n)%ac72(t)%SoilLayer(1)%sat = SAT(AC72_struc(n)%ac72(t)%soiltype) * 100
                AC72_struc(n)%ac72(t)%SoilLayer(1)%fc = FC(AC72_struc(n)%ac72(t)%soiltype) * 100
                AC72_struc(n)%ac72(t)%SoilLayer(1)%InfRate = INFRATE(AC72_struc(n)%ac72(t)%soiltype)
                sd_tmp = sd(AC72_struc(n)%ac72(t)%soiltype)*100
                cl_tmp = cl(AC72_struc(n)%ac72(t)%soiltype)*100
                si_tmp = si(AC72_struc(n)%ac72(t)%soiltype)*100
                    
                ! define default CN
                if (AC72_struc(n)%ac72(t)%SoilLayer(1)%InfRate>864) then
                    AC72_struc(n)%ac72(t)%Soil%CNvalue = 46
                elseif (AC72_struc(n)%ac72(t)%SoilLayer(1)%InfRate>=347) then
                    AC72_struc(n)%ac72(t)%Soil%CNvalue = 61
                elseif (AC72_struc(n)%ac72(t)%SoilLayer(1)%InfRate>=36) then
                    AC72_struc(n)%ac72(t)%Soil%CNvalue = 72
                else
                    AC72_struc(n)%ac72(t)%Soil%CNvalue = 77
                endif
                
                ! define default REW (only top layer will be used)
                Z_surf = 0.04
                REW=nint(10.0*(AC72_struc(n)%ac72(t)%SoilLayer(1)%fc-(AC72_struc(n)%ac72(t)%SoilLayer(1)%wp/2.0))*Z_surf)
                if (REW < 0) REW = 0
                if (REW > 15) REW = 15
                AC72_struc(n)%ac72(t)%Soil%REW = REW

                !  associate soil class with USDA soil type for soil description
                if ((1.5 * cl_tmp + si_tmp) < 15.0) then
                    descr = 0
                elseif (((1.5 * cl_tmp + si_tmp) >= 15) .and. ((2 * cl_tmp + si_tmp) <= 30)) then
                    descr = 1
                elseif ((cl_tmp >= 7 .and. cl_tmp < 20 .and. sd_tmp >= 52) .and. (2 * cl_tmp + si_tmp >= 30)) then
                    descr = 2
                elseif ((cl_tmp < 7) .and. (si_tmp < 50) .and. (2 * cl_tmp + si_tmp >= 30)) then
                    descr = 2
                elseif ((cl_tmp >= 7) .and. (cl_tmp < 27) .and. (si_tmp >= 28) .and. (si_tmp < 50) .and. (sd_tmp < 52)) then
                    descr = 3
                elseif ((si_tmp >= 50) .and. (cl_tmp >= 12) .and. (cl_tmp < 27)) then
                    descr = 4
                elseif ((si_tmp >= 50) .and. (si_tmp < 80) .and. (cl_tmp < 12)) then
                    descr = 4
                elseif ((si_tmp >= 80) .and. (cl_tmp < 12)) then
                    descr = 5
                elseif ((cl_tmp >= 20) .and. (cl_tmp < 35) .and. (si_tmp < 28) .and. (sd_tmp > 45)) then
                    descr = 6
                elseif ((cl_tmp >= 27) .and. (cl_tmp < 40) .and. (sd_tmp >= 20) .and. (sd_tmp < 45)) then
                    descr = 7
                elseif ((cl_tmp >= 27) .and. (cl_tmp < 40) .and. (sd_tmp < 20)) then
                    descr = 8
                elseif ((cl_tmp >= 35) .and. (cl_tmp < 55) .and. (sd_tmp >= 45) .and. (sd_tmp < 65)) then
                    descr = 9
                elseif ((cl_tmp >= 40) .and. (si_tmp >= 40)) then
                    descr = 10
                elseif ((cl_tmp >= 40) .and. (sd_tmp < 45) .and. (si_tmp < 40)) then
                    descr = 11
                else
                   write(LIS_logunit, *) 'no soil texture found'
                   write(LIS_logunit, *) 'program stopping ...'
                   call LIS_endrun
                end if
                InfRate_tmp = AC72_struc(n)%ac72(t)%SoilLayer(1)%InfRate
                if ((descr == 0) .or. (descr == 1) .or. (descr == 2)) then
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRa = -0.3112-10**(-5)*InfRate_tmp
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRb = -1.4936+0.2416*log(InfRate_tmp)
                elseif ((descr == 3) .or. (descr == 4) .or. (descr == 5)) then
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRa = -0.4986-9*10**(-5)*InfRate_tmp
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRb = -2.1320+0.4778*log(InfRate_tmp)
                elseif ((descr == 6) .or. (descr == 7) .or. (descr == 9)) then
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRa = -0.5677-4*10**(-5)*InfRate_tmp
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRb = -3.7189+0.5922*log(InfRate_tmp)
                elseif ((descr == 8) .or. (descr == 10) .or. (descr == 11)) then
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRa = -0.6366+8*10**(-4)*InfRate_tmp
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%CRb = -1.9165+0.7163*log(InfRate_tmp)
                endif
                ! Set GravelMass and Penetrability
                AC72_struc(n)%ac72(t)%SoilLayer(1)%Penetrability = 100.0
                AC72_struc(n)%ac72(t)%SoilLayer(1)%GravelMass = 10.0
                ! Set GravelVol
                AC72_struc(n)%ac72(t)%SoilLayer(1)%GravelVol = &
                    FromGravelMassToGravelVolume(AC72_struc(n)%ac72(t)%SoilLayer(1)%sat,&
                    AC72_struc(n)%ac72(t)%SoilLayer(1)%GravelMass)

                AC72_struc(n)%ac72(t)%SoilLayer(1)%Description = 'soil type from LIS'


                AC72_struc(n)%ac72(t)%SoilLayer(1)%Thickness = AC72_struc(n)%Thickness(1)

                ! loop if more than 2 layers
                if (AC72_struc(n)%NrSoilLayers.gt.1) then
                    do l = 2, AC72_struc(n)%NrSoilLayers
                        AC72_struc(n)%ac72(t)%SoilLayer(l) = AC72_struc(n)%ac72(t)%SoilLayer(1)
                        AC72_struc(n)%ac72(t)%SoilLayer(l)%Thickness = AC72_struc(n)%Thickness(l)
                        AC72_struc(n)%ac72(t)%SoilLayer(l)%GravelMass = 0.0
                        AC72_struc(n)%ac72(t)%SoilLayer(l)%GravelVol = 0.0
                    enddo
                endif

                if (AC72_struc(n)%NrSoilLayers.lt.5) then !repeat the last layer
                    do l = AC72_struc(n)%NrSoilLayers+1, 5
                        AC72_struc(n)%ac72(t)%SoilLayer(l) = AC72_struc(n)%ac72(t)%SoilLayer(AC72_struc(n)%NrSoilLayers)
                    enddo
                endif
                
                AC72_struc(n)%ac72(t)%Soil%NrSoilLayers = AC72_struc(n)%NrSoilLayers
                AC72_struc(n)%ac72(t)%NrCompartments = AC72_struc(n)%max_No_compartments

                ! Set soil global variables 
                call SetSoilLayer(AC72_struc(n)%ac72(t)%SoilLayer)
                call SetSoil(AC72_struc(n)%ac72(t)%Soil)
                call SetNrCompartments(AC72_struc(n)%ac72(t)%NrCompartments)

                call SetPathNameProg('')
                !
                call SetMultipleProjectDescription('undefined')
                ProgramParametersAvailable = .false.
                ! Running with default program parameters
                call SetFullFileNameProgramParameters("(None)")
                call LoadProgramParametersProjectPlugIn(&
                      GetFullFileNameProgramParameters(), &
                             ProgramParametersAvailable)
                call SetSimulation_MultipleRun(.true.)
                call SetSimulation_NrRuns(TotalSimRuns)

                ! InitializeSimulation (year)
                AC72_struc(n)%ac72(t)%irun = 1

                call SetClimRecord_DataType(0_int8)
                call SetClimRecord_fromd(0)
                call SetClimRecord_fromdaynr(ProjectInput(1)%Simulation_DayNr1)
                call SetClimRecord_fromm(0)
                call SetClimRecord_fromstring("")
                call SetClimRecord_fromy(LIS_rc%syr)
                call SetClimRecord_NrObs(999)
                call SetClimRecord_tod(0)
                call SetClimRecord_todaynr(ProjectInput(GetSimulation_NrRuns())%Simulation_DayNrN)
                call SetClimRecord_tom(0)
                call SetClimRecord_tostring("")
                call SetClimRecord_toy(0)
                call SetClimFile('(External)')

                AC72_struc(n)%ac72(t)%WPi = 0.

                ! Set Global variable to pass T record to AquaCrop
                call SetTminRun(AC72_struc(n)%ac72(t)%Tmin_record)    
                call SetTmaxRun(AC72_struc(n)%ac72(t)%Tmax_record)
                call SetTmin(AC72_struc(n)%ac72(t)%Tmin_record(1))    
                call SetTmax(AC72_struc(n)%ac72(t)%Tmax_record(1))

                ! Set Tmin and Tmax reference to compute the stress realtions
                call SetTminTnxReference12MonthsRun(AC72_struc(n)%ac72(t)%tmincli_monthly(:))
                call SetTmaxTnxReference12MonthsRun(AC72_struc(n)%ac72(t)%tmaxcli_monthly(:))

                ! Set reference year for CO2 for stress functions
                call SetTnxReferenceYear(AC72_struc(n)%tempcli_refyr)
                call SetTnxReferenceFile('(External)')

                ! InitializeRunPart1
                call InitializeRunPart1(int(AC72_struc(n)%ac72(t)%irun, kind=int8), AC72_struc(n)%ac72(t)%TheProjectType)
                call InitializeSimulationRunPart2()
                AC72_struc(n)%ac72(t)%InitializeRun = 0
                AC72_struc(n)%ac72(t)%read_Trecord = 0
                ! Check if enough GDDays to complete cycle, if not, turn on flag to warn the user
                AC72_struc(n)%AC72(t)%crop = GetCrop()
                if(GetCrop_ModeCycle().eq.ModeCycle_GDDays)then
                    if (((GetCrop_Day1()+GetCrop_DaysToHarvest()).gt.GetSimulation_ToDayNr()) &
                        .or.(GetCrop_DaysToHarvest()<1)) then
                        AC72_struc(n)%ac72(t)%cycle_complete = 0
                    else
                        AC72_struc(n)%ac72(t)%cycle_complete = 1
                    endif
                endif
                ! Close irrigation file after Run Initialization
                ! Note: only 2 irrigation records can be passed
                if((GetIrriMode().eq.IrriMode_Generate)&
                   .or.(GetIrriMode().eq.IrriMode_Manual)) then
                    call fIrri_close()
                endif      

                ! Set AC72_struc after Initialization
                AC72_struc(n)%AC72(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
                AC72_struc(n)%AC72(t)%RootZoneWC_FC = GetRootZoneWC_FC()
                AC72_struc(n)%AC72(t)%RootZoneWC_WP = GetRootZoneWC_WP()
                AC72_struc(n)%AC72(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
                AC72_struc(n)%AC72(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
                AC72_struc(n)%AC72(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
                AC72_struc(n)%AC72(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
                AC72_struc(n)%AC72(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
                AC72_struc(n)%AC72(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
                AC72_struc(n)%AC72(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
                AC72_struc(n)%AC72(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
                AC72_struc(n)%AC72(t)%Compartment = GetCompartment()
                AC72_struc(n)%AC72(t)%TotalSaltContent = GetTotalSaltContent()
                AC72_struc(n)%AC72(t)%TotalWaterContent = GetTotalWaterContent()
                AC72_struc(n)%AC72(t)%effectiverain = Geteffectiverain()
                AC72_struc(n)%AC72(t)%SumWaBal = GetSumWaBal()
                AC72_struc(n)%AC72(t)%RootZoneSalt = GetRootZoneSalt()
                AC72_struc(n)%AC72(t)%Simulation = GetSimulation()
                AC72_struc(n)%AC72(t)%IrriInterval = GetIrriInterval()
                AC72_struc(n)%AC72(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
                AC72_struc(n)%AC72(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
                AC72_struc(n)%AC72(t)%Irrigation = GetIrrigation()
                AC72_struc(n)%AC72(t)%IrriBeforeSeason = GetIrriBeforeSeason()
                AC72_struc(n)%AC72(t)%IrriAfterSeason = GetIrriAfterSeason()
                AC72_struc(n)%AC72(t)%SoilLayer = GetSoilLayer()
                AC72_struc(n)%AC72(t)%daynri = GetDayNri()
                AC72_struc(n)%AC72(t)%NrCompartments = GetNrCompartments()
                do l=1, AC72_struc(n)%AC72(t)%NrCompartments
                     AC72_struc(n)%AC72(t)%smc(l) = GetCompartment_theta(l)
                enddo
                AC72_struc(n)%AC72(t)%IrriECw = GetIrriECw()
                AC72_struc(n)%AC72(t)%Management = GetManagement()
                AC72_struc(n)%AC72(t)%PerennialPeriod = GetPerennialPeriod()
                AC72_struc(n)%AC72(t)%simulparam = GetSimulParam()
                AC72_struc(n)%AC72(t)%Cuttings = GetManagement_Cuttings()
                AC72_struc(n)%AC72(t)%onset = GetOnset()
                AC72_struc(n)%AC72(t)%endseason = GetEndSeason()
                AC72_struc(n)%AC72(t)%Soil = GetSoil()
                AC72_struc(n)%AC72(t)%TemperatureRecord = GetTemperatureRecord()
                AC72_struc(n)%AC72(t)%ClimRecord = GetClimRecord()
                AC72_struc(n)%AC72(t)%RainRecord = GetRainRecord()
                AC72_struc(n)%AC72(t)%EToRecord = GetEToRecord()

                AC72_struc(n)%AC72(t)%GenerateTimeMode = GetGenerateTimeMode()
                AC72_struc(n)%AC72(t)%GenerateDepthMode = GetGenerateDepthMode()
                AC72_struc(n)%AC72(t)%IrriMode = GetIrriMode()
                AC72_struc(n)%AC72(t)%IrriMethod = GetIrriMethod()
                AC72_struc(n)%AC72(t)%DaySubmerged = GetDaySubmerged()
                AC72_struc(n)%AC72(t)%MaxPlotNew = GetMaxPlotNew()
                AC72_struc(n)%AC72(t)%IrriFirstDayNr = GetIrriFirstDayNr()
                AC72_struc(n)%AC72(t)%ZiAqua = GetZiAqua()
                AC72_struc(n)%AC72(t)%IniPercTAW = GetIniPercTAW()
                AC72_struc(n)%AC72(t)%MaxPlotTr = GetMaxPlotTr()
                AC72_struc(n)%AC72(t)%OutputAggregate = GetOutputAggregate()

                AC72_struc(n)%AC72(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
                AC72_struc(n)%AC72(t)%PreDay = GetPreDay()
                AC72_struc(n)%AC72(t)%OutDaily = GetOutDaily()
                AC72_struc(n)%AC72(t)%Out1Wabal = GetOut1Wabal()
                AC72_struc(n)%AC72(t)%Out2Crop = GetOut2Crop()
                AC72_struc(n)%AC72(t)%Out3Prof = GetOut3Prof()
                AC72_struc(n)%AC72(t)%Out4Salt = GetOut4Salt()
                AC72_struc(n)%AC72(t)%Out5CompWC = GetOut5CompWC()
                AC72_struc(n)%AC72(t)%Out6CompEC = GetOut6CompEC()
                AC72_struc(n)%AC72(t)%Out7Clim = GetOut7Clim()
                AC72_struc(n)%AC72(t)%Out8Irri = GetOut8Irri()
                AC72_struc(n)%AC72(t)%Part1Mult = GetPart1Mult()
                AC72_struc(n)%AC72(t)%Part2Eval = GetPart2Eval()

                !
                AC72_struc(n)%AC72(t)%CCiActual = GetCCiActual()
                AC72_struc(n)%AC72(t)%CCiprev = GetCCiprev()
                AC72_struc(n)%AC72(t)%CCiTopEarlySen = GetCCiTopEarlySen()
                AC72_struc(n)%AC72(t)%CRsalt = GetCRsalt ()
                AC72_struc(n)%AC72(t)%CRwater = GetCRwater()
                AC72_struc(n)%AC72(t)%ECdrain = GetECdrain()
                AC72_struc(n)%AC72(t)%ECiAqua = GetECiAqua()
                AC72_struc(n)%AC72(t)%ECstorage = GetECstorage()
                AC72_struc(n)%AC72(t)%Eact = GetEact()
                AC72_struc(n)%AC72(t)%Epot = GetEpot()
                AC72_struc(n)%AC72(t)%ETo = GetETo()
                AC72_struc(n)%AC72(t)%Drain = GetDrain()
                AC72_struc(n)%AC72(t)%Infiltrated = GetInfiltrated()
                AC72_struc(n)%AC72(t)%prcp = GetRain()
                AC72_struc(n)%AC72(t)%RootingDepth = GetRootingDepth()
                AC72_struc(n)%AC72(t)%Runoff = GetRunoff()
                AC72_struc(n)%AC72(t)%SaltInfiltr = GetSaltInfiltr()
                AC72_struc(n)%AC72(t)%Surf0 = GetSurf0()
                AC72_struc(n)%AC72(t)%SurfaceStorage = GetSurfaceStorage()
                AC72_struc(n)%AC72(t)%Tact = GetTact()
                AC72_struc(n)%AC72(t)%Tpot = GetTpot()
                AC72_struc(n)%AC72(t)%TactWeedInfested = 0. !not ini in AC GetTactWeedInfested()
                AC72_struc(n)%AC72(t)%Tmax = GetTmax()
                AC72_struc(n)%AC72(t)%Tmin = GetTmin()


                AC72_struc(n)%AC72(t)%GwTable = GetGwTable()
                AC72_struc(n)%AC72(t)%PlotVarCrop = GetPlotVarCrop()
                AC72_struc(n)%AC72(t)%StressTot = GetStressTot()
                AC72_struc(n)%AC72(t)%CutInfoRecord1 = GetCutInfoRecord1()
                AC72_struc(n)%AC72(t)%CutInfoRecord2 = GetCutInfoRecord2()
                AC72_struc(n)%AC72(t)%Transfer = GetTransfer()
                AC72_struc(n)%AC72(t)%PreviousSum = GetPreviousSum()
                AC72_struc(n)%AC72(t)%Tadj = GetTadj()
                AC72_struc(n)%AC72(t)%GDDTadj = GetGDDTadj()
                AC72_struc(n)%AC72(t)%DayLastCut = GetDayLastCut()
                AC72_struc(n)%AC72(t)%NrCut = GetNrCut()
                AC72_struc(n)%AC72(t)%SumInterval = GetSumInterval()
                AC72_struc(n)%AC72(t)%PreviousStressLevel = GetPreviousStressLevel()
                AC72_struc(n)%AC72(t)%StressSFadjNEW = GetStressSFadjNEW()
                AC72_struc(n)%AC72(t)%Bin = GetBin()
                AC72_struc(n)%AC72(t)%Bout = GetBout()
                AC72_struc(n)%AC72(t)%GDDayi = GetGDDayi()
                AC72_struc(n)%AC72(t)%CO2i = GetCO2i()
                AC72_struc(n)%AC72(t)%FracBiomassPotSF = GetFracBiomassPotSF()
                AC72_struc(n)%AC72(t)%SumETo = GetSumETo()
                AC72_struc(n)%AC72(t)%SumGDD = GetSumGDD()
                AC72_struc(n)%AC72(t)%Ziprev = GetZiprev()
                AC72_struc(n)%AC72(t)%SumGDDPrev = GetSumGDDPrev()
                AC72_struc(n)%AC72(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
                AC72_struc(n)%AC72(t)%Coeffb0 = GetCoeffb0()
                AC72_struc(n)%AC72(t)%Coeffb1 = GetCoeffb1()
                AC72_struc(n)%AC72(t)%Coeffb2 = GetCoeffb2()
                AC72_struc(n)%AC72(t)%Coeffb0Salt = GetCoeffb0Salt()
                AC72_struc(n)%AC72(t)%Coeffb1Salt = GetCoeffb1Salt()
                AC72_struc(n)%AC72(t)%Coeffb2Salt = GetCoeffb2Salt()
                AC72_struc(n)%AC72(t)%StressLeaf = GetStressLeaf()
                AC72_struc(n)%AC72(t)%StressSenescence = GetStressSenescence()
                AC72_struc(n)%AC72(t)%DayFraction = GetDayFraction()
                AC72_struc(n)%AC72(t)%GDDayFraction = GetGDDayFraction()
                AC72_struc(n)%AC72(t)%CGCref = GetCGCref()
                AC72_struc(n)%AC72(t)%GDDCGCref = GetGDDCGCref()
                AC72_struc(n)%AC72(t)%TimeSenescence = GetTimeSenescence()
                AC72_struc(n)%AC72(t)%SumKcTop = GetSumKcTop()
                AC72_struc(n)%AC72(t)%SumKcTopStress = GetSumKcTopStress()
                AC72_struc(n)%AC72(t)%SumKci = GetSumKci()
                AC72_struc(n)%AC72(t)%CCoTotal = GetCCoTotal()
                AC72_struc(n)%AC72(t)%CCxTotal = GetCCxTotal()
                AC72_struc(n)%AC72(t)%CDCTotal = GetCDCTotal()
                AC72_struc(n)%AC72(t)%GDDCDCTotal = GetGDDCDCTotal()
                AC72_struc(n)%AC72(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
                AC72_struc(n)%AC72(t)%WeedRCi = GetWeedRCi()
                AC72_struc(n)%AC72(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
                AC72_struc(n)%AC72(t)%fWeedNoS = GetfWeedNoS()
                AC72_struc(n)%AC72(t)%Zeval = GetZeval()
                AC72_struc(n)%AC72(t)%BprevSum = GetBprevSum()
                AC72_struc(n)%AC72(t)%YprevSum = GetYprevSum()
                AC72_struc(n)%AC72(t)%SumGDDcuts = GetSumGDDcuts()
                AC72_struc(n)%AC72(t)%HItimesBEF = GetHItimesBEF()
                AC72_struc(n)%AC72(t)%ScorAT1 = GetScorAT1()
                AC72_struc(n)%AC72(t)%ScorAT2 = GetScorAT2()
                AC72_struc(n)%AC72(t)%HItimesAT1 = GetHItimesAT1()
                AC72_struc(n)%AC72(t)%HItimesAT2 = GetHItimesAT2()
                AC72_struc(n)%AC72(t)%HItimesAT = GetHItimesAT()
                AC72_struc(n)%AC72(t)%alfaHI = GetalfaHI()
                AC72_struc(n)%AC72(t)%alfaHIAdj = GetalfaHIAdj()
                AC72_struc(n)%AC72(t)%NextSimFromDayNr = GetNextSimFromDayNr()
                AC72_struc(n)%AC72(t)%DayNr1Eval = GetDayNr1Eval()
                AC72_struc(n)%AC72(t)%DayNrEval = GetDayNrEval()
                AC72_struc(n)%AC72(t)%LineNrEval = GetLineNrEval()
                AC72_struc(n)%AC72(t)%PreviousSumETo = GetPreviousSumETo()
                AC72_struc(n)%AC72(t)%PreviousSumGDD = GetPreviousSumGDD()
                AC72_struc(n)%AC72(t)%PreviousBmob = GetPreviousBmob()
                AC72_struc(n)%AC72(t)%PreviousBsto = GetPreviousBsto()
                AC72_struc(n)%AC72(t)%StageCode = GetStageCode()
                AC72_struc(n)%AC72(t)%PreviousDayNr = GetPreviousDayNr()
                AC72_struc(n)%AC72(t)%NoYear = GetNoYear()
                AC72_struc(n)%AC72(t)%WaterTableInProfile = GetWaterTableInProfile()
                AC72_struc(n)%AC72(t)%StartMode = GetStartMode()
                AC72_struc(n)%AC72(t)%NrRuns = GetSimulation_NrRuns()
                AC72_struc(n)%AC72(t)%TheProjectType = TheProjectType

                ! Check for irrigation (irrigation file management)
                if(AC72_struc(n)%ac72(t)%IrriMode.eq.IrriMode_Manual)then
                    if(AC72_struc(n)%ac72(t)%IrriInfoRecord1%NoMoreInfo)then
                        AC72_struc(n)%ac72(t)%irri_lnr = 9
                    else
                        AC72_struc(n)%ac72(t)%irri_lnr = 10
                    endif
                elseif(AC72_struc(n)%ac72(t)%IrriMode.eq.IrriMode_Generate)then
                    if(AC72_struc(n)%ac72(t)%IrriInfoRecord1%NoMoreInfo)then
                        AC72_struc(n)%ac72(t)%irri_lnr = 11
                    else
                        AC72_struc(n)%ac72(t)%irri_lnr = 12
                    endif
                else ! no irrigation, set to 0
                    AC72_struc(n)%ac72(t)%irri_lnr = 0
                endif

                ! Check if we need to start a new sim period
                call LIS_get_julhr(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                                  0,0,0,time1julhours)
                time1days = (time1julhours - timerefjulhours)/24
                ! If we restart on the first day of simulation
                ! Do not read Trecord in main but initialize run
                if ((AC72_struc(n)%Sim_AnnualStartMonth.eq.LIS_rc%smo) &
                    .and.(AC72_struc(n)%Sim_AnnualStartDay.eq.LIS_rc%sda)) then
                    AC72_struc(n)%ac72(t)%InitializeRun = 1
                endif
                
        enddo ! do t = 1, LIS_rc%npatch(n, mtype)
    enddo
end subroutine AC72_setup

!BOP
!
! !ROUTINE: AC72_read_MULTILEVEL_param
!  \label{AC72_read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  04 NOV Louise Busschaert; Initial implementation
!
! !INTERFACE:
subroutine AC72_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
! !USES:
    use netcdf
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_localPet,   &   
                            LIS_ews_halo_ind, LIS_ewe_halo_ind, &
                            LIS_nss_halo_ind, LIS_nse_halo_ind   
    use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    implicit none
! !ARGUMENTS: 
    integer, intent(in)          :: n
    integer, intent(in)          :: level
    character(len=*), intent(in) :: ncvar_name 
    real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
! !DESCRIPTION:
!  This subroutine reads MULTILEVEL parameters from the LIS
!  NetCDF parameter data file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[level]
!    level index (month, quarter, soil layer, snow layer) of the data to be read
!   \item[array]
!    array containing returned values
!   \end{description}
!
!EOP      

    integer       :: ios1
    integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
    integer       :: nc, nr, t, nlevel, k
    real, pointer :: level_data(:, :, :)
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in AC72_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AC72_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AC72_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC72_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC72_read_MULTILEVEL_param')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

        ! inquire the length of the level dimension
        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

        ! allocate memory
        allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in AC72_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in AC72_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) 'MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) 'program stopping ...'
        call LIS_endrun
    endif
 end subroutine AC72_read_MULTILEVEL_param

!BOP
!
! !ROUTINE: AC72_read_croptype
!  \label{AC72_read_croptype}
!
! !REVISION HISTORY:
!  04 NOV 2024; Louise Busschaert, initial implementation
!
! !INTERFACE:
subroutine ac72_read_croptype(n)
! !USES
    use ESMF
    use LIS_fileIOMod
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_coreMod
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use AC72_lsmMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    implicit none
! !ARGUMENTS
    integer, intent(in)     :: n

! !DESCRIPTION:
!  This subroutine checks the crop classifcation and assigns a
!  a crop (corresponding to a crop file .CRO) to each tile.
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \end{description}
!
!EOP   
    character(len=LIS_CONST_PATH_LEN) :: crop_path
    integer                 :: ftn
    integer                 :: t,j,col,row,IINDEX,CT,n_crops
    integer                 :: nid, ios, status, croptypeId
    integer                 :: rc
    logical                 :: file_exists
    character(len=100),allocatable :: croptypes(:)
    real, allocatable       :: l_croptype(:,:)
    real, allocatable       :: glb_croptype(:,:)
    real, allocatable       :: glb_croptype1(:,:)
    character*128 :: mess
! __________________________________________________________________________
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    call ESMF_ConfigFindLabel(LIS_config,"AquaCrop.7.2 crop library directory:",rc=rc)
    call ESMF_ConfigGetAttribute(LIS_config,crop_path,rc=rc)
    call LIS_verify(rc,'AquaCrop.7.2 crop library directory: not specified')

 !- Set path to crop files
      AC72_struc(n)%PathCropFiles = trim(crop_path)//"AC_Crop.Files/"

 !- Read in crop table
      open(19, FILE=trim(crop_path)//"AC_Crop.Inventory",FORM='FORMATTED',STATUS='OLD',IOSTAT=rc)
      call LIS_verify(rc,'AC72_setup.F: failure opening AC_Crop.Inventory')

      read(19,*) n_crops

      allocate(character(len=100) :: croptypes(n_crops))

      write( mess , *) 'AC_Crop.Inventory contains ', n_crops, ' types' 
      if (LIS_masterproc) then
            write(LIS_logunit, *) trim(mess)
      endif
      read(19,*)
      do CT=1,n_crops
        read(19,*) IINDEX, croptypes(CT)
      enddo
      close(19)

    allocate(l_croptype(LIS_rc%lnc(n),LIS_rc%lnr(n)))


 !- Read in crop type map file (specified in LIS parameter input file)
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       ios = nf90_open(path=LIS_rc%paramfile(n),&
                       mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) "[INFO] Reading in the crop type field ... "
       
       allocate(glb_croptype(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       
       ios = nf90_inq_varid(nid,'AC_CROPT',croptypeId)
       call LIS_verify(ios,'nf90_inq_varid failed for AC_CROPT')
       
       ios = nf90_get_var(nid, croptypeId, glb_croptype)
       call LIS_verify(ios,'nf90_get_var failed for AC_CROPT')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in AC72_setup')
       
       l_croptype(:,:) = glb_croptype(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1))
       deallocate(glb_croptype)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          
          AC72_struc(n)%ac72(t)%cropt = trim(croptypes(nint(l_croptype(col,row))))
       enddo
       deallocate(croptypes)

    else
       write(LIS_logunit,*) "[ERR] The croptype inventory: ",&
             LIS_rc%paramfile(n)," does not exist."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun
    endif
    deallocate(l_croptype)
#endif

  end subroutine ac72_read_croptype
                                          

