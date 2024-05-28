!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: Ac71_setup
! \label{Ac71_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   2023; Michel Bechtold, Initial implementation
!   06 MAR 2024; Louise Busschaert, Clean-up
!
! !INTERFACE:
subroutine Ac71_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface, LIS_masterproc
    use LIS_timeMgrMod
    use LIS_mpiMod,    only: LIS_mpi_comm

    use module_sf_aclsm_71, only: &
           OC, WP, SAT, FC, INFRATE, SD, CL, SI 
    use Ac71_lsmMod
    use ac71_prep_f,    only: ac71_read_Trecord

    use ac_project_input, only: set_project_input, &
                                allocate_project_input
    use ac_global, only:    GetCCiActual,&
                            GetCCiprev,&
                            GetCCiTopEarlySen,&
                            GetClimRecord,&
                            GetCompartment,&
                            GetCrop,&
                            GetDaySubmerged,&
                            GetDrain,&  
                            GetEact,& 
                            GetECstorage,& 
                            GetFullFileNameProgramParameters, &
                            GetIrriAfterSeason,&
                            GetIrriBeforeSeason,&
                            GetIrriECw,&
                            GetIrrigation,&
                            GetManagement,&
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
                            GetSimulation,&
                            GetSimulation_SumGDD,&
                            GetSimulation_ToDayNr, &
                            GetSoil,&
                            GetSoilLayer,&
                            GetSumWaBal,&
                            GetSurfaceStorage,&
                            GetTact,&
                            GetTactWeedInfested,&
                            GetTpot,&
                            SetCCiActual,&
                            SetClimFile,&
                            SetClimRecord_DataType, &
                            SetClimRecord_fromd, &
                            SetClimRecord_FromDayNr,&
                            SetClimRecord_fromm, &
                            SetClimRecord_FromString,&
                            SetClimRecord_fromy, &
                            SetClimRecord_NrObs, &
                            SetClimRecord_tod, &
                            SetClimRecord_ToDayNr,&
                            SetClimRecord_tom, &
                            SetClimRecord_ToString,&
                            SetClimRecord_toy,&
                            SetCO2Description,&
                            SetCompartment,&
                            SetETo,& 
                            SetEToRecord, &
                            SetFullFileNameProgramParameters, &
                            SetMultipleProjectDescription, &
                            SetNrCompartments,&
                            SetPathNameList, &
                            SetPathNameOutp, &
                            SetPathNameParam, &
                            SetPathNameProg, &
                            SetPathNameSimul, &
                            SetProfFilefull,&
                            SetRain,& 
                            SetRainRecord,&
                            SetSimulation,&
                            SetSimulation_MultipleRun, &
                            SetSimulation_NrRuns, &
                            SetSimulation_Storage_CropString,&
                            SetSoil,&
                            SetSoilLayer,&
                            SetTemperatureRecord, &
                            SetTmax,& 
                            SetTmin,&
                            typeproject_typeprm
    use ac_project_input, only: ProjectInput 
    use ac_run, only:   GetalfaHI,&
                        GetalfaHIAdj,&
                        GetBin,&
                        GetBout,&
                        GetCCiActualWeedInfested,&
                        GetCCoTotal,&
                        GetCCxCropWeedsNoSFstress,&
                        GetCCxTotal,&
                        GetCCxWitheredTpotNoS,&
                        GetCDCTotal,&
                        GetCGCref,&
                        GetCoeffb0,&
                        GetCoeffb0Salt,&
                        GetCoeffb1,&
                        GetCoeffb1Salt,&
                        GetCoeffb2,&
                        GetCoeffb2Salt,&
                        GetDayFraction,&
                        GetDayNri,&
                        GetGDDayFraction,&
                        GetGDDayi,&
                        GetGDDCDCTotal,&
                        GetGDDCGCref ,&
                        GetGDDTadj,&
                        GetHItimesAT,&
                        GetHItimesAT1,&
                        GetHItimesAT2,&
                        GetHItimesBEF,&
                        GetIrriInfoRecord1,&
                        GetIrriInfoRecord2,&
                        GetNoMoreCrop,&
                        GetPreviousStressLevel,&
                        GetScorAT1,&
                        GetScorAT2,&
                        GetSimulation_ToDayNr, &
                        GetStressLeaf,&
                        GetStressSenescence ,&
                        GetStressSFadjNEW,&
                        GetStressTot,&
                        GetSumGDDcuts,&
                        GetSumInterval,&
                        GetSumKci,&
                        GetSumKcTop,&
                        GetSumKcTopStress,&
                        GetTadj,&
                        GetTimeSenescence ,&
                        GetWeedRCi,&
                        GetYprevSum,&
                        GetZiprev,&
                        InitializeRunPart1 ,&
                        InitializeSimulation ,&
                        InitializeSimulationRunPart2, &
                        SetETo, &
                        SetRain,& 
                        SetTmax,& 
                        SetTmin
        use ac_kinds, only: intEnum, &
                            int32, &
                            int8, &
                            dp,&
                            sp                           
        use ac_startunit, only: GetSimulation_NrRuns, &
                                GetTimeAggregationResults, &
                                GetRequestDailyResults, &
                                GetRequestParticularResults, &
                                LoadProgramParametersProjectPlugIn
        use ac_initialsettings, only: InitializeSettings
    !
    ! !DESCRIPTION:
    !
    !  This routine is the entry point to set up the parameters
    !  required for Ac71.  
    !  The routines invoked are:
    !  \begin{description}
    !  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
    !    retrieves LIS parameter data from NetCDF file
    !  \item[AC71\_read\_MULTILEVEL\_param](\ref{AC71_read_MULTILEVEL_param}) \newline
    !    retrieves MULTILEVEL spatial parameter from NetCDF file
    !  \item[AC71\_read\_croptype](\ref{AC71_read_croptype}) \newline
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

        character(256) :: irr_dir, man_dir

        mtype = LIS_rc%lsm_index

        
        do n=1, LIS_rc%nnest
            ! allocate memory for place holder for #n nest
            allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))

            ! Check for AC simulation period start at LIS start for coldstart
            if (trim(LIS_rc%startcode) .eq. "coldstart") then
                if ((AC71_struc(n)%Sim_AnnualStartMonth.eq.LIS_rc%smo) &
                    .and.(AC71_struc(n)%Sim_AnnualStartDay.eq.LIS_rc%sda)) then
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

            call ac71_read_croptype(n)

            ! soiltype takes value from the LIS built-in parameter soilt
            if(LIS_rc%usetexturemap(n) .ne. 'none') then
                write(LIS_logunit,*) "Ac71: retrieve parameter SOILTYPE from LIS"
                do t=1, LIS_rc%npatch(n, mtype)
                    AC71_struc(n)%ac71(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
                enddo
            else 
                ! read: soiltype
                write(LIS_logunit,*) "Ac71: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
                call LIS_read_param(n, trim(AC71_struc(n)%LDT_ncvar_soiltype), placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    AC71_struc(n)%ac71(t)%soiltype = placeholder(col, row)
                enddo 
            endif
            deallocate(placeholder)
            ! Read soil table
            call SOIL_PARM_71(AC71_struc(n)%soil_tbl_name)
            TheProjectType = typeproject_typeprm

            ! Create AquaCrop 'Run' Structure
            call LIS_get_julhr(1901,1,1,0,0,0,timerefjulhours)
            TotalSimRuns = LIS_rc%eyr - LIS_rc%syr + 1
            call allocate_project_input(TotalSimRuns)
            do l=1, TotalSimRuns  ! TotalSimRuns
                !! Check irrigation and management
                ! Check if irrigation ON
                if(LIS_rc%irrigation_type.ne."none") then
                    if(trim(LIS_rc%irrigation_type).eq."Sprinkler") then
                        AC71_struc(n)%Irrigation_Filename = "sprinkler.IRR"
                        irr_dir = trim(AC71_struc(n)%PathNameSimul)
                    else ! Other options can be implemented later
                        write(LIS_logunit, *) trim(LIS_rc%irrigation_type), &
                        " irrigation type not compatible with AquaCrop.7.1"
                        call LIS_endrun()
                    endif
                else
                    AC71_struc(n)%Irrigation_Filename = "(None)"
                    irr_dir = "(None)"
                endif

                ! Check if management file
                if(trim(AC71_struc(n)%Management_Filename).ne."(None)") then
                    man_dir = trim(AC71_struc(n)%PathNameSimul)
                else
                    man_dir = "(None)"
                endif
                !! End check irrigation and management

                ! MB: for current generic crop this is fixed to 1
                call set_project_input(l, 'Simulation_YearSeason', 1_int8)
                ! Simulation
                call LIS_get_julhr(LIS_rc%syr+(l-1), AC71_struc(n)%Sim_AnnualStartMonth, &
                                   AC71_struc(n)%Sim_AnnualStartDay,0,0,0,time1julhours)
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
                call LIS_get_julhr(LIS_rc%syr+(l-1),AC71_struc(n)%Crop_AnnualStartMonth, &
                                   AC71_struc(n)%Crop_AnnualStartDay,0,0,0,time1julhours)
                time1days = (time1julhours - timerefjulhours)/24 + 1
                !Note: end of cropping period is defined by crop params
                call set_project_input(l, 'Crop_Day1', time1days)
                call set_project_input(l, 'Crop_DayN', time2days)
                ! Note: '(External)' input sources can vary spatially (e.g. soil, crop, meteo)
                call set_project_input(l, 'Description', ' LIS ')
                call set_project_input(l, 'Climate_Info', '(External)')
                call set_project_input(l, 'Climate_Filename', '(External)')
                call set_project_input(l, 'Climate_Directory', '(None)')
                call set_project_input(l, 'VersionNr', 7.1_dp)
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
                call set_project_input(l, 'CO2_Filename', trim(AC71_struc(n)%CO2_Filename))
                call set_project_input(l, 'CO2_Directory', trim(AC71_struc(n)%PathNameSimul))
                call set_project_input(l, 'Calendar_Info', '(None)')
                call set_project_input(l, 'Calendar_Filename', '(None)')
                call set_project_input(l, 'Calendar_Directory', '(None)')
                call set_project_input(l, 'Crop_Info', '(None)')
                call set_project_input(l, 'Crop_Filename', '(External)')
                call set_project_input(l, 'Crop_Directory', trim(AC71_struc(n)%PathCropFiles))
                call set_project_input(l, 'Irrigation_Info', ' LIS ')
                call set_project_input(l, 'Irrigation_Filename', trim(AC71_struc(n)%Irrigation_Filename))
                call set_project_input(l, 'Irrigation_Directory', trim(irr_dir))
                call set_project_input(l, 'Management_Info', ' LIS ')
                call set_project_input(l, 'Management_Filename',  trim(AC71_struc(n)%Management_Filename))
                call set_project_input(l, 'Management_Directory', trim(man_dir))
                call set_project_input(l, 'GroundWater_Info', '(None)')
                call set_project_input(l, 'GroundWater_Filename', '(None)')
                call set_project_input(l, 'GroundWater_Directory', '(None)')
                call set_project_input(l, 'Soil_Info', '(External)')
                call set_project_input(l, 'Soil_Filename', '(External)')
                call set_project_input(l, 'Soil_Directory', '(External)')
                call set_project_input(l, 'SWCIni_Info', '(None)')
                if (l == 1) then
                    call set_project_input(l, 'SWCIni_Filename', '(None)') !Initial conditions are deifned in config
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
            ! Get end day and month
            call LIS_julhr_date((time2days-1)*24 + timerefjulhours,yr2,mo2,da2,hr2)
            AC71_struc(n)%Sim_AnnualEndMonth = mo2
            AC71_struc(n)%Sim_AnnualEndDay = da2

            ! Read annual temperature record if GDD_Mode
            if(LIS_masterproc.and.(AC71_struc(n)%GDD_Mode.eq.1)) then 
                call ac71_read_Trecord(n)
            endif
            ! Other processes wait for the master
            write(LIS_logunit, *) 'Waiting for temeprature reading for AC71 simul period...'
            call mpi_barrier(LIS_mpi_comm, ierr)

            do t = 1, LIS_rc%npatch(n, mtype)
                
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row

                call SetPathNameSimul(trim(AC71_struc(n)%PathNameSimul)) ! needed for InitializeSettings
                ! Initializes all settings to default, later overwritten in main
                call InitializeSettings(use_default_soil_file=.false.,use_default_crop_file=.false.)
                ! It uses the defaults anyway
                !!! Start the Program AC71
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
                call SetSimulation_Storage_CropString('None') ! perennial crops for later

                ! Set AC71 soil parameters based on soiltype
                ! Note: 2 soil layers with the same texture
                AC71_struc(n)%ac71(t)%SoilLayer(1)%wp = WP(AC71_struc(n)%ac71(t)%soiltype) * 100
                AC71_struc(n)%ac71(t)%SoilLayer(1)%sat = SAT(AC71_struc(n)%ac71(t)%soiltype) * 100
                AC71_struc(n)%ac71(t)%SoilLayer(1)%fc = FC(AC71_struc(n)%ac71(t)%soiltype) * 100
                AC71_struc(n)%ac71(t)%SoilLayer(1)%InfRate = INFRATE(AC71_struc(n)%ac71(t)%soiltype) * 86400000
                sd_tmp = sd(AC71_struc(n)%ac71(t)%soiltype)*100
                cl_tmp = cl(AC71_struc(n)%ac71(t)%soiltype)*100
                si_tmp = si(AC71_struc(n)%ac71(t)%soiltype)*100
                    
                ! define default CN
                if (AC71_struc(n)%ac71(t)%SoilLayer(1)%InfRate>864) then
                    AC71_struc(n)%ac71(t)%Soil%CNvalue = 46
                elseif (AC71_struc(n)%ac71(t)%SoilLayer(1)%InfRate>=347) then
                    AC71_struc(n)%ac71(t)%Soil%CNvalue = 61
                elseif (AC71_struc(n)%ac71(t)%SoilLayer(1)%InfRate>=36) then
                    AC71_struc(n)%ac71(t)%Soil%CNvalue = 72
                else
                    AC71_struc(n)%ac71(t)%Soil%CNvalue = 77
                endif
                
                ! define default REW (only top layer will be used)
                Z_surf = 0.04_dp
                REW=nint(10.0_dp*(AC71_struc(n)%ac71(t)%SoilLayer(1)%fc-(AC71_struc(n)%ac71(t)%SoilLayer(1)%wp/2.0))*Z_surf)
                if (REW < 0) REW = 0
                if (REW > 15) REW = 15
                AC71_struc(n)%ac71(t)%Soil%REW = REW

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
                InfRate_tmp = AC71_struc(n)%ac71(t)%SoilLayer(1)%InfRate
                if ((descr == 0) .or. (descr == 1) .or. (descr == 2)) then
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRa = -0.3112-10**(-5)*InfRate_tmp
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRb = -1.4936+0.2416*log(InfRate_tmp)
                elseif ((descr == 3) .or. (descr == 4) .or. (descr == 5)) then
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRa = -0.4986-9*10**(-5)*InfRate_tmp
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRb = -2.1320+0.4778*log(InfRate_tmp)
                elseif ((descr == 6) .or. (descr == 7) .or. (descr == 9)) then
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRa = -0.5677-4*10**(-5)*InfRate_tmp
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRb = -3.7189+0.5922*log(InfRate_tmp)
                elseif ((descr == 8) .or. (descr == 10) .or. (descr == 11)) then
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRa = -0.6366+8*10**(-4)*InfRate_tmp
                    AC71_struc(n)%ac71(t)%SoilLayer(1)%CRb = -1.9165+0.7163*log(InfRate_tmp)
                endif
                ! Set GravelMass and Penetrability
                AC71_struc(n)%ac71(t)%SoilLayer(1)%Penetrability = 100.0
                AC71_struc(n)%ac71(t)%SoilLayer(1)%GravelMass = 10.0
                AC71_struc(n)%ac71(t)%SoilLayer(1)%Description = 'soil type from LIS'


                AC71_struc(n)%ac71(t)%SoilLayer(1)%Thickness = AC71_struc(n)%Thickness(1)

                ! loop if more than 2 layers
                if (AC71_struc(n)%NrSoilLayers.gt.1) then
                    do l = 2, AC71_struc(n)%NrSoilLayers
                        AC71_struc(n)%ac71(t)%SoilLayer(l) = AC71_struc(n)%ac71(t)%SoilLayer(1)
                        AC71_struc(n)%ac71(t)%SoilLayer(l)%Thickness = AC71_struc(n)%Thickness(l)
                        AC71_struc(n)%ac71(t)%SoilLayer(l)%GravelMass = 0.0
                    enddo
                endif

                if (AC71_struc(n)%NrSoilLayers.lt.5) then !repeat the last layer
                    do l = AC71_struc(n)%NrSoilLayers+1, 5
                        AC71_struc(n)%ac71(t)%SoilLayer(l) = AC71_struc(n)%ac71(t)%SoilLayer(AC71_struc(n)%NrSoilLayers)
                    enddo
                endif
                
                AC71_struc(n)%ac71(t)%Soil%NrSoilLayers = AC71_struc(n)%NrSoilLayers
                AC71_struc(n)%ac71(t)%NrCompartments = AC71_struc(n)%max_No_compartments

                ! Set soil global variables 
                call SetSoilLayer(AC71_struc(n)%ac71(t)%SoilLayer)
                call SetSoil(AC71_struc(n)%ac71(t)%Soil)
                call SetNrCompartments(AC71_struc(n)%ac71(t)%NrCompartments)

                call SetPathNameProg('')
                !
                call SetMultipleProjectDescription('undefined')
                ProgramParametersAvailable = .false.
                call SetFullFileNameProgramParameters("(None)") ! Running with default rogram parameters
                call LoadProgramParametersProjectPlugIn(&
                      GetFullFileNameProgramParameters(), &
                             ProgramParametersAvailable)
                call SetSimulation_MultipleRun(.true.)
                call SetSimulation_NrRuns(TotalSimRuns)

                ! InitializeSimulation (year)
                AC71_struc(n)%ac71(t)%irun = 1

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

                AC71_struc(n)%ac71(t)%WPi = 0._dp

                ! Set crop file (crop parameters are read when calling InitializeRunPart1)
                call set_project_input(int(AC71_struc(n)%ac71(t)%irun, kind=int32), &
                                       'Crop_Filename', &
                                        trim(AC71_struc(n)%ac71(t)%cropt)//'.CRO')

                ! Set Global variable to pass T record to AquaCrop
                call SetTminRun(AC71_struc(n)%Trecord(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)%Tmin_record)    
                call SetTmaxRun(AC71_struc(n)%Trecord(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)%Tmax_record)
                ! InitializeRunPart1    
                call InitializeRunPart1(int(AC71_struc(n)%ac71(t)%irun, kind=int8), AC71_struc(n)%ac71(t)%TheProjectType)
                call InitializeSimulationRunPart2()
                AC71_struc(n)%ac71(t)%HarvestNow = .false.
                AC71_struc(n)%ac71(t)%InitializeRun = 0

                ! Set AC71_struc after Initialization
                AC71_struc(n)%ac71(t)%Bin = GetBin()
                AC71_struc(n)%ac71(t)%Bout = GetBout()
                AC71_struc(n)%ac71(t)%CCiActual = GetCCiActual()
                AC71_struc(n)%ac71(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
                AC71_struc(n)%ac71(t)%CCiTopEarlySen = GetCCiTopEarlySen()
                AC71_struc(n)%ac71(t)%CCiprev = GetCCiprev()
                AC71_struc(n)%ac71(t)%CCoTotal = GetCCoTotal()
                AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
                AC71_struc(n)%ac71(t)%CCxTotal = GetCCxTotal()
                AC71_struc(n)%ac71(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
                AC71_struc(n)%ac71(t)%CDCTotal = GetCDCTotal()
                AC71_struc(n)%ac71(t)%CGCref = GetCGCref()
                AC71_struc(n)%ac71(t)%Coeffb0 = GetCoeffb0()
                AC71_struc(n)%ac71(t)%Coeffb0Salt = GetCoeffb0Salt()
                AC71_struc(n)%ac71(t)%Coeffb1 = GetCoeffb1()
                AC71_struc(n)%ac71(t)%Coeffb1Salt = GetCoeffb1Salt()
                AC71_struc(n)%ac71(t)%Coeffb2 = GetCoeffb2()
                AC71_struc(n)%ac71(t)%Coeffb2Salt = GetCoeffb2Salt()
                AC71_struc(n)%ac71(t)%Compartment = GetCompartment()
                AC71_struc(n)%ac71(t)%DayFraction = GetDayFraction()
                AC71_struc(n)%ac71(t)%DaySubmerged = GetDaySubmerged()
                AC71_struc(n)%ac71(t)%Drain = GetDrain()  
                AC71_struc(n)%ac71(t)%ECstorage = GetECstorage()
                AC71_struc(n)%ac71(t)%Eact = GetEact()
                AC71_struc(n)%ac71(t)%GDDCDCTotal = GetGDDCDCTotal()
                AC71_struc(n)%ac71(t)%GDDCGCref = GetGDDCGCref ()
                AC71_struc(n)%ac71(t)%GDDTadj = GetGDDTadj()
                AC71_struc(n)%ac71(t)%GDDayFraction = GetGDDayFraction()
                AC71_struc(n)%ac71(t)%GDDayi = GetGDDayi()
                AC71_struc(n)%ac71(t)%HItimesAT = GetHItimesAT()
                AC71_struc(n)%ac71(t)%HItimesAT1 = GetHItimesAT1()
                AC71_struc(n)%ac71(t)%HItimesAT2 = GetHItimesAT2()
                AC71_struc(n)%ac71(t)%HItimesBEF = GetHItimesBEF()
                AC71_struc(n)%ac71(t)%IrriAfterSeason = GetIrriAfterSeason()
                AC71_struc(n)%ac71(t)%IrriBeforeSeason = GetIrriBeforeSeason()
                AC71_struc(n)%ac71(t)%IrriECw = GetIrriECw()
                AC71_struc(n)%ac71(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
                AC71_struc(n)%ac71(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
                AC71_struc(n)%ac71(t)%Irrigation = GetIrrigation()
                AC71_struc(n)%ac71(t)%Management = GetManagement()
                AC71_struc(n)%ac71(t)%NoMoreCrop = GetNoMoreCrop()
                AC71_struc(n)%ac71(t)%PreviousStressLevel = GetPreviousStressLevel()
                AC71_struc(n)%ac71(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
                AC71_struc(n)%ac71(t)%RootZoneWC_FC = GetRootZoneWC_FC()
                AC71_struc(n)%ac71(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
                AC71_struc(n)%ac71(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
                AC71_struc(n)%ac71(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
                AC71_struc(n)%ac71(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
                AC71_struc(n)%ac71(t)%RootZoneWC_WP = GetRootZoneWC_WP()
                AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
                AC71_struc(n)%ac71(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
                AC71_struc(n)%ac71(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
                AC71_struc(n)%ac71(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
                AC71_struc(n)%ac71(t)%Runoff = GetRunoff()
                AC71_struc(n)%ac71(t)%ScorAT1 = GetScorAT1()
                AC71_struc(n)%ac71(t)%ScorAT2 = GetScorAT2()
                AC71_struc(n)%ac71(t)%Simulation = GetSimulation()
                AC71_struc(n)%ac71(t)%Soil = GetSoil()
                AC71_struc(n)%ac71(t)%SoilLayer = GetSoilLayer()
                AC71_struc(n)%ac71(t)%StressLeaf = GetStressLeaf()
                AC71_struc(n)%ac71(t)%StressSFadjNEW = GetStressSFadjNEW()
                AC71_struc(n)%ac71(t)%StressSenescence = GetStressSenescence ()
                AC71_struc(n)%ac71(t)%StressTot = GetStressTot()
                AC71_struc(n)%ac71(t)%SumGDDcuts = GetSumGDDcuts()
                AC71_struc(n)%ac71(t)%SumInterval = GetSumInterval()
                AC71_struc(n)%ac71(t)%SumKcTop = GetSumKcTop()
                AC71_struc(n)%ac71(t)%SumKcTopStress = GetSumKcTopStress()
                AC71_struc(n)%ac71(t)%SumKci = GetSumKci()
                AC71_struc(n)%ac71(t)%SumWaBal = GetSumWaBal()
                AC71_struc(n)%ac71(t)%SurfaceStorage = GetSurfaceStorage()
                AC71_struc(n)%ac71(t)%Tact = GetTact()
                AC71_struc(n)%ac71(t)%TactWeedInfested = GetTactWeedInfested()
                AC71_struc(n)%ac71(t)%Tadj = GetTadj()
                AC71_struc(n)%ac71(t)%TheProjectType = TheProjectType
                AC71_struc(n)%ac71(t)%TimeSenescence = GetTimeSenescence ()
                AC71_struc(n)%ac71(t)%Tpot = GetTpot()
                AC71_struc(n)%ac71(t)%WeedRCi = GetWeedRCi()
                AC71_struc(n)%ac71(t)%YprevSum = GetYprevSum()
                AC71_struc(n)%ac71(t)%Ziprev = GetZiprev()
                AC71_struc(n)%ac71(t)%alfaHI = GetalfaHI()
                AC71_struc(n)%ac71(t)%alfaHIAdj = GetalfaHIAdj()
                AC71_struc(n)%ac71(t)%crop = GetCrop()
                AC71_struc(n)%ac71(t)%daynri = GetDayNri()

                if ((LIS_rc%mo .eq. AC71_struc(n)%Sim_AnnualEndMonth) &
                    .and.(LIS_rc%da .eq. AC71_struc(n)%Sim_AnnualEndDay)) then
                    AC71_struc(n)%ac71(t)%irun = 2 ! Means that we need to start a new sim
                    AC71_struc(n)%ac71(t)%InitializeRun = 1
                endif
        enddo ! do t = 1, LIS_rc%npatch(n, mtype)
    enddo
end subroutine Ac71_setup

!BOP
!
! !ROUTINE: AC71_read_MULTILEVEL_param
!  \label{AC71_read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
!  30 Oct  2013: Shugong Wang; Generalization for reading MULTILEVEL spatial parameter
!
! !INTERFACE:
subroutine AC71_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
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
        call LIS_verify(ios, 'Error in nf90_open in AC71_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AC71_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AC71_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC71_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC71_read_MULTILEVEL_param')

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
        call LIS_verify(ios, 'Error in nf90_get_var in AC71_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in AC71_read_MULTILEVEL_param')

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
 end subroutine AC71_read_MULTILEVEL_param

!BOP
!
! !ROUTINE: Ac71_read_croptype
!  \label{Ac71_read_croptype}
!
! !REVISION HISTORY:
!  06 MAR 2024; Louise Busschaert, initial implementation
!
! !INTERFACE:
subroutine ac71_read_croptype(n)
! !USES
    use ESMF
    use LIS_fileIOMod
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_coreMod
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use Ac71_lsmMod
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

    call ESMF_ConfigFindLabel(LIS_config,"AquaCrop.7.1 crop library directory:",rc=rc)
    call ESMF_ConfigGetAttribute(LIS_config,crop_path,rc=rc)
    call LIS_verify(rc,'AquaCrop.7.1 crop library directory: not specified')

 !- Set path to crop files
      AC71_struc(n)%PathCropFiles = trim(crop_path)//"AC_Crop.Files/"

 !- Read in crop table
      open(19, FILE=trim(crop_path)//"AC_Crop.Inventory",FORM='FORMATTED',STATUS='OLD',IOSTAT=rc)
      call LIS_verify(rc,'Ac71_setup.F: failure opening AC_Crop.Inventory')

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
       call LIS_verify(ios,'nf90_close failed in AC71_setup')
       
       l_croptype(:,:) = glb_croptype(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1))
       deallocate(glb_croptype)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          
          AC71_struc(n)%ac71(t)%cropt = trim(croptypes(nint(l_croptype(col,row))))
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

  end subroutine ac71_read_croptype
                                          

