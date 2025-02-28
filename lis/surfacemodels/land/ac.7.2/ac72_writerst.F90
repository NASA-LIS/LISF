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
! !ROUTINE: AC72_writerst
! \label{AC72_writerst}
!
! !REVISION HISTORY:
!   04 NOV 2024: Louise Busschaert; initial implementation for AC72
!
! !INTERFACE:
subroutine AC72_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use AC72_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for AC72.
!  This includes all relevant AC72 variables to restart.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[AC72\_dump\_restart](\ref{AC72_dump_restart}) \newline
!   writes the AC72 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck, alarmCheck_sf
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "AC72 restart alarm")
    alarmCheck_sf = LIS_isAlarmRinging(LIS_rc, "AC72 model alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(AC72_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1) .or. &
      ((AC72_struc(n)%ac72(1)%InitializeRun.eq.1).and.alarmCheck_sf)) then
    ! Writes restart file at the end of AquaCrop simulation period
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "AC72",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in AC72_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in AC72_writerst")
#endif
             endif
        endif
    
        call AC72_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in AC72_writerst")
#endif
            endif
            write(LIS_logunit, *) "AC72 archive restart written: ", trim(filen)
        endif
    endif
end subroutine AC72_writerst

!BOP
!
! !ROUTINE: AC72_dump_restart
! \label{AC72_dump_restart}
!
! !REVISION HISTORY:
!  04 NOV 2024: Louise Busschaert; initial implementation for AC72
! !INTERFACE:
subroutine AC72_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use AC72_lsmMod

    implicit none

    integer, intent(in) :: ftn
    integer, intent(in) :: n
    character(len=*), intent(in) :: wformat
!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!   \item[wformat]
!    restart file format (binary/netcdf)
!  \end{description}
!
!
! The routines invoked are:
! \begin{description}
!   \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart}) \newline
!      writes the global header information
!   \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart}) \newline
!      writes the header information for a variable
!   \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart}) \newline
!      close the header
!   \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart}) \newline
!      writes a variable to the restart file
! \end{description}
! 
!EOP 
               
    integer :: l, t 
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: tmptilen_int(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: dimID(11)

    integer :: SMC_ID

    !! reals
    integer :: alfaHI_ID
    integer :: alfaHIAdj_ID
    integer :: Bin_ID
    integer :: Bout_ID
    integer :: CCiActual_ID
    integer :: CCiActualWeedInfested_ID
    integer :: CCiPrev_ID
    integer :: CCiTopEarlySen_ID
    integer :: CCxWitheredTpotNoS_ID
    integer :: Crop_CCoAdjusted_ID
    integer :: Crop_CCxAdjusted_ID
    integer :: Crop_CCxWithered_ID
    integer :: Crop_pActStom_ID
    integer :: Crop_pLeafAct_ID
    integer :: Crop_pSenAct_ID
    integer :: DayFraction_ID
    integer :: ECstorage_ID
    integer :: HItimesAT_ID
    integer :: HItimesAT1_ID
    integer :: HItimesAT2_ID
    integer :: HItimesBEF_ID
    integer :: PlotVarCrop_ActVal_ID
    integer :: PlotVarCrop_PotVal_ID
    integer :: RootZoneWC_Actual_ID
    integer :: RootZoneWC_FC_ID
    integer :: RootZoneWC_Leaf_ID
    integer :: RootZoneWC_SAT_ID
    integer :: RootZoneWC_Sen_ID
    integer :: RootZoneWC_Thresh_ID
    integer :: RootZoneWC_WP_ID
    integer :: RootZoneWC_ZtopAct_ID
    integer :: RootZoneWC_ZtopFC_ID
    integer :: RootZoneWC_ZtopThresh_ID
    integer :: RootZoneWC_ZtopWP_ID
    integer :: ScorAT1_ID
    integer :: ScorAT2_ID
    integer :: Simulation_EffectStress_CDecline_ID
    integer :: Simulation_EvapWCSurf_ID
    integer :: Simulation_SumEToStress_ID
    integer :: Simulation_Scor_ID
    integer :: Simulation_SumGDD_ID
    integer :: Simulation_SumGDDfromDay1_ID
    integer :: StressLeaf_ID
    integer :: StressSenescence_ID
    integer :: SumGDDcuts_ID
    integer :: SumGDDPrev_ID
    integer :: SumKci_ID
    integer :: SumKcTopStress_ID
    integer :: SumWaBal_Biomass_ID
    integer :: SumWaBal_BiomassPot_ID
    integer :: SumWaBal_BiomassTot_ID
    integer :: SumWaBal_BiomassUnlim_ID
    integer :: SumWaBal_Tact_ID
    integer :: SumWaBal_Tpot_ID
    integer :: SumWaBal_YieldPart_ID
    integer :: SurfaceStorage_ID
    integer :: Tact_ID
    integer :: TactWeedInfested_ID
    integer :: Tadj_ID
    integer :: TimeSenescence_ID
    integer :: Tpot_ID
    integer :: WeedRCi_ID
    integer :: WPi_ID
    integer :: Ziprev_ID

    !! integers
    integer :: Crop_DaysToFullCanopySF_ID
    integer :: Crop_GDDaysToFullCanopySF_ID
    integer :: DayLastCut_ID
    integer :: DayNri_ID
    integer :: DaySubmerged_ID
    integer :: Management_WeedDeltaRC_ID
    integer :: PreviousStressLevel_ID
    integer :: Simulation_DayAnaero_ID
    integer :: Simulation_DelayedDays_ID
    integer :: Simulation_EffectStress_RedCGC_ID
    integer :: Simulation_EffectStress_RedCCx_ID
    integer :: Simulation_EffectStress_RedWP_ID
    integer :: Simulation_EffectStress_RedKsSto_ID
    integer :: Simulation_EvapStartStg2_ID
    integer :: Simulation_HIfinal_ID
    integer :: StressSFadjNEW_ID
    integer :: SumInterval_ID

    !! logicals
    integer :: NoMoreCrop_ID
    integer :: Simulation_EvapLimitON_ID
    integer :: Simulation_Germinate_ID
    integer :: Simulation_ProtectedSeedling_ID
    integer :: Simulation_SWCtopSoilConsidered_ID

    !! From derived types
    !!! Compartment (reals)
    !integer :: Compartment_Salt_ID
    !integer :: Compartment_Depo_ID
    integer :: Compartment_fluxout_ID
    integer :: Compartment_Smax_ID
    integer :: Compartment_FCadj_ID
    integer :: Compartment_WFactor_ID
    !!! Compartment (integers)
    integer :: Compartment_DayAnaero_ID

    !!! StressTot (reals)
    integer :: StressTot_Salt_ID
    integer :: StressTot_Temp_ID
    integer :: StressTot_Exp_ID
    integer :: StressTot_Sto_ID
    integer :: StressTot_Weed_ID
    !!! StressTot (integers)
    integer :: StressTot_NrD_ID


    ! write the header of the restart file
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
         "AC72", &
         dim1 = AC72_struc(n)%max_No_Compartments, &
         dimID = dimID, &
         output_format = trim(wformat))

    ! write the header for state variable smc
    call LIS_writeHeader_restart(ftn, n, dimID, smc_ID, "SMC", &
                                 "volumtric soil moisture", &
                                 "m3/m3", vlevels=AC72_struc(n)%max_No_Compartments , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
    !! reals
    ! write the header for state variable alfaHI
    call LIS_writeHeader_restart(ftn, n, dimID, alfaHI_ID, "alfaHI", &
                            "alfaHI at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable alfaHIAdj
    call LIS_writeHeader_restart(ftn, n, dimID, alfaHIAdj_ID, "alfaHIAdj", &
                            "alfaHIAdj at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Bin
    call LIS_writeHeader_restart(ftn, n, dimID, Bin_ID, "Bin", &
                            "Bin at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Bout
    call LIS_writeHeader_restart(ftn, n, dimID, Bout_ID, "Bout", &
                            "Bout at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable CCiActual
    call LIS_writeHeader_restart(ftn, n, dimID, CCiActual_ID, "CCiActual", &
                            "CCiActual at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable CCiActualWeedInfested
    call LIS_writeHeader_restart(ftn, n, dimID, CCiActualWeedInfested_ID, "CCiActualWeedInfested", &
                            "CCiActualWeedInfested at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable CCiPrev
    call LIS_writeHeader_restart(ftn, n, dimID, CCiPrev_ID, "CCiPrev", &
                            "CCiPrev at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable CCiTopEarlySen
    call LIS_writeHeader_restart(ftn, n, dimID, CCiTopEarlySen_ID, "CCiTopEarlySen", &
                            "CCiTopEarlySen at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable CCxWitheredTpotNoS
    call LIS_writeHeader_restart(ftn, n, dimID, CCxWitheredTpotNoS_ID, "CCxWitheredTpotNoS", &
                            "CCxWitheredTpotNoS at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Crop_CCoAdjusted
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_CCoAdjusted_ID, "Crop_CCoAdjusted", &
                            "Crop_CCoAdjusted at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Crop_CCxAdjusted
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_CCxAdjusted_ID, "Crop_CCxAdjusted", &
                            "Crop_CCxAdjusted at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Crop_CCxWithered
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_CCxWithered_ID, "Crop_CCxWithered", &
                            "Crop_CCxWithered at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Crop_pActStom
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_pActStom_ID, "Crop_pActStom", &
                            "Crop_pActStom at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Crop_pLeafAct
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_pLeafAct_ID, "Crop_pLeafAct", &
                            "Crop_pLeafAct at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Crop_pSenAct
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_pSenAct_ID, "Crop_pSenAct", &
                            "Crop_pSenAct at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable DayFraction
    call LIS_writeHeader_restart(ftn, n, dimID, DayFraction_ID, "DayFraction", &
                            "DayFraction at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable ECstorage
    call LIS_writeHeader_restart(ftn, n, dimID, ECstorage_ID, "ECstorage", &
                            "ECstorage at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable HItimesAT
    call LIS_writeHeader_restart(ftn, n, dimID, HItimesAT_ID, "HItimesAT", &
                            "HItimesAT at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable HItimesAT1
    call LIS_writeHeader_restart(ftn, n, dimID, HItimesAT1_ID, "HItimesAT1", &
                            "HItimesAT1 at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable HItimesAT2
    call LIS_writeHeader_restart(ftn, n, dimID, HItimesAT2_ID, "HItimesAT2", &
                            "HItimesAT2 at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable HItimesBEF
    call LIS_writeHeader_restart(ftn, n, dimID, HItimesBEF_ID, "HItimesBEF", &
                            "HItimesBEF at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable PlotVarCrop_ActVal
    call LIS_writeHeader_restart(ftn, n, dimID, PlotVarCrop_ActVal_ID, "PlotVarCrop_ActVal", &
                            "PlotVarCrop_ActVal at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable PlotVarCrop_PotVal
    call LIS_writeHeader_restart(ftn, n, dimID, PlotVarCrop_PotVal_ID, "PlotVarCrop_PotVal", &
                            "PlotVarCrop_PotVal at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_Actual
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_Actual_ID, "RootZoneWC_Actual", &
                            "RootZoneWC_Actual at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_FC
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_FC_ID, "RootZoneWC_FC", &
                            "RootZoneWC_FC at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_Leaf
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_Leaf_ID, "RootZoneWC_Leaf", &
                            "RootZoneWC_Leaf at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_SAT
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_SAT_ID, "RootZoneWC_SAT", &
                            "RootZoneWC_SAT at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_Sen
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_Sen_ID, "RootZoneWC_Sen", &
                            "RootZoneWC_Sen at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_Thresh
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_Thresh_ID, "RootZoneWC_Thresh", &
                            "RootZoneWC_Thresh at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_WP
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_WP_ID, "RootZoneWC_WP", &
                            "RootZoneWC_WP at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_ZtopAct
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_ZtopAct_ID, "RootZoneWC_ZtopAct", &
                            "RootZoneWC_ZtopAct at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_ZtopFC
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_ZtopFC_ID, "RootZoneWC_ZtopFC", &
                            "RootZoneWC_ZtopFC at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_ZtopThresh
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_ZtopThresh_ID, "RootZoneWC_ZtopThresh", &
                            "RootZoneWC_ZtopThresh at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RootZoneWC_ZtopWP
    call LIS_writeHeader_restart(ftn, n, dimID, RootZoneWC_ZtopWP_ID, "RootZoneWC_ZtopWP", &
                            "RootZoneWC_ZtopWP at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable ScorAT1
    call LIS_writeHeader_restart(ftn, n, dimID, ScorAT1_ID, "ScorAT1", &
                            "ScorAT1 at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable ScorAT2
    call LIS_writeHeader_restart(ftn, n, dimID, ScorAT2_ID, "ScorAT2", &
                            "ScorAT2 at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EffectStress_CDecline
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EffectStress_CDecline_ID, "Simulation_EffectStress_CDecline", &
                           "Simulation_EffectStress_CDecline at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EvapWCSurf
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EvapWCSurf_ID, "Simulation_EvapWCSurf", &
                           "Simulation_EvapWCSurf at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_SCor
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_SCor_ID, "Simulation_SCor", &
                           "Simulation_SCor at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_SumEToStress
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_SumEToStress_ID, "Simulation_SumEToStress", &
                           "Simulation_SumEToStress at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EvapLimitON
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EvapLimitON_ID, "Simulation_EvapLimitON", &
                           "Simulation_EvapLimitON at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_Germinate
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_Germinate_ID, "Simulation_Germinate", &
                           "Simulation_Germinate at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_ProtectedSeedling
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_ProtectedSeedling_ID, "Simulation_ProtectedSeedling", &
                           "Simulation_ProtectedSeedling at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_SumGDD
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_SumGDD_ID, "Simulation_SumGDD", &
                           "Simulation_SumGDD at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_SumGDDFromDay1
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_SumGDDFromDay1_ID, "Simulation_SumGDDFromDay1", &
                           "Simulation_SumGDDFromDay1 at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_SWCtopSoilConsidered
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_SWCtopSoilConsidered_ID, "Simulation_SWCtopSoilConsidered", &
                           "Simulation_SWCtopSoilConsidered at last time step", &
                           "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable StressLeaf
    call LIS_writeHeader_restart(ftn, n, dimID, StressLeaf_ID, "StressLeaf", &
                            "StressLeaf at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable StressSenescence
    call LIS_writeHeader_restart(ftn, n, dimID, StressSenescence_ID, "StressSenescence", &
                            "StressSenescence at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumGDDcuts
    call LIS_writeHeader_restart(ftn, n, dimID, SumGDDcuts_ID, "SumGDDcuts", &
                            "SumGDDcuts at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumGDDPrev
    call LIS_writeHeader_restart(ftn, n, dimID, SumGDDPrev_ID, "SumGDDPrev", &
                            "SumGDDPrev at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumKci
    call LIS_writeHeader_restart(ftn, n, dimID, SumKci_ID, "SumKci", &
                            "SumKci at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumKcTopStress
    call LIS_writeHeader_restart(ftn, n, dimID, SumKcTopStress_ID, "SumKcTopStress", &
                            "SumKcTopStress at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumWaBal_Biomass
    call LIS_writeHeader_restart(ftn, n, dimID, SumWaBal_Biomass_ID, "SumWaBal_Biomass", &
                            "SumWaBal_Biomass at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumWaBal_BiomassPot
    call LIS_writeHeader_restart(ftn, n, dimID, SumWaBal_BiomassPot_ID, "SumWaBal_BiomassPot", &
                            "SumWaBal_BiomassPot at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumWaBal_BiomassTot
    call LIS_writeHeader_restart(ftn, n, dimID, SumWaBal_BiomassTot_ID, "SumWaBal_BiomassTot", &
                            "SumWaBal_BiomassTot at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumWaBal_BiomassUnlim
    call LIS_writeHeader_restart(ftn, n, dimID, SumWaBal_BiomassUnlim_ID, "SumWaBal_BiomassUnlim", &
                            "SumWaBal_BiomassUnlim at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumWaBal_Tact
    call LIS_writeHeader_restart(ftn, n, dimID, SumWaBal_Tact_ID, "SumWaBal_Tact", &
                            "SumWaBal_Tact at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumWaBal_Tpot
    call LIS_writeHeader_restart(ftn, n, dimID, SumWaBal_Tpot_ID, "SumWaBal_Tpot", &
                            "SumWaBal_Tpot at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumWaBal_YieldPart
    call LIS_writeHeader_restart(ftn, n, dimID, SumWaBal_YieldPart_ID, "SumWaBal_YieldPart", &
                            "SumWaBal_YieldPart at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SurfaceStorage
    call LIS_writeHeader_restart(ftn, n, dimID, SurfaceStorage_ID, "SurfaceStorage", &
                            "SurfaceStorage at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Tact
    call LIS_writeHeader_restart(ftn, n, dimID, Tact_ID, "Tact", &
                            "Tact at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable TactWeedInfested
    call LIS_writeHeader_restart(ftn, n, dimID, TactWeedInfested_ID, "TactWeedInfested", &
                            "TactWeedInfested at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Tadj
    call LIS_writeHeader_restart(ftn, n, dimID, Tadj_ID, "Tadj", &
                            "Tadj at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable TimeSenescence
    call LIS_writeHeader_restart(ftn, n, dimID, TimeSenescence_ID, "TimeSenescence", &
                            "TimeSenescence at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Tpot
    call LIS_writeHeader_restart(ftn, n, dimID, Tpot_ID, "Tpot", &
                            "Tpot at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable WeedRCi
    call LIS_writeHeader_restart(ftn, n, dimID, WeedRCi_ID, "WeedRCi", &
                            "WeedRCi at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable WPi
    call LIS_writeHeader_restart(ftn, n, dimID, WPi_ID, "WPi", &
                            "WPi at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Ziprev
    call LIS_writeHeader_restart(ftn, n, dimID, Ziprev_ID, "Ziprev", &
                            "Ziprev at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    

    !! integers
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_DaysToFullCanopySF_ID, "Crop_DaysToFullCanopySF", &
                            "Crop_DaysToFullCanopySF at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    call LIS_writeHeader_restart(ftn, n, dimID, Crop_GDDaysToFullCanopySF_ID, "Crop_GDDaysToFullCanopySF", &
                            "Crop_GDDaysToFullCanopySF at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable DayLastCut
    call LIS_writeHeader_restart(ftn, n, dimID, DayLastCut_ID, "DayLastCut", &
                            "DayLastCut at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable DayNri
    call LIS_writeHeader_restart(ftn, n, dimID, DayNri_ID, "DayNri", &
                            "DayNri at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable DaySubmerged
    call LIS_writeHeader_restart(ftn, n, dimID, DaySubmerged_ID, "DaySubmerged", &
                            "DaySubmerged at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable PreviousStressLevel
    call LIS_writeHeader_restart(ftn, n, dimID, PreviousStressLevel_ID, "PreviousStressLevel", &
                            "PreviousStressLevel at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable StressSFadjNEW
    call LIS_writeHeader_restart(ftn, n, dimID, StressSFadjNEW_ID, "StressSFadjNEW", &
                            "StressSFadjNEW at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SumInterval
    call LIS_writeHeader_restart(ftn, n, dimID, SumInterval_ID, "SumInterval", &
                            "SumInterval at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Management_WeedDeltaRC
    call LIS_writeHeader_restart(ftn, n, dimID, Management_WeedDeltaRC_ID, "Management_WeedDeltaRC", &
                            "Management_WeedDeltaRC at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_DayAnaero
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_DayAnaero_ID, "Simulation_DayAnaero", &
                            "Simulation_DayAnaero at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_DelayedDays
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_DelayedDays_ID, "Simulation_DelayedDays", &
                            "Simulation_DelayedDays at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EffectStress_RedCGC
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EffectStress_RedCGC_ID, "Simulation_EffectStress_RedCGC", &
                            "Simulation_EffectStress_RedCGC at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EffectStress_RedCCx
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EffectStress_RedCCx_ID, "Simulation_EffectStress_RedCCx", &
                            "Simulation_EffectStress_RedCCx at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EffectStress_RedWP
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EffectStress_RedWP_ID, "Simulation_EffectStress_RedWP", &
                            "Simulation_EffectStress_RedWP at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EffectStress_RedKsSto
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EffectStress_RedKsSto_ID, "Simulation_EffectStress_RedKsSto", &
                            "Simulation_EffectStress_RedKsSto at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_EvapStartStg2
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EvapStartStg2_ID, "Simulation_EvapStartStg2", &
                            "Simulation_EvapStartStg2 at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable Simulation_HIfinal
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_HIfinal_ID, "Simulation_HIfinal", &
                            "Simulation_HIfinal at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    
    !! logicals
    ! write the header for state variable NoMoreCrop
    call LIS_writeHeader_restart(ftn, n, dimID, NoMoreCrop_ID, "NoMoreCrop", &
                            "NoMoreCrop at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    !! from Compartment

    ! Restart not implemented for salinity
    ! write the header for state variable Compartment_Salt
    ! write the header for state variable Compartment_Depo

    ! write the header for state variable Compartment_fluxout
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_fluxout_ID, "Compartment_fluxout", &
                            "Compartment_fluxout at last time step", &
                            "-", vlevels=AC72_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_Smax
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_Smax_ID, "Compartment_Smax", &
                            "Compartment_Smax at last time step", &
                            "-", vlevels=AC72_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_FCadj
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_FCadj_ID, "Compartment_FCadj", &
                            "Compartment_FCadj at last time step", &
                            "-", vlevels=AC72_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_WFactor
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_WFactor_ID, "Compartment_WFactor", &
                            "Compartment_WFactor at last time step", &
                            "-", vlevels=AC72_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_DayAnaero
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_DayAnaero_ID, "Compartment_DayAnaero", &
                            "Compartment_DayAnaero at last time step", &
                            "-", vlevels=AC72_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! From StressTot
    ! write the header for state variable StressTot_Salt
    call LIS_writeHeader_restart(ftn, n, dimID, StressTot_Salt_ID, "StressTot_Salt", &
                            "StressTot_Salt at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable StressTot_Temp
    call LIS_writeHeader_restart(ftn, n, dimID, StressTot_Temp_ID, "StressTot_Temp", &
                            "StressTot_Temp at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable StressTot_Exp
    call LIS_writeHeader_restart(ftn, n, dimID, StressTot_Exp_ID, "StressTot_Exp", &
                            "StressTot_Exp at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable StressTot_Sto
    call LIS_writeHeader_restart(ftn, n, dimID, StressTot_Sto_ID, "StressTot_Sto", &
                            "StressTot_Sto at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable StressTot_Weed
    call LIS_writeHeader_restart(ftn, n, dimID, StressTot_Weed_ID, "StressTot_Weed", &
                            "StressTot_Weed at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable StressTot_NrD
    call LIS_writeHeader_restart(ftn, n, dimID, StressTot_NrD_ID, "StressTot_NrD", &
                            "StressTot_NrD at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, AC72_struc(n)%rstInterval)

    ! write state variables into restart file
    ! volumtric soil moisture
    do l=1, AC72_struc(n)%max_No_Compartments
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AC72_struc(n)%ac72(t)%smc(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=smc_ID, dim=l, wformat=wformat)
    enddo

    !! From Compartment
    ! Compartment_fluxout
    do l=1, AC72_struc(n)%max_No_Compartments
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AC72_struc(n)%ac72(t)%Compartment(l)%fluxout
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=Compartment_fluxout_ID, dim=l, wformat=wformat)
    enddo

    ! Compartment_Smax
    do l=1, AC72_struc(n)%max_No_Compartments
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AC72_struc(n)%ac72(t)%Compartment(l)%Smax
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=Compartment_Smax_ID, dim=l, wformat=wformat)
    enddo

    ! Compartment_FCadj
    do l=1, AC72_struc(n)%max_No_Compartments
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AC72_struc(n)%ac72(t)%Compartment(l)%FCadj
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=Compartment_FCadj_ID, dim=l, wformat=wformat)
    enddo

    ! Compartment_WFactor
    do l=1, AC72_struc(n)%max_No_Compartments
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AC72_struc(n)%ac72(t)%Compartment(l)%WFactor
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=Compartment_WFactor_ID, dim=l, wformat=wformat)
    enddo

    ! Compartment_DayAnaero
    do l=1, AC72_struc(n)%max_No_Compartments
        tmptilen_int = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen_int(t) = AC72_struc(n)%ac72(t)%Compartment(l)%DayAnaero
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                  varid=Compartment_DayAnaero_ID, dim=l, wformat=wformat)
    enddo



    !! reals
    ! alfaHI
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%alfaHI, &
                            varid=alfaHI_ID, dim=1, wformat=wformat)

    ! alfaHIAdj
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%alfaHIAdj, &
                            varid=alfaHIAdj_ID, dim=1, wformat=wformat)

    ! Bin 
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Bin, &
                            varid=Bin_ID, dim=1, wformat=wformat)

    ! Bout
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Bout, &
                            varid=Bout_ID, dim=1, wformat=wformat)

    ! CCiActual
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%CCiActual, &
                            varid=CCiActual_ID, dim=1, wformat=wformat)

    ! CCiActualWeedInfested
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%CCiActualWeedInfested, &
                            varid=CCiActualWeedInfested_ID, dim=1, wformat=wformat)

    ! CCiPrev
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%CCiPrev, &
                            varid=CCiPrev_ID, dim=1, wformat=wformat)

    ! CCiTopEarlySen
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%CCiTopEarlySen, &
                            varid=CCiTopEarlySen_ID, dim=1, wformat=wformat)

    ! CCxWitheredTpotNoS
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%CCxWitheredTpotNoS, &
                            varid=CCxWitheredTpotNoS_ID, dim=1, wformat=wformat)

    ! DayFraction
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%DayFraction, &
                            varid=DayFraction_ID, dim=1, wformat=wformat)

    ! ECstorage
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%ECstorage, &
                            varid=ECstorage_ID, dim=1, wformat=wformat)

    ! HItimesAT
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%HItimesAT, &
                            varid=HItimesAT_ID, dim=1, wformat=wformat)

    ! HItimesAT1
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%HItimesAT1, &
                            varid=HItimesAT1_ID, dim=1, wformat=wformat)

    ! HItimesAT2
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%HItimesAT2, &
                            varid=HItimesAT2_ID, dim=1, wformat=wformat)

    ! HItimesBEF
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%HItimesBEF, &
                            varid=HItimesBEF_ID, dim=1, wformat=wformat)

    ! RootZoneWC_Actual
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_Actual, &
                            varid=RootZoneWC_Actual_ID, dim=1, wformat=wformat)

    ! RootZoneWC_FC
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_FC, &
                            varid=RootZoneWC_FC_ID, dim=1, wformat=wformat)

    ! RootZoneWC_Leaf
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_Leaf, &
                            varid=RootZoneWC_Leaf_ID, dim=1, wformat=wformat)

    ! RootZoneWC_SAT
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_SAT, &
                            varid=RootZoneWC_SAT_ID, dim=1, wformat=wformat)

    ! RootZoneWC_Sen
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_Sen, &
                            varid=RootZoneWC_Sen_ID, dim=1, wformat=wformat)

    ! RootZoneWC_Thresh
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_Thresh, &
                            varid=RootZoneWC_Thresh_ID, dim=1, wformat=wformat)

    ! RootZoneWC_WP
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_WP, &
                            varid=RootZoneWC_WP_ID, dim=1, wformat=wformat)

    ! RootZoneWC_ZtopAct
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_ZtopAct, &
                            varid=RootZoneWC_ZtopAct_ID, dim=1, wformat=wformat)

    ! RootZoneWC_ZtopFC
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_ZtopFC, &
                            varid=RootZoneWC_ZtopFC_ID, dim=1, wformat=wformat)

    ! RootZoneWC_ZtopThresh
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_ZtopThresh, &
                            varid=RootZoneWC_ZtopThresh_ID, dim=1, wformat=wformat)

    ! RootZoneWC_ZtopWP
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%RootZoneWC_ZtopWP, &
                            varid=RootZoneWC_ZtopWP_ID, dim=1, wformat=wformat)

    ! ScorAT1
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%ScorAT1, &
                            varid=ScorAT1_ID, dim=1, wformat=wformat)

    ! ScorAT2
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%ScorAT2, &
                            varid=ScorAT2_ID, dim=1, wformat=wformat)

    ! StressLeaf
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%StressLeaf, &
                            varid=StressLeaf_ID, dim=1, wformat=wformat)

    ! StressSenescence
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%StressSenescence, &
                            varid=StressSenescence_ID, dim=1, wformat=wformat)

    ! SumGDDcuts
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%SumGDDcuts, &
                            varid=SumGDDcuts_ID, dim=1, wformat=wformat)

    ! SumGDDPrev
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%SumGDDPrev, &
                            varid=SumGDDPrev_ID, dim=1, wformat=wformat)

    ! SumKci
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%SumKci, &
                            varid=SumKci_ID, dim=1, wformat=wformat)

    ! SumKcTopStress
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%SumKcTopStress, &
                            varid=SumKcTopStress_ID, dim=1, wformat=wformat)

    ! SurfaceStorage
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%SurfaceStorage, &
                            varid=SurfaceStorage_ID, dim=1, wformat=wformat)

    ! Tact
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Tact, &
                            varid=Tact_ID, dim=1, wformat=wformat)

    ! TactWeedInfested
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%TactWeedInfested, &
                            varid=TactWeedInfested_ID, dim=1, wformat=wformat)

    ! Tadj
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Tadj, &
                            varid=Tadj_ID, dim=1, wformat=wformat)

    ! TimeSenescence
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%TimeSenescence, &
                            varid=TimeSenescence_ID, dim=1, wformat=wformat)

    ! Tpot
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Tpot, &
                            varid=Tpot_ID, dim=1, wformat=wformat)

    ! WeedRCi
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%WeedRCi, &
                            varid=WeedRCi_ID, dim=1, wformat=wformat)

    ! WPi
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%WPi, &
                            varid=WPi_ID, dim=1, wformat=wformat)

    ! Ziprev
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Ziprev, &
                            varid=Ziprev_ID, dim=1, wformat=wformat)

    
    !! integers
    ! Crop_DaysToFullCanopySF
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Crop%DaysToFullCanopySF, &
                            varid=Crop_DaysToFullCanopySF_ID, dim=1, wformat=wformat)

    ! Crop_GDDaysToFullCanopySF
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Crop%GDDaysToFullCanopySF, &
                            varid=Crop_GDDaysToFullCanopySF_ID, dim=1, wformat=wformat)
    ! DayLastCut
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%DayLastCut, &
                            varid=DayLastCut_ID, dim=1, wformat=wformat)

    ! DayNri
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%DayNri, &
                            varid=DayNri_ID, dim=1, wformat=wformat)

    ! DaySubmerged
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%DaySubmerged, &
                            varid=DaySubmerged_ID, dim=1, wformat=wformat)

    ! PreviousStressLevel
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%PreviousStressLevel, &
                            varid=PreviousStressLevel_ID, dim=1, wformat=wformat)

    ! StressSFadjNEW
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%StressSFadjNEW, &
                            varid=StressSFadjNEW_ID, dim=1, wformat=wformat)

    ! SumInterval
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%SumInterval, &
                            varid=SumInterval_ID, dim=1, wformat=wformat)

    ! Management_WeedDeltaRC
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AC72_struc(n)%ac72%Management%WeedDeltaRC, &
                            varid=Management_WeedDeltaRC_ID, dim=1, wformat=wformat)


    !! logicals (convert to integer)
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        if (AC72_struc(n)%ac72(t)%NoMoreCrop) then
            tmptilen_int(t) = 1
        endif
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                            varid=NoMoreCrop_ID, dim=1, wformat=wformat)


    !! From Management
    ! Management_WeedDeltaRC
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Management%WeedDeltaRC
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Management_WeedDeltaRC_ID, dim=1, wformat=wformat)

    !! From Crop
    ! Crop_CCoAdjusted
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Crop%CCoAdjusted
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Crop_CCoAdjusted_ID, dim=1, wformat=wformat)

    ! Crop_CCxAdjusted
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Crop%CCxAdjusted
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Crop_CCxAdjusted_ID, dim=1, wformat=wformat)

    ! Crop_CCxWithered
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Crop%CCxWithered
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Crop_CCxWithered_ID, dim=1, wformat=wformat)

    ! Crop_pActStom
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Crop%pActStom
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Crop_pActStom_ID, dim=1, wformat=wformat)

    ! Crop_pSenAct
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Crop%pSenAct
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Crop_pSenAct_ID, dim=1, wformat=wformat)

    ! Crop_pLeafAct
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Crop%pLeafAct
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Crop_pLeafAct_ID, dim=1, wformat=wformat)

    !! From Simulation
    ! Simulation_EvapLimitON
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        if (AC72_struc(n)%ac72(t)%Simulation%EvapLimitON) then
            tmptilen_int(t) = 1
        endif
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_EvapLimitON_ID, dim=1, wformat=wformat)

    ! Simulation_Germinate
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        if (AC72_struc(n)%ac72(t)%Simulation%Germinate) then
            tmptilen_int(t) = 1
        endif
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_Germinate_ID, dim=1, wformat=wformat)

    ! Simulation_ProtectedSeedling
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        if (AC72_struc(n)%ac72(t)%Simulation%ProtectedSeedling) then
            tmptilen_int(t) = 1
        endif
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_ProtectedSeedling_ID, dim=1, wformat=wformat)

    ! Simulation_SWCtopSoilConsidered
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        if (AC72_struc(n)%ac72(t)%Simulation%SWCtopSoilConsidered) then
            tmptilen_int(t) = 1
        endif
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_SWCtopSoilConsidered_ID, dim=1, wformat=wformat)

    ! Simulation_EvapWCSurf
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Simulation%EvapWCsurf
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Simulation_EvapWCSurf_ID, dim=1, wformat=wformat)

    ! Simulation_SCor
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Simulation%SCor
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Simulation_SCor_ID, dim=1, wformat=wformat)

    ! Simulation_SumEToStress
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Simulation%SumEToStress
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Simulation_SumEToStress_ID, dim=1, wformat=wformat)

    ! Simulation_SumGDD
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Simulation%SumGDD
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Simulation_SumGDD_ID, dim=1, wformat=wformat)

    ! Simulation_SumGDDFromDay1
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Simulation%SumGDDFromDay1
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Simulation_SumGDDFromDay1_ID, dim=1, wformat=wformat)

    !! From Simulation%EffectStress
    ! Simulation_EffectStress_CDecline
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%Simulation%EffectStress%CDecline
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=Simulation_EffectStress_CDecline_ID, dim=1, wformat=wformat)

    ! Simulation_DayAnaero
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%DayAnaero
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_DayAnaero_ID, dim=1, wformat=wformat)

    ! Simulation_DelayedDays
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%DelayedDays
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_DelayedDays_ID, dim=1, wformat=wformat)

    ! Simulation_EffectStress_RedCGC
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%EffectStress%RedCGC
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_EffectStress_RedCGC_ID, dim=1, wformat=wformat)

    ! Simulation_EffectStress_RedCCx
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%EffectStress%RedCCx
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_EffectStress_RedCCx_ID, dim=1, wformat=wformat)

    ! Simulation_EffectStress_RedWP
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%EffectStress%RedWP
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_EffectStress_RedWP_ID, dim=1, wformat=wformat)

    ! Simulation_EffectStress_RedKsSto
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%EffectStress%RedKsSto
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_EffectStress_RedKsSto_ID, dim=1, wformat=wformat)

    ! Simulation_EvapStartStg2
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%EvapStartStg2
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_EvapStartStg2_ID, dim=1, wformat=wformat)


    ! Simulation_HIfinal
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%Simulation%HIfinal
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=Simulation_HIfinal_ID, dim=1, wformat=wformat)

    !! From PlotVarCrop
    ! PlotVarCrop_ActVal
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%PlotVarCrop%ActVal
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=PlotVarCrop_ActVal_ID, dim=1, wformat=wformat)

    ! PlotVarCrop_PotVal
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%PlotVarCrop%PotVal
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=PlotVarCrop_PotVal_ID, dim=1, wformat=wformat)


    !! From StressTot
    ! StressTot_Salt
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%StressTot%Salt
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=StressTot_Salt_ID, dim=1, wformat=wformat)

    ! StressTot_Temp
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%StressTot%Temp
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=StressTot_Temp_ID, dim=1, wformat=wformat)

    ! StressTot_Exp
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%StressTot%Exp
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=StressTot_Exp_ID, dim=1, wformat=wformat)
                
    ! StressTot_Sto
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%StressTot%Sto
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=StressTot_Sto_ID, dim=1, wformat=wformat)

    ! StressTot_Weed
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%StressTot%Weed
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=StressTot_Weed_ID, dim=1, wformat=wformat)

    ! StressTot_NrD
    tmptilen_int = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen_int(t) = AC72_struc(n)%ac72(t)%StressTot%NrD
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                                varid=StressTot_NrD_ID, dim=1, wformat=wformat)

    !! From SumWaBal
    ! SumWaBal_Biomass
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%SumWaBal%Biomass
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=SumWaBal_Biomass_ID, dim=1, wformat=wformat)

    ! SumWaBal_BiomassPot
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%SumWaBal%BiomassPot
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=SumWaBal_BiomassPot_ID, dim=1, wformat=wformat)

    ! SumWaBal_BiomassTot
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%SumWaBal%BiomassTot
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=SumWaBal_BiomassTot_ID, dim=1, wformat=wformat)

    ! SumWaBal_BiomassUnlim
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%SumWaBal%BiomassUnlim
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=SumWaBal_BiomassUnlim_ID, dim=1, wformat=wformat)

    ! SumWaBal_Tact
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%SumWaBal%Tact
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=SumWaBal_Tact_ID, dim=1, wformat=wformat)

    ! SumWaBal_Tpot
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%SumWaBal%Tpot
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=SumWaBal_Tpot_ID, dim=1, wformat=wformat)

    ! SumWaBal_YieldPart
    tmptilen = 0
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tmptilen(t) = AC72_struc(n)%ac72(t)%SumWaBal%YieldPart
    enddo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                varid=SumWaBal_YieldPart_ID, dim=1, wformat=wformat)



end subroutine AC72_dump_restart