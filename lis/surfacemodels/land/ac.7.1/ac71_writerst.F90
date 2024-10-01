!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: Ac71_writerst
! \label{Ac71_writerst}
!
! !REVISION HISTORY:
!   26 FEB 2024: Louise Busschaert; initial implementation for AC71
!
! !INTERFACE:
subroutine Ac71_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use Ac71_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for Ac71.
!  This includes all relevant AC71 variables to restart.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[Ac71\_dump\_restart](\ref{Ac71_dump_restart}) \newline
!   writes the Ac71 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "Ac71 restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(AC71_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "AC71",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in Ac71_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in Ac71_writerst")
#endif
             endif
        endif
    
        call Ac71_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in Ac71_writerst")
#endif
            endif
            write(LIS_logunit, *) "Ac71 archive restart written: ", trim(filen)
        endif
    endif
end subroutine Ac71_writerst

!BOP
!
! !ROUTINE: Ac71_dump_restart
! \label{Ac71_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  9/4/14: Shugong Wang, initial implementation for LIS 7 and Ac71
! !INTERFACE:
subroutine Ac71_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use Ac71_lsmMod

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
    integer :: StressLeaf_ID
    integer :: StressSenescence_ID
    integer :: SumGDDcuts_ID
    integer :: SumKci_ID
    integer :: SumKcTopStress_ID
    integer :: SumWaBal_Biomass_ID
    integer :: SumWaBal_BiomassPot_ID
    integer :: SumWaBal_BiomassTot_ID
    integer :: SumWaBal_BiomassUnlim_ID
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
    integer :: DayNri_ID
    integer :: DaySubmerged_ID
    integer :: Management_WeedDeltaRC_ID
    integer :: PreviousStressLevel_ID
    integer :: Simulation_DayAnaero_ID
    integer :: Simulation_EffectStress_RedCGC_ID
    integer :: Simulation_EffectStress_RedCCx_ID
    integer :: Simulation_EffectStress_RedWP_ID
    integer :: Simulation_EffectStress_RedKsSto_ID
    integer :: Simulation_EvapStartStg2_ID
    integer :: StressSFadjNEW_ID
    integer :: SumInterval_ID

    !! logicals
    integer :: NoMoreCrop_ID
    integer :: Simulation_EvapLimitON_ID
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
         "AC71", &
         dim1 = AC71_struc(n)%max_No_Compartments, &
         ! dim2 = AC71_struc(n)%max_No_Compartments*11, & ! number of cells in AC7.1
         dimID = dimID, &
         output_format = trim(wformat))

    ! write the header for state variable smc
    call LIS_writeHeader_restart(ftn, n, dimID, smc_ID, "SMC", &
                                 "volumtric soil moisture", &
                                 "m3/m3", vlevels=AC71_struc(n)%max_No_Compartments , valid_min=-99999.0, valid_max=99999.0, &
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
    ! write the header for state variable Simulation_EvapLimitON
    call LIS_writeHeader_restart(ftn, n, dimID, Simulation_EvapLimitON_ID, "Simulation_EvapLimitON", &
                           "Simulation_EvapLimitON at last time step", &
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
    
    !! logicals
    ! write the header for state variable NoMoreCrop
    call LIS_writeHeader_restart(ftn, n, dimID, NoMoreCrop_ID, "NoMoreCrop", &
                            "NoMoreCrop at last time step", &
                            "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    !! from Compartment
    ! write the header for state variable Compartment_Salt

    ! write the header for state variable Compartment_Depo

    ! write the header for state variable Compartment_fluxout
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_fluxout_ID, "Compartment_fluxout", &
                            "Compartment_fluxout at last time step", &
                            "-", vlevels=AC71_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_Smax
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_Smax_ID, "Compartment_Smax", &
                            "Compartment_Smax at last time step", &
                            "-", vlevels=AC71_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_FCadj
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_FCadj_ID, "Compartment_FCadj", &
                            "Compartment_FCadj at last time step", &
                            "-", vlevels=AC71_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_WFactor
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_WFactor_ID, "Compartment_WFactor", &
                            "Compartment_WFactor at last time step", &
                            "-", vlevels=AC71_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
                            var_flag = "dim1")

    ! write the header for state variable Compartment_DayAnaero
    call LIS_writeHeader_restart(ftn, n, dimID, Compartment_DayAnaero_ID, "Compartment_DayAnaero", &
                            "Compartment_DayAnaero at last time step", &
                            "-", vlevels=AC71_struc(n)%max_No_Compartments, valid_min=-99999.0, valid_max=99999.0, &
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
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, AC71_struc(n)%rstInterval)


end subroutine Ac71_dump_restart
