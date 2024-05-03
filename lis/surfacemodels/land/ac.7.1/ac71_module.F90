!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module Ac71_module
!BOP
!
! !MODULE: Ac71_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the Ac71 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[landuse\_tbl\_name]
!     Noah model landuse parameter table. unit: -
!   !!!!!LB provide documentation and units
!   \end{description}
!
! !REVISION HISTORY:
!  06 MAR 2024; Louise Busschaert, initial implementation
!
!EOP
  
  use ac_global, only: &
          CompartmentIndividual,&
          SoilLayerIndividual, &
          rep_Crop,&
          rep_DayEventInt,&
          rep_IrriECw,&
          rep_Manag,&
          rep_sim,&
          rep_soil,&
          rep_sum

  use ac_run, only: &
          repIrriInfoRecord,&
          rep_StressTot

  use ac_kinds, only: dp,&
                int8,&
                int32,&
                intEnum,&
                sp
  implicit none
  private
  type, public :: ac71dec
    !-------------------------------------------------------------------------
    ! forcing
    !-------------------------------------------------------------------------
    real               :: tair
    real               :: tmax
    real               :: tmin
    real               :: tdew
    real               :: psurf
    real               :: wndspd
    real               :: swdown
    real               :: prcp
    real               :: eto
    !-------------------------------------------------------------------------
    ! spatial parameter
    !-------------------------------------------------------------------------
    character(len=100) :: cropt
    integer            :: soiltype
    !-------------------------------------------------------------------------
    ! multilevel spatial parameter
    !-------------------------------------------------------------------------
    real, pointer      :: smceq(:)
    !-------------------------------------------------------------------------
    ! state
    !-------------------------------------------------------------------------
    real, pointer      :: smc(:)
    !-------------------------------------------------------------------------
    ! AC specific
    !-------------------------------------------------------------------------
    integer        :: daynri
    integer(int32) :: DaySubmerged
    integer(int32) :: GDDTadj
    integer(int32) :: InitializeRun
    integer(int32) :: NrCompartments
    integer(int32) :: SumInterval
    integer(int32) :: Tadj
    integer        :: PreviousStressLevel
    integer        :: StressSFadjNEW
    integer       :: irun
    integer       :: NoMoreCrop
    integer(intEnum) :: TheProjectType
    logical  :: HarvestNow
    logical :: NoYear
    real    :: RootZoneWC_Actual
    real    :: RootZoneWC_FC
    real    :: RootZoneWC_Leaf
    real    :: RootZoneWC_SAT
    real    :: RootZoneWC_Sen
    real    :: RootZoneWC_Thresh
    real    :: RootZoneWC_WP
    real    :: RootZoneWC_ZtopAct
    real    :: RootZoneWC_ZtopFC
    real    :: RootZoneWC_ZtopThresh
    real    :: RootZoneWC_ZtopWP
    real :: alfaHI
    real :: alfaHIAdj
    real :: Bin
    real :: Bout
    real :: CCiActual
    real :: CCiActualWeedInfested
    real :: CCiprev
    real :: CCiTopEarlySen
    real :: CCoTotal
    real :: CCxCropWeedsNoSFstress
    real :: CCxTotal
    real :: CCxWitheredTpotNoS
    real :: CDCTotal
    real :: CGCref,GDDCGCref 
    real :: Coeffb0
    real :: Coeffb0Salt
    real :: Coeffb1
    real :: Coeffb1Salt
    real :: Coeffb2
    real :: Coeffb2Salt
    real :: DayFraction
    real :: Drain  
    real :: Eact
    real :: ECstorage
    real :: Epot 
    real :: GDDayFraction
    real :: GDDayi
    real :: GDDCDCTotal
    real :: HItimesAT
    real :: HItimesAT1
    real :: HItimesAT2
    real :: HItimesBEF
    real :: Irrigation 
    real :: Runoff  
    real :: ScorAT1
    real :: ScorAT2
    real :: StressLeaf
    real :: StressSenescence
    real :: SumGDDcuts
    real :: SumKci
    real :: SumKcTop
    real :: SumKcTopStress
    real :: SurfaceStorage
    real :: Tact 
    real :: TactWeedInfested
    real :: TimeSenescence
    real :: Tpot 
    real :: WeedRCi
    real :: WPi
    real :: YprevSum
    real :: Ziprev
    type(CompartmentIndividual), dimension(12) :: Compartment
    type(repIrriInfoRecord) :: IrriInfoRecord1
    type(repIrriInfoRecord) :: IrriInfoRecord2
    type(rep_Crop) :: Crop
    type(rep_IrriECw) :: IrriECw
    type(rep_DayEventInt), dimension(5) :: IrriAfterSeason
    type(rep_DayEventInt), dimension(5) :: IrriBeforeSeason
    type(rep_Manag) :: Management
    type(rep_sim) :: Simulation
    type(rep_soil) :: Soil
    type(rep_StressTot) :: StressTot
    type(rep_sum) :: SumWaBal
    type(SoilLayerIndividual), dimension(5) :: soillayer

  end type ac71dec
end module Ac71_module
