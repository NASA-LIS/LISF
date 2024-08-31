!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module geowrsi2_module
!BOP
!
! !MODULE: geowrsi2_module.F90
!
! !REVISION HISTORY:
!
! 31 Jul 2011: Brad Wind; Initial Definition
! 11 Jan 2013: KR Arsenault; Cleaned up comments/documentation
! 13 Feb 2013: KR Arsenault; Added new single-input parameter declarations
!                             (to replace the old "e1v" routines)
! 25 Oct 2013: KR Arsenault;  Added GeoWRSI2.0 model to LIS-7
!
! !DESCRIPTION:
!  WRSI LSM option (WRSI) variables
!
! !USES:        
  use geowrsi2_physics_module, only : gSOScalcNtStepsOfForcingUsed
!EOP

  implicit none

!-------------------------------------------------------------------------
! WRSI Variables
!-------------------------------------------------------------------------
  type geowrsi2dec

   logical   :: Photosensitive
   integer*2 :: isPermWilted
   integer*2 :: PET                       ! Total potential evaporation [kg m-2]
   integer*2 :: PPT                       ! Total rainfall precip [kg m-2]
   integer*2 :: WHC                       ! Water holding capacity; parameter [mm]
   integer*2 :: Mask                      ! FEWSNET-African region landmask [-]
  ! integer*4 :: LGP                       ! Length of growing period; parameter [dekad]
              ! Deprecated this var in favor of using only LGP_TIMESTEPS
   integer*4 :: LGP_TIMESTEPS
   integer*4 :: SOS                       ! Start-of-season [in dekads]
              ! LIS can refresh (mid-run) the SOS, and initial & final ts for current LIS_rc%yr
  ! integer*4 :: SOS_TIMESTEP             ! Deprecated this var in favor of using only SOS
   integer*4 :: SOSa                      ! Start-of-season Anomaly [in dekads]
   integer*4 :: SOSClimFromDataFile       ! Added for convenience (see note in geowrsi2_readInputSettings)
   integer*4 :: SOSClim                   ! SOS Climatology [+/- dekad]
   integer*2 :: WRSI                      ! Water requirement satisfaction index [ratio]
   real*4    :: KF                        ! Fraction of Season [-];  fracSeasonData(x,y)
   integer*2 :: KF2                       ! Percent of growing season [%];   percentSeasonData(x,y)
   real*4    :: SWAT                      ! Soil water content [mm]
   integer*2 :: SWI                       ! Soil Water Index [%]
   real*8    :: SumWR                     ! Sum of Water Requirement [mm], over season
   real*8    :: SumET                     ! Sum of Evapotranspiration [mm], over season
   integer*2 :: WRSIphoto                 ! WRSI * (flrng + SOSa - 100) / flrng 
   real*8    :: TotalSurplusWater         ! Total surplus water ~ different stages of crop growth [mm]
   real*8    :: MaxSurplusWater           ! Max surplus water experienced in 1 dekad [mm]
   real*8    :: TotalWaterDeficit         ! Total water deficit ~ different stages of crop growth [mm]
   real*8    :: MaxWaterDeficit           ! Max water deficit experienced in 1 dekad [mm]
   real*4    :: TotalAETInitial           ! Actual evapotranspiration ~ Initial stage [mm]
   real*4    :: TotalWRInitial            ! Water requirement         ~ Initial stage [mm]
   real*8    :: TotalSurplusWaterInitial  ! Surplus water      ~ Initial stage [mm]
   real*8    :: TotalWaterDeficitInitial  ! Water deficit      ~ Initial stage [mm]
   real*4    :: TotalAETVeg               ! Actual evapotransp ~ Vegetative stage [mm]
   real*4    :: TotalWRVeg                ! Water requirement  ~ Vegetative stage [mm]
   real*8    :: TotalSurplusWaterVeg      ! Surplus water      ~ Vegetative stage [mm]
   real*8    :: TotalWaterDeficitVeg      ! Water deficit      ~ Vegetative stage [mm]
   real*4    :: TotalAETFlower            ! Actual evapotransp ~ Flowering stage [mm]
   real*4    :: TotalWRFlower             ! Water requirement  ~ Flowering stage [mm]
   real*8    :: TotalSurplusWaterFlower   ! Surplus water      ~ Flowering stage [mm]
   real*8    :: TotalWaterDeficitFlower   ! Water deficit      ~ Flowering stage [mm]
   real*4    :: TotalAETRipe              ! Actual evapotransp ~ Ripening stage [mm]
   real*4    :: TotalWRRipe               ! Water requirement  ~ Ripening stage [mm]
   real*8    :: TotalSurplusWaterRipe     ! Surplus water      ~ Ripening stage [mm]
   real*8    :: TotalWaterDeficitRipe     ! Water deficit      ~ Ripening stage [mm]

   integer*4 :: PermWiltDate              ! Permanent wilting date [dekad]
   integer*4 :: Wilting1                  ! First wilting date [dekad]
   integer*4 :: Wilting2                  ! Second wilting date [dekad]

   integer*2 :: WRSIa                     ! WRSI anomaly [-]; uses WRSI climatology parameter file
   integer*2 :: WRSIClim                  ! WRSI climatology parameter file (read-in from config file)

   integer*2 :: WRSI_TimeStep             ! WRSI per one dekad (+SY)
   real*8    :: WR_TimeStep               ! Water requirement (WR) per one dekad 
   real*8    :: AET_TimeStep              ! AET per one dekad 
   real*8    :: SurplusWater_TimeStep     ! Surplus Water per one dekad 

  ! Short description about this group of variables and also those variables 
  ! prefixed with 'LastCurrent' in their varnames below.  The 'Initial' and 'Final' 
  ! variables are values that are not calculated (e.g. not sos [start-of-season] 
  ! for the 'Initial', nor SOS+LGP [length-of-growing-period] for the 'Final' terms) but
  ! rather represent stored values read from disk.  Therefore, to some extent these 
  ! likely represent at best regional or local estimates, perhaps climatological. 
  ! The original GeoWRSI model used these data in exactly the same manner as is done 
  ! here in LIS, as is the case with all of the variables translated over.  The other 
  ! variables in this grouping constitute quantities derived directly from these data. +BW

   integer*2 :: InitialYear           ! Initial year of run specified from input file
   integer*4 :: InitialTStep          ! Initial timestep of growing season (based on region)

!- Added to support multiple growing seasons (+BW)
   integer*2 :: physics_activation_yr      ! Used for multi-season capability
   integer*4 :: physics_activation_ts      ! Used for multi-season capability
   integer*2 :: FinalYearFromDataFile      ! BW distinguishs here FinalYear read from data file (?)
   integer   :: growing_season             ! Added for multi-season run capability in LIS

 ! Here are the *actual* end-of-season (EOS) variables as-computed 
 !  (by adding SOS+LGP) by SetFinalTStep function:
   integer*2 :: FinalYear
   integer*4 :: FinalTStep

   logical   :: PW_EnablePermanentWilting
   logical   :: PermanentWiltingIsFromSWI
   logical   :: PermanentWiltingIsFromWRSI
   integer*4 :: PermanentWiltingThresholdValue

   logical   :: SOS_CALCULATE
   integer*4 :: SOS_offset
   logical   :: RunRegardlessOfSufficientData

 ! Crop coefficients (read in from crop type files):
   real*4    :: Crop_c1
   real*4    :: Crop_c2
   real*4    :: Crop_F1
   real*4    :: Crop_F2
   real*4    :: Crop_F3
   real*4    :: Crop_K1
   real*4    :: Crop_K2
   real*4    :: Crop_K3
   real*4    :: Crop_kp
   real*4    :: Crop_rini
   real*4    :: Crop_r
   real*4    :: PreSeason_Kc

 ! These next 2 were introduced during integration of CalcSOS:
   integer*2, dimension(gSOScalcNtStepsOfForcingUsed) :: PPTcacheForSOScalc
   integer*2, dimension(gSOScalcNtStepsOfForcingUsed) :: PETcacheForSOScalc

!- Accumulated PPT and PET variables for each WRSI timestep:
   real*4    :: ACC_PPT
   real*4    :: ACC_PET

!- SOScalc Inputs:
   logical   :: SOScalcUseRainThreshold
   integer*4 :: SOScalcTStep1Threshold
   integer*4 :: SOScalcTStep2Plus3Threshold
   logical   :: SOScalcUseSWIThreshold
   integer*4 :: SOScalcSoilMoistureThreshold
   logical   :: SOScalcUseNearestToDesiredPracticalPlantingTStep
   integer*4 :: SOScalcDesiredPlantingTStep
   logical   :: SOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart
   integer*4 :: SOScalcRestartThresholdWithWRSI
   integer*4 :: SOScalcRestartThresholdWithSWI
   logical   :: SOScalcRestartCropWhenFailed
   integer*4 :: SOScalcCropCanStillRestartWithinPctLGPThreshold
   logical   :: SOScalcIgnoreClimatology
   integer*4 :: SOScalcMaxTStepsLate
   integer*4 :: SOScalcMaxTStepsEarly
   logical   :: SOScalcExcludeIncompleteAreasFromSOS
   integer*4 :: SOScalcAcceptablePercentOfTotalSeasonalGrowth

   real, pointer :: sos_write(:)
   real, pointer :: sosa_write(:)
!   real :: sos_write
!   real :: sosa_write

  end type geowrsi2dec

! __________________________________________________________

!- WRSI User Input Settings - Parameter/variable declarations:

!= Integer:
   integer*4 :: InitialYear               ! Array: InitialYear
   integer*4 :: SimulatedCurrentDekad     ! Array: SimulatedCurrentDekad
   integer*4 :: SOSoffset                 ! Array: SOS_offset
!   integer*4 :: LGPnumdekads              ! 0
!   integer*4 :: SOSdekad                  ! 1
!   integer*4 :: WHCnum_mils               ! 0

  ! SOScalcs:
   integer*4 :: Dekad1Threshold           ! Array: SOScalcTStep1Threshold
   integer*4 :: Dekad2plus3Threshold      ! Array: SOScalcTStep2Plus3Threshold
   integer*4 :: SoilMoistureThreshold     ! Array: SOScalcSoilMoistureThreshold
   integer*4 :: maxDeksEarly              ! Array: SOScalcMaxTStepsEarly
   integer*4 :: maxDeksLate               ! Array: SOScalcMaxTStepsLate
   integer*4 :: restartThresholdWithWRSI  ! Array: SOScalcRestartThresholdWithWRSI
   integer*4 :: restartThresholdWithSWI   ! Array: SOScalcRestartThresholdWithSWI
   integer*4 :: CropCanStillRestartWithinPctLGPthreshold ! Array: SOScalcCropCanStillRestartWithinPctLGPThreshold
   integer*4 :: AcceptablePercentOfTotalSeasonalGrowth   ! Array: SOScalcAcceptablePercentOfTotalSeasonalGrowth
   integer*4 :: DesiredPlantingDekad      ! Array: SOScalcDesiredPlantingTStep
  ! EOScalcs:
!   integer*4 :: drynessThreshold !10

  ! PermWiltCalcs:
   integer*4 :: PWThresholdValue          ! Array: PermanentWiltingThresholdValue

!= Logical:
   logical :: EnablePermanentWilting      ! Array: PW_EnablePermanentWilting
   logical :: UseSimulatedCurrentDekad    ! Array: UseSimulatedCurrentDekad
!   logical :: UseDefaultMask    !True
!   logical :: Go_Beyond_Dek36   !False

 ! SOScalcs:
   logical :: calcsUseRainThreshold          ! Array: SOScalcUseRainThreshold
   logical :: calcsUseSWIthreshold           ! Array: SOScalcUseSWIThreshold
   logical :: IgnoreClimatology              ! Array: SOScalcIgnoreClimatology
   logical :: RestartCropWhenFailed          ! Array: SOScalcRestartCropWhenFailed
   logical :: ExcludeIncompleteAreasFromSOS  ! Array: SOScalcExcludeIncompleteAreasFromSOS
   logical :: useWRSIthresholdInsteadofSWIthresholdToRestart ! Array: SOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart
   logical :: UseNearestToDesiredPracticalPlantingDekad ! Array: SOScalcUseNearestToDesiredPracticalPlantingTStep

 ! PermWiltCalcs:
   logical :: PermanentWiltingIsFromSWI      ! Array: PermanentWiltingIsFromSWI
   logical :: PermanentWiltingIsFromWRSI     ! Array: PermanentWiltingIsFromWRSI
   logical :: PermanentWiltingIsFromPlantPhysiology !False  ?? Array: Photosensitive ??

!= Character String:
   character(40) :: AnalysisRegion ! "East Africa (Oct-Feb)"
   character(40) :: CropType       ! "maize"

!= Crop file inputs:
   real*4    :: crop_F1
   real*4    :: crop_F2
   real*4    :: crop_F3
   real*4    :: crop_K1
   real*4    :: crop_K2
   real*4    :: crop_K3
   real*4    :: crop_kp
   real*4    :: crop_rini
   real*4    :: crop_r
   real*4    :: crop_c1
   real*4    :: crop_c2

end module geowrsi2_module
