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
! !ROUTINE: geowrsi2_main
! \label{geowrsi2_main}
!
! !REVISION HISTORY:
! 31 Jul 2011: Brad Wind; Initial Definition
! 19 Dec 2012: KR Arsenault; Removed test code
! 19 Feb 2013: KR Arsenault; Removed param array for multi season component
! 25 Oct 2013: KR Arsenault; Added GeoWRSI2.0 model to LIS-7
! 
! !INTERFACE:
subroutine geowrsi2_main(n)

! !USES:
  use LIS_coreMod,   only : LIS_rc, LIS_domain
  use LIS_logMod,    only : LIS_logunit, LIS_endrun
  use LIS_timeMgrMod
  use LIS_histDataMod
  use geowrsi2_module
  use geowrsi2_lsmMod, only : geowrsi2_struc, geowrsi2_CalcSOSlsmRunMode, &
                              num_growing_seasons, geowrsi2_udef, outVar

  use geowrsi2_physics_module, only :  &
      offsetTStepByNumTSteps, CalcWRSI, CalcWRSIAnomalies,      &
      CalcSOS, resetWBdataFor1point, gTimeStepsPerYear,         &
      gCurrentTStep, gCurrentYear, gInitialTStep, gInitialYear, &
      gPW_EnablePermanentWilting, gPermanentWiltingIsFromSWI,   &
      gPermanentWiltingIsFromWRSI, gPermanentWiltingThresholdValue, &
      gCrop_c1, gCrop_c2, gCrop_F1, gCrop_F2, gCrop_F3, gCrop_K1,   &
      gCrop_K2, gCrop_K3, gCrop_kp, gCrop_rini, gCrop_r, gPreSeason_Kc, &
      gLastCurrentTStep, gLastCurrentYear, &
      gSOScalcUseRainThreshold, gSOScalcTStep1Threshold,     &
      gSOScalcTStep2Plus3Threshold, gSOScalcUseSWIThreshold, &
      gSOScalcSoilMoistureThreshold, gSOScalcDesiredPlantingTStep,  &
      gSOScalcUseNearestToDesiredPracticalPlantingTStep,      &
      gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart, &
      gSOScalcRestartThresholdWithWRSI, gSOScalcRestartThresholdWithSWI,& 
      gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart, &
      gSOScalcRestartCropWhenFailed, &
      gSOScalcCropCanStillRestartWithinPctLGPThreshold,   &
      gFinalYear, gFinalTStep, gSOScalcIgnoreClimatology, &
      gSOScalcMaxTStepsLate, gSOScalcMaxTStepsEarly,      &
      gSOScalcExcludeIncompleteAreasFromSOS,              &
      gSOScalcAcceptablePercentOfTotalSeasonalGrowth,     &
      gMASK_INC, gWRSI_NA, gKF2_SEASONTOSTARTDANDORACCUM, &
      DifferenceOf2TSteps

  use fbil_module, only : CInt2, CInt4

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
!
! !DESCRIPTION:
! 
!  Calls the run routines for the forcing-only option (geowrsi2)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  integer             :: t, s, Tm
  integer             :: offset
  logical             :: isTimeToRunCheck
  integer*4           :: timestep
  logical             :: tmpPW_EnablePermanentWilting
  logical             :: do_cycle
  logical             :: attemptPostSOScalcSOSreset
  logical             :: resetWBdata
  type(geowrsi2dec), pointer :: geowrsi2Pt

  integer*2, allocatable, dimension(:) :: total_precip_in
  integer*2, allocatable, dimension(:) :: total_pet_in
  character*3         :: fnest
! ______________________________________________________________________


  if( LIS_rc%npatch(n,LIS_rc%lsm_index) == 0 ) then
!        write(LIS_logunit,*) " MSG: THIS PROCESSING ELEMENT,",LIS_localPet,&
!       " HAS ZERO-LAND TILES WITHIN DOMAIN FOR NEST, ",n,".  SKIPPING ..."
     return
  endif

  write(fnest,'(i3.3)') n    

  offset = 0  

  isTimeToRunCheck = LIS_isAlarmRinging(LIS_rc, "GEOWRSI2 model alarm "//trim(fnest), "dekad")

!  if(isTimeToRunCheck == .true. ) print *, "Geowrsi runtime check: ", isTimeToRunCheck
  if(isTimeToRunCheck .neqv. .true.) return

! Estimate dekad-based timestep in timing with the starting dekad ...
  geowrsi2_struc(n)%mon_dekad = LIS_getDekad( LIS_rc, offset=offset )

  if( offset > 0 .and. geowrsi2_struc(n)%mon_dekad == 3 ) then
    timestep = ((LIS_rc%mo-1)*3) 
  else
    timestep = ((LIS_rc%mo-1)*3) + geowrsi2_struc(n)%mon_dekad
  endif
  if( timestep < 1 ) timestep = gTimeStepsPerYear    ! 36, for dekads

  geowrsi2_struc(n)%year_dekad = timestep

  write(LIS_logunit,'(a50,3(1x,i4,a1))')">> LIS-WRSI dekad tstep, month, dekad of month:",&
                     geowrsi2_struc(n)%year_dekad,",",LIS_rc%mo,",",geowrsi2_struc(n)%mon_dekad
!  write(*,'(a50,4(1x,i4,a1))')">> LIS-WRSI year, dekad tstep, month, dekad of month:",&
!             LIS_rc%yr,",",geowrsi2_struc(n)%year_dekad,",",LIS_rc%mo,",",geowrsi2_struc(n)%mon_dekad


! Allocate total Precip and PET arrays:
  allocate(total_precip_in(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(total_pet_in(LIS_rc%npatch(n,LIS_rc%lsm_index)))

!!! BEGIN MAIN WRSI PHYSICS RUN LOOP !!!

  do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)

     geowrsi2_struc(n)%wrsi(t)%PPT = CInt2(real(geowrsi2_struc(n)%wrsi(t)%ACC_PPT, 8))
     geowrsi2_struc(n)%wrsi(t)%PET = CInt2(real(geowrsi2_struc(n)%wrsi(t)%ACC_PET, 8))
     total_precip_in(t) = geowrsi2_struc(n)%wrsi(t)%PPT   ! save PPT for output
     total_pet_in(t)    = geowrsi2_struc(n)%wrsi(t)%PET   ! save PET for output

     geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)

     geowrsi2Pt%WR_TimeStep   = LIS_rc%udef        ! +SY
     geowrsi2Pt%AET_TimeStep  = LIS_rc%udef               
     geowrsi2Pt%WRSI_TimeStep = LIS_rc%udef               
     geowrsi2Pt%SurplusWater_TimeStep = LIS_rc%udef       


  !- Skip WRSI computation for the present grid cell if we are before 
  !  the start of the physics activation for the grid cell.
     if( 0 < DifferenceOf2TSteps(                    &
               geowrsi2Pt%physics_activation_yr,     &
               geowrsi2Pt%physics_activation_ts,     &
               CInt2(real(LIS_rc%yr, 8)), timestep,  &
               gTimeStepsPerYear)                    &
       ) cycle 

  !- WRSI Physics Are Activated Here: ----
     if( (geowrsi2Pt%physics_activation_yr == LIS_rc%yr) .and. &
         (geowrsi2Pt%physics_activation_ts == timestep) ) then

     ! Reset all water balance quantities at start of season:
       if( geowrsi2Pt%Mask == gMASK_INC ) then  ! If mask = 1
          call resetWBdataFor1point(               &
              .false.,                             &
              geowrsi2Pt%isPermWilted,             &
              geowrsi2Pt%WRSI,                     &
              geowrsi2Pt%KF,                       &
              geowrsi2Pt%KF2,                      &
              geowrsi2Pt%SWAT,                     &
              geowrsi2Pt%SumWR,                    &
              geowrsi2Pt%SumET,                    &
              geowrsi2Pt%SWI,                      &
              geowrsi2Pt%MaxSurplusWater,          & 
              geowrsi2Pt%MaxWaterDeficit,          &
              geowrsi2Pt%TotalSurplusWaterInitial, &
              geowrsi2Pt%TotalWaterDeficitInitial, & 
              geowrsi2Pt%TotalSurplusWaterVeg,     &
              geowrsi2Pt%TotalWaterDeficitVeg,     &
              geowrsi2Pt%TotalSurplusWaterFlower,  &
              geowrsi2Pt%TotalWaterDeficitFlower,  &
              geowrsi2Pt%TotalSurplusWaterRipe,    &
              geowrsi2Pt%TotalWaterDeficitRipe,    &
              geowrsi2Pt%PermWiltDate,             &
              geowrsi2Pt%TotalSurplusWater,        &
              geowrsi2Pt%TotalWaterDeficit,        &
              geowrsi2Pt%TotalAETInitial,          &
              geowrsi2Pt%TotalWRInitial,           &
              geowrsi2Pt%TotalAETVeg,              &
              geowrsi2Pt%TotalWRVeg,               &
              geowrsi2Pt%TotalAETFlower,           &
              geowrsi2Pt%TotalWRFlower,            &
              geowrsi2Pt%TotalAETRipe,             &
              geowrsi2Pt%TotalWRRipe,              &
              geowrsi2Pt%WRSIphoto,                &
              geowrsi2Pt%SOSa )

          geowrsi2Pt%SWAT = 0

       else  ! Exclude where mask = 0 or undef
          geowrsi2Pt%SWAT = -99
          geowrsi2Pt%SWI  = 255
       endif
       geowrsi2Pt%WRSI = gWRSI_NA
       geowrsi2Pt%KF   = -1
       geowrsi2Pt%KF2  = gKF2_SEASONTOSTARTDANDORACCUM
       if(geowrsi2Pt%WHC < 1) geowrsi2Pt%SWI = 255  ! N/A

     ! Advance to next physics activation year (if multiple growing seasons):
       if( t == LIS_rc%npatch(n,LIS_rc%lsm_index) ) then
          geowrsi2_struc(n)%next_physics_act_year = geowrsi2Pt%physics_activation_yr + 1
       endif

     endif   ! End Physics Activation Check for Water Balance Reset

  !- SOS Run-mode:  Store PPT and PET
     if( geowrsi2_CalcSOSlsmRunMode .eqv. .true. ) then
       if( geowrsi2Pt%SOS_CALCULATE .neqv. .true. ) &
         cycle ! Do nothing;
               ! Existing SOS value read during parms setup 
               !  will be the value retained, unchanged

     ! Make sure SOS calculation is not made past the final timestep
     !  set for a given region (found in region file):
       if( (timestep-3) == geowrsi2Pt%FinalTStep ) then  ! +KRA (added 7.12.2013)
          geowrsi2_struc(n)%end_soscalc = .true.
       endif

     ! SOS calculation is performed in-sync but out-of-step with 
     !  real-time by using these cache buffers, the contents of 
     !  which must always shift with the timestep as it marches forward.
       if( geowrsi2_struc(n)%end_soscalc ) then   ! Last SOS calc timestep
         geowrsi2Pt%PPTcacheForSOScalc = 0.

         do Tm = 1, (gSOScalcNtStepsOfForcingUsed-1)
            geowrsi2Pt%PETcacheForSOScalc(Tm) = geowrsi2Pt%PETcacheForSOScalc(Tm+1)
         end do
         geowrsi2Pt%PETcacheForSOScalc(gSOScalcNtStepsOfForcingUsed) = geowrsi2Pt%PET
         geowrsi2Pt%PET = geowrsi2Pt%PETcacheForSOScalc(1)
         geowrsi2Pt%PPT = geowrsi2Pt%PPTcacheForSOScalc(1)

       else  ! Store Precipitation Inputs for SOS calculation:
         do Tm = 1, (gSOScalcNtStepsOfForcingUsed-1)
            geowrsi2Pt%PPTcacheForSOScalc(Tm) = geowrsi2Pt%PPTcacheForSOScalc(Tm+1)
            geowrsi2Pt%PETcacheForSOScalc(Tm) = geowrsi2Pt%PETcacheForSOScalc(Tm+1)
         end do
         geowrsi2Pt%PPTcacheForSOScalc(gSOScalcNtStepsOfForcingUsed) = geowrsi2Pt%PPT
         geowrsi2Pt%PETcacheForSOScalc(gSOScalcNtStepsOfForcingUsed) = geowrsi2Pt%PET
         geowrsi2Pt%PPT = geowrsi2Pt%PPTcacheForSOScalc(1)
         geowrsi2Pt%PET = geowrsi2Pt%PETcacheForSOScalc(1)
       end if

     endif   ! End SOS Run-mode condition


   ! First check whether any of the quantities required for CalcWRSI are undefined.
     if( (geowrsi2Pt%InitialYear  == CInt2(real(geowrsi2_udef, 8))) .or. &
         (geowrsi2Pt%InitialTStep == CInt4(real(geowrsi2_udef, 8))) .or. &
         (geowrsi2Pt%Crop_c1      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_c2      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_F1      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_F2      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_F3      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_K1      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_K2      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_K3      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_kp      == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_rini    == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%Crop_r       == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%PreSeason_Kc == geowrsi2_udef)                 .or. &
         (geowrsi2Pt%SOSClim      == CInt4(real(geowrsi2_udef, 8)))      &
       ) cycle

     if( (geowrsi2Pt%PPT == CInt2(real(geowrsi2_udef, 8))) .or. &
         (geowrsi2Pt%PET == CInt2(real(geowrsi2_udef, 8)))      &
       ) cycle

     if( geowrsi2Pt%SOS == CInt4(real(geowrsi2_udef, 8)) ) cycle


  !! Align the time with LIS-master time manager ...
     gCurrentYear  = LIS_rc%yr
     gCurrentTStep = timestep

     if( geowrsi2_CalcSOSlsmRunMode .eqv. .true. ) then   ! SOS Run-mode
         call offsetTStepByNumTSteps(gCurrentTStep, (0-2), gCurrentYear)
     endif

  !! Align all other variables in the model physics module code for CalcWRSI:
     gInitialTStep = geowrsi2Pt%InitialTStep
     gInitialYear  = geowrsi2Pt%InitialYear

     gPW_EnablePermanentWilting      = geowrsi2Pt%PW_EnablePermanentWilting
     gPermanentWiltingIsFromSWI      = geowrsi2Pt%PermanentWiltingIsFromSWI
     gPermanentWiltingIsFromWRSI     = geowrsi2Pt%PermanentWiltingIsFromWRSI
     gPermanentWiltingThresholdValue = geowrsi2Pt%PermanentWiltingThresholdValue

     gCrop_c1   = geowrsi2Pt%Crop_c1
     gCrop_c2   = geowrsi2Pt%Crop_c2
     gCrop_F1   = geowrsi2Pt%Crop_F1
     gCrop_F2   = geowrsi2Pt%Crop_F2
     gCrop_F3   = geowrsi2Pt%Crop_F3
     gCrop_K1   = geowrsi2Pt%Crop_K1
     gCrop_K2   = geowrsi2Pt%Crop_K2
     gCrop_K3   = geowrsi2Pt%Crop_K3
     gCrop_kp   = geowrsi2Pt%Crop_kp
     gCrop_rini = geowrsi2Pt%Crop_rini
     gCrop_r    = geowrsi2Pt%Crop_r
     gPreSeason_Kc = geowrsi2Pt%PreSeason_Kc

     tmpPW_EnablePermanentWilting = gPW_EnablePermanentWilting

     if( geowrsi2_CalcSOSlsmRunMode .eqv. .true. ) then
        gPW_EnablePermanentWilting = .false.
     endif
     if( (gPW_EnablePermanentWilting .eqv. .true.) .and. &
         (geowrsi2Pt%PermanentWiltingThresholdValue == CInt4(real(geowrsi2_udef, 8))) &
       ) then
        gPW_EnablePermanentWilting = tmpPW_EnablePermanentWilting
        cycle
     end if

   ! Note: the following physics call will have gCalcSOS mode set 
   !  (see geowrsi2_lsmMod and wrsi_lsm_ini) when running in that mode
!     print *, " About to call CalcWRSI from geowrsi2_main ... "

   ! Calculate WRSI and other model variables:
     call CalcWRSI(                            &
          geowrsi2Pt%photosensitive,           &
          geowrsi2Pt%isPermWilted,             &
          geowrsi2Pt%PET,                      &
          geowrsi2Pt%PPT,                      &
          geowrsi2Pt%LGP_TIMESTEPS,            &
          geowrsi2Pt%WHC,                      &
          geowrsi2Pt%Mask,                     &
          geowrsi2Pt%SOS,                      &
          geowrsi2Pt%WRSI,                     &
          geowrsi2Pt%SOSClim,                  &
          geowrsi2Pt%KF,                       &
          geowrsi2Pt%KF2,                      &
          geowrsi2Pt%SWAT,                     &
          geowrsi2Pt%SumWR,                    &
          geowrsi2Pt%SumET,                    &
          geowrsi2Pt%SWI,                      &
          geowrsi2Pt%WRSIphoto,                &
          geowrsi2Pt%SOSa,                     &
          geowrsi2Pt%TotalSurplusWater,        &
          geowrsi2Pt%MaxSurplusWater,          &
          geowrsi2Pt%TotalWaterDeficit,        &
          geowrsi2Pt%MaxWaterDeficit,          &
          geowrsi2Pt%TotalAETInitial,          &
          geowrsi2Pt%TotalWRInitial,           &
          geowrsi2Pt%TotalSurplusWaterInitial, &
          geowrsi2Pt%TotalWaterDeficitInitial, &
          geowrsi2Pt%TotalAETVeg,              &
          geowrsi2Pt%TotalWRVeg,               &
          geowrsi2Pt%TotalSurplusWaterVeg,     &
          geowrsi2Pt%TotalWaterDeficitVeg,     &
          geowrsi2Pt%TotalAETFlower,           &
          geowrsi2Pt%TotalWRFlower,            &
          geowrsi2Pt%TotalSurplusWaterFlower,  &
          geowrsi2Pt%TotalWaterDeficitFlower,  &
          geowrsi2Pt%TotalAETRipe,             &
          geowrsi2Pt%TotalWRRipe,              &
          geowrsi2Pt%TotalSurplusWaterRipe,    &
          geowrsi2Pt%TotalWaterDeficitRipe,    &
          geowrsi2Pt%PermWiltDate,             &
          geowrsi2Pt%WR_TimeStep,              &
          geowrsi2Pt%AET_TimeStep,             &
          geowrsi2Pt%WRSI_TimeStep,            &
          geowrsi2Pt%SurplusWater_TimeStep ) 


  !- WRSI-runmode: Calculate WRSI anomalies
     if( geowrsi2_CalcSOSlsmRunMode .neqv. .true. ) then

    ! The following check ensures that we compute anomalies only once the season 
    !  starts and never before (same as done in the original GeoWRSI model):
        if( 0 >= DifferenceOf2TSteps(geowrsi2Pt%InitialYear,    &
                                     geowrsi2Pt%InitialTStep,   &
                                     CInt2(real(LIS_rc%yr, 8)), & 
                                     timestep, gTimeStepsPerYear) ) then
           call CalcWRSIAnomalies(       &
               geowrsi2Pt%Mask,          &
               geowrsi2Pt%LGP_TIMESTEPS, &
               geowrsi2Pt%WRSIa,         &
               geowrsi2Pt%WRSI,          &
               geowrsi2Pt%SWI,           & 
               geowrsi2Pt%KF2,           &
               geowrsi2Pt%WRSIClim )
             cycle
        endif
     endif

  !! CalcSOS LIS Runmode - Specific code begin !!
     gPW_EnablePermanentWilting = tmpPW_EnablePermanentWilting

     do_cycle = .false.
     do Tm = 1, gSOScalcNtStepsOfForcingUsed
        if( geowrsi2Pt%PPTcacheForSOScalc(Tm) == CInt2(real(geowrsi2_udef, 8)) ) then
           do_cycle = .true.
           exit
        endif
     end do

     if( do_cycle .eqv. .true.                                                      .or. &
        (geowrsi2Pt%SWI == CInt2(real(geowrsi2_udef, 8)))                           .or. &
        (geowrsi2Pt%SOScalcTStep1Threshold      == CInt4(real(geowrsi2_udef, 8)))   .or. &
        (geowrsi2Pt%SOScalcTStep2Plus3Threshold == CInt4(real(geowrsi2_udef, 8)))   .or. &
        (geowrsi2Pt%SOScalcSoilMoistureThreshold== CInt4(real(geowrsi2_udef, 8)))   .or. &
        (geowrsi2Pt%SOScalcDesiredPlantingTStep == CInt4(real(geowrsi2_udef, 8)))   .or. &
        (geowrsi2Pt%SOScalcRestartThresholdWithSWI == CInt4(real(geowrsi2_udef, 8))).or. &
        (geowrsi2Pt%SOScalcCropCanStillRestartWithinPctLGPThreshold ==               &
                    CInt4(real(geowrsi2_udef, 8)))                                  .or. &
        (geowrsi2Pt%FinalYearFromDataFile == CInt2(real(geowrsi2_udef, 8)))         .or. &
        (geowrsi2Pt%FinalTStep            == CInt4(real(geowrsi2_udef, 8)))         .or. &
        (geowrsi2Pt%SOScalcMaxTStepsLate  == CInt4(real(geowrsi2_udef, 8)))         .or. &
        (geowrsi2Pt%SOScalcMaxTStepsEarly == CInt4(real(geowrsi2_udef, 8)))         .or. &
        (geowrsi2Pt%SOScalcAcceptablePercentOfTotalSeasonalGrowth ==                 &
           CInt4(real(geowrsi2_udef, 8)))                                            &
      ) cycle

   ! Align all other variables in the model physics module for CalcSOS:
   ! LastCurrentYear is the end year of the run:
     gLastCurrentYear = geowrsi2_struc(n)%LastCurrentYear

   ! LastCurrentTStep is the timestep-of-year (e.g., dekad-of-year in the range of [1,36]) 
   !      of the end of the run, and where typically this date represents the end of the 
   !      very last growing season of the run.
     if( geowrsi2_struc(n)%lastSOScalcOfSeasonTStep == 1 ) then 
        gLastCurrentTStep = gTimeStepsPerYear  ! 36, for dekad
     elseif( geowrsi2_struc(n)%lastSOScalcOfSeasonTStep > 1  .and. &
             geowrsi2_struc(n)%lastSOScalcOfSeasonTStep <= gTimeStepsPerYear ) then  ! <=36 
        gLastCurrentTStep = geowrsi2_struc(n)%lastSOScalcOfSeasonTStep - 1
     endif

   ! Set remaining variables/terms to grid points:
     gFinalYear                    = geowrsi2Pt%FinalYearFromDataFile
     gFinalTStep                   = geowrsi2Pt%FinalTStep
     gSOScalcMaxTStepsLate         = geowrsi2Pt%SOScalcMaxTStepsLate
     gSOScalcMaxTStepsEarly        = geowrsi2Pt%SOScalcMaxTStepsEarly
     gSOScalcUseRainThreshold      = geowrsi2Pt%SOScalcUseRainThreshold
     gSOScalcTStep1Threshold       = geowrsi2Pt%SOScalcTStep1Threshold
     gSOScalcTStep2Plus3Threshold  = geowrsi2Pt%SOScalcTStep2Plus3Threshold
     gSOScalcUseSWIThreshold       = geowrsi2Pt%SOScalcUseSWIThreshold
     gSOScalcSoilMoistureThreshold = geowrsi2Pt%SOScalcSoilMoistureThreshold
     gSOScalcUseNearestToDesiredPracticalPlantingTStep = &
                geowrsi2Pt%SOScalcUseNearestToDesiredPracticalPlantingTStep
     gSOScalcDesiredPlantingTStep  = geowrsi2Pt%SOScalcDesiredPlantingTStep
     gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart = &
                geowrsi2Pt%SOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart
     gSOScalcRestartThresholdWithWRSI = geowrsi2Pt%SOScalcRestartThresholdWithWRSI
     gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart = &
                geowrsi2Pt%SOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart
     gSOScalcRestartThresholdWithSWI  = geowrsi2Pt%SOScalcRestartThresholdWithSWI
     gSOScalcRestartCropWhenFailed    = geowrsi2Pt%SOScalcRestartCropWhenFailed
     gSOScalcCropCanStillRestartWithinPctLGPThreshold = &
                geowrsi2Pt%SOScalcCropCanStillRestartWithinPctLGPThreshold
     gSOScalcIgnoreClimatology     = geowrsi2Pt%SOScalcIgnoreClimatology
     gSOScalcExcludeIncompleteAreasFromSOS = geowrsi2Pt%SOScalcExcludeIncompleteAreasFromSOS
     gSOScalcAcceptablePercentOfTotalSeasonalGrowth = &
                geowrsi2Pt%SOScalcAcceptablePercentOfTotalSeasonalGrowth

   ! Post-SOS calculation and SOS reset parameter:
     attemptPostSOScalcSOSreset = .false.
     if( (LIS_rc%yr == geowrsi2_struc(n)%lastSOScalcOfSeasonYr ) .and. &
         (timestep  == geowrsi2_struc(n)%lastSOScalcOfSeasonTStep) ) then
       attemptPostSOScalcSOSreset = .true. 
     endif

   ! Attempt post-season SOS reset (to occur under certain specific 
   !  conditions given in that procedure; geowrsi2Pt%SOSClim=60)

     call CalcSOS(                         &
          geowrsi2Pt%Mask,                 &
          geowrsi2Pt%SOSClim,              &
          geowrsi2Pt%PPTcacheForSOScalc,   &
          geowrsi2Pt%WRSI,                 &
          geowrsi2Pt%KF2,                  &
          geowrsi2Pt%SWI,                  &
          geowrsi2Pt%Wilting1,             &
          geowrsi2Pt%Wilting2,             &
          resetWBdata,                     &
          geowrsi2Pt%SOS,                  &
          geowrsi2Pt%SOSa,                 &
          attemptPostSOScalcSOSreset )

     if( resetWBdata .eqv. .true. ) then 
        call resetWBdataFor1point(              &
           .true.,                              &
           geowrsi2Pt%isPermWilted,             &
           geowrsi2Pt%WRSI,                     &
           geowrsi2Pt%KF,                       &
           geowrsi2Pt%KF2,                      &
           geowrsi2Pt%SWAT,                     &
           geowrsi2Pt%SumWR,                    &
           geowrsi2Pt%SumET,                    &
           geowrsi2Pt%SWI,                      &
           geowrsi2Pt%MaxSurplusWater,          & 
           geowrsi2Pt%MaxWaterDeficit,          &
           geowrsi2Pt%TotalSurplusWaterInitial, &
           geowrsi2Pt%TotalWaterDeficitInitial, & 
           geowrsi2Pt%TotalSurplusWaterVeg,     &
           geowrsi2Pt%TotalWaterDeficitVeg,     &
           geowrsi2Pt%TotalSurplusWaterFlower,  &
           geowrsi2Pt%TotalWaterDeficitFlower,  &
           geowrsi2Pt%TotalSurplusWaterRipe,    &
           geowrsi2Pt%TotalWaterDeficitRipe,    &
           geowrsi2Pt%PermWiltDate              &
        )
     endif

  end do  ! end tile loop

!!! END MAIN WRSI PHYSICS RUN LOOP !!!


!!! WRSI MAIN MODEL VARIABLE OUTPUT DIAGNOSIS !!!

   do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
      geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)

   ! LIS standard output:
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF, &
          value=real(total_precip_in(t)),vlevel=1,unit="kg m-2",direction="DN", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP, &
          value=real(total_pet_in(t)),vlevel=1,unit="kg m-2",direction="UP", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOS, &
          value=real(geowrsi2Pt%SOS),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WRSI, &
          value=real(geowrsi2Pt%WRSI),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_KF2, &
          value=real(geowrsi2Pt%KF2),vlevel=1,unit="%",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWI, &
          value=real(geowrsi2Pt%SWI),vlevel=1,unit="%",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOSa, &
          value=real(geowrsi2Pt%SOSa),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_PermWiltDate, &
          value=real(geowrsi2Pt%PermWiltDate),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_Wilting1, &
          value=real(geowrsi2Pt%Wilting1),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_Wilting2, &
          value=real(geowrsi2Pt%Wilting2),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WRSIa, &
          value=real(geowrsi2Pt%WRSIa),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_growing_season, &
          value=real(geowrsi2Pt%growing_season),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
 
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SumWR, &
          value=real(outVar(geowrsi2Pt%SumWR,geowrsi2Pt%WRSI)), vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SumET, &
          value=real(outVar(geowrsi2Pt%SumET,geowrsi2Pt%WRSI)), vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalSurplusWater, &
          value=real(outVar(geowrsi2Pt%TotalSurplusWater,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_MaxSurplusWater, &
          value=real(outVar(geowrsi2Pt%MaxSurplusWater,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWaterDeficit, &
          value=real(outVar(geowrsi2Pt%TotalWaterDeficit,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_MaxWaterDeficit, &
          value=real(outVar(geowrsi2Pt%MaxWaterDeficit,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalAETInitial, &
          value=real(outVar(real(geowrsi2Pt%TotalAETInitial,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWRInitial, &
          value=real(outVar(real(geowrsi2Pt%TotalWRInitial,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalSurplusWaterInitial, &
          value=real(outVar(geowrsi2Pt%TotalSurplusWaterInitial,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWaterDeficitInitial, &
          value=real(outVar(geowrsi2Pt%TotalWaterDeficitInitial,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalAETVeg, &
          value=real(outVar(real(geowrsi2Pt%TotalAETVeg,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWRVeg, &
          value=real(outVar(real(geowrsi2Pt%TotalWRVeg,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalSurplusWaterVeg, &
          value=real(outVar(geowrsi2Pt%TotalSurplusWaterVeg,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWaterDeficitVeg, &
          value=real(outVar(geowrsi2Pt%TotalWaterDeficitVeg,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalAETFlower, &
          value=real(outVar(real(geowrsi2Pt%TotalAETFlower,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWRFlower, &
          value=real(outVar(real(geowrsi2Pt%TotalWRFlower,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalSurplusWaterFlower, &
          value=real(outVar(geowrsi2Pt%TotalSurplusWaterFlower,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWaterDeficitFlower, &
          value=real(outVar(geowrsi2Pt%TotalWaterDeficitFlower,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalAETRipe, &
          value=real(outVar(real(geowrsi2Pt%TotalAETRipe,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWRRipe, &
          value=real(outVar(real(geowrsi2Pt%TotalWRRipe,8),geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalSurplusWaterRipe, &
          value=real(outVar(geowrsi2Pt%TotalSurplusWaterRipe,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TotalWaterDeficitRipe, &
          value=real(outVar(geowrsi2Pt%TotalWaterDeficitRipe,geowrsi2Pt%WRSI)),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WHC, &
          value=real(geowrsi2Pt%WHC),vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LGP, &
          value=real(geowrsi2Pt%LGP_TIMESTEPS),vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)

     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WR_TimeStep, &
          value=real(geowrsi2Pt%WR_TimeStep), vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_AET_TimeStep, &
          value=real(geowrsi2Pt%AET_TimeStep), vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WRSI_TimeStep, &
          value=real(geowrsi2Pt%WRSI_TimeStep), vlevel=1,unit="-",direction="-", &
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SurplusWater_TimeStep, &
          value=real(geowrsi2Pt%SurplusWater_TimeStep), vlevel=1,unit="kg m-2",direction="-", &
          surface_type=LIS_rc%lsm_index)

   end do  ! end tile loop


! -------------------------------------

!- Assign SOS and SOS-anomaly values for each growing season and tile:
  if( geowrsi2_CalcSOSlsmRunMode ) then
    do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)
       do s = 1, num_growing_seasons
!       s = geowrsi2_struc(n)%season_count
!       if( geowrsi2_struc(n)%wrsi(t)%growing_season == geowrsi2_struc(n)%season_count .and. &
          if( geowrsi2_struc(n)%wrsi(t)%growing_season == s .and. &
            (geowrsi2_struc(n)%wrsi(t)%KF2 == 100. .or. &   ! end-of-growing season
             geowrsi2_struc(n)%wrsi(t)%KF2 == 200. .or. &   ! season is over
             geowrsi2_struc(n)%wrsi(t)%KF2 == 210.)) then   ! season has not started

            geowrsi2_struc(n)%wrsi(t)%sos_write(s) = geowrsi2_struc(n)%wrsi(t)%SOS
            geowrsi2_struc(n)%wrsi(t)%sosa_write(s)= geowrsi2_struc(n)%wrsi(t)%SOSa

          end if
       end do
    end do
  endif

! Advance seasonal parameters as applicable for the relevant time period; 
!  a.k.a. the Start of Season(SOS) and associated (time-checking, 
!  for final year and time step) parameters (cont'd):

  if( (LIS_rc%yr == geowrsi2_struc(n)%lastSOScalcOfSeasonYr ) .and. &
      (timestep  == geowrsi2_struc(n)%lastSOScalcOfSeasonTStep) ) then

  ! Trigger end-of-season (EOS) alarm:
    geowrsi2_struc(n)%eos_alarm = .true.

  ! Advance Last SOS Calculation Year:
    geowrsi2_struc(n)%lastSOScalcOfSeasonYr = geowrsi2_struc(n)%lastSOScalcOfSeasonYr + 1

  ! Advance growing season count for grid point:
    geowrsi2_struc(n)%season_count = geowrsi2_struc(n)%season_count + 1

    do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)

      ! Advance Initial and Final Year:
        geowrsi2Pt%InitialYear = geowrsi2Pt%InitialYear + 1
        geowrsi2Pt%FinalYear   = geowrsi2Pt%FinalYear + 1
        geowrsi2Pt%FinalYearFromDataFile = geowrsi2Pt%FinalYearFromDataFile + 1

      ! Advance to next physics activation year (if multiple growing seasons):
        geowrsi2Pt%physics_activation_yr = geowrsi2_struc(n)%next_physics_act_year

      ! Advance growing season count for grid point:
        geowrsi2Pt%growing_season = geowrsi2Pt%growing_season + 1   

      ! Reset Flag to determining when to end calculating SOS values: 
        geowrsi2_struc(n)%end_soscalc = .false.

     enddo
  endif

! -------------------------------------

!- Physics have been run, reset both the total precipitation accumulator
!  and the total evapotranspiration accumulator:

   geowrsi2_struc(n)%wrsi(:)%ACC_PPT = 0.0
   geowrsi2_struc(n)%wrsi(:)%ACC_PET = 0.0


end subroutine geowrsi2_main

