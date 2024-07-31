!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module geowrsi2_physics_module
!BOP
!
! !MODULE: geowrsi2_physics_module
! \label{geowrsi2_physics_module}
!
! !REVISION HISTORY:
! 31 Jul 2011: Brad Wind; Initial Definition
! 20 Dec 2012: KR Arsenault; Cleaned-up code
! 25 Oct 2013: KR Arsenault; Added GeoWRSI2.0 model to LIS-7
!
! !USES:
  use geowrsi2_arraymgmt_module
  use fbil_module
!EOP

implicit none
!
! Variables preceded with 'g'   :  global variables
! Comments annotated with 'cgrs':  any variables to be implemented as 
!        Crop Growing Region-Specific (CGRS or cgrs) in LIS-WRSI 
!  -- (+BW; April, 2011)
!
! Higher numerical prec for the following CALCSOS variables (integer assigned short):
  integer*4 :: gCurrentTStep   
  integer*2 :: gCurrentYear    
  integer*4 :: gInitialTStep   
  integer*2 :: gInitialYear  
! Note: Based on T. Magadire's April 27, 2011 e-mail and GeoWRSI in-program Help,
!  (see 'Configuring The GeoWRSI Water Balance'), gInitialYear represents the 
!  year containing the time step when the crop's growing season actually starts.  
!  For multi-crop runs this will be the year containing the time step of the start 
!  of the first growing season, for whichever crop growing region that happens to be.

! The following introduced together with introduction of CalcSOS ("cgrs"):
  integer*4 :: gFinalTStep           ! Short   
  integer*2 :: gFinalYear            ! Short   
  integer*2 :: gLastCurrentYear      ! Integer 
  integer*4 :: gLastCurrentTStep     ! Integer 

  logical   :: gPhotosensitive                 ! Boolean 
  logical   :: gPW_EnablePermanentWilting      ! Boolean 
  logical   :: gPermanentWiltingIsFromSWI      ! Boolean 
  logical   :: gPermanentWiltingIsFromWRSI     ! Boolean 
  integer*4 :: gPermanentWiltingThresholdValue ! Integer ThresholdValue 

! All crop-related gVars (coeff's, etc) are Single in GeoWRSI (VB), "cgrs":
  real*4 :: gCrop_c1      
  real*4 :: gCrop_c2      
  real*4 :: gCrop_F1      
  real*4 :: gCrop_F2      
  real*4 :: gCrop_F3     
  real*4 :: gCrop_K1      
  real*4 :: gCrop_K2      
  real*4 :: gCrop_K3      
  real*4 :: gCrop_kp      
  real*4 :: gCrop_rini    
  real*4 :: gCrop_r       
  real*4 :: gPreSeason_Kc 

! New gVars introduced for flexibility in design for WRSI and SOScalc:
  integer*4, parameter :: gTimeStepsPerYear         = 36   ! Number of dekads per year
  integer*4, parameter :: gTStepsBeforeSeasonStarts = 6    ! Number of TSteps before season starts
                                                           !  This is considered the "sowing" period
  integer*4, parameter :: gTStepsAfterSeasonEnds    = 10   ! Number of Tsteps after season ends

! The following for CalcSOS to be able to call CalcWRSI with minimal redundancy
  logical :: gCalcSOSmode = .false.

! The following for CalcSOS gSOScalcTStep2Plus3Threshold and associated cache array(s) sizing:
  integer*4, parameter :: gSOScalcNtStepsOfForcingUsed = 3  ! hard-wired to 3 for dekad time schema

! Global Variable ENCODED_VALUES Section
! Values are assigned here for those gVars;
! Otherwise, their values are set at run-time (e.g, from configuration information)
  integer*2, parameter :: gISPERMWILTED_NOTWILTED = 0
  integer*2, parameter :: gISPERMWILTED_WILTED    = 1
  integer*2, parameter :: gPPT_MISSINGDATABELOWTHISQTY = 0
  integer*2, parameter :: gPET_MISSINGDATABELOWTHISQTY = 0
  integer*4, parameter :: gLGP_MISSINGDATA = 0
  integer*2, parameter :: gWHC_MISSINGDATA = 0
  integer*4, parameter :: gLGP_PROCWRSIABOVETHISQTY = 0
  integer*2, parameter :: gWHC_PROCWRSIABOVETHISQTY = 0
  integer*2, parameter :: gPET_PROCWRSIABOVETHISQTY = 0
  integer*2, parameter :: gPPT_PROCWRSIGTETHISQTY   = 0
  integer*2, parameter :: gKF2_SEASONTOSTARTDANDORACCUM = 210
  integer*2, parameter :: gKF2_SEASONISOVER           = 200
  integer*2, parameter :: gKF2_SEASONHASPOSSIBLESTART = 0
  real*8, parameter    :: gSUMWR_EXITWRSIBELOWTHISQTY = 0
  real*8, parameter    :: gSUMWR_EXITWRSIABOVETHISQTY = 10000
  real*8, parameter    :: gSUMET_EXITWRSIABOVETHISQTY = 10000
  real*4, parameter    :: gKF_INIBELOWTHISPCT  = 0.16
  real*4, parameter    :: gKF_VEGGTETHISPCT    = 0.16
  real*4, parameter    :: gKF_VEGBELOWTHISPCT  = 0.43
  real*4, parameter    :: gKF_FLWRGTETHISPCT   = 0.43
  real*4, parameter    :: gKF_FLWRBELOWTHISPCT = 0.75
  real*4, parameter    :: gKF_RIPEGTETHISPCT   = 0.75
  integer*2, parameter :: gWRSI_NOSTART = 253 
                        ! No start of season *(normally occurs after normal start of season occurred) 
  integer*2, parameter :: gWRSI_SEASONYETTOSTART = 254    
                        ! *Season is yet to start
  integer*2, parameter :: gWRSI_NA = 255 
                        ! *N/A - masked out areas, where the WRSI is not relevant
  integer*2, parameter :: gWRSI_MAXANOMALYPCT    = 200    ! Hard-wired to 200 for dekad tstep
  integer*2, parameter :: gWRSI_MAXNUMDENOMANOMALY = 252  ! Hard-wired to 252 for dekad tstep
  integer*4, parameter :: gSOSPREOK_DIFLSTCURCURGTETHIS = 0
!*Note: These starred explanations provided by T. Magadzire, Dec 2010.

! ====== Start potentially timestep size-dependent ENCODED_VALUES ======
  integer*4, parameter :: gSOSPREOK_DIFLSTCURCURBELOWTHIS = 2    ! Hard-wired to 2 for dekad tstep +BW
  integer*4, parameter :: gSOS_NOSTART = 60                      ! Hard-wired to 60 for dekad tstep +BW
  integer*4, parameter :: gSOS_POSSIBLESTART1 = 50               ! Hard-wired to 50 for dekad tstep +BW
  integer*4, parameter :: gSOSCLIM_NOSTART = 60                  ! Hard-wired to 60 for dekad tstep +BW
  integer*4, parameter :: gSOSCALCUSENRSTTODSRDPRACPLNTTSTEPTHRESHOLD = 10  ! Hard-wired to 10 for dekad tstep +BW
  integer*4, parameter :: gPPTDATALIMITEDORMISSINGTHRESHOLD = 10 ! Hard-wired to 10 for dekad tstep +BW
! ====== End potentially timestep size-dependent ENCODED_VALUES ======

  integer*4, parameter :: gSOS_NA      = 0
  integer*4, parameter :: gSOSCLIM_NA  = 0
  integer*2, parameter :: gMASK_EXC    = 0
  integer*2, parameter :: gMASK_INC    = 1    ! gMaskNulVal
  integer*4, parameter :: gSOS_INITIAL = 0
  integer*4, parameter :: gWILTING1_DEFAULT = 0
  integer*4, parameter :: gWILTING2_DEFAULT = 0

! The following introduced together with introduction of CalcSOS ("cgrs"):
  logical   :: gSOScalcUseRainThreshold               ! Boolean 
  integer*4 :: gSOScalcTStep1Threshold                ! Integer 
  integer*4 :: gSOScalcTStep2Plus3Threshold           ! Integer 
  logical   :: gSOScalcUseSWIThreshold                ! Boolean  
  integer*4 :: gSOScalcSoilMoistureThreshold          ! Integer 
  logical   :: gSOScalcUseNearestToDesiredPracticalPlantingTStep  ! Boolean  
  integer*4 :: gSOScalcDesiredPlantingTStep           ! Integer
  logical   :: gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart ! Boolean 
  integer*4 :: gSOScalcRestartThresholdWithWRSI       ! Integer 
  integer*4 :: gSOScalcRestartThresholdWithSWI        ! Integer 
  logical   :: gSOScalcRestartCropWhenFailed          ! Boolean 
  integer*4 :: gSOScalcCropCanStillRestartWithinPctLGPThreshold   ! Integer 
  logical   :: gSOScalcIgnoreClimatology              ! Boolean 
  integer*4 :: gSOScalcMaxTStepsLate                  ! Integer 
  integer*4 :: gSOScalcMaxTStepsEarly                 ! Integer 
  logical   :: gSOScalcExcludeIncompleteAreasFromSOS  ! Boolean 
  integer*4 :: gSOScalcAcceptablePercentOfTotalSeasonalGrowth     ! Integer 

!_____________________________________________________________

 CONTAINS
!________


! "point" or "pixel-level" code of the WRSI computation

 subroutine CalcWRSI( photosensitive, isPermWilted, PET, PPT, LGP, WHC,              &
                     Mask, SOS, WRSI, SOSClim, KF, KF2, SWAT, SumWR, SumET,          &
                     SWI, WRSIphoto, SOSa, TotalSurplusWater, MaxSurplusWater,       &
                     TotalWaterDeficit, MaxWaterDeficit, TotalAETInitial,            &
                     TotalWRInitial, TotalSurplusWaterInitial,                       &
                     TotalWaterDeficitInitial, TotalAETVeg, TotalWRVeg,              &
                     TotalSurplusWaterVeg, TotalWaterDeficitVeg, TotalAETFlower,     &
                     TotalWRFlower, TotalSurplusWaterFlower, TotalWaterDeficitFlower,&
                     TotalAETRipe, TotalWRRipe, TotalSurplusWaterRipe,               &
                     TotalWaterDeficitRipe, PermWiltDate, WR_TimeStep, AET_TimeStep, &
                     WRSI_TimeStep, SurplusWater_TimeStep )

!
!  Routine:  CalcWRSI
!
!  Description:
!

 ! GeoWRSI (VB): Data types are retained and when variable used; (BW 12/15/2010)
   logical, intent(in)      :: photosensitive    ! Boolean photosensitive ! Determined atop main control loop, here checked lower-down
   integer*2, intent(inout) :: isPermWilted      ! Byte    isPermWilted(x,y) ! 0=not wilted, 1=wilted
   integer*2, intent(in)    :: PET               ! Short   petCube(x,y,deki) 
   integer*2, intent(in)    :: PPT               ! Short   pptCube(x,y,deki)
 ! Increased numerical prc of following variable to allow flexibility of design by bringing everything back to seconds
   integer*4, intent(inout) :: LGP               ! Byte    lgpData(x,y)
   integer*2, intent(inout) :: WHC               ! Short   WHC = whcData(x, y) ! WHC = whcData(x,y) * crop_r - WHICH IS RIGHT?
   integer*2, intent(inout) :: Mask              ! Byte    maskData(x,y)
 ! Increased numerical prec of following variables for calcsos which assigned an int to shrt
   integer*4, intent(in)    :: SOS               ! Short   sosData(x,y)
   integer*2, intent(inout) :: WRSI              ! Byte    geowrsi2Data(x,y)
   integer*2, intent(inout) :: WRSI_TimeStep     ! added by SY
   integer*4, intent(in)    :: SOSClim           ! Byte    sosClim(x,y)
   real*4, intent(inout)    :: KF                ! Single  fracSeasonData(x,y)
   integer*2, intent(inout) :: KF2               ! Byte    percentSeasonData(x,y)
   real*4, intent(inout)    :: SWAT              ! Single  swatData(x,y)
   real*8, intent(inout)    :: SumWR             ! Double  swrData(x,y)
   real*8, intent(inout)    :: WR_TimeStep       ! added by SY
   real*8, intent(inout)    :: SumET             ! Double  setData(x,y)
   real*8, intent(inout)    :: AET_TimeStep      ! added by SY
   integer*2, intent(inout) :: SWI               ! Byte    swiData(x,y)
   integer*2, intent(inout) :: WRSIphoto         ! Byte    geowrsi2data_photo(x,y)
   integer*4, intent(in)    :: SOSa              ! Byte    sosaData(x,y)
   real*8, intent(inout)    :: TotalSurplusWater ! Single  totalSurplusWater(x,y)
   real*8, intent(inout)    :: SurplusWater_TimeStep ! added by SY
   real*8, intent(inout)    :: MaxSurplusWater   ! Single  maxSurplusWater(x,y)
   real*8, intent(inout)    :: TotalWaterDeficit ! Single  totalWaterDeficit(x,y)
   real*8, intent(inout)    :: MaxWaterDeficit   ! Single  maxWaterDeficit(x,y)
   real*4, intent(inout)    :: TotalAETInitial, TotalWRInitial
                                                 ! Single totalAETinitial(x,y)
                                                 ! Single totalWRinitial(x,y)
   real*8, intent(inout)    :: TotalSurplusWaterInitial, TotalWaterDeficitInitial
                                                 ! Single totalSurplusWaterInitial(x,y)
                                                 ! Single totalWaterDeficitInitial(x,y)
   real*4, intent(inout)    :: TotalAETVeg, TotalWRVeg
                                                 ! Single totalAETveg(x,y)
                                                 ! Single totalWRveg(x,y)
   real*8, intent(inout)    :: TotalSurplusWaterVeg, TotalWaterDeficitVeg
                                                 ! Single totalSurplusWaterVeg(x,y)
                                                 ! Single totalWaterDeficitVeg(x,y)
   real*4, intent(inout)    :: TotalAETFlower, TotalWRFlower
                                                 ! Single totalAETflower(x,y)
                                                 ! Single totalWRflower(x,y)
   real*8, intent(inout)    :: TotalSurplusWaterFlower, TotalWaterDeficitFlower
                                                 ! Single totalSurplusWaterFlower(x,y)
                                                 ! Single totalWaterDeficitFlower(x,y)
   real*4, intent(inout)    :: TotalAETRipe, TotalWRRipe
                                                 ! Single totalAETripe(x,y)
                                                 ! Single totalWRripe(x,y)
   real*8, intent(inout)    :: TotalSurplusWaterRipe, TotalWaterDeficitRipe
                                                 ! Single totalSurplusWaterRipe(x,y)
                                                 ! Single totalWaterDeficitRipe(x,y)
   integer*4, intent(inout) :: PermWiltDate      ! Byte   PermWiltDate(x,y)

   logical   :: proceed
   real*4    :: Kc
   real*4    :: MAD
   real*4    :: SW
   real*4    :: ET
   real*8    :: SurplusWater
   real*8    :: flrng
   real*8    :: ETminWR
   integer*2 :: WRSIphoto_TimeStep   ! +SY
! ____________________________________________________________________________

   if( (gPW_EnablePermanentWilting .eqv. .true.) .and. &
       (isPermWilted /= gISPERMWILTED_NOTWILTED) ) return

 ! Eliminate cell from analysis if it has missing data:
   if( (PPT < gPPT_MISSINGDATABELOWTHISQTY) .or. &
       (PET < gPET_MISSINGDATABELOWTHISQTY) ) then
      LGP  = gLGP_MISSINGDATA
      WHC  = gWHC_MISSINGDATA
      Mask = gMASK_EXC
   endif

!- Select SOS value type case:
   select case (SOS)

      case (gSOS_NA)
         WRSI = gWRSI_NA
         if (.not. (gCalcSOSmode .eqv. .true.)) return
      case (gSOS_NOSTART)
         WRSI = gWRSI_NOSTART ! No start
                              ! Usually occurs after normal start of season has occurred
       ! Analyze/adjust further only those areas where SOSClim itself has a start 
       !  (i.e. does not have a "no-start"):
         if ( (SOSClim /= gSOSCLIM_NA) .and. (SOSClim /= gSOSCLIM_NOSTART) ) then
            if (.not. (sosReached(SOSClim,  gInitialTStep,  &
                                  gCurrentYear, gCurrentTStep, gInitialYear) &
                             .eqv. .true.) ) then
                  WRSI = gWRSI_SEASONYETTOSTART  ! season is yet to start
            !else WRSI = gWRSI_NOSTART (i.e. remains no start)
            endif
         endif
         if (.not. (gCalcSOSmode .eqv. .true.)) return  ! If WRSI runmode, return
      case default

   end select

 ! Only perform computation if these conditions are satisfied:
   proceed =                     &
     (  (LGP > gLGP_PROCWRSIABOVETHISQTY)   .and. &
        (WHC > gWHC_PROCWRSIABOVETHISQTY)   .and. &
        (PET > gPET_PROCWRSIABOVETHISQTY)   .and. &
        ((gCalcSOSmode .eqv. .true.)   .or.       &
         (PPT >= gPPT_PROCWRSIGTETHISQTY))        &
      )

   if (.not. (proceed .eqv. .true.) ) then ! Do nothing
   else
     ! percentSeason:
     ! 0         If season has possible start,
     ! 210       If
     !              season has not started but we want to accumulate, or
     !              season has not yet begun
     ! 200       If season is over
     ! % 0-100   Otherwise
      KF = percentSeason(SOS, int(gInitialYear, 4), gInitialTStep, &
                         int(gCurrentYear, 4), gCurrentTStep, LGP)

      KF2 = 3
      if ( (KF==(0-1)) .or. (KF==(0-2))             ) KF2 = gKF2_SEASONTOSTARTDANDORACCUM
      if (                              (KF==(0-3)) ) KF2 = gKF2_SEASONISOVER
      if (KF2==3) then ! KF is not -1,-2, or -3
         ! Just so that the percent is 0 if negative, 
         ! 1 rather than 0 when the season is exactly at the beginning, 
         ! and the percent otherwise
         if (KF <  0 )  KF2 = gKF2_SEASONHASPOSSIBLESTART
         if (KF == 0 )  KF2 = 1
         if (KF >  0 )  KF2 = CInt2(real((KF * 100),8))
      endif

      if (KF < (0-1)) then  ! Not in growing season nor possible start: Do Nothing
      else
         if (KF < 0) then
         ! Just accumulate ...

            SWAT = SWAT + PPT - (PET * KcFn(real(0,4)))
            SumWR = 0
            SumET = 0
            if (SWAT > WHC) then
               SWAT = WHC
            elseif (SWAT < 0) then
               SWAT = 0
            endif
            SWI = CInt2(real((100 * SWAT / WHC),8))

         else  ! if (KF >= 0)
         ! Now - seasonal calculation:

            Kc = KcFn(KF)

            MAD = WHC * rootfunction(KF) * gCrop_kp
            SW = SWAT + PPT

            ET = AETcalcUsingPETandSoilMoisture(SW, MAD, real(PET, 4), Kc)

            SW = SWAT + PPT - ET   ! SW at the end of the simulation period...

          ! Reset SW to max or min when it is outside bucket size:
            SurplusWater = 0
            if (SW > WHC) then
               SurplusWater = SW - WHC
               SW = WHC
            elseif (SW < 0) then
               SW = 0
            endif
            SurplusWater_TimeStep = SurplusWater ! SY

          ! Store SW back in our array
            SWAT = SW
            SWI = CInt2(real((100 * SWAT / WHC),8))
            SumWR = SumWR + real(PET, 8) * Kc
            WR_TimeStep = real(PET, 8) * Kc    ! SY
            SumET = SumET + ET
            AET_TimeStep = ET     ! SY

          ! At first BW was using a minimum SumWR of 1, but he saw that some PET can be even 0, 
          !  so he decided to limit SumWR to 0. But on the other hand, note that if your PET is 0, 
          !  the crop will probably die of frost (because temperatures will be too low). 
          ! There should be a way to check for this ... 
            if ( (SumWR < gSUMWR_EXITWRSIBELOWTHISQTY) .or. &
               (SumWR > gSUMWR_EXITWRSIABOVETHISQTY)   .or. &
               (SumET > gSUMET_EXITWRSIABOVETHISQTY) ) then  !-9997 to -9999 are program flags

               if (.not. (gCalcSOSmode .eqv. .true.)) return

            else

               if (SumET > SumWR) then
                  WRSI = 100
               else
                  WRSI = CInt2(100 * SumET / SumWR)
                  if( photosensitive .eqv. .true. ) then
                     flrng = LGP * real(gCrop_F2, 8)

                   ! Original VB code had the following troubling direct assignment of 
                   ! double calculation; changed for correctness as follows - BW
                   !  WRSIphoto = WRSI * (flrng + SOSa - 100) / flrng 
                     WRSIphoto = CInt2(WRSI * (flrng + SOSa - 100) / flrng)

                     if (WRSIphoto > 100) WRSIphoto = 100
                     if (WRSIphoto < 0)   WRSIphoto = 0
                     WRSI = WRSIphoto
                  endif
               endif

             ! SY: Begin
             ! Estimate WRSI per timestep (e.g., for each dekad): 
               if (AET_TimeStep > WR_TimeStep) then
                  WRSI_TimeStep = 100
               else
                  WRSI_TimeStep = CInt2(100 * AET_TimeStep / WR_TimeStep)
                  if( photosensitive .eqv. .true. ) then
                     flrng = LGP * real(gCrop_F2, 8)
                     WRSIphoto_TimeStep = CInt2(WRSI_TimeStep * (flrng + SOSa - 100) / flrng)

                     if (WRSIphoto_TimeStep > 100) WRSIphoto_TimeStep = 100
                     if (WRSIphoto_TimeStep < 0)   WRSIphoto_TimeStep = 0
                     WRSI_TimeStep = WRSIphoto_TimeStep
                  endif
               endif
             ! SY: End

            endif

          ! Conditionally set total and maximum surplus water:
            if (SurplusWater > 0) then
               TotalSurplusWater = TotalSurplusWater + SurplusWater
               if (SurplusWater > MaxSurplusWater) MaxSurplusWater = SurplusWater
            endif

            ETminWR = (ET - (PET * real(Kc, 8)))

         !- For WRSI runmode only !
            if (.not. (gCalcSOSmode .eqv. .true.)) then

             ! Conditionally set total and maximum water deficit:
               if (ETminWR < 0) then
                  TotalWaterDeficit = TotalWaterDeficit - ETminWR
                  if ( ((0-1) * ETminWR) > MaxWaterDeficit ) MaxWaterDeficit = ((0-1) * ETminWR)
               endif

             ! Initial Total: AET, WR, surplus water, water deficit
               if (KF < gKF_INIBELOWTHISPCT) &
                  call setTotalAETwrSURPLUSWATERwaterdeficit( &
                          TotalAETInitial, TotalWRInitial,    &
                          TotalSurplusWaterInitial, TotalWaterDeficitInitial, &
                          ET, real(PET, 4), Kc, SurplusWater, ETminWR)

             ! Veg Total: AET, WR, surplus water, water deficit
               if ( (KF >= gKF_VEGGTETHISPCT) .and. (KF < gKF_VEGBELOWTHISPCT) ) &
                  call setTotalAETwrSURPLUSWATERwaterdeficit(&
                          TotalAETVeg, TotalWRVeg, TotalSurplusWaterVeg, &
                          TotalWaterDeficitVeg, ET, real(PET, 4), Kc,    &
                          SurplusWater, ETminWR)

             ! Flower Total: AET, WR, surplus water, water deficit
               if ( (KF >= gKF_FLWRGTETHISPCT) .and. (KF < gKF_FLWRBELOWTHISPCT) ) &
                  call setTotalAETwrSURPLUSWATERwaterdeficit( &
                          TotalAETFlower, TotalWRFlower, TotalSurplusWaterFlower,  &
                          TotalWaterDeficitFlower, ET, real(PET, 4), Kc,   &
                          SurplusWater, ETminWR)

             ! Ripe Total: AET, WR, surplus water, water deficit
               if (KF >= gKF_RIPEGTETHISPCT) &
                  call setTotalAETwrSURPLUSWATERwaterdeficit( &
                          TotalAETRipe, TotalWRRipe, TotalSurplusWaterRipe, &
                          TotalWaterDeficitRipe, ET, real(PET, 4), Kc,      &
                          SurplusWater, ETminWR)

            end if  ! end WRSI run mode only

            ! The following is within 'if KF>=0' in order to avoid checking failure
            ! before crop starts growing or is well established:
            if ( (gPW_EnablePermanentWilting .eqv. .true. ) .and. &
                 (isPermWilted == gISPERMWILTED_NOTWILTED) ) then

               if (gPermanentWiltingIsFromSWI .eqv. .true. ) then
                  if (SWI < gPermanentWiltingThresholdValue) then
                     IsPermWilted = gISPERMWILTED_WILTED
                     PermWiltDate = gCurrentTStep
                     WRSI = 1 ! so that crop does not die and remain with a high WRSI
                     WRSI_TimeStep = 1 ! SY
                  endif
               elseif (gPermanentWiltingIsFromWRSI .eqv. .true.) Then
                  if (WRSI < gPermanentWiltingThresholdValue) then
                        IsPermWilted = gISPERMWILTED_WILTED
                        PermWiltDate = gCurrentTStep
                  endif
               endif
            endif

         endif ! ends if KF<0 else KF>=0
      endif    ! ends if KF>=(0-1)

   endif ! end: only perform computation if these conditions are satisfied:
         !    LGP > 0 .and. WHC > 0 .and. PET > 0 .and. PPT >= 0 

   if (SOS == gSOS_POSSIBLESTART1) then  ! 50 is listed as Possible Start in 
                                         ! '\GeoWRSI\colors\SOSCENA1.CLR'
      WRSI = gWRSI_NOSTART ! No start
      return
   endif
                    
   if (WRSI == gWRSI_NA) then
      if (.not.( sosReached(SOS, gInitialTStep, &
                    gCurrentYear, gCurrentTStep, gInitialYear).eqv. .true.) ) then
            WRSI = gWRSI_SEASONYETTOSTART ! season is yet to start
      else
            WRSI = gWRSI_NOSTART          ! season never started (KRA)
      endif
   endif

 end subroutine CalcWRSI

!
! Function: stepnFromTStep 
! Description: Returns a count of timesteps, 0 if initial timestep; 
!     0+ if initial timestep greater than that
!
 function stepnFromTStep(year, tstep, begYear, begTStep, offset)

      integer*4 :: stepnFromTStep

      integer*4, intent(in) :: year
      integer*4, intent(in) :: tstep
      integer*4, intent(in) :: begYear
      integer*4, intent(in) :: begTStep
      integer*4, intent(in) :: offset

      if (year <  begYear) stepnFromTStep = tstep - gTimeStepsPerYear - begTStep
      if (year == begYear) stepnFromTStep = tstep - begTStep
      if (year >  begYear) stepnFromTStep = gTimeStepsPerYear + tstep - begTStep

      ! stepnFromTStep = stepnFromTStep + 1   
      ! TTM changed this subroutine so that it doesn't give distance + 1. 

      stepnFromTStep = stepnFromTStep + offset

 end function stepnFromTStep


 function sosReached(SOS,   initialTStep,   currentYear, currentTStep, initialYear)

   logical :: sosReached

   integer*4, intent(in) :: SOS
   integer*4, intent(in) :: initialTStep
   integer*2, intent(in) :: currentYear
   integer*4, intent(in) :: currentTStep
   integer*2, intent(in) :: initialYear

   integer*4 :: addYear
   integer*4 :: dist
! _____________________________________________

   addYear = 0

  ! Check whether the 'no start' data are really still due to start
  ! i.e., the SOS tstep has not yet been reached and mark these
  !  pixels also as yet to start.
   if (SOS < initialTStep) then
                          ! Note: this can only work
                          ! if your initial tstep occurs quite late.
                          ! I.e., if it occurs early in the year
                          ! this condition will fail
      addYear = 1
   endif

   dist = stepnFromTStep( int(currentYear,4), currentTStep, &
          (initialYear + addYear), SOS, 0 )
   if ( dist < 0 ) then    ! <== This was 1 before I changed ndeksfromdekad. TTM
      sosReached = .false. ! Start of season has not yet been reached
      return
   endif
   sosReached = .true.     ! Start of season has been reached

 end function sosReached


 function percentSeason(SOS, initialYear, initialTStep, currentYear, currentTStep, LGP)

   real*4 :: percentSeason

   integer*4, intent(in) :: SOS
   integer*4, intent(in) :: initialYear
   integer*4, intent(in) :: initialTStep
   integer*4, intent(in) :: currentYear
   integer*4, intent(in) :: currentTStep
   integer*4, intent(in) :: LGP

   integer*4, parameter :: sos_PossibleStart0 = gTimeStepsPerYear + 1
   integer*4 :: tStep_n
   integer*4 :: iYear

   select case (SOS)
      case (gSOS_NA)
         percentSeason = (0.0-2.0)
         return
      case (sos_PossibleStart0:gSOS_POSSIBLESTART1)
         percentSeason = (0.0-1.0)
         return
      case (gSOS_NOSTART)
         percentSeason = (0.0-3.0)
       ! The following will tell the program that it has not started but should accumulate
         if (gCalcSOSmode .eqv. .true.) percentSeason = (0.0-1.0)
         return
      case default
   end select

   iYear = initialYear ! if (SOS >= initialTStep)
   if (SOS < initialTStep) iYear = initialYear + 1

   tStep_n = stepnFromTStep(currentYear, currentTStep, iYear, SOS, 0)   !# timesteps from SOS
   select case (tStep_n)
      case ( (0-gTStepsBeforeSeasonStarts):(0-1) )  ! we accumulate SWAT
         percentSeason = (0.0-1.0)
      case ( :(0-(gTStepsBeforeSeasonStarts+1)) )
         percentSeason = (0.0-2.0)
      case default
         if ( tStep_n > LGP ) then
            percentSeason = (0.0-3.0)
         else
            percentSeason = (real(tStep_n, 4) + 1.0) / LGP
            if ( percentSeason > 1.0 ) percentSeason = 1.0
         endif
   end select

 end function percentSeason


 function KcFn(KF)

   real*4 :: KcFn
   real*4, intent(in) :: KF

!- Set crop coefficient:
   if     (KF < 0.0)         then
      KcFn = gPreSeason_Kc

   elseif (KF <= gCrop_F1)   then
      KcFn = gCrop_K1

   elseif (KF <= gCrop_F2)   then
      KcFn = gCrop_K1 + (gCrop_K2 - gCrop_K1) / (gCrop_F2 - gCrop_F1) * (KF - gCrop_F1)

   elseif (KF <= gCrop_F3)   then
      KcFn = gCrop_K2

   elseif (KF <= 1.0)        then
      KcFn = gCrop_K2 + (gCrop_K3 - gCrop_K2) / (1.0 - gCrop_F3) * (KF - gCrop_F3)

   else
      KcFn = gCrop_K3

   endif

 end function KcFn


 function rootfunction(KFi)

   real*4 :: rootfunction

   real*4, intent(in) :: KFi

   real*4 :: KF

   KF=KFi
   if (KF < 0.0) KF = 0.0

    ! rinit is initial root depth at emergence ~ 0.1 for most crops
      if (KF <= gCrop_F2) then
        !the following changes in the algorithm by Senay in July 2004
         rootfunction = gCrop_rini + (1.0 - gCrop_rini) * KF / gCrop_F2
         !rootfunction = gCrop_rini + (gCrop_r - gCrop_rini) * KF / gCrop_F2
      else
         rootfunction = 1.0
        !rootfunction = gCrop_r
     endif

 end function rootfunction


 function AETcalcUsingPETandSoilMoisture(SWi, MAD, PET, Kc)

   real*4 :: AETcalcUsingPETandSoilMoisture

   real*4, intent(in) :: SWi
   real*4, intent(in) :: MAD
   real*4, intent(in) :: PET
   real*4, intent(in) :: Kc

   real*4 :: ET1
   real*4 :: ET2   
   real*4 :: SW

   ET2 = 0   ! added by BW
   SW = SWi

   if (SW >= MAD) then
      ET1 = PET * Kc
      if (ET1 > SW) ET1 = SW
      SW = SW - ET1
      if (SW >= MAD) then
         ET2 = PET * Kc
         if (ET2 > SW) ET2 = SW
      elseif ( (SW < MAD) .and. (SW > 0.0) ) then
         ET2 = SW / MAD * PET * Kc
         if (ET2 > SW) ET2 = SW
      else
         ET2 = 0.0
      endif
   else
      ET1 = SW / MAD * PET * Kc
      if (ET1 > SW) ET1 = SW
      SW = SW - ET1
      if (SW > 0.0) then
         ET2 = SW / MAD * PET * Kc
         if (ET2 > SW) then
            ET2 = SW
         else
            ET2 = 0.0
         endif
      endif
   endif

   AETcalcUsingPETandSoilMoisture = &
      gCrop_c1 * ET1 + gCrop_c2 * ET2  ! c1/c2 0.75/.25 or 0/1 for rice
   ! c1 is found at the initialization subroutine
   ! ET calculated

 end function AETcalcUsingPETandSoilMoisture


 subroutine setTotalAETwrSURPLUSWATERwaterdeficit(      &
               tAET, tWR, tSurplusWater, tWaterDeficit, &
               ET, PET, Kc, SurplusWater, ETminWR  )

   real*4, intent(inout) :: tAET
   real*4, intent(inout) :: tWR
   real*8, intent(inout) :: tSurplusWater
   real*8, intent(inout) :: tWaterDeficit
   real*4, intent(in) :: ET
   real*4, intent(in) :: PET
   real*4, intent(in) :: Kc
   real*8, intent(in) :: SurplusWater
   real*8, intent(in) :: ETminWR

   tAET = tAET + ET
   tWR = tWR + PET * Kc
   if (SurplusWater > 0) tSurplusWater = tSurplusWater + SurplusWater
   if (ETminWR < 0)      tWaterDeficit = tWaterDeficit - ETminWR

 end subroutine setTotalAETwrSURPLUSWATERwaterdeficit


 subroutine CalcWRSIAnomalies( Mask, LGP, WRSIa, WRSI, SWI, KF2, WRSIClim )

! "Point" or "pixel-level" code of the WRSI anomalies computation.  Called 
! together with CalcWRSI, immediately following the call to it (i.e., this 
! subroutine must execute immediately following execution of the WRSI 
! computation subroutine).

! "Calculate anomalies as % of geowrsi2 climalogies.  This subroutine also 
!  takes care of other parameters to give them the gWRSI_NOSTART(e.g. 253) 
!  and gWRSI_SEASONYETTOSTART(e.g. 254) coding."

   integer*2, intent(in)    :: Mask              ! Byte    maskData(x,y)
 ! Increased numerical prec of following variables to allow flexibility 
 !  of design by bringing everything back to seconds
   integer*4, intent(in)    :: LGP               ! Byte    lgpData(x,y)
   integer*2, intent(inout) :: WRSIa             ! Byte    geowrsi2aData(x,y)
   integer*2, intent(in)    :: WRSI              ! Byte    geowrsi2Data(x,y)
   integer*2, intent(inout) :: SWI               ! Byte    swiData(x,y)
   integer*2, intent(inout) :: KF2               ! Byte    percentSeasonData(x,y)
   integer*2, intent(in)    :: WRSIClim          ! Byte    geowrsi2Clim(x,y)
   
   real*4 :: t1, t2, t3
! _______________________________________

   if ( (Mask /= gMASK_INC) .or. (LGP < 1) ) then
!#ifndef MULTIPLE_CROPS_
      if (Mask /= gMASK_INC) WRSIa = gWRSI_NA
!#endif
      if (LGP < 1)           WRSIa = gWRSI_NA
      return
   endif 

 ! gWRSI_NOSTART:  no start (late) ; e.g. 253 for dekad tStep
 ! gWRSI_SEASONYETTOSTART: yet to start; e.g. 254 for dekad tStep
 ! gWRSI_NA:  N/A (no data value for WRSI); e.g. 255 for dekad tStep
   select case (WRSI)
      case (gWRSI_NOSTART:gWRSI_SEASONYETTOSTART)
         WRSIa = WRSI
         SWI   = WRSI
         KF2   = WRSI
      case (gWRSI_NA)
         WRSIa = WRSI
         SWI   = WRSI
         KF2   = 0 
      case default
         t1 = real(WRSI,4)
         t2 = real(WRSIClim,4)
         if (t2 == 0) t2 = gWRSI_NA
         t3 = 100.0 * t1 / t2
         if (t3 > gWRSI_MAXANOMALYPCT) t3 = gWRSI_MAXANOMALYPCT 
         if ( (t2 > gWRSI_MAXNUMDENOMANOMALY) .or. &
              (t1 > gWRSI_MAXNUMDENOMANOMAlY) ) t3 = gWRSI_NA
         WRSIa = CInt2(real(t3,8))
   end select

 end subroutine CalcWRSIAnomalies

 function DifferenceOf2TSteps(minuendYear, minuendTStep, subtrahendYear, &
                              subtrahendTStep, timestepsPerYear)

! Subtracts one date specified as Year/TStep (subtrahendYear/TStep) 
!  from another (minuendYear/TStep) the resulting difference returned is in tsteps.

   implicit none

   integer*4 :: DifferenceOf2TSteps

   integer*2, intent(in)  :: minuendYear
   integer*4, intent(in)  :: minuendTStep
   integer*2, intent(in)  :: subtrahendYear
   integer*4, intent(in)  :: subtrahendTStep
   integer*4, intent(in)  :: timeStepsPerYear

   DifferenceOf2TSteps = (minuendYear - subtrahendYear) * timeStepsPerYear &
                       + (minuendTStep - subtrahendTStep)

 end function DifferenceOf2TSteps


!== Main CALC SOS Routines ==

 subroutine offsetTStepByNumTSteps (TStep, offsetNumTSteps, Year)

! Adds offsetNumTSteps to TStep, adjusts Year in accordance, if necessary

   integer*4, intent(inout) :: TStep
   integer*4, intent(in)    :: offsetNumTSteps
   integer*2, intent(inout), optional :: Year

   TStep = TStep + offsetNumTSteps

   if (TStep <= 0) then
      TStep = TStep + gTimeStepsPerYear
      if(present(Year)) Year = Year - 1
   endif
   if (TStep > gTimeStepsPerYear) then
      TStep = TStep - gTimeStepsPerYear
      if(present(Year)) Year = Year + 1
   endif

 end subroutine offsetTStepByNumTSteps


 subroutine setWiltings(wiltingTStep, Wilting1, Wilting2)

   integer*4, intent(in)    :: wiltingTStep
   integer*4, intent(inout) :: Wilting1
   integer*4, intent(inout) :: Wilting2

   if (Wilting1 == gWILTING1_DEFAULT) then   ! if there has been no wilting before:
      Wilting1 = wiltingTStep
   else   ! If there has already been a wilting once before:
      if (Wilting2 == gWILTING2_DEFAULT) Wilting2 = wiltingTStep
   endif

 end subroutine setWiltings


 subroutine CalcSOS( Mask, SOSClim, PPT, WRSI, KF2,   &
                SWI, Wilting1, Wilting2, resetWBdata, &
                SOS, SOSa, attemptPostSOScalcSOSreset )
!
!  Routine:  CalcSOS
!
!  Description:
!   Typically, precipitation estimates are used to identify the period of 
!    onset of rains (or start-of-season) for the agricultural growing season. 
!    Cumulative rainfall accounting criteria are applied to the rainfall 
!    estimates on a per pixel basis.  Beginning several periods in advance 
!    of the usual onset of rains, each pixel is tested to identify the first 
!    10-day period (e.g., in which at least 25 mm) of rain has fallen. 
!    Onset of rains is then identified by the following criteria: 
!
!       1st 10-day period in which at least 25 mm of rain falls, 
!       followed by at least 20 mm of rain in the next two dekadal updates. 
!
!    If rainfall for the latter two updates totals at least 20 mm of rain, 
!    the 1st 10-day period of 25 mm rainfall is identified as the period of 
!    onset of rains, and the start of the agricultural growing season. 
!    If the latter two update periods do not total at least 20 mm of rain, 
!    the 1st 10-day period is declared a "failed planting", and the search 
!    for onset of rain, or planting period, continues with the next 10-day 
!    accumulation period. 
!
!   The onset of rains algorithm applied is from an algorithm originally 
!     developed for a dekadal (~10-day) time step for the Sahel region of 
!     West Africa (AGRHYMET, 1996).  
!
!  References:
!    http://earlywarning.cr.usgs.gov/fews/ccm/centralamerica/web/readme.php?symbol=p3
!   
!    AGRHYMET, 1996. Methodologie de suivi des zones a risque. AGRHYMET FLASH, 
!     Bulletin de Suivi de la Campagne Agricole au Sahel, Centre Regional AGRHYMET, 
!     B.P. 11011, Niamey, Niger, Vol. 2, No. 0/96, 2 pages.
!
! _________________________________________________________________________________

   integer*2, intent(in)               :: Mask       ! maskData(x,y)
   integer*4, intent(in)               :: SOSClim    ! sosClim(x,y)
   integer*2, dimension(:), intent(in) :: PPT        ! pptCube(x,y)
   integer*2, intent(inout)            :: WRSI       ! geowrsi2Data(x,y)
   integer*2, intent(inout)            :: KF2        ! percentSeasonData(x,y)
   integer*2, intent(in)               :: SWI        ! swiData(x,y)
 ! Increased numerical prec of the following two variables which were short 
 !  to which what needed to be single, per soscalc, was assigned:
   integer*4, intent(inout)            :: Wilting1   ! Wilting1Data(x,y) 
   integer*4, intent(inout)            :: Wilting2   ! Wilting2Data(x,y)
   logical,   intent(inout)            :: resetWBdata
   integer*4, intent(inout)            :: SOS        ! sosData(x,y)
   integer*4, intent(inout)            :: SOSa       ! sosaData(x,y)
   logical,   intent(in)               :: attemptPostSOScalcSOSreset

   integer*4 :: no_start_val
   logical   :: SOSrequirementSatisfied 
   integer*4 :: diffFromLastCurr
   logical   :: allRainSoFarIsCurrent

   integer*4 :: aclimsosTStep   ! the date(expressed as the timestep) of the clim. SOS
   integer*2 :: ansosclimYear   ! the year of the climatological SOS 
   integer*4 :: climDiff
! ________________________________________________________________________

   resetWBdata = .false.
   no_start_val = gSOS_NOSTART

 ! This following default value was assigned not in this routine but outside of it 
 ! and outside of the loop that runs it; i.e. once for all (spatial) points. +BW

   SOSrequirementSatisfied = .false.

   if (Mask /= gMASK_INC) return

   diffFromLastCurr = DifferenceOf2TSteps( gLastCurrentYear, gLastCurrentTStep, &
                                           gCurrentYear, gCurrentTStep, gTimeStepsPerYear )

 ! Check if the water (rainfall, swi, etc) data satisfies the conditions for an SOS:
   if( gSOScalcUseRainThreshold .eqv. .true. ) then
       SOSrequirementSatisfied = .false.

      if( (diffFromLastCurr >= 2 )            .and. &
          (PPT(1) >= gSOScalcTStep1Threshold) .and. &
          ((PPT(2) + PPT(3)) >= gSOScalcTStep2Plus3Threshold) ) then
         SOSrequirementSatisfied = .true.
      endif

    ! The following set of statements is added to allow SOS to be activated 
    !   even before it's confirmed.
    ! ==== Start SOS preliminary ok modification ====
      if (.not. (SOSrequirementSatisfied .eqv. .true.) ) then
         if ((diffFromLastCurr >= gSOSPREOK_DIFLSTCURCURGTETHIS)   .and. &
             (diffFromLastCurr <  gSOSPREOK_DIFLSTCURCURBELOWTHIS) .and. &
             (PPT(1) >= gSOScalcTStep1Threshold) ) then
            SOSrequirementSatisfied = .true.
         endif
      endif
    ! ==== End SOS preliminary ok modification ====

   elseif (gSOScalcUseSWIThreshold .eqv. .true.) then

      allRainSoFarIsCurrent = .true.

!!! BW LastCurrent edit begin !!!
#if 0
    ! Here, '2' initializes the index being in fortran +BW; 
      do i = 2, timeToProcess
         if (.not. (PPTisAnnual(i) .eqv. .true.)) allRainSoFarIsCurrent = .false. 
      end do
#endif
    ! Although I may not understand it fully, I took the tack here of re-
    ! implementing the logic as-is w/ exclusively 'last current' information +BW

      call offsetTStepByNumTSteps( gInitialTstep, &
                 (0-gTStepsBeforeSeasonStarts), gInitialYear )

      if( DifferenceOf2TSteps(gLastCurrentYear, gLastCurrentTStep, &
                              gInitialYear, gInitialTStep, &
                              gTimeStepsPerYear) /= 0 ) then
         if( diffFromLastCurr <= 0 ) then
            allRainSoFarIsCurrent = .false.
         endif
      endif
!!! BW LastCurrent edit end   !!

      if (allRainSoFarIsCurrent .eqv. .true.) then
         SOSrequirementSatisfied = .false.

         if (SWI >= gSOScalcSoilMoistureThreshold) then
            if (.not. (gSOScalcUseNearestToDesiredPracticalPlantingTStep .eqv. .true.)) then
               SOSrequirementSatisfied = .true.
            else
             ! The following desired planting time check would seem to be most complete if it 
             ! incorporated not just timestep but also year as well (??); but this is how the 
             ! original was coded and it is assumed that the original was in fact correct. +BW
               if ( (gCurrentTStep >= gSOScalcDesiredPlantingTStep) .and. &
                  ( (gCurrentTStep - gSOScalcDesiredPlantingTStep)  <     &
                     gSOSCALCUSENRSTTODSRDPRACPLNTTSTEPTHRESHOLD) ) then
                  SOSrequirementSatisfied = .true.
               endif
            endif
         endif
      endif

   endif

 ! If the time-of-season (i.e. >= gTStepsBeforeSeasonStarts of PPT) has started,
 !  check if the crop has failed, if necessary:
   if( SOS /= no_start_val ) then
     if( &
        ((gCurrentYear > gInitialYear) .or. (gCurrentTStep >= gInitialTStep) )  .and. &
        (((gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart .eqv. .true.) .and. &
          (WRSI <= gSOScalcRestartThresholdWithWRSI) ) &
         .or. &
         ((.not. (gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart .eqv. .true.)) .and. &
          (SWI <= gSOScalcRestartThresholdWithSWI) &
         ) & 
        ) &
       ) then 

       call setWiltings(gCurrentTStep, Wilting1, Wilting2)

       if( (gSOScalcRestartCropWhenFailed .eqv. .true.)           .and. &
           (KF2 <= gSOScalcCropCanStillRestartWithinPctLGPThreshold)    &
         ) then    ! Note: 'KF2' is 'percentSeason'
          resetWBdata = .true.  ! Set to have reset occur after execution of this routine
          WRSI = 100
          KF2  = 0 ! As of here it is determined that KF2 is to be reset and it is 
                   !  important to do so here b/c it is potentially referenced further 
                   !  down during same call to this soscalc routine; same for WRSI
          SOS = no_start_val
       endif
    endif !!!

   endif

 ! Determine the SOS climatology date:
 ! Do the calculation only if all the tsteps have realtime data (next line):
   if( (SOSClim == gSOSCLIM_NOSTART) .or. (SOSClim == gSOSCLIM_NA) ) then   
     ! SOSClim is climatological SOS.
     ! So the SOS can get a value even when the clim SOS has "no planting" indicated.
     ! Therefore, if the pixel has a climatological SOS value of 60, a fake SOS 
     !  anomaly will be assigned equal to average.

#if 0
!8/2*!
      aclimsosTStep = stepnFromTStep(int(gCurrentYear, 4), gCurrentTStep, &
                                     int(gInitialYear, 4), gInitialTStep, &
                                     gTStepsBeforeSeasonStarts)
#endif
      ansosclimYear = gCurrentYear 
      aclimsosTStep = gCurrentTStep 

   else   ! If SOSClim /= 0 then

    ! This is the default setting, and will work whether or not the season crosses years:
      ansosclimYear = gInitialYear   
      if( gInitialTStep > gFinalTStep ) then   ! The WRSI run crosses Dec. 31st
         if( SOSClim < gInitialTStep ) ansosclimYear = gFinalYear
      endif
      !8/2*!aclimsosTStep = stepnFromTStep(int(ansosclimYear, 4), SOSClim, &
      !                        int(gInitialYear, 4), gInitialTStep, gTStepsBeforeSeasonStarts)
      aclimsosTStep = SOSClim 

   endif

  !8/2*!climDiff = DifferenceOf2TSteps(gCurrentYear, gCurrentTStep, &
  !                                    ansosclimYear, aclimsosTStep, gTimeStepsPerYear)
   climDiff = DifferenceOf2TSteps(gCurrentYear, gCurrentTStep, ansosclimYear, &
                                  aclimsosTStep, gTimeStepsPerYear) 

  ! if(diffFromLastCurr == 0) SOSrequirementSatisfied = .false.  ! print*, "...at the last ts!!" 
  !  print*, 'SOSrequirementSatisfied=', SOSrequirementSatisfied, 'gCurrentTStep=',gCurrentTStep

   if (.not. (SOSrequirementSatisfied .eqv. .true.)) then

      if ( (SOSa == 0) .and. &
           (diffFromLastCurr == 0) ) then 
         SOSa = 100 - climDiff
         if (climDiff <= 0) then
            SOSa = 150
         endif
      endif

   else
   ! if (SOSrequirementSatisfied .eqv. .true.) then

     !  logical   :: gSOScalcRestartCropWhenFailed = .true.
     !  logical   :: gSOScalcUseWRSIThresholdInsteadofSWIThresholdToRestart = .true.

     ! The following threshold allows you to determine whether the crop has failed so 
     !  that you can restart using normal methods.
     !  integer*4 :: gSOScalcRestartThresholdWithWRSI = 30   

      if( (gCurrentYear > gInitialYear) .or. (gCurrentTStep >= gInitialTStep) ) then
       ! Season potential start time reached:
         if (gSOScalcIgnoreClimatology .eqv. .true.) then
          ! ... so that it doesn't put a new value if it has already started.
            if (SOS == no_start_val) then   
               SOS = gCurrentTStep
               if (SOS > gTimeStepsPerYear) SOS = SOS - gTimeStepsPerYear 
               SOSa = 100   ! SOSa is SOS-Anomaly data
            endif
         else
            if ( (climDiff <= gSOScalcMaxTStepsLate) .and. &
                 (climDiff >= ((0-1) * gSOScalcMaxTStepsEarly)) .and. &
                 (SOS == no_start_val) ) then   
             ! So that it doesn't put a new value if it's already started.
               SOS = gCurrentTStep
               SOSa = 100 - climDiff
               if (SOSClim == gSOSCLIM_NOSTART) SOSa = 100
            endif
         endif
      endif

   endif

   if(attemptPostSOScalcSOSreset .eqv. .true.) then
      if (gSOScalcExcludeIncompleteAreasFromSOS .eqv. .true.) then
         if (KF2 < gSOScalcAcceptablePercentOfTotalSeasonalGrowth) then
          ! Note: 'KF2' is 'percentSeason'
            SOS  = no_start_val
            SOSa = no_start_val
          !print*,'resetting SOS gCurrentTStep=',gCurrentTStep
         endif
      endif

    ! The condition below guards against the possibility of the SOS-anomaly continuously 
    !  having "yet to start" if the SOS-clim value is 0 where mask is 1. (T.Magadzire, March 17, 2011)
      if ((SOS == no_start_val) .and. (SOSClim == 0)) then
         SOSa = 0
      endif
   endif

 end subroutine CalcSOS


! The following notes were for if/when CalcSOS was part of init of same process -Jun
! After output for a crop for a timestep as desired (complete being the default).
! For resetWBdata, call following in masked spatial loop w/ exceptSoilWater = .false.

 subroutine resetWBdataFor1point( &
      exceptSoilWater,                                      &
      isPermWilted, WRSI, KF, KF2, SWAT, SumWR, SumET, SWI, &
      MaxSurplusWater, MaxWaterDeficit,                     &
      TotalSurplusWaterInitial, TotalWaterDeficitInitial,   &
      TotalSurplusWaterVeg, TotalWaterDeficitVeg,           &
      TotalSurplusWaterFlower, TotalWaterDeficitFlower,     &
      TotalSurplusWaterRipe, TotalWaterDeficitRipe,         &
      PermWiltDate,                                         &
      TotalSurplusWater,                                    &
      TotalWaterDeficit,                                    &
      TotalAETInitial,                                      &
      TotalWRInitial,                                       &
      TotalAETVeg,                                          &
      TotalWRVeg,                                           &
      TotalAETFlower,                                       &
      TotalWRFlower,                                        &
      TotalAETRipe,                                         &
      TotalWRRipe,                                          &
      WRSIphoto,                                            &
      SOSa ) 

   logical,  intent(in)    :: exceptSoilWater
   integer*2,intent(inout) :: isPermWilted      ! Byte    isPermWilted(x,y) ! 0=not wilted, 1=wilted
   integer*2,intent(inout) :: WRSI              ! Byte    geowrsi2Data(x,y)
   real*4,   intent(inout) :: KF                ! Single  fracSeasonData(x,y)
   integer*2,intent(inout) :: KF2               ! Byte    percentSeasonData(x,y)
   real*4,   intent(inout) :: SWAT              ! Single  swatData(x,y)
   real*8,   intent(inout) :: SumWR             ! Double  swrData(x,y)
   real*8,   intent(inout) :: SumET             ! Double  setData(x,y)
   integer*2,intent(inout) :: SWI               ! Byte    swiData(x,y)
   real*8,   intent(inout) :: MaxSurplusWater   ! Single  maxSurplusWater(x,y)
   real*8,   intent(inout) :: MaxWaterDeficit   ! Single  maxWaterDeficit(x,y)
   real*8,   intent(inout) :: TotalSurplusWaterInitial, TotalWaterDeficitInitial
                                                 ! Single totalSurplusWaterInitial(x,y)
                                                 ! Single totalWaterDeficitInitial(x,y)
   real*8,   intent(inout) :: TotalSurplusWaterVeg, TotalWaterDeficitVeg
                                                 ! Single totalSurplusWaterVeg(x,y)
                                                 ! Single totalWaterDeficitVeg(x,y)
   real*8,   intent(inout) :: TotalSurplusWaterFlower, TotalWaterDeficitFlower
                                                 ! Single totalSurplusWaterFlower(x,y)
                                                 ! Single totalWaterDeficitFlower(x,y)
   real*8,   intent(inout) :: TotalSurplusWaterRipe, TotalWaterDeficitRipe
                                                 ! Single totalSurplusWaterRipe(x,y)
                                                 ! Single totalWaterDeficitRipe(x,y)
   integer*4, intent(inout) :: PermWiltDate      ! Byte   PermWiltDate(x,y)

 ! The following were not in the original routine but added for completeness, 
 ! especially anticipating multi-season runs.
   real*8, intent(inout), optional    :: TotalSurplusWater
   real*8, intent(inout), optional    :: TotalWaterDeficit
   real*4, intent(inout), optional    :: TotalAETInitial, TotalWRInitial
   real*4, intent(inout), optional    :: TotalAETVeg, TotalWRVeg
   real*4, intent(inout), optional    :: TotalAETFlower, TotalWRFlower
   real*4, intent(inout), optional    :: TotalAETRipe, TotalWRRipe
   integer*2, intent(inout), optional :: WRSIphoto
   integer*4, intent(inout), optional :: SOSa

! ___________________________________________________________________________________________

   IsPermWilted = 0
   WRSI = 100
   KF = 0
   KF2 = 0

 ! SWAT is not [normally] reset because it is the result
 !  of several previous timesteps of monitoring.
   if (.not. (exceptSoilWater .eqv. .true.)) then
      SWAT = 0
      SWI = 0

    ! The following weren't in the original routine but added for completeness, 
    !  especially anticipating multi-season runs.
      if(present(TotalSurplusWater)) TotalSurplusWater = 0
      if(present(TotalWaterDeficit)) TotalWaterDeficit = 0
      if(present(TotalAETInitial))   TotalAETInitial = 0
      if(present(TotalWRInitial))    TotalWRInitial = 0
      if(present(TotalAETVeg))       TotalAETVeg = 0
      if(present(TotalWRVeg))        TotalWRVeg = 0
      if(present(TotalAETFlower))    TotalAETFlower = 0
      if(present(TotalWRFlower))     TotalWRFlower = 0
      if(present(TotalAETRipe))      TotalAETRipe = 0
      if(present(TotalWRRipe))       TotalWRRipe = 0
      if(present(WRSIphoto))         WRSIphoto = 0
!      if(present(SOSa)) SOSa = 0   ! KRA: To enable writing of SOSa in WRSI-mode
   endif

   SumWR = 0
   SumET = 0
   MaxSurplusWater = 0
   MaxWaterDeficit = 0
   TotalSurplusWaterInitial = 0
   TotalWaterDeficitInitial = 0
   TotalSurplusWaterVeg = 0
   TotalWaterDeficitVeg = 0
   TotalSurplusWaterFlower = 0
   TotalWaterDeficitFlower = 0
   TotalSurplusWaterRipe = 0
   TotalWaterDeficitRipe = 0
   PermWiltDate = 0

end subroutine resetWBdataFor1point

! End CALCSOS set of code

!___________________

end module geowrsi2_physics_module
