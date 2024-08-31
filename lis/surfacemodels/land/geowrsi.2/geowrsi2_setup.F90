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
! !ROUTINE: geowrsi2_setup
! \label{geowrsi2_setup}
!
! !REVISION HISTORY:
!
! 31 Jul 2011: Brad Wind; Initial Definition
! 21 Dec 2012: KR Arsenault; Cleaned up code and comments
! 25 Oct 2013: KR Arsenault; Added GeoWRSI2.0 model to LIS-7
! 
! !INTERFACE:
subroutine geowrsi2_setup()

! !USES:
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_logunit, LIS_endrun

  use geowrsi2_module
  use geowrsi2_lsmMod
  use geowrsi2_physics_module, only : offsetTStepByNumTSteps, &
           gSOScalcNtStepsOfForcingUsed, &
           gPPTDATALIMITEDORMISSINGTHRESHOLD, gTimeStepsPerYear, &
           gWILTING1_DEFAULT, gWILTING2_DEFAULT, gMASK_EXC, &
           gSOS_INITIAL, gSOS_NOSTART, gCurrentYear, gCurrentTStep, &
           DifferenceOf2TSteps

  use fbil_module

!
! !DESCRIPTION:
!
!  This routine is the entry point to complete the set up parameters
!  required for the FEWSNET WRSI LSM. These include everything from crop phenology 
!  water requirement parameterizations, to soil water holding capacities, 
!  to start and end of crop growing season, to analyst options and state variables 
!  initializations in WRSI.
! 
!EOP
  implicit none

  integer   :: n, t, Tm
  integer*2 :: earliestLastYear
  integer*4 :: earliestLastTstep
  type(geowrsi2dec), pointer :: geowrsi2Pt
! _____________________________________

  write(LIS_logunit,*) "Initializing WRSI LSM variables (geowrsi2_setup)"

  do n = 1, LIS_rc%nnest
    do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)

       geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)

    !- Initialize PPT and PET fields for SOS calculation:
       do Tm = 1, gSOScalcNtStepsOfForcingUsed
          geowrsi2Pt%PPTcacheForSOScalc(Tm) = CInt2(real(geowrsi2_udef, 8))
          geowrsi2Pt%PETcacheForSOScalc(Tm) = CInt2(real(geowrsi2_udef, 8))
       end do

    !- Initialize output variables:
       geowrsi2Pt%WR_TimeStep   = LIS_rc%udef          ! "TimeStep" terms +SY
       geowrsi2Pt%AET_TimeStep  = LIS_rc%udef          
       geowrsi2Pt%WRSI_TimeStep = LIS_rc%udef          
       geowrsi2Pt%SurplusWater_TimeStep = LIS_rc%udef  
       geowrsi2Pt%PPT     = 0
       geowrsi2Pt%PET     = 0
       geowrsi2Pt%ACC_PPT = 0.0
       geowrsi2Pt%ACC_PET = 0.0

     ! Initialize LIS's output variables prior to the start of season. +BW (Jan. 2012)
       geowrsi2Pt%SumWR   = 0
       geowrsi2Pt%SumET   = 0
       geowrsi2Pt%TotalSurplusWater = 0
       geowrsi2Pt%MaxSurplusWater   = 0
       geowrsi2Pt%TotalWaterDeficit = 0
       geowrsi2Pt%MaxWaterDeficit   = 0
       geowrsi2Pt%TotalAETInitial   = 0
       geowrsi2Pt%TotalWRInitial    = 0
       geowrsi2Pt%TotalSurplusWaterInitial = 0
       geowrsi2Pt%TotalWaterDeficitInitial = 0
       geowrsi2Pt%TotalAETVeg = 0
       geowrsi2Pt%TotalWRVeg  = 0
       geowrsi2Pt%TotalSurplusWaterVeg = 0
       geowrsi2Pt%TotalWaterDeficitVeg = 0
       geowrsi2Pt%TotalAETFlower = 0
       geowrsi2Pt%TotalWRFlower  = 0
       geowrsi2Pt%TotalSurplusWaterFlower = 0
       geowrsi2Pt%TotalWaterDeficitFlower = 0
       geowrsi2Pt%TotalAETRipe   = 0
       geowrsi2Pt%TotalWRRipe    = 0
       geowrsi2Pt%TotalSurplusWaterRipe = 0
       geowrsi2Pt%TotalWaterDeficitRipe = 0
       geowrsi2Pt%SWI  = 0
       geowrsi2Pt%isPermWilted  = 0
       geowrsi2Pt%LGP_TIMESTEPS = 0
       geowrsi2Pt%WHC  = 0
       geowrsi2Pt%Mask = 0
       geowrsi2Pt%WRSI = 0
       geowrsi2Pt%KF   = 0
       geowrsi2Pt%KF2  = 0
       geowrsi2Pt%SWAT = 0
       geowrsi2Pt%WRSIphoto = 0
       geowrsi2Pt%SOSa = 0
       geowrsi2Pt%PermWiltDate = 0
       geowrsi2Pt%Wilting1 = 0
       geowrsi2Pt%Wilting2 = 0
       geowrsi2Pt%wrsia = 0
    end do
  end do

!- Set up reading in WRSI User Input Settings (File): 
   call geowrsi2_readInputSettings()


!- SOS Calc Runmode:
  if( geowrsi2_CalcSOSlsmRunMode .eqv. .true. ) then
    write(LIS_logunit,*) "Initializing WRSI SOS Calc Mode Options (geowrsi2_setup)"

    do n = 1, LIS_rc%nnest
      do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
         geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)
       ! Calculate SOS only over unmasked areas, omitting areas where there is a mask.
         if( (geowrsi2Pt%Mask == CInt2(real(geowrsi2_udef, 8))) .or. &
             (geowrsi2Pt%Mask == gMASK_EXC) ) then   ! 0
            geowrsi2Pt%SOS = gSOS_INITIAL   ! 0
         else
            geowrsi2Pt%SOS = gSOS_NOSTART   ! 60
         endif
      end do
    end do

  endif  ! End SOS field initialization

end subroutine geowrsi2_setup
