!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! BOP
! 
! !ROUTINE: irrigation_model
! \label{irrigation_model}

! !DESCRIPTION:        
!
! Irrigation model for a single or combination of irrigation types 
! (Sprinkler, Drip, Flood, and Rice paddy) that supports calls from 
! different LSMs.  
! The model consists of four main components with multiple configurable options:
! 1. Crop growing season determination (vegetation seasonality, crop calendar)
! 2. Irrigation trigger (soil moisture deficit, schedule)
! 3. Irrigation rate (dynamic, prescribed)
! 4. Irrigation application method (precipitation, soil moisture content)
!
! Growing season is determined by greenness fraction or LAI seasonality using
! a threshold given the annual range of max/min values.  Alternatively, crop
! specific planting/harvesting dates can be obtained from a crop calendar.
!
! For all irrigation types, the irrigation season starts once it is determined 
! to be in the crop’s growing season, and when root zone soil moisture 
! availability falls below a threshold (e.g. 50 %) of the field capacity 
! (reference soil moisture).  These deciding factors for when to irrigate is 
! called "irrigation trigger". The water requirement in the root zone is based 
! on the actual maximum root depth rather than the root zone parameter in Noah.
!
! Irrigation trigger is checked at a specified time (e.g. 6 am LST) for the 
! deficit option, or at all model time steps for the schedule option.
! Irrigation amount is determined once per irrigation cycle, at the check time
! for the dynamic option, or can be prescribed. Crop specific water requirement
! is based on the actual root zone soil moisture avaiability in Sprinkler and 
! Drip, and the amount to bring the surface/column soil to saturation in Flood 
! and Paddy. Subsequent conversion from amount (mm) to rate (mm/s)
! is controlled by the irrigation duration (e.g. between 6-10 am LST), varied 
! by irrigation methods. 
! 
! Irrigation rate is scaled to grid total crop fraction when irrigation 
! fraction is less than the crop fraction.  Optionally, 
! efficiency correction is applied to account for field loss that can vary by
! application method (e.g. sprinkler 75%, drip 85%). "efcor" is defined as 
! water loss through inefficiency in the system where crop_water_deficit is 
! multiplied by 100/(100-efcor). For example, efcor is 25% for the sprinkler 
! system with 75% efficiency.
! Aligning the irrigation fraction with the crop fraction is done in LDT for 
! crop tile option, and for no-tiling or no-crop tiling options, a scale is 
! predetermined in the alltypes_irrigationMod module (i.e. IM%IrrigScale).  
! So, IrrigScale should be 1 for a crop tile if irrgation fraction is equal
! to the crop tile fraction, but in cases it is not, it is a part of 
! soil moisture availability formula and irrigation
! application for flood and drip in the model routines.  
!
! The irrigation rate is passed back to respective LSM interface routine and 
! added directly to soil moisture fields (drip, flood, and paddy) or 
! added to precipitation (sprinkler) in the public update routine.
!
! Note, the irrigation amount is a diagnostic field that needs to be subtracted
! from the precipitation for a modified water budget when sources of irrigation 
! water is not considered (i.e. P + IRR = E + R + dS).
!
! Schedule turns on at a frequency specified in the config file (eg. every 2.5 days). 
! It keeps going even when it rains unless the shutoff is triggerred--
! (using 90% RZSM for now, should test a)root zone fully saturated,
! b) surface layer fully saturated). Sprinkler schedule keeps irrigating over the days to 
! mimic revolution time of pivot system. On the other hand, Flood and Drip schedule is
! the frequency of irrigation event (eg. once 10 days) and irrigated over the duration hours.
! 
! Flood assumes furrow field where every other ditch is flooded, so the amount of irrigation
! is reduced by half.
!
! REVISION HISTORY:
!
! Dec 2021: Sarith Mahanama; Initial implementation, created a generic module 
!                            structure separated from model dependent routine 
!                            (eg. noah33_getirrigstates).
! Feb 2022: Hiroko Beaudoing; Modified subroutines.
! Oct 2022: Hiroko Beaudoing; Added schedule option for drip and flood. Flood has
!                             a special option for prescribed rate per soil type.
! Feb 2023: Hiroko Beaudoing; Fixed schedule option for flood and drip to be frequency of 
!                             irrigation event (eg once every 10-days). Removed the
!                             special option for prescribed rate per soil type.
!                             Flood assumes furrow field, thus rate is cut by a half.
!                             Pass TileNo to trigger and rate subtourines for easier
!                             debugging.
! May 2023: Hiroko Beaudoing; Fixed a bug in the schedule option for sprinkler

!EOP

MODULE IRRIGATION_MODULE

  use ESMF
  use LIS_coreMod
  use LIS_irrigationMod
  use LIS_logMod
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_timeMgrMod
  
  IMPLICIT NONE

  PRIVATE

  type, public :: irrig_state
     
     real,  pointer :: irrigRate(:)
     real,  pointer :: irrigAppRate(:)
     real,  pointer :: irrigFrac(:)
     real,  pointer :: irrigType(:)
     real,  pointer :: irrigScale(:)
     real,  pointer :: irrigRootDepth(:)
     real,  pointer :: irriggwratio(:)
     real*8,  pointer :: irrig_schedule_start(:)
     real,  pointer :: irrig_schedule_timer(:)
     
  end type irrig_state
     
  type, public, extends (irrig_state) :: irrigation_model
     
   contains
     
     ! public
     procedure, public :: get_irrigstate
     procedure, public :: update_irrigrate
     procedure, public :: get_irrig_vegindex
          
  end type irrigation_model

contains

  SUBROUTINE get_irrigstate (IM,nest,irrigState)

    implicit none
    
    class (irrigation_model), intent(inout) :: IM
    integer, intent(in)                     :: nest
    type(ESMF_State)                        :: irrigState
    type(ESMF_Field)                        :: irrigRateField,irrigFracField,irrigTypeField,prcpField
    type(ESMF_Field)                        :: irrigRootDepthField,irrigScaleField,irriggwratioField
    type(ESMF_Field)                        :: irrigAppRateField
    type(ESMF_Field)                        :: irrigScheduleTimerField
    type(ESMF_Field)                        :: irrigScheduleStartField
    integer                                 :: rc
    
    call ESMF_StateGet(irrigState, "Irrigation rate",irrigRateField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation rate')    
    call ESMF_FieldGet(irrigRateField, localDE=0,farrayPtr=IM%irrigRate,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation rate')
    
    call ESMF_StateGet(irrigState, "Applied Irrigation rate",irrigAppRateField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Applied Irrigation rate')    
    call ESMF_FieldGet(irrigAppRateField, localDE=0,farrayPtr=IM%irrigAppRate,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Applied Irrigation rate')
    
    call ESMF_StateGet(irrigState, "Irrigation frac",&
         irrigFracField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation frac')    
    call ESMF_FieldGet(irrigFracField, localDE=0,&
         farrayPtr=IM%irrigFrac,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation frac')
    
    call ESMF_StateGet(irrigState, "Irrigation max root depth",&
         irrigRootDepthField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation max root depth')    
    call ESMF_FieldGet(irrigRootDepthField, localDE=0,&
         farrayPtr=IM%irrigRootDepth,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation root depth')
    
    call ESMF_StateGet(irrigState, "Irrigation scale",&
         irrigScaleField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation scale')    
    call ESMF_FieldGet(irrigScaleField, localDE=0,&
         farrayPtr=IM%irrigScale,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation scale')  

    call ESMF_StateGet(irrigState, "Irrigation type",&
         irrigTypeField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation type')    
    call ESMF_FieldGet(irrigTypeField, localDE=0,&
         farrayPtr=IM%irrigType,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation type')

    if (LIS_irrig_struc(nest)%irrigation_GWabstraction == 1) then
       call ESMF_StateGet(irrigState, "Groundwater irrigation ratio",&
            irriggwratioField,rc=rc)
       call LIS_verify(rc,'ESMF_StateGet failed for Groundwater irrigation ratio')
       call ESMF_FieldGet(irriggwratioField, localDE=0,&
            farrayPtr=IM%irriggwratio,rc=rc)
       call LIS_verify(rc,'ESMF_FieldGet failed for Groundwater irrigation ratio')
    endif
    
    call ESMF_StateGet(irrigState, "Irrigation schedule start time",&
         irrigScheduleStartField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrig schedule start time')    
    call ESMF_FieldGet(irrigScheduleStartField, localDE=0,&
         farrayPtr=IM%irrig_schedule_start,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrig schedule start time')

    call ESMF_StateGet(irrigState, "Irrigation schedule timer",&
         irrigScheduleTimerField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrig schedule timer')    
    call ESMF_FieldGet(irrigScheduleTimerField, localDE=0,&
         farrayPtr=IM%irrig_schedule_timer,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrig schedule timer')

  END SUBROUTINE get_irrigstate
  
  ! ----------------------------------------------------------------------------

  SUBROUTINE update_irrigrate (IM, nest, TileNo, croptype, longitude, veg_trigger, &
      veg_thresh, SMWP, SMSAT, SMREF, SMCNT, RDPTH, soiltype)
    
    ! INPUTS:
    ! -------
    ! NEST         : LIS grid nest identifier
    ! TileNo       : Tile ID
    ! croptype     : MIRCA1-26, CROPMAP1-19, etc
    ! LONGITUDE    : Tile longitude
    ! VEG_TRIGGER  : Current vegetation trigger value
    ! VEG_THRESH   : vegetation threshold to turn the trigger on
    ! SMWP         : soil moisture content at wilting point [m^3/m^3]
    ! SMSAT        : soil moisture content at saturation [m^3/m^3]
    ! SMREF        : ~soil field capacity - the upper limit of water content 
    !                that soil can hold for plants [m^3/m^3]
    ! SMCNT(layers): soil moisture content in soil layers where plant roots 
    !                are active [m^3/m^3]
    ! RDPTH(layers): thicknesses of active soil layers [m]
    ! soiltype     : soil texture class index

    ! OUTPUT (internal ESMF states updated)
    ! ------
    ! irrigRate       : irrigation rate  [kg m-2 s-1]
    ! irrigAppRate    : irrigation rate applied to soil and rain [kg m-2 s-1]
    ! irrigStartTime  : start time of a scheduled irrigation application 
    !                   cycle [real*8 LIS time unit]
    ! irrigTimer      : scheduled irrigation application timer 
    !                   [seconds since the start, -1 nonactive]
    
    implicit none
    class (irrigation_model), intent(inout) :: IM
    integer, intent (in)                    :: nest, TileNo, croptype, soiltype
    real, intent (in)                       :: longitude, veg_trigger, veg_thresh, SMWP, SMSAT, SMREF, SMCNT(:), RDPTH(:)
 
    ! locals
    real     :: HC, T1, T2, ma, asmc, tsmcwlt,tsmcref
    integer  :: layer, season
    logical  :: season_active
    logical  :: irrigOn, irrigStart
    real*8   :: curtime

    asmc    = 0.0
    tsmcwlt = 0.0
    tsmcref = 0.0
    ma      = 0.0

    !-------------------------------------------------------------
    !     Compute the root zone accumlative soil moisture [mm], 
    !     field capacity [mm], and wilting point [mm] 
    !-------------------------------------------------------------
    
    SOIL_LAYERS :do layer = 1, SIZE (RDPTH)
       asmc = asmc + SMCNT(layer)*RDPTH(layer)*1000.
       tsmcwlt = tsmcwlt + smwp * rdpth(layer)*1000.0
       tsmcref = tsmcref + smref* rdpth(layer)*1000.0
    end do SOIL_LAYERS
    
    !---------------------------------------------------------------
    !     Get the root zone moisture availability to the plant
    !---------------------------------------------------------------
    
    if(tsmcref > tsmcwlt)then
       !ma = (asmc - tsmcwlt) /(tsmcref - tsmcwlt) / IM%IrrigScale(TileNo)
       !HKB: test
       ma = (asmc - tsmcwlt) /(tsmcref - tsmcwlt)
    else
       ma = -1.     
       ! HKB: keep running with warning, do not stop
       if ( LIS_rc%tscount(nest) == 1) then
       write(LIS_logunit,* )'WARN: negative ma should not happen!! ', &
        TileNo,tsmcref,tsmcwlt,IM%IrrigScale(TileNo),size(rdpth),rdpth(:)
       !call LIS_endrun()
       endif
    endif
    LIS_irrig_struc(nest)%irrigma(TileNo) = ma
    if(croptype == LIS_rc%ricecrop .and. &
                 IM%irrigType(TileNo) == 3 ) then
          LIS_irrig_struc(nest)%paddyf(TileNo) = 1.0   ! save paddy flag
    else
          LIS_irrig_struc(nest)%paddyf(TileNo) = 0.0   ! save paddy flag
    endif

    ! --------
    ! Set time
    ! --------
    
    HC = real(LIS_rc%hr) + real(LIS_rc%mn)/60. + real(LIS_rc%ss)/3600. + &
         12.* longitude/180.
    IF (HC >= 24.) HC = HC - 24.
    IF (HC <   0.) HC = HC + 24.
    T1 = CEILING (HC)     - real(LIS_rc%ts)/3600.
    T2 = FLOOR   (HC + 1) + real(LIS_rc%ts)/3600.
    
    if((HC >= T1).and.(HC < T2))then
       HC = real(NINT(HC))
    end if
    
    season_active = .false.
           
    CALENDAR: if ( LIS_irrig_struc(nest)%cropcalendar .eq. "none" ) then
              
       ! -----------------------------
       ! Vegetation Trigger
       ! -----------------------------
              
       if(veg_trigger >= veg_thresh) season_active = .true.

    else
              
       ! -----------------------------
       ! Crop calendar
       ! -----------------------------
       
       NOF_SEASONS: do season = 1, SIZE(LIS_irrig_struc(nest)%plantDay(TileNo,:))
          IF(IS_WITHIN_SEASON(LIS_rc%doy,NINT(LIS_irrig_struc(nest)%PLANTDAY(TileNo, season)), &
               NINT(LIS_irrig_struc(nest)%harvestDay(TileNo, season)))) season_active = .true.
          ! HKB: drain rice field 10 days before harvest (assume harvestDay > plantDay)
          IF(LIS_irrig_struc(nest)%paddyf(TileNo) .eq. 1.0 ) THEN
             IF(NINT(LIS_irrig_struc(nest)%PLANTDAY(TileNo, season)) .ne. -9999 .and. &
                NINT(LIS_irrig_struc(nest)%harvestDay(TileNo, season)) .ne. -9999 ) THEN
                IF(LIS_rc%doy >= 1 .and.  &
                   LIS_rc%doy > (NINT(LIS_irrig_struc(nest)%harvestDay(TileNo, season))-10) ) &
                   season_active = .false.
             ENDIF
          ENDIF
       END do NOF_SEASONS
       
    endif CALENDAR
    
    ! ----------------------------------------------------------
    ! Run irrigation model if the crop growing season is active
    ! ----------------------------------------------------------
    
    CROP_GROWING_SEASON: if (season_active) then

       curtime = LIS_rc%time   

       PADDY: if(croptype == LIS_rc%ricecrop .and. &
                 IM%irrigType(TileNo) == 3 ) then
          ! PADDY--continually add water to keep the surface soil layer at 
          ! saturation, which would maximize soil evaporation and recharge, 
          ! as would be the case with ponded water during the growing season
          IM%irrigRate(TileNo) = (smsat - SMCNT(1))*rdpth(1)*1000.0 * &
                        (100.0/(100.0-LIS_irrig_struc(nest)%flood_efcor))
          ! convert from mm to mm/s
          IM%irrigRate(TileNo) = IM%irrigRate(TileNo)*IM%IrrigScale(TileNo)/LIS_rc%ts

       else  
          ! SPRINKLER, DRIP, or FLOOD

          call irrig_trigger ( nest, TileNo, HC, ma, IM%irrigType(TileNo),  &
                               irrigOn, irrigStart, &
                               curtime, IM%irrig_schedule_start(TileNo), &
                               IM%irrig_schedule_timer(TileNo) )

          if ( irrigStart ) then
           call get_irrig_rate ( nest, TileNo, HC, ma, SMREF, SMSAT, SMCNT, &
                                 RDPTH, soiltype, &
                                 IM%irrigScale(TileNo),IM%irrigType(TileNo), &
                                 IM%irrigRate(TileNo) )
          endif

          if ( irrigOn ) then  
             ! during irrigation hours
             ! keep IM%irrigRate(TileNo) unchanged unless flood
             if ( IM%irrigType(TileNo) == 3 ) then
               ! make sure it starts at irrig_start time
               if ( IM%irrigRate(TileNo).gt.0 ) then
               call get_irrig_rate ( nest, TileNo, HC, ma, SMREF, SMSAT, &
                                 SMCNT, RDPTH, soiltype, &
                                 IM%irrigScale(TileNo),IM%irrigType(TileNo), &
                                 IM%irrigRate(TileNo) )
               endif 
             endif
          else
             ! outside of irrigation hour or shut off triggerred
             IM%irrigRate(TileNo) = 0.
          endif

       endif PADDY

    else
       !  Outside the season
       IM%irrigRate(TileNo) = 0.
       IM%irrig_schedule_start(TileNo) = -1.0
       IM%irrig_schedule_timer(TileNo) = -1.0
       
    endif CROP_GROWING_SEASON

    
  END SUBROUTINE update_irrigrate
  
  ! ----------------------------------------------------------------------------

  REAL FUNCTION crop_water_deficit (SMCNT, RDPTH, SMREF)

    implicit none

    real, intent (in)                    :: SMCNT(:), RDPTH(:), SMREF
    real                                 :: twater
    integer                              :: layer

    !---------------------------------------------------------------
    !     Get the moisture availability to the plant
    !---------------------------------------------------------------    

    crop_water_deficit = 0.
    
    do layer = 1,SIZE (SMCNT)
       crop_water_deficit = crop_water_deficit + (smref -smcnt(layer))*rdpth(layer)*1000.0
    enddo
    
    
  END FUNCTION crop_water_deficit
  
  ! ----------------------------------------------------------------------------

  logical FUNCTION IS_WITHIN_SEASON (DOY,DP, DH)

    implicit none

    integer, intent (in)                 :: DOY,DP, DH

    IS_WITHIN_SEASON = .false.
    if(DH > DP) then
       if((DOY >= DP).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    elseif (DH < DP) then
       if((DOY >= DP).AND.(DOY <= 366)) IS_WITHIN_SEASON = .true.
       if((DOY >=  1).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    endif
      
  end FUNCTION IS_WITHIN_SEASON
  ! ----------------------------------------------------------------------------
  
  SUBROUTINE irrig_trigger (nest, tileNo, HC, ma, irrigType, irrigOn, irrigStart, &
                            curtime, irrigStartTime, timer)
! Checks the irrigation trigger by soil moisture deficit or schedule option
! Returns on and off switch and indicator for time to update the irrigation rate
! at H1 for the given irrigation Type
! Updates time information for schedule irrigation 

    implicit none

    INTEGER, intent (in)                    :: nest, tileNo
    REAL, intent (in)                       :: HC, ma, irrigType
    LOGICAL, intent (out)                   :: irrigOn, irrigStart
    REAL*8, intent(in)                      :: curtime
    REAL*8, intent(inout)                   :: irrigStartTime
    REAL, intent(inout)                     :: timer
    REAL                                    :: H1, H2
    REAL                                    :: IT
    REAL*8                                  :: time2
    type(ESMF_Time)                         :: irrigTime
    LOGICAL                                 :: sprinklerOn, dripOn, floodOn
    INTEGER                                 :: doy,yr,mo,da,hr,mn,ss
    REAL                                    :: gmt
    REAL                                    :: sprinklerFreq
    REAL                                    :: dripFreq
    REAL                                    :: floodFreq
    REAL                                    :: duration
    INTEGER                                 :: irrDays

    irrigOn = .false.

    SPRINKLER: if( irrigType == 1. ) then
      IT = LIS_irrig_struc(nest)%sprinkler_thresh

      SM_DEFICIT_OR_SCHEDULE1: if ( LIS_irrig_struc(nest)%sprinkler_schedule == 0 ) then  
         ! SM_DEFICIT-- valid only during H1 <= HC < H2
         H1 = LIS_irrig_struc(nest)%sprinkler_start
         H2 = LIS_irrig_struc(nest)%sprinkler_start + LIS_irrig_struc(nest)%sprinkler_duration
         if ((HC >= H1).AND.(HC < H2)) then
            ! check soil moisture availability at H1 daily
            ! and keep it on for the duration (H1 <= HC < H2).
            if((ma <= IT).AND.(H1 == HC)) then
             sprinklerOn = .true.
             irrigStart = .true.
            else
             sprinklerOn = .true.
             irrigStart = .false.
            endif
         else
            sprinklerOn = .false.
            irrigStart = .false.
         endif
      else  
         ! SCHEDULE-- can trigger any time
         sprinklerFreq = LIS_irrig_struc(nest)%sprinkler_frequency * 86400.0
         ! set ending irrigation hour: irrigStart + duration in LIS time unites
         if ( irrigStartTime == -1 ) then   ! start a new cycle
            yr=LIS_rc%yr    !now
            mo=LIS_rc%mo
            da=LIS_rc%da
            hr=LIS_rc%hr
            mn=LIS_rc%mn
            ss=0
            call LIS_tick(time2,doy,gmt,yr,mo,da,hr,mn,ss,sprinklerFreq)
         endif
         ! check rootzone soil moisture each time
         if( ma <= IT ) then
            ! is it first DOY irrigation or after shutoff?
            if ( timer == -1 ) then   ! yes
              sprinklerOn = .true.
              irrigStart = .true.
              timer = LIS_rc%ts
              irrigStartTime = curtime
            else                      ! no
              irrigStart = .false.
              ! is it during the revolution?
              if ( timer <= sprinklerFreq ) then
                timer = timer + LIS_rc%ts
                  sprinklerOn = .true.
              else
                timer = -1
                irrigStartTime = -1
                sprinklerOn = .false.
              endif
            endif   ! timer
         !check for saturation_shutoff:90% reference RZSM
         ! HKB: NEED TO TEST OTHER SHUTOFF CONDITIONS!
         elseif ( ma > 0.90 ) then  
            sprinklerOn = .false.
            irrigStart = .false.
            irrigStartTime = -1
            timer = -1
         else  ! 0.5 < ma <= 0.9
            irrigStart = .false.
            ! is it during the revolution?
            ! ensure the timer hasn't started
            if ( timer <= sprinklerFreq .and. timer .gt. -1.0 ) then
              timer = timer + LIS_rc%ts
              ! keep irrigStartTime unchanged
                sprinklerOn = .true.
            else
              timer = -1
              irrigStartTime = -1
              sprinklerOn = .false.
            endif
         endif  ! ma
       
      endif SM_DEFICIT_OR_SCHEDULE1
      irrigOn = sprinklerOn
    endif SPRINKLER

    DRIP: if( irrigType == 2. ) then
      IT = LIS_irrig_struc(nest)%drip_thresh

      SM_DEFICIT_OR_SCHEDULE2: if ( LIS_irrig_struc(nest)%drip_schedule == 0 ) then  
         ! SM_DEFICIT-- 
         H1 = LIS_irrig_struc(nest)%drip_start
         H2 = LIS_irrig_struc(nest)%drip_start + LIS_irrig_struc(nest)%drip_duration

         if ((HC >= H1).AND.(HC < H2)) then
            ! check soil moisture availability at H1 daily and 
            ! irrigate until root zone becomes field capacity or 
            ! for the duration (H1 <= HC < H2).
           if( ma <= IT ) then
             dripOn = .true.
             if (H1 == HC) then
               irrigStart = .true.
             else
               irrigStart = .false.
             endif
!           elseif( ma > 0.99 ) then ! at field capacity, shut off
           elseif( ma > 0.90 ) then ! at field capacity, shut off
             dripOn = .false.
             irrigStart = .false.
           else
             dripOn = .true.
             irrigStart = .false.
           endif
         else
           dripOn = .false.
           irrigStart = .false.
         endif
      else  
         ! SCHEDULE-- can trigger any time but applied only during the duration
         ! drip schedule is frequency of irrigation event (eg every 3-day)
         dripFreq = LIS_irrig_struc(nest)%drip_frequency * 86400.0
         if ( irrigStartTime == -1 ) then   ! start a new cycle
            yr=LIS_rc%yr    !now
            mo=LIS_rc%mo
            da=LIS_rc%da
            hr=LIS_rc%hr
            mn=LIS_rc%mn
            ss=0
            duration = LIS_irrig_struc(nest)%drip_duration * 3600.
         else
         ! set ending irrigation hour: irrigStart + duration in LIS time unites
            call LIS_time2date(irrigStartTime,doy,gmt,yr,mo,da,hr,mn)
            irrDays = int(timer/86400.)
            ! duration in sec
            duration = LIS_irrig_struc(nest)%drip_duration * 3600. + &
                       irrDays * 86400.
            ss=0
         endif
         call LIS_tick(time2,doy,gmt,yr,mo,da,hr,mn,ss,duration)
         ! check rootzone soil moisture each time
         if( ma <= IT ) then
            ! is it first DOY irrigation or after shutoff?
            if ( timer == -1 ) then   ! yes
              dripOn = .true.
              irrigStart = .true.
              timer = LIS_rc%ts
              irrigStartTime = curtime
            else                      ! no
              irrigStart = .false.
              ! is it during the irrigation event frequency?
              if ( timer < dripFreq ) then
                timer = timer + LIS_rc%ts
                ! irrigate only during the hours under scheduled dates
                if ( curtime <= time2 ) then
                  dripOn = .true.
                else
                  dripOn = .false.
                endif
              else
                timer = -1
                irrigStartTime = -1
                dripOn = .false.
              endif
            endif   ! timer
         !check for saturation_shutoff:90% reference RZSM
         ! HKB: NEED TO TEST OTHER SHUTOFF CONDITIONS!
         elseif ( ma > 0.90 ) then
            dripOn = .false.
            irrigStart = .false.
            irrigStartTime = -1
            timer = -1
         else  ! 0.5 < ma <= 0.9
            irrigStart = .false.
            ! is it during the irrigation event frequency?
            ! ensure the timer hasn't started
            if ( timer < dripFreq .and. timer > -1.0 ) then
              timer = timer + LIS_rc%ts
              ! irrigate only during the hours under scheduled dates
              if ( curtime <= time2 ) then
                dripOn = .true.
              else
                dripOn = .false.
              endif
            else
              timer = -1
              irrigStartTime = -1
              dripOn = .false.
            endif
         endif  ! ma
      endif SM_DEFICIT_OR_SCHEDULE2
      irrigOn = dripOn
    endif DRIP

    FLOOD: if( irrigType == 3. ) then
      IT = LIS_irrig_struc(nest)%flood_thresh

      SM_DEFICIT_OR_SCHEDULE3: if ( LIS_irrig_struc(nest)%flood_schedule == 0 ) then  
         ! SM_DEFICIT-- 
         H1 = LIS_irrig_struc(nest)%flood_start
         H2 = LIS_irrig_struc(nest)%flood_start + LIS_irrig_struc(nest)%flood_duration

         if ((HC >= H1).AND.(HC < H2)) then
            ! check soil moisture availability at H1 daily
            ! and compute rate each time to keep soil saturated for 
            ! the duration (H1 <= HC < H2).
          if((ma <= IT).AND.(H1 == HC)) then
             floodOn = .true.
             irrigStart = .true.
          else
             floodOn = .true.
             irrigStart = .false.
          endif
         else
          floodOn = .false.
          irrigStart = .false.
         endif
      else  
         ! SCHEDULE-- can trigger any time but applied only during the duration
         ! flood schedule is frequency of irrigation event (eg every 10-day)
         floodFreq = LIS_irrig_struc(nest)%flood_frequency * 86400.0
         if ( irrigStartTime == -1 ) then   ! start a new cycle
            yr=LIS_rc%yr    !now
            mo=LIS_rc%mo
            da=LIS_rc%da
            hr=LIS_rc%hr
            mn=LIS_rc%mn
            ss=0
            duration = LIS_irrig_struc(nest)%flood_duration * 3600. 
         else
         ! set ending irrigation time: irrigStart + duration in LIS time unites
            call LIS_time2date(irrigStartTime,doy,gmt,yr,mo,da,hr,mn)
            irrDays = int(timer/86400.)
            ! duration in sec
            duration = LIS_irrig_struc(nest)%flood_duration * 3600. + &
                       irrDays * 86400.
            ss=0
         endif
         call LIS_tick(time2,doy,gmt,yr,mo,da,hr,mn,ss,duration)
         ! check rootzone soil moisture each time
         if( ma <= IT ) then
            ! is it first DOY irrigation or after shutoff?
            if ( timer == -1 ) then   ! yes
              floodOn = .true.
              irrigStart = .true.
              timer = LIS_rc%ts
              irrigStartTime = curtime
            else                      ! no
              irrigStart = .false.
              ! is it during the irrigation event frequency?
              if ( timer < floodFreq ) then
                timer = timer + LIS_rc%ts
                ! irrigate only during the hours under scheduled dates
                if ( curtime <= time2 ) then
                  floodOn = .true.
                else
                  floodOn = .false.
                endif
              else
                timer = -1
                irrigStartTime = -1
                floodOn = .false.
              endif
            endif   ! timer
         !check for saturation_shutoff:90% reference RZSM
         ! HKB: NEED TO TEST OTHER SHUTOFF CONDITIONS!
         elseif ( ma > 0.90 ) then
            floodOn = .false.
            irrigStart = .false.
            irrigStartTime = -1
            timer = -1
         else  ! 0.5 < ma <= 0.9
            irrigStart = .false.
            ! is it during the irrigation event frequency?
            ! ensure the timer hasn't started
            if ( timer < floodFreq .and. timer > -1.0 ) then
              timer = timer + LIS_rc%ts
              ! irrigate only during the hours under scheduled dates
              if ( curtime <= time2 ) then
                floodOn = .true.
              else
                floodOn = .false.
              endif
            else
              timer = -1
              irrigStartTime = -1
              floodOn = .false.
            endif
         endif  ! ma
      endif SM_DEFICIT_OR_SCHEDULE3
      irrigOn = floodOn
    endif FLOOD

  END SUBROUTINE irrig_trigger
  ! ----------------------------------------------------------------------------

  SUBROUTINE get_irrig_rate (nest, tileNo, HC, ma, SMREF, SMSAT, SMCNT, RDPTH, &
                             soiltype, IrrigScale, irrigType, IRATE)
! This subroutine computes irrigation rate for Sprinkler, Drip, or Flood.
! The irrigation rate can be computed dunamically based on the water deficit,
! or prescribed (NEW!). This routine is called only when irrig_trigger routine
! determined it is the time to compute the rate via irrigStart flag.
! The irrigation amount is computed at H1 (start of 
! irrigation, eg. 6am local time) and kept during the irrigation hours for
! Sprinkler and Drip.  For flood, irrigation amount is checked and determined
! each time step to keep surface soil layer saturated for the fixed duration.
! The rate is zero outside of the irrigation hours.

    implicit none

    INTEGER, intent (in)                    :: nest, soiltype, tileNo
    REAL, intent (in)                       :: HC, ma, SMREF, SMSAT, SMCNT(:), RDPTH(:)
    REAL, intent (in)                       :: IrrigScale, irrigType
    REAL                                    :: H1, H2, IT
    REAL, intent (out)                      :: IRATE
    REAL                                    :: SRATE, DRATE, FRATE
    REAL                                    :: maxinfrate

   SPRINKLER: if( irrigType == 1. ) then
    ! The rate and duration are fixed
       H1 = LIS_irrig_struc(nest)%sprinkler_start
       H2 = LIS_irrig_struc(nest)%sprinkler_start + LIS_irrig_struc(nest)%sprinkler_duration
       IT = LIS_irrig_struc(nest)%sprinkler_thresh

    DYNAMIC_OR_PRESCRIBED1: if( LIS_irrig_struc(nest)%sprinkler_rate .gt. 0 ) then
       ! PRESCRIBED
     SCHEDULE1P: if ( LIS_irrig_struc(nest)%sprinkler_schedule .gt. 0 ) then
       ! ON SCHEDULE, get rate anytime
               SRATE = LIS_irrig_struc(nest)%sprinkler_rate*(100.0/(100.0-LIS_irrig_struc(nest)%sprinkler_efcor))*IrrigScale/3600.
     else   
       ! not on SCHEDULE, compute rate only during irrigation hours
       if ((HC >= H1).AND.(HC < H2)) then
          ! get the prescribed rate in mm/hr and convert to mm/s at H1 to compute the
          ! rate for the day and maintain the same rate through out the irrigation
          ! duration (H1 <= HC < H2).
          if(H1 == HC) &
               SRATE = LIS_irrig_struc(nest)%sprinkler_rate*(100.0/(100.0-LIS_irrig_struc(nest)%sprinkler_efcor))*IrrigScale/3600.
       else
          SRATE = 0.
       endif
     endif SCHEDULE1P
    else  ! DYNAMIC

     SCHEDULE1D: if ( LIS_irrig_struc(nest)%sprinkler_schedule .gt. 0 ) then
       ! ON SCHEDULE, get rate anytime
       ! Note: if sprinkler is on schedule, the rate is normally prescribed
       ! however, the option to compute the demand based rate is provided
       ! below, but needs refinement (i.e. rate is for H1/H2 duration).
               SRATE = crop_water_deficit (SMCNT, RDPTH, SMREF)*(100.0/(100.0-LIS_irrig_struc(nest)%sprinkler_efcor))*IrrigScale/(H2 - H1)/3600.

     else   
       ! not on SCHEDULE, compute rate only during irrigation hours
       if ((HC >= H1).AND.(HC < H2)) then
          ! The model uses rootzone soil moisture state at H1 to compute irrigation
          ! rates for the day and maintains the same rate through out the irrigation
          ! duration (H1 <= HC < H2).
          if((ma <= IT).AND.(H1 == HC)) &
               SRATE = crop_water_deficit (SMCNT, RDPTH, SMREF)*(100.0/(100.0-LIS_irrig_struc(nest)%sprinkler_efcor))*IrrigScale/(H2 - H1)/3600.
       else
          SRATE = 0.
       endif
     endif SCHEDULE1D
    endif DYNAMIC_OR_PRESCRIBED1
    IRATE = SRATE
   endif SPRINKLER

   DRIP: if( irrigType == 2. ) then
    ! The rate is fixed but duration is variable.
    H1 = LIS_irrig_struc(nest)%drip_start
    H2 = LIS_irrig_struc(nest)%drip_start + LIS_irrig_struc(nest)%drip_duration
    IT = LIS_irrig_struc(nest)%drip_thresh

    DYNAMIC_OR_PRESCRIBED2: if( LIS_irrig_struc(nest)%drip_rate .gt. 0 ) then
       ! PRESCRIBED
     SCHEDULE2P: if ( LIS_irrig_struc(nest)%drip_schedule .gt. 0 ) then
       ! ON SCHEDULE, get rate anytime
       DRATE = LIS_irrig_struc(nest)%drip_rate * &
               (100.0/(100.0-LIS_irrig_struc(nest)%drip_efcor))*IrrigScale/3600.
     else   
       ! not on SCHEDULE, compute rate only during irrigation hours
       if ((HC >= H1).AND.(ma < 1.0)) then
          ! get the prescribed rate in mm/hr and convert to mm/s at H1 to compute the
          ! rate for the day and maintain the same rate through out the irrigation
          ! until rootzone becomes field capacity
          ! note to myself: when applying amount to surface layer, make sure it doesn't exceed SMSAT?
          if(H1 == HC) &
               DRATE = LIS_irrig_struc(nest)%drip_rate * &
                      (100.0/(100.0-LIS_irrig_struc(nest)%drip_efcor))*IrrigScale/3600.
       else  
          DRATE = 0.
       endif
     endif SCHEDULE2P
    else ! DYNAMIC

     SCHEDULE2D: if ( LIS_irrig_struc(nest)%drip_schedule .gt. 0 ) then
       ! ON SCHEDULE, get rate anytime
       ! Note: if drip is on schedule, the rate is normally prescribed
       ! however, the option to compute the demand based rate is provided
       ! below, but needs testing.
       DRATE = crop_water_deficit (SMCNT, RDPTH, SMREF) *  &
               (100.0/(100.0-LIS_irrig_struc(nest)%drip_efcor))*IrrigScale/(H2 - H1)/3600.
     else   
       ! not on SCHEDULE, compute rate only during irrigation hours
       if ((HC >= H1).AND.(ma < 1.0)) then
          ! Check root zone soil moisture at H1 and irrigate until root zone becomes 
          ! field capacity (ma = 1)
          ! To make drip rate lower than sprinkler with deficit option, irrigation duration 
          ! may need to be increased.
          ! Alternatively, deficit can be computed for surface layer and multiply by
          ! LIS_rc%ts like in flood.  
          if((ma <= IT).AND.(H1 == HC)) then
               DRATE = crop_water_deficit (SMCNT, RDPTH, SMREF) *  &
                       (100.0/(100.0-LIS_irrig_struc(nest)%drip_efcor))*IrrigScale/(H2 - H1)/3600.
          endif
       else
          DRATE = 0.
       endif
     endif SCHEDULE2D
    endif DYNAMIC_OR_PRESCRIBED2
    IRATE = DRATE
   endif  DRIP

   FLOOD: if( irrigType == 3. ) then
    ! The duration is fixed but rate is variable 
    H1 = LIS_irrig_struc(nest)%flood_start
    H2 = LIS_irrig_struc(nest)%flood_start + LIS_irrig_struc(nest)%flood_duration
    IT = LIS_irrig_struc(nest)%flood_thresh

    DYNAMIC_OR_PRESCRIBED3: if( LIS_irrig_struc(nest)%flood_rate .gt. 0 ) then
       ! PRESCRIBED
     SCHEDULE3P: if ( LIS_irrig_struc(nest)%flood_schedule .gt. 0 ) then
       ! ON SCHEDULE, get rate anytime
       FRATE = LIS_irrig_struc(nest)%flood_rate* &
           (100.0/(100.0-LIS_irrig_struc(nest)%flood_efcor))*IrrigScale/3600.
     else   
       ! not on SCHEDULE, compute rate only during irrigation hours
       if ((HC >= H1).AND.(HC < H2)) then
          ! get the prescribed rate in mm/hr and convert to mm/s 
          FRATE = LIS_irrig_struc(nest)%flood_rate* &
            (100.0/(100.0-LIS_irrig_struc(nest)%flood_efcor))*IrrigScale/3600.
       else  ! reset
          FRATE = 0.
       endif
     endif SCHEDULE3P
    else  ! DYNAMIC
     ! use SMSAT instead of SMREF
     ! toplayer => LIS_rc%irrigation_mxsoildpth = 1
     ! but can be expanded to entire column 
     ! note: should we check the rate is greater than infiltration rate?
     SCHEDULE3D: if ( LIS_irrig_struc(nest)%flood_schedule .gt. 0 ) then
       ! ON SCHEDULE, get rate anytime
       ! Note: if flood is on schedule, the rate is normally prescribed
       ! however, the option to compute the demand based rate is provided
       ! below, but not tested.
       FRATE = crop_water_deficit (                                 &
               SMCNT(1:LIS_irrig_struc(nest)%irrigation_mxsoildpth),&
               RDPTH(1:LIS_irrig_struc(nest)%irrigation_mxsoildpth),&
               SMSAT)*                                              &
               (100.0/(100.0- LIS_irrig_struc(nest)%flood_efcor))*  &
               IrrigScale/LIS_rc%ts
     else   
       ! not on SCHEDULE, compute rate only during irrigation hours
       if ((HC >= H1).AND.(HC < H2)) then
          ! check rootzone soil moisture and saturate top soil layer 
          ! during H1 <= HC < H2 (the rate is variable). 
          ! new: do not check for ma, for FRACE to be variable each time step 
          ! whereas normally ma is checked at irrig_start_hour 
          ! if( ma <= IT ) then
            FRATE = crop_water_deficit (                                 &
                    SMCNT(1:LIS_irrig_struc(nest)%irrigation_mxsoildpth),&
                    RDPTH(1:LIS_irrig_struc(nest)%irrigation_mxsoildpth),&
                    SMSAT)*                                              &
                    (100.0/(100.0- LIS_irrig_struc(nest)%flood_efcor))*  &
                    IrrigScale/LIS_rc%ts
          !endif
       else
          FRATE = 0.
       endif
     endif SCHEDULE3D
    endif  DYNAMIC_OR_PRESCRIBED3
    ! Assume furrow field, every other ditch is filled, so 50% of field is
    ! saturated.  Apply this assumption for both demand & prescribed/schedule.
    IRATE = FRATE*0.5
   endif  FLOOD

  END SUBROUTINE get_irrig_rate
    ! ----------------------------------------------------------------

  subroutine get_irrig_vegindex(IM, veg_index1, veg_index2, nlctypes)

 implicit none
  class (irrigation_model), intent(inout) :: IM
  integer, intent(inout)    :: veg_index1,veg_index2,nlctypes

  if(LIS_rc%lcscheme.eq."UMD") then !UMD
     veg_index1 = 6
     veg_index2 = 11
     nlctypes   = 13
  elseif(LIS_rc%lcscheme.eq."UMD+MIRCA") then !UMD+MIRCA (Temporary, KRA)
     veg_index1 = 6
     veg_index2 = 40
     nlctypes   = 14
  elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") then
     veg_index1 = 6
     veg_index2 = 14
     nlctypes   = 20
  elseif(LIS_rc%lcscheme.eq."IGBPNCEP+MIRCA") then
     veg_index1 = 6
     veg_index2 = 46
     nlctypes   = 20
  elseif(LIS_rc%lcscheme.eq."USGS") then
     veg_index1 = 2
     veg_index2 = 10
     nlctypes   = 24
  elseif(LIS_rc%lcscheme.eq."UMDCROPMAP") then
     veg_index1 = 5
     veg_index2 = 32
     nlctypes   = 13
  else
     write(LIS_logunit,*) '[ERR] The landcover scheme ',trim(LIS_rc%lcscheme)
     write(LIS_logunit,*) '[ERR] is not supported for the Noah.3.3 irrigation module.'
     call LIS_endrun()
  endif

  end subroutine get_irrig_vegindex
 
  subroutine get_soiltex_infilrate(n,soiltype,maxinfrate) 
!  Map out soiltype into four hydrologic types and assign infiltration rate 
!  Four Hydrologic soil groups A-D correspond to general soil types that
!  range in the basic infiltration rate in mm/hr:
!  A) sand and sandy loam, B) loam, C) clay loam, and D) Clay
!  Group: Infiltration rate [mm/hr], (soil texture classes)  
!  A: < 30, 20-30 (gravel, sandy gravel, silty gravels, gravelly sands, sand
!                  sand, loamy sand, sandy loam)
!  B: 10 - 20     (silty sands, loam, silt loam)
!  C: 5 - 10      (sandy clay loam, silts)
!  D: 1 - 5       (clay loam, silty clay loam, sandy clay, silty clay, clay)
!   
!  Values based on FAO Annex 2 Infiltration rate and infiltration test website.
!  Take the max value for testing the upper end of flood irrigation.
! 
!  Note: this routine is not used but keeping it for a reference
!
   implicit none
    integer, intent(in)    :: n,soiltype
    real, intent(out)      :: maxinfrate

    maxinfrate = 0.
    if (LIS_irrig_struc(n)%soiltextscheme.eq."STATSGO")then     !NOAH or NoahMP(19)
       select case (soiltype)
       case (1,2,3,13,19)  !A
        !sand,loamy sand,sandy loam,organic material,white sand
         maxinfrate = 30.
       case (4,6)  !B
        !silt loam,loam
         maxinfrate = 20.
       case (5,7)  !C
        !silt,sandy clay loam,
         maxinfrate = 10.
       case (8,9,10,11,12,17)  !D
        !slity clay loam,clay loam,sandy clay,silty clay,clay,plava
         maxinfrate = 5.
       case default  ! 14,15,16,18
        !water,bedrock,other(land-ice),lava
         write(LIS_logunit,*) '[ERR] infiltration--Should not have STATSGO type:',&
                              LIS_irrig_struc(n)%soiltextscheme 
       end select 
    elseif (LIS_irrig_struc(n)%soiltextscheme.eq."Zobler") then !NOAH or NoahMP
       write(LIS_logunit,*) '[ERR] Noah Zobler soil tex mapping not yet implemented' 
    elseif (LIS_irrig_struc(n)%soiltextscheme.eq."Soil texture not selected") then
     if (LIS_rc%lsm.eq."CLSM F2.5") then
       write(LIS_logunit,*) '[ERR] CLSM soil tex mapping not yet implemented' 
     else  ! VIC
       write(LIS_logunit,*) '[ERR] VIC soil tex mapping not yet implemented' 
     endif
    else
     write(LIS_logunit,*) '[ERR] The soil texture scheme ',trim(LIS_irrig_struc(n)%soiltextscheme)
     write(LIS_logunit,*) '[ERR] is not supported in get_soiltex_infilrate'
     call LIS_endrun()
    endif

  end subroutine get_soiltex_infilrate
 
END MODULE IRRIGATION_MODULE
