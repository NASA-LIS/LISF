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
! !ROUTINE: noah33_getirrigationstates
! \label{noah33_getirrigationstates}
! 
! !INTERFACE:
subroutine noah33_getirrigationstates(n,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use noah33_lsmMod
  use LIS_irrigationMod

! !DESCRIPTION:        
!
! Calculate water requirement and apply the amount for Sprinkler, Drip, or
! Flood irrigation types accordingly. The water requirement is checked once
! a day at specified local time. Irrigation check time, duration, and threshold
! are specifield for each irrigation type in the lis.config. 
!
! Sprinker: water requirement is based on the actual root zone soil moisture
!           availability. The irrigation threshold determines when to irrigate
!           (i.e. when root zone soil moisture falls below a percentage of 
!           the fields capacity). Apply the amount to precipitation  
!           during the specified hours, eg 6-10 am LST.
! Drip: water requirement is based on the transpiration stress.  Apply the 
!       amount to the top layer soil moisture during the specified hours.
! Flood: water requirement is based on the amount to bring the entire column
!        to saturation and applied to all soil layers at once.
!
! The root zone is actual maximum root depth rather than NOAH's vegetation 
! parameter "NROOT". The root depth grows following the GVF seasonality to 
! reflect the water demand increase/decrease for crops. (note: need age of
! trees for orchard).
! 
! Irrigation amount is scaled to grid total crop fraction when intensity
! is less than the fraction.  Irrigation is expanded to non-crop, non-forest,
! non-baresoil/urban tiles if intensity exceeds grid total crop fraction.
! In the latter case, scaled irrigation is applied to grassland first, 
! then further applied over the rest of tiles equally if the intensity 
! exceeds grassland fraction as well. The irrigation scale and alighning 
! intensity with the crops or land cover types, as well as irrigation type
! determination for the tile are done before this routine is called.   
!
! Optionally efficiency correction is applied to account for field loss. 
!
! This version includes modifications and updates as follows:
! 1) Change growing season threshold (i.e. 40% of GFRAC range) in lis.config
! 2) Use Crop Plant/Harvest dates in addition to GVF for trigger
! 3) Irrigation check time, duration, threshold, and efficiency are set in 
!    lis.config and can be different for irrigation types 
!
! Note: need to add rice paddy in flood 
!       add option to use crop calender for growing season check
!
! REVISION HISTORY:
!
! Aug 2008: Hiroko Kato; Initial code
! Nov 2012: Sujay Kumar, Incorporated into LIS
! Jun 2014: Ben Zaitchik; Added flood scheme
! Feb 2020: Jessica Erlingis; Fix sprinkler irrigation winodw 
! Dec 2020: Hiroko Beaudoing; Updated things based on old LIS/Noah and 
!                             incorporated Sarith's concurrent irrigation types
!                             and Wanshu/Ben's modifications.
!                             
!EOP
  implicit none
! moved to lis.config
  ! Sprinkler parameters
!  real, parameter      :: otimess = 6.0 ! local trigger check start time [hour]
!  real, parameter      :: irrhrs = 4.   ! duration of irrigation hours 
  ! Drip parameters (not currently implemented)
!  real, parameter      :: otimeds = 6.0 ! local trigger check start time [hour]
!  real, parameter      :: irrhrd = 12.0   ! duration of irrigation hours 
 ! Flood parameters
!  real, parameter      :: otimefs = 6.0 ! local trigger check start time [hour]
!  real, parameter      :: irrhrf = 1.0   ! duration of irrigation hours 
  !!!real, parameter      :: ffreq = 0.0 ! frequency of flood irrig [days] set to 0.0 to use thresh instead
  
!  real, parameter      :: efcor = 0.0      ! Efficiency Correction (%)

  integer              :: n
  integer              :: rc
  integer              :: t,k,gid,vegt,l
  type(ESMF_State)     :: irrigState
  type(ESMF_Field)     :: irrigRateField,irrigFracField,irrigTypeField
  type(ESMF_Field)     :: irrigRootDepthField,irrigScaleField
  
  real,  pointer       :: irrigRate(:), irrigFrac(:), irrigType(:)
  real,  pointer       :: irrigRootDepth(:), irrigScale(:)
  integer              :: chhr, lhr
  real                 :: asmc, tsmcwlt, tsmcref, ma
  real                 :: sldpth(noah33_struc(n)%nslay)
  real                 :: rdpth(noah33_struc(n)%nslay)
  real                 :: zdpth(noah33_struc(n)%nslay)
  real                 :: water(noah33_struc(n)%nslay)
  real                 :: twater, twater1, twater2
  real                 :: ippix, crootd
  real                 :: smcmax, smcref, smcwlt
  real                 :: shdfac
  integer              :: lroot,veg_index1,veg_index2
  real                 :: gsthresh, ltime
  logical              :: irrig_check_frozen_soil
  real                 :: timestep, shift_otimess, shift_otimese
  real                 :: shift_otimeds,shift_otimede
  real                 :: shift_otimefs,shift_otimefe
! _______________________________________________________

  call ESMF_StateGet(irrigState, "Irrigation rate",irrigRateField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation rate')    
  call ESMF_FieldGet(irrigRateField, localDE=0,farrayPtr=irrigRate,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation rate')

  call ESMF_StateGet(irrigState, "Irrigation frac",&
       irrigFracField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation frac')    
  call ESMF_FieldGet(irrigFracField, localDE=0,&
       farrayPtr=irrigFrac,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation frac')

  call ESMF_StateGet(irrigState, "Irrigation max root depth",&
       irrigRootDepthField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation max root depth')    
  call ESMF_FieldGet(irrigRootDepthField, localDE=0,&
       farrayPtr=irrigRootDepth,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation root depth')

  call ESMF_StateGet(irrigState, "Irrigation scale",&
       irrigScaleField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation scale')    
  call ESMF_FieldGet(irrigScaleField, localDE=0,&
       farrayPtr=irrigScale,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation scale')  

  call ESMF_StateGet(irrigState, "Irrigation type",&
       irrigTypeField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation type')    
  call ESMF_FieldGet(irrigTypeField, localDE=0,&
       farrayPtr=irrigType,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation type')

!JE  irrigRate = 0.0  

!----------------------------------------------------------------------
! Set start and end times for each irrigation type
!----------------------------------------------------------------------
  timestep = LIS_rc%ts
  ! Adjust bounds by timestep to account for the fact that LIS_rc%hr, etc.
  ! will represents the END of the integration timestep window
  ! Sprinkler
  shift_otimess = LIS_irrig_struc(n)%sprinkler_start + (timestep/3600.)
  shift_otimese = (LIS_irrig_struc(n)%sprinkler_start + &
                   LIS_irrig_struc(n)%sprinkler_duration) + (timestep/3600.)
  ! Drip
  shift_otimeds = LIS_irrig_struc(n)%drip_start + (timestep/3600.)
  shift_otimede = (LIS_irrig_struc(n)%drip_start + &
                   LIS_irrig_struc(n)%drip_duration) + (timestep/3600.)
  ! Flood
  shift_otimefs = LIS_irrig_struc(n)%flood_start + (timestep/3600.)
  shift_otimefe = (LIS_irrig_struc(n)%flood_start + &
                   LIS_irrig_struc(n)%flood_duration) + (timestep/3600.)

! Set vegetation type index to be irrigated 
  if(LIS_rc%lcscheme.eq."UMD") then !UMD
      veg_index1 = 6
      veg_index2 = 11
  elseif(LIS_rc%lcscheme.eq."UMD+MIRCA") then !UMD+MIRCA (Temporary, KRA)
      veg_index1 = 6
      veg_index2 = 40
  elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") then
      veg_index1 = 6
      veg_index2 = 14
  elseif(LIS_rc%lcscheme.eq."IGBPNCEP+MIRCA") then  
      veg_index1 = 6
      veg_index2 = 46
  elseif(LIS_rc%lcscheme.eq."USGS") then 
      veg_index1 = 2
      veg_index2 = 10
  elseif(LIS_rc%lcscheme.eq."UMDCROPMAP") then 
      veg_index1 = 5
      veg_index2 = 32
  else
      write(LIS_logunit,*) '[ERR] The landcover scheme ',trim(LIS_rc%lcscheme)
      write(LIS_logunit,*) '[ERR] is not supported for the Noah.3.3 irrigation module.'
      call LIS_endrun()
  endif

! Set global parameters
  sldpth(1) = noah33_struc(n)%lyrthk(1)         ! Soil layer thicknesses (m)
  sldpth(2) = noah33_struc(n)%lyrthk(2)
  sldpth(3) = noah33_struc(n)%lyrthk(3)
  sldpth(4) = noah33_struc(n)%lyrthk(4)
  zdpth(1) = sldpth(1)         ! Soil layer depth from top (m)
  zdpth(2) = sldpth(1) + sldpth(2)
  zdpth(3) = sldpth(1) + sldpth(2) + sldpth(3)
  zdpth(4) = sldpth(1) + sldpth(2) + sldpth(3) + sldpth(4)
 
!-----------------------------------------------------
! Main tile loop
!-----------------------------------------------------
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

   irrig_check_frozen_soil = .false.
     
   if((noah33_struc(n)%noah(t)%smc(1) - &
        noah33_struc(n)%noah(t)%sh2o(1)).gt.0001) then 
      irrig_check_frozen_soil = .true. 
   elseif((noah33_struc(n)%noah(t)%smc(2) - &
        noah33_struc(n)%noah(t)%sh2o(2)).gt.0001) then 
      irrig_check_frozen_soil = .true. 
   elseif((noah33_struc(n)%noah(t)%smc(3) - &
        noah33_struc(n)%noah(t)%sh2o(3)).gt.0001) then 
      irrig_check_frozen_soil = .true. 
   elseif((noah33_struc(n)%noah(t)%smc(4) - &
        noah33_struc(n)%noah(t)%sh2o(4)).gt.0001) then 
      irrig_check_frozen_soil = .true. 
   elseif(noah33_struc(n)%noah(t)%stc(2).le.LIS_CONST_TKFRZ) then
      irrig_check_frozen_soil = .true. 
   elseif(noah33_struc(n)%noah(t)%stc(3).le.LIS_CONST_TKFRZ) then
      irrig_check_frozen_soil = .true. 
   elseif(noah33_struc(n)%noah(t)%stc(4).le.LIS_CONST_TKFRZ) then
      irrig_check_frozen_soil = .true. 
   endif

   if(.not.irrig_check_frozen_soil) then 
     twater  = 0.0
     water   = 0.0
     asmc    = 0.0
     tsmcwlt = 0.0
     tsmcref = 0.0
     ma      = 0.0
     crootd  = 0.0
     lroot   = 0

     smcmax =  noah33_struc(n)%noah(t)%smcmax
     smcref = noah33_struc(n)%noah(t)%smcref
     smcwlt = noah33_struc(n)%noah(t)%smcwlt
     shdfac =  noah33_struc(n)%noah(t)%shdfac

     ! Process only irrigated tiles
     if(irrigFrac(t).gt.0) then 
       ippix = irrigFrac(t)*0.01
           
       ! Determine the amount of irrigation to apply if irrigated tile
       if( IrrigScale(t).gt.0.0 ) then ! irrigated tile

       ! Proceed if it is non-forest, non-baresoil, non-urban
         vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
         if(vegt.ge.veg_index1.and.vegt.le.veg_index2&
             .and.vegt.ne.LIS_rc%bareclass.and.&
             vegt.ne.LIS_rc%urbanclass) then 
                 
       ! Check GVF growing season 
          gsthresh = noah33_struc(n)%noah(t)%shdmin + & 
                     0.40 * (noah33_struc(n)%noah(t)%shdmax - &
                     noah33_struc(n)%noah(t)%shdmin)
       ! HKB ==> Change to below after benchmarking
       ! let gsthresh be a function of the range, which means the larger
       ! the range is, the higher GVF threshold will be for this grid. (WN)
       !   gsthresh = noah33_struc(n)%noah(t)%shdmin + & 
       !           (LIS_rc%irrigation_GVFparam1 + LIS_rc%irrigation_GVFparam2*&
       !     (noah33_struc(n)%noah(t)%shdmax - noah33_struc(n)%noah(t)%shdmin))*&
       !     (noah33_struc(n)%noah(t)%shdmax - noah33_struc(n)%noah(t)%shdmin)
                     
                 
          if(shdfac .ge. gsthresh) then    
       ! Calculate vegetation and root depth parameters
            crootd = irrigRootdepth(t)*shdfac
            if(crootd.gt.0.and.crootd.lt.zdpth(1)) then 
              lroot = 1
              rdpth(1) = crootd
            elseif(crootd .ge. zdpth(1).and.crootd .lt. zdpth(2) ) then
              lroot = 2
              rdpth(1) = sldpth(1)
              rdpth(2) = crootd - zdpth(1)
            elseif ( crootd.ge.zdpth(2).and.crootd .lt. zdpth(3) ) then
              lroot = 3
              rdpth(1) = sldpth(1)
              rdpth(2) = sldpth(2)
              rdpth(3) = crootd - zdpth(2)
            elseif ( crootd.ge.zdpth(3).and.crootd .lt. zdpth(4) ) then
              lroot = 4
              rdpth(1) = sldpth(1)
              rdpth(2) = sldpth(2)
              rdpth(3) = sldpth(3)
              rdpth(4) = crootd - zdpth(3)
            endif
       
            gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
            chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
            if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
               (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
               chhr = chhr -1
            lhr = LIS_rc%hr +chhr
            if(lhr.ge.24) lhr = lhr-24
            if(lhr.lt.0) lhr = lhr+24
                
            ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0

   

!!!!! SPRINKLER IRRIGATION
            if(irrigType(t).eq.1) then
             ! If we are outside of the irrigation window, set rate to 0
               if ((ltime.ge.shift_otimese).or.(ltime.lt.shift_otimess)) then
                  irrigRate(t) = 0.0
               endif
!----------------------------------------------------------------------       
!    Set the irrigation rate at start time; keep the value till next day
!----------------------------------------------------------------------       
               if(ltime.eq.shift_otimess) then 
                  irrigRate(t) = 0.0
!-------------------------------------------------------------
!     Compute the root zone accumlative soil moisture [mm], 
!     field capacity [mm], and wilting point [mm] 
!-------------------------------------------------------------
                  if(lroot.gt.0) then 
                    do k=1,lroot
                      asmc = asmc + noah33_struc(n)%noah(t)%smc(k)*&
                             rdpth(k)*1000.0
                      tsmcwlt = tsmcwlt + smcwlt * rdpth(k)*1000.0
                      tsmcref = tsmcref + smcref * rdpth(k)*1000.0
                     enddo
!---------------------------------------------------------------
!     Get the root zone moisture availability to the plant
!--------------------------------------------------------------- 
                     ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)
                     if(ma.le.LIS_irrig_struc(n)%sprinkler_thresh) then 
                        do k=1,lroot
                           water(k) = &
                                  (smcref-noah33_struc(n)%noah(t)%smc(k))*&
                                        rdpth(k)*1000.0
                           twater = twater + water(k)
                        enddo
                                
!-----------------------------------------------------------------------------
!     Scale the irrigation amount for thr irrigated fraction (intensity).
!-----------------------------------------------------------------------------
                        twater1 = twater
                        twater = twater * irrigScale(t)
                                
!-----------------------------------------------------------------------------
!     Apply efficiency correction
!-----------------------------------------------------------------------------
                        twater2 = twater
                        twater = twater*(100.0/(100.0-LIS_irrig_struc(n)%sprinkler_efcor))
!-----------------------------------------------------------------------------
!     Compute irrigation rate
!-----------------------------------------------------------------------------
                        irrigRate(t) = twater/(LIS_irrig_struc(n)%sprinkler_duration*3600.0)
!HKB check                        write(LIS_logunit,*) '[INFO] Irrigating SP ', &
!                                     t,ltime, twater
                     endif   ! ma < thresh
                  else
                    write(LIS_logunit,*) '[ERR] lroot should be > 0!'
                    call LIS_endrun()
                  endif   ! lroot > 0
               endif   ! ltime = shift_otimes
!!!!! DRIP IRRIGATION (NOT CURRENTLY IMPLEMENTED)
            elseif(irrigType(t).eq.2) then
! Need to get crop coefficient so that we can caculate unstressed Transp
!       RC=RSMIN/(XLAI*RCS*RCT*RCQ)
!       PCIRR=(RR+DELTA)/(RR*(1.+RC*CH)+DELTA)
! CALL TRANSP (with PCIRR)

! Then add enough water to get from actual Transp to unstressed Transp
                       twater = 0.0
!-----------------------------------------------------------------------------
!     Apply efficiency correction
!-----------------------------------------------------------------------------
                       twater2 = twater
                       twater = twater*(100.0/(100.0-LIS_irrig_struc(n)%drip_efcor))
!-----------------------------------------------------------------------------
!     Compute irrigation rate
!-----------------------------------------------------------------------------
                       irrigRate(t) = twater  ! for drip calculation, twater is a rate [kg m-2/s]
                       noah33_struc(n)%noah(t)%smc(1) = &
                            noah33_struc(n)%noah(t)%smc(1) + (twater-twater2)/(sldpth(1)*1000.0) !! check this with Sujay
!HKB check                       write(LIS_logunit,*) '[INFO] Irrigating DP ', &
!                                     t,ltime, twater
                     
!!!!! FLOOD IRRIGATION
            elseif(irrigType(t).eq.3) then
!               if(ltime.eq.shift_otimefs) then 
               if((ltime.ge.shift_otimefs).and.(ltime.lt.shift_otimefe)) then 
!-------------------------------------------------------------
!     Compute the root zone accumlative soil moisture [mm], 
!     field capacity [mm], and wilting point [mm] 
!-------------------------------------------------------------
                 if(lroot.gt.0) then 
                   do k=1,lroot
                     asmc = asmc + noah33_struc(n)%noah(t)%smc(k)*&
                            rdpth(k)*1000.0
                     tsmcwlt = tsmcwlt + smcwlt * rdpth(k)*1000.0
                     tsmcref = tsmcref + smcref * rdpth(k)*1000.0
                   enddo
!---------------------------------------------------------------
!     Get the root zone moisture availability to the plant
!--------------------------------------------------------------- 
!                  ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)   ! Original
                   ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)/IrrigScale(t) ! BZ added IrrigScale
                        
                   if( ma .le. LIS_irrig_struc(n)%flood_thresh ) then
                     do l = 1, LIS_rc%irrigation_mxsoildpth
                        if( l == 1 ) then
                          twater = (smcmax - noah33_struc(n)%noah(t)%smc(l))*sldpth(l)*1000.0
                        else
                        ! BZ modification 4/2/2015 to saturate entire column and apply ippix 
                          twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(l))*sldpth(l)*1000.0
                        endif
                      end do
                           
!-----------------------------------------------------------------------------
!     Scale the irrigation amount for thr irrigated fraction (intensity).
!-----------------------------------------------------------------------------
                      twater1 = twater
                      twater = twater * irrigScale(t)
!-----------------------------------------------------------------------------
!     Apply efficiency correction
!-----------------------------------------------------------------------------
                      twater2 = twater
                      twater = twater*(100.0/(100.0-LIS_irrig_struc(n)%flood_efcor))
!-----------------------------------------------------------------------------
!     Compute irrigation rate
!-----------------------------------------------------------------------------
                      irrigRate(t) = twater/LIS_rc%ts
                           
! BZ modification 4/2/2015 to account for ippix and all soil layers:
                      do l = 1, LIS_rc%irrigation_mxsoildpth
                        noah33_struc(n)%noah(t)%smc(l) = IrrigScale(t)*smcmax + &
                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(l)
                      end do
!HKB check                      write(LIS_logunit,*) '[INFO] Irrigating FD ', &
!                                     t,ltime, twater
                    endif   ! ma <= thresh
                  else
                    write(LIS_logunit,*) '[ERR] lroot should be > 0!!!'
                    call LIS_endrun()
                  endif  ! lroot > 0
               end if    ! flood irrigation hours
                       
            endif  ! irrigation type
           endif  ! growing season
         end if   ! vegclass
       end if  ! irrigation scale > 0
     end if  ! irrigated tile
   endif  ! not frozen soil
  enddo  ! t
  
end subroutine noah33_getirrigationstates
