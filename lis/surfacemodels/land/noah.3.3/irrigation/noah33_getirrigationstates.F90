!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
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
subroutine noah33_getirrigationstates(nest,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use noah33_lsmMod
  use LIS_irrigationMod
  use IRRIGATION_MODULE
  
! !DESCRIPTION:        
!
! This module is a LSM interface that collects NOAH LSM specific states and 
! parameters needed for the irrigation model.  Root zone soil moisture, 
! greenness/LAI states, soil parameters, and vegetation parameters of a given 
! tile are used to calculate water requirement of a crop, based on the 
! irrigation type predetermined for the crop tile.  The module calls irrigation
! routines that return irrigation rate, then applies the water via 
! precipitation (Sprinkler) or directly modify soil moisture content 
! (Drip, Flood, and Paddy).  
!
! The root zone is actual maximum root depth rather than NOAH's vegetation 
! parameter "NROOT". The root depth grows following the GVF seasonality to 
! reflect the water demand increase/decrease for crops. (note: need age of
! trees for orchard).
! 
! Growing season threshold (i.e. 40% of GFRAC range) can be changed in 
! lis.config.
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
! Dec 2021: Sarith Mahanama; Separated into general irrigation model and LSM 
!                            dependant components.
! Feb 2022: Hiroko Beaudoing; Added irrigation application section.
!                             
!EOP
  implicit none

  integer, intent(in)  :: nest
  integer              :: rc
  integer              :: TileNo,tid,gid,vegt,l
  type(ESMF_State)     :: irrigState
  real                 :: sldpth(noah33_struc(nest)%nslay)
  real                 :: rdpth(noah33_struc(nest)%nslay)
  real                 :: zdpth(noah33_struc(nest)%nslay)
  real                 :: crootd
  integer              :: lroot,veg_index1,veg_index2, nlctypes
  integer              :: croptype
  real                 :: gsthresh
  logical              :: irrig_check_frozen_soil
  real                 :: amount
  integer              :: checktid
  type(irrigation_model) :: IM 

  ! _______________________________________________________

  ! Get irrigation state variables
  ! -----------------------------------------
  call IM%get_irrigstate (nest, irrigState)

  ! Set vegetation type index to be irrigated
  ! -----------------------------------------
  call IM%get_irrig_vegindex (veg_index1, veg_index2, nlctypes)
  
  ! Set global soil  parameters
  sldpth(1) = noah33_struc(nest)%lyrthk(1)         ! Soil layer thicknesses (m)
  sldpth(2) = noah33_struc(nest)%lyrthk(2)
  sldpth(3) = noah33_struc(nest)%lyrthk(3)
  sldpth(4) = noah33_struc(nest)%lyrthk(4)
  zdpth(1) = sldpth(1)         ! Soil layer depth from top (m)
  zdpth(2) = sldpth(1) + sldpth(2)
  zdpth(3) = sldpth(1) + sldpth(2) + sldpth(3)
  zdpth(4) = sldpth(1) + sldpth(2) + sldpth(3) + sldpth(4)
  
  !---------------------------------------------------------------
  ! Main tile loop
  !---------------------------------------------------------------
  
  TILE_LOOP: do TileNo = 1,LIS_rc%npatch(nest,LIS_rc%lsm_index)
     
     gid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%index
     tid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%tile_id

     ! frozen tile check
     irrig_check_frozen_soil = .false.
     
     if((noah33_struc(nest)%noah(TileNo)%smc(1) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(1)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif((noah33_struc(nest)%noah(TileNo)%smc(2) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(2)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif((noah33_struc(nest)%noah(TileNo)%smc(3) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(3)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif((noah33_struc(nest)%noah(TileNo)%smc(4) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(4)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif(noah33_struc(nest)%noah(TileNo)%stc(2).le.LIS_CONST_TKFRZ) then
        irrig_check_frozen_soil = .true. 
     elseif(noah33_struc(nest)%noah(TileNo)%stc(3).le.LIS_CONST_TKFRZ) then
        irrig_check_frozen_soil = .true. 
     elseif(noah33_struc(nest)%noah(TileNo)%stc(4).le.LIS_CONST_TKFRZ) then
        irrig_check_frozen_soil = .true. 
     endif
     
     ! Process only non-frozen tiles
     FROZEN: if(.not.irrig_check_frozen_soil) then
        
        ! Process only irrigated tiles
        IRRF: if(IM%irrigFrac(TileNo).gt.0) then
           
           ! Determine the amount of irrigation to apply if irrigated tile
           IRRS: if( IM%IrrigScale(TileNo).gt.0.0 ) then ! irrigated tile
           
              ! Proceed if it is non-forest, non-baresoil, non-urban
              vegt = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%vegt  
              VEGIF: if(vegt.ge.veg_index1.and.vegt.le.veg_index2     &
                   .and.vegt.ne.LIS_rc%bareclass                      &
                   .and.vegt.ne.LIS_rc%urbanclass) then

                 ! Calculate active root depth
                 ! ---------------------------
                 ! make sure greenness data is not zero
                 if ( noah33_struc(nest)%noah(TileNo)%shdfac.eq.0.0) then 
                   crootd = IM%irrigRootdepth(TileNo)
                 else
                   crootd = IM%irrigRootdepth(TileNo)*noah33_struc(nest)%noah(TileNo)%shdfac
                 endif
                 lroot  = 0
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

                 if (lroot == 0 .and. LIS_rc%tscount(nest) == 1 ) then
                    ! HKB: keep running with warning, do not stop 
                    write(LIS_logunit,*) '[WARN] lroot should be > 0!'
                    write(LIS_logunit,*) TileNo,gid,vegt,crootd,&
                     IM%irrigRootdepth(TileNo), &
                     noah33_struc(nest)%noah(TileNo)%shdfac, &
                     LIS_domain(nest)%grid(LIS_domain(nest)%gindex( &
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,&
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row))%lat, &
                     LIS_domain(nest)%grid(LIS_domain(nest)%gindex( &
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,&
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row))%lon
                     !call LIS_endrun()
                 endif

                 ! compute vegetation threshold for the trigger
                 ! --------------------------------------------
                 
                 !gsthresh = noah33_struc(nest)%noah(TileNo)%shdmin +   & 
                 !     LIS_irrig_struc(nest)%veg_thresh * (noah33_struc(nest)%noah(TileNo)%shdmax - &
                 !     noah33_struc(nest)%noah(TileNo)%shdmin)
                 ! let gsthresh be a function of the range, which means the larger
                 ! the range is, the higher GVF threshold will be for this grid.  
                  gsthresh = noah33_struc(nest)%noah(TileNo)%shdmin + &
                      (LIS_irrig_struc(nest)%irrigation_GVFparam1 + &
                       LIS_irrig_struc(nest)%irrigation_GVFparam2* &
                       (noah33_struc(nest)%noah(TileNo)%shdmax- &
                        noah33_struc(nest)%noah(TileNo)%shdmin)) * &
                       (noah33_struc(nest)%noah(TileNo)%shdmax - &
                        noah33_struc(nest)%noah(TileNo)%shdmin)

                 ! get irrigation rates from the irrigation model
                 ! ----------------------------------------------
           
                 croptype = vegt - nlctypes

                 call IM%update_irrigrate (                                             &
                      nest,TileNo, croptype, LIS_domain(nest)%grid(gid)%lon,     &
                      noah33_struc(nest)%noah(TileNo)%shdfac,gsthresh,                  &
                      noah33_struc(nest)%noah(TileNo)%smcwlt,                           &
                      noah33_struc(nest)%noah(TileNo)%smcmax,                           &
                      noah33_struc(nest)%noah(TileNo)%smcref,                           &
                      noah33_struc(nest)%noah(TileNo)%smc(:lroot),                      &
                      rdpth(:lroot),noah33_struc(nest)%noah(TileNo)%soiltype)
                                   
              endif VEGIF
           endif IRRS
        endif IRRF
     endif FROZEN
  end do TILE_LOOP

  ! Update land surface model's moisture state
  ! ------------------------------------------
  TILE_LOOP2: do TileNo = 1,LIS_rc%npatch(nest,LIS_rc%lsm_index)
        vegt = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%vegt  
        croptype = vegt - nlctypes
        ! Process only irrigated tiles
        IRRF2: if(IM%irrigFrac(TileNo).gt.0) then
           IRRS2: if( IM%IrrigScale(TileNo).gt.0.0 ) then ! irrigated tile
              ! Apply irrigation rate-- zero outside of irrigation hour, 
              ! not in season, or not in deficit
              IRRON: if ( IM%irrigRate(TileNo).gt.0 ) then

                !SPRINKLER application is in alltypes_irrigation_updates 

                DRIP: if ( IM%irrigType(TileNo) == 2 ) then
                ! convert rate (mm/s) to accumulation (mm or kg m-2) and 
                ! then volumetric and apply to surface layer
                ! if root zone moisture becomes field capacity at the next 
                ! time step, it will be shut off 
                   amount = IM%irrigRate(TileNo)*LIS_rc%ts/sldpth(1)/1000.
                   noah33_struc(nest)%noah(TileNo)%smc(1) =  &
                        amount + noah33_struc(nest)%noah(TileNo)%smc(1)      
                endif DRIP

                FLOOD: if ( IM%irrigType(TileNo) == 3 ) then
                   if ( croptype == LIS_rc%ricecrop ) then   ! rice
                   ! PADDY-- surface layer only
                        noah33_struc(nest)%noah(TileNo)%smc(1) =  &
                             IM%IrrigScale(TileNo)*noah33_struc(nest)%noah(TileNo)%smcmax + &
                             (1-IM%IrrigScale(TileNo))*noah33_struc(nest)%noah(TileNo)%smc(1)
                   else
                   ! NON-PADDY -- surface only or  multiple layers
                   ! BZ modification 4/2/2015 to account for ippix and all soil layers:
                   ! raise SM to saturation instantly and keep it saturated
                   ! throughout the irrigation duration 
                      do l = 1, LIS_irrig_struc(nest)%irrigation_mxsoildpth
                        noah33_struc(nest)%noah(TileNo)%smc(l) =  &
                             IM%IrrigScale(TileNo)*noah33_struc(nest)%noah(TileNo)%smcmax + &
                             (1-IM%IrrigScale(TileNo))*noah33_struc(nest)%noah(TileNo)%smc(l)
                      end do
                   endif
                endif FLOOD

              endif   IRRON
           endif IRRS2
        endif IRRF2
  end do TILE_LOOP2
    
end subroutine noah33_getirrigationstates
