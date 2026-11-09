!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
!
! !ROUTINE: noahmp36_getirrigationstates
! \label{noahmp36_getirrigationstates}
!
! !INTERFACE:
subroutine noahmp36_getirrigationstates(nest,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP36_lsmMod
  use MODULE_SF_NOAHMPLSM_36, only: MAXSMC, REFSMC, WLTSMC
  use LIS_vegDataMod, only: LIS_read_shdmin, LIS_read_shdmax
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_irrigationMod
  use IRRIGATION_MODULE
! use MODULE_SF_NOAHMPLSM_36, only: SMCMAX, PSISAT, DKSAT, BEXP

! !DESCRIPTION:
!
! This module is a LSM interface that collects NOAHMP LSM specific states and
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
! Following irrigation settings per Sprinkler, Drip, and Flood are set in
! lis.config:
! Local trigger check start time,
! Duration of irrigation hours
! Frequency of irrigation [days]
! Efficiency Correction (%)
! Irrigation trigger threshold
! Growing season threshold (i.e. 40% of GFRAC range)
!
! REVISION HISTORY:
!
! Aug 2008: Hiroko Kato; Initial code
! Nov 2012: Sujay Kumar, Incorporated into LIS
! Jun 2014: Ben Zaitchik; Added flood scheme
! Aug 2016: Wanshu Nie; Incorporated into NoahMP
! May 2018: Wanshu Nie; Add temperature check for GRACE-DA purpose
! May 2019: Jessica Erlingis; Incorporate W. Nie's updates into LIS
!                             and add optional flag for groundwater abstraction
! Feb 2020: Jessica Erlingis; Correct sprinkler scheme so that it checks moisture
!                             at otimess and applies constant rate for irrhrs
! Feb 2026: Hiroko Beaudoing; Updated to the irrigation scheme of Rodell et al 2026
! Aug 2026: Hiroko Beaudoing: Do not irrigate if soil is frozen. Enforce consitent
!                             smc and sh2o in furrow.

!EOP
  implicit none

  integer              :: nest
  integer              :: rc
  integer              :: TileNo,gid,vegt,l,tid
  type(ESMF_State)     :: irrigState
  real                 :: sldpth(NOAHMP36_struc(nest)%nsoil)
  real                 :: rdpth(NOAHMP36_struc(nest)%nsoil)
  real                 :: zdpth(NOAHMP36_struc(nest)%nsoil)
  real                 :: crootd
  integer              :: lroot,veg_index1,veg_index2,nlctypes
  integer              :: croptype
  real                 :: gsthresh
  real                 :: shdfac, shdmin, shdmax
! _______________wanshu  add GW extraction_______________________________
  real                 :: AWS
  real                 :: Dtime

!---------------wanshu readin shdmax and shdmin from map----------------

  real, allocatable    :: placeshdmax(:,:), placeshdmin(:,:)
!--------------wanshu-----add temp check-------
  real                 :: sfctemp, tempcheck
!--------------new irrigation scheme
  real                 :: amount,added,totamount
  integer              :: i, m
  real                 :: sfctemp_avg
  real                 :: shdfac_avg
  real                 :: tile_sat_target, SICE
  type(ESMF_Field)   :: tmpField
  real,pointer       :: tmp(:)
  integer            :: status
  type(irrigation_model) :: IM
  ! _______________________________________________________

  ! Get irrigation state variables
  ! -----------------------------------------
  call IM%get_irrigstate (nest, irrigState)

  ! Set vegetation type index to be irrigated
  ! -----------------------------------------
  call IM%get_irrig_vegindex (veg_index1, veg_index2, nlctypes)

 !-------------wanshu-----call LIS_read_shdmax and min------------
  allocate(placeshdmax(LIS_rc%lnc(nest),LIS_rc%lnr(nest)))
  allocate(placeshdmin(LIS_rc%lnc(nest),LIS_rc%lnr(nest)))

  call LIS_read_shdmax(nest,placeshdmax)
  call LIS_read_shdmin(nest,placeshdmin)

  ! Set global soil  parameters
  sldpth(1) = NOAHMP36_struc(nest)%sldpth(1)
  sldpth(2) = NOAHMP36_struc(nest)%sldpth(2)
  sldpth(3) = NOAHMP36_struc(nest)%sldpth(3)
  sldpth(4) = NOAHMP36_struc(nest)%sldpth(4)
  zdpth(1) = sldpth(1)         ! Soil layer depth from top (m)
  zdpth(2) = sldpth(1) + sldpth(2)
  zdpth(3) = sldpth(1) + sldpth(2) + sldpth(3)
  zdpth(4) = sldpth(1) + sldpth(2) + sldpth(3) + sldpth(4)

 ! Set frozen soil temperature threshold (currently experimental at 2.5 C)
  tempcheck = 273.16 + 2.5

  call ESMF_StateGet(LIS_FORC_State(nest),(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status,'noahmp36_getirrigationstates: error getting Tair')
  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr= tmp,rc=status)
  call LIS_verify(status,'noahmp36_getirrigationstates: error retrieving tmp')

  !---------------------------------------------------------------
  ! Main tile loop
  !---------------------------------------------------------------

  TILE_LOOP: do i=1,LIS_rc%npatch(nest,LIS_rc%lsm_index)/LIS_rc%nensem(nest)

     sfctemp_avg = 0.
     shdfac_avg  = 0.

     do m=1,LIS_rc%nensem(nest)

        TileNo=(i-1)*LIS_rc%nensem(nest)+m

        if ( tmp(TileNo) .gt. 330 ) then
        !HKB: TODO--commented out below because forc_count is zero before calling noahmp36_f2t.
        !tair needs to be divided by forc_count if using more than one forcing sources (i.e. tair is a sum)
        !sfctmp = NOAHMP36_struc(n)%noahmp36(t)%tair/NOAHMP36_struc(n)%forc_count
         write(LIS_logunit, *) "[ERR] need to account for forc_count stopping"
                call LIS_endrun()
        endif
        sfctemp_avg=sfctemp_avg+tmp(TileNo)
        shdfac_avg=shdfac_avg+NOAHMP36_struc(nest)%noahmp36(TileNo)%shdfac_monthly(LIS_rc%mo)

     end do

     sfctemp_avg=sfctemp_avg/LIS_rc%nensem(nest)
     shdfac_avg=shdfac_avg/LIS_rc%nensem(nest)

     ENS_LOOP: do m=1,LIS_rc%nensem(nest)

       TileNo = (i-1)*LIS_rc%nensem(nest)+m
       gid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%index
       tid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%tile_id

       ! Process only non-frozen tiles-- air temp can be above zero but soil frozen
       FROZEN: if(sfctemp_avg.gt.tempcheck) then

        ! Process only irrigated tiles
        IRRF: if(IM%irrigFrac(TileNo).gt.0) then

           ! Determine the amount of irrigation to apply if irrigated tile
           IRRS: if( IM%IrrigScale(TileNo).gt.0.0 ) then ! irrigated tile

              ! Proceed if it is non-forest, non-baresoil, non-urban
              vegt = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%vegt
              shdfac =  NOAHMP36_struc(nest)%noahmp36(TileNo)%shdfac_monthly(LIS_rc%mo)
              VEGIF: if(vegt.ge.veg_index1.and.vegt.le.veg_index2     &
                   .and.vegt.ne.LIS_rc%bareclass                      &
                   .and.vegt.ne.LIS_rc%urbanclass) then

                 ! Calculate active root depth
                 ! ---------------------------
                 ! make sure greenness data is not zero
                 if ( shdfac.eq.0.0) then
                   crootd = IM%irrigRootdepth(TileNo)
                 else
                   crootd = IM%irrigRootdepth(TileNo)*shdfac
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
                    write(LIS_logunit,*) '[WARN] lroot should be > 0!'
                    write(LIS_logunit,*) TileNo,gid,vegt,crootd,&
                     IM%irrigRootdepth(TileNo), &
                     shdfac, &
                     LIS_domain(nest)%grid(LIS_domain(nest)%gindex( &
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,&
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row))%lat, &
                     LIS_domain(nest)%grid(LIS_domain(nest)%gindex( &
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,&
                      LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row))%lon
                     call LIS_endrun()
                 endif

                 ! compute vegetation threshold for the trigger
                 ! --------------------------------------------
                 shdmin =placeshdmin(LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,     &
                                     LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row)
                 shdmax =placeshdmax(LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,     &
                                     LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row)

                 ! write(*,*) 'test for greenness fraction value:', shdfac, shdmin, shdmax
                 ! let gsthresh be a function of the range, which means the larger
                 ! the range is, the higher GVF threshold will be for this grid.
                 ! JE Gsthresh is a GVF threshold used to identify a growing season for each
                 ! pixel and allow irrigation during that time
                 gsthresh = shdmin + &
                      (LIS_irrig_struc(nest)%irrigation_GVFparam1 + LIS_irrig_struc(nest)%irrigation_GVFparam2*&
                       (shdmax-shdmin)) * (shdmax - shdmin)

                 ! get irrigation rates from the irrigation model
                 ! ----------------------------------------------

                 croptype = vegt - nlctypes

                 ! make sure surface soil is not frozen
                 SICE = NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(1)-NOAHMP36_struc(nest)%noahmp36(TileNo)%sh2o(1)
                 ! distinguish paddy from others
                 ! ----------------------------------------------
                 PADDY: if ( IM%irrigType(TileNo) == 3 .and. croptype == LIS_rc%ricecrop) then
                   ! PADDY-- use sh2o and soil layer depths for irrigateinput
                    ! make sure surface soil is not frozen
                    if ( SICE .eq. 0 ) then
                     call IM%update_irrigrate (                                  &
                      nest,TileNo, croptype, LIS_domain(nest)%grid(gid)%lon,     &
                      shdfac,gsthresh,                  &
                      WLTSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype),                           &
                      MAXSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype),                           &
                      REFSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype),                           &
                      NOAHMP36_struc(nest)%noahmp36(TileNo)%sh2o(:lroot),                      &
                      sldpth(:lroot),NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype)
                    else
                     IM%irrigRate(TileNo) = 0.0
                    endif
                 else
                    ! make sure surface soil is not froze
                    if ( SICE .eq. 0 ) then
                     call IM%update_irrigrate (                                  &
                      nest,TileNo, croptype, LIS_domain(nest)%grid(gid)%lon,     &
                      shdfac,gsthresh,                  &
                      WLTSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype),                           &
                      MAXSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype),                           &
                      REFSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype),                           &
                      NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(:lroot),                      &
                      rdpth(:lroot),NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype)
                    else
                     IM%irrigRate(TileNo) = 0.0
                    endif
                 endif PADDY

              endif VEGIF
           endif IRRS
        endif IRRF
       endif FROZEN
    end do ENS_LOOP
  end do TILE_LOOP

  ! Update land surface model's moisture state
  ! ------------------------------------------
  TILE_LOOP2: do i=1,LIS_rc%npatch(nest,LIS_rc%lsm_index)/LIS_rc%nensem(nest)
   ENS_LOOP2: do m=1,LIS_rc%nensem(nest)
        TileNo = (i-1)*LIS_rc%nensem(nest)+m
        gid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%index
        vegt = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%vegt
        croptype = vegt - nlctypes
        ! Process only irrigated tiles
        IRRF2: if(IM%irrigFrac(TileNo).gt.0) then
           IRRS2: if( IM%IrrigScale(TileNo).gt.0.0 ) then ! irrigated tile
              ! Apply irrigation rate-- zero outside of irrigation hour,
              ! not in season, or not in deficit
              IRRON: if ( IM%irrigRate(TileNo).gt.0 ) then

                !SPRINKLER application is in alltypes_irrigation_updates
                SPRINKLER: if ( IM%irrigType(TileNo) == 1 ) then
                  IM%irrigAppRate(TileNo) = IM%irrigRate(TileNo)
                endif SPRINKLER

                DRIP: if ( IM%irrigType(TileNo) == 2 ) then
                ! convert rate (mm/s) to accumulation (mm or kg m-2) and
                ! then volumetric and apply to surface layer
                ! if root zone moisture becomes field capacity at the next
                ! time step, it will be shut off
                   amount = IM%irrigRate(TileNo)*LIS_rc%ts/sldpth(1)/1000.
                ! HKB: application need to weight by irrigated fraction of tile
                   NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(1) =  &
                        IM%IrrigScale(TileNo)*(amount + NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(1)) &
                        + (1-IM%IrrigScale(TileNo))*NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(1)
                   NOAHMP36_struc(nest)%noahmp36(TileNo)%sh2o(1) =  &
                        IM%IrrigScale(TileNo)*(amount + NOAHMP36_struc(nest)%noahmp36(TileNo)%sh2o(1)) &
                        + (1-IM%IrrigScale(TileNo))*NOAHMP36_struc(nest)%noahmp36(TileNo)%sh2o(1)
                   IM%irrigAppRate(TileNo) = IM%irrigRate(TileNo)
                endif DRIP

                FLOOD: if ( IM%irrigType(TileNo) == 3 ) then
                   if ( croptype == LIS_rc%ricecrop ) then   ! rice
                   ! PADDY-- surface layer only
                   ! Bug fix--compute applied amount first!!!!! qri, qritest2
                   ! NEW PADDY-- if smc is below target, irrigRate > 0. Reset smc to target.
                      tile_sat_target = IM%IrrigScale(TileNo)*MAXSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype)+&
                                  (1.0-IM%IrrigScale(TileNo))*REFSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype)

                      NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(1) = tile_sat_target
                      NOAHMP36_struc(nest)%noahmp36(TileNo)%sh2o(1) = tile_sat_target

                   else
                   ! NON-PADDY -- surface only or  multiple layers
                   ! BZ modification 4/2/2015 to account for ippix and all soil layers:
                   ! raise SM to saturation instantly and keep it saturated
                   ! throughout the irrigation duration
                   ! HKB: make sure the rate is high enough to bring it to
                   ! saturation in volmetric
                     totamount = 0.0
                     do l = 1, LIS_irrig_struc(nest)%irrigation_mxsoildpth
                      amount = IM%irrigRate(TileNo)*LIS_rc%ts/sldpth(l)/1000.
                      added = amount + NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(l)
                      if ( added >= MAXSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype) ) then
                        !add the amount
                        totamount = totamount + amount
                        NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(l) =  &
                           IM%IrrigScale(TileNo)*added + &
                           (1-IM%IrrigScale(TileNo))*NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(l)
                      else
                        !bring to saturation
                        totamount = totamount + &
                           IM%IrrigScale(TileNo)* &
                           (MAXSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype) - NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(l)) &
                           *sldpth(l)*1000.       ! mm
                        NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(l) =  &
                           IM%IrrigScale(TileNo)*MAXSMC(NOAHMP36_struc(nest)%noahmp36(TileNo)%soiltype) + &
                           (1-IM%IrrigScale(TileNo))*NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(l)
                      endif
                      !HKB: make smc and sh2o are consistent 8/7/2026
                      NOAHMP36_struc(nest)%noahmp36(TileNo)%sh2o(l) = &
                           NOAHMP36_struc(nest)%noahmp36(TileNo)%smc(l)
                     end do
                     IM%irrigAppRate(TileNo) = totamount/LIS_rc%ts
                   endif
                endif FLOOD

              else
                IM%irrigAppRate(TileNo) = 0.0
              endif   IRRON

              ! Remove irrigated water from groundwater
              !JE Add in flag to turn groundwater abstraction on/off
              if (LIS_irrig_struc(nest)%irrigation_GWabstraction.eq.1) then
                 AWS = NOAHMP36_struc(nest)%noahmp36(TileNo)%wa
                 Dtime = NOAHMP36_struc(nest)%dt
                 NOAHMP36_struc(nest)%noahmp36(TileNo)%wa = AWS - IM%irrigAppRate(TileNo)*Dtime
              end if

           endif IRRS2
        endif IRRF2

   end do ENS_LOOP2
  end do TILE_LOOP2
  end subroutine noahmp36_getirrigationstates
