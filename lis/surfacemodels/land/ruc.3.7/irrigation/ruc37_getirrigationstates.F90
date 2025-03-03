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
! !ROUTINE: RUC37_getirrigationstates
! \label{RUC37_getirrigationstates}
! 
! !INTERFACE:
subroutine RUC37_getirrigationstates(n,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use RUC37_lsmMod

! !DESCRIPTION:        
!
! Calculate water requirement and apply the amount to precipitation.
!
! Irrigate when root zone soil moisture falls below 50 % of 
! the field capacity (reference soil moiture) at 6 am LST.  
! The root zone is actual maximum root depth rather than NOAH root zone.
! Method of irrigation is by precipitation between 6-10 am LST.
!
! Irrigation amount is scaled to grid total crop fraction when intensity
! is less than the fraction.  Irrigation is expanded to non-crop, non-forest,
! non-baresoil/urban tiles if intensity exceeds grid total crop fraction.
! In latter case, scaled irrigation is applied to grassland first, 
! then further applied over the rest of tiles equally if the intensity 
! exceeds grassland fraction as well.   
!
! Optionally efficiency correction is applied to account for field loss. 
!
! Optionally outputs amount of water put into the system to a text file. 
!
! This version includes modifications to irr4 as follows:
! 1) Use location specific growing season threshold (40% of GFRAC range)
! 2) Allow irrigation in non-crop/non-forest tiles when irrigation 
!    intensity exceeds total crop fraction
!
! REVISION HISTORY:
!
! Aug 2008: Hiroko Kato; Initial code
! Nov 2012: Sujay Kumar, Incorporated into LIS
! Jun 2014: Ben Zaitchik; Added flood scheme
!EOP
  implicit none
  ! Sprinkler parameters
  real, parameter      :: otimess = 6.0 ! local trigger check start time [hour]
  real, parameter      :: irrhrs = 4.   ! duration of irrigation hours 
  ! Drip parameters (not currently implemented)
  real, parameter      :: otimeds = 6.0 ! local trigger check start time [hour]
  real, parameter      :: irrhrd = 12.0   ! duration of irrigation hours 
 ! Flood parameters
  real, parameter      :: otimefs = 6.0 ! local trigger check start time [hour]
  real, parameter      :: irrhrf = 1.0   ! duration of irrigation hours 
  !!!real, parameter      :: ffreq = 0.0 ! frequency of flood irrig [days] set to 0.0 to use thresh instead
  
  real, parameter      :: efcor = 0.0      ! Efficiency Correction (%)
  integer, parameter   :: nsoil = 9

  integer              :: n
  integer              :: rc
  integer              :: t,k,gid,vegt,l
  type(ESMF_State)     :: irrigState
  type(ESMF_Field)     :: irrigRateField,irrigFracField
  type(ESMF_Field)     :: irrigRootDepthField,irrigScaleField
  
  real,  pointer       :: irrigRate(:), irrigFrac(:)
  real,  pointer       :: irrigRootDepth(:), irrigScale(:)
  integer              :: chhr, lhr
  real                 :: asmc, tsmcwlt, tsmcref, ma, otimes, otimee, irrhr
  real                 :: sldpth(nsoil)
  real                 :: rdpth(nsoil)
  real                 :: zdpth(nsoil)
  real                 :: water(nsoil)
  real                 :: twater, twater1, twater2
  real                 :: ippix, crootd
  real                 :: smcmax, smcref, smcwlt, psisat, dksat
  real                 :: smcref1, smcwlt1,shdfac
  real                 :: bexp, smhigh, smlow
  integer              :: lroot,veg_index1,veg_index2
  real                 :: gsthresh, ltime

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

  irrigRate = 0.0  

!----------------------------------------------------------------------
! Set start and end times for selected irrigation type
!----------------------------------------------------------------------
  if(LIS_rc%irrigation_type.eq."Sprinkler") then
     otimes = otimess
     irrhr = irrhrs
     otimee = otimess + irrhrs
  elseif(LIS_rc%irrigation_type.eq."Drip") then
     otimes = otimeds
     irrhr = irrhrd
     otimee = otimeds + irrhrd
  elseif(LIS_rc%irrigation_type.eq."Flood") then
     otimes = otimefs
     irrhr = irrhrf
     otimee = otimefs + irrhrf
  endif
  
 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     twater  = 0.0
     water   = 0.0
     asmc    = 0.0
     tsmcwlt = 0.0
     tsmcref = 0.0
     ma      = 0.0
     crootd  = 0.0
     lroot   = 0

     sldpth(1) = 0.01         ! Soil layer thicknesses (m)
     sldpth(2) = 0.03
     sldpth(3) = 0.06
     sldpth(4) = 0.2
     sldpth(5) = 0.3         ! Soil layer thicknesses (m)
     sldpth(6) = 0.4
     sldpth(7) = 0.6
     sldpth(8) = 1.4
     sldpth(9) = 1.

     zdpth(1) = 0.         ! Soil layer depth from top (m)
     zdpth(2) = sldpth(1)
     zdpth(3) = zdpth(2) + sldpth(2)
     zdpth(4) = zdpth(3) + sldpth(3)
     zdpth(5) = zdpth(4) + sldpth(4)
     zdpth(6) = zdpth(5) + sldpth(5)
     zdpth(7) = zdpth(6) + sldpth(6)
     zdpth(8) = zdpth(7) + sldpth(7)
     zdpth(9) = zdpth(8) + sldpth(8)

     smcmax =  RUC37_struc(n)%ruc37(t)%qmax
     psisat =  RUC37_struc(n)%ruc37(t)%psis
     dksat  =  RUC37_struc(n)%ruc37(t)%ksat
     bexp   =  RUC37_struc(n)%ruc37(t)%bclh
     smlow = 0.5
     smhigh = 6.0
     smcref1 = smcmax*(5.79e-9/dksat)**(1.0/(2.0*bexp+3.0))
     smcref = smcref1 + (smcmax-smcref1) / smhigh
     smcwlt1 = smcmax * (200.0/psisat)**(-1.0/bexp)
     smcwlt = smcwlt1 - smlow * smcwlt1

     gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
     chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
     if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
          (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
          chhr = chhr -1
     lhr = LIS_rc%hr +chhr
     if(lhr.ge.24) lhr = lhr-24
     if(lhr.lt.0) lhr = lhr+24
                
     ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
      
     shdfac =  RUC37_struc(n)%ruc37(t)%shdfac
       
! Calculate vegetation and root depth parameters
   
     if((ltime.ge.otimes).and.(ltime.lt.otimee)) then 
        vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
       !----------------------------------------------------------------------       
       !    Proceed if it is non-forest, non-baresoil, non-urban
       !----------------------------------------------------------------------       
        if(LIS_rc%lcscheme.eq."UMD") then !UMD
           veg_index1 = 6
           veg_index2 = 11
        elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") then 
           veg_index1 = 6
           veg_index2 = 14
        elseif(LIS_rc%lcscheme.eq."USGS") then !UMD
           veg_index1 = 2
           veg_index2 = 10
        else
           write(LIS_logunit,*) 'The landcover scheme ',trim(LIS_rc%lcscheme)
           write(LIS_logunit,*) 'is not supported for irrigation '
           call LIS_endrun()
        endif
        
        if(vegt.ge.veg_index1.and.vegt.le.veg_index2&
             .and.vegt.ne.LIS_rc%bareclass.and.&
             vegt.ne.LIS_rc%urbanclass) then 
           if(irrigFrac(t).gt.0) then 
              ippix = irrigFrac(t)*0.01
              
              ! Determine the amount of irrigation to apply if irrigated tile
              if( IrrigScale(t).gt.0.0 ) then ! irrigated tile
!                if(ippix.gt.0.0) then  ! irrigated tile
                 
                 gsthresh = RUC37_struc(n)%ruc37(t)%shdmin + & 
                      0.40 * (RUC37_struc(n)%ruc37(t)%shdmax - &
                      RUC37_struc(n)%ruc37(t)%shdmin)
                 
                 if(shdfac .ge. gsthresh) then 
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
                      do k = 1,lroot-1
                       rdpth(k) = sldpth(k)
                      enddo
                       rdpth(lroot) = crootd - zdpth(lroot-1)
                    elseif ( crootd.ge.zdpth(4).and.crootd .lt. zdpth(5) ) then
                       lroot = 5
                      do k = 1,lroot-1
                       rdpth(k) = sldpth(k)
                      enddo
                       rdpth(lroot) = crootd - zdpth(lroot-1)
                    elseif ( crootd.ge.zdpth(5).and.crootd .lt. zdpth(6) ) then
                       lroot = 6
                      do k = 1,lroot-1
                       rdpth(k) = sldpth(k)
                      enddo
                       rdpth(lroot) = crootd - zdpth(lroot-1)
                    elseif ( crootd.ge.zdpth(6).and.crootd .lt. zdpth(7) ) then
                       lroot = 7
                      do k = 1,lroot-1
                       rdpth(k) = sldpth(k)
                      enddo
                       rdpth(lroot) = crootd - zdpth(lroot-1)
                    elseif ( crootd.ge.zdpth(7).and.crootd .lt. zdpth(8) ) then
                       lroot = 8
                      do k = 1,lroot-1
                       rdpth(k) = sldpth(k)
                      enddo
                       rdpth(lroot) = crootd - zdpth(lroot-1)

                    elseif ( crootd.ge.zdpth(8).and.crootd .lt. zdpth(9) ) then
                       lroot = 9
                      do k = 1,lroot-1
                       rdpth(k) = sldpth(k)
                      enddo
                       rdpth(lroot) = crootd - zdpth(lroot-1)

!                      else
!                         print*,'error getting root depth'
!                         stop
                    endif
       
                !!!!! SPRINKLER IRRIGATION
                    if(LIS_rc%irrigation_type.eq."Sprinkler") then
                       !----------------------------------------------------------------------       
                       !    Set the irrigation rate at start time; keep the value till next day
                       !    If local time at the tile fall in the irrigation check
                       !    hour then check the root zone average soil moisture
                       !----------------------------------------------------------------------       
                       if(ltime.eq.otimes) then 
                         !-------------------------------------------------------------
                         !     Compute the root zone accumlative soil moisture [mm], 
                         !     field capacity [mm], and wilting point [mm] 
                         !-------------------------------------------------------------
                          if(lroot.gt.0) then 
                             do k=1,lroot
                                asmc = asmc + RUC37_struc(n)%ruc37(t)%smc(k)*&
                                     rdpth(k)*1000.0
                                tsmcwlt = tsmcwlt + smcwlt * rdpth(k)*1000.0
                                tsmcref = tsmcref + smcref * rdpth(k)*1000.0
                             enddo
                         !---------------------------------------------------------------
                         !     Get the root zone moisture availability to the plant
                         !--------------------------------------------------------------- 
                             ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)
                             if(ma.le.LIS_rc%irrigation_thresh) then 
                                do k=1,lroot
                                   water(k) = &
                                        (smcref-RUC37_struc(n)%ruc37(t)%smc(k))*&
                                        rdpth(k)*1000.0
                                   twater = twater + water(k)
                                enddo
                                
                             !-----------------------------------------------------------------------------
                             !     Scale the irrigation intensity to the crop % when intensity < crop%.
                             !     Expand irrigation for non-crop, non-forest when intensity > crop %
                             !     in preference order of grassland first then rest.
                             !     *scale is pre-computed for each tile in getirrpmapetc module in a way
                             !     that is transparent to every tile irrigated or non-irrigated
                             !-----------------------------------------------------------------------------
                                twater1 = twater
                                twater = twater * irrigScale(t)
                                
                             !-----------------------------------------------------------------------------
                             !     Apply efficiency correction
                             !-----------------------------------------------------------------------------
                                twater2 = twater
                                twater = twater*(100.0/(100.0-efcor))
                             !-----------------------------------------------------------------------------
                             !     Compute irrigation rate
                             !-----------------------------------------------------------------------------
                                irrigRate(t) = twater/(irrhr*3600.0)
                           endif
                        endif
                     endif
                       !!!!! DRIP IRRIGATION (NOT CURRENTLY IMPLEMENTED)
                  elseif(LIS_rc%irrigation_type.eq."Drip") then
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
                         twater = twater*(100.0/(100.0-efcor))
                         !-----------------------------------------------------------------------------
                         !     Compute irrigation rate
                         !-----------------------------------------------------------------------------
                         irrigRate(t) = twater  ! for drip calculation, twater is a rate [kg m-2/s]
                         RUC37_struc(n)%ruc37(t)%smc(1) = &
                         RUC37_struc(n)%ruc37(t)%smc(1) + (twater-twater2)/(sldpth(1)*1000.0) !! check this with Sujay

                       !!!!! FLOOD IRRIGATION
                         elseif(LIS_rc%irrigation_type.eq."Flood") then
                         !-------------------------------------------------------------
                         !     Compute the root zone accumlative soil moisture [mm], 
                         !     field capacity [mm], and wilting point [mm] 
                         !-------------------------------------------------------------
                         if(lroot.gt.0) then 
                            do k=1,lroot
                               asmc = asmc + RUC37_struc(n)%ruc37(t)%smc(k)*&
                                    rdpth(k)*1000.0
                               tsmcwlt = tsmcwlt + smcwlt * rdpth(k)*1000.0
                               tsmcref = tsmcref + smcref * rdpth(k)*1000.0
                            enddo
                         !---------------------------------------------------------------
                         !     Get the root zone moisture availability to the plant
                         !--------------------------------------------------------------- 
!                            ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)   ! Original
                            ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)/IrrigScale(t) ! BZ added IrrigScale

                            if( ma .le. LIS_rc%irrigation_thresh ) then
                              do l = 1, LIS_rc%irrigation_mxsoildpth
                                 if( l == 1 ) then
                                   twater = (smcmax - RUC37_struc(n)%ruc37(t)%smc(l))*sldpth(l)*1000.0
                                 else
                                 ! BZ modification 4/2/2015 to saturate entire column and apply ippix 
                                   twater = twater + (smcmax - RUC37_struc(n)%ruc37(t)%smc(l))*sldpth(l)*1000.0
                                 endif
                              end do

                              !-----------------------------------------------------------------------------
                              !     Scale the irrigation intensity to the crop % when intensity < crop%.
                              !     Expand irrigation for non-crop, non-forest when intensity > crop %
                              !     in preference order of grassland first then rest.
                              !     *scale is pre-computed for each tile in getirrpmapetc module in a way
                              !     that is transparent to every tile irrigated or non-irrigated
                              !-----------------------------------------------------------------------------
                              twater1 = twater
                              twater = twater * irrigScale(t)
                              !-----------------------------------------------------------------------------
                              !     Apply efficiency correction
                              !-----------------------------------------------------------------------------
                              twater2 = twater
                              twater = twater*(100.0/(100.0-efcor))
                              !-----------------------------------------------------------------------------
                              !     Compute irrigation rate
                              !-----------------------------------------------------------------------------
                              irrigRate(t) = twater/LIS_rc%ts

                            ! BZ modification 4/2/2015 to account for ippix and all soil layers:
                               do l = 1, LIS_rc%irrigation_mxsoildpth
                                  RUC37_struc(n)%ruc37(t)%smc(l) = IrrigScale(t)*smcmax + &
                                                 (1-IrrigScale(t))*RUC37_struc(n)%ruc37(t)%smc(l)
                               end do
                            endif  
                         endif

                      endif

                   endif
                end if
             end if
          end if
       end if
    enddo
  end subroutine RUC37_getirrigationstates
