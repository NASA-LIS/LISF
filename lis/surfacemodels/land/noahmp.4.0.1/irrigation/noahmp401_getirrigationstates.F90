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
! !ROUTINE: noahmp401_getirrigationstates
! \label{noahmp401_getirrigationstates}
! 
! !INTERFACE:
subroutine noahmp401_getirrigationstates(n,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP401_lsmMod
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE, SMCWLT_TABLE, SMCREF_TABLE
  use LIS_vegDataMod, only: LIS_read_shdmin, LIS_read_shdmax
  use MODULE_SF_NOAHMPLSM_401, only: DVEG

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
! Aug 2016: Wanshu Nie; Incorporated into NoahMP
! May 2018: Wanshu Nie; Add temperature check for GRACE-DA purpose
! May 2019: Jessica Erlingis; Incorporate W. Nie's updates into LIS
!                             and add optional flag for groundwater abstraction
! Feb 2020: Jessica Erlingis; Correct sprinkler scheme so that it checks moisture
!                             at otimess and applies constant rate for irrhrs
! March 2020: Jessica Erlingis; Add to Noah-MP 4.0.1
! Apr 2021: Wanshu Nie; Add option to interact with DVEG
! May 2021: Wanshu Nie; update irrigation using ensemble mean when runing with DA.
!
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
  integer, parameter   :: nsoil = 4

  integer              :: n
  integer              :: rc
  integer              :: t,k,gid,vegt,l
  type(ESMF_State)     :: irrigState
  type(ESMF_Field)     :: irrigRateField,irrigFracField
  type(ESMF_Field)     :: irrigRootDepthField,irrigScaleField
  
  real,  pointer       :: irrigRate(:), irrigFrac(:)
  real,  pointer       :: irrigRootDepth(:), irrigScale(:)
  integer              :: chhr, lhr
  integer                :: soiltyp           ! soil type index [-]
  real                 :: asmc, tsmcwlt, tsmcref, ma, otimes, otimee, irrhr
  real                 :: sldpth(nsoil)
  real                 :: rdpth(nsoil)
  real                 :: zdpth(nsoil)
  real                 :: water(nsoil)
  real                 :: twater, twater1, twater2
  real                 :: ippix, crootd
  real                 :: smcmax, smcref, smcwlt
  !real                 :: smcref1, smcwlt1, shdfac
  real                 :: smhigh, smlow
  integer              :: lroot,veg_index1,veg_index2
  real                 :: gsthresh, ltime
  real                 :: shdfac, shdmin, shdmax
  real                 :: timestep, shift_otimes, shift_otimee
  real                 :: AWS
  real                 :: Dtime
  real, allocatable    :: placeshdmax(:,:), placeshdmin(:,:)
  real                 :: sfctemp, tempcheck

  type(ESMF_Field)     :: irriggwratioField
  real,  pointer       :: irriggwratio(:)

  integer              :: i, m
  real                 :: sfctemp_avg
  real                 :: shdfac_avg
  real                 :: smc_avg(nsoil)

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

  call ESMF_StateGet(irrigState, "Groundwater irrigation ratio",&
       irriggwratioField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Groundwater irrigation ratio')
  call ESMF_FieldGet(irriggwratioField, localDE=0,&
       farrayPtr=irriggwratio,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Groundwater irrigation ratio')

  allocate(placeshdmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(placeshdmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  call LIS_read_shdmax(n,placeshdmax)
  call LIS_read_shdmin(n,placeshdmin)

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
  


  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)

     sfctemp_avg = 0.
     shdfac_avg  = 0.
     smc_avg     = 0.

     do m=1,LIS_rc%nensem(n)

        t=(i-1)*LIS_rc%nensem(n)+m
        
        sfctemp_avg=sfctemp_avg+NOAHMP401_struc(n)%noahmp401(t)%sfctmp        
        shdfac_avg=shdfac_avg+NOAHMP401_struc(n)%noahmp401(t)%fveg            
        
        do k=1,nsoil
           
           smc_avg(k)=smc_avg(k)+NOAHMP401_struc(n)%noahmp401(t)%smc(k)
        
        end do
  
     end do

     sfctemp_avg=sfctemp_avg/LIS_rc%nensem(n)
     shdfac_avg=shdfac_avg/LIS_rc%nensem(n)

     do k=1,nsoil
  
        smc_avg(k)=smc_avg(k)/LIS_rc%nensem(n)
    
     end do


     do m=1,LIS_rc%nensem(n)

        t=(i-1)*LIS_rc%nensem(n)+m

     timestep = NOAHMP401_struc(n)%dt
     soiltyp = noahmp401_struc(n)%noahmp401(t)%soiltype
    
 
     ! Adjust bounds by timestep to account for the fact that LIS_rc%hr, etc. 
     ! will represents the END of the integration timestep window


     shift_otimes = otimes + (timestep/3600.)
     shift_otimee = otimee + (timestep/3600.)

     twater  = 0.0
     water   = 0.0
     asmc    = 0.0
     tsmcwlt = 0.0
     tsmcref = 0.0
     ma      = 0.0
     crootd  = 0.0
     lroot   = 0

     !JE this code block will need to be changed to account for variable
     ! soil layers in Noah-MP

     sldpth(1) = 0.1         ! Soil layer thicknesses (m)
     sldpth(2) = 0.3
     sldpth(3) = 0.6
     sldpth(4) = 1.0
     zdpth(1) = sldpth(1)         ! Soil layer depth from top (m)
     zdpth(2) = sldpth(1) + sldpth(2)
     zdpth(3) = sldpth(1) + sldpth(2) + sldpth(3)
     zdpth(4) = sldpth(1) + sldpth(2) + sldpth(3) + sldpth(4)

     smcmax = SMCMAX_TABLE(soiltyp)
     smcref = SMCREF_TABLE(soiltyp)
     smcwlt = SMCWLT_TABLE(soiltyp)

  !   sfctemp = NOAHMP401_struc(n)%noahmp401(t)%sfctmp
     tempcheck = 273.16 + 2.5

     gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
     chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
     if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
          (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
          chhr = chhr -1
     lhr = LIS_rc%hr +chhr
     if(lhr.ge.24) lhr = lhr-24
     if(lhr.lt.0) lhr = lhr+24
                
     ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
    
     if((DVEG == 2 .OR. DVEG == 5 .OR. DVEG == 6) .AND. LIS_rc%irrigation_dveg == 1) then
        shdfac = shdfac_avg

     else
        shdfac =  NOAHMP401_struc(n)%noahmp401(t)%shdfac_monthly(LIS_rc%mo)
     end if


   ! If we are outside of the irrigation window, set rate to 0
     if ((ltime.gt.shift_otimee).or.(ltime.lt.shift_otimes)) then
       irrigRate(t) = 0.0
     endif

! Calculate vegetation and root depth parameters
   
! JE This temperature check avoids irrigating at temperatures near or below 0C
     if((ltime.ge.shift_otimes).and.(ltime.le.shift_otimee).and. &
         (sfctemp_avg.gt.tempcheck)) then 
!------------------------------------
        vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
       !----------------------------------------------------------------------       
       !    Proceed if it is non-forest, non-baresoil, non-urban
       !----------------------------------------------------------------------       
        if(LIS_rc%lcscheme.eq."UMD") then !UMD
           veg_index1 = 6
           veg_index2 = 11
        elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") &
             then 
           veg_index1 = 6
           veg_index2 = 14
        elseif(LIS_rc%lcscheme.eq."USGS") then !UMD
           veg_index1 = 2
           veg_index2 = 10
        else
           write(LIS_logunit,*) '[ERR] The landcover scheme ',&
                trim(LIS_rc%lcscheme),' is not supported for irrigation '
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
                 
                 !shdmin = minval(NOAHMP36_struc(n)%noahmp36(t)%shdfac_monthly)
                 !shdmax = maxval(NOAHMP36_struc(n)%noahmp36(t)%shdfac_monthly)
                 shdmin =placeshdmin(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,     &
                                                       LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
                 shdmax =placeshdmax(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,     &
                                        LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
                 
               ! let gsthresh be a function of the range, which means the larger
               ! the range is, the higher GVF threshold will be for this grid.                           
               ! JE Gsthresh is a GVF threshold used to identify a growing season for each
               ! pixel and allow irrigation during that time
                  gsthresh = shdmin + & 
                      (LIS_rc%irrigation_GVFparam1 + LIS_rc%irrigation_GVFparam2*&
                       (shdmax-shdmin)) * (shdmax - shdmin)


                 !JE Changes needed to this code block to account for variable soil layers
                 ! in Noah-MP
                 
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
                       rdpth(1) = sldpth(1)
                       rdpth(2) = sldpth(2)
                       rdpth(3) = sldpth(3)
                       rdpth(4) = crootd - zdpth(3)
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
                       if(ltime.eq.shift_otimes) then !Check moisture availability at otimes only
                         !-------------------------------------------------------------
                         !     Compute the root zone accumlative soil moisture [mm], 
                         !     field capacity [mm], and wilting point [mm] 
                         !-------------------------------------------------------------
                          if(lroot.gt.0) then 
                             do k=1,lroot
                                asmc = asmc + smc_avg(k)*rdpth(k)*1000.0
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
                                        (smcref-smc_avg(k))*rdpth(k)*1000.0
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
                         irrigRate(t) = twater  ! for drip calculation, twater is a rate [kg/m2/s]
                         NOAHMP401_struc(n)%noahmp401(t)%smc(1) = &
                         NOAHMP401_struc(n)%noahmp401(t)%smc(1) + (twater-twater2)/(sldpth(1)*1000.0) !! check this with Sujay

                       !!!!! FLOOD IRRIGATION
                         elseif(LIS_rc%irrigation_type.eq."Flood") then
                         !-------------------------------------------------------------
                         !     Compute the root zone accumlative soil moisture [mm], 
                         !     field capacity [mm], and wilting point [mm] 
                         !-------------------------------------------------------------
                         if(lroot.gt.0) then 
                            do k=1,lroot
                               asmc = asmc + NOAHMP401_struc(n)%noahmp401(t)%smc(k)*&
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
                                   twater = (SMCMAX - NOAHMP401_struc(n)%noahmp401(t)%smc(l))*sldpth(l)*1000.0
                                 else
                                 ! BZ modification 4/2/2015 to saturate entire column and apply ippix 
                                   twater = twater + (smcmax - NOAHMP401_struc(n)%noahmp401(t)%smc(l))*sldpth(l)*1000.0
!                                 twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(2))*sldpth(2)*1000.0
!                                 twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(3))*sldpth(3)*1000.0
!                                 twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(4))*sldpth(4)*1000.0
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
                                  NOAHMP401_struc(n)%noahmp401(t)%smc(l) = IrrigScale(t)*smcmax + &
                                                 (1-IrrigScale(t))*NOAHMP401_struc(n)%noahmp401(t)%smc(l)
                               end do
                            endif  
                         endif

                      endif

                   endif
                end if
             end if
          end if
  
           ! Remove irrigated water from groundwater 
           !JE Add in flag to turn groundwater abstraction on/off
           if (LIS_rc%irrigation_GWabstraction.eq.1) then
              AWS = NOAHMP401_struc(n)%noahmp401(t)%wa
              Dtime = NOAHMP401_struc(n)%ts
              if (LIS_rc%irrigation_SourcePartition.eq.1) then
                  if(irriggwratio(t).gt.0) then
                  NOAHMP401_struc(n)%noahmp401(t)%wa = AWS - irrigRate(t)*Dtime*irriggwratio(t)/100

                  end if
              else
                 NOAHMP401_struc(n)%noahmp401(t)%wa = AWS - irrigRate(t)*Dtime
              end if
           end if
       end if

    end do
  end do

  end subroutine noahmp401_getirrigationstates
