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
! !ROUTINE: clsmf25_getirrigationstates
! \label{clsmf25_getirrigationstates}
! 
! !INTERFACE:
subroutine clsmf25_getirrigationstates(n,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use clsmf25_lsmMod
  use LIS_vegDataMod

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
! Aug 2008: Hiroko Kato; Initial code for Noah LSM
! Feb 2014: Sujay Kumar; Implemented in LIS based on the work of
!           John Bolten and student. 
! Jul 2014: Ben Zaitchik; added flood routine

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
  
  real, parameter      :: efcor = 76.0      ! Efficiency Correction (%)
  integer, parameter   :: nsoil = 4

  integer              :: n
  integer              :: rc
  integer              :: t,k,gid,tid,vegt
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
  real                 :: smcref1, smcwlt1,laifac,lai
  real                 :: bexp, smhigh, smlow
  integer              :: lroot,veg_index1,veg_index2
  real                 :: laithresh,ltime
  real,      allocatable   :: laimax(:,:),laimin(:,:)
  
  if(clsmf25_struc(n)%modelStart) then 
     clsmf25_struc(n)%modelStart = .false. 

     allocate(laimax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     allocate(laimin(LIS_rc%lnc(n),LIS_rc%lnr(n)))

     call LIS_read_laimax(n,laimax)
     call LIS_read_laimin(n,laimin)

     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if (laimax(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,       &
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row).ne.-9999.00) then
           clsmf25_struc(n)%cat_param(t)%laimax =            &
                laimax(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,     &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
        else
           clsmf25_struc(n)%cat_param(t)%laimax = 1.0
        endif
        if (laimin(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,        &
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row).ne.-9999.00) then
           clsmf25_struc(n)%cat_param(t)%laimin =            &
                laimin(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,     &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
        else
           clsmf25_struc(n)%cat_param(t)%laimin = 0.0
        endif
     enddo
     
     deallocate(laimax)
     deallocate(laimin)

  endif

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

 ! irrigRate = 0.0  
  
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

     gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id

     chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
     if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
          (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
          chhr = chhr -1
     lhr = LIS_rc%hr +chhr
     if(lhr.ge.24) lhr = lhr-24
     if(lhr.lt.0) lhr = lhr+24
     
     ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
      
!SVK: Check 
     lai = LIS_lai(n)%tlai(tid)

! Calculate vegetation and root depth parameters
      
     if((ltime.ge.otimes).and.(ltime.lt.otimee)) then 
        vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
!----------------------------------------------------------------------       
!    Proceed if it is non-forest, non-baresoil, non-urban
!----------------------------------------------------------------------       
      select case ( LIS_rc%lcscheme )
        case( "UMD" )
           veg_index1 = 6
           veg_index2 = 11
        case( "IGBP", "IGBPNCEP", "MODIS" )
           veg_index1 = 6
           veg_index2 = 14
        case( "USGS" )
           veg_index1 = 2
           veg_index2 = 10
        case default
          write(LIS_logunit,*) "The landcover scheme, ",trim(LIS_rc%lcscheme),","
          write(LIS_logunit,*) "is not supported for irrigation. Stopping program ... "
          call LIS_endrun()
      end select

!        if(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col.eq.214.and.&
!             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row.eq.128) then 
!             print*, t, LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
!                  LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
!                  LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row))%lat, &
!                  LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
!                  LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
!                  LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row))%lon
!             stop
!          endif

        if(vegt.ge.veg_index1.and.vegt.le.veg_index2&
             .and.vegt.ne.LIS_rc%bareclass.and.&
             vegt.ne.LIS_rc%urbanclass) then 
           if(irrigFrac(t).gt.0) then 
              ippix = irrigFrac(t)*0.01
              
              smcwlt = clsmf25_struc(n)%cat_param(t)%wpwet* & 
                   clsmf25_struc(n)%cat_param(t)%poros
!              smcref = clsmf25_struc(n)%cat_diagn(t)%ar2*&
!                   clsmf25_struc(n)%cat_param(t)%poros

              smcref = (clsmf25_struc(n)%cat_param(t)%wpwet + & 
                   0.333* (1-clsmf25_struc(n)%cat_param(t)%wpwet))*&
                   clsmf25_struc(n)%cat_param(t)%poros
!              smcref = clsmf25_struc(n)%cat_param(t)%wpwet + & 
!                   0.01* (1-clsmf25_struc(n)%cat_param(t)%wpwet)     

              smcmax = clsmf25_struc(n)%cat_param(t)%poros
	      	                    
             !Determine the amount of irrigation to apply if irrigated tile
              if(ippix.gt.0.0) then !irrigated tile
                 laithresh =   clsmf25_struc(n)%cat_param(t)%laimin + & 
                      0.60 *  (clsmf25_struc(n)%cat_param(t)%laimax - &
                      clsmf25_struc(n)%cat_param(t)%laimin)
                 
                 if(clsmf25_struc(n)%cat_param(t)%laimax.ne.&
                      clsmf25_struc(n)%cat_param(t)%laimin) then 
                    laifac = (lai-clsmf25_struc(n)%cat_param(t)%laimin)/&
                         (clsmf25_struc(n)%cat_param(t)%laimax- &
                         clsmf25_struc(n)%cat_param(t)%laimin)
                 else
                    laifac = 0.0
                 endif
                 crootd = irrigRootdepth(t)*laifac                 
                 if(lai .ge. laithresh.and.crootd.gt.0) then 
		 
		 !!! SPRINKLER IRRIGATION
		    if(LIS_rc%irrigation_type.eq."Sprinkler") then		 
		    !----------------------------------------------------------------------       
		    !    Set the irrigation rate at 6 am; keep the value till next day
		    !    If local time at the tile fall in the irrigation check
		    !    hour then check the root zone average soil moisture
		    !---------------------------------------------------------------------- 
                      if(ltime.eq.otimes) then 	
                         irrigRate(t) = 0.0
		      !-----------------------------------------------------------------------------
		      !     Compute the root zone accumlative soil moisture [mm], 
		      !     field capacity [mm], and wilting point [mm] 
		      !-----------------------------------------------------------------------------
                	  asmc = clsmf25_struc(n)%cat_diagn(t)%rzmc* &
                               crootd*1000.0
                	  tsmcwlt = smcwlt * 1000.0
                	  tsmcref = smcref * 1000.0

		 !-----------------------------------------------------------------------------
		 !     Get the root zone moisture availability to the plant
		 !----------------------------------------------------------------------------- 
                	  if(tsmcref.ge.tsmcwlt) then 
                	     ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt) 
                	  else
                	     ma = -1
                	  endif

                	  if(ma .le.LIS_rc%irrigation_thresh.and. ma.ge.0) then 
                	     twater = tsmcref - asmc
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
		 !                       print*, 'irrig ',t, tsmcref, asmc, irrigRate(t)
		 !                       if(t.eq.35378) then 
		 !                          print*, LIS_rc%da, LIS_rc%hr, asmc, tsmcwlt, tsmcref, irrigRate(t)
                          endif !threshold
		      endif ! otimes	  
		    
		    
		    elseif(LIS_rc%irrigation_type.eq."Flood") then
		      !-----------------------------------------------------------------------------
		      !     Compute the root zone accumlative soil moisture [mm], 
		      !     field capacity [mm], and wilting point [mm] 
		      !-----------------------------------------------------------------------------
                	  asmc = clsmf25_struc(n)%cat_diagn(t)%rzmc* &
                               crootd*1000.0
                	  tsmcwlt = smcwlt * 1000.0
                	  tsmcref = smcref * 1000.0		      
		 !-----------------------------------------------------------------------------
		 !     Get the root zone moisture availability to the plant
		 !----------------------------------------------------------------------------- 
                	  if(tsmcref.ge.tsmcwlt) then 
                	     ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt) 
                	  else
                	     ma = -1
                	  endif

                	  if(ma .le.LIS_rc%irrigation_thresh.and. ma.ge.0) then 		      
		           twater = (smcmax - clsmf25_struc(n)%cat_diagn(t)%sfmc) &
			             *clsmf25_struc(n)%cat_param(t)%dzsf 
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
		 !                       print*, 'irrig ',t, tsmcref, asmc, irrigRate(t)
		 !                       if(t.eq.35378) then 
		 !                          print*, LIS_rc%da, LIS_rc%hr, asmc, tsmcwlt, tsmcref, irrigRate(t)
                           ! clsmf25_struc(n)%cat_diagn(t)%sfmc = smcmax
			   clsmf25_struc(n)%cat_progn(t)%srfexc =  &
			      clsmf25_struc(n)%cat_progn(t)%srfexc + twater
			  endif !threshold		      
		    else
                        ! No Irrigation
                        irrigRate(t) = 0.0   
                    endif ! irrigation type
		    
                 endif
              endif
           end if
        end if
     end if
  enddo
end subroutine clsmf25_getirrigationstates
