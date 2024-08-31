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
! !ROUTINE: cable_driver
! \label{cable_driver}
!
! !REVISION HISTORY:
!  28 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  14 Dec 2011: Claire Carouge (ccc), CABLE LSM improvements
!  20 Aug 2012: Claire Carouge (ccc), CABLE LSM improvements
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
  subroutine cable_driver(n)
! !USES:
    use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_surface
    use LIS_logMod,         only : LIS_logunit, LIS_endrun
    use LIS_timeMgrMod,     only : LIS_date2time, LIS_tick, &
         LIS_isAlarmRinging
    use LIS_constantsMod,  only : LIS_CONST_RHOFW
    use LIS_histDataMod
    use cable_dimensions,   only : ms,msn,ncp,ncs,mf
    use cable_physical_constants, only : sboltz,emleaf,emsoil
    use cable_radiation,   only : sinbet
    use cable_main_module,  only : cable_main
    use cable_arrays
    use cable_lsmMod
    use cable_listime,      only : calc_localtime

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This subroutine is the entry point for calling the CABLE LSM physics.
!  The "cable\_driver" routine is called, which calls additional routines
!  in the proper order to solve the CABLE water/energy/carbon equations.
!  For additional documentation on CABLE, see the webpage at CAWCR:
!      http://www.cawcr.gov.au/projects/access/cable/
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!  \end{description}
!EOP
    integer :: t,l,index1
    integer :: tstep,kstart,kend,nvegt,nsoilt
    integer :: local_hour,change
    integer :: locdoy,locyr,locmo,locda,lochr,locmn,locss
    real*8  :: loctime
    real    :: locgmt
    real    :: dels
    real    :: landmask   ! To output a better landmask when water points included
    logical :: prin
    real    :: soilm_prev, soilm_new
    real    :: eps = 1.0  ! Epsilon to find the water points
    real    :: avail_smoist(ms)
    real    :: avail_sm_mm
    real    :: soilwet
    logical :: alarmCheck
    character*3   :: fnest
  
    write(fnest,'(i3.3)') n
  
    alarmCheck = LIS_isAlarmRinging(LIS_rc,"CABLE model alarm "//trim(fnest))
    if(alarmCheck) then 
       tstep = cable_struc(n)%ktau
    ! kstart and kend are only used in the carbon routine to initialize,
    ! sum up, and then print carbon fluxes.  These sums are not done
    ! anymore, so kstart and kend are not used anymore.  ccc
       kstart = 1
       kend   = 1

       nvegt = cable_struc(n)%nvegt
       nsoilt = cable_struc(n)%nsoilt

    ! temporary - dmm
       rough%za = cable_struc(n)%refheight

       do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          dels = LIS_rc%ts

       ! For water points for now. (ccc)
          if (abs(cable_struc(n)%cable(t)%lai+9999) > eps) then
             
             prin = .false.
             if (((cable_struc(n)%tileprint.eq.t).or.                      &
                  (cable_struc(n)%tileprint.eq.0)).and.                    &
                  (cable_struc(n)%verbose)) prin = .true.

          ! Grid cell index in case subtiling for vegetation.
             index1 = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index

          ! Calculate local time - CABLE uses local time instead of GMT - dmm
             call calc_localtime(LIS_rc%hr,LIS_domain(n)%grid(index1)%lon, &
                  local_hour,change)

             change = change * 3600
             locyr = LIS_rc%yr
             locmo = LIS_rc%mo
             locda = LIS_rc%da
             lochr = LIS_rc%hr
             locmn = LIS_rc%mn
             locss = LIS_rc%ss
             
             call LIS_date2time(loctime,locdoy,locgmt,locyr,locmo,         &
                  locda,lochr,locmn,locss)
             call LIS_tick(loctime,locdoy,locgmt,locyr,locmo,              &
                  locda,lochr,locmn,locss,real(change))
             
             ! Set point attributes (longitude not needed in physics) - dmm
             rad%latitude = LIS_domain(n)%grid(index1)%lat
             
             if (prin) then
                write(LIS_logunit,*) 'Latitude: ',rad%latitude
                write(LIS_logunit,*) 'Longitude: ',LIS_domain(n)%grid(index1)%lon
                write(LIS_logunit,*) 'Canopy Flag: ',model_structure%canopy
                write(LIS_logunit,*) 'Photosynthesis Flag: ',model_structure%photosynthesis
                write(LIS_logunit,*) 'Soil Flag: ',model_structure%soil
                if (any(model_structure%soil.eq.'sli')) then
                   write(LIS_logunit,*) 'sli soil litter Flag: ',model_structure%sli_litter
                   write(LIS_logunit,*) 'sli soil isotope Flag: ',model_structure%sli_isotope
                   write(LIS_logunit,*) 'sli soil coupled Flag: ',model_structure%sli_coupled
                endif
             endif
             
          ! Set the "met_type" values for this timestep and point
             met%year = LIS_rc%yr
             met%moy = LIS_rc%mo
             met%doy = float(locdoy)
             met%hod = float(lochr) + (float(locmn)/60.0) +                &
                  (float(locss)/3600.0)
             met%ca = cable_struc(n)%cable(t)%co2
             met%fsd = cable_struc(n)%cable(t)%swdown
             met%fld = cable_struc(n)%cable(t)%lwdown
             ! Convert units to mm/timestep
             met%precip = (cable_struc(n)%cable(t)%rainf +                 &
                  cable_struc(n)%cable(t)%snowf) * dels
             met%precip_s = cable_struc(n)%cable(t)%snowf * dels
             met%tc = cable_struc(n)%cable(t)%tair - 273.16
             met%tk = cable_struc(n)%cable(t)%tair
             ! tvair and tvrad set to tk before entering physics timestep - dmm
             met%tvair = met%tk
             met%tvrad = met%tk
             met%pmb = cable_struc(n)%cable(t)%psurf
             met%ua = (cable_struc(n)%cable(t)%uwind *                     &
                  cable_struc(n)%cable(t)%uwind) +                    &
                  (cable_struc(n)%cable(t)%vwind *                     &
                  cable_struc(n)%cable(t)%vwind)
             met%ua = sqrt(met%ua)
             met%qv = cable_struc(n)%cable(t)%qair
             ! qvair, da, & dva set in canopy routines - dmm
             met%qvair = 0.0
             met%da = 0.0
             met%dva = 0.0
             ! Set cosine of zenith angle (may be provided by atmospheric model when online):
             met%coszen = sinbet(met%doy, rad%latitude, met%hod)
             
             ! q2sat is given by WRF if coupled
#if (defined COUPLED)
             air%qsat(1) = cable_struc(n)%cable(t)%q2sat
#endif
             
          ! tk_old is only used by "sli" - ccc
             if (cable_struc(n)%soilflag.eq.'sli') then
                if (cable_struc(n)%ktau.eq.1) then
                   met%tk_old = met%tk
                else
                   met%tk_old = cable_struc(n)%cable(t)%tk_old
                endif
             endif
          ! Set solid precip based on temp (EK nov2007)
             met%precip_s = 0.0
             if (met%tc(1).le.0.0) met%precip_s = met%precip
             
             if (prin) then
                write(LIS_logunit,*) 'Year: ',met%year
                write(LIS_logunit,*) 'Month: ',met%moy
                write(LIS_logunit,*) 'Day of year: ',met%doy
                write(LIS_logunit,*) 'Hour: ',met%hod
                write(LIS_logunit,*) 'CO2: ',met%ca
                write(LIS_logunit,*) 'SWdown: ',met%fsd
                write(LIS_logunit,*) 'LWdown: ',met%fld
                write(LIS_logunit,*) 'TotalPrecip: ',met%precip
                write(LIS_logunit,*) 'Snowf: ',met%precip_s
                write(LIS_logunit,*) 'Tc: ',met%tc
                write(LIS_logunit,*) 'Tk: ',met%tk
                write(LIS_logunit,*) 'Pmb: ',met%pmb
                write(LIS_logunit,*) 'Ua: ',met%ua
                write(LIS_logunit,*) 'Qv: ',met%qv
                write(LIS_logunit,*) 'Coszen: ',met%coszen
             endif
             
             ! Set "soil_parameter_type" values for this timestep and point
             soil%isoilm  = cable_struc(n)%cable(t)%soiltype
             soil%albsoil = cable_struc(n)%cable(t)%albsoil
             soil%silt    = cable_struc(n)%cable(t)%silt
             soil%clay    = cable_struc(n)%cable(t)%clay
             soil%sand    = cable_struc(n)%cable(t)%sand
             soil%swilt   = cable_struc(n)%cable(t)%swilt
             soil%sfc     = cable_struc(n)%cable(t)%sfc
             soil%ssat    = cable_struc(n)%cable(t)%ssat
             soil%bch     = cable_struc(n)%cable(t)%bch
             soil%hyds    = cable_struc(n)%cable(t)%hyds
             soil%sucs    = cable_struc(n)%cable(t)%sucs
             soil%rhosoil = cable_struc(n)%cable(t)%rhosoil
             soil%css     = cable_struc(n)%cable(t)%css
             soil%rs20    = cable_struc(n)%cable(t)%rs20
             soil%zse     = cable_struc(n)%cable(t)%zse
             soil%zshh    = cable_struc(n)%cable(t)%zshh
             ! Calculations from subroutine derived_parameters
             soil%cnsd    = soil%sand*0.3 + soil%clay*0.25 + soil%silt*0.265
             soil%hsbh    = soil%hyds*abs(soil%sucs)*soil%bch ! difsat*etasat
             soil%ibp2    = nint(soil%bch)+2
             soil%i2bp3   = 2*nint(soil%bch)+3
             if (any(model_structure%soil=='sli')) then
                ! use 2 soil horizons globally
                ! For now just set B horizon parameters to be the same as A
                soil%nhorizons = 2
                soil%bchB      = soil%bch
                soil%clayB     = soil%clay
                soil%sandB     = soil%sand
                soil%siltB     = soil%silt
                soil%cssB      = soil%css
                soil%hydsB     = soil%hyds
                soil%rhosoilB  = soil%rhosoil
                soil%sfcB      = soil%sfc
                soil%ssatB     = soil%ssat
                soil%sucsB     = soil%sucs
                soil%swiltB    = soil%swilt
                soil%swilt_vec = spread(soil%swilt,2,ms)
                ! Arbitrarily set A horiz depth to be first three layers
                soil%depthA    = sum(soil%zse(1:3))
                soil%depthB    = sum(soil%zse(4:6))
                soil%clitt     = 5.0    ! (tC / ha)
                soil%cnsdB     = soil%sandB*0.3 + soil%clayB*0.25 +        &
                     soil%siltB*0.265
             endif
             if (any(model_structure%canopy=='canopy_vh')) then
                soil%swilt_vec = spread(soil%swilt,2,ms)
             endif
             
             if (prin) then
                write(LIS_logunit,*) 'Soil type: ',soil%isoilm
                write(LIS_logunit,*) 'albsoil: ',soil%albsoil
                write(LIS_logunit,*) 'silt: ',soil%silt
                write(LIS_logunit,*) 'clay: ',soil%clay
                write(LIS_logunit,*) 'sand: ',soil%sand
                write(LIS_logunit,*) 'swilt: ',soil%swilt
                write(LIS_logunit,*) 'sfc: ',soil%sfc
                write(LIS_logunit,*) 'ssat: ',soil%ssat
                write(LIS_logunit,*) 'bch: ',soil%bch
                write(LIS_logunit,*) 'hyds: ',soil%hyds
                write(LIS_logunit,*) 'sucs: ',soil%sucs
                write(LIS_logunit,*) 'rhosoil: ',soil%rhosoil
                write(LIS_logunit,*) 'css: ',soil%css
                write(LIS_logunit,*) 'cnsd: ',soil%cnsd
                write(LIS_logunit,*) 'hsbh: ',soil%hsbh
                write(LIS_logunit,*) 'ibp2: ',soil%ibp2
                write(LIS_logunit,*) 'i2bp3: ',soil%i2bp3
                write(LIS_logunit,*) 'rs20: ',soil%rs20
                write(LIS_logunit,*) 'zse: ',soil%zse
                write(LIS_logunit,*) 'zshh: ',soil%zshh
                if (any(model_structure%soil=='sli')) then
                   write(LIS_logunit,*) 'nhorizons: ',soil%nhorizons
                   write(LIS_logunit,*) 'bchB: ',soil%bchB
                   write(LIS_logunit,*) 'clayB: ',soil%clayB
                   write(LIS_logunit,*) 'sandB: ',soil%sand
                   write(LIS_logunit,*) 'siltB: ',soil%silt
                   write(LIS_logunit,*) 'cssB: ',soil%cssB
                   write(LIS_logunit,*) 'hydsB: ',soil%hydsB
                   write(LIS_logunit,*) 'rhosoilB: ',soil%rhosoilB
                   write(LIS_logunit,*) 'sfcB: ',soil%sfcB
                   write(LIS_logunit,*) 'ssatB: ',soil%ssatB
                   write(LIS_logunit,*) 'sucsB: ',soil%sucsB
                   write(LIS_logunit,*) 'cnsdB: ',soil%cnsdB
                   write(LIS_logunit,*) 'swiltB: ',soil%swiltB
                   write(LIS_logunit,*) 'swilt_vec: ',soil%swilt_vec
                   write(LIS_logunit,*) 'depthA: ',soil%depthA
                   write(LIS_logunit,*) 'depthB: ',soil%depthB
                   write(LIS_logunit,*) 'clitt: ',soil%clitt
                endif
             endif
             
             ! Set "veg_parameter_type" values for this timestep and point
             veg%iveg       = cable_struc(n)%cable(t)%vegtype
             veg%canst1     = cable_struc(n)%cable(t)%canst1
             veg%dleaf      = cable_struc(n)%cable(t)%dleaf
             veg%ejmax      = cable_struc(n)%cable(t)%ejmax
             veg%extkn      = cable_struc(n)%cable(t)%extkn
             veg%hc         = cable_struc(n)%cable(t)%hc
             veg%rp20       = cable_struc(n)%cable(t)%rp20
             veg%rpcoef     = cable_struc(n)%cable(t)%rpcoef
             veg%shelrb     = cable_struc(n)%cable(t)%shelrb
             veg%tmaxvj     = cable_struc(n)%cable(t)%tmaxvj
             veg%tminvj     = cable_struc(n)%cable(t)%tminvj
             veg%vbeta      = cable_struc(n)%cable(t)%vbeta
             veg%vegcf      = cable_struc(n)%cable(t)%vegcf
             veg%vcmax      = cable_struc(n)%cable(t)%vcmax
             veg%wai        = cable_struc(n)%cable(t)%wai
             veg%xfang      = cable_struc(n)%cable(t)%xfang
             veg%vlai       = cable_struc(n)%cable(t)%lai
             veg%frac4      = cable_struc(n)%cable(t)%frac4
             veg%froot(1,:) = cable_struc(n)%cable(t)%froot(:)
             if (any(model_structure%soil=='sli')) then
                ! New veg parameters:
                veg%F10     = 0.85
                veg%ZR      = 5.0
                veg%gamma   = 1.e-2
             endif
             if (any(model_structure%canopy=='canopy_vh')) then
                veg%gamma   = 1.e-2
                veg%gamma   = 1.0e-5
             endif
             veg%deciduous = .false.
             ! Check this code against CASA/UMD 13 vegetation types - dmm
             if (LIS_rc%lcscheme.eq."UMD") then
                where((veg%iveg.eq.2).or.(veg%iveg.eq.5))                  &
                     veg%deciduous = .true.
                !check if this true. Earlier it was set to 6
             elseif (LIS_rc%lcscheme.eq."MODIS") then
                where((veg%iveg.eq.3).or.(veg%iveg.eq.4))                  &
                     veg%deciduous = .true.
             endif
             veg%meth = cable_struc(n)%cable(t)%meth
             
             if (prin) then
                write(LIS_logunit,*) 'Vegetation type: ',veg%iveg
                write(LIS_logunit,*) 'canst1: ',veg%canst1
                write(LIS_logunit,*) 'dleaf: ',veg%dleaf
                write(LIS_logunit,*) 'ejmax: ',veg%ejmax
                write(LIS_logunit,*) 'extkn: ',veg%extkn
                write(LIS_logunit,*) 'hc: ',veg%hc
                write(LIS_logunit,*) 'rp20: ',veg%rp20
                write(LIS_logunit,*) 'rpcoef: ',veg%rpcoef
                write(LIS_logunit,*) 'shelrb: ',veg%shelrb
                write(LIS_logunit,*) 'tmaxvj: ',veg%tmaxvj
                write(LIS_logunit,*) 'tminvj: ',veg%tminvj
                write(LIS_logunit,*) 'vbeta: ',veg%vbeta
                write(LIS_logunit,*) 'vegcf: ',veg%vegcf
                write(LIS_logunit,*) 'vcmax: ',veg%vcmax
                write(LIS_logunit,*) 'wai: ',veg%wai
                write(LIS_logunit,*) 'xfang: ',veg%xfang
                write(LIS_logunit,*) 'vlai: ',veg%vlai
                write(LIS_logunit,*) 'frac4: ',veg%frac4
                write(LIS_logunit,*) 'froot: ',veg%froot
                if (any(model_structure%soil=='sli')) then
                   write(LIS_logunit,*) 'F10: ',veg%F10
                   write(LIS_logunit,*) 'ZR: ',veg%ZR
                   write(LIS_logunit,*) 'gamma: ',veg%gamma
                endif
                if (any(model_structure%canopy=='canopy_vh')) then
                   write(LIS_logunit,*) 'gamma: ',veg%gamma
                endif
                write(LIS_logunit,*) 'deciduous: ',veg%deciduous
                write(LIS_logunit,*) 'meth: ',veg%meth
             endif
             
             ! Set "canopy_type" values for this timestep and point
             canopy%cansto = cable_struc(n)%cable(t)%cansto
             
             ! Set "soil_snow_type" values for this timestep and point
             ssoil%albsoilsn(1,:) = cable_struc(n)%cable(t)%albsoilsn(:)
             ssoil%isflag = cable_struc(n)%cable(t)%isflag
             ssoil%osnowd = cable_struc(n)%cable(t)%osnowd
             ssoil%rtsoil = cable_struc(n)%cable(t)%rtsoil
             ssoil%smass(1,:) = cable_struc(n)%cable(t)%smass(:)
             ssoil%snage = cable_struc(n)%cable(t)%snage
             ssoil%snowd = cable_struc(n)%cable(t)%snowd
             ssoil%ssdn(1,:) = cable_struc(n)%cable(t)%ssdn(:)
             ssoil%ssdnn = cable_struc(n)%cable(t)%ssdnn
             ssoil%tgg(1,:) = cable_struc(n)%cable(t)%tgg(:)
             ssoil%tggsn(1,:) = cable_struc(n)%cable(t)%tggsn(:)
             ssoil%wb(1,:) = cable_struc(n)%cable(t)%wb(:)
             ssoil%wbice(1,:) = cable_struc(n)%cable(t)%wbice(:)
             !ccc
             ssoil%pwb_min(1) = cable_struc(n)%cable(t)%pwb_min
             ssoil%gammzz(1,:) = cable_struc(n)%cable(t)%gammzz
             ssoil%wbtot(1) = cable_struc(n)%cable(t)%wbtot
             
             if (prin) then
                write(LIS_logunit,*) 'cansto: ',canopy%cansto
                write(LIS_logunit,*) 'albsoilsn: ',ssoil%albsoilsn
                write(LIS_logunit,*) 'isflag: ',ssoil%isflag
                write(LIS_logunit,*) 'osnowd: ',ssoil%osnowd
                write(LIS_logunit,*) 'rtsoil: ',ssoil%rtsoil
                write(LIS_logunit,*) 'smass: ',ssoil%smass
                write(LIS_logunit,*) 'snage: ',ssoil%snage
                write(LIS_logunit,*) 'snowd: ',ssoil%snowd
                write(LIS_logunit,*) 'ssdn: ',ssoil%ssdn
                write(LIS_logunit,*) 'ssdnn: ',ssoil%ssdnn
                write(LIS_logunit,*) 'tgg: ',ssoil%tgg
                write(LIS_logunit,*) 'tggsn: ',ssoil%tggsn
                write(LIS_logunit,*) 'wb: ',ssoil%wb
                write(LIS_logunit,*) 'wbice: ',ssoil%wbice
                !ccc
                write(LIS_logunit,*) 'pwb_min: ',ssoil%pwb_min
                write(LIS_logunit,*) 'gammzz: ',ssoil%gammzz
                write(LIS_logunit,*) 'wbtot:  ',ssoil%wbtot
             endif
             
             ! Set "bge_pool_type" values for this timestep and point
             bgc%cplant(1,:) = cable_struc(n)%cable(t)%cplant(:)
             bgc%csoil(1,:) = cable_struc(n)%cable(t)%csoil(:)
            
             soilm_prev =  0 
             do l=1,ms
                soilm_prev = soilm_prev + & 
                     cable_struc(n)%cable(t)%zse(l)*LIS_CONST_RHOFW*& 
                     cable_struc(n)%cable(t)%wb(l)
             enddo
             ! CALL land surface scheme for this timestep, all grid points:
             ! NOTE - calling "cable_main", but all grid points here is
             !        defined as a single tile in the LIS structure - dmm
             call cable_main(tstep, kstart, kend, dels, air, bgc, canopy,  &
                  met, bal, rad, rough, soil, ssoil, sum_flux, veg, nvegt,    &
                  nsoilt, model_structure, prin)
             
             ! After the call to cable_main, return the state variables at this
             ! point to the LIS arrays for output and the next timestep. - dmm
             cable_struc(n)%cable(t)%cansto = canopy%cansto(1)
             cable_struc(n)%cable(t)%albsoilsn(:) = ssoil%albsoilsn(1,:)
             cable_struc(n)%cable(t)%isflag = ssoil%isflag(1)
             cable_struc(n)%cable(t)%osnowd = ssoil%osnowd(1)
             cable_struc(n)%cable(t)%rtsoil = ssoil%rtsoil(1)
             cable_struc(n)%cable(t)%smass(:) = ssoil%smass(1,:)
             cable_struc(n)%cable(t)%snage = ssoil%snage(1)
             cable_struc(n)%cable(t)%snowd = ssoil%snowd(1)
             cable_struc(n)%cable(t)%ssdn(:) = ssoil%ssdn(1,:)
             cable_struc(n)%cable(t)%tgg(:) = ssoil%tgg(1,:)
             cable_struc(n)%cable(t)%tggsn(:) = ssoil%tggsn(1,:)
             cable_struc(n)%cable(t)%wb(:) = ssoil%wb(1,:)
             cable_struc(n)%cable(t)%wbice(:) = ssoil%wbice(1,:)
             cable_struc(n)%cable(t)%cplant(:) = bgc%cplant(1,:)
             cable_struc(n)%cable(t)%csoil(:) = bgc%csoil(1,:)
             !ccc
             cable_struc(n)%cable(t)%pwb_min = ssoil%pwb_min(1)
             cable_struc(n)%cable(t)%gammzz  = ssoil%gammzz(1,:)
             cable_struc(n)%cable(t)%wbtot   = ssoil%wbtot(1)
             
             ! Added for coupling to WRF (ccc)
             cable_struc(n)%cable(t)%z0m = rough%z0m(1)
             cable_struc(n)%cable(t)%smelt = ssoil%smelt(1)
             cable_struc(n)%cable(t)%trad = rad%trad(1)
             
             ! Prepare diagnostics values (ccc)
             cable_struc(n)%cable(t)%swnet = sum(rad%qcan(1,:,1),1) +      &
                  sum(rad%qcan(1,:,2),1) +      &
                  rad%qssabs(1)
             cable_struc(n)%cable(t)%lwnet = met%fld(1) - sboltz*emleaf*   &
                  canopy%tv(1)**4*(1-rad%transd(1)) - &
                  rad%flws(1)*rad%transd(1)
             cable_struc(n)%cable(t)%rnet = cable_struc(n)%cable(t)%swnet+ &
                  cable_struc(n)%cable(t)%lwnet
             cable_struc(n)%cable(t)%qle = canopy%fe(1)
             cable_struc(n)%cable(t)%qh = canopy%fh(1)
             cable_struc(n)%cable(t)%qg = canopy%ga(1)
             cable_struc(n)%cable(t)%qs = ssoil%rnof1(1)/dels
             cable_struc(n)%cable(t)%qsb = ssoil%rnof2(1)/dels
             cable_struc(n)%cable(t)%snowf = met%precip_s(1)/dels
             cable_struc(n)%cable(t)%rainf = met%precip(1)/dels
             cable_struc(n)%cable(t)%evap = canopy%fe(1)/air%rlam(1)
             cable_struc(n)%cable(t)%vegt = canopy%tv(1)
             cable_struc(n)%cable(t)%radt = (((1.0-rad%transd(1))*emleaf*  &
                  sboltz*canopy%tv(1)**4         &
                  + rad%transd(1)*emsoil*sboltz*   &
                  ((1-ssoil%isflag(1))*ssoil%tgg(1,1) &
                  + ssoil%isflag(1)*ssoil%tggsn(1,1)) ** 4) &
                  / sboltz)**0.25
             cable_struc(n)%cable(t)%albedo = 0.5*(rad%albedo(1,1)+rad%albedo(1,2))
             cable_struc(n)%cable(t)%sdepth = SUM(ssoil%sdepth(1,:),1)
             cable_struc(n)%cable(t)%ecanop = canopy%fevw(1)/air%rlam(1)
             cable_struc(n)%cable(t)%tveg = canopy%fevc(1) / air%rlam(1)
             cable_struc(n)%cable(t)%esoil = canopy%fes(1)/air%rlam(1)
             cable_struc(n)%cable(t)%autoresp=(canopy%frp(1)+canopy%frday(1))/  &
                  12.01E-6
             cable_struc(n)%cable(t)%heteroresp = canopy%frs(1)/12.01E-6
             cable_struc(n)%cable(t)%leafresp = canopy%frday(1)/12.01E-6
             cable_struc(n)%cable(t)%npp = (-1.0*canopy%fpn(1)-canopy%frp(1))/ &
                  12.01E-6
             cable_struc(n)%cable(t)%gpp = (canopy%frday(1)-canopy%fpn(1))/ &
                  12.01E-6
             cable_struc(n)%cable(t)%nee = canopy%fnee(1)/12.01E-6

             !Jatin Added output for Potential Evapotranspiration, potev
             cable_struc(n)%cable(t)%potev = ssoil%potev(1)
             landmask = 1.0

          else ! To have diagnostics for water points

             cable_struc(n)%cable(t)%swnet = -9999
             cable_struc(n)%cable(t)%lwnet = -9999
             cable_struc(n)%cable(t)%rnet = -9999
             cable_struc(n)%cable(t)%qle = -9999
             cable_struc(n)%cable(t)%qh = -9999
             cable_struc(n)%cable(t)%qg = -9999
             cable_struc(n)%cable(t)%qs = -9999
             cable_struc(n)%cable(t)%qsb = -9999
             cable_struc(n)%cable(t)%snowf = -9999
             cable_struc(n)%cable(t)%rainf = -9999
             cable_struc(n)%cable(t)%evap = -9999
             cable_struc(n)%cable(t)%tggsn(1) = -9999
             cable_struc(n)%cable(t)%vegt = -9999
             cable_struc(n)%cable(t)%tgg(1) = -9999
             cable_struc(n)%cable(t)%radt = -9999
             cable_struc(n)%cable(t)%albedo = -9999
             cable_struc(n)%cable(t)%snowd = -9999
             cable_struc(n)%cable(t)%sdepth = -9999
             cable_struc(n)%cable(t)%vegt = -9999
             do l = 1,ms
                cable_struc(n)%cable(t)%tgg(l) = -9999
                cable_struc(n)%cable(t)%wb(l) = -9999
                cable_struc(n)%cable(t)%wbice(l) = -9999
             enddo
             cable_struc(n)%cable(t)%cansto = -9999
             cable_struc(n)%cable(t)%ecanop = -9999
             cable_struc(n)%cable(t)%tveg = -9999
             cable_struc(n)%cable(t)%esoil = -9999
             cable_struc(n)%cable(t)%autoresp = -9999
             cable_struc(n)%cable(t)%heteroresp = -9999
             cable_struc(n)%cable(t)%leafresp = -9999
             cable_struc(n)%cable(t)%npp = -9999
             cable_struc(n)%cable(t)%gpp = -9999
             cable_struc(n)%cable(t)%nee = -9999
             !Jatin
             cable_struc(n)%cable(t)%potev = 0
             landmask = 0.0

          endif

          ! General energy balance components output
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWNET,value=           &
               cable_struc(n)%cable(t)%swnet,                     &
               vlevel=1,unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LWNET,value=           &
               cable_struc(n)%cable(t)%lwnet,                     &
               vlevel=1,unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RNET,value=            &
               cable_struc(n)%cable(t)%rnet,                      &
               vlevel=1,unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QLE,value=             &
               cable_struc(n)%cable(t)%qle,                       &
               vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QH,value=              &
               cable_struc(n)%cable(t)%qh,                        &
               vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QG,value=              &
               cable_struc(n)%cable(t)%qg,                        &
               vlevel=1,unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)

          ! Water balance output
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QS,value=              &
               cable_struc(n)%cable(t)%qs,                        &
               vlevel=1,unit="kg m-2 s-1",direction="OUT",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSB,value=             &
               cable_struc(n)%cable(t)%qsb,                       &
               vlevel=1,unit="kg m-2 s-1",direction="OUT",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=           &
               cable_struc(n)%cable(t)%snowf,                     &
               vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=           &
               cable_struc(n)%cable(t)%rainf,                     &
               vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EVAP,value=            &
               cable_struc(n)%cable(t)%evap,                      &
               vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
          soilm_new = 0 
          do l=1,ms
             soilm_new = soilm_new + & 
                  cable_struc(n)%cable(t)%zse(l)*LIS_CONST_RHOFW*& 
                  cable_struc(n)%cable(t)%wb(l)
          enddo
          
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSOILMOIST, &
               value=(soilm_new - soilm_prev),vlevel=1,unit="kg m-2",&
               direction="INC",surface_type=LIS_rc%lsm_index)
          ! Surface state variables output
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWT,value=           &
               cable_struc(n)%cable(t)%tggsn(1),                  &
               vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_VEGT,value=            &
               cable_struc(n)%cable(t)%vegt,                      &
               vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_BARESOILT,value=       &
               cable_struc(n)%cable(t)%tgg(1),                    &
               vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RADT,value=            &
               cable_struc(n)%cable(t)%radt,                      &
               vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ALBEDO,value=          &
               cable_struc(n)%cable(t)%albedo,                    &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE,value=             &
               cable_struc(n)%cable(t)%snowd,                     &
               vlevel=1,unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWDEPTH,value=       &
               cable_struc(n)%cable(t)%sdepth,                    &
               vlevel=1,unit="m",direction="-",surface_type=LIS_rc%lsm_index)

          ! Subsurface state variables output
          do l = 1,ms
             call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILTEMP,value=     &
                  cable_struc(n)%cable(t)%tgg(l),                 &
                  vlevel=l,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
             call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,value=    &
                  cable_struc(n)%cable(t)%wb(l),                  &
                  vlevel=l,unit="m^3 m-3",direction="-",&
                  surface_type=LIS_rc%lsm_index)
             call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,value= &
                  cable_struc(n)%cable(t)%wb(l)*LIS_CONST_RHOFW*&
                  cable_struc(n)%cable(t)%zse(l),                  &
                  vlevel=l,unit="kg m-2",direction="-",&
                  surface_type=LIS_rc%lsm_index)
             ! Putting "wbice" into the ALMA "SmFrozFrac" variable for now
             ! This required a change in LIS_histDataMod to add this option
             ! Can output in m^3 m-3 or fraction of the total - dmm
             call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SMFROZFRAC,value=   &
                  cable_struc(n)%cable(t)%wbice(l),               &
                  vlevel=l,unit="m^3 m-3",direction="-",surface_type=LIS_rc%lsm_index)
             call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SMFROZFRAC,value=   &
                  (cable_struc(n)%cable(t)%wbice(l) /             &
                  cable_struc(n)%cable(t)%wb(l)),                 &
                  vlevel=l,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          enddo

!soil wetness
          do l=1,ms
             avail_smoist(l) = max(real(cable_struc(n)%cable(t)%wb(l) &
                  -cable_struc(n)%cable(t)%wbice(l)) - &
                  cable_struc(n)%cable(t)%swilt, 0.0)
             avail_sm_mm = avail_sm_mm + avail_smoist(l)*1000.0*&
                  cable_struc(n)%cable(t)%zse(l)
          enddo
          soilwet = avail_sm_mm/((cable_struc(n)%cable(t)%ssat-&
               cable_struc(n)%cable(t)%swilt)*&
               SUM(cable_struc(n)%cable(t)%zse(:)))
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILWET,value=soilwet,   &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)

          ! Evaporation components output
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CANOPINT,value=        &
               cable_struc(n)%cable(t)%cansto,                    &
               vlevel=1,unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,value=          &
               cable_struc(n)%cable(t)%ecanop,                    &
               vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,value=            &
               cable_struc(n)%cable(t)%tveg,                      &
               vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,value=           &
               cable_struc(n)%cable(t)%esoil,                     &
               vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
! Jatin adding call for potev 
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,value=         &
              cable_struc(n)%cable(t)%potev,                     &
               vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)

          ! Carbon output
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_AUTORESP,value=        &
               cable_struc(n)%cable(t)%autoresp,                  &
               vlevel=1,unit="umol m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_HETERORESP,value=      &
               cable_struc(n)%cable(t)%heteroresp,                &
               vlevel=1,unit="umol m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LEAFRESP,value=        &
               cable_struc(n)%cable(t)%leafresp,                  &
               vlevel=1,unit="umol m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_NPP,value=             &
               cable_struc(n)%cable(t)%npp,                       &
               vlevel=1,unit="umol m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GPP,value=             &
               cable_struc(n)%cable(t)%gpp,                       &
               vlevel=1,unit="umol m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_NEE,value=             &
               cable_struc(n)%cable(t)%nee,                       &
               vlevel=1,unit="umol m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)

          ! Parameter output
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LANDMASK,value=        &
               landmask,                                          &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LANDCOVER,value=       &
               float(cable_struc(n)%cable(t)%vegtype),            &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILTYPE,value=        &
               float(cable_struc(n)%cable(t)%soiltype),           &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SANDFRAC,value=        &
               cable_struc(n)%cable(t)%sand,                      &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SILTFRAC,value=        &
               cable_struc(n)%cable(t)%silt,                      &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CLAYFRAC,value=        &
               cable_struc(n)%cable(t)%clay,                      &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POROSITY,value=        &
               cable_struc(n)%cable(t)%ssat,                      &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LAI,value=             &
               cable_struc(n)%cable(t)%lai,                       &
               vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)

       enddo

       cable_struc(n)%ktau = cable_struc(n)%ktau + 1
    end if
  end subroutine cable_driver
