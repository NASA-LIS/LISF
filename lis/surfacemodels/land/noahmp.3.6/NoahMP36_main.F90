!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: NoahMP36_main
! \label{NoahMP36_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   9/4/14: Shugong Wang; initial implementation for NoahMP36 with LIS-7
!   2/7/18: Soni Yatheendradas; code added for OPTUE to work
!
! !INTERFACE:
subroutine NoahMP36_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_constantsMod,  only : LIS_CONST_RHOFW
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use NoahMP36_lsmMod
   !use other modules
    use ESMF
    use LIS_routingMod, only : LIS_runoff_state
  
    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    real                 :: dt
    real                 :: lat, lon
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    logical              :: alarmCheck

!
! !DESCRIPTION:
!  This is the entry point for calling the NoahMP36 physics.
!  This routine calls the {\tt noahmp\_driver\_36 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for NoahMP36
    character(len=256)   :: tmp_landuse_tbl_name   ! Noah model landuse parameter table [-]
    character(len=256)   :: tmp_soil_tbl_name      ! Noah model soil parameter table [-]
    character(len=256)   :: tmp_gen_tbl_name       ! Noah model general parameter table [-]
    character(len=256)   :: tmp_noahmp_tbl_name    ! NoahMP parameter table [-]
    character(len=256)   :: tmp_landuse_scheme_name ! landuse classification scheme [-]
    character(len=256)   :: tmp_soil_scheme_name   ! soil classification scheme [-]
    integer              :: tmp_dveg_opt           ! vegetation model [-]
    integer              :: tmp_crs_opt            ! canopy stomatal resistance [-]
    integer              :: tmp_btr_opt            ! soil moisture factor for stomatal resistance [-]
    integer              :: tmp_run_opt            ! runoff and groundwater [-]
    integer              :: tmp_sfc_opt            ! surface layer drag coefficients (CH & CM) [-]
    integer              :: tmp_frz_opt            ! supercooled liquid water [-]
    integer              :: tmp_inf_opt            ! frozen soil permeability [-]
    integer              :: tmp_rad_opt            ! radiation transfer [-]
    integer              :: tmp_alb_opt            ! snow surface albedo [-]
    integer              :: tmp_snf_opt            ! rainfall & snowfall [-]
    integer              :: tmp_tbot_opt           ! lower boundary of soil temperature [-]
    integer              :: tmp_stc_opt            ! snow/soil temperature time scheme [-]
    integer              :: tmp_nslcats            ! the number of total soil types in parameter table [-]
    integer              :: tmp_nlucats            ! the number of total land cover types in parameter table [-]
    integer              :: tmp_nslpcats           ! the number of total slope category for Noah baseflow [-]
    real                 :: tmp_latitude           ! latitude in decimal degree [-]
    real                 :: tmp_longitude          ! longitude in decimal degree [-]
    integer              :: tmp_year               ! year of the current time step [-]
    integer              :: tmp_month              ! month of the current time step [-]
    integer              :: tmp_day                ! day of the current time step [-]
    integer              :: tmp_hour               ! hour of the current time step [-]
    integer              :: tmp_minute             ! minute of the current time step [-]
    real                 :: tmp_dt                 ! time step in seconds [s]
    real                 :: tmp_albd(2)
    real                 :: tmp_albi(2)
    integer              :: tmp_nsoil              ! number of soil layers [-]
    real, allocatable    :: tmp_sldpth(:)          ! thickness of soil layers [-]
    integer              :: tmp_nsnow              ! maximum number of snow layers [-]
    real, allocatable    :: tmp_shdfac_monthly(:)  ! monthly values for green vegetation fraction [-]
    integer              :: tmp_vegetype           ! land cover type index [-]
    integer              :: tmp_soiltype           ! soil type index [-]
    integer              :: tmp_slopetype          ! slope type for Noah baseflow [-]
    integer              :: tmp_urban_vegetype     ! urban land cover type index [-]
    integer              :: tmp_ice_flag           ! ice flag: 0 = no ice, 1 = ice [-]
    integer              :: tmp_st_flag            ! surface type 1=soil, 2=lake [-]
    integer              :: tmp_sc_idx             ! soil color type [-]
    integer              :: tmp_iz0tlnd            ! option of Chen adjustment of Czil [-]
    real, allocatable    :: tmp_smceq(:)           ! equilibrium soil water content [m^3 m-3]
    real                 :: tmp_tair               ! air temperature [K]
    real                 :: tmp_psurf              ! air pressure [Pa]
    real                 :: tmp_wind_e             ! eastward wind speed [m s-1]
    real                 :: tmp_wind_n             ! northward wind speed [m s-1]
    real                 :: tmp_qair               ! near Surface Specific Humidity [kg kg-1]
    real                 :: tmp_swdown             ! downward solar radiation [w/m2]
    real                 :: tmp_lwdown             ! downward longwave radiation [w/m2]
    real                 :: tmp_prcp               ! total precip (rainfall+snowfall) Rate [kg m-2 s-1]   
    !real                 :: tmp_rainf              ! rainfall Rate [kg m-2 s-1]
    !real                 :: tmp_snowf              ! snowfall Rate [kg m-2 s-1]
    real                 :: tmp_tbot               ! deep-layer soil temperature [K]
    real                 :: tmp_pblh               ! planetary boundary layer height [m]
    real                 :: tmp_zlvl               ! reference height of temperature and humidity [m]
    real                 :: tmp_albold             ! snow albedo at last time step [-]
    real                 :: tmp_sneqvo             ! snow mass at the last time step [mm]
    real, allocatable    :: tmp_sstc(:)            ! snow/soil temperature [K]
    real, allocatable    :: tmp_sh2o(:)            ! volumetric liquid soil moisture [m^3 m-3]
    real, allocatable    :: tmp_smc(:)             ! volumetric soil moisture, ice + liquid [m^3 m-3]
    real                 :: tmp_tah                ! canopy air temperature [K]
    real                 :: tmp_eah                ! canopy air vapor pressure [Pa]
    real                 :: tmp_fwet               ! wetted or snowed fraction of canopy [-]
    real                 :: tmp_canliq             ! intercepted liquid water [mm]
    real                 :: tmp_canice             ! intercepted ice mass [mm]
    real                 :: tmp_tv                 ! vegetation temperature [K]
    real                 :: tmp_tg                 ! ground temperature [K]
    real                 :: tmp_qsnow              ! snowfall on the ground [mm s-1]
    integer              :: tmp_isnow              ! actual number of snow layers [-]
    real, allocatable    :: tmp_zss(:)             ! snow/soil layer-bottom depth from snow surface [m]
    real                 :: tmp_snowh              ! snow height [m]
    real                 :: tmp_sneqv              ! snow water equivalent [mm]
    real, allocatable    :: tmp_snowice(:)         ! snow-layer ice [mm]
    real, allocatable    :: tmp_snowliq(:)         ! snow-layer liquid water [mm]
    real                 :: tmp_zwt                ! depth to water table [m]
    real                 :: tmp_wa                 ! water storage in aquifer [mm]
    real                 :: tmp_wt                 ! water in aquifer and saturated soil [mm]
    real                 :: tmp_wslake             ! lake water storage [mm]
    real                 :: tmp_lfmass             ! leaf mass [g/m2]
    real                 :: tmp_rtmass             ! mass of fine roots [g/m2]
    real                 :: tmp_stmass             ! stem mass [g/m2]
    real                 :: tmp_wood               ! mass of wood including woody roots [g/m2]
    real                 :: tmp_stblcp             ! stable carbon in deep soil [g/m2]
    real                 :: tmp_fastcp             ! short-lived carbon in shallow soil [g/m2]
    real                 :: tmp_lai                ! leaf area index [-]
    real                 :: tmp_sai                ! stem area index [-]
    real                 :: tmp_cm                 ! momentum drag coefficient [m s-1]
    real                 :: tmp_ch                 ! sensible heat exchange coefficient [m s-1]
    real                 :: tmp_tauss              ! snow aging term [-]
    real                 :: tmp_smcwtd             ! soil water content between bottom of the soil and water table [m^3 m-3]
    real                 :: tmp_deeprech           ! recharge to or from the water table when deep [m]
    real                 :: tmp_rech               ! recharge to or from the water table when shallow [m]
    real                 :: tmp_fsa                ! total absorbed solar radiation [W m-2]
    real                 :: tmp_fsr                ! total reflected solar radiation [W m-2]
    real                 :: tmp_fira               ! total net longwave radiation to atmosphere [W m-2]
    real                 :: tmp_fsh                ! total sensible heat to atmosphere [W m-2]
    real                 :: tmp_ssoil              ! ground heat flux to soil [W m-2]
    real                 :: tmp_fcev               ! canopy evaporative heat to atmosphere [W m-2]
    real                 :: tmp_fgev               ! ground evaporative heat to atmosphere [W m-2]
    real                 :: tmp_fctr               ! transpiration heat to atmosphere [W m-2]
    real                 :: tmp_ecan               ! evaporation rate of canopy water [kg m-2 s-1]
    real                 :: tmp_etran              ! transpiration rate [kg m-2 s-1]
    real                 :: tmp_edir               ! direct evaporation rate from surface [kg m-2 s-1]
    real                 :: tmp_trad               ! surface radiative temperature [K]
    real                 :: tmp_subsnow            ! snow sublimation rate [kg m-2 s-1]
    real                 :: tmp_tgb                ! ground temperature [K]
    real                 :: tmp_tgv                ! ground surface temperature [K]
    real                 :: tmp_t2mv               ! 2-m air temperature over vegetated part [K]
    real                 :: tmp_t2mb               ! 2-m height air temperature [K]
    real                 :: tmp_q2v                ! 2-m specific humidity over vegetation [kg kg-1]
    real                 :: tmp_q2b                ! 2-m air specific humidity [kg kg-1]
    real                 :: tmp_runsrf             ! surface runoff [kg m-2 s-1]
    real                 :: tmp_runsub             ! baseflow (saturation excess) [kg m-2 s-1]
    real                 :: tmp_apar               ! photosynthesis active energy by canopy [W m-2]
    real                 :: tmp_psn                ! total photosynthesis of CO2 [umol m-2 s-1]
    real                 :: tmp_sav                ! solar radiation absorbed by vegetation [W m-2]
    real                 :: tmp_sag                ! solar radiation absorbed by ground [W m-2]
    real                 :: tmp_fsno               ! snow-cover fraction on the ground [-]
    real                 :: tmp_nee                ! net ecosystem exchange of CO2 [g/m2s]
    real                 :: tmp_gpp                ! net instantaneous assimilation of carbon [g/m2s]
    real                 :: tmp_npp                ! net primary productivity of carbon [g/m2s]
    real                 :: tmp_fveg               ! green vegetation fraction [-]
    real                 :: tmp_albedo             ! surface albedo [-]
    real                 :: tmp_qsnbot             ! melting water out of snow bottom [kg m-2 s-1]
    real                 :: tmp_ponding            ! surface ponding [mm]
    real                 :: tmp_ponding1           ! surface ponding1 [mm]
    real                 :: tmp_ponding2           ! surface ponding2 [mm]
    real                 :: tmp_rssun              ! sunlit stomatal resistance [s/m]
    real                 :: tmp_rssha              ! shaded stomatal resistance [s/m]
    real                 :: tmp_bgap               ! between canopy gap fraction for beam [-]
    real                 :: tmp_wgap               ! within canopy gap fraction for beam [-]
    real                 :: tmp_chv                ! sensible heat exchange coefficient over vegetated fraction [m s-1]
    real                 :: tmp_chb                ! sensible heat exchange coefficient over bare-ground fraction [m s-1]
    real                 :: tmp_emissi             ! surface emissivity [-]
    real                 :: tmp_shg                ! ground sensible heat [W m-2]
    real                 :: tmp_shc                ! canopy sensible heat [W m-2]
    real                 :: tmp_shb                ! bare ground sensible heat [W m-2]
    real                 :: tmp_evg                ! ground evaporation heat [W m-2]
    real                 :: tmp_evb                ! bare ground evaporation heat [W m-2]
    real                 :: tmp_ghv                ! ground heat flux [W m-2]
    real                 :: tmp_ghb                ! bare ground heat flux [W m-2]
    real                 :: tmp_irg                ! ground net long wave radiation [W m-2]
    real                 :: tmp_irc                ! canopy net long wave radiation [W m-2]
    real                 :: tmp_irb                ! bare ground net long wave radiation [W m-2]
    real                 :: tmp_tr                 ! transpiration heat [W m-2]
    real                 :: tmp_evc                ! canopy evaporation heat [W m-2]
    real                 :: tmp_chleaf             ! leaf exchange coefficient [-]
    real                 :: tmp_chuc               ! under canopy exchange coefficient [-]
    real                 :: tmp_chv2               ! sensible heat exchange coefficient over vegetated fraction [-]
    real                 :: tmp_chb2               ! sensible heat exchange coefficient over bare-ground [-]
    real                 :: tmp_fpice              ! snow fraction in precipitation [-]
    real                 :: tmp_sfcheadrt          ! extra output for WRF-HYDRO [m]
    
    ! SY: Begin for enabling OPTUE
    ! SY: Begin corresponding to REDPRM
    ! SY: Begin SOIL PARAMETERS 
    real                 :: tmp_csoil          ! vol. soil heat capacity [j/m3/K] 
    real                 :: tmp_bexp           ! B parameter 
    real                 :: tmp_dksat          ! saturated soil hydraulic conductivity
    real                 :: tmp_dwsat          ! saturated soil hydraulic diffusivity
    real                 :: tmp_psisat         ! saturated soil matric potential 
    real                 :: tmp_quartz         ! soil quartz content 
    real                 :: tmp_smcmax         ! porosity, saturated value of soil moisture (volumetric) 
    real                 :: tmp_smcref         ! reference soil moisture (field capacity) 
    real                 :: tmp_smcwlt         ! wilting point soil moisture (volumetric) 
    ! SY: End SOIL PARAMETERS 
    ! SY: Begin UNIVERSAL PARAMETERS
    real                 :: tmp_czil           ! Calculate roughness length of heat 
    real                 :: tmp_frzk           ! frozen ground parameter 
    real                 :: tmp_refdk          ! parameters in the surface runoff parameteriz.
    real                 :: tmp_refkdt         ! parameters in the surface runoff parameteriz. 
    real                 :: tmp_slope          ! slope index (0 - 1)
    ! SY: End UNIVERSAL PARAMETERS
    ! SY: Begin PARAMETERS for prescribed VEGETATION evolution in time
    real                 :: tmp_topt           ! optimum transpiration air temperature 
    real                 :: tmp_rgl            ! parameter used in radiation stress function
    real                 :: tmp_rsmax          ! maximum stomatal resistance
    real                 :: tmp_rsmin          ! minimum Canopy Resistance [s/m]
    real                 :: tmp_hs             ! parameter used in vapor pressure deficit function 
    real                 :: tmp_nroot
    ! SY: End PARAMETERS for prescribed VEGETATION evolution in time 
    ! SY: End corresponding to REDPRM
    ! SY: Begin corresponding to read_mp_veg_parameters
    real                 :: tmp_CH2OP          !     maximum intercepted h2o per unit lai+sai [mm]
    real                 :: tmp_DLEAF          !     characteristic leaf dimension [m]
    real                 :: tmp_Z0MVT          !     momentum roughness length [m]
    real                 :: tmp_HVT            !     top of canopy [m]
    real                 :: tmp_HVB            !     bottom of canopy [m]
    real                 :: tmp_RC             !     tree crown radius [m]
    real                 :: tmp_RHOL1          !     leaf reflectance (1=vis)
    real                 :: tmp_RHOL2          !     leaf reflectance (2=nir)
    real                 :: tmp_RHOS1          !     stem reflectance (1=vis)
    real                 :: tmp_RHOS2          !     stem reflectance (2=nir)
    real                 :: tmp_TAUL1          !     leaf transmittance (1=vis)
    real                 :: tmp_TAUL2          !     leaf transmittance (2=nir)
    real                 :: tmp_TAUS1          !     stem transmittance (1=vis)
    real                 :: tmp_TAUS2          !     stem transmittance (2=nir)
    real                 :: tmp_XL             !     leaf/stem orientation index
    real                 :: tmp_CWPVT          !     empirical canopy wind parameter
    real                 :: tmp_C3PSN          !     photosynthetic pathway (0. = c4, 1. = c3)
    real                 :: tmp_KC25           !     co2 michaelis-menten constant at 25c [pa]
    real                 :: tmp_AKC            !     q10 for kc25
    real                 :: tmp_KO25           !     o2 michaelis-menten constant at 25c [pa]
    real                 :: tmp_AKO            !     q10 for ko25
    real                 :: tmp_AVCMX          !     q10 for vcmx25
    real                 :: tmp_AQE            !     q10 for qe25
    real                 :: tmp_LTOVRC         !     leaf turnover [1/s]
    real                 :: tmp_DILEFC         !     coeficient for leaf stress death [1/s]
    real                 :: tmp_DILEFW         !     coeficient for leaf stress death [1/s]
    real                 :: tmp_RMF25          !     leaf maintenance respiration at 25c [umol co2/m**2/s]
    real                 :: tmp_SLA            !     single-side leaf area per Kg [m2/kg]
    real                 :: tmp_FRAGR          !     fraction of growth respiration. original was 0.3
    real                 :: tmp_TMIN           !     minimum temperature for photosynthesis [k]
    real                 :: tmp_VCMX25         !     maximum rate of carboxylation at 25c. unit [umol co2/m**2/s]
    real                 :: tmp_TDLEF          !     characteristic T for leaf freezing [K]
    real                 :: tmp_BP             !     minimum leaf conductance [umol/m**2/s]
    real                 :: tmp_MP             !     slope of conductance-to-photosynthesis relationship
    real                 :: tmp_QE25           !     quantum efficiency at 25c [umol co2 / umol photon]
    real                 :: tmp_RMS25          !     stem maintenance respiration at 25c [umol co2/kg bio/s]
    real                 :: tmp_RMR25          !     root maintenance respiration at 25c [umol co2/kg bio/s]
    real                 :: tmp_ARM            !     q10 for maintenance respiration
    real                 :: tmp_FOLNMX         !     foliage nitrogen concentration when f(n)=1 [%]
    real                 :: tmp_WDPOOL         !     wood pool (switch 1 or 0) depending on woody or not [-]
    real                 :: tmp_WRRAT          !     wood to non-wood ratio
    real                 :: tmp_MRP            !     microbial respiration parameter [umol co2 /kg c/ s]
    
    ! SY: End corresponding to read_mp_veg_parameters
    ! SY: End for enabling OPTUE

    ! Code added by Shugong Wang 09/12/2014
    real                 :: soil_temp(NOAHMP36_struc(n)%nsoil)   
    real                 :: snow_temp(NOAHMP36_struc(n)%nsnow)
    real                 :: AvgSurfT_out, Qle_out, Evap_out
    ! Code added by David Mocko
    real                 :: TWS_out, Rainf_out, Snowf_out
    real                 :: bdsno
    ! Code added by Rhae Sung Kim 
    real                 :: layersd

    !ag (12Sep2019)
    real                 :: tmp_rivsto         !     MNB ADD LATER
    real                 :: tmp_fldsto         !     MNB ADD LATER
    real                 :: tmp_fldfrc         !     MNB ADD LATER

    !ag (18Sep2019)
    real,   allocatable   :: rivsto(:)
    real,   allocatable   :: fldsto(:)
    real,   allocatable   :: fldfrc(:)
    real,   allocatable   :: tmp_nensem(:,:,:)

    integer               :: status
    integer               :: c,r
    integer               :: ios, nid,rivid,fldid

    integer            :: tid

    allocate( tmp_sldpth( NOAHMP36_struc(n)%nsoil ) )
    allocate( tmp_shdfac_monthly( 12 ) )
    allocate( tmp_smceq( NOAHMP36_struc(n)%nsoil ) )
    allocate( tmp_sstc( NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow) )
    allocate( tmp_sh2o( NOAHMP36_struc(n)%nsoil ) )
    allocate( tmp_smc( NOAHMP36_struc(n)%nsoil ) )
    allocate( tmp_zss( NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow) )
    allocate( tmp_snowice( NOAHMP36_struc(n)%nsnow ) )
    allocate( tmp_snowliq( NOAHMP36_struc(n)%nsnow ) )
    
    ! check NoahMP36 alarm. If alarm is ring, run model. 
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "NoahMP36 model alarm")
    if (alarmCheck) Then

       do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          dt = LIS_rc%ts
          row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
          col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
          lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
          lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

            ! retrieve forcing data from NOAHMP36_struc(n)%noahmp36(t) and assign to local variables
            ! tair: air temperature
            tmp_tair       = NOAHMP36_struc(n)%noahmp36(t)%tair   / NOAHMP36_struc(n)%forc_count
            NOAHMP36_struc(n)%noahmp36(t)%sfctmp = tmp_tair
            ! psurf: air pressure
            tmp_psurf      = NOAHMP36_struc(n)%noahmp36(t)%psurf  / NOAHMP36_struc(n)%forc_count
 
            ! wind_e: eastward wind speed
            tmp_wind_e     = NOAHMP36_struc(n)%noahmp36(t)%wind_e / NOAHMP36_struc(n)%forc_count
 
            ! wind_n: northward wind speed
            tmp_wind_n     = NOAHMP36_struc(n)%noahmp36(t)%wind_n / NOAHMP36_struc(n)%forc_count
 
            ! qair: near Surface Specific Humidity
            tmp_qair       = NOAHMP36_struc(n)%noahmp36(t)%qair   / NOAHMP36_struc(n)%forc_count
 
            ! swdown: downward solar radiation
            tmp_swdown     = NOAHMP36_struc(n)%noahmp36(t)%swdown / NOAHMP36_struc(n)%forc_count
 
            ! lwdown: downward longwave radiation
            tmp_lwdown     = NOAHMP36_struc(n)%noahmp36(t)%lwdown / NOAHMP36_struc(n)%forc_count
 
            ! prcp: precipitation Rate
            tmp_prcp       = NOAHMP36_struc(n)%noahmp36(t)%prcp   / NOAHMP36_struc(n)%forc_count
            
            !ag(18Sep2019)
            ! rivsto/fldsto: River storage and flood storage
            ! NOAHMP36_struc(n)%noahmp36(t)%rivsto and NOAHMP36_struc(n)%noahmp36(t)%fldsto
            ! are updated in noahmp36_getsws_hymap2.F90
            tmp_rivsto = NOAHMP36_struc(n)%noahmp36(t)%rivsto
            tmp_fldsto = NOAHMP36_struc(n)%noahmp36(t)%fldsto
            tmp_fldfrc = NOAHMP36_struc(n)%noahmp36(t)%fldfrc

!            if(t.eq.13859) write(111,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2,8E14.3)') &
!                 LIS_rc%yr,LIS_rc%mo,LIS_rc%da, LIS_rc%hr, LIS_rc%mn, &
!                 tmp_tair, tmp_qair, tmp_swdown, tmp_lwdown,tmp_wind_n,tmp_wind_e,tmp_psurf,&
!                 tmp_prcp
            ! 
            ! check validity of tair
            if(tmp_tair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable tair in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of psurf
            if(tmp_psurf .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable psurf in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_e
            if(tmp_wind_e .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable wind_e in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_n
            if(tmp_wind_n .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable wind_n in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of qair
            if(tmp_qair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable qair in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of swdown
            if(tmp_swdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable swdown in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of lwdown
            if(tmp_lwdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable lwdown in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of prcp
            if(tmp_prcp .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable prcp in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            !ag (23Sep2019)
            ! check validity of rivsto
            if(tmp_rivsto .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable rivsto in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of fldsto
            if(tmp_fldsto .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable fldsto in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of fldfrc
            if(tmp_fldfrc .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable fldfrc in NoahMP36"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! 
            tmp_latitude  = lat
            tmp_longitude = lon
            tmp_year   = LIS_rc%yr
            tmp_month  = LIS_rc%mo
            tmp_day    = LIS_rc%da
            tmp_hour   = LIS_rc%hr
            tmp_minute = LIS_rc%mn

            ! get parameters
            tmp_landuse_tbl_name                    = NOAHMP36_struc(n)%landuse_tbl_name
            tmp_soil_tbl_name                       = NOAHMP36_struc(n)%soil_tbl_name
            tmp_gen_tbl_name                        = NOAHMP36_struc(n)%gen_tbl_name
            tmp_noahmp_tbl_name                     = NOAHMP36_struc(n)%noahmp_tbl_name
            tmp_landuse_scheme_name                 = NOAHMP36_struc(n)%landuse_scheme_name
            tmp_soil_scheme_name                    = NOAHMP36_struc(n)%soil_scheme_name
            tmp_dveg_opt                            = NOAHMP36_struc(n)%dveg_opt
            tmp_crs_opt                             = NOAHMP36_struc(n)%crs_opt
            tmp_btr_opt                             = NOAHMP36_struc(n)%btr_opt
            tmp_run_opt                             = NOAHMP36_struc(n)%run_opt
            tmp_sfc_opt                             = NOAHMP36_struc(n)%sfc_opt
            tmp_frz_opt                             = NOAHMP36_struc(n)%frz_opt
            tmp_inf_opt                             = NOAHMP36_struc(n)%inf_opt
            tmp_rad_opt                             = NOAHMP36_struc(n)%rad_opt
            tmp_alb_opt                             = NOAHMP36_struc(n)%alb_opt
            tmp_snf_opt                             = NOAHMP36_struc(n)%snf_opt
            tmp_tbot_opt                            = NOAHMP36_struc(n)%tbot_opt
            tmp_stc_opt                             = NOAHMP36_struc(n)%stc_opt
            tmp_nslcats                             = NOAHMP36_struc(n)%nslcats
            tmp_nlucats                             = NOAHMP36_struc(n)%nlucats
            tmp_nslpcats                            = NOAHMP36_struc(n)%nslpcats
            tmp_dt                                  = NOAHMP36_struc(n)%dt
            tmp_nsoil                               = NOAHMP36_struc(n)%nsoil
            tmp_sldpth(:)                           = NOAHMP36_struc(n)%sldpth(:)
            tmp_nsnow                               = NOAHMP36_struc(n)%nsnow
            tmp_shdfac_monthly(:)                   = NOAHMP36_struc(n)%noahmp36(t)%shdfac_monthly(:)
            tmp_vegetype                            = NOAHMP36_struc(n)%noahmp36(t)%vegetype
            tmp_soiltype                            = NOAHMP36_struc(n)%noahmp36(t)%soiltype
            tmp_slopetype                           = NOAHMP36_struc(n)%noahmp36(t)%slopetype
            tmp_urban_vegetype                      = NOAHMP36_struc(n)%urban_vegetype
            tmp_ice_flag                            = NOAHMP36_struc(n)%ice_flag
            tmp_st_flag                             = NOAHMP36_struc(n)%st_flag
            tmp_sc_idx                              = NOAHMP36_struc(n)%sc_idx
            tmp_iz0tlnd                             = NOAHMP36_struc(n)%iz0tlnd
            tmp_smceq(:)                            = NOAHMP36_struc(n)%noahmp36(t)%smceq(:)
            tmp_tbot                                = NOAHMP36_struc(n)%noahmp36(t)%tbot
            tmp_pblh                                = NOAHMP36_struc(n)%noahmp36(t)%pblh
            !tmp_zlvl                                = NOAHMP36_struc(n)%zlvl

            ! get state variables
            tmp_albold      = NOAHMP36_struc(n)%noahmp36(t)%albold
            tmp_sneqvo      = NOAHMP36_struc(n)%noahmp36(t)%sneqvo
            tmp_sstc(:)     = NOAHMP36_struc(n)%noahmp36(t)%sstc(:)
            tmp_sh2o(:)     = NOAHMP36_struc(n)%noahmp36(t)%sh2o(:)
            tmp_smc(:)      = NOAHMP36_struc(n)%noahmp36(t)%smc(:)
            tmp_tah         = NOAHMP36_struc(n)%noahmp36(t)%tah
            tmp_eah         = NOAHMP36_struc(n)%noahmp36(t)%eah
            tmp_fwet        = NOAHMP36_struc(n)%noahmp36(t)%fwet
            tmp_canliq      = NOAHMP36_struc(n)%noahmp36(t)%canliq
            tmp_canice      = NOAHMP36_struc(n)%noahmp36(t)%canice
            tmp_tv          = NOAHMP36_struc(n)%noahmp36(t)%tv
            tmp_tg          = NOAHMP36_struc(n)%noahmp36(t)%tg
            tmp_qsnow       = NOAHMP36_struc(n)%noahmp36(t)%qsnow
            tmp_isnow       = NOAHMP36_struc(n)%noahmp36(t)%isnow
            tmp_zss(:)      = NOAHMP36_struc(n)%noahmp36(t)%zss(:)
            tmp_snowh       = NOAHMP36_struc(n)%noahmp36(t)%snowh
            tmp_sneqv       = NOAHMP36_struc(n)%noahmp36(t)%sneqv
            tmp_snowice(:)  = NOAHMP36_struc(n)%noahmp36(t)%snowice(:)
            tmp_snowliq(:)  = NOAHMP36_struc(n)%noahmp36(t)%snowliq(:)
            tmp_zwt         = NOAHMP36_struc(n)%noahmp36(t)%zwt
            tmp_wa          = NOAHMP36_struc(n)%noahmp36(t)%wa
            tmp_wt          = NOAHMP36_struc(n)%noahmp36(t)%wt
            tmp_wslake      = NOAHMP36_struc(n)%noahmp36(t)%wslake
            tmp_lfmass      = NOAHMP36_struc(n)%noahmp36(t)%lfmass
            tmp_rtmass      = NOAHMP36_struc(n)%noahmp36(t)%rtmass
            tmp_stmass      = NOAHMP36_struc(n)%noahmp36(t)%stmass
            tmp_wood        = NOAHMP36_struc(n)%noahmp36(t)%wood
            tmp_stblcp      = NOAHMP36_struc(n)%noahmp36(t)%stblcp
            tmp_fastcp      = NOAHMP36_struc(n)%noahmp36(t)%fastcp
            tmp_lai         = NOAHMP36_struc(n)%noahmp36(t)%lai
            tmp_sai         = NOAHMP36_struc(n)%noahmp36(t)%sai
            tmp_cm          = NOAHMP36_struc(n)%noahmp36(t)%cm
            tmp_ch          = NOAHMP36_struc(n)%noahmp36(t)%ch
            tmp_tauss       = NOAHMP36_struc(n)%noahmp36(t)%tauss
            tmp_smcwtd      = NOAHMP36_struc(n)%noahmp36(t)%smcwtd
            tmp_deeprech    = NOAHMP36_struc(n)%noahmp36(t)%deeprech
            tmp_rech        = NOAHMP36_struc(n)%noahmp36(t)%rech
            tmp_zlvl        = NOAHMP36_struc(n)%noahmp36(t)%zlvl
            tmp_albd      = NOAHMP36_struc(n)%noahmp36(t)%albd
            tmp_albi      = NOAHMP36_struc(n)%noahmp36(t)%albi

            ! SY: Begin for enabling OPTUE: get calibratable parameters
            ! SY: Begin corresponding to REDPRM
            tmp_csoil       = NOAHMP36_struc(n)%noahmp36(t)%csoil
            tmp_bexp        = NOAHMP36_struc(n)%noahmp36(t)%bexp
            tmp_dksat       = NOAHMP36_struc(n)%noahmp36(t)%dksat
            tmp_dwsat       = NOAHMP36_struc(n)%noahmp36(t)%dwsat
            tmp_psisat      = NOAHMP36_struc(n)%noahmp36(t)%psisat
            tmp_quartz      = NOAHMP36_struc(n)%noahmp36(t)%quartz
            tmp_smcmax      = NOAHMP36_struc(n)%noahmp36(t)%smcmax
            tmp_smcref      = NOAHMP36_struc(n)%noahmp36(t)%smcref
            tmp_smcwlt      = NOAHMP36_struc(n)%noahmp36(t)%smcwlt
            tmp_czil        = NOAHMP36_struc(n)%noahmp36(t)%czil
            tmp_frzk        = NOAHMP36_struc(n)%noahmp36(t)%frzk
            tmp_refdk       = NOAHMP36_struc(n)%noahmp36(t)%refdk
            tmp_refkdt      = NOAHMP36_struc(n)%noahmp36(t)%refkdt
            tmp_slope       = NOAHMP36_struc(n)%noahmp36(t)%slope
            tmp_topt        = NOAHMP36_struc(n)%noahmp36(t)%topt
            tmp_rgl         = NOAHMP36_struc(n)%noahmp36(t)%rgl
            tmp_rsmax       = NOAHMP36_struc(n)%noahmp36(t)%rsmax
            tmp_rsmin       = NOAHMP36_struc(n)%noahmp36(t)%rsmin
            tmp_hs          = NOAHMP36_struc(n)%noahmp36(t)%hs
            tmp_nroot       = NOAHMP36_struc(n)%noahmp36(t)%nroot
            ! SY: End corresponding to REDPRM
            ! SY: Begin corresponding to read_mp_veg_parameters
            tmp_CH2OP       = NOAHMP36_struc(n)%noahmp36(t)%CH2OP
            tmp_DLEAF       = NOAHMP36_struc(n)%noahmp36(t)%DLEAF
            tmp_Z0MVT       = NOAHMP36_struc(n)%noahmp36(t)%Z0MVT
            tmp_HVT         = NOAHMP36_struc(n)%noahmp36(t)%HVT
            tmp_HVB         = NOAHMP36_struc(n)%noahmp36(t)%HVB
            tmp_RC          = NOAHMP36_struc(n)%noahmp36(t)%RC
            tmp_RHOL1       = NOAHMP36_struc(n)%noahmp36(t)%RHOL1
            tmp_RHOL2       = NOAHMP36_struc(n)%noahmp36(t)%RHOL2
            tmp_RHOS1       = NOAHMP36_struc(n)%noahmp36(t)%RHOS1
            tmp_RHOS2       = NOAHMP36_struc(n)%noahmp36(t)%RHOS2
            tmp_TAUL1       = NOAHMP36_struc(n)%noahmp36(t)%TAUL1
            tmp_TAUL2       = NOAHMP36_struc(n)%noahmp36(t)%TAUL2
            tmp_TAUS1       = NOAHMP36_struc(n)%noahmp36(t)%TAUS1
            tmp_TAUS2       = NOAHMP36_struc(n)%noahmp36(t)%TAUS2
            tmp_XL          = NOAHMP36_struc(n)%noahmp36(t)%XL
            tmp_CWPVT       = NOAHMP36_struc(n)%noahmp36(t)%CWPVT
            tmp_C3PSN       = NOAHMP36_struc(n)%noahmp36(t)%C3PSN
            tmp_KC25        = NOAHMP36_struc(n)%noahmp36(t)%KC25
            tmp_AKC         = NOAHMP36_struc(n)%noahmp36(t)%AKC
            tmp_KO25        = NOAHMP36_struc(n)%noahmp36(t)%KO25
            tmp_AKO         = NOAHMP36_struc(n)%noahmp36(t)%AKO
            tmp_AVCMX       = NOAHMP36_struc(n)%noahmp36(t)%AVCMX
            tmp_AQE         = NOAHMP36_struc(n)%noahmp36(t)%AQE
            tmp_LTOVRC      = NOAHMP36_struc(n)%noahmp36(t)%LTOVRC
            tmp_DILEFC      = NOAHMP36_struc(n)%noahmp36(t)%DILEFC
            tmp_DILEFW      = NOAHMP36_struc(n)%noahmp36(t)%DILEFW
            tmp_RMF25       = NOAHMP36_struc(n)%noahmp36(t)%RMF25
            tmp_SLA         = NOAHMP36_struc(n)%noahmp36(t)%SLA
            tmp_FRAGR       = NOAHMP36_struc(n)%noahmp36(t)%FRAGR
            tmp_TMIN        = NOAHMP36_struc(n)%noahmp36(t)%TMIN
            tmp_VCMX25      = NOAHMP36_struc(n)%noahmp36(t)%VCMX25
            tmp_TDLEF       = NOAHMP36_struc(n)%noahmp36(t)%TDLEF
            tmp_BP          = NOAHMP36_struc(n)%noahmp36(t)%BP
            tmp_MP          = NOAHMP36_struc(n)%noahmp36(t)%MP
            tmp_QE25        = NOAHMP36_struc(n)%noahmp36(t)%QE25
            tmp_RMS25       = NOAHMP36_struc(n)%noahmp36(t)%RMS25
            tmp_RMR25       = NOAHMP36_struc(n)%noahmp36(t)%RMR25
            tmp_ARM         = NOAHMP36_struc(n)%noahmp36(t)%ARM
            tmp_FOLNMX      = NOAHMP36_struc(n)%noahmp36(t)%FOLNMX
            tmp_WDPOOL      = NOAHMP36_struc(n)%noahmp36(t)%WDPOOL
            tmp_WRRAT       = NOAHMP36_struc(n)%noahmp36(t)%WRRAT
            tmp_MRP         = NOAHMP36_struc(n)%noahmp36(t)%MRP
            ! SY: End corresponding to read_mp_veg_parameters
            ! SY: End for enabling OPTUE: get calibratable parameters
            call noahmp_driver_36(LIS_localPet, t,tmp_landuse_tbl_name  , & ! in    - Noah model landuse parameter table [-]
                                  tmp_soil_tbl_name     , & ! in    - Noah model soil parameter table [-]
                                  tmp_gen_tbl_name      , & ! in    - Noah model general parameter table [-]
                                  tmp_noahmp_tbl_name   , & ! in    - NoahMP parameter table [-]
                                  tmp_landuse_scheme_name, & ! in    - landuse classification scheme [-]
                                  tmp_soil_scheme_name  , & ! in    - soil classification scheme [-]
                                  tmp_dveg_opt          , & ! in    - vegetation model [-]
                                  tmp_crs_opt           , & ! in    - canopy stomatal resistance [-]
                                  tmp_btr_opt           , & ! in    - soil moisture factor for stomatal resistance [-]
                                  tmp_run_opt           , & ! in    - runoff and groundwater [-]
                                  tmp_sfc_opt           , & ! in    - surface layer drag coefficients (CH & CM) [-]
                                  tmp_frz_opt           , & ! in    - supercooled liquid water [-]
                                  tmp_inf_opt           , & ! in    - frozen soil permeability [-]
                                  tmp_rad_opt           , & ! in    - radiation transfer [-]
                                  tmp_alb_opt           , & ! in    - snow surface albedo [-]
                                  tmp_snf_opt           , & ! in    - rainfall & snowfall [-]
                                  tmp_tbot_opt          , & ! in    - lower boundary of soil temperature [-]
                                  tmp_stc_opt           , & ! in    - snow/soil temperature time scheme [-]
                                  tmp_nslcats           , & ! in    - the number of total soil types in parameter table [-]
                                  tmp_nlucats           , & ! in    - the number of total land cover types in parameter table [-]
                                  tmp_nslpcats          , & ! in    - the number of total slope category for Noah baseflow [-]
                                  tmp_latitude          , & ! in    - latitude in decimal degree [-]
                                  tmp_longitude         , & ! in    - longitude in decimal degree [-]
                                  tmp_year              , & ! in    - year of the current time step [-]
                                  tmp_month             , & ! in    - month of the current time step [-]
                                  tmp_day               , & ! in    - day of the current time step [-]
                                  tmp_hour              , & ! in    - hour of the current time step [-]
                                  tmp_minute            , & ! in    - minute of the current time step [-]
                                  tmp_dt                , & ! in    - time step in seconds [s]
                                  NOAHMP36_struc(n)%noahmp36(t)%alb_upd_flag,&
                                  tmp_albd            , & 
                                  tmp_albi            , & 
                                  tmp_nsoil             , & ! in    - number of soil layers [-]
                                  tmp_sldpth            , & ! in    - thickness of soil layers [-]
                                  tmp_nsnow             , & ! in    - maximum number of snow layers [-]
                                  tmp_shdfac_monthly    , & ! in    - monthly values for green vegetation fraction [-]
                                  tmp_vegetype          , & ! in    - land cover type index [-]
                                  tmp_soiltype          , & ! in    - soil type index [-]
                                  tmp_slopetype         , & ! in    - slope type for Noah baseflow [-]
                                  tmp_urban_vegetype    , & ! in    - urban land cover type index [-]
                                  tmp_ice_flag          , & ! in    - ice flag: 0 = no ice, 1 = ice [-]
                                  tmp_st_flag           , & ! in    - surface type 1=soil, 2=lake [-]
                                  tmp_sc_idx            , & ! in    - soil color type [-]
                                  tmp_iz0tlnd           , & ! in    - option of Chen adjustment of Czil [-]
                                  tmp_smceq             , & ! in    - equilibrium soil water content [m^3 m-3]
                                  tmp_tair              , & ! in    - air temperature [K]
                                  tmp_psurf             , & ! in    - air pressure [Pa]
                                  tmp_wind_e            , & ! in    - eastward wind speed [m s-1]
                                  tmp_wind_n            , & ! in    - northward wind speed [m s-1]
                                  tmp_qair              , & ! in    - near Surface Specific Humidity [kg kg-1]
                                  tmp_swdown            , & ! in    - downward solar radiation [w/m2]
                                  tmp_lwdown            , & ! in    - downward longwave radiation [w/m2]
                                  tmp_prcp              , & ! in    - total precip (rainfall+snowfall) Rate [kg m-2 s-1]   
                                  tmp_tbot              , & ! in    - deep-layer soil temperature [K]
                                  tmp_pblh              , & ! in    - planetary boundary layer height [m]
                                  tmp_zlvl              , & ! in    - reference height of temperature and humidity [m]
                                  tmp_csoil             , & ! in    - vol. soil heat capacity [j/m3/K] 
                                  tmp_bexp              , & ! in    - B parameter 
                                  tmp_dksat             , & ! in    - saturated soil hydraulic conductivity
                                  tmp_dwsat             , & ! in    - saturated soil hydraulic diffusivity
                                  tmp_psisat            , & ! in    - saturated soil matric potential 
                                  tmp_quartz            , & ! in    - soil quartz content 
                                  tmp_smcmax            , & ! in    - porosity, saturated value of soil moisture (volumetric) 
                                  tmp_smcref            , & ! in    - reference soil moisture (field capacity) 
                                  tmp_smcwlt            , & ! in    - wilting point soil moisture (volumetric) 
                                  tmp_czil              , & ! in    - Calculate roughness length of heat 
                                  tmp_frzk              , & ! in    - frozen ground parameter 
                                  tmp_refdk             , & ! in    - parameters in the surface runoff parameteriz.
                                  tmp_refkdt            , & ! in    - parameters in the surface runoff parameteriz. 
                                  tmp_slope             , & ! in    - slope index (0 - 1)
                                  tmp_topt              , & ! in    - optimum transpiration air temperature 
                                  tmp_rgl               , & ! in    - parameter used in radiation stress function
                                  tmp_rsmax             , & ! in    - maximum stomatal resistance
                                  tmp_rsmin             , & ! in    - minimum Canopy Resistance [s/m]
                                  tmp_hs                , & ! in    - parameter used in vapor pressure deficit function 
                                  tmp_nroot             , & ! in    - 
                                  tmp_CH2OP             , & ! in    - maximum intercepted h2o per unit lai+sai [mm]
                                  tmp_DLEAF             , & ! in    - characteristic leaf dimension [m]
                                  tmp_Z0MVT             , & ! in    - momentum roughness length [m]
                                  tmp_HVT               , & ! in    - top of canopy [m]
                                  tmp_HVB               , & ! in    - bottom of canopy [m]
                                  tmp_RC                , & ! in    - tree crown radius [m]
                                  tmp_RHOL1             , & ! in    - leaf reflectance (1=vis)
                                  tmp_RHOL2             , & ! in    - leaf reflectance (2=nir)
                                  tmp_RHOS1             , & ! in    - stem reflectance (1=vis)
                                  tmp_RHOS2             , & ! in    - stem reflectance (2=nir)
                                  tmp_TAUL1             , & ! in    - leaf transmittance (1=vis)
                                  tmp_TAUL2             , & ! in    - leaf transmittance (2=nir)
                                  tmp_TAUS1             , & ! in    - stem transmittance (1=vis)
                                  tmp_TAUS2             , & ! in    - stem transmittance (2=nir)
                                  tmp_XL                , & ! in    - leaf/stem orientation index
                                  tmp_CWPVT             , & ! in    - empirical canopy wind parameter
                                  tmp_C3PSN             , & ! in    - photosynthetic pathway (0. = c4, 1. = c3)
                                  tmp_KC25              , & ! in    - co2 michaelis-menten constant at 25c [pa]
                                  tmp_AKC               , & ! in    - q10 for kc25
                                  tmp_KO25              , & ! in    - o2 michaelis-menten constant at 25c [pa]
                                  tmp_AKO               , & ! in    - q10 for ko25
                                  tmp_AVCMX             , & ! in    - q10 for vcmx25
                                  tmp_AQE               , & ! in    - q10 for qe25
                                  tmp_LTOVRC            , & ! in    - leaf turnover [1/s]
                                  tmp_DILEFC            , & ! in    - coeficient for leaf stress death [1/s]
                                  tmp_DILEFW            , & ! in    - coeficient for leaf stress death [1/s]
                                  tmp_RMF25             , & ! in    - leaf maintenance respiration at 25c [umol co2/m**2/s]
                                  tmp_SLA               , & ! in    - single-side leaf area per Kg [m2/kg]
                                  tmp_FRAGR             , & ! in    - fraction of growth respiration. original was 0.3
                                  tmp_TMIN              , & ! in    - minimum temperature for photosynthesis [k]
                                  tmp_VCMX25            , & ! in    - maximum rate of carboxylation at 25c. unit [umol co2/m**2/s]
                                  tmp_TDLEF             , & ! in    - characteristic T for leaf freezing [K]
                                  tmp_BP                , & ! in    - minimum leaf conductance [umol/m**2/s]
                                  tmp_MP                , & ! in    - slope of conductance-to-photosynthesis relationship
                                  tmp_QE25              , & ! in    - quantum efficiency at 25c [umol co2 / umol photon]
                                  tmp_RMS25             , & ! in    - stem maintenance respiration at 25c [umol co2/kg bio/s]
                                  tmp_RMR25             , & ! in    - root maintenance respiration at 25c [umol co2/kg bio/s]
                                  tmp_ARM               , & ! in    - q10 for maintenance respiration
                                  tmp_FOLNMX            , & ! in    - foliage nitrogen concentration when f(n)=1 [%]
                                  tmp_WDPOOL            , & ! in    - wood pool (switch 1 or 0) depending on woody or not [-]
                                  tmp_WRRAT             , & ! in    - wood to non-wood ratio
                                  tmp_MRP               , & ! in    - microbial respiration parameter [umol co2 /kg c/ s]
                                  tmp_albold            , & ! inout - snow albedo at last time step [-]
                                  tmp_sneqvo            , & ! inout - snow mass at the last time step [mm]
                                  tmp_sstc              , & ! inout - snow/soil temperature [K]
                                  tmp_sh2o              , & ! inout - volumetric liquid soil moisture [m^3 m-3]
                                  tmp_smc               , & ! inout - volumetric soil moisture, ice + liquid [m^3 m-3]
                                  tmp_tah               , & ! inout - canopy air temperature [K]
                                  tmp_eah               , & ! inout - canopy air vapor pressure [Pa]
                                  tmp_fwet              , & ! inout - wetted or snowed fraction of canopy [-]
                                  tmp_canliq            , & ! inout - intercepted liquid water [mm]
                                  tmp_canice            , & ! inout - intercepted ice mass [mm]
                                  tmp_tv                , & ! inout - vegetation temperature [K]
                                  tmp_tg                , & ! inout - ground temperature [K]
                                  tmp_qsnow             , & ! inout - snowfall on the ground [mm s-1]
                                  tmp_isnow             , & ! inout - actual number of snow layers [-]
                                  tmp_zss               , & ! inout - snow/soil layer-bottom depth from snow surface [m]
                                  tmp_snowh             , & ! inout - snow height [m]
                                  tmp_sneqv             , & ! inout - snow water equivalent [mm]
                                  tmp_snowice           , & ! inout - snow-layer ice [mm]
                                  tmp_snowliq           , & ! inout - snow-layer liquid water [mm]
                                  tmp_zwt               , & ! inout - depth to water table [m]
                                  tmp_wa                , & ! inout - water storage in aquifer [mm]
                                  tmp_wt                , & ! inout - water in aquifer and saturated soil [mm]
                                  tmp_wslake            , & ! inout - lake water storage [mm]
                                  tmp_lfmass            , & ! inout - leaf mass [g/m2]
                                  tmp_rtmass            , & ! inout - mass of fine roots [g/m2]
                                  tmp_stmass            , & ! inout - stem mass [g/m2]
                                  tmp_wood              , & ! inout - mass of wood including woody roots [g/m2]
                                  tmp_stblcp            , & ! inout - stable carbon in deep soil [g/m2]
                                  tmp_fastcp            , & ! inout - short-lived carbon in shallow soil [g/m2]
                                  tmp_lai               , & ! inout - leaf area index [-]
                                  tmp_sai               , & ! inout - stem area index [-]
                                  tmp_cm                , & ! inout - momentum drag coefficient [m s-1]
                                  tmp_ch                , & ! inout - sensible heat exchange coefficient [m s-1]
                                  tmp_tauss             , & ! inout - snow aging term [-]
                                  tmp_smcwtd            , & ! inout - soil water content between bottom of the soil and water table [m^3 m-3]
                                  tmp_deeprech          , & ! inout - recharge to or from the water table when deep [m]
                                  tmp_rech              , & ! inout - recharge to or from the water table when shallow [m]
                                  tmp_fsa               , & ! out   - total absorbed solar radiation [W m-2]
                                  tmp_fsr               , & ! out   - total reflected solar radiation [W m-2]
                                  tmp_fira              , & ! out   - total net longwave radiation to atmosphere [W m-2]
                                  tmp_fsh               , & ! out   - total sensible heat to atmosphere [W m-2]
                                  tmp_ssoil             , & ! out   - ground heat flux to soil [W m-2]
                                  tmp_fcev              , & ! out   - canopy evaporative heat to atmosphere [W m-2]
                                  tmp_fgev              , & ! out   - ground evaporative heat to atmosphere [W m-2]
                                  tmp_fctr              , & ! out   - transpiration heat to atmosphere [W m-2]
                                  tmp_ecan              , & ! out   - evaporation rate of canopy water [kg m-2 s-1]
                                  tmp_etran             , & ! out   - transpiration rate [kg m-2 s-1]
                                  tmp_edir              , & ! out   - direct evaporation rate from surface [kg m-2 s-1]
                                  tmp_trad              , & ! out   - surface radiative temperature [K]
                                  tmp_subsnow           , & ! out   - snow sublimation rate [kg m-2 s-1]
                                  tmp_tgb               , & ! out   - ground temperature [K]
                                  tmp_tgv               , & ! out   - ground surface temperature [K]
                                  tmp_t2mv              , & ! out   - 2-m air temperature over vegetated part [K]
                                  tmp_t2mb              , & ! out   - 2-m height air temperature [K]
                                  tmp_q2v               , & ! out   - 2-m specific humidity over vegetation [kg kg-1]
                                  tmp_q2b               , & ! out   - 2-m air specific humidity [kg kg-1]
                                  tmp_runsrf            , & ! out   - surface runoff [kg m-2 s-1]
                                  tmp_runsub            , & ! out   - baseflow (saturation excess) [kg m-2 s-1]
                                  tmp_apar              , & ! out   - photosynthesis active energy by canopy [W m-2]
                                  tmp_psn               , & ! out   - total photosynthesis of CO2 [umol m-2 s-1]
                                  tmp_sav               , & ! out   - solar radiation absorbed by vegetation [W m-2]
                                  tmp_sag               , & ! out   - solar radiation absorbed by ground [W m-2]
                                  tmp_fsno              , & ! out   - snow-cover fraction on the ground [-]
                                  tmp_nee               , & ! out   - net ecosystem exchange of CO2 [g/m2s]
                                  tmp_gpp               , & ! out   - net instantaneous assimilation of carbon [g/m2s]
                                  tmp_npp               , & ! out   - net primary productivity of carbon [g/m2s]
                                  tmp_fveg              , & ! out   - green vegetation fraction [-]
                                  tmp_albedo            , & ! out   - surface albedo [-]
                                  tmp_qsnbot            , & ! out   - melting water out of snow bottom [kg m-2 s-1]
                                  tmp_ponding           , & ! out   - surface ponding [mm]
                                  tmp_ponding1          , & ! out   - surface ponding1 [mm]
                                  tmp_ponding2          , & ! out   - surface ponding2 [mm]
                                  tmp_rssun             , & ! out   - sunlit stomatal resistance [s/m]
                                  tmp_rssha             , & ! out   - shaded stomatal resistance [s/m]
                                  tmp_bgap              , & ! out   - between canopy gap fraction for beam [-]
                                  tmp_wgap              , & ! out   - within canopy gap fraction for beam [-]
                                  tmp_chv               , & ! out   - sensible heat exchange coefficient over vegetated fraction [m s-1]
                                  tmp_chb               , & ! out   - sensible heat exchange coefficient over bare-ground fraction [m s-1]
                                  tmp_emissi            , & ! out   - surface emissivity [-]
                                  tmp_shg               , & ! out   - ground sensible heat [W m-2]
                                  tmp_shc               , & ! out   - canopy sensible heat [W m-2]
                                  tmp_shb               , & ! out   - bare ground sensible heat [W m-2]
                                  tmp_evg               , & ! out   - ground evaporation heat [W m-2]
                                  tmp_evb               , & ! out   - bare ground evaporation heat [W m-2]
                                  tmp_ghv               , & ! out   - ground heat flux [W m-2]
                                  tmp_ghb               , & ! out   - bare ground heat flux [W m-2]
                                  tmp_irg               , & ! out   - ground net long wave radiation [W m-2]
                                  tmp_irc               , & ! out   - canopy net long wave radiation [W m-2]
                                  tmp_irb               , & ! out   - bare ground net long wave radiation [W m-2]
                                  tmp_tr                , & ! out   - transpiration heat [W m-2]
                                  tmp_evc               , & ! out   - canopy evaporation heat [W m-2]
                                  tmp_chleaf            , & ! out   - leaf exchange coefficient [-]
                                  tmp_chuc              , & ! out   - under canopy exchange coefficient [-]
                                  tmp_chv2              , & ! out   - sensible heat exchange coefficient over vegetated fraction [-]
                                  tmp_chb2              , & ! out   - sensible heat exchange coefficient over bare-ground [-]
                                  tmp_fpice             , & ! out   - snow fraction in precipitation [-]
                                  !ag (12Sep2019)
                                  tmp_rivsto            , & ! in   - river storage [m/s] 
                                  tmp_fldsto            , & ! in   - flood storage [m/s]
                                  tmp_fldfrc            , & ! in   - flood storage [m/s]
                                  
                                  tmp_sfcheadrt         )   ! out   - extra output for WRF-HYDRO [m]
            
           ! save state variables from local variables to global variables
            NOAHMP36_struc(n)%noahmp36(t)%albold      = tmp_albold
            NOAHMP36_struc(n)%noahmp36(t)%sneqvo      = tmp_sneqvo
            NOAHMP36_struc(n)%noahmp36(t)%sstc(:)     = tmp_sstc(:)
            NOAHMP36_struc(n)%noahmp36(t)%sh2o(:)     = tmp_sh2o(:)
            NOAHMP36_struc(n)%noahmp36(t)%smc(:)      = tmp_smc(:)
            NOAHMP36_struc(n)%noahmp36(t)%tah         = tmp_tah
            NOAHMP36_struc(n)%noahmp36(t)%eah         = tmp_eah
            NOAHMP36_struc(n)%noahmp36(t)%fwet        = tmp_fwet
            NOAHMP36_struc(n)%noahmp36(t)%canliq      = tmp_canliq
            NOAHMP36_struc(n)%noahmp36(t)%canice      = tmp_canice
            NOAHMP36_struc(n)%noahmp36(t)%tv          = tmp_tv
            NOAHMP36_struc(n)%noahmp36(t)%tg          = tmp_tg
            NOAHMP36_struc(n)%noahmp36(t)%qsnow       = tmp_qsnow
            NOAHMP36_struc(n)%noahmp36(t)%isnow       = tmp_isnow
            NOAHMP36_struc(n)%noahmp36(t)%zss(:)      = tmp_zss(:)
            NOAHMP36_struc(n)%noahmp36(t)%snowh       = tmp_snowh
            NOAHMP36_struc(n)%noahmp36(t)%sneqv       = tmp_sneqv
            NOAHMP36_struc(n)%noahmp36(t)%snowice(:)  = tmp_snowice(:)
            NOAHMP36_struc(n)%noahmp36(t)%snowliq(:)  = tmp_snowliq(:)
            NOAHMP36_struc(n)%noahmp36(t)%zwt         = tmp_zwt
            NOAHMP36_struc(n)%noahmp36(t)%wa          = tmp_wa
            NOAHMP36_struc(n)%noahmp36(t)%wt          = tmp_wt
            NOAHMP36_struc(n)%noahmp36(t)%wslake      = tmp_wslake
            NOAHMP36_struc(n)%noahmp36(t)%lfmass      = tmp_lfmass
            NOAHMP36_struc(n)%noahmp36(t)%rtmass      = tmp_rtmass
            NOAHMP36_struc(n)%noahmp36(t)%stmass      = tmp_stmass
            NOAHMP36_struc(n)%noahmp36(t)%wood        = tmp_wood
            NOAHMP36_struc(n)%noahmp36(t)%stblcp      = tmp_stblcp
            NOAHMP36_struc(n)%noahmp36(t)%fastcp      = tmp_fastcp
            NOAHMP36_struc(n)%noahmp36(t)%lai         = tmp_lai
            NOAHMP36_struc(n)%noahmp36(t)%sai         = tmp_sai
            NOAHMP36_struc(n)%noahmp36(t)%cm          = tmp_cm
            NOAHMP36_struc(n)%noahmp36(t)%ch          = tmp_ch
            NOAHMP36_struc(n)%noahmp36(t)%tauss       = tmp_tauss
            NOAHMP36_struc(n)%noahmp36(t)%smcwtd      = tmp_smcwtd
            NOAHMP36_struc(n)%noahmp36(t)%deeprech    = tmp_deeprech
            NOAHMP36_struc(n)%noahmp36(t)%rech        = tmp_rech

            ! save output variables from local variables to global variables
            NOAHMP36_struc(n)%noahmp36(t)%fsa          = tmp_fsa
            NOAHMP36_struc(n)%noahmp36(t)%fsr          = tmp_fsr
            NOAHMP36_struc(n)%noahmp36(t)%fira         = tmp_fira
            NOAHMP36_struc(n)%noahmp36(t)%fsh          = tmp_fsh
            NOAHMP36_struc(n)%noahmp36(t)%ssoil        = tmp_ssoil
            NOAHMP36_struc(n)%noahmp36(t)%fcev         = tmp_fcev
            NOAHMP36_struc(n)%noahmp36(t)%fgev         = tmp_fgev
            NOAHMP36_struc(n)%noahmp36(t)%fctr         = tmp_fctr
            NOAHMP36_struc(n)%noahmp36(t)%ecan         = tmp_ecan
            NOAHMP36_struc(n)%noahmp36(t)%etran        = tmp_etran
            NOAHMP36_struc(n)%noahmp36(t)%edir         = tmp_edir
            NOAHMP36_struc(n)%noahmp36(t)%trad         = tmp_trad
            NOAHMP36_struc(n)%noahmp36(t)%tgb          = tmp_tgb
            NOAHMP36_struc(n)%noahmp36(t)%tgv          = tmp_tgv
            NOAHMP36_struc(n)%noahmp36(t)%t2mv         = tmp_t2mv
            NOAHMP36_struc(n)%noahmp36(t)%t2mb         = tmp_t2mb
            NOAHMP36_struc(n)%noahmp36(t)%q2v          = tmp_q2v
            NOAHMP36_struc(n)%noahmp36(t)%q2b          = tmp_q2b
            NOAHMP36_struc(n)%noahmp36(t)%runsrf       = tmp_runsrf
            NOAHMP36_struc(n)%noahmp36(t)%runsub       = tmp_runsub
            NOAHMP36_struc(n)%noahmp36(t)%apar         = tmp_apar
            NOAHMP36_struc(n)%noahmp36(t)%psn          = tmp_psn
            NOAHMP36_struc(n)%noahmp36(t)%sav          = tmp_sav
            NOAHMP36_struc(n)%noahmp36(t)%sag          = tmp_sag
            NOAHMP36_struc(n)%noahmp36(t)%fsno         = tmp_fsno
            NOAHMP36_struc(n)%noahmp36(t)%nee          = tmp_nee
            NOAHMP36_struc(n)%noahmp36(t)%gpp          = tmp_gpp
            NOAHMP36_struc(n)%noahmp36(t)%npp          = tmp_npp
            NOAHMP36_struc(n)%noahmp36(t)%fveg         = tmp_fveg
            NOAHMP36_struc(n)%noahmp36(t)%albedo       = tmp_albedo
            NOAHMP36_struc(n)%noahmp36(t)%qsnbot       = tmp_qsnbot
            NOAHMP36_struc(n)%noahmp36(t)%ponding      = tmp_ponding
            NOAHMP36_struc(n)%noahmp36(t)%ponding1     = tmp_ponding1
            NOAHMP36_struc(n)%noahmp36(t)%ponding2     = tmp_ponding2
            NOAHMP36_struc(n)%noahmp36(t)%rssun        = tmp_rssun
            NOAHMP36_struc(n)%noahmp36(t)%rssha        = tmp_rssha
            NOAHMP36_struc(n)%noahmp36(t)%bgap         = tmp_bgap
            NOAHMP36_struc(n)%noahmp36(t)%wgap         = tmp_wgap
            NOAHMP36_struc(n)%noahmp36(t)%chv          = tmp_chv
            NOAHMP36_struc(n)%noahmp36(t)%chb          = tmp_chb
            NOAHMP36_struc(n)%noahmp36(t)%emissi       = tmp_emissi
            NOAHMP36_struc(n)%noahmp36(t)%shg          = tmp_shg
            NOAHMP36_struc(n)%noahmp36(t)%shc          = tmp_shc
            NOAHMP36_struc(n)%noahmp36(t)%shb          = tmp_shb
            NOAHMP36_struc(n)%noahmp36(t)%evg          = tmp_evg
            NOAHMP36_struc(n)%noahmp36(t)%evb          = tmp_evb
            NOAHMP36_struc(n)%noahmp36(t)%ghv          = tmp_ghv
            NOAHMP36_struc(n)%noahmp36(t)%ghb          = tmp_ghb
            NOAHMP36_struc(n)%noahmp36(t)%irg          = tmp_irg
            NOAHMP36_struc(n)%noahmp36(t)%irc          = tmp_irc
            NOAHMP36_struc(n)%noahmp36(t)%irb          = tmp_irb
            NOAHMP36_struc(n)%noahmp36(t)%tr           = tmp_tr
            NOAHMP36_struc(n)%noahmp36(t)%evc          = tmp_evc
            NOAHMP36_struc(n)%noahmp36(t)%chleaf       = tmp_chleaf
            NOAHMP36_struc(n)%noahmp36(t)%chuc         = tmp_chuc
            NOAHMP36_struc(n)%noahmp36(t)%chv2         = tmp_chv2
            NOAHMP36_struc(n)%noahmp36(t)%chb2         = tmp_chb2
            NOAHMP36_struc(n)%noahmp36(t)%fpice        = tmp_fpice
            NOAHMP36_struc(n)%noahmp36(t)%sfcheadrt    = tmp_sfcheadrt
            NOAHMP36_struc(n)%noahmp36(t)%albd       = tmp_albd
            NOAHMP36_struc(n)%noahmp36(t)%albi       = tmp_albi  

            ![ 1] output variable: soil_temp (unit=K). ***  soil layer temperature
            soil_temp(1:NOAHMP36_struc(n)%nsoil) = NOAHMP36_struc(n)%noahmp36(t)%sstc(NOAHMP36_struc(n)%nsnow+1 : NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow)
            do i=1, NOAHMP36_struc(n)%nsoil
               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILTEMP, value = soil_temp(i), &
                    vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 2] output variable: snow_temp (unit=K). ***  snow layer temperature
            snow_temp(1:NOAHMP36_struc(n)%nsnow) = NOAHMP36_struc(n)%noahmp36(t)%sstc(1:NOAHMP36_struc(n)%nsnow) 
            do i=1, NOAHMP36_struc(n)%nsnow
                ! Test code to reset snow temperature to undefined
                ! when there is no corresponding snow layer - Mocko
                if ((i + abs(NOAHMP36_struc(n)%noahmp36(t)%isnow))     &
                    .le.NOAHMP36_struc(n)%nsnow) then
                   snow_temp(i) = LIS_rc%udef
                endif
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTPROF, value = snow_temp(i), &
                     vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
             end do
            
            ![ 3] output variable: sh2o (unit=m^3 m-3). ***  volumetric liquid soil moisture 
             do i=1, NOAHMP36_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMLIQFRAC, value = NOAHMP36_struc(n)%noahmp36(t)%sh2o(i), &
                     vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
             end do
            
            ![ 4] output variable: smc (unit=m^3 m-3 ). ***  volumetric soil moisture, ice + liquid 
            do i=1, NOAHMP36_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = NOAHMP36_struc(n)%noahmp36(t)%smc(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = NOAHMP36_struc(n)%noahmp36(t)%smc(i)*tmp_sldpth(i)*LIS_CONST_RHOFW,  &
                                                  vlevel=i, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 5] output variable: tah (unit=K  ). ***  canopy air temperature 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_TEMP, value = NOAHMP36_struc(n)%noahmp36(t)%tah, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 6] output variable: eah (unit=Pa  ). ***  canopy air vapor pressure 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_VP, value = NOAHMP36_struc(n)%noahmp36(t)%eah, &
                                              vlevel=1, unit="Pa", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 7] output variable: fwet (unit=-  ). ***  wetted or snowed fraction of canopy 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_WF, value = NOAHMP36_struc(n)%noahmp36(t)%fwet, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 8] output variable: canliq (unit=mm). ***  intercepted liquid water 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_INTL, value = NOAHMP36_struc(n)%noahmp36(t)%canliq, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 9] output variable: canice (unit=mm). ***  intercepted ice mass 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWEVEG, value = NOAHMP36_struc(n)%noahmp36(t)%canice, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 10] output variable: tv (unit=K ). ***  vegetation temperature 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGT, value = NOAHMP36_struc(n)%noahmp36(t)%tv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 11] output variable: tg (unit=K). ***  ground averaged temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GROUNDAVGT, value = NOAHMP36_struc(n)%noahmp36(t)%tg, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 12] output variable: isnow (unit=-). ***  actual number of snow layers 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOWN_NLAYER, value = -1.0*NOAHMP36_struc(n)%noahmp36(t)%isnow, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 13] output variable: z_snow (unit=m). ***  snow layer-bottom depth from snow surface
            do i=1, NOAHMP36_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW_LBDFSS, value = NOAHMP36_struc(n)%noahmp36(t)%zss(i), &
                                                  vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 14] output variable: z_soil (unit=m). ***  soil layer-bottom depth from snow surface
            do i=1, NOAHMP36_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOIL_LBDFSS, value = NOAHMP36_struc(n)%noahmp36(t)%zss(i+NOAHMP36_struc(n)%nsnow), &
                                                  vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 15] output variable: snowh (unit=m ). ***  snow height 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, value = NOAHMP36_struc(n)%noahmp36(t)%snowh, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 16] output variable: sneqv (unit=mm ). ***  snow water equivalent 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWE, value = NOAHMP36_struc(n)%noahmp36(t)%sneqv, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 17] output variable: snowice (unit=mm ). ***  snow-layer ice 
            do i=1, NOAHMP36_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWICE, value = NOAHMP36_struc(n)%noahmp36(t)%snowice(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 18] output variable: snowliq (unit=mm ). ***  snow-layer liquid water 
            do i=1, NOAHMP36_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQ, value = NOAHMP36_struc(n)%noahmp36(t)%snowliq(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            !snow density for each layer
            do i=1, NOAHMP36_struc(n)%nsnow

               bdsno = LIS_rc%udef
              
               if(-1.0*NOAHMP36_struc(n)%noahmp36(t)%zss(i).gt.0) then 
                  if(abs(NOAHMP36_struc(n)%noahmp36(t)%isnow).gt.1.and.i.gt.1) then 
                     bdsno = (NOAHMP36_struc(n)%noahmp36(t)%snowliq(i) + NOAHMP36_struc(n)%noahmp36(t)%snowice(i))/& 
                          (-1.0*(NOAHMP36_struc(n)%noahmp36(t)%zss(i) - & 
                          NOAHMP36_struc(n)%noahmp36(t)%zss(i-1)))
                  else
                     bdsno = (NOAHMP36_struc(n)%noahmp36(t)%snowliq(i) + NOAHMP36_struc(n)%noahmp36(t)%snowice(i))/& 
                          (-1.0*NOAHMP36_struc(n)%noahmp36(t)%zss(i))
                  endif
               endif

               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAYERSNOWDENSITY, value = bdsno,&
                    vlevel=i, unit="kg m-3", direction="-", surface_type = LIS_rc%lsm_index)
            enddo

            !snow depth for each layer, Rhae Sung 04/27/2018
            do i=1, NOAHMP36_struc(n)%nsnow

               layersd = LIS_rc%udef


               if(-1.0*NOAHMP36_struc(n)%noahmp36(t)%zss(i).gt.0) then
                  if(abs(NOAHMP36_struc(n)%noahmp36(t)%isnow).gt.1.and.i.gt.1) then
                     layersd =  (-1.0*(NOAHMP36_struc(n)%noahmp36(t)%zss(i) - &
                                       NOAHMP36_struc(n)%noahmp36(t)%zss(i-1)))
                  else
                     layersd =  (-1.0*NOAHMP36_struc(n)%noahmp36(t)%zss(i))
                  endif
               endif

               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAYERSNOWDEPTH, value = layersd,&
                    vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            enddo



            ![ 19] output variable: zwt (unit=m). ***  depth to water table 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WATERTABLED, value = NOAHMP36_struc(n)%noahmp36(t)%zwt, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 20] output variable: wa (unit=mm). ***  water storage in aquifer 
            !water storage in aquifer = ground water storage - David Mocko
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GWS, value = NOAHMP36_struc(n)%noahmp36(t)%wa, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 21] output variable: wt (unit=mm). ***  water in aquifer and saturated soil 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WT_AQUI_SATSOIL, value = NOAHMP36_struc(n)%noahmp36(t)%wt, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ! TWS should be SWE + CanopInt + Soil moisture + WA - David Mocko
            TWS_out = 0.0
            do i = 1,NOAHMP36_struc(n)%nsoil
               TWS_out = TWS_out +                                     &
                      (NOAHMP36_struc(n)%noahmp36(t)%smc(i)  *         &
                       tmp_sldpth(i)*LIS_CONST_RHOFW)
            enddo
            TWS_out =  NOAHMP36_struc(n)%noahmp36(t)%sneqv   +         &
                      (NOAHMP36_struc(n)%noahmp36(t)%canliq  +         &
                       NOAHMP36_struc(n)%noahmp36(t)%canice) +         &
                       NOAHMP36_struc(n)%noahmp36(t)%wa      + TWS_out
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TWS, value = TWS_out, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 22] output variable: wslake (unit=mm). ***  lake water storage 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKEWATER, value = NOAHMP36_struc(n)%noahmp36(t)%wslake, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 23] output variable: lfmass (unit=g/m2). ***  leaf mass 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LEAFMASS, value = NOAHMP36_struc(n)%noahmp36(t)%lfmass, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 24] output variable: rtmass (unit=g/m2 ). ***  mass of fine roots 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ROOTMASS, value = NOAHMP36_struc(n)%noahmp36(t)%rtmass, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 25] output variable: stmass (unit=g/m2 ). ***  stem mass 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_STEMMASS, value = NOAHMP36_struc(n)%noahmp36(t)%stmass, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 26] output variable: wood (unit=g/m2). ***  mass of wood including woody roots 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WOODMASS, value = NOAHMP36_struc(n)%noahmp36(t)%wood, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 27] output variable: stblcp (unit=g/m2). ***  stable carbon in deep soil 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CARBON_DEEPSOIL, value = NOAHMP36_struc(n)%noahmp36(t)%stblcp, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 28] output variable: fastcp (unit=g/m2 ). ***  short-lived carbon in shallow soil 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CARBON_SHALLOWSOIL, value = NOAHMP36_struc(n)%noahmp36(t)%fastcp, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 29] output variable: lai (unit=-). ***  leaf area index 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAI, value = NOAHMP36_struc(n)%noahmp36(t)%lai, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 30] output variable: sai (unit=- ). ***  stem area index 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAI, value = NOAHMP36_struc(n)%noahmp36(t)%sai, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 31] output variable: cm (unit=m s-1 ). ***  momentum drag coefficient 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CM, value = NOAHMP36_struc(n)%noahmp36(t)%cm, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 32] output variable: ch (unit=m s-1 ). ***  sensible heat exchange coefficient 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CH, value = NOAHMP36_struc(n)%noahmp36(t)%ch, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 33] output variable: tauss (unit=- ). ***  snow aging term 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWAGE, value = NOAHMP36_struc(n)%noahmp36(t)%tauss, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 34] output variable: smcwtd (unit=m^3 m-3). ***  soil water content between bottom of the soil and water table 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BETWEENWATER, value = NOAHMP36_struc(n)%noahmp36(t)%smcwtd, &
                                              vlevel=1, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 35] output variable: deeprech (unit=m). ***  recharge to the water table when groundwater is deep 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QRECTOGW, value = NOAHMP36_struc(n)%noahmp36(t)%deeprech, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 36] output variable: rech (unit=m). ***  recharge from the water table when groundwater is shallow 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QRECFROMGW, value = NOAHMP36_struc(n)%noahmp36(t)%rech, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 37] output variable: fsa (unit=W m-2). ***  total absorbed solar radiation 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWNET, value = NOAHMP36_struc(n)%noahmp36(t)%fsa, &
                                              vlevel=1, unit="W m-2", direction="DN", surface_type = LIS_rc%lsm_index)
            
            ![ 38] output variable: fsr (unit=W m-2 ). ***  total reflected solar radiation 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FSR, value = NOAHMP36_struc(n)%noahmp36(t)%fsr, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 39] output variable: fira (unit=W m-2 ). ***  total net longwave radiation to atmosphere 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWUP, value = NOAHMP36_struc(n)%noahmp36(t)%fira, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 40] output variable: fsh (unit=W m-2). ***  total sensible heat to atmosphere 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QH, value = NOAHMP36_struc(n)%noahmp36(t)%fsh, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 41] output variable: ssoil (unit=W m-2). ***  ground heat flux to soil 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QG, value = NOAHMP36_struc(n)%noahmp36(t)%ssoil, &
                                              vlevel=1, unit="W m-2", direction="DN", surface_type = LIS_rc%lsm_index)
            
            ![ 42] output variable: fcev (unit=W m-2). ***  canopy evaporative heat to atmosphere 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FCEV, value = NOAHMP36_struc(n)%noahmp36(t)%fcev, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 43] output variable: fgev (unit=W m-2  ). ***  ground evaporative heat to atmosphere 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FGEV, value = NOAHMP36_struc(n)%noahmp36(t)%fgev, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 44] output variable: fctr (unit=W m-2). ***  transpiration heat to atmosphere 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FCTR, value = NOAHMP36_struc(n)%noahmp36(t)%fctr, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 45] output variable: ecan (unit=kg m-2 s-1 ). ***  evaporation rate of canopy water 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ECANOP, value = NOAHMP36_struc(n)%noahmp36(t)%ecan, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 46] output variable: etran (unit=kg m-2 s-1 ). ***  transpiration rate 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TVEG, value = NOAHMP36_struc(n)%noahmp36(t)%etran, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 47] output variable: edir (unit=kg m-2 s-1 ). ***  direct evaporation rate from surface 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ESOIL, value = NOAHMP36_struc(n)%noahmp36(t)%edir, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 48] output variable: trad (unit=K  ). ***  surface radiative temperature 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RADT, value = NOAHMP36_struc(n)%noahmp36(t)%trad, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 49] output variable: tgb (unit=K). ***  bare ground surface temperature 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARESOILT, value = NOAHMP36_struc(n)%noahmp36(t)%tgb, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 50] output variable: tgv (unit=K). ***  vegetated ground surface temperature 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GROUNDVEGT, value = NOAHMP36_struc(n)%noahmp36(t)%tgv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 51] output variable: t2mv (unit=K). ***  2-m air temperature over vegetated part 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGE2MT, value = NOAHMP36_struc(n)%noahmp36(t)%t2mv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 52] output variable: t2mb (unit=K ). ***  2-m height air temperature 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARE2MT, value = NOAHMP36_struc(n)%noahmp36(t)%t2mb, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 53] output variable: q2v (unit=kg kg-1). ***  2-m specific humidity over vegetation 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGE2MQ2, value = NOAHMP36_struc(n)%noahmp36(t)%q2v, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 54] output variable: q2b (unit=kg kg-1). ***  2-m air specific humidity 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARE2MQ2, value = NOAHMP36_struc(n)%noahmp36(t)%q2b, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 55] output variable: runsrf (unit=kg m-2 s-1). ***  surface runoff 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS, value = NOAHMP36_struc(n)%noahmp36(t)%runsrf, &
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 56] output variable: runsub (unit=kg m-2 s-1 ). ***  baseflow (saturation excess) 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSB, value = NOAHMP36_struc(n)%noahmp36(t)%runsub, &
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 57] output variable: apar (unit=W m-2). ***  photosynthesis active energy by canopy 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_APAR, value = NOAHMP36_struc(n)%noahmp36(t)%apar, &
                                              vlevel=1, unit="W m-2", direction="IN", surface_type = LIS_rc%lsm_index)
            
            ![ 58] output variable: psn (unit=umol m-2 s-1 ). ***  total photosynthesis of CO2 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PSCO2, value = NOAHMP36_struc(n)%noahmp36(t)%psn, &
                                              vlevel=1, unit="umol m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)
            
            ![ 59] output variable: sav (unit=W m-2 ). ***  solar radiation absorbed by vegetation 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAV, value = NOAHMP36_struc(n)%noahmp36(t)%sav, &
                                              vlevel=1, unit="W m-2", direction="IN", surface_type = LIS_rc%lsm_index)
            
            ![ 60] output variable: sag (unit=W m-2 ). ***  solar radiation absorbed by ground 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAG, value = NOAHMP36_struc(n)%noahmp36(t)%sag, &
                                              vlevel=1, unit="W m-2", direction="IN", surface_type = LIS_rc%lsm_index)
            
            ![ 61] output variable: fsno (unit=-). ***  snow-cover fraction on the ground 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER, value = NOAHMP36_struc(n)%noahmp36(t)%fsno, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 62] output variable: nee (unit=g m-2 s-1 ). ***  net ecosystem exchange of CO2 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NEE, value = NOAHMP36_struc(n)%noahmp36(t)%nee, &
                                              vlevel=1, unit="g m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 63] output variable: gpp (unit=g m-2 s-1 ). ***  net instantaneous assimilation of carbon 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GPP, value = NOAHMP36_struc(n)%noahmp36(t)%gpp, &
                                              vlevel=1, unit="g m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)
            
            ![ 64] output variable: npp (unit=g m-2 s-1). ***  net primary productivity of carbon 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NPP, value = NOAHMP36_struc(n)%noahmp36(t)%npp, &
                                              vlevel=1, unit="g m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 65] output variable: fveg (unit=-). ***  green vegetation fraction 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GREENNESS, value = NOAHMP36_struc(n)%noahmp36(t)%fveg, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GREENNESS, value = NOAHMP36_struc(n)%noahmp36(t)%fveg*100.0, &
                                              vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
            
            if(tmp_albedo.lt.0) then 
               tmp_albedo = LIS_rc%udef
            endif
            ![ 66] output variable: albedo (unit=- ). ***  surface albedo 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO, &
                 value = tmp_albedo,&
                 vlevel=1, unit="-", direction="-", &
                 surface_type = LIS_rc%lsm_index)
            
            ![ 67] output variable: qsnbot (unit=kg m-2 s-1). ***  melting water out of snow bottom 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSM, value = NOAHMP36_struc(n)%noahmp36(t)%qsnbot, &
                                              vlevel=1, unit="kg m-2 s-1", direction="S2L", surface_type = LIS_rc%lsm_index)
            
            ![ 68] output variable: ponding (unit=mm). ***  surface ponding 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PONDING, value = NOAHMP36_struc(n)%noahmp36(t)%ponding, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 69] output variable: ponding1 (unit=mm). ***  surface ponding1 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PONDING1, value = NOAHMP36_struc(n)%noahmp36(t)%ponding1, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 70] output variable: ponding2 (unit=mm ). ***  surface ponding2 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PONDING2, value = NOAHMP36_struc(n)%noahmp36(t)%ponding2, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 71] output variable: rssun (unit=s m-1). ***  sunlit stomatal resistance 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RSSUN, value = NOAHMP36_struc(n)%noahmp36(t)%rssun, &
                                              vlevel=1, unit="s m-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 72] output variable: rssha (unit=s m-1). ***  shaded stomatal resistance 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RSSHA, value = NOAHMP36_struc(n)%noahmp36(t)%rssha, &
                                              vlevel=1, unit="s m-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 73] output variable: bgap (unit=-). ***  between-canopy gap fraction for beam 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BGAP, value = NOAHMP36_struc(n)%noahmp36(t)%bgap, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 74] output variable: wgap (unit=- ). ***  within-canopy gap fraction for beam 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WGAP, value = NOAHMP36_struc(n)%noahmp36(t)%wgap, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 75] output variable: chv (unit=m s-1). ***  sensible heat exchange coefficient over vegetated fraction 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHV, value = NOAHMP36_struc(n)%noahmp36(t)%chv, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 76] output variable: chb (unit=m s-1). ***  sensible heat exchange coefficient over bare-ground fraction 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHB, value = NOAHMP36_struc(n)%noahmp36(t)%chb, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 77] output variable: emissi (unit=- ). ***  surface emissivity 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EMISSFORC, value = NOAHMP36_struc(n)%noahmp36(t)%emissi, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 78] output variable: shg (unit=W m-2     ). ***  ground sensible heat  
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHG, value = NOAHMP36_struc(n)%noahmp36(t)%shg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 79] output variable: shc (unit=W m-2   ). ***  canopy sensible heat
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHC, value = NOAHMP36_struc(n)%noahmp36(t)%shc, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 80] output variable: shb (unit=W m-2     ). ***  bare ground sensible heat
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHB, value = NOAHMP36_struc(n)%noahmp36(t)%shb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 81] output variable: evg (unit=W m-2  ). ***  ground evaporation heat
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVG, value = NOAHMP36_struc(n)%noahmp36(t)%evg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 82] output variable: evb (unit=W m-2  ). ***  bare ground evaporation heat
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVB, value = NOAHMP36_struc(n)%noahmp36(t)%evb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 83] output variable: ghv (unit=W m-2 ). ***  vegetated ground heat flux
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GHV, value = NOAHMP36_struc(n)%noahmp36(t)%ghv, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 84] output variable: ghb (unit=W m-2 ). ***  bare ground heat flux 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GHB, value = NOAHMP36_struc(n)%noahmp36(t)%ghb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 85] output variable: irg (unit=W m-2 ). ***  vegeted ground net long wave radiation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRV, value = NOAHMP36_struc(n)%noahmp36(t)%irg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 86] output variable: irc (unit=W m-2 ). ***  canopy net long wave radiation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRC, value = NOAHMP36_struc(n)%noahmp36(t)%irc, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 87] output variable: irb (unit=W m-2 ). ***  bare ground net long wave radiation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRB, value = NOAHMP36_struc(n)%noahmp36(t)%irb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 88] output variable: tr (unit=W m-2 ). ***  transpiration heat
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_HTR, value = NOAHMP36_struc(n)%noahmp36(t)%tr, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 89] output variable: evc (unit=W m-2 ). ***  canopy evaporation heat
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_HEVC, value = NOAHMP36_struc(n)%noahmp36(t)%evc, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 90] output variable: chleaf (unit=m s-1). ***  leaf exchange coefficient 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHLEAF, value = NOAHMP36_struc(n)%noahmp36(t)%chleaf, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 91] output variable: chuc (unit=m s-1). ***  under canopy exchange coefficient  
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHUC, value = NOAHMP36_struc(n)%noahmp36(t)%chuc, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 92] output variable: chv2 (unit=m s-1). ***  sensible heat exchange coefficient over vegetated fraction 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHV2, value = NOAHMP36_struc(n)%noahmp36(t)%chv2, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 93] output variable: chb2 (unit=m s-1). ***  sensible heat exchange coefficient over bare-ground 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHB2, value = NOAHMP36_struc(n)%noahmp36(t)%chb2, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 94] output variable: fpice (unit=- ). ***  snow fraction in precipitation 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FPICE, value = NOAHMP36_struc(n)%noahmp36(t)%fpice, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            ![95] swnet
            !call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_SWNET,value=tmp_swdown * &
            !      (1.0-NOAHMP36_struc(n)%noahmp36(t)%albedo ),vlevel=1,unit="W m-2",&
            !      direction="DN",surface_type=LIS_rc%lsm_index)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWNET,value=NOAHMP36_struc(n)%noahmp36(t)%fsa,&
            !      vlevel=1,unit="W m-2", direction="DN",surface_type=LIS_rc%lsm_index)
            
            ![96] lwnet
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWNET,vlevel=1,  &
                  value=(-1.0 * NOAHMP36_struc(n)%noahmp36(t)%fira), &
                  unit="W m-2", direction="DN", surface_type=LIS_rc%lsm_index)
            
            ![97] average surface temperature 
            AvgSurfT_out = NOAHMP36_struc(n)%noahmp36(t)%fveg * NOAHMP36_struc(n)%noahmp36(t)%tv + &
                (1.0-NOAHMP36_struc(n)%noahmp36(t)%fveg) * NOAHMP36_struc(n)%noahmp36(t)%tgb
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AVGSURFT, value = AvgSurfT_out, &
                                              vlevel=1, unit="K", direction="-",  &
                                              surface_type = LIS_rc%lsm_index)
            ![98] latent heat flux  
            Qle_out = NOAHMP36_struc(n)%noahmp36(t)%fgev + NOAHMP36_struc(n)%noahmp36(t)%fcev + NOAHMP36_struc(n)%noahmp36(t)%fctr 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QLE, value=Qle_out, vlevel=1,unit="W m-2",&
                  direction="UP",surface_type=LIS_rc%lsm_index)
    
            ![99] Total evapotranspiration
            Evap_out = NOAHMP36_struc(n)%noahmp36(t)%ecan + NOAHMP36_struc(n)%noahmp36(t)%etran + NOAHMP36_struc(n)%noahmp36(t)%edir 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVAP, value=Evap_out, vlevel=1,unit="kg m-2 s-1",&
                  direction="UP",surface_type=LIS_rc%lsm_index)
            
            ! rainf/snowf added by David Mocko
            ! %fpice is the snow fraction of precipitation in Noah-MP
            !      and varies as a function of "snf" option
            Rainf_out = (tmp_prcp) * (1.0 - NOAHMP36_struc(n)%noahmp36(t)%fpice)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF, value=Rainf_out, &
                  vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF, value=Rainf_out*dt, &
                  vlevel=1,unit="kg m-2", direction="DN",surface_type=LIS_rc%lsm_index)
            Snowf_out = (tmp_prcp) * NOAHMP36_struc(n)%noahmp36(t)%fpice 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWF, value=Snowf_out, &
                  vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWF, value=Snowf_out*dt, &
                  vlevel=1,unit="kg m-2", direction="DN",surface_type=LIS_rc%lsm_index)

            ! canopint added by David Mocko
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPINT, value = &
                    (NOAHMP36_struc(n)%noahmp36(t)%canliq +                   &
                     NOAHMP36_struc(n)%noahmp36(t)%canice),                   &
                 vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ! subsnow added by David Mocko
            ! Note that sublimation is already included in the
            !   %edir term for evaporation from the surface
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SUBSNOW, value =  &
                 tmp_subsnow, vlevel=1,         &
                 unit="kg m-2 s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ! reset forcing variables to 0 for accumulation 
            NOAHMP36_struc(n)%noahmp36(t)%tair = 0.0
            NOAHMP36_struc(n)%noahmp36(t)%psurf = 0.0
            NOAHMP36_struc(n)%noahmp36(t)%wind_e = 0.0
            NOAHMP36_struc(n)%noahmp36(t)%wind_n = 0.0
            NOAHMP36_struc(n)%noahmp36(t)%qair = 0.0
            NOAHMP36_struc(n)%noahmp36(t)%swdown = 0.0
            NOAHMP36_struc(n)%noahmp36(t)%lwdown = 0.0
            NOAHMP36_struc(n)%noahmp36(t)%prcp = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        NOAHMP36_struc(n)%forc_count = 0 
    endif ! end of alarmCheck loop 

    deallocate( tmp_sldpth )
    deallocate( tmp_shdfac_monthly )
    deallocate( tmp_smceq )
    deallocate( tmp_sstc )
    deallocate( tmp_sh2o )
    deallocate( tmp_smc )
    deallocate( tmp_zss )
    deallocate( tmp_snowice )
    deallocate( tmp_snowliq )
end subroutine NoahMP36_main
